"""
Local finish-up script: resume DNABERT-2 from epoch 9 checkpoint (downloaded from
Colab/Drive) and complete epoch 10 on the local GPU. Already run on 2026-08-16;
epoch_10.pt and the final DNABERT-2 metrics in this repo came from this script.

ClassifierHead below (mean pooling, `encoder`/`classifier` state_dict keys) mirrors
the DNABERT-2 training cell as it existed on 2026-08-16 -- not the current
notebooks/train_transformers_colab.ipynb, which was rewritten on 2026-09-06 to
share src/training/transformer_utils.py's SequenceWrapper (CLS-token pooling,
`m`/`h` keys) for both models. That rewrite happened after DNABERT-2 finished and
was only ever exercised for Nucleotide Transformer, so DNABERT-2's checkpoint
lineage never used CLS pooling. Re-running DNABERT-2 through today's notebook would
produce a differently-keyed checkpoint incompatible with epoch_01-10.pt here; use
ClassifierHead (as below), not SequenceWrapper, if resuming this lineage again. See
REPRODUCTION_MANIFEST.yaml's checkpoints.dnabert2.pooling_provenance for details.

Before running:
1. Place the downloaded checkpoint at:
   results/checkpoints/dnabert2/epoch_09.pt
2. Run from the repo root: python scripts/finish_dnabert2_locally.py
"""

import copy
import json
import os
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.model_selection import train_test_split
from torch import nn
from torch.utils.data import DataLoader
from transformers import AutoConfig, AutoModel, AutoTokenizer

REPO_ROOT = Path(__file__).resolve().parent.parent  # adjust if script isn't in scripts/
sys.path.insert(0, str(REPO_ROOT))
from src.training.run_experiment import PairSet, metrics, seed

WINDOW = 100
MODEL = "dnabert2"
CHECKPOINT_NAME = "zhihan1996/DNABERT-2-117M"

seed(42)
DATA = REPO_ROOT / "data/processed" / f"variants_{WINDOW}bp.parquet"
OUT = REPO_ROOT

df = pd.read_parquet(DATA)
train_df, val_df = train_test_split(
    df, test_size=0.2, random_state=42, stratify=df.label
)

tokenizer = AutoTokenizer.from_pretrained(CHECKPOINT_NAME, trust_remote_code=True)
if tokenizer.pad_token is None:
    tokenizer.pad_token = tokenizer.eos_token if tokenizer.eos_token else "[PAD]"

config = AutoConfig.from_pretrained(CHECKPOINT_NAME, trust_remote_code=True)
encoder = AutoModel.from_pretrained(
    CHECKPOINT_NAME,
    config=config,
    trust_remote_code=True,
    low_cpu_mem_usage=False,
)


class ClassifierHead(nn.Module):
    def __init__(self, encoder, hidden_size):
        super().__init__()
        self.encoder = encoder
        self.classifier = nn.Linear(hidden_size, 2)

    def forward(self, input_ids, attention_mask, **kwargs):
        out = self.encoder(input_ids=input_ids, attention_mask=attention_mask)
        hidden = out[0] if isinstance(out, tuple) else out.last_hidden_state
        mask = attention_mask.unsqueeze(-1).to(hidden.dtype)
        pooled = (hidden * mask).sum(1) / mask.sum(1).clamp(min=1e-6)
        return self.classifier(pooled)


hidden_size = getattr(config, "hidden_size", None) or getattr(config, "d_model", None)
model = ClassifierHead(encoder, hidden_size)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model.to(device)
print("device:", device)

# same Triton flash-attention patch as the Colab notebook
patched = 0
for _name, _mod in list(sys.modules.items()):
    if _mod is not None and ("flash_attn_triton" in _name or "bert_layers" in _name):
        if hasattr(_mod, "flash_attn_qkvpacked_func"):
            _mod.flash_attn_qkvpacked_func = None
            patched += 1
print(
    f"Patched {patched} module(s) to disable incompatible Triton flash-attention kernel."
)

maxlen = 2 * WINDOW
train_loader = DataLoader(
    PairSet(train_df, "transformer", tokenizer, maxlen), batch_size=16, shuffle=True
)
val_loader = DataLoader(
    PairSet(val_df, "transformer", tokenizer, maxlen), batch_size=16
)

counts = np.bincount(train_df.label, minlength=2)
weights = torch.tensor(
    len(train_df) / (2 * np.maximum(counts, 1)), dtype=torch.float, device=device
)
loss_fn = nn.CrossEntropyLoss(weight=weights)
optimizer = torch.optim.AdamW(model.parameters(), lr=2e-5, weight_decay=0.01)
use_amp = device.type == "cuda"
scaler = torch.amp.GradScaler("cuda", enabled=use_amp)

CKPT_DIR = OUT / "results/checkpoints" / MODEL
LOG = OUT / "results/metrics" / f"{MODEL}_{WINDOW}bp_epoch_times.csv"


def save_checkpoint(epoch, best_score, best_state):
    p = CKPT_DIR / f"epoch_{epoch:02d}.pt"
    torch.save(
        {
            "epoch": epoch,
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "scaler": scaler.state_dict(),
            "best_score": best_score,
            "best_state": best_state,
        },
        p,
    )
    return p


# --- resume from the checkpoint you downloaded from Drive ---
resume_path = CKPT_DIR / "epoch_09.pt"
assert resume_path.exists(), (
    f"Expected checkpoint at {resume_path} — did you download and place it?"
)
ck = torch.load(resume_path, map_location=device)
model.load_state_dict(ck["model"])
optimizer.load_state_dict(ck["optimizer"])
scaler.load_state_dict(ck["scaler"])
start_epoch = ck["epoch"]  # 9
best_score = ck["best_score"]  # 0.8580...
best_state = ck.get("best_state")
print(f"Resumed after epoch {start_epoch}, best_score so far: {best_score:.4f}")

# --- run only the remaining epoch(s) ---
for epoch in range(start_epoch, 10):
    epoch_start = time.time()
    model.train()
    for x, y in train_loader:
        y = y.to(device)
        optimizer.zero_grad(set_to_none=True)
        with torch.autocast(device_type="cuda", dtype=torch.float16, enabled=use_amp):
            loss = loss_fn(model(**{k: v.to(device) for k, v in x.items()}), y)
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
        scaler.step(optimizer)
        scaler.update()

    model.eval()
    ys = []
    ps = []
    probs = []
    with torch.no_grad():
        for x, y in val_loader:
            q = (
                torch.softmax(model(**{k: v.to(device) for k, v in x.items()}), 1)[:, 1]
                .cpu()
                .numpy()
            )
            probs.extend(q)
            ps.extend((q >= 0.5).astype(int))
            ys.extend(y.numpy())
    m = metrics(ys, ps, probs)
    score = m["weighted_f1"]
    if score > best_score:
        best_score = score
        best_state = copy.deepcopy(model.state_dict())

    ckpt = save_checkpoint(epoch + 1, best_score, best_state)
    end = time.time()
    pd.DataFrame(
        [
            {
                "model": MODEL,
                "epoch": epoch + 1,
                "start_time": time.ctime(epoch_start),
                "end_time": time.ctime(end),
                "duration_seconds": end - epoch_start,
            }
        ]
    ).to_csv(LOG, index=False, mode="a", header=not LOG.exists())
    print(
        json.dumps({"epoch": epoch + 1, "weighted_f1": score, "checkpoint": str(ckpt)})
    )

# --- final eval + save (same output format as Colab notebook) ---
model.load_state_dict(best_state)
model.eval()
ys = []
ps = []
probs = []
with torch.no_grad():
    for x, y in val_loader:
        q = (
            torch.softmax(model(**{k: v.to(device) for k, v in x.items()}), 1)[:, 1]
            .cpu()
            .numpy()
        )
        probs.extend(q)
        ps.extend((q >= 0.5).astype(int))
        ys.extend(y.numpy())
result = metrics(ys, ps, probs)
pred = OUT / "results/predictions"
pred.mkdir(parents=True, exist_ok=True)
pd.DataFrame(
    {"label": ys, "prediction": ps, "probability_pathogenic": probs}
).to_parquet(pred / f"{MODEL}_{WINDOW}bp.parquet", index=False)
(OUT / "results/metrics" / f"{MODEL}_{WINDOW}bp_final.json").write_text(
    json.dumps({"model": MODEL, "window": WINDOW, "metrics": result}, indent=2)
)
print("FINAL RESULT:", result)
