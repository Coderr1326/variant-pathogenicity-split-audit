#!/usr/bin/env python
"""Resumable local transformer training (dnabert2/nt) for unattended multi-hour runs.

run_experiment.py has NO mid-training checkpointing or resume support — it trains
straight through all epochs in memory and only writes one final best-state file at
the very end. An interruption (crash, OOM, power loss, laptop suspend gone wrong)
loses all progress. This script mirrors train_transformers_colab.ipynb's per-epoch
checkpoint+resume logic (originally built for Colab's session cutoff) so a local
unattended run is crash-safe too, using the exact same shared src/ functions
(load_transformer, compute_max_length, PairSet, metrics, seed) already exercised by
the compute gate — so the architecture/tokenization match what was gate-tested.

Usage:
    python scripts/train_transformer_local.py --model nt --data data/processed/variants_100bp.parquet

Resuming: just re-run the same command. It auto-detects and resumes from the latest
results/checkpoints/<model>/epoch_NN.pt if one exists. Pass --resume <path> to
override, or --restart to ignore existing checkpoints and start over from epoch 0.
"""
from __future__ import annotations
import argparse, copy, json, sys, time
from pathlib import Path

import numpy as np, pandas as pd, torch
from sklearn.model_selection import train_test_split
from torch import nn
from torch.utils.data import DataLoader
from transformers import get_linear_schedule_with_warmup

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
from src.training.run_experiment import PairSet, metrics, seed
from src.training.transformer_utils import load_transformer, compute_max_length
from src.data.splits import load_manifest, apply_split

CHECKPOINT_NAMES = {
    "dnabert2": "zhihan1996/DNABERT-2-117M",
    "nt": "InstaDeepAI/nucleotide-transformer-v2-100m-multi-species",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=["dnabert2", "nt"], required=True)
    ap.add_argument("--data", required=True)
    ap.add_argument("--window", type=int, default=100)
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--resume", default=None, help="checkpoint path to resume from (default: auto-detect latest)")
    ap.add_argument("--restart", action="store_true", help="ignore existing checkpoints, start from epoch 0")
    ap.add_argument("--split-manifest", default=None, help="path to a frozen split manifest (e.g. data/splits/split_manifest_b.json); "
                                                             "default omits this and uses the original random 80/20 split (SPLIT-A / S1), unchanged")
    ap.add_argument("--warmup-ratio", type=float, default=0.1,
                     help="fraction of total training steps spent on linear LR warmup from 0 to the base LR "
                          "(transformers.get_linear_schedule_with_warmup, stepped once per batch). "
                          "Set to 0 to disable and reproduce the old fixed-LR behavior exactly (pre-warmup-fix "
                          "SPLIT-A/SPLIT-C runs used no scheduler at all).")
    ap.add_argument("--fraction", type=float, default=1.0, help="subsample this fraction of train+val data (stratified by label); for smoke tests")
    ap.add_argument("--log-every", type=int, default=50, help="print per-step training loss every N steps")
    ap.add_argument("--tag", default="", help="extra suffix for checkpoint/metrics/prediction paths, to isolate one-off runs "
                                                "(e.g. smoke tests) from the real checkpoint history for a model/split")
    ap.add_argument("--seed", type=int, default=42,
                     help="RNG seed for shuffle order, classifier-head init and dropout. Does NOT affect the train/val split "
                          "(frozen manifest) or --fraction sampling. A non-default seed requires --tag to contain 'seed<N>' "
                          "so it cannot overwrite the seed-42 outputs.")
    args = ap.parse_args()
    if args.seed != 42 and f"seed{args.seed}" not in args.tag:
        ap.error(f"--seed {args.seed} requires --tag to contain 'seed{args.seed}' (e.g. --tag _seed{args.seed}) to avoid overwriting seed-42 outputs")
    seed(args.seed)

    df = pd.read_parquet(args.data)
    full_pairs = df.ref_sequence + df.mut_sequence
    full_longest = full_pairs.loc[full_pairs.str.len().nlargest(20).index].tolist()
    if args.fraction < 1.0:
        df = df.groupby("label", group_keys=False).sample(frac=args.fraction, random_state=42)
    split_tag = ""
    if args.split_manifest:
        manifest = load_manifest(args.split_manifest)
        split_tag = "_split" + manifest["split_id"].upper()
        train_df, val_df = apply_split(df, manifest)
    else:
        train_df, val_df = train_test_split(df, test_size=0.2, random_state=42, stratify=df.label)
    split_tag = split_tag + args.tag

    tokenizer, model = load_transformer(CHECKPOINT_NAMES[args.model])
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model.to(device)

    maxlen = compute_max_length(tokenizer, full_longest) if args.model == "nt" else 2 * args.window
    print("max_length used:", maxlen, file=sys.stderr)
    print("device:", device, "train:", len(train_df), "val:", len(val_df), "batch_size:", args.batch_size, file=sys.stderr)

    train_loader = DataLoader(PairSet(train_df, "transformer", tokenizer, maxlen), batch_size=args.batch_size, shuffle=True)
    val_loader = DataLoader(PairSet(val_df, "transformer", tokenizer, maxlen), batch_size=args.batch_size)

    counts = np.bincount(train_df.label, minlength=2)
    weights = torch.tensor(len(train_df) / (2 * np.maximum(counts, 1)), dtype=torch.float, device=device)
    loss_fn = nn.CrossEntropyLoss(weight=weights)
    optimizer = torch.optim.AdamW(model.parameters(), lr=2e-5, weight_decay=0.01)
    use_amp = device.type == "cuda"
    scaler = torch.amp.GradScaler("cuda", enabled=use_amp)

    total_steps = len(train_loader) * args.epochs
    scheduler = None
    if args.warmup_ratio > 0:
        warmup_steps = int(args.warmup_ratio * total_steps)
        scheduler = get_linear_schedule_with_warmup(optimizer, num_warmup_steps=warmup_steps, num_training_steps=total_steps)
        print(f"LR warmup: {warmup_steps}/{total_steps} steps ({args.warmup_ratio:.0%}), linear 0 -> 2e-5", file=sys.stderr)

    ckpt_dir = REPO_ROOT / "results/checkpoints" / (args.model + split_tag.lower())
    ckpt_dir.mkdir(parents=True, exist_ok=True)
    log_path = REPO_ROOT / "results/metrics" / f"{args.model}_{args.window}bp{split_tag}_epoch_times.csv"
    log_path.parent.mkdir(parents=True, exist_ok=True)

    def save_checkpoint(epoch, best_score, best_state):
        p = ckpt_dir / f"epoch_{epoch:02d}.pt"
        torch.save(
            {
                "epoch": epoch,
                "model": model.state_dict(),
                "optimizer": optimizer.state_dict(),
                "scaler": scaler.state_dict(),
                "scheduler": scheduler.state_dict() if scheduler is not None else None,
                "best_score": best_score,
                "best_state": best_state,
            },
            p,
        )
        return p

    start_epoch, best_score, best_state = 0, -1.0, None
    resume_path = args.resume
    if resume_path is None and not args.restart:
        existing = sorted(ckpt_dir.glob("epoch_*.pt"))
        if existing:
            resume_path = str(existing[-1])
    if resume_path:
        ck = torch.load(resume_path, map_location=device)
        model.load_state_dict(ck["model"])
        optimizer.load_state_dict(ck["optimizer"])
        scaler.load_state_dict(ck["scaler"])
        if scheduler is not None and ck.get("scheduler") is not None:
            scheduler.load_state_dict(ck["scheduler"])
        start_epoch, best_score, best_state = ck["epoch"], ck["best_score"], ck.get("best_state")
        print(f"Resumed from {resume_path} after epoch {start_epoch}, best_score so far {best_score:.4f}", file=sys.stderr)

    if start_epoch >= args.epochs:
        print(f"start_epoch ({start_epoch}) >= --epochs ({args.epochs}); nothing left to train.", file=sys.stderr)

    for epoch in range(start_epoch, args.epochs):
        epoch_start = time.time()
        model.train()
        for step, (x, y) in enumerate(train_loader, start=1):
            y = y.to(device)
            optimizer.zero_grad(set_to_none=True)
            with torch.autocast(device_type="cuda", dtype=torch.float16, enabled=use_amp):
                loss = loss_fn(model(**{k: v.to(device) for k, v in x.items()}), y)
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            scaler.step(optimizer)
            scaler.update()
            if scheduler is not None:
                scheduler.step()
            if step % args.log_every == 0 or step == 1:
                print(json.dumps({
                    "epoch": epoch + 1, "step": step, "loss": loss.item(),
                    "lr": optimizer.param_groups[0]["lr"],
                }))
                sys.stdout.flush()

        model.eval()
        ys, ps, probs = [], [], []
        with torch.no_grad():
            for x, y in val_loader:
                q = torch.softmax(model(**{k: v.to(device) for k, v in x.items()}), 1)[:, 1].cpu().numpy()
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
            [{"model": args.model, "epoch": epoch + 1, "start_time": time.ctime(epoch_start), "end_time": time.ctime(end), "duration_seconds": end - epoch_start}]
        ).to_csv(log_path, index=False, mode="a", header=not log_path.exists())
        print(json.dumps({"epoch": epoch + 1, "weighted_f1": score, "checkpoint": str(ckpt)}))
        sys.stdout.flush()

    # final eval + save, same output format as the notebook / run_experiment.py
    model.load_state_dict(best_state)
    model.eval()
    ys, ps, probs = [], [], []
    with torch.no_grad():
        for x, y in val_loader:
            q = torch.softmax(model(**{k: v.to(device) for k, v in x.items()}), 1)[:, 1].cpu().numpy()
            probs.extend(q)
            ps.extend((q >= 0.5).astype(int))
            ys.extend(y.numpy())
    result = metrics(ys, ps, probs)
    pred_dir = REPO_ROOT / "results/predictions"
    pred_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"label": ys, "prediction": ps, "probability_pathogenic": probs}).to_parquet(
        pred_dir / f"{args.model}_{args.window}bp{split_tag}.parquet", index=False
    )
    (REPO_ROOT / "results/metrics" / f"{args.model}_{args.window}bp{split_tag}_final.json").write_text(
        json.dumps({"model": args.model, "window": args.window, "metrics": result,
                    "seed": args.seed, "warmup_ratio": args.warmup_ratio}, indent=2)
    )
    print(json.dumps({"final_result": result}))


if __name__ == "__main__":
    main()
