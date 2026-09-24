#!/usr/bin/env python
"""Validated join of saved test-set predictions to genes.

Saved prediction files hold (label, prediction, probability_pathogenic) with NO variant key, so genes are attached by
row order to the test rows reconstructed from the frozen manifest + data/processed/variants_100bp.parquet. A file is
only returned if its label column matches the reconstructed labels for EVERY row (AssertionError otherwise), and
its 'prediction' column equals (probability >= 0.5).

Row-order provenance: every B/C file and every M3 file is in manifest (apply_split) order. The transformers' SPLIT-A
files were written by the inline train_test_split(test_size=0.2, random_state=42, stratify=label) val order, which is a
shuffled order, not manifest order; that order is reconstructed the same way and checked against the manifest membership.

No model inference happens here. Classical models (CNN, BiLSTM, CNN+BiLSTM, Ensemble) have no saved probabilities and are
listed in SKIPPED_MODELS. CPU only.

    python scripts/pred_loader.py      # prints the join-check table and the SPLIT-B gene check
"""
from __future__ import annotations

import os

os.environ["CUDA_VISIBLE_DEVICES"] = ""

import json
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from src.data.splits import SPLITS_DIR, apply_split, load_manifest, variant_ids  # noqa: E402

DATA = ROOT / "data" / "processed" / "variants_100bp.parquet"
MANIFEST_FILE = {"A": "split_manifest_a.json", "B": "split_manifest_b.json", "C": "split_manifest_c.json"}
SUFFIX = {"A": "", "B": "_splitB", "C": "_splitC"}
EXPECTED_N = {"A": 33_786, "B": 33_780, "C": 33_536}
MODELS = ["m3", "nt", "dnabert2"]
MODEL_LABEL = {"m3": "M3 (gene-only)", "nt": "NT", "dnabert2": "DNABERT-2"}
SKIPPED_MODELS = {m: "no saved probabilities (only weights + history); predictions are not regenerated"
                  for m in ["cnn", "bilstm", "cnn_bilstm", "ensemble"]}
DNABERT2_B_SEEDS = (42, 1, 2, 3)  # seed 42 = original run; 1/2/3 = warmup 0.1. The warmup-0 run (collapsed) is never exposed.


def _paths(model: str, split: str, seed: int | None) -> tuple[Path, Path]:
    """(predictions parquet, final-result json) for a run."""
    if model == "m3":
        d = ROOT / "results" / "m3_gene_only" / f"split_{split}"
        return d / "predictions.parquet", d / "final_result.json"
    stem = f"{model}_100bp{SUFFIX[split]}"
    if model == "dnabert2" and seed not in (None, 42):
        if split != "B" or seed not in DNABERT2_B_SEEDS:
            raise FileNotFoundError(f"no DNABERT-2 seed-{seed} run exists for SPLIT-{split}")
        stem += f"_seed{seed}"
    elif seed not in (None, 42) and model != "dnabert2":
        raise FileNotFoundError(f"{model} has no seed runs")
    return ROOT / "results" / "predictions" / f"{stem}.parquet", ROOT / "results" / "metrics" / f"{stem}_final.json"


@lru_cache(maxsize=1)
def _data() -> pd.DataFrame:
    df = pd.read_parquet(DATA, columns=["chrom", "pos", "ref", "alt", "gene", "label"])
    df["vid"] = variant_ids(df)
    return df


@lru_cache(maxsize=None)
def _manifest_split(split: str):
    manifest = load_manifest(SPLITS_DIR / MANIFEST_FILE[split])  # raises if the content hash does not match
    train, test = apply_split(_data(), manifest)
    return manifest, train.reset_index(drop=True), test.reset_index(drop=True)


def _order_kind(model: str, split: str) -> str:
    return "inline" if (split == "A" and model in ("nt", "dnabert2")) else "manifest"


@lru_cache(maxsize=None)
def test_frame(split: str, order: str = "manifest") -> pd.DataFrame:
    """Reconstructed test rows (vid, gene, label) in the given row order."""
    _, _, test = _manifest_split(split)
    if order == "manifest":
        return test[["vid", "gene", "label"]].copy()
    if order == "inline" and split == "A":
        df = _data()
        _, val = train_test_split(df, test_size=0.2, random_state=42, stratify=df.label)
        if set(val.vid) != set(test.vid):
            raise AssertionError("SPLIT-A: inline train_test_split val set differs from the frozen manifest's val set")
        return val.reset_index(drop=True)[["vid", "gene", "label"]].copy()
    raise ValueError(f"order '{order}' is not defined for SPLIT-{split}")


def load_predictions(model: str, split: str, seed: int | None = None) -> pd.DataFrame:
    """DataFrame(gene, label, score) for one run, validated row by row. seed None/42 = the original run."""
    if model in SKIPPED_MODELS:
        raise FileNotFoundError(f"{model}: {SKIPPED_MODELS[model]}")
    pred_path, _ = _paths(model, split, seed)
    if not pred_path.exists():
        raise FileNotFoundError(f"missing prediction file: {pred_path.relative_to(ROOT)}")
    pred = pd.read_parquet(pred_path)
    frame = test_frame(split, _order_kind(model, split))
    tag = f"{model} SPLIT-{split} seed={seed}"
    if len(pred) != len(frame) or len(pred) != EXPECTED_N[split]:
        raise AssertionError(f"{tag}: {len(pred)} prediction rows vs {len(frame)} reconstructed test rows (expected {EXPECTED_N[split]})")
    bad = int((pred["label"].to_numpy(dtype=int) != frame["label"].to_numpy(dtype=int)).sum())
    if bad:
        raise AssertionError(f"{tag}: label mismatch in {bad} of {len(pred)} rows -- row order does not match the reconstructed test set")
    score = pred["probability_pathogenic"].to_numpy(dtype=float)
    if int(((score >= 0.5).astype(int) != pred["prediction"].to_numpy(dtype=int)).sum()):
        raise AssertionError(f"{tag}: 'prediction' column is not (probability >= 0.5)")
    return pd.DataFrame({"gene": frame["gene"].to_numpy(), "label": frame["label"].to_numpy(dtype=int), "score": score})


def available_runs() -> list[tuple[str, str, int | None]]:
    """Every (model, split, seed) with a prediction file that this loader exposes."""
    runs = []
    for model in MODELS:
        for split in "ABC":
            if model == "dnabert2" and split == "B":
                seeds: list[int | None] = list(DNABERT2_B_SEEDS)
            else:
                seeds = [None]
            for seed in seeds:
                if _paths(model, split, seed)[0].exists():
                    runs.append((model, split, seed))
    return runs


def run_name(model: str, seed: int | None) -> str:
    return MODEL_LABEL[model] + (f" (seed {seed})" if (model == "dnabert2" and seed is not None) else "")


def join_check_table() -> pd.DataFrame:
    """Per run: validated join + pooled AUROC recomputed from the joined data vs. roc_auc in the run's final JSON (4 decimals)."""
    rows = []
    for model, split, seed in available_runs():
        row = {"run": run_name(model, seed), "split": split, "n_rows": np.nan, "labels_match": "", "pooled_auroc_recomputed": np.nan,
               "auroc_in_final_json": np.nan, "match_4dp": "", "status": ""}
        try:
            d = load_predictions(model, split, seed)
            _, final_path = _paths(model, split, seed)
            saved = float(json.loads(final_path.read_text())["metrics"]["roc_auc"])
            got = float(roc_auc_score(d["label"], d["score"]))
            row.update(n_rows=len(d), labels_match="all rows", pooled_auroc_recomputed=got, auroc_in_final_json=saved,
                       match_4dp="yes" if round(got, 4) == round(saved, 4) else "NO", status="OK" if round(got, 4) == round(saved, 4) else "FAILED (AUROC)")
        except (AssertionError, FileNotFoundError, ValueError) as e:
            row.update(labels_match="NO", status=f"FAILED: {e}")
        rows.append(row)
    return pd.DataFrame(rows)


def splitb_gene_check() -> dict:
    manifest, train, test = _manifest_split("B")
    test_genes, train_genes = sorted(set(test.gene)), set(train.gene)
    seen = sorted(set(test_genes) & train_genes)
    if len(test_genes) != 14 or seen:
        raise AssertionError(f"SPLIT-B gene check failed: {len(test_genes)} test genes, {len(seen)} seen in train {seen}")
    return {"test_genes": test_genes, "n_test_genes": len(test_genes), "seen_in_train": seen, "manifest_content_hash": manifest["content_hash"]}


if __name__ == "__main__":
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 120)
    print("skipped (no saved probabilities):", ", ".join(SKIPPED_MODELS), "\n")
    t = join_check_table()
    print(t.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    g = splitb_gene_check()
    print(f"\nSPLIT-B: {g['n_test_genes']} test genes, {len(g['seen_in_train'])} seen in train; genes: {', '.join(g['test_genes'])}")
    if (t["status"] != "OK").any():
        sys.exit("join check FAILED for at least one file (see table)")
