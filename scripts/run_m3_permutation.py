#!/usr/bin/env python
"""D4: gene-label permutation test for the M3 gene-only baseline (SPLIT-A and SPLIT-C).

For each permutation, the gene labels are shuffled among TRAINING variants only (np.random.default_rng(seed),
seeds 0..n_perm-1), the same smoothed per-gene rate is refit (fit_smoothed_rate / score_smoothed imported from
run_m3_gene_only.py, not reimplemented), and AUROC is evaluated on the REAL, unshuffled test set. SPLIT-B is
skipped: it is a constant prior by construction. CPU only.

Empirical p = (1 + count(null >= real)) / (1 + n_perm).
"""
from __future__ import annotations

import os

os.environ["CUDA_VISIBLE_DEVICES"] = ""

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from run_m3_gene_only import SPLIT_CFG, fit_smoothed_rate, score_smoothed  # noqa: E402  (reuse, do not reimplement)
from src.data.splits import SPLITS_DIR, apply_split, load_manifest  # noqa: E402
from src.evaluation.generate_results import FIG_DIR  # noqa: E402
from src.evaluation.generate_split_comparison import SPLIT_COLORS, SPLITS  # noqa: E402

M3_DIR = ROOT / "results" / "m3_gene_only"
OUT_DIR = M3_DIR / "d4_permutation"
SPLIT_LABEL = {sid: label for sid, _, label in SPLITS}


def log(msg: str) -> None:
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def run_split(sid: str, df: pd.DataFrame, alpha: float, n_perm: int) -> dict:
    t0 = time.time()
    manifest = load_manifest(SPLITS_DIR / SPLIT_CFG[sid][0])
    train, test = apply_split(df, manifest)
    y = test.label.to_numpy(dtype=int)

    rates, g = fit_smoothed_rate(train, alpha)
    real = float(roc_auc_score(y, score_smoothed(test, rates, g)))
    saved = json.loads((M3_DIR / f"split_{sid}" / "final_result.json").read_text())["metrics"]["roc_auc"]
    if abs(real - saved) > 1e-12:
        raise RuntimeError(f"SPLIT-{sid}: recomputed real AUROC {real} != saved M3 AUROC {saved}")
    log(f"SPLIT-{sid}: real M3 AUROC {real:.6f} (matches saved final_result.json); running {n_perm} permutations")

    null = np.empty(n_perm)
    base_genes = train.gene.to_numpy()
    for seed in range(n_perm):
        shuffled = train.copy()
        shuffled["gene"] = np.random.default_rng(seed).permutation(base_genes)  # labels stay put, gene identities shuffled
        r, gr = fit_smoothed_rate(shuffled, alpha)  # gr == g (label rate unchanged by shuffling genes)
        null[seed] = roc_auc_score(y, score_smoothed(test, r, gr))
        if (seed + 1) % 50 == 0:
            log(f"SPLIT-{sid}: {seed + 1}/{n_perm} permutations done")

    p_emp = float((1 + np.sum(null >= real)) / (1 + n_perm))
    out = {
        "split": sid, "alpha": alpha, "n_perm": n_perm, "seeds": [0, n_perm - 1], "real_auroc": real,
        "null_mean": float(null.mean()), "null_std": float(null.std(ddof=1)), "null_min": float(null.min()), "null_max": float(null.max()),
        "n_null_ge_real": int((null >= real).sum()), "empirical_p": p_emp, "null_aurocs": [float(v) for v in null],
        "manifest_content_hash": manifest["content_hash"],
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    (OUT_DIR / f"permutation_split_{sid}.json").write_text(json.dumps(out, indent=2))
    log(f"SPLIT-{sid}: null mean={out['null_mean']:.4f} std={out['null_std']:.4f} min={out['null_min']:.4f} max={out['null_max']:.4f} | "
        f"real={real:.4f} | null>=real: {out['n_null_ge_real']} | empirical p={p_emp:.5f} | {time.time() - t0:.1f}s")
    return out


def make_figure(results: dict[str, dict]) -> Path:
    fig, axes = plt.subplots(1, len(results), figsize=(6.6 * len(results), 5.0), squeeze=False)
    for ax, (sid, r) in zip(axes[0], results.items()):
        ax.hist(r["null_aurocs"], bins=25, color=SPLIT_COLORS[sid], alpha=0.85, edgecolor="white")
        ax.axvline(r["real_auroc"], color="black", linestyle="--", linewidth=1.8, label=f"Real M3 AUROC = {r['real_auroc']:.4f}")
        ax.axvline(0.5, color="dimgray", linestyle=":", linewidth=1.1, label="AUROC 0.5")
        ax.set_title(f"D4 Gene-Label Permutation Null — M3 ({SPLIT_LABEL[sid]})\n"
                     f"null mean {r['null_mean']:.3f}, sd {r['null_std']:.3f}, range {r['null_min']:.3f}–{r['null_max']:.3f}; "
                     f"empirical p = {r['empirical_p']:.4f}", fontsize=10)
        ax.set_xlabel("ROC-AUC on real test set (gene labels shuffled among train variants)")
        ax.set_ylabel("Permutations")
        ax.grid(axis="y", alpha=0.25)
        ax.legend(fontsize=8.5)
    fig.tight_layout()
    path = FIG_DIR / "d4_permutation_null_m3.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--alpha", type=float, default=1.0)
    ap.add_argument("--splits", nargs="+", choices=["A", "C"], default=["A", "C"])
    args = ap.parse_args()
    df = pd.read_parquet(ROOT / "data" / "processed" / "variants_100bp.parquet", columns=["chrom", "pos", "ref", "alt", "gene", "label"])
    log(f"load: {len(df):,} variants")
    # Expectation (comment only): the null should sit near 0.50. What actually comes out is printed below.
    results = {sid: run_split(sid, df, args.alpha, args.n_perm) for sid in args.splits}
    (OUT_DIR / "d4_summary.json").write_text(json.dumps(
        {sid: {k: v for k, v in r.items() if k != "null_aurocs"} for sid, r in results.items()}, indent=2))
    fig = make_figure(results)
    log(f"save: {fig.relative_to(ROOT)}, {OUT_DIR.relative_to(ROOT)}/")
    rows = [{"split": s, "real AUROC": r["real_auroc"], "null mean": r["null_mean"], "null std": r["null_std"],
             "null min": r["null_min"], "null max": r["null_max"], "n null>=real": r["n_null_ge_real"], "empirical p": r["empirical_p"]}
            for s, r in results.items()]
    print("\nD4 permutation summary (actual output):", flush=True)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.4f}"), flush=True)


if __name__ == "__main__":
    main()
