#!/usr/bin/env python
"""D2 / M3: gene-identity-only baseline. Sees NO sequence -- only the variant's gene.

Model: smoothed per-gene pathogenic rate fitted on TRAIN ONLY,
    p_gene = (pathogenic_count + alpha * global_train_rate) / (n_variants + alpha)
Genes unseen in train fall back to the global train rate. A one-hot-gene
LogisticRegression (train only) is fitted as a sanity check and its AUROC is logged next to
the smoothed-rate AUROC. Predict pathogenic at p >= 0.5, same threshold as the other models.

Splits come from the frozen manifests (src.data.splits); metrics come from
run_experiment.metrics(), so numbers are computed exactly like the other models'. Output JSON
follows the existing final_result schema {model, window, metrics}; ``window`` is kept only for
schema parity (the model ignores sequence). CPU-only.

Usage:
    PYTHONUNBUFFERED=1 python scripts/run_m3_gene_only.py [--splits A B C] [--alpha 1.0]
"""
from __future__ import annotations

import os

os.environ["CUDA_VISIBLE_DEVICES"] = ""  # M3 is CPU-only; never touch the GPU (a transformer job may be running)

import argparse
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder

REPO_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_ROOT))
from src.data.splits import SPLITS_DIR, apply_split, load_manifest, variant_ids  # noqa: E402
from src.training.run_experiment import metrics, seed  # noqa: E402

DATA = REPO_ROOT / "data" / "processed" / "variants_100bp.parquet"
SPLIT_CFG = {  # split id -> (manifest file, results-file suffix used by the other models)
    "A": ("split_manifest_a.json", ""),
    "B": ("split_manifest_b.json", "_splitB"),
    "C": ("split_manifest_c.json", "_splitC"),
}
COMPARE_MODELS = ["dnabert2", "nt"]
AUROC_AGREEMENT_TOL = 0.01


def log(msg: str) -> None:
    print(f"[{datetime.now():%H:%M:%S}] {msg}", flush=True)


def fit_smoothed_rate(train: pd.DataFrame, alpha: float) -> tuple[pd.Series, float]:
    global_rate = float(train.label.mean())
    g = train.groupby("gene").label.agg(["sum", "count"])
    rates = (g["sum"] + alpha * global_rate) / (g["count"] + alpha)
    return rates, global_rate


def score_smoothed(df: pd.DataFrame, rates: pd.Series, global_rate: float) -> np.ndarray:
    return df.gene.map(rates).fillna(global_rate).to_numpy(dtype=float)


def score_onehot_lr(train: pd.DataFrame, test: pd.DataFrame) -> np.ndarray:
    enc = OneHotEncoder(handle_unknown="ignore")  # unseen gene -> all-zero row -> intercept-only prediction
    x_train = enc.fit_transform(train[["gene"]])
    clf = LogisticRegression(max_iter=1000, random_state=42).fit(x_train, train.label)
    return clf.predict_proba(enc.transform(test[["gene"]]))[:, 1]


def load_reference(model: str, suffix: str) -> dict | None:
    path = REPO_ROOT / "results" / "metrics" / f"{model}_100bp{suffix}_final.json"
    return json.loads(path.read_text())["metrics"] if path.exists() else None


def run_split(sid: str, df: pd.DataFrame, alpha: float, out_root: Path) -> dict:
    manifest_file, suffix = SPLIT_CFG[sid]
    t0 = time.time()
    log(f"SPLIT-{sid}: load -- reading manifest {manifest_file}")
    manifest = load_manifest(SPLITS_DIR / manifest_file)  # raises on content_hash mismatch
    train, test = apply_split(df, manifest)
    if (len(train), len(test)) != (manifest["n_train"], manifest["n_val"]):
        raise RuntimeError(f"SPLIT-{sid}: split sizes {(len(train), len(test))} != manifest {(manifest['n_train'], manifest['n_val'])}")
    if sid == "A":
        # SPLIT-A's manifest is an audit freeze; run_experiment.py recomputes the split inline. Prove they agree.
        inline_tr, _ = train_test_split(df, test_size=0.2, random_state=42, stratify=df.label)
        if set(variant_ids(inline_tr)) != set(variant_ids(train)):
            raise RuntimeError("SPLIT-A: frozen manifest differs from run_experiment.py's inline train_test_split")
        log("SPLIT-A: frozen manifest matches run_experiment.py's inline split (seed 42, stratified)")
    log(f"SPLIT-{sid}: manifest hash {manifest['content_hash']} | train={len(train):,} test={len(test):,}")

    train_genes, test_genes = set(train.gene), set(test.gene)
    seen = test_genes & train_genes
    seen_frac = float(test.gene.isin(train_genes).mean())
    diag = {
        "n_train_genes": len(train_genes), "n_test_genes": len(test_genes),
        "n_test_genes_seen_in_train": len(seen), "frac_test_variants_gene_seen_in_train": seen_frac,
    }
    log(f"SPLIT-{sid}: test genes={len(test_genes)} seen-in-train={len(seen)} | fraction of test variants with seen gene={seen_frac:.6f}")
    if sid == "B" and (len(seen) != 0 or seen_frac != 0.0):
        raise RuntimeError(f"SPLIT-B is not gene-disjoint: {len(seen)} test genes appear in train ({sorted(seen)}); manifest is invalid")

    log(f"SPLIT-{sid}: fit -- smoothed gene rate (alpha={alpha}) + one-hot LogisticRegression, train only")
    rates, global_rate = fit_smoothed_rate(train, alpha)
    log(f"SPLIT-{sid}: global train pathogenic rate={global_rate:.6f} | genes fitted={len(rates)}")

    log(f"SPLIT-{sid}: evaluate")
    y = test.label.to_numpy(dtype=int)
    score = score_smoothed(test, rates, global_rate)
    pred = (score >= 0.5).astype(int)
    m = metrics(y.tolist(), pred.tolist(), score.tolist())
    lr_score = score_onehot_lr(train, test)
    lr_auroc = float(roc_auc_score(y, lr_score))
    delta = abs(lr_auroc - m["roc_auc"])
    log(f"SPLIT-{sid}: acc={m['accuracy']:.4f} wF1={m['weighted_f1']:.4f} AUROC={m['roc_auc']:.4f} "
        f"path_recall={m['pathogenic_recall']:.4f} | tn={m['tn']} fp={m['fp']} fn={m['fn']} tp={m['tp']}")
    log(f"SPLIT-{sid}: AUROC smoothed-rate={m['roc_auc']:.6f} vs one-hot LR={lr_auroc:.6f} (|diff|={delta:.6f})"
        + ("" if delta <= AUROC_AGREEMENT_TOL else f"  ** WARNING: exceeds {AUROC_AGREEMENT_TOL} **"))
    chance_floor = bool(np.unique(score).size == 1)
    if chance_floor:
        majority = max(float((y == 0).mean()), float((y == 1).mean()))
        log(f"SPLIT-{sid}: CHANCE FLOOR -- every test score is the constant global prior {score[0]:.6f}; "
            f"AUROC={m['roc_auc']} (expected exactly 0.5), accuracy={m['accuracy']:.6f} vs majority-class rate {majority:.6f}")

    log(f"SPLIT-{sid}: save")
    out_dir = out_root / f"split_{sid}"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "final_result.json").write_text(json.dumps({"model": "m3_gene_only", "window": 100, "metrics": m}, indent=2))
    pd.DataFrame({"label": y, "prediction": pred, "probability_pathogenic": score}).to_parquet(out_dir / "predictions.parquet", index=False)
    diag.update({
        "manifest_content_hash": manifest["content_hash"], "n_train": len(train), "n_test": len(test),
        "global_train_rate": global_rate, "alpha": alpha, "chance_floor": chance_floor,
        "lr_onehot_auroc": lr_auroc, "auroc_abs_diff_smoothed_vs_lr": delta,
    })
    (out_dir / "diagnostics.json").write_text(json.dumps(diag, indent=2))
    log(f"SPLIT-{sid}: done in {time.time() - t0:.1f}s -> {out_dir.relative_to(REPO_ROOT)}")
    return {"metrics": m, "diagnostics": diag}


def verify_saved(sid: str, out_root: Path) -> None:
    """Re-read what was written and re-derive everything independently of the in-memory objects."""
    out_dir = out_root / f"split_{sid}"
    saved = json.loads((out_dir / "final_result.json").read_text())["metrics"]
    p = pd.read_parquet(out_dir / "predictions.parquet")
    y, pr = p.label.astype(int), p.prediction.astype(int)
    cm = {"tn": int(((y == 0) & (pr == 0)).sum()), "fp": int(((y == 0) & (pr == 1)).sum()),
          "fn": int(((y == 1) & (pr == 0)).sum()), "tp": int(((y == 1) & (pr == 1)).sum())}
    for k, v in cm.items():
        if saved[k] != v:
            raise RuntimeError(f"SPLIT-{sid}: saved JSON {k}={saved[k]} != recomputed from predictions {v}")
    if abs(roc_auc_score(y, p.probability_pathogenic) - saved["roc_auc"]) > 1e-12:
        raise RuntimeError(f"SPLIT-{sid}: saved AUROC does not match recomputed AUROC")
    if sid == "B":
        d = json.loads((out_dir / "diagnostics.json").read_text())
        if saved["roc_auc"] != 0.5 or d["frac_test_variants_gene_seen_in_train"] != 0.0:
            raise RuntimeError(f"SPLIT-B chance-floor check failed: AUROC={saved['roc_auc']}, seen fraction={d['frac_test_variants_gene_seen_in_train']}")
    log(f"SPLIT-{sid}: verify OK -- confusion counts, AUROC re-derived from saved predictions match saved JSON")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--splits", nargs="+", choices=list(SPLIT_CFG), default=list(SPLIT_CFG))
    ap.add_argument("--alpha", type=float, default=1.0, help="smoothing strength toward the global train rate")
    ap.add_argument("--data", default=str(DATA))
    ap.add_argument("--out-dir", default=str(REPO_ROOT / "results" / "m3_gene_only"))
    args = ap.parse_args()
    seed(42)
    out_root = Path(args.out_dir)

    log(f"load: reading {Path(args.data).relative_to(REPO_ROOT) if Path(args.data).is_relative_to(REPO_ROOT) else args.data} (metadata columns only; no sequence)")
    df = pd.read_parquet(args.data, columns=["chrom", "pos", "ref", "alt", "gene", "label"])
    log(f"load: {len(df):,} variants, {df.gene.nunique()} genes, pathogenic (label=1) rate {df.label.mean():.4f}")

    results = {sid: run_split(sid, df, args.alpha, out_root) for sid in args.splits}
    for sid in args.splits:
        verify_saved(sid, out_root)

    rows = []
    summary = {"model": "m3_gene_only", "alpha": args.alpha, "seed": 42, "splits": {}}
    for sid in args.splits:
        m = results[sid]["metrics"]
        ref = {name: load_reference(name, SPLIT_CFG[sid][1]) for name in COMPARE_MODELS}
        rows.append({"split": sid, "accuracy": m["accuracy"], "weighted F1": m["weighted_f1"], "ROC-AUC": m["roc_auc"],
                     "pathogenic recall": m["pathogenic_recall"],
                     "DNABERT-2 AUROC": ref["dnabert2"]["roc_auc"] if ref["dnabert2"] else float("nan"),
                     "NT AUROC": ref["nt"]["roc_auc"] if ref["nt"] else float("nan")})
        summary["splits"][sid] = {**results[sid], "reference_auroc": {k: (v["roc_auc"] if v else None) for k, v in ref.items()}}
    out_root.mkdir(parents=True, exist_ok=True)
    (out_root / "summary.json").write_text(json.dumps(summary, indent=2))
    log(f"save: summary -> {(out_root / 'summary.json').relative_to(REPO_ROOT) if out_root.is_relative_to(REPO_ROOT) else out_root / 'summary.json'}")
    print("\nM3 gene-only baseline (DNABERT-2 / NT AUROC read from existing final_result JSONs):", flush=True)
    print(pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.4f}"), flush=True)


if __name__ == "__main__":
    main()
