#!/usr/bin/env python
"""Fixes to the M3 outputs plus the all-model split-shift delta table. Add-only companion to make_m3_figures.py.

make_m3_figures.py is deliberately left untouched (re-running it would restore the old formats); this script
imports its helpers and regenerates only:
  * local_results_summary_m3.{csv,md,png}   -- metric columns rounded to 4 decimals in md/png (csv stays numeric)
  * confusion_matrix_m3{,_splitB,_splitC}.png -- vmin=0 and one common vmax across the three splits
and creates new: split_shift_delta_all_models.{csv,md,png}. CPU only; every value is read from result JSONs.
"""
from __future__ import annotations

import os

os.environ["CUDA_VISIBLE_DEVICES"] = ""

import sys
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import make_m3_figures as base  # noqa: E402
from src.evaluation.generate_results import FIG_DIR, MODEL_LABELS, read_results, write_table  # noqa: E402

DELTA_MODELS = ["m3", "cnn", "bilstm", "cnn_bilstm", "ensemble", "dnabert2", "nt"]
DELTA_METRICS = [("accuracy", "Acc"), ("weighted_f1", "F1"), ("roc_auc", "AUC"), ("pathogenic_recall", "PathRecall")]
M3_MARK = "‡"


def results_table_rounded(m3: dict) -> list[Path]:
    rows = []
    for sid in base.SUFFIX:
        m = m3[sid]["metrics"]
        prf = lambda p: base._shorten_prf(f"{m[p + '_precision']:.6f}/{m[p + '_recall']:.6f}/{m[p + '_f1']:.6f}")
        rows.append({"Model": base.M3_LABEL, "Split": sid, "Accuracy": m["accuracy"], "Weighted Precision": m["weighted_precision"],
                     "Weighted Recall": m["weighted_recall"], "Weighted F1": m["weighted_f1"], "ROC-AUC": m["roc_auc"],
                     "Benign P/R/F1": prf("benign"), "Pathogenic P/R/F1": prf("pathogenic")})
    note = "SPLIT-B: " + base.chance_floor_text(m3["B"]["diag"]) + "; AUROC = 0.5 exactly and tp = 0 by construction (expected, not a bug)."
    return write_table(pd.DataFrame(rows), "local_results_summary_m3", "M3 Gene-Only Baseline Results at 100 bp (no sequence input)",
                       footnote=note, display_transform=base.fixed4)


def confusion_matrices_common_scale(m3: dict) -> list[Path]:
    vmax = max(m3[s]["metrics"][k] for s in base.SUFFIX for k in ("tn", "fp", "fn", "tp"))  # from the JSONs
    paths = []
    for sid, suffix, label in base.SPLITS:
        m, diag = m3[sid]["metrics"], m3[sid]["diag"]
        cm = [[m["tn"], m["fp"]], [m["fn"], m["tp"]]]
        fig, ax = plt.subplots(figsize=(6.6, 5.0))
        sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=True, vmin=0, vmax=vmax,
                    cbar_kws={"label": "Variants (common scale across SPLIT-A/B/C)"},
                    xticklabels=["Benign", "Pathogenic"], yticklabels=["Benign", "Pathogenic"], ax=ax)
        ax.set_xlabel("Predicted class")
        ax.set_ylabel("Actual class")
        ax.set_title(f"Confusion Matrix — M3 Gene-Only ({label})\nNo sequence input, N={sum(map(sum, cm)):,}", fontsize=11)
        fig.tight_layout()
        if sid == "B":
            note = base.chance_floor_text(diag) + f". Expected behaviour, not a bug: AUROC = {m['roc_auc']:.4f}, tp = {m['tp']}."
            fig.text(0.5, 0.0, textwrap.fill(note, 70), ha="center", va="top", fontsize=7.5, style="italic", color="dimgray")
        path = FIG_DIR / f"confusion_matrix_m3{suffix}.png"
        fig.savefig(path, dpi=220, bbox_inches="tight" if sid == "B" else None)
        plt.close(fig)
        paths.append(path)
    print(f"confusion-matrix colour scale: vmin=0, vmax={vmax} (max cell over the three splits' JSONs)")
    return paths


def collect_split_metrics(m3: dict) -> tuple[dict, list[str]]:
    """vals[model][split][metric]; anything not found is listed in `missing`, never filled in."""
    vals = {m: {} for m in DELTA_MODELS}
    missing = []
    for sid in base.SUFFIX:
        vals["m3"][sid] = {k: m3[sid]["metrics"].get(k) for k, _ in DELTA_METRICS}
    for sid, suffix, _ in base.SPLITS:
        res = read_results(suffix)
        for model in DELTA_MODELS[1:]:
            best = res.get(model, {}).get("best")
            if best is None:
                missing.append(f"{model} SPLIT-{sid}: no result found")
                continue
            vals[model][sid] = {k: best.get(k) for k, _ in DELTA_METRICS}
    for model, per in vals.items():
        for sid, d in per.items():
            missing += [f"{model} SPLIT-{sid}: metric {k} missing" for k, v in d.items() if v is None]
    return vals, missing


def split_shift_table(m3: dict) -> list[Path]:
    vals, missing = collect_split_metrics(m3)
    rows = []
    for model in DELTA_MODELS:
        row = {"Model": base.M3_LABEL if model == "m3" else MODEL_LABELS[model]}
        for x in ("B", "C"):
            for key, short in DELTA_METRICS:
                a, b = vals[model].get("A", {}).get(key), vals[model].get(x, {}).get(key)
                row[f"A->{x} {short}"] = round(a - b, 4) if (a is not None and b is not None) else float("nan")
        rows.append(row)
    df = pd.DataFrame(rows)

    seeds = base.dnabert2_splitb_seed_aurocs()
    s = pd.Series(seeds)
    a_auc = vals["dnabert2"]["A"]["roc_auc"]
    seed_txt = "/".join(str(k) for k in sorted(seeds, key=lambda k: (k != 42, k)))
    note = (f"{base.WARMUP_NOTE} DNABERT-2 SPLIT-B AUROC over seeds {seed_txt} (warmup 0.1): mean {s.mean():.4f}, min-max {s.min():.4f}-{s.max():.4f} "
            f"(A->B AUROC delta over these seeds: mean {a_auc - s.mean():.4f}, min-max {a_auc - s.max():.4f}-{a_auc - s.min():.4f}); "
            f"the table row uses seed 42. {M3_MARK} M3 SPLIT-B is a chance floor by construction (constant global prior; "
            f"{m3['B']['diag']['n_test_genes_seen_in_train']} of {m3['B']['diag']['n_test_genes']} test genes seen in train), so its A->B deltas "
            f"measure loss of gene signal, not model degradation.")

    def mark(d: pd.DataFrame) -> pd.DataFrame:
        d = base.fixed4(d)
        d["Model"] = d["Model"].replace({MODEL_LABELS["dnabert2"]: MODEL_LABELS["dnabert2"] + base.WARMUP_MARK,
                                         base.M3_LABEL: base.M3_LABEL + M3_MARK})
        return d.replace("nan", "N/A")

    if missing:
        print("VALUES NOT FOUND (left as N/A):", *missing, sep="\n  ")
    else:
        print("all 7 models x 3 splits x 4 metrics found in the result JSONs")
    return write_table(df, "split_shift_delta_all_models", "Split-Shift Delta (SPLIT-A minus SPLIT-B/C)", footnote=note, display_transform=mark)


def main() -> None:
    m3 = base.load_m3()
    created = []
    created += results_table_rounded(m3)
    created += confusion_matrices_common_scale(m3)
    created += split_shift_table(m3)
    print("\nWritten:")
    for p in created:
        print(f"  {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
