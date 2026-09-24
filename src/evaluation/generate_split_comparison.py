#!/usr/bin/env python
"""Generate SPLIT-B/C figures+tables and a cross-split A/B/C comparison (classical models + DNABERT-2 + NT).

Reuses generate_results.py's read_results / make_confusion_matrices / make_local_comparison /
local_summary / write_table as-is (now split-suffix-parameterized there) -- this file adds no new
metric-reading or best-epoch-selection logic, only a cross-split comparison layer that doesn't
exist for SPLIT-A alone. All six models now have SPLIT-A/B/C results and are picked up automatically. Tables iterate
MODEL_ORDER (all six), so a model/split missing its results would show "pending" rather than
being silently dropped.

By default, per-split figures are (re)generated only for the transformer models, so the
already-verified classical confusion matrices / local_model_comparison charts are left untouched;
pass --regenerate-classical to rebuild those too.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from src.evaluation.generate_results import (
    ROOT, FIG_DIR, TABLE_DIR, MODEL_ORDER, MODEL_LABELS,
    read_results, make_confusion_matrices, make_local_comparison, local_summary, write_table,
)

SPLITS = [
    ("A", "", "SPLIT-A (random 80/20)"),
    ("B", "_splitB", "SPLIT-B (gene-disjoint)"),
    ("C", "_splitC", "SPLIT-C (temporal)"),
]
CLASSICAL_MODELS = ["cnn", "bilstm", "cnn_bilstm", "ensemble"]
TRANSFORMER_MODELS = ["dnabert2", "nt"]
CHART_MODELS = CLASSICAL_MODELS + ["dnabert2", "nt"]
# Research-transparency note: required wherever DNABERT-2's SPLIT-B number appears.
WARMUP_MARK = "\u2020"
WARMUP_NOTE = (f"{WARMUP_MARK} DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced "
               "representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule.")
CROSS_METRICS = ["accuracy", "weighted_f1", "roc_auc", "pathogenic_recall"]
CROSS_METRIC_LABELS = {
    "accuracy": "Accuracy", "weighted_f1": "Weighted F1",
    "roc_auc": "ROC-AUC", "pathogenic_recall": "Pathogenic Recall",
}
# Short forms for the leakage-delta table only -- 8 metric columns plus Model needs
# tighter headers than the per-split tables to avoid overlapping in write_table()'s
# fixed-width rendering.
SHORT_METRIC_LABELS = {
    "accuracy": "Acc", "weighted_f1": "F1",
    "roc_auc": "AUC", "pathogenic_recall": "PathRecall",
}
SPLIT_COLORS = {"A": "#2171b5", "B": "#fd8d3c", "C": "#a50f15"}


def _shorten_prf(value: str, digits: int = 3) -> str:
    """"0.596123/0.782456/0.676789" -> "0.596/0.782/0.677": local_results_summary's
    "Benign P/R/F1" and "Pathogenic P/R/F1" columns carry three full-precision floats
    each, which overflow write_table()'s fixed per-column width. Only these two
    columns are affected -- Accuracy/Weighted F1/ROC-AUC etc. stay full precision,
    matching SPLIT-A's local_results_summary.png untouched."""
    return "/".join(f"{float(p):.{digits}f}" for p in value.split("/"))


def local_summary_for_table(results: dict[str, dict]) -> pd.DataFrame:
    df = local_summary(results, models=MODEL_ORDER)
    for col in ("Benign P/R/F1", "Pathogenic P/R/F1"):
        if col in df.columns:
            df[col] = df[col].apply(_shorten_prf)
    return df


def build_results_by_split() -> dict[str, dict]:
    return {sid: read_results(suffix) for sid, suffix, _label in SPLITS}


def cross_split_table(results_by_split: dict[str, dict]) -> pd.DataFrame:
    """All 6 models x all 3 splits x [accuracy, weighted_f1, roc_auc, pathogenic_recall].
    A model/split not yet trained (currently: dnabert2/nt on B/C) reads "pending"."""
    rows = []
    for model in MODEL_ORDER:
        for sid, _suffix, _label in SPLITS:
            best = results_by_split[sid].get(model, {}).get("best")
            row = {"Model": MODEL_LABELS[model], "Split": sid,
                   "Best Epoch": best["epoch"] if best else "pending"}
            for metric in CROSS_METRICS:
                row[CROSS_METRIC_LABELS[metric]] = round(best[metric], 4) if best else "pending"
            rows.append(row)
    return pd.DataFrame(rows)


def leakage_delta_table(results_by_split: dict[str, dict]) -> pd.DataFrame:
    """Per model: (SPLIT-A score - SPLIT-B score) and (SPLIT-A score - SPLIT-C score), the D1
    leakage-gap diagnostic. "pending" wherever either side of the subtraction isn't trained yet."""
    rows = []
    res_a = results_by_split["A"]
    for model in MODEL_ORDER:
        best_a = res_a.get(model, {}).get("best")
        row = {"Model": MODEL_LABELS[model]}
        for sid in ["B", "C"]:
            best_x = results_by_split[sid].get(model, {}).get("best")
            for metric in CROSS_METRICS:
                col = f"D1 A-{sid} {SHORT_METRIC_LABELS[metric]}"
                row[col] = round(best_a[metric] - best_x[metric], 4) if (best_a and best_x) else "pending"
        rows.append(row)
    return pd.DataFrame(rows)


def make_cross_split_chart(results_by_split: dict[str, dict]) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(20, 6), sharey=True)
    x = np.arange(len(CHART_MODELS))
    width = 0.25
    for ax, metric in zip(axes, ["accuracy", "weighted_f1", "roc_auc"]):
        for i, (sid, _suffix, _label) in enumerate(SPLITS):
            vals = []
            for m in CHART_MODELS:
                best = results_by_split[sid].get(m, {}).get("best")
                vals.append(best[metric] if best else np.nan)
            bars = ax.bar(x + (i - 1) * width, vals, width, label=f"SPLIT-{sid}", color=SPLIT_COLORS[sid])
            labels = [f"{v:.3f}" + (WARMUP_MARK if (sid == "B" and m == "dnabert2") else "") for v, m in zip(vals, CHART_MODELS)]
            ax.bar_label(bars, labels=labels, fontsize=7, padding=2)
        ax.set_title(CROSS_METRIC_LABELS[metric])
        ax.set_xticks(x, [MODEL_LABELS[m] for m in CHART_MODELS], rotation=20, ha="right")
        ax.set_ylim(0, 1.0)
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Validation score")
    axes[0].legend()
    fig.suptitle("All Models Across Splits A -> B -> C (Leakage Delta D1)", fontsize=15, weight="bold")
    fig.tight_layout(rect=(0, 0.04, 1, 0.93))
    fig.text(0.5, 0.012, WARMUP_NOTE, ha="center", fontsize=9, style="italic", color="dimgray")
    path = FIG_DIR / "cross_split_comparison_chart.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def _mark_dnabert2_model(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Model"] = df["Model"].where(df["Model"] != MODEL_LABELS["dnabert2"], MODEL_LABELS["dnabert2"] + WARMUP_MARK)
    return df


def _mark_dnabert2_splitb_row(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    hit = (df["Model"] == MODEL_LABELS["dnabert2"]) & (df["Split"] == "B")
    df["Split"] = df["Split"].astype(object)
    df.loc[hit, "Split"] = "B" + WARMUP_MARK
    return df


def write_methodology_note(results_by_split: dict[str, dict]) -> Path:
    b = results_by_split["B"].get("dnabert2", {}).get("best")
    path = TABLE_DIR / "dnabert2_splitB_methodology_note.md"
    lines = [
        "# DNABERT-2 SPLIT-B methodology note", "",
        WARMUP_NOTE, "",
        "- The first SPLIT-B attempt (fixed LR 2e-5, no warmup) collapsed: accuracy 0.5377, ROC-AUC 0.4999, "
        "tp=0 (every validation example predicted benign). It was discarded.",
        "- The reported SPLIT-B result is the rerun with linear LR warmup (`--warmup-ratio 0.1`, warmup over the "
        "first 10% of total steps to 2e-5, then linear decay; `scripts/train_transformer_local.py`).",
        "- NT needed no warmup: NT SPLIT-A/B/C all used the fixed-LR schedule (`--warmup-ratio 0`), so NT rows carry no marker.",
        "- SPLIT-A and SPLIT-C DNABERT-2 runs used the original fixed-LR schedule (no scheduler). "
        "The DNABERT-2 SPLIT-A/C vs SPLIT-B comparison therefore also differs in LR schedule, not only in split.",
    ]
    if b:
        lines.append(f"- Reported SPLIT-B numbers (best epoch {b['epoch']}): accuracy {b['accuracy']:.4f}, weighted F1 {b['weighted_f1']:.4f}, "
                     f"ROC-AUC {b['roc_auc']:.4f}, pathogenic recall {b['pathogenic_recall']:.4f}.")
    lines += [
        "- `results/metrics/dnabert2_100bp_splitB_epoch_times.csv` is append-only and contains BOTH runs: the first 10 rows "
        "(2026-09-18) are the discarded collapsed run, the last 10 rows (2026-09-19) are the reported rerun.", "",
        f"Where marked ({WARMUP_MARK}): confusion_matrix_dnabert2_splitB.png, local_results_summary_splitB, cross_split_comparison, "
        "leakage_delta, cross_split_comparison_chart.png. The .csv files carry no marker (kept machine-clean); this note applies to them.",
    ]
    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--per-split-models", nargs="+", default=TRANSFORMER_MODELS, choices=MODEL_ORDER,
                    help="models to (re)draw per-split confusion matrices for (default: both transformers)")
    ap.add_argument("--regenerate-classical", action="store_true",
                    help="also rebuild the classical models' per-split confusion matrices and local_model_comparison charts "
                         "(default: leave those already-verified files untouched)")
    args = ap.parse_args()
    results_by_split = build_results_by_split()
    created: list[Path] = []

    # Per-split figures/table for B and C. SPLIT-A's own outputs are
    # generate_results.py's job and are not touched or regenerated here.
    for sid, suffix, label in SPLITS:
        if sid == "A":
            continue
        results = results_by_split[sid]
        cm_models = None if args.regenerate_classical else args.per_split_models
        cm_notes = {"dnabert2": WARMUP_NOTE} if sid == "B" else None
        created.extend(make_confusion_matrices(results, output_suffix=suffix, split_label=label,
                                               models=cm_models, footnotes=cm_notes))
        if args.regenerate_classical:
            created.append(make_local_comparison(results, output_suffix=suffix, split_label=label))
        created.extend(write_table(local_summary_for_table(results), f"local_results_summary{suffix}",
                                    f"Local Model Results at 100 bp -- {label}",
                                    footnote=WARMUP_NOTE if sid == "B" else None,
                                    display_transform=_mark_dnabert2_model if sid == "B" else None))

    # Cross-split comparison (new; doesn't exist for SPLIT-A alone).
    cross_df = cross_split_table(results_by_split)
    created.extend(write_table(cross_df, "cross_split_comparison", "Model Results Across Splits A/B/C",
                                footnote=WARMUP_NOTE, display_transform=_mark_dnabert2_splitb_row))

    delta_df = leakage_delta_table(results_by_split)
    created.extend(write_table(delta_df, "leakage_delta", "Leakage Delta (D1): SPLIT-A minus SPLIT-B/C",
                                footnote=WARMUP_NOTE, display_transform=_mark_dnabert2_model))
    created.append(write_methodology_note(results_by_split))

    created.append(make_cross_split_chart(results_by_split))

    print("Trained models found per split:")
    for sid, _suffix, label in SPLITS:
        found = [MODEL_LABELS[m] for m in MODEL_ORDER if m in results_by_split[sid]]
        missing = [MODEL_LABELS[m] for m in MODEL_ORDER if m not in results_by_split[sid]]
        print(f"  {label}: {found}" + (f"  (pending: {missing})" if missing else ""))

    print("\nCreated files:")
    for path in created:
        print(f"  {path.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
