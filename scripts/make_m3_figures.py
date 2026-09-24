#!/usr/bin/env python
"""M3 (gene-identity-only baseline) figures and tables, styled like the existing 6-model pipeline.

Adds NEW files only (results/figures/*_m3*, results/tables/*m3*); never rewrites an existing figure or
table. Reuses write_table / SPLIT_COLORS / the dagger note / the P/R-F1 3-decimal shortening from
src.evaluation. Every number is read from result JSONs (M3 final_result.json, the other models'
final/history JSONs, the DNABERT-2 SPLIT-B seed JSONs), nothing is hardcoded. CPU only.
"""
from __future__ import annotations

import os

os.environ["CUDA_VISIBLE_DEVICES"] = ""

import json
import sys
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from src.evaluation.generate_results import FIG_DIR, METRICS_DIR, MODEL_LABELS, TABLE_DIR, read_results, write_table  # noqa: E402
from src.evaluation.generate_split_comparison import (  # noqa: E402
    CROSS_METRICS, SHORT_METRIC_LABELS, SPLIT_COLORS, SPLITS, WARMUP_MARK, WARMUP_NOTE, _shorten_prf,
)

M3_DIR = ROOT / "results" / "m3_gene_only"
M3_LABEL = "M3 (gene-only)"
CHART_MODELS = ["m3", "cnn", "bilstm", "cnn_bilstm", "ensemble", "dnabert2", "nt"]
SUFFIX = {sid: suffix for sid, suffix, _ in SPLITS}
SPLIT_LABEL = {sid: label for sid, _, label in SPLITS}
DNABERT2_MAIN_WARMUP = 0.1


def load_m3() -> dict[str, dict]:
    out = {}
    for sid in SUFFIX:
        d = M3_DIR / f"split_{sid}"
        out[sid] = {"metrics": json.loads((d / "final_result.json").read_text())["metrics"],
                    "diag": json.loads((d / "diagnostics.json").read_text())}
    return out


def chance_floor_text(diag: dict) -> str:
    return (f"chance floor: constant global prior ({diag['n_test_genes_seen_in_train']} of {diag['n_test_genes']} "
            f"test genes seen in train)")


def m3_confusion_matrices(m3: dict) -> list[Path]:
    paths = []
    for sid, suffix, label in SPLITS:
        m, diag = m3[sid]["metrics"], m3[sid]["diag"]
        cm = [[m["tn"], m["fp"]], [m["fn"], m["tp"]]]
        fig, ax = plt.subplots(figsize=(6.6, 5.0))
        sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=True,
                    xticklabels=["Benign", "Pathogenic"], yticklabels=["Benign", "Pathogenic"], ax=ax)
        ax.set_xlabel("Predicted class")
        ax.set_ylabel("Actual class")
        ax.set_title(f"Confusion Matrix — M3 Gene-Only ({label})\nNo sequence input, N={sum(map(sum, cm)):,}", fontsize=11)
        fig.tight_layout()
        if sid == "B":
            note = (chance_floor_text(diag) + f". Expected behaviour, not a bug: AUROC = {m['roc_auc']:.4f}, tp = {m['tp']}.")
            fig.text(0.5, 0.0, textwrap.fill(note, 70), ha="center", va="top", fontsize=7.5, style="italic", color="dimgray")
        path = FIG_DIR / f"confusion_matrix_m3{suffix}.png"
        fig.savefig(path, dpi=220, bbox_inches="tight" if sid == "B" else None)
        plt.close(fig)
        paths.append(path)
    return paths


def m3_results_table(m3: dict) -> list[Path]:
    rows = []
    for sid in SUFFIX:
        m = m3[sid]["metrics"]
        prf = lambda p: _shorten_prf(f"{m[p + '_precision']:.6f}/{m[p + '_recall']:.6f}/{m[p + '_f1']:.6f}")
        rows.append({"Model": M3_LABEL, "Split": sid, "Accuracy": m["accuracy"], "Weighted Precision": m["weighted_precision"],
                     "Weighted Recall": m["weighted_recall"], "Weighted F1": m["weighted_f1"], "ROC-AUC": m["roc_auc"],
                     "Benign P/R/F1": prf("benign"), "Pathogenic P/R/F1": prf("pathogenic")})
    note = "SPLIT-B: " + chance_floor_text(m3["B"]["diag"]) + "; AUROC = 0.5 exactly and tp = 0 by construction (expected, not a bug)."
    return write_table(pd.DataFrame(rows), "local_results_summary_m3", "M3 Gene-Only Baseline Results at 100 bp (no sequence input)", footnote=note)


def dnabert2_splitb_seed_aurocs() -> dict[int, float]:
    seeds = {}
    for p in sorted(METRICS_DIR.glob("dnabert2_100bp_splitB*_final.json")):
        j = json.loads(p.read_text())
        if abs(j.get("warmup_ratio", DNABERT2_MAIN_WARMUP) - DNABERT2_MAIN_WARMUP) > 1e-12:
            continue  # the warmup-0 secondary check is not part of the seed range
        seeds[int(j.get("seed", 42))] = float(j["metrics"]["roc_auc"])
    return seeds


def collect_auroc(m3: dict) -> dict[str, dict[str, float]]:
    """auroc[model][split]; classical/transformer values come from the same read_results() the existing figures use."""
    auroc = {"m3": {sid: float(m3[sid]["metrics"]["roc_auc"]) for sid in SUFFIX}}
    for sid, suffix, _ in SPLITS:
        for model, r in read_results(suffix).items():
            auroc.setdefault(model, {})[sid] = float(r["best"]["roc_auc"])
    return auroc


def m3_comparison_chart(auroc: dict, seeds: dict[int, float], m3: dict) -> Path:
    fig, ax = plt.subplots(figsize=(16, 6.8))
    x = np.arange(len(CHART_MODELS))
    width = 0.26
    lo, hi = (min(seeds.values()), max(seeds.values())) if seeds else (None, None)
    for i, (sid, _suffix, _label) in enumerate(SPLITS):
        vals = np.array([auroc[m].get(sid, np.nan) for m in CHART_MODELS])
        pos = x + (i - 1) * width
        bars = ax.bar(pos, vals, width, color=SPLIT_COLORS[sid])
        for bar, model, v in zip(bars, CHART_MODELS, vals):
            top = v
            if model == "dnabert2" and sid == "B" and seeds:
                err = [[max(0.0, v - lo)], [max(0.0, hi - v)]]
                ax.errorbar(bar.get_x() + width / 2, v, yerr=err, fmt="none", ecolor="black", elinewidth=1.4, capsize=4, zorder=5)
                top = max(v, hi)
            if model == "m3" and sid == "B":
                bar.set_hatch("//")
                bar.set_edgecolor("black")
                bar.set_linewidth(0.8)
            mark = WARMUP_MARK if (model == "dnabert2" and sid == "B") else ""
            ax.text(bar.get_x() + width / 2, top + 0.008, f"{v:.3f}{mark}", ha="center", va="bottom", fontsize=7.5)
    ax.axhline(0.5, color="dimgray", linestyle="--", linewidth=1.1, zorder=1)
    ax.annotate("chance floor\nby construction", xy=(x[0], 0.55), xytext=(x[0], 0.86), ha="center", fontsize=8.5,
                arrowprops=dict(arrowstyle="->", color="black"))
    if seeds:
        d2 = CHART_MODELS.index("dnabert2")
        ax.annotate(f"{len(seeds)}-seed range\n{lo:.3f}–{hi:.3f}", xy=(x[d2], hi + 0.055), xytext=(x[d2], 1.005),
                    ha="center", va="center", fontsize=8, bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="gray", lw=0.8),
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.8), zorder=6)
    ax.set_xticks(x, [M3_LABEL if m == "m3" else MODEL_LABELS[m] for m in CHART_MODELS], rotation=15, ha="right")
    ax.set_ylabel("ROC-AUC (validation)")
    ax.set_ylim(0, 1.08)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.grid(axis="y", alpha=0.25)
    ax.set_title("ROC-AUC Across Splits A -> B -> C: M3 Gene-Only Baseline vs. All Models", fontsize=14, weight="bold")
    handles = [Patch(facecolor=SPLIT_COLORS[s], label=f"SPLIT-{s}") for s in SUFFIX]
    handles += [Line2D([0], [0], color="dimgray", linestyle="--", label="AUROC 0.5 (chance)"),
                Line2D([0], [0], color="black", marker="_", linestyle="none", markersize=10, label="DNABERT-2 SPLIT-B seed range")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=5, fontsize=9, frameon=False)
    seed_txt = ", ".join(str(s) for s in sorted(seeds, key=lambda s: (s != 42, s))) if seeds else "n/a"
    foot = [WARMUP_NOTE,
            f"M3 SPLIT-B: {chance_floor_text(m3['B']['diag'])}; AUROC exactly 0.5 by construction. "
            f"DNABERT-2 SPLIT-B bar = seed 42; error bar = min-max over seeds {seed_txt} (all warmup 0.1)."]
    fig.tight_layout(rect=(0, 0.11, 1, 1))
    fig.text(0.5, 0.035, "\n".join(textwrap.fill(t, 190) for t in foot), ha="center", va="center", fontsize=8.5, style="italic", color="dimgray")
    path = FIG_DIR / "cross_split_comparison_chart_m3.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def fixed4(d: pd.DataFrame) -> pd.DataFrame:
    """Display-only: every float column as fixed 4 decimals (keeps 0.5000 / 0.6740 instead of 0.5 / 0.674)."""
    d = d.copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: f"{v:.4f}")
    return d


def m3_vs_transformers_table(auroc: dict) -> list[Path]:
    rows = []
    for sid in SUFFIX:
        m3v, d2, nt = auroc["m3"][sid], auroc["dnabert2"][sid], auroc["nt"][sid]
        rows.append({"Split": sid, "M3 AUROC": round(m3v, 4), "DNABERT-2 AUROC": round(d2, 4), "NT AUROC": round(nt, 4),
                     "DNABERT-2 minus M3": round(d2 - m3v, 4), "NT minus M3": round(nt - m3v, 4)})
    df = pd.DataFrame(rows)

    def mark(d: pd.DataFrame) -> pd.DataFrame:
        d = fixed4(d)
        d["Split"] = d["Split"].astype(object)
        d.loc[d["Split"] == "B", "Split"] = "B" + WARMUP_MARK
        return d

    note = WARMUP_NOTE + " M3 SPLIT-B AUROC is 0.5 by construction (constant global prior; no test gene seen in train)."
    return write_table(df, "m3_vs_transformers_auroc", "M3 Gene-Only Baseline vs. Transformers (ROC-AUC)", footnote=note, display_transform=mark)


def m3_leakage_delta_table(m3: dict) -> list[Path]:
    a = m3["A"]["metrics"]
    row = {"Model": M3_LABEL}
    for sid in ["B", "C"]:
        x = m3[sid]["metrics"]
        for metric in CROSS_METRICS:
            row[f"D1 A-{sid} {SHORT_METRIC_LABELS[metric]}"] = round(a[metric] - x[metric], 4)
    note = ("M3 SPLIT-B is a chance floor by construction (constant global prior), so its A-B deltas measure the loss of gene-identity "
            "signal, not model degradation.")
    return write_table(pd.DataFrame([row]), "leakage_delta_m3", "Leakage Delta (D1): SPLIT-A minus SPLIT-B/C -- M3 Gene-Only Baseline", footnote=note,
                       display_transform=fixed4)


def main() -> None:
    m3 = load_m3()
    created: list[Path] = []
    created += m3_confusion_matrices(m3)
    created += m3_results_table(m3)
    auroc = collect_auroc(m3)
    seeds = dnabert2_splitb_seed_aurocs()
    print("DNABERT-2 SPLIT-B warmup-0.1 seed AUROCs found:", {s: round(v, 4) for s, v in sorted(seeds.items())})
    created.append(m3_comparison_chart(auroc, seeds, m3))
    pd.DataFrame([{"model": m, "split": s, "roc_auc": v} for m in CHART_MODELS for s, v in auroc.get(m, {}).items()]
                 ).to_csv(TABLE_DIR / "cross_split_auroc_with_m3.csv", index=False)
    pd.DataFrame([{"seed": s, "warmup_ratio": DNABERT2_MAIN_WARMUP, "roc_auc": v} for s, v in sorted(seeds.items())]
                 ).to_csv(TABLE_DIR / "dnabert2_splitB_seed_auroc.csv", index=False)
    created += [TABLE_DIR / "cross_split_auroc_with_m3.csv", TABLE_DIR / "dnabert2_splitB_seed_auroc.csv"]
    created += m3_vs_transformers_table(auroc)
    created += m3_leakage_delta_table(m3)
    print("\nCreated files:")
    for p in created:
        print(f"  {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
