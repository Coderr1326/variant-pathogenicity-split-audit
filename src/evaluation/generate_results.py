#!/usr/bin/env python
"""Generate paper-ready result figures and tables from metric histories.

All reproduction values are read from results/metrics/*.json.  The only values
defined in this file are the paper's body-text reference values, which are not
available in machine-readable result files.
"""
from __future__ import annotations

import json
import math
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

ROOT = Path(__file__).resolve().parents[2]
METRICS_DIR = ROOT / "results" / "metrics"
PREDICTIONS_DIR = ROOT / "results" / "predictions"
FIG_DIR = ROOT / "results" / "figures"
TABLE_DIR = ROOT / "results" / "tables"
FIG_DIR.mkdir(parents=True, exist_ok=True)
TABLE_DIR.mkdir(parents=True, exist_ok=True)
DATASET_TOTAL = 168_927
VALIDATION_TOTAL = 33_786

# Paper body-text reference values, sourced from results/reports/paper_specification.md
# under "Results reported by the paper". These are not reproduction metrics.
PAPER_REFERENCE = {
    "cnn": {"weighted_f1": 0.702, "roc_auc": 0.823},
    "dnabert2": {"accuracy": 0.89, "weighted_f1": 0.89, "roc_auc": 0.95},
    "nt": {"accuracy": 0.93, "weighted_f1": 0.93, "roc_auc": 0.97},
}

MODEL_ORDER = ["cnn", "bilstm", "cnn_bilstm", "ensemble", "dnabert2", "nt"]
MODEL_LABELS = {
    "cnn": "CNN",
    "bilstm": "BiLSTM",
    "cnn_bilstm": "CNN+BiLSTM",
    "ensemble": "Ensemble",
    "dnabert2": "DNABERT-2",
    "nt": "NT",
}
LOCAL_MODELS = MODEL_ORDER[:5]
METRICS = ["accuracy", "weighted_f1", "roc_auc"]
METRIC_LABELS = {"accuracy": "Accuracy", "weighted_f1": "Weighted F1", "roc_auc": "ROC-AUC"}


def checkpoint_best_epoch(model: str, final_score: float | None = None, split_suffix: str = "") -> tuple[int | None, str | None]:
    """Resolve the epoch whose state produced the saved best result.

    Checkpoints store cumulative ``best_score`` and ``best_state``. The final
    checkpoint can therefore be later than the epoch that last improved the
    selected score. This rule is model-agnostic and takes precedence over a
    simple last-epoch assumption whenever checkpoint provenance exists.

    ``split_suffix`` (e.g. "_splitB") looks under the split-tagged checkpoint dir
    that scripts/train_transformer_local.py and the Colab notebook write to when
    given --split-manifest; empty (default) is SPLIT-A's untagged directory.
    """
    import torch

    checkpoint_dir = ROOT / "results" / "checkpoints" / (model + split_suffix.lower())
    records = []
    for path in sorted(checkpoint_dir.glob("epoch_*.pt")):
        try:
            checkpoint = torch.load(path, map_location="cpu")
            epoch = int(checkpoint.get("epoch", int(path.stem.split("_")[-1])))
            score = checkpoint.get("best_score")
            if score is not None:
                records.append((epoch, float(score)))
        except (OSError, ValueError, TypeError, RuntimeError, EOFError):
            continue
    if not records:
        return None, None
    records.sort()
    previous = -float("inf")
    last_improvement = None
    for epoch, score in records:
        if score > previous + 1e-12:
            last_improvement = epoch
            previous = score
    if final_score is not None and abs(previous - final_score) > 1e-8:
        return None, f"checkpoint best_score {previous:.10f} does not match final score {final_score:.10f}"
    return last_improvement, "checkpoint cumulative best_score improvement"


def select_best_row(model: str, history: list[dict], final_score: float | None = None, split_suffix: str = "") -> tuple[dict, int, str]:
    """Select metrics using checkpoint provenance, then history as fallback."""
    history_best = max(history, key=lambda row: row["weighted_f1"]) if history else None
    checkpoint_epoch, checkpoint_source = checkpoint_best_epoch(model, final_score, split_suffix)
    if checkpoint_epoch is not None:
        for row in history:
            if int(row.get("epoch")) == checkpoint_epoch:
                return row, checkpoint_epoch, checkpoint_source or "checkpoint provenance"
    if history_best is not None:
        return history_best, int(history_best["epoch"]), "history maximum weighted_f1"
    raise ValueError(f"No metric history or resolvable checkpoint epoch for {model}")


def read_results(split_suffix: str = "") -> dict[str, dict]:
    """Load best-epoch results per model. ``split_suffix`` (e.g. "_splitB") reads
    SPLIT-B/C's tagged files instead of SPLIT-A's untagged ones; a model/split
    combination that hasn't been trained yet (e.g. transformers on B/C) is simply
    absent from the returned dict, not an error -- callers should treat a missing
    key as "not trained yet" rather than fail.
    """
    results = {}
    # Classical models retain their full per-epoch histories.
    for model in ["cnn", "bilstm", "cnn_bilstm", "ensemble"]:
        path = METRICS_DIR / f"{model}_100bp{split_suffix}.json"
        try:
            data = json.loads(path.read_text())
            model = data["model"]
            history = data["history"]
            if not history:
                continue
            best, best_epoch, best_source = select_best_row(model, history, split_suffix=split_suffix)
            results[model] = {"path": path, "data": data, "best": best,
                              "best_epoch": best_epoch, "best_epoch_source": best_source,
                              "curve_history": history, "confusion_source": "metrics history"}
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue

    # Transformer models' final output uses a flat metrics object. The older
    # {model}_100bp.json is a one-epoch sanity run and is deliberately ignored.
    for model in ["dnabert2", "nt"]:
        final_path = METRICS_DIR / f"{model}_100bp{split_suffix}_final.json"
        pred_path = PREDICTIONS_DIR / f"{model}_100bp{split_suffix}.parquet"
        if not (final_path.exists() and pred_path.exists()):
            continue
        try:
            final = json.loads(final_path.read_text())
            metrics = dict(final["metrics"])
            predictions = pd.read_parquet(pred_path)
            required = {"label", "prediction"}
            if not required.issubset(predictions.columns):
                raise ValueError(f"{model} predictions missing columns: {required - set(predictions.columns)}")
            labels = predictions["label"].astype(int)
            predicted = predictions["prediction"].astype(int)
            tn = int(((labels == 0) & (predicted == 0)).sum())
            fp = int(((labels == 0) & (predicted == 1)).sum())
            fn = int(((labels == 1) & (predicted == 0)).sum())
            tp = int(((labels == 1) & (predicted == 1)).sum())
            cm = {"tn": tn, "fp": fp, "fn": fn, "tp": tp}
            if sum(cm.values()) != len(predictions):
                raise ValueError(f"{model} prediction count {len(predictions)} does not match confusion-matrix total")
            for key, value in cm.items():
                if int(metrics[key]) != value:
                    raise ValueError(f"{model} {key} mismatch: final JSON={metrics[key]}, predictions={value}")
            metrics.update(cm)
            best_epoch, best_source = checkpoint_best_epoch(model, metrics["weighted_f1"], split_suffix)
            if best_epoch is None:
                raise ValueError(f"Could not resolve {model} best epoch from checkpoints")
            metrics["epoch"] = best_epoch
            results[model] = {
                "path": final_path, "data": {"model": model, "window": 100, "metrics": metrics},
                "best": metrics, "best_epoch": best_epoch, "best_epoch_source": best_source,
                "curve_history": None,
                "confusion_source": str(pred_path.relative_to(ROOT)),
            }
        except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError):
            continue
    return results


def fmt(value) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "N/A"
    return f"{value:.6f}"


def write_table(df: pd.DataFrame, stem: str, title: str, footnote: str | None = None,
                display_transform=None) -> list[Path]:
    """``footnote`` / ``display_transform`` affect only the .md and .png renderings (the
    transform lets callers add a footnote marker to a cell); the .csv stays unmarked."""
    csv_path = TABLE_DIR / f"{stem}.csv"
    md_path = TABLE_DIR / f"{stem}.md"
    png_path = TABLE_DIR / f"{stem}.png"
    df.to_csv(csv_path, index=False)
    shown_df = display_transform(df) if display_transform else df
    try:
        md = shown_df.to_markdown(index=False)
    except (ImportError, ModuleNotFoundError):
        md = shown_df.to_csv(index=False)
    md_path.write_text(f"# {title}\n\n{md}\n" + (f"\n{footnote}\n" if footnote else ""))

    # Render a readable paper/PPT-ready table image.
    shown = shown_df.astype(object).where(pd.notna(shown_df), "N/A")
    fig_h = max(2.4, 0.38 * (len(shown) + 2))
    fig, ax = plt.subplots(figsize=(max(12, len(df.columns) * 1.55), fig_h))
    ax.axis("off")
    table = ax.table(cellText=shown.values, colLabels=shown.columns, loc="center", cellLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(8.5)
    table.scale(1, 1.5)
    for (row, col), cell in table.get_celld().items():
        if row == 0:
            cell.set_text_props(weight="bold", color="white")
            cell.set_facecolor("#274c77")
        elif row % 2 == 0:
            cell.set_facecolor("#edf2f7")
    ax.set_title(title, pad=18, fontsize=13, weight="bold")
    fig.tight_layout()
    if footnote:
        fig.text(0.5, 0.0, textwrap.fill(footnote, 150), ha="center", va="top", fontsize=8, style="italic", color="dimgray")
    fig.savefig(png_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return [csv_path, md_path, png_path]


def make_training_curves(results: dict[str, dict]) -> Path:
    models = [m for m in MODEL_ORDER if m in results and results[m].get("curve_history")]
    ncols = 3 if len(models) > 4 else 2
    nrows = math.ceil(len(models) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(5.8 * ncols, 4.2 * nrows), squeeze=False)
    for ax, model in zip(axes.flat, models):
        history = results[model]["curve_history"]
        epochs = [row["epoch"] for row in history]
        ax.plot(epochs, [row["weighted_f1"] for row in history], marker="o", label="Weighted F1")
        ax.plot(epochs, [row["roc_auc"] for row in history], marker="s", label="ROC-AUC")
        best = results[model]["best"]["epoch"]
        ax.axvline(best, color="gray", linestyle="--", alpha=0.6, label=f"Best F1 epoch ({best})")
        ax.set_title(MODEL_LABELS.get(model, model))
        ax.set_xlabel("Epoch")
        ax.set_ylabel("Validation score")
        ax.set_ylim(0, 1.02)
        ax.grid(alpha=0.25)
        ax.legend(fontsize=8)
    for ax in axes.flat[len(models):]:
        ax.remove()
    fig.suptitle("Validation Training Curves by Model at 100 bp", fontsize=16, weight="bold")
    fig.text(0.5, 0.01, "DNABERT-2 omitted: final output has no per-epoch metric history; timing log contains duration only.",
             ha="center", fontsize=9, color="dimgray")
    fig.tight_layout(rect=(0, 0.05, 1, 0.96))
    path = FIG_DIR / "local_training_curves.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_confusion_matrices(results: dict[str, dict], output_suffix: str = "", split_label: str = "",
                            models: list[str] | None = None, footnotes: dict[str, str] | None = None) -> list[Path]:
    """``models`` restricts which models get a figure (default: all present);
    ``footnotes`` maps model key -> caption text drawn under that model's figure."""
    paths = []
    for model in MODEL_ORDER:
        if model not in results or (models is not None and model not in models):
            continue
        row = results[model]["best"]
        cm = [[row["tn"], row["fp"]], [row["fn"], row["tp"]]]
        fig, ax = plt.subplots(figsize=(6.6, 5.0))
        sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", cbar=True,
                    xticklabels=["Benign", "Pathogenic"],
                    yticklabels=["Benign", "Pathogenic"], ax=ax)
        ax.set_xlabel("Predicted class")
        ax.set_ylabel("Actual class")
        title = f"Confusion Matrix — {MODEL_LABELS.get(model, model)}"
        if split_label:
            title += f" ({split_label})"
        title += f"\nBest Epoch {row['epoch']}, N={sum(map(sum, cm)):,}"
        ax.set_title(title, fontsize=11)
        fig.tight_layout()
        note = (footnotes or {}).get(model)
        if note:
            fig.text(0.5, 0.0, textwrap.fill(note, 70), ha="center", va="top", fontsize=7.5, style="italic", color="dimgray")
        path = FIG_DIR / f"confusion_matrix_{model}{output_suffix}.png"
        fig.savefig(path, dpi=220, bbox_inches="tight" if note else None)
        plt.close(fig)
        paths.append(path)
    return paths


def local_summary(results: dict[str, dict], models: list[str] | None = None) -> pd.DataFrame:
    rows = []
    for model in (models if models is not None else LOCAL_MODELS):
        if model not in results:
            continue
        row = results[model]["best"]
        rows.append({
            "Model": MODEL_LABELS[model], "Best Epoch": row["epoch"],
            "Accuracy": row["accuracy"], "Weighted Precision": row["weighted_precision"],
            "Weighted Recall": row["weighted_recall"], "Weighted F1": row["weighted_f1"],
            "ROC-AUC": row["roc_auc"],
            "Benign P/R/F1": f"{row['benign_precision']:.6f}/{row['benign_recall']:.6f}/{row['benign_f1']:.6f}",
            "Pathogenic P/R/F1": f"{row['pathogenic_precision']:.6f}/{row['pathogenic_recall']:.6f}/{row['pathogenic_f1']:.6f}",
        })
    return pd.DataFrame(rows)


def make_local_comparison(results: dict[str, dict], output_suffix: str = "", split_label: str = "", n_total: int | None = None) -> Path:
    models = [m for m in LOCAL_MODELS if m in results]
    x = np.arange(len(models))
    width = 0.24
    fig, ax = plt.subplots(figsize=(11, 6))
    for i, metric in enumerate(METRICS):
        vals = [results[m]["best"][metric] for m in models]
        ax.bar(x + (i - 1) * width, vals, width, label=METRIC_LABELS[metric])
    ax.set_xticks(x, [MODEL_LABELS[m] for m in models])
    ax.set_ylabel("Validation score")
    ax.set_xlabel("Model")
    ax.set_ylim(0, 1.05)
    title = f"Local Model Comparison at 100 bp"
    if split_label:
        title += f" — {split_label}"
    title += f" (N={n_total if n_total is not None else DATASET_TOTAL:,})"
    ax.set_title(title)
    ax.legend()
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    path = FIG_DIR / f"local_model_comparison{output_suffix}.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def paper_comparison(results: dict[str, dict]) -> pd.DataFrame:
    rows = []
    for model in MODEL_ORDER:
        best = results.get(model, {}).get("best")
        for metric in METRICS:
            paper = PAPER_REFERENCE.get(model, {}).get(metric)
            reproduction = best.get(metric) if best else None
            rows.append({
                "Model": MODEL_LABELS[model], "Metric": METRIC_LABELS[metric],
                "Paper value": paper, "Reproduction value": reproduction,
                "Absolute difference": None if paper is None or reproduction is None else abs(reproduction - paper),
                "Relative difference": None if paper in (None, 0) or reproduction is None else abs(reproduction - paper) / abs(paper),
            })
    return pd.DataFrame(rows)


def make_paper_comparison_chart(df: pd.DataFrame) -> Path:
    fig, axes = plt.subplots(1, 3, figsize=(19, 6), sharey=True)
    labels = [MODEL_LABELS[m] for m in MODEL_ORDER]
    x = np.arange(len(labels))
    for ax, metric_label in zip(axes, [METRIC_LABELS[m] for m in METRICS]):
        sub = df[df["Metric"] == metric_label].set_index("Model").reindex(labels)
        width = 0.36
        paper = pd.to_numeric(sub["Paper value"], errors="coerce").to_numpy(dtype=float)
        repro = pd.to_numeric(sub["Reproduction value"], errors="coerce").to_numpy(dtype=float)
        ax.bar(x - width / 2, paper, width, label="Paper", color="#6baed6")
        ax.bar(x + width / 2, repro, width, label="Reproduction", color="#2171b5")
        ax.set_title(metric_label)
        ax.set_xticks(x, labels, rotation=35, ha="right")
        ax.set_xlabel("Model")
        ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.25)
        ax.legend()
    axes[0].set_ylabel("Score")
    fig.suptitle(f"Paper Body Values vs. Reproduction at 100 bp (N={DATASET_TOTAL:,})", fontsize=16, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    path = FIG_DIR / "paper_vs_reproduction_comparison.png"
    fig.savefig(path, dpi=220)
    plt.close(fig)
    return path


def make_discrepancy_summary(df: pd.DataFrame) -> Path:
    lines = ["# Automatically Generated Discrepancy Summary", ""]
    for model in [MODEL_LABELS[m] for m in MODEL_ORDER]:
        sub = df[df["Model"] == model]
        parts = []
        for _, row in sub.iterrows():
            if pd.isna(row["Paper value"]) or pd.isna(row["Reproduction value"]):
                continue
            diff = row["Reproduction value"] - row["Paper value"]
            direction = "above" if diff > 0 else "below" if diff < 0 else "equal to"
            parts.append(f"{row['Metric']} reproduction {row['Reproduction value']:.3f} vs. paper {row['Paper value']:.3f}, gap {diff:+.3f} ({abs(diff):.3f} {direction})")
        if parts:
            lines.append(f"- **{model}:** " + "; ".join(parts) + ".")
        else:
            lines.append(f"- **{model}:** no paper reference values are reported for this model.")
    path = TABLE_DIR / "discrepancy_summary.md"
    path.write_text("\n".join(lines) + "\n")
    return path


def main() -> None:
    results = read_results()
    created = []
    created.append(make_training_curves(results))
    created.extend(make_confusion_matrices(results))
    created.append(make_local_comparison(results))
    created.extend(write_table(local_summary(results), "local_results_summary", "Local Model Results at 100 bp"))
    comparison = paper_comparison(results)
    created.append(make_paper_comparison_chart(comparison))
    created.extend(write_table(comparison, "paper_vs_reproduction", "Paper Body Values vs. Reproduction"))
    created.append(make_discrepancy_summary(comparison))

    windows = sorted({int(p.stem.rsplit("_", 1)[1].replace("bp", "")) for p in METRICS_DIR.glob("*_*.json") if "_100bp" not in p.name and p.name.endswith(".json")})
    log = TABLE_DIR / "plot_generation.log"
    log.write_text(
        "Generated from metric histories/final outputs using best validation weighted-F1 selection where history exists.\n"
        f"Canonical dataset size: {DATASET_TOTAL:,} variants; DNABERT-2 validation predictions: {VALIDATION_TOTAL:,} (20% split).\n"
        + ("Window-size comparison skipped: no 30/50/200bp metric files found.\n" if not windows else f"Window-size results detected: {windows}; no window table generator was requested in this pass.\n")
        + ("DNABERT-2 included from results/metrics/dnabert2_100bp_final.json and confusion counts recomputed from results/predictions/dnabert2_100bp.parquet; best epoch resolved from checkpoint cumulative best_score provenance.\n" if "dnabert2" in results else "DNABERT-2 skipped: final metrics/prediction files unavailable or inconsistent.\n")
        + ("Nucleotide Transformer included from results/metrics/nt_100bp_final.json and confusion counts recomputed from results/predictions/nt_100bp.parquet; best epoch resolved from checkpoint cumulative best_score provenance.\n" if "nt" in results else "Nucleotide Transformer skipped: final metrics/prediction files unavailable or inconsistent, so its rows/bars are N/A.\n")
    )
    created.append(log)

    print("Best epochs selected programmatically:")
    for model in MODEL_ORDER:
        if model in results:
            row = results[model]["best"]
            print(f"  {MODEL_LABELS[model]}: epoch {row['epoch']} (weighted F1 {row['weighted_f1']:.6f}; source: {results[model].get('best_epoch_source', 'unknown')})")
    print("\nCreated files:")
    for path in created:
        print(f"  {path.relative_to(ROOT)}")
    print("\nSkipped: window-size comparison (no 30/50/200bp results).")
    if "nt" not in results:
        print("Pending: Nucleotide Transformer (no trained result included).")


if __name__ == "__main__":
    main()
