#!/usr/bin/env python
"""D5: per-gene breakdown on SPLIT-B (14 test genes) and pooled vs within-gene vs between-gene AUROC decomposition.

Uses scripts/pred_loader.py (validated join of saved predictions to genes). No model inference. Runs available:
M3, NT, DNABERT-2 (seeds 42/1/2/3, warmup 0.1). CNN/BiLSTM/CNN+BiLSTM/Ensemble have no saved probabilities and are skipped.
The DNABERT-2 warmup-0 run (collapsed) is not used. Outputs go to results/d5_per_gene/ (add-only). CPU only.
Numbers only; no interpretation is written by this script.
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
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import pred_loader as pl  # noqa: E402
import src.evaluation.generate_results as gr  # noqa: E402
from src.evaluation.generate_split_comparison import WARMUP_MARK, WARMUP_NOTE  # noqa: E402

OUT = ROOT / "results" / "d5_per_gene"
OUT.mkdir(parents=True, exist_ok=True)
gr.TABLE_DIR = OUT  # write_table resolves TABLE_DIR at call time: md/png/csv views land in results/d5_per_gene/

MIN_PER_CLASS, N_BOOT = 30, 1000
SEEDS = list(pl.DNABERT2_B_SEEDS)
SKIP_NOTE = "Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities)."
INSUFF_NOTE = f"insufficient = fewer than {MIN_PER_CLASS} pathogenic or benign variants in the gene."
f4 = lambda v: "N/A" if pd.isna(v) else f"{v:.4f}"  # noqa: E731


def auroc(y, s) -> float:
    y = np.asarray(y)
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(s)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def per_gene(d: pd.DataFrame, genes: list[str]) -> pd.DataFrame:
    rows = []
    for gene in genes:
        g = d[d.gene == gene]
        y, s = g.label.to_numpy(), g.score.to_numpy()
        n, n1 = len(y), int(y.sum())
        n0 = n - n1
        valid = n1 >= MIN_PER_CLASS and n0 >= MIN_PER_CLASS
        row = {"gene": gene, "n": n, "n_pathogenic": n1, "n_benign": n0, "valid": valid, "auroc": np.nan, "ci_lo": np.nan, "ci_hi": np.nan,
               "accuracy": float(((s >= 0.5).astype(int) == y).mean()) if n else np.nan, "pred_pos_frac": float((s >= 0.5).mean()) if n else np.nan}
        if valid:
            row["auroc"] = auroc(y, s)
            idx = np.random.default_rng(0).integers(0, n, size=(N_BOOT, n))
            vals = np.array([auroc(y[i], s[i]) for i in idx])
            row["ci_lo"], row["ci_hi"] = (float(v) for v in np.nanpercentile(vals, [2.5, 97.5]))
        rows.append(row)
    return pd.DataFrame(rows)


def decompose(d: pd.DataFrame, pg: pd.DataFrame) -> dict:
    valid = pg[pg.valid]
    pooled = auroc(d.label, d.score)
    if abs(pooled - roc_auc_score(d.label, d.score)) > 1e-12:
        raise RuntimeError("rank-based AUROC disagrees with sklearn")
    between = auroc(d.label, d.groupby("gene").score.transform("mean"))
    return {"pooled_auroc": pooled, "macro_within_gene_auroc": float(valid.auroc.mean()),
            "n_weighted_within_gene_auroc": float((valid.n * valid.auroc).sum() / valid.n.sum()), "between_gene_auroc": between,
            "n_valid_genes": int(len(valid)), "frac_variants_in_valid_genes": float(valid.n.sum() / pg.n.sum())}


def main() -> None:
    # ---- Step 2 outputs: join check + gene check ------------------------------------------------------------------
    jc = pl.join_check_table()
    print("Join check:\n", jc.to_string(index=False, float_format=lambda v: f"{v:.4f}"), flush=True)
    if (jc.status != "OK").any():
        raise SystemExit("join check FAILED; refusing to continue")
    jc_view = jc.rename(columns={"run": "Run", "split": "Split", "n_rows": "Rows", "labels_match": "Labels match", "pooled_auroc_recomputed": "AUROC recomputed",
                                 "auroc_in_final_json": "AUROC in final JSON", "match_4dp": "Match (4 dp)", "status": "Status"})
    for c in ("AUROC recomputed", "AUROC in final JSON"):
        jc_view[c] = jc_view[c].map(f4)
    gr.write_table(jc_view, "join_check", "Join Check: Saved Predictions vs. Reconstructed Test Set",
                   footnote="Genes are attached by row order; every row's label was asserted equal to the reconstructed test label. "
                            "The DNABERT-2 warmup-0 run is not used. " + SKIP_NOTE)
    gene_info = pl.splitb_gene_check()
    print(f"\nSPLIT-B: {gene_info['n_test_genes']} test genes, {len(gene_info['seen_in_train'])} seen in train: {', '.join(gene_info['test_genes'])}", flush=True)

    # ---- Step 3: per-gene metrics on SPLIT-B ----------------------------------------------------------------------
    runs_b = [("M3", "m3", None), ("NT", "nt", None)] + [(f"DNABERT-2 seed {s}", "dnabert2", s) for s in SEEDS]
    data_b = {name: pl.load_predictions(model, "B", seed) for name, model, seed in runs_b}
    any_run = data_b["M3"]
    order = any_run.groupby("gene").size().sort_values(ascending=False, kind="stable").index.tolist()
    pd.DataFrame({"gene": order, "n_test": [int((any_run.gene == g).sum()) for g in order],
                  "n_pathogenic": [int(any_run[any_run.gene == g].label.sum()) for g in order],
                  "seen_in_train": 0}).to_csv(OUT / "splitB_test_genes.csv", index=False)

    long_b = []
    for name, model, seed in runs_b:
        pg = per_gene(data_b[name], order)
        pg.insert(0, "run", name)
        pg.insert(0, "model", model)
        long_b.append(pg)
    long_b = pd.concat(long_b, ignore_index=True)
    long_b.to_csv(OUT / "per_gene_metrics_splitB.csv", index=False)

    def one(name: str) -> pd.DataFrame:
        return long_b[long_b.run == name].set_index("gene").loc[order]

    def view_single(name: str, stem: str, title: str) -> None:
        t = one(name)
        df = pd.DataFrame({"Gene": order, "n": t.n.to_numpy(), "n pathogenic": t.n_pathogenic.to_numpy(), "n benign": t.n_benign.to_numpy(),
                           "AUROC": [f4(a) if v else "insufficient" for a, v in zip(t.auroc, t.valid)],
                           "95% CI (bootstrap)": [f"{f4(lo)}-{f4(hi)}" if v else "-" for lo, hi, v in zip(t.ci_lo, t.ci_hi, t.valid)],
                           "Accuracy@0.5": t.accuracy.map(f4).to_numpy(), "Pred. pos. frac": t.pred_pos_frac.map(f4).to_numpy()})
        gr.write_table(df, stem, title, footnote=f"{INSUFF_NOTE} 95% CI: {N_BOOT} bootstrap resamples of the gene's variants, seed 0.")

    view_single("M3", "view_per_gene_splitB_m3", "M3 (gene-only): Per-Gene Metrics on SPLIT-B")
    view_single("NT", "view_per_gene_splitB_nt", "NT: Per-Gene Metrics on SPLIT-B")

    seeds_tbl = {s: one(f"DNABERT-2 seed {s}") for s in SEEDS}
    auc_mat = pd.DataFrame({s: seeds_tbl[s].auroc for s in SEEDS}).loc[order]
    acc_mat = pd.DataFrame({s: seeds_tbl[s].accuracy for s in SEEDS}).loc[order]
    ppf_mat = pd.DataFrame({s: seeds_tbl[s].pred_pos_frac for s in SEEDS}).loc[order]
    valid = seeds_tbl[42].valid.loc[order]
    d2 = pd.DataFrame({"gene": order, "n": seeds_tbl[42].n.loc[order].to_numpy(), "valid": valid.to_numpy(),
                       "auroc_mean": auc_mat.mean(axis=1).to_numpy(), "auroc_min": auc_mat.min(axis=1).to_numpy(), "auroc_max": auc_mat.max(axis=1).to_numpy(),
                       "accuracy_mean": acc_mat.mean(axis=1).to_numpy(), "accuracy_min": acc_mat.min(axis=1).to_numpy(), "accuracy_max": acc_mat.max(axis=1).to_numpy(),
                       "pred_pos_mean": ppf_mat.mean(axis=1).to_numpy(), "pred_pos_min": ppf_mat.min(axis=1).to_numpy(), "pred_pos_max": ppf_mat.max(axis=1).to_numpy()})
    d2.to_csv(OUT / "per_gene_dnabert2_seed_summary_splitB.csv", index=False)
    s42 = seeds_tbl[42]
    df = pd.DataFrame({"Gene": order, "n": s42.n.loc[order].to_numpy(), "n pathogenic": s42.n_pathogenic.loc[order].to_numpy(), "n benign": s42.n_benign.loc[order].to_numpy()})
    for s in SEEDS:
        df[f"AUROC s{s}"] = [f4(a) if v else "insufficient" for a, v in zip(auc_mat[s], valid)]
    df["AUROC mean"] = [f4(a) if v else "-" for a, v in zip(d2.auroc_mean, valid)]
    df["AUROC min"] = [f4(a) if v else "-" for a, v in zip(d2.auroc_min, valid)]
    df["AUROC max"] = [f4(a) if v else "-" for a, v in zip(d2.auroc_max, valid)]
    df["95% CI (s42)"] = [f"{f4(lo)}-{f4(hi)}" if v else "-" for lo, hi, v in zip(s42.ci_lo.loc[order], s42.ci_hi.loc[order], valid)]
    df["Acc@0.5 mean"] = d2.accuracy_mean.map(f4).to_numpy()
    df["Pred. pos. mean"] = d2.pred_pos_mean.map(f4).to_numpy()
    gr.write_table(df, "view_per_gene_splitB_dnabert2", f"DNABERT-2{WARMUP_MARK}: Per-Gene Metrics on SPLIT-B (seeds 42/1/2/3, warmup 0.1)",
                   footnote=f"{WARMUP_NOTE} {INSUFF_NOTE} Per-seed CIs are in per_gene_metrics_splitB.csv ({N_BOOT} bootstrap resamples, seed 0).")

    # ---- Step 4: decomposition ------------------------------------------------------------------------------------
    dec = []
    for name, _, _ in runs_b:
        dec.append({"run": name, **decompose(data_b[name], long_b[long_b.run == name].reset_index(drop=True))})
    dec = pd.DataFrame(dec)
    d2_rows = dec[dec.run.str.startswith("DNABERT-2")].drop(columns="run")
    agg = [{"run": f"DNABERT-2 {stat} over seeds", **getattr(d2_rows, stat)().to_dict()} for stat in ("mean", "min", "max")]
    dec_all = pd.concat([dec, pd.DataFrame(agg)], ignore_index=True)
    dec_all.to_csv(OUT / "decomposition_splitB.csv", index=False)
    short = lambda r: (r.replace("DNABERT-2 seed", f"DNABERT-2{WARMUP_MARK} seed") if r.startswith("DNABERT-2 seed")
                       else f"DNABERT-2{WARMUP_MARK} " + r.split()[1]  if r.startswith("DNABERT-2") else r)
    dv = pd.DataFrame({"Run": [short(r) for r in dec_all.run],
                       "Pooled": dec_all.pooled_auroc.map(f4), "Within-gene (macro)": dec_all.macro_within_gene_auroc.map(f4),
                       "Within-gene (n-wtd)": dec_all.n_weighted_within_gene_auroc.map(f4), "Between-gene": dec_all.between_gene_auroc.map(f4),
                       "Valid genes": [f"{v:.0f}" if isinstance(v, float) and not float(v).is_integer() else str(int(v)) for v in dec_all.n_valid_genes],
                       "Variants in valid genes": dec_all.frac_variants_in_valid_genes.map(lambda v: f"{100 * v:.1f}%")})
    gr.write_table(dv, "view_decomposition_splitB", "SPLIT-B: Pooled vs. Within-Gene vs. Between-Gene AUROC",
                   footnote=f"{WARMUP_NOTE} All columns except the last two are AUROC; mean/min/max are over the 4 seeds; macro = unweighted mean over valid genes, n-wtd = weighted by gene size. Between-gene: each variant's score replaced by its gene's mean score, AUROC over all variants. "
                            f"Within-gene metrics use genes with >= {MIN_PER_CLASS} variants in each class. {SKIP_NOTE}")

    # ---- Step 5: same per-gene AUROC on SPLIT-A for the same 14 genes ---------------------------------------------
    runs_a = [("M3", "m3"), ("NT", "nt"), ("DNABERT-2", "dnabert2")]
    long_a = []
    for name, model in runs_a:
        d = pl.load_predictions(model, "A", None)
        d = d[d.gene.isin(order)]
        pg = per_gene(d, order)
        pg.insert(0, "run", name)
        long_a.append(pg)
    long_a = pd.concat(long_a, ignore_index=True)
    long_a.to_csv(OUT / "per_gene_metrics_splitA_14genes.csv", index=False)
    ref = long_a[long_a.run == "M3"].set_index("gene").loc[order]
    va = pd.DataFrame({"Gene": order, "n (SPLIT-A test)": ref.n.to_numpy(), "n pathogenic": ref.n_pathogenic.to_numpy(), "n benign": ref.n_benign.to_numpy()})
    for name, _ in runs_a:
        t = long_a[long_a.run == name].set_index("gene").loc[order]
        va[f"{name} AUROC"] = [f4(a) if v else "insufficient" for a, v in zip(t.auroc, t.valid)]
    gr.write_table(va, "view_per_gene_splitA_14genes", "Per-Gene AUROC on SPLIT-A Test Set (the same 14 genes as SPLIT-B)",
                   footnote=f"{INSUFF_NOTE} SPLIT-A has no per-seed DNABERT-2 runs; CIs are in per_gene_metrics_splitA_14genes.csv. {SKIP_NOTE}")

    # ---- Step 6: figures ------------------------------------------------------------------------------------------
    figure1(order, one, d2, valid)
    figure2(order, long_a, d2)

    # ---- printed summary ------------------------------------------------------------------------------------------
    print("\nSPLIT-B decomposition (AUROC):", flush=True)
    print(dec_all.to_string(index=False, float_format=lambda v: f"{v:.4f}"), flush=True)
    comp = pd.DataFrame({"gene": order, "n_B": one("M3").n.loc[order].to_numpy(), "n_A": ref.n.to_numpy()})
    comp["B: M3"], comp["B: NT"] = [one("M3").auroc.loc[order].to_numpy(), one("NT").auroc.loc[order].to_numpy()]
    comp["B: DNABERT-2 mean"], comp["B: DNABERT-2 min"], comp["B: DNABERT-2 max"] = d2.auroc_mean, d2.auroc_min, d2.auroc_max
    for name, _ in runs_a:
        comp[f"A: {name}"] = long_a[long_a.run == name].set_index("gene").loc[order].auroc.to_numpy()
    print("\nper-gene AUROC (NaN = insufficient); n_B / n_A = test variants of the gene in SPLIT-B / SPLIT-A:", flush=True)
    print(comp.to_string(index=False, float_format=lambda v: f"{v:.4f}"), flush=True)
    print("\nFiles in", OUT.relative_to(ROOT), ":", *sorted(p.name for p in OUT.iterdir()), sep="\n  ")


def figure1(order, one, d2, valid) -> None:
    cols = ["M3", "NT", f"DNABERT-2{WARMUP_MARK}\nmean (min-max), 4 seeds"]
    mat = np.full((len(order), 3), np.nan)
    ann = np.full((len(order), 3), "", dtype=object)
    for j, name in enumerate(["M3", "NT"]):
        t = one(name)
        for i, g in enumerate(order):
            if t.valid.loc[g]:
                mat[i, j], ann[i, j] = t.auroc.loc[g], f"{t.auroc.loc[g]:.3f}"
    for i, g in enumerate(order):
        if valid.iloc[i]:
            mat[i, 2] = d2.auroc_mean.iloc[i]
            ann[i, 2] = f"{d2.auroc_mean.iloc[i]:.3f}\n({d2.auroc_min.iloc[i]:.3f}-{d2.auroc_max.iloc[i]:.3f})"
    n_by_gene = one("M3").n
    labels = [f"{g}  (n={n_by_gene.loc[g]:,})" for g in order]
    fig, ax = plt.subplots(figsize=(8.8, 8.2))
    sns.heatmap(mat, annot=ann, fmt="", cmap="RdBu", vmin=0.3, vmax=1.0, center=0.5, mask=np.isnan(mat), linewidths=0.6, linecolor="white",
                cbar_kws={"label": "Within-gene AUROC (0.5 = chance)"}, xticklabels=cols, yticklabels=labels, ax=ax, annot_kws={"fontsize": 9})
    for i, g in enumerate(order):
        for j in range(3):
            if np.isnan(mat[i, j]):
                ax.add_patch(plt.Rectangle((j, i), 1, 1, facecolor="#e6e6e6", edgecolor="white", lw=0.6))
                ax.text(j + 0.5, i + 0.5, "insufficient", ha="center", va="center", fontsize=8, color="dimgray", style="italic")
    ax.set_title("Within-Gene AUROC on SPLIT-B (14 test genes)", fontsize=13, weight="bold")
    ax.tick_params(axis="y", rotation=0)
    fig.tight_layout(rect=(0, 0.10, 1, 1))
    foot = f"{WARMUP_NOTE}\n{INSUFF_NOTE} {SKIP_NOTE}"
    fig.text(0.5, 0.045, "\n".join(textwrap.fill(p, 120) for p in foot.split("\n")), ha="center", va="center", fontsize=8, style="italic", color="dimgray")
    fig.savefig(OUT / "fig1_within_gene_auroc_heatmap_splitB.png", dpi=220)
    plt.close(fig)


def place_labels(fig, ax, names, xs, ys, sizes, ylo=None, yhi=None, fontsize: float = 8.5) -> None:
    """Greedy label placement (in points): biggest markers first; avoid other labels, other markers and the axes edge."""
    k = 72 / fig.dpi
    P = ax.transData.transform(np.c_[xs, ys]) * k
    R = np.sqrt(np.asarray(sizes) / np.pi)  # marker radius in points (scatter s = area in pt^2)
    ex = ax.get_window_extent()
    bounds = (ex.x0 * k, ex.y0 * k, ex.x1 * k, ex.y1 * k)
    overlap = lambda a, b: max(0.0, min(a[2], b[2]) - max(a[0], b[0])) * max(0.0, min(a[3], b[3]) - max(a[1], b[1]))  # noqa: E731
    err_boxes = []
    if ylo is not None:
        lo_pt = ax.transData.transform(np.c_[xs, ylo]) * k
        hi_pt = ax.transData.transform(np.c_[xs, yhi]) * k
        err_boxes = [(P[j, 0] - 3, lo_pt[j, 1] - 2, P[j, 0] + 3, hi_pt[j, 1] + 2) for j in range(len(names))]  # error bar + caps
    placed: list[tuple] = []
    dirs = [(1, 1), (1, -1), (-1, 1), (-1, -1), (0, 1), (0, -1), (1, 0), (-1, 0)]
    for i in sorted(range(len(names)), key=lambda i: -sizes[i]):
        w, h = 0.62 * fontsize * len(names[i]), 1.3 * fontsize
        best = None
        for extra in (2, 14, 28, 44):
            for dx, dy in dirs:
                cx = P[i, 0] + dx * (R[i] + extra + w / 2) * (1 if dy or True else 0)
                cy = P[i, 1] + dy * (R[i] + extra + h / 2)
                box = (cx - w / 2, cy - h / 2, cx + w / 2, cy + h / 2)
                score = sum(1000 + overlap(box, b) for b in placed if overlap(box, b) > 0)
                for j in range(len(names)):
                    c = (P[j, 0] - R[j], P[j, 1] - R[j], P[j, 0] + R[j], P[j, 1] + R[j])
                    if overlap(box, c) > 0:
                        score += 200 + overlap(box, c)
                score += sum(300 + overlap(box, e) for e in err_boxes if overlap(box, e) > 0)
                if box[0] < bounds[0] or box[2] > bounds[2] or box[1] < bounds[1] or box[3] > bounds[3]:
                    score += 5000
                score += extra * 2 + (6 if (dx and dy) else 0)  # prefer side/above/below placement over diagonals
                if best is None or score < best[0]:
                    best = (score, cx, cy, box, extra)
            if best[0] < 60:  # a clean spot at this distance: stop searching farther out
                break
        _, cx, cy, box, extra = best
        placed.append(box)
        ax.annotate(names[i], (xs[i], ys[i]), xytext=(cx - P[i, 0], cy - P[i, 1]), textcoords="offset points", ha="center", va="center", fontsize=fontsize,
                    arrowprops=dict(arrowstyle="-", color="dimgray", lw=0.9, shrinkA=0, shrinkB=R[i]))


def figure2(order, long_a, d2) -> None:
    b_n = d2.set_index("gene").n
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.8))
    omitted, label_jobs = {}, []
    for ax, name in zip(axes, ["DNABERT-2", "NT"]):
        a = long_a[long_a.run == name].set_index("gene")
        if name == "NT":
            bt = pd.read_csv(OUT / "per_gene_metrics_splitB.csv")
            bt = bt[bt.run == "NT"].set_index("gene")
            y = bt.auroc
            lo = hi = None
            b_valid = bt.valid
        else:
            d = d2.set_index("gene")
            y, lo, hi, b_valid = d.auroc_mean, d.auroc_min, d.auroc_max, d.valid
        pts = [g for g in order if a.valid.loc[g] and b_valid.loc[g]]
        omitted[name] = [g for g in order if g not in pts]
        xs, ys = a.auroc.loc[pts].to_numpy(), y.loc[pts].to_numpy()
        sizes = 30 + 650 * (b_n.loc[pts].to_numpy() / b_n.max())
        lims = (min(xs.min(), ys.min()) - 0.06, 1.0)
        ax.plot([lims[0], 1.0], [lims[0], 1.0], color="dimgray", linestyle="--", linewidth=1.1, zorder=1, label="y = x")
        ax.axhline(0.5, color="lightgray", linewidth=0.8, zorder=0)
        ax.axvline(0.5, color="lightgray", linewidth=0.8, zorder=0)
        if lo is not None:
            ax.errorbar(xs, ys, yerr=[ys - lo.loc[pts].to_numpy(), hi.loc[pts].to_numpy() - ys], fmt="none", ecolor="black", elinewidth=1.1, capsize=3, zorder=3)
        ax.scatter(xs, ys, s=sizes, color="#2171b5" if name == "DNABERT-2" else "#a50f15", alpha=0.75, edgecolor="white", linewidth=0.8, zorder=2)
        ax.set_xlim(lims)
        ax.set_ylim(lims)
        label_jobs.append((ax, pts, xs, ys, sizes, None if lo is None else lo.loc[pts].to_numpy(), None if hi is None else hi.loc[pts].to_numpy()))
        ax.set_xlabel("Within-gene AUROC on SPLIT-A (same genes)")
        ax.set_ylabel("Within-gene AUROC on SPLIT-B" + (" (mean over seeds; bars = min-max)" if name == "DNABERT-2" else ""))
        ax.set_title(f"{name}{WARMUP_MARK if name == 'DNABERT-2' else ''}: SPLIT-A vs. SPLIT-B", fontsize=12, weight="bold")
        ax.grid(alpha=0.25)
        ax.legend(loc="lower right", fontsize=9)
    fig.suptitle("Per-Gene AUROC: SPLIT-A vs. SPLIT-B (marker area ~ SPLIT-B variants per gene)", fontsize=14, weight="bold")
    om = "; ".join(f"{k}: {', '.join(v) if v else 'none'}" for k, v in omitted.items())
    foot = f"{WARMUP_NOTE}\nGenes omitted (insufficient on SPLIT-A or SPLIT-B): {om}. {SKIP_NOTE}"
    fig.tight_layout(rect=(0, 0.09, 1, 0.94))
    fig.text(0.5, 0.04, "\n".join(textwrap.fill(p, 170) for p in foot.split("\n")), ha="center", va="center", fontsize=8.5, style="italic", color="dimgray")
    fig.canvas.draw()
    for job in label_jobs:
        place_labels(fig, *job)
    fig.savefig(OUT / "fig2_gene_auroc_splitA_vs_splitB.png", dpi=220)
    plt.close(fig)
    print("Figure 2 omitted genes:", omitted, flush=True)


if __name__ == "__main__":
    main()
