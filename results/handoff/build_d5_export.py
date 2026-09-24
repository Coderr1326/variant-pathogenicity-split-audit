#!/usr/bin/env python
"""Package the D5 per-gene numbers into ONE markdown file (results/handoff/d5_export.md).

Numbers only, no interpretation. Values are copied from results/d5_per_gene/*.csv (full precision) and re-verified at
export time against an independent sklearn recomputation (per-gene AUROC, n, class counts, decomposition). Bootstrap CIs are
copied, not recomputed. No model inference. Add-only: writes only results/handoff/d5_export.md.
"""
from __future__ import annotations

import os

os.environ["CUDA_VISIBLE_DEVICES"] = ""

import csv
import io
import json
import platform
import re
import subprocess
import sys
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from scipy.stats import kendalltau, spearmanr
from sklearn.metrics import roc_auc_score

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import pred_loader as pl  # noqa: E402

D5 = ROOT / "results" / "d5_per_gene"
OUT = ROOT / "results" / "handoff" / "d5_export.md"
MIN_PER_CLASS = 30


def fmt(v) -> str:
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return ""
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, (int, np.integer)):
        return str(int(v))
    if isinstance(v, (float, np.floating)):
        return f"{float(v):.6f}"
    return str(v)


def block(title: str, description: str, header: list[str], rows: list[list]) -> str:
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(header)
    for r in rows:
        assert len(r) == len(header), (title, r)
        w.writerow([fmt(x) for x in r])
    return f"## {title}\n\n{description}\n\n```csv\n{buf.getvalue()}```\n"


def reason(n1: int, n0: int) -> str:
    parts = []
    if n1 < MIN_PER_CLASS:
        parts.append(f"n_pathogenic={n1}<{MIN_PER_CLASS}")
    if n0 < MIN_PER_CLASS:
        parts.append(f"n_benign={n0}<{MIN_PER_CLASS}")
    return "; ".join(parts)


def main() -> None:
    sections = []
    b = pd.read_csv(D5 / "per_gene_metrics_splitB.csv")
    a = pd.read_csv(D5 / "per_gene_metrics_splitA_14genes.csv")
    seed_sum = pd.read_csv(D5 / "per_gene_dnabert2_seed_summary_splitB.csv")
    dec = pd.read_csv(D5 / "decomposition_splitB.csv")
    order = pd.read_csv(D5 / "splitB_test_genes.csv").gene.tolist()
    run_order_b = ["M3", "NT", "DNABERT-2 seed 42", "DNABERT-2 seed 1", "DNABERT-2 seed 2", "DNABERT-2 seed 3"]

    # ---- independent re-verification against the loader + sklearn -------------------------------------------------
    loaders = {"M3": ("m3", None), "NT": ("nt", None), "DNABERT-2 seed 42": ("dnabert2", 42), "DNABERT-2 seed 1": ("dnabert2", 1),
               "DNABERT-2 seed 2": ("dnabert2", 2), "DNABERT-2 seed 3": ("dnabert2", 3)}
    data = {run: pl.load_predictions(m, "B", s) for run, (m, s) in loaders.items()}
    n_checked = 0
    for run, d in data.items():
        rows = b[b.run == run].set_index("gene")
        pooled = roc_auc_score(d.label, d.score)
        valid_aucs, valid_n = [], []
        for g in order:
            gd = d[d.gene == g]
            r = rows.loc[g]
            assert (len(gd), int(gd.label.sum())) == (r.n, r.n_pathogenic), (run, g)
            if bool(r.valid):
                auc = roc_auc_score(gd.label, gd.score)
                assert abs(auc - r.auroc) < 1e-9, (run, g, auc, r.auroc)
                valid_aucs.append(auc)
                valid_n.append(len(gd))
                n_checked += 1
        dr = dec[dec.run == run].iloc[0]
        between = roc_auc_score(d.label, d.groupby("gene").score.transform("mean"))
        w = np.array(valid_n)
        assert abs(pooled - dr.pooled_auroc) < 1e-9 and abs(np.mean(valid_aucs) - dr.macro_within_gene_auroc) < 1e-9
        assert abs((w * np.array(valid_aucs)).sum() / w.sum() - dr.n_weighted_within_gene_auroc) < 1e-9 and abs(between - dr.between_gene_auroc) < 1e-9
    print(f"independent re-verification OK: {n_checked} per-gene AUROCs + 6 decompositions", flush=True)

    # ---- 1. join_check --------------------------------------------------------------------------------------------
    warmup_note = {("dnabert2", "B"): "0.1"}
    jrows = []
    for model, split, seed in pl.available_runs():
        pred_path, final_path = pl._paths(model, split, seed)
        pred = pd.read_parquet(pred_path)
        frame = pl.test_frame(split, pl._order_kind(model, split))
        assert len(pred) == len(frame)
        mism = int((pred.label.astype(int).to_numpy() != frame.label.astype(int).to_numpy()).sum())
        recomputed = float(roc_auc_score(pred.label, pred.probability_pathogenic))
        saved = float(json.loads(final_path.read_text())["metrics"]["roc_auc"])
        fj = json.loads(final_path.read_text())
        if model == "m3":
            seed_v, warm = "n/a", "n/a"
        else:
            seed_v = fj.get("seed", 42 if seed in (None, 42) else seed)
            warm = fj.get("warmup_ratio", warmup_note.get((model, split), "none"))
        jrows.append([{"m3": "M3", "nt": "NT", "dnabert2": "DNABERT-2"}[model], split, seed_v, warm, len(pred), mism, recomputed, saved, abs(recomputed - saved)])
    sections.append(block("1. join_check",
                          "One row per usable run. label_mismatches = rows where the saved label differs from the label reconstructed from the frozen manifest "
                          "(row-order join). warmup: 0.1 = linear LR warmup ratio; none = fixed LR, no scheduler; n/a = not applicable (M3 has no training). "
                          "DNABERT-2 SPLIT-B seed 42 predates the warmup_ratio JSON key (its checkpoint scheduler state confirms warmup 0.1).",
                          ["model", "split", "seed", "warmup", "n_rows", "label_mismatches", "pooled_auroc_recomputed", "pooled_auroc_from_final_json", "abs_diff"], jrows))

    # ---- 2. decomposition_splitB ----------------------------------------------------------------------------------
    drows = [[r.run, r.pooled_auroc, r.macro_within_gene_auroc, r.n_weighted_within_gene_auroc, r.between_gene_auroc, int(r.n_valid_genes), r.frac_variants_in_valid_genes]
             for r in dec.itertuples()]
    sections.append(block("2. decomposition_splitB",
                          "AUROC on the SPLIT-B test set. between_gene = AUROC after replacing each variant's score by its gene's mean score (all variants). "
                          "Within-gene metrics use genes with >= 30 pathogenic and >= 30 benign variants (n_valid_genes of 14). "
                          "mean/min/max rows are over DNABERT-2 seeds 42/1/2/3 (warmup 0.1).",
                          ["run", "pooled_auroc", "macro_within_gene", "n_weighted_within_gene", "between_gene", "n_valid_genes", "frac_variants_in_valid_genes"], drows))

    # ---- 3. per_gene_splitB ---------------------------------------------------------------------------------------
    hdr = ["run", "gene", "n", "n_pathogenic", "n_benign", "pathogenic_fraction", "auroc", "auroc_blank_reason", "ci_low", "ci_high",
           "accuracy_at_0.5", "predicted_positive_fraction"]

    def gene_rows(df: pd.DataFrame, runs: list[str]) -> list[list]:
        out = []
        for run in runs:
            t = df[df.run == run].set_index("gene")
            for g in order:
                r = t.loc[g]
                ok = bool(r.valid)
                out.append([run, g, int(r.n), int(r.n_pathogenic), int(r.n_benign), r.n_pathogenic / r.n, r.auroc if ok else None,
                            "" if ok else reason(int(r.n_pathogenic), int(r.n_benign)), r.ci_lo if ok else None, r.ci_hi if ok else None,
                            r.accuracy, r.pred_pos_frac])
        return out

    rows3 = gene_rows(b, run_order_b)
    s42 = b[b.run == "DNABERT-2 seed 42"].set_index("gene")
    ss = seed_sum.set_index("gene")
    for stat in ("mean", "min", "max"):
        for g in order:
            r = s42.loc[g]
            ok = bool(ss.loc[g, "valid"])
            rows3.append([f"DNABERT-2 {stat} over seeds", g, int(r.n), int(r.n_pathogenic), int(r.n_benign), r.n_pathogenic / r.n,
                          ss.loc[g, f"auroc_{stat}"] if ok else None, "" if ok else reason(int(r.n_pathogenic), int(r.n_benign)), None, None,
                          ss.loc[g, f"accuracy_{stat}"], ss.loc[g, f"pred_pos_{'mean' if stat == 'mean' else stat}"]])
    sections.append(block("3. per_gene_splitB",
                          "SPLIT-B test set, one row per (run, gene); 14 genes ordered by n descending. auroc is blank when the gene has fewer than 30 variants in either class "
                          "(reason given). ci_low/ci_high: 95% percentile bootstrap, 1000 resamples of the gene's variants, seed 0, per run. mean/min/max rows aggregate "
                          "DNABERT-2 seeds 42/1/2/3 metric-by-metric (min and max are taken independently per column); ci is blank for those rows.",
                          hdr, rows3))

    # ---- 4. per_gene_splitA_same_genes ----------------------------------------------------------------------------
    rows4 = gene_rows(a, ["M3", "NT", "DNABERT-2"])
    sections.append(block("4. per_gene_splitA_same_genes",
                          "SPLIT-A test set, the same 14 genes. SPLIT-A has one run per model (DNABERT-2 seed 42, no warmup). Same columns and conventions as section 3.",
                          hdr, rows4))

    # ---- 5. per_gene_comparison_splitB ----------------------------------------------------------------------------
    nt = b[b.run == "NT"].set_index("gene")
    rows5 = []
    for g in order:
        ok = bool(nt.loc[g].valid)
        rows5.append([g, int(nt.loc[g].n), nt.loc[g].auroc if ok else None, ss.loc[g, "auroc_mean"] if ok else None, ss.loc[g, "auroc_min"] if ok else None,
                      ss.loc[g, "auroc_max"] if ok else None, (nt.loc[g].auroc - ss.loc[g, "auroc_mean"]) if ok else None])
    sections.append(block("5. per_gene_comparison_splitB",
                          "Per-gene SPLIT-B AUROC. dnabert2_* columns are over seeds 42/1/2/3. Blank = insufficient class counts (see section 3).",
                          ["gene", "n", "nt_auroc", "dnabert2_mean_auroc", "dnabert2_min_auroc", "dnabert2_max_auroc", "nt_minus_dnabert2_mean"], rows5))

    # ---- 6. rank_agreement ----------------------------------------------------------------------------------------
    valid_genes = [g for g in order if bool(nt.loc[g].valid)]
    x = nt.loc[valid_genes].auroc.to_numpy()
    rows6 = []
    for label, y in [("NT vs DNABERT-2 mean over seeds", ss.loc[valid_genes, "auroc_mean"].to_numpy())] + \
                    [(f"NT vs DNABERT-2 seed {s}", b[b.run == f"DNABERT-2 seed {s}"].set_index("gene").loc[valid_genes].auroc.to_numpy()) for s in (42, 1, 2, 3)]:
        rows6.append([label, len(valid_genes), float(spearmanr(x, y).statistic), float(kendalltau(x, y).statistic)])
    sections.append(block("6. rank_agreement",
                          f"Correlation of per-gene SPLIT-B AUROC between NT and DNABERT-2 over the {len(valid_genes)} valid genes (excluded: "
                          f"{', '.join(g for g in order if g not in valid_genes)}). Spearman rho and Kendall tau-b (scipy defaults).",
                          ["comparison", "n_genes", "spearman_rho", "kendall_tau_b"], rows6))

    # ---- 7. gene_size_context -------------------------------------------------------------------------------------
    cap = int(re.search(r"cap_per_gene:\s*(\d+)", (ROOT / "REPRODUCTION_MANIFEST.yaml").read_text()).group(1))
    full = pd.read_parquet(ROOT / "data" / "processed" / "variants_100bp.parquet", columns=["gene"]).gene.value_counts()
    rows7 = []
    for g in order:
        n_b = int(nt.loc[g].n)
        assert int(full.loc[g]) == n_b, (g, full.loc[g], n_b)  # whole genes sit on one side of SPLIT-B
        rows7.append([g, n_b, int(full.loc[g]), cap, n_b == cap])
    sections.append(block("7. gene_size_context",
                          f"n_splitB_test = variants of the gene in the SPLIT-B test set; n_dataset_total = variants of the gene in the full 168,927-variant dataset "
                          f"(equal, because SPLIT-B assigns whole genes to one side). cap_per_gene = {cap} (REPRODUCTION_MANIFEST.yaml); hit_cap = n_dataset_total == cap.",
                          ["gene", "n_splitB_test", "n_dataset_total", "cap_per_gene", "hit_cap"], rows7))

    # ---- 8. excluded_runs -----------------------------------------------------------------------------------------
    wu0 = json.loads((ROOT / "results" / "metrics" / "dnabert2_100bp_splitB_seed1_warmup0_final.json").read_text())["metrics"]["roc_auc"]
    rows8 = []
    for m in ("CNN", "BiLSTM", "CNN+BiLSTM", "Ensemble"):
        rows8.append([m, "A;B;C", "no saved probabilities: only model weights (results/<model>_100bp*.pt) and per-epoch metric history exist; predictions were not regenerated (no-inference rule)"])
    rows8.append(["DNABERT-2 seed 1, warmup 0", "B", f"collapsed run (pooled AUROC {wu0:.4f}); excluded by instruction; its prediction file exists but pred_loader does not expose it"])
    rows8.append(["DNABERT-2 first SPLIT-B attempt (fixed LR, seed 42)", "B", "collapsed and discarded (accuracy 0.5377, AUROC 0.4999, tp=0 per results/tables/dnabert2_splitB_methodology_note.md); its files were overwritten by the warmup rerun, so no per-gene data exists"])
    rows8.append(["DNABERT-2 seeds 1, 2, 3", "A", "never run: seed runs exist for SPLIT-B only"])
    rows8.append(["M3; NT; DNABERT-2", "C", "outside the D5 scope (SPLIT-B and SPLIT-A only); prediction files exist"])
    sections.append(block("8. excluded_runs", "Every run that is not in sections 1-7, and why.", ["run", "splits", "reason"], rows8))

    # ---- 9. environment -------------------------------------------------------------------------------------------
    def git(*args: str) -> str:
        return subprocess.run(["git", "-C", str(ROOT.parent), *args], capture_output=True, text=True).stdout.strip()

    env = [["date", datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z")], ["python", platform.python_version()], ["platform", platform.platform()]]
    for pkg in ("numpy", "pandas", "scikit-learn", "scipy", "matplotlib", "seaborn", "torch", "transformers"):
        env.append([pkg, version(pkg)])
    env += [["git_head", git("rev-parse", "HEAD")], ["git_head_subject", git("log", "-1", "--format=%s")],
            ["git_note", "Test/ is largely untracked; results, scripts and figures are not covered by this commit hash"],
            ["bootstrap", "1000 resamples, numpy default_rng(0) per (run, gene)"], ["min_class_count_for_auroc", MIN_PER_CLASS]]
    sections.append(block("9. environment", "Software and repository state at export time.", ["key", "value"], env))

    header = ("# D5 per-gene export\n\n"
              f"Generated {datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %z')}. Numbers only; no interpretation. Floats have 6 decimals; blank = not defined.\n"
              "Source: results/d5_per_gene/*.csv (full precision), re-verified at export time against an independent sklearn recomputation "
              "(per-gene AUROC, n, class counts, pooled/macro/n-weighted/between-gene AUROC). Bootstrap CIs are copied, not recomputed.\n\n"
              "Sections: 1 join_check, 2 decomposition_splitB, 3 per_gene_splitB, 4 per_gene_splitA_same_genes, 5 per_gene_comparison_splitB, "
              "6 rank_agreement, 7 gene_size_context, 8 excluded_runs, 9 environment.\n\n"
              "DNABERT-2 SPLIT-B results (seeds 42/1/2/3) all used LR warmup 0.1 after the initial fixed-LR SPLIT-B attempt collapsed; DNABERT-2 SPLIT-A and SPLIT-C used the "
              "original fixed-LR schedule.\n\n")
    OUT.write_text(header + "\n".join(sections))
    print(f"wrote {OUT.relative_to(ROOT)}", flush=True)


if __name__ == "__main__":
    main()
