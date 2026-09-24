#!/usr/bin/env python
"""D7: how much does each model's score depend on train/test locus overlap?

Locus definition (reused exactly from the ad-hoc computation behind the 36.1% figure in docs/terminology.md):
  locus = chrom + "_" + pos; a test variant is LOCUS-SHARED when its locus also occurs among that split's TRAIN variants
  (any alt allele), otherwise LOCUS-NOVEL. SPLIT-A must give 12,186 / 33,786 = 36.1% (asserted), SPLIT-B 0.

Uses scripts/pred_loader.py (validated join of saved predictions to genes). No model inference. CPU only. Add-only.
Models: M3, NT, DNABERT-2 seeds 42/1/2/3 (warmup 0.1 on SPLIT-B; A/C single run, fixed LR). CNN/BiLSTM/CNN+BiLSTM/Ensemble have no
saved probabilities and are skipped. The DNABERT-2 warmup-0 run is never loaded. Numbers only; no mechanism is claimed.

    python scripts/run_d7_locus_overlap.py --before-sha <sha256 manifest> --before-ckpt <checkpoint size+mtime manifest>
"""
from __future__ import annotations

import os

os.environ["CUDA_VISIBLE_DEVICES"] = ""

import argparse
import csv
import hashlib
import io
import platform
import sys
import textwrap
from datetime import datetime
from pathlib import Path

sys.dont_write_bytecode = True
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import pred_loader as pl  # noqa: E402
import src.evaluation.generate_results as gr  # noqa: E402
from src.evaluation.generate_split_comparison import SPLIT_COLORS, WARMUP_MARK, WARMUP_NOTE  # noqa: E402

# ---- parameters, all fixed before any result was looked at --------------------------------------------------------
MIN_PER_CLASS = 30                       # per-class count for a gene to count as valid (Step 4: in BOTH subsets)
N_BOOT = 1000                            # bootstrap resamples
BOOT_SEED = 0                            # numpy default_rng(0), re-created per (run, split, subset) cell
CI_PCT = (2.5, 97.5)                     # 95% percentile interval
TOP_GENES = 5
EXPECTED_LOCUS = {"A": (12186, 33786, 36.1), "B": (0, 33780, 0.0), "C": (11495, 33536, 34.3)}  # (shared, n_test, pct)
SUBSETS = {"A": ["full", "shared", "novel"], "B": ["full"], "C": ["full", "shared", "novel"]}
RUNS = [("M3", "m3", None, "ABC"), ("NT", "nt", None, "ABC"), ("DNABERT-2 seed 42", "dnabert2", 42, "ABC"),
        ("DNABERT-2 seed 1", "dnabert2", 1, "B"), ("DNABERT-2 seed 2", "dnabert2", 2, "B"), ("DNABERT-2 seed 3", "dnabert2", 3, "B")]
CONTROL_RUNS = ["M3", "NT", "DNABERT-2 seed 42"]
SKIPPED = ["CNN", "BiLSTM", "CNN+BiLSTM", "Ensemble"]
SKIP_NOTE = "Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities)."
DESC_NOTE = "Subsets are not random samples, so differences are descriptive."
OUT = ROOT / "results" / "d7_locus_overlap"
HANDOFF = ROOT / "results" / "handoff"
gr.TABLE_DIR = OUT  # write_table resolves TABLE_DIR at call time
f4 = lambda v: "-" if v is None or pd.isna(v) else f"{v:.4f}"  # noqa: E731


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


# ---- locus labelling ----------------------------------------------------------------------------------------------
def _locus_of_vid() -> pd.Series:
    d = pl._data()
    return pd.Series((d.chrom.astype(str) + "_" + d.pos.astype(str)).to_numpy(), index=d.vid.to_numpy())


LOCUS = _locus_of_vid()


def train_loci(split: str) -> set:
    _, train, _ = pl._manifest_split(split)
    return set(train.chrom.astype(str) + "_" + train.pos.astype(str))


def shared_flags(frame: pd.DataFrame, split: str) -> np.ndarray:
    return LOCUS.loc[frame.vid].isin(train_loci(split)).to_numpy()


def manifest_frame(split: str) -> pd.DataFrame:
    fr = pl.test_frame(split, "manifest").copy()
    fr["shared"] = shared_flags(fr, split)
    return fr


# ---- AUROC and bootstrap ------------------------------------------------------------------------------------------
def auroc(y, s) -> float:
    y = np.asarray(y)
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(s)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def auroc_ci(y, s) -> tuple[float, float, float]:
    y, s = np.asarray(y), np.asarray(s)
    n = len(y)
    rng = np.random.default_rng(BOOT_SEED)
    vals = np.empty(N_BOOT)
    for i in range(N_BOOT):
        idx = rng.integers(0, n, n)
        vals[i] = auroc(y[idx], s[idx])
    lo, hi = np.nanpercentile(vals, CI_PCT)
    return auroc(y, s), float(lo), float(hi)


def load_run(model: str, split: str, seed) -> pd.DataFrame:
    d = pl.load_predictions(model, split, seed)
    fr = pl.test_frame(split, pl._order_kind(model, split))
    assert (d.gene.to_numpy() == fr.gene.to_numpy()).all() and (d.label.to_numpy() == fr.label.to_numpy()).all()
    d["shared"] = shared_flags(fr, split)
    return d


def subset_mask(d: pd.DataFrame, subset: str) -> np.ndarray:
    return {"full": np.ones(len(d), bool), "shared": d.shared.to_numpy(), "novel": ~d.shared.to_numpy()}[subset]


# ---- integrity check ----------------------------------------------------------------------------------------------
def integrity(before_sha: Path | None, before_ckpt: Path | None) -> list[list]:
    if before_sha is None or before_ckpt is None:
        return [["integrity_check", "not run (no before-manifest given)"]]

    def sha(p: Path) -> str:
        h = hashlib.sha256()
        with open(p, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 22), b""):
                h.update(chunk)
        return h.hexdigest()

    before = {l.rstrip("\n").split("  ", 1)[1]: l.split("  ", 1)[0] for l in open(before_sha) if l.strip()}
    now = {}
    for top in ("results", "scripts", "docs"):
        for p in sorted((ROOT / top).rglob("*")):
            rel = p.relative_to(ROOT).as_posix()
            if p.is_file() and p.suffix != ".pyc" and not rel.startswith("results/checkpoints/"):
                now[rel] = sha(p)
    modified = sorted(k for k in before if k in now and now[k] != before[k])
    removed = sorted(k for k in before if k not in now)
    added = sorted(k for k in now if k not in before)
    cb = {}
    for l in open(before_ckpt):
        if l.strip():
            path, size, mtime = l.strip().rsplit(" ", 2)
            cb[path] = (int(size), float(mtime))
    ca = {}
    for p in (ROOT / "results" / "checkpoints").rglob("*"):
        if p.is_file():
            st = p.stat()
            ca[p.relative_to(ROOT).as_posix()] = (st.st_size, st.st_mtime)
    c_changed = sorted(k for k in cb if k in ca and (ca[k][0] != cb[k][0] or abs(ca[k][1] - cb[k][1]) > 1e-3))
    c_added, c_removed = sorted(set(ca) - set(cb)), sorted(set(cb) - set(ca))
    ok = not (modified or removed or c_changed or c_added or c_removed)
    return [["files_hashed_sha256_before", len(before)], ["files_hashed_sha256_now", len(now)],
            ["modified_preexisting_files", "none" if not modified else ";".join(modified)], ["removed_files", "none" if not removed else ";".join(removed)],
            ["added_files", ";".join(added)], ["checkpoints_size_mtime_before", len(cb)], ["checkpoints_size_mtime_now", len(ca)],
            ["checkpoints_changed_added_removed", "none" if not (c_changed or c_added or c_removed) else f"{c_changed};{c_added};{c_removed}"],
            ["result", "PASS: no pre-existing file changed" if ok else "FAIL"]]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--before-sha", type=Path)
    ap.add_argument("--before-ckpt", type=Path)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    HANDOFF.mkdir(parents=True, exist_ok=True)

    # ---- Step 1: reproduce and extend the locus-overlap fraction ----------------------------------------------------
    frames = {s: manifest_frame(s) for s in "ABC"}
    ov_rows = []
    for s in "ABC":
        n, k = len(frames[s]), int(frames[s].shared.sum())
        exp_k, exp_n, exp_pct = EXPECTED_LOCUS[s]
        pct = round(100 * k / n, 1)
        if (k, n, pct) != (exp_k, exp_n, exp_pct):
            raise AssertionError(f"SPLIT-{s}: locus-shared {k}/{n} = {pct}% differs from the recorded {exp_k}/{exp_n} = {exp_pct}%")
        ov_rows.append([s, n, k, n - k, k / n, 100 * k / n, exp_pct, exp_k, True])
    ov = pd.DataFrame(ov_rows, columns=["split", "n_test", "n_locus_shared", "n_locus_novel", "frac_shared", "pct_shared", "expected_pct", "expected_n_shared", "matches_expected"])
    ov.to_csv(OUT / "locus_overlap_fraction.csv", index=False)
    gr.write_table(pd.DataFrame({"Split": ov.split, "Test variants": ov.n_test, "Locus-shared": ov.n_locus_shared, "Locus-novel": ov.n_locus_novel,
                                 "Fraction shared": ov.frac_shared.map(f4), "Recorded in docs (%)": ov.expected_pct.map(lambda v: f"{v:.1f}")}),
                   "view_locus_overlap_fraction", "Locus Overlap: Test Variants Sharing a Locus with a Training Variant",
                   footnote="locus = chrom_pos; shared = the same locus occurs among the split's training variants (any alt allele). Reproduces the 36.1% / 34.3% figures in docs/terminology.md.")

    # ---- Step 2: subset composition ---------------------------------------------------------------------------------
    comp_rows, comp_view = [], []
    for s in "ABC":
        fr = frames[s]
        for sub in SUBSETS[s]:
            m = fr[subset_mask(fr, sub)]
            top = m.gene.value_counts().head(TOP_GENES)
            comp_rows.append([s, sub, len(m), int(m.label.sum()), m.label.mean(), m.gene.nunique(), ";".join(f"{g}:{c}" for g, c in top.items())])
            comp_view.append({"Split": s, "Subset": sub, "n": len(m), "n pathogenic": int(m.label.sum()), "Pathogenic frac.": f4(m.label.mean()), "Genes": m.gene.nunique(),
                              **{f"Top {i + 1}": (f"{g} ({c:,})" if i < len(top) else "-") for i, (g, c) in enumerate(list(top.items()) + [(None, None)] * (TOP_GENES - len(top)))}})
    comp = pd.DataFrame(comp_rows, columns=["split", "subset", "n", "n_pathogenic", "pathogenic_fraction", "n_genes", "top5_genes_count"])
    comp.to_csv(OUT / "subset_composition.csv", index=False)
    gr.write_table(pd.DataFrame(comp_view), "view_subset_composition", "Test-Set Subset Composition (locus-shared vs. locus-novel)",
                   footnote="SPLIT-B is gene-disjoint: every test variant is locus-novel, so only the full set is shown. " + DESC_NOTE)

    # ---- Step 3: AUROC with bootstrap CI per model and subset -------------------------------------------------------
    cells, data_cache = {}, {}
    for name, model, seed, splits in RUNS:
        for s in splits:
            d = load_run(model, s, seed)
            data_cache[(name, s)] = d
            for sub in SUBSETS[s]:
                m = subset_mask(d, sub)
                a, lo, hi = auroc_ci(d.label.to_numpy()[m], d.score.to_numpy()[m])
                cells[(name, s, sub)] = {"n": int(m.sum()), "n_pos": int(d.label.to_numpy()[m].sum()), "auroc": a, "lo": lo, "hi": hi}
    long_rows = [[n, s, sub, c["n"], c["n_pos"], c["auroc"], c["lo"], c["hi"]] for (n, s, sub), c in cells.items()]
    d2_b = [cells[(f"DNABERT-2 seed {k}", "B", "full")] for k in (42, 1, 2, 3)]
    for stat, fn in (("mean", np.mean), ("min", np.min), ("max", np.max)):
        long_rows.append([f"DNABERT-2 {stat} over seeds", "B", "full", d2_b[0]["n"], d2_b[0]["n_pos"], float(fn([c["auroc"] for c in d2_b])), None, None])
    pd.DataFrame(long_rows, columns=["run", "split", "subset", "n", "n_pathogenic", "auroc", "ci_low", "ci_high"]).to_csv(OUT / "auroc_by_subset_long.csv", index=False)

    cols = [("A", "full"), ("A", "shared"), ("A", "novel"), ("B", "full"), ("C", "full"), ("C", "shared"), ("C", "novel")]
    run_names = [r[0] for r in RUNS] + [f"DNABERT-2 {x} over seeds" for x in ("mean", "min", "max")]
    lookup = {(r[0], r[1], r[2]): r for r in long_rows}

    def label(r: str) -> str:
        return r.replace(" over seeds", "").replace("DNABERT-2", "DNABERT-2" + WARMUP_MARK) if r.startswith("DNABERT-2") else r

    pt, ci, diffs = [], [], []
    for r in run_names:
        row_pt, row_ci = {"Run": label(r)}, {"Run": label(r)}
        for s, sub in cols:
            hit = lookup.get((r, s, sub))
            row_pt[f"{s} {sub}"] = f4(hit[5]) if hit else "-"
            row_ci[f"{s} {sub}"] = f"{f4(hit[6])}-{f4(hit[7])}" if (hit and hit[6] is not None) else "-"
        for s in "AC":
            sh, no = lookup.get((r, s, "shared")), lookup.get((r, s, "novel"))
            row_pt[f"{s} shared-novel"] = f4(sh[5] - no[5]) if (sh and no) else "-"
            if sh and no:
                diffs.append([r, s, sh[5], no[5], sh[5] - no[5]])
        pt.append(row_pt)
        ci.append(row_ci)
    order = ["Run", "A full", "A shared", "A novel", "A shared-novel", "B full", "C full", "C shared", "C novel", "C shared-novel"]
    gr.write_table(pd.DataFrame(pt)[order], "view_auroc_by_subset", "ROC-AUC by Locus-Overlap Subset",
                   footnote=f"{WARMUP_NOTE} DNABERT-2 SPLIT-A and SPLIT-C have one run each (seed 42, fixed LR); seeds 1-3 exist for SPLIT-B only, so mean/min/max are over the four SPLIT-B seeds. "
                            f"'shared-novel' is the point difference. {DESC_NOTE} {SKIP_NOTE}")
    gr.write_table(pd.DataFrame(ci[:len(RUNS)])[["Run"] + order[1:4] + order[5:9]], "view_auroc_ci_by_subset", "ROC-AUC 95% Bootstrap CI by Locus-Overlap Subset",
                   footnote=f"CI: {N_BOOT} bootstrap resamples of the subset's variants, seed {BOOT_SEED}, percentile interval. {WARMUP_NOTE} {DESC_NOTE}")
    pd.DataFrame(diffs, columns=["run", "split", "auroc_shared", "auroc_novel", "difference"]).to_csv(OUT / "shared_minus_novel.csv", index=False)

    # ---- Step 4: gene-composition control ---------------------------------------------------------------------------
    gene_rows, detail_rows, ctl_rows, ctl_view, genes_used = [], [], [], [], {}
    for s in "AC":
        fr = frames[s]
        g_used = []
        for g, gd in fr.groupby("gene"):
            cnt = {sub: (int(gd.label[m].sum()), int((1 - gd.label[m]).sum())) for sub, m in (("shared", gd.shared.to_numpy()), ("novel", ~gd.shared.to_numpy()))}
            if all(v >= MIN_PER_CLASS for pair in cnt.values() for v in pair):
                g_used.append(g)
                gene_rows.append([s, g, sum(cnt["shared"]), cnt["shared"][0], cnt["shared"][1], sum(cnt["novel"]), cnt["novel"][0], cnt["novel"][1]])
        genes_used[s] = sorted(g_used)
    for run in CONTROL_RUNS:
        for s in "AC":
            d = data_cache[(run, s)]
            per = {"shared": [], "novel": []}
            for g in genes_used[s]:
                for sub in ("shared", "novel"):
                    m = (d.gene.to_numpy() == g) & subset_mask(d, sub)
                    a = auroc(d.label.to_numpy()[m], d.score.to_numpy()[m])
                    per[sub].append(a)
                    detail_rows.append([run, s, g, sub, int(m.sum()), a])
            ms, mn = float(np.mean(per["shared"])), float(np.mean(per["novel"]))
            ctl_rows.append([run, s, len(genes_used[s]), ms, mn, ms - mn])
            ctl_view.append({"Run": run, "Split": s, "Genes used": len(genes_used[s]), "Macro shared": f4(ms), "Macro novel": f4(mn), "Shared-novel": f4(ms - mn)})
    pd.DataFrame(gene_rows, columns=["split", "gene", "n_shared", "n_shared_pathogenic", "n_shared_benign", "n_novel", "n_novel_pathogenic", "n_novel_benign"]).to_csv(OUT / "within_gene_control_genes.csv", index=False)
    pd.DataFrame(detail_rows, columns=["run", "split", "gene", "subset", "n", "auroc"]).to_csv(OUT / "within_gene_control_per_gene_auroc.csv", index=False)
    pd.DataFrame(ctl_rows, columns=["run", "split", "n_genes_used", "macro_within_gene_auroc_shared", "macro_within_gene_auroc_novel", "difference"]).to_csv(OUT / "within_gene_control.csv", index=False)
    gr.write_table(pd.DataFrame(ctl_view), "view_within_gene_control", "Gene-Composition Control: Macro Within-Gene AUROC, Locus-Shared vs. Locus-Novel",
                   footnote=f"Genes with >= {MIN_PER_CLASS} pathogenic and >= {MIN_PER_CLASS} benign test variants in BOTH subsets. M3 is constant within a gene, so its within-gene AUROC is 0.5. "
                            f"DNABERT-2 SPLIT-A/C: single run, seed 42, fixed LR. {DESC_NOTE} {SKIP_NOTE}")
    gv = pd.DataFrame(gene_rows, columns=["Split", "Gene", "n shared", "shared path.", "shared benign", "n novel", "novel path.", "novel benign"])
    for s_ in "AC":
        gr.write_table(gv[gv.Split == s_].reset_index(drop=True), f"view_within_gene_control_genes_split{s_}", f"Genes Used in the Gene-Composition Control: SPLIT-{s_}",
                       footnote=f"Eligibility: >= {MIN_PER_CLASS} pathogenic and >= {MIN_PER_CLASS} benign variants in both subsets. " + DESC_NOTE)

    # ---- Figure: models x [A shared, A novel, B] --------------------------------------------------------------------
    figure(cells, d2_b)

    # ---- parameters + export ----------------------------------------------------------------------------------------
    params = [["min_class_count", MIN_PER_CLASS], ["n_bootstrap", N_BOOT], ["bootstrap_seed", BOOT_SEED], ["ci", "95% percentile (2.5, 97.5)"], ["top_genes_listed", TOP_GENES],
              ["locus_definition", "locus = chrom + '_' + pos; shared = locus occurs among the split's training variants (any alt allele); novel = otherwise"],
              ["models", "M3; NT; DNABERT-2 seeds 42/1/2/3 (SPLIT-B warmup 0.1); DNABERT-2 SPLIT-A/C single run seed 42, fixed LR"], ["python", platform.python_version()],
              ["generated", datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z")]]
    pd.DataFrame(params, columns=["key", "value"]).to_csv(OUT / "parameters.csv", index=False)
    integ = integrity(args.before_sha, args.before_ckpt)
    write_export(ov, comp, long_rows, diffs, ctl_rows, gene_rows, detail_rows, params, integ)

    a = lambda r, s, sub: cells[(r, s, sub)]["auroc"]  # noqa: E731
    print(f"report: {(HANDOFF / 'd7_export.md').relative_to(ROOT)}  (tables/figure: {OUT.relative_to(ROOT)}/)")
    print(f"1 locus check: A {ov.pct_shared[0]:.1f}% ({ov.n_locus_shared[0]:,}/{ov.n_test[0]:,}) matched; B {ov.n_locus_shared[1]}; C {ov.pct_shared[2]:.1f}%; skipped: {', '.join(SKIPPED)}")
    for r in ("M3", "NT", "DNABERT-2 seed 42"):
        print(f"  {r}: A shared {a(r, 'A', 'shared'):.4f} / novel {a(r, 'A', 'novel'):.4f} | C shared {a(r, 'C', 'shared'):.4f} / novel {a(r, 'C', 'novel'):.4f}")
    b = [c["auroc"] for c in d2_b]
    print(f"  B full: M3 {a('M3', 'B', 'full'):.4f}, NT {a('NT', 'B', 'full'):.4f}, DNABERT-2 seeds mean {np.mean(b):.4f} (min-max {min(b):.4f}-{max(b):.4f})")
    for r in ("NT", "DNABERT-2 seed 42"):
        cc = {s: next(x for x in ctl_rows if x[0] == r and x[1] == s) for s in "AC"}
        print(f"  control {r}: A ({cc['A'][2]} genes) shared {cc['A'][3]:.4f} / novel {cc['A'][4]:.4f}; C ({cc['C'][2]} genes) shared {cc['C'][3]:.4f} / novel {cc['C'][4]:.4f}")
    print(f"integrity: {[r for r in integ if r[0] == 'result'][0][1] if any(r[0] == 'result' for r in integ) else integ[0][1]}")


def figure(cells: dict, d2_b: list) -> None:
    models = [("M3 (gene-only)", "M3"), (f"NT", "NT"), (f"DNABERT-2{WARMUP_MARK}", "DNABERT-2 seed 42")]
    series = [("SPLIT-A locus-shared", "A", "shared", "#08519c"), ("SPLIT-A locus-novel", "A", "novel", "#9ecae1"), ("SPLIT-B (all locus-novel)", "B", "full", SPLIT_COLORS["B"])]
    x, width = np.arange(len(models)), 0.26
    fig, ax = plt.subplots(figsize=(10.5, 6.6))
    for i, (lab, s, sub, color) in enumerate(series):
        for j, (_, run) in enumerate(models):
            c = cells[(run, s, sub)]
            xp = x[j] + (i - 1) * width
            bar = ax.bar(xp, c["auroc"], width, color=color, edgecolor="black" if (run == "M3" and s == "B") else "none", hatch="//" if (run == "M3" and s == "B") else None, linewidth=0.8)
            ax.errorbar(xp, c["auroc"], yerr=[[c["auroc"] - c["lo"]], [c["hi"] - c["auroc"]]], fmt="none", ecolor="black", elinewidth=1.1, capsize=3, zorder=5)
            mark = WARMUP_MARK if (run.startswith("DNABERT-2") and s == "B") else ""
            ax.text(xp, c["hi"] + 0.008, f"{c['auroc']:.3f}{mark}", ha="center", va="bottom", fontsize=8)
    ax.axhline(0.5, color="dimgray", linestyle="--", linewidth=1.1, zorder=1)
    ax.set_xticks(x, [m[0] for m in models])
    ax.set_ylabel("ROC-AUC (validation)")
    ax.set_ylim(0, 1.08)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    ax.grid(axis="y", alpha=0.25)
    ax.set_title("ROC-AUC by Locus Overlap: SPLIT-A Locus-Shared vs. Locus-Novel, and SPLIT-B", fontsize=13, weight="bold")
    handles = [Patch(facecolor=c, label=l) for l, _, _, c in series] + [Line2D([0], [0], color="dimgray", linestyle="--", label="AUROC 0.5 (chance)")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.07), ncol=4, fontsize=9, frameon=False)
    b = [c["auroc"] for c in d2_b]
    foot = (f"{WARMUP_NOTE} Error bars: bootstrap 95% CI within each subset ({N_BOOT} resamples, seed {BOOT_SEED}). DNABERT-2 SPLIT-B bar = seed 42; seeds 42/1/2/3 range {min(b):.3f}-{max(b):.3f}. "
            f"M3 SPLIT-B is a chance floor by construction (no test gene seen in train). {DESC_NOTE} {SKIP_NOTE}")
    fig.tight_layout(rect=(0, 0.13, 1, 1))
    fig.text(0.5, 0.06, "\n".join(textwrap.fill(p, 150) for p in foot.split("\n")), ha="center", va="center", fontsize=8, style="italic", color="dimgray")
    fig.savefig(OUT / "d7_grouped_bar_auroc_subsets.png", dpi=220)
    plt.close(fig)


def write_export(ov, comp, long_rows, diffs, ctl_rows, gene_rows, detail_rows, params, integ) -> None:
    secs = [block("1. locus_reproduction_check",
                  "Locus-shared test variants per split; expected values are those recorded when 36.1% / 34.3% were first computed (docs/terminology.md). The script asserts every match.",
                  ["split", "n_test", "n_locus_shared", "n_locus_novel", "frac_shared", "pct_shared", "expected_pct", "expected_n_shared", "matches_expected"], ov.values.tolist()),
            block("2. skipped_models", "Models without saved probabilities; no predictions were regenerated (no-inference rule). The DNABERT-2 warmup-0 run is never loaded.",
                  ["model", "reason"], [[m, "no saved probabilities (only weights and per-epoch history)"] for m in SKIPPED]),
            block("3. subset_composition", "Test-set subsets. top5_genes_count = gene:variant count, ordered by count.",
                  ["split", "subset", "n", "n_pathogenic", "pathogenic_fraction", "n_genes", "top5_genes_count"], comp.values.tolist()),
            block("4. auroc_by_subset", "AUROC within each subset; ci = 95% percentile bootstrap (1000 resamples of the subset's variants, seed 0). DNABERT-2 SPLIT-A/C: one run (seed 42, fixed LR); "
                  "seeds 1/2/3 exist for SPLIT-B only (warmup 0.1). mean/min/max rows are over the four SPLIT-B seeds (no CI).",
                  ["run", "split", "subset", "n", "n_pathogenic", "auroc", "ci_low", "ci_high"], long_rows),
            block("5. shared_minus_novel", "Point difference in AUROC, locus-shared minus locus-novel (descriptive; subsets are not random samples).",
                  ["run", "split", "auroc_shared", "auroc_novel", "difference"], diffs),
            block("6. within_gene_control", f"Macro-mean within-gene AUROC over genes with >= {MIN_PER_CLASS} pathogenic and >= {MIN_PER_CLASS} benign test variants in BOTH subsets (gene lists in section 7).",
                  ["run", "split", "n_genes_used", "macro_within_gene_auroc_shared", "macro_within_gene_auroc_novel", "difference"], ctl_rows),
            block("7. within_gene_control_genes", "Genes used in section 6, with per-subset counts.",
                  ["split", "gene", "n_shared", "n_shared_pathogenic", "n_shared_benign", "n_novel", "n_novel_pathogenic", "n_novel_benign"], gene_rows),
            block("8. within_gene_control_per_gene_auroc", "Per-gene AUROC behind section 6.", ["run", "split", "gene", "subset", "n", "auroc"], detail_rows),
            block("9. parameters", "Thresholds and resample counts, fixed in the script before any result was inspected.", ["key", "value"], params),
            block("10. integrity_check", "sha256 of every file under results/, scripts/, docs/ before vs. after (results/checkpoints: size + mtime only). Added files are this task's outputs.", ["check", "value"], integ)]
    header = ("# D7 locus-overlap export\n\n"
              f"Generated {datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %z')}. Numbers only; no mechanism is claimed. Floats have 6 decimals; blank = not defined. "
              "Subsets are not random samples, so differences are descriptive.\n\n"
              "DNABERT-2 SPLIT-B results (seeds 42/1/2/3) used LR warmup 0.1 after the initial fixed-LR SPLIT-B attempt collapsed; DNABERT-2 SPLIT-A and SPLIT-C used the original fixed-LR schedule.\n\n")
    (HANDOFF / "d7_export.md").write_text(header + "\n".join(secs))


if __name__ == "__main__":
    main()
