#!/usr/bin/env python
"""D6: sequence similarity between SPLIT-B test windows and train windows (different genes, by construction of the split).

Window = ref_sequence (the 100-bp reference window oriented to the gene strand). Step 2: exact duplicates (forward and reverse complement)
of test windows among train windows, for ref_sequence and mut_sequence. Step 3: each window -> the set of its distinct 25-mers, encoded
exactly as 50-bit integers (a collision-free 'hash'); an inverted index over train windows (sorted k-mer array + searchsorted, no all-pairs
comparison) gives, per test window, the maximum fraction of its distinct 25-mers shared with any single train window. Two variants:
'forward' (k-mers as written; primary) and 'canonical' (min of k-mer and its reverse complement; strand-independent). K-mers containing N are
skipped. Read-only on existing files, CPU only, add-only. Numbers only; no interpretation.

    python scripts/run_d6_cross_gene_similarity.py --before-sha <sha256 manifest> --before-ckpt <checkpoint size+mtime manifest>
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
from collections import defaultdict
from datetime import datetime
from pathlib import Path

sys.dont_write_bytecode = True
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from src.data.splits import SPLITS_DIR, apply_split, load_manifest, variant_ids  # noqa: E402
from src.evaluation.generate_split_comparison import SPLIT_COLORS  # noqa: E402

# ---- parameters, fixed before any result was looked at ------------------------------------------------------------
K = 25
THRESHOLDS = (0.5, 0.8, 0.95)
TOP_PAIRS = 15
QUANTILES = (0.0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1.0)
DATA = ROOT / "data" / "processed" / "variants_100bp.parquet"
MANIFEST = SPLITS_DIR / "split_manifest_b.json"
SEQ_COLS = ["ref_sequence", "mut_sequence", "genomic_ref_sequence", "genomic_mut_sequence"]
WINDOW_COL = "ref_sequence"
OUT = ROOT / "results" / "d6_cross_gene_similarity"
HANDOFF = ROOT / "results" / "handoff"
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


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


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


# ---- k-mer machinery ----------------------------------------------------------------------------------------------
LUT = np.full(256, 4, dtype=np.uint8)
for _i, _ch in enumerate("ACGT"):
    LUT[ord(_ch)] = _i
    LUT[ord(_ch.lower())] = _i
SENTINEL = np.uint64(0xFFFFFFFFFFFFFFFF)


def encode(seqs: pd.Series) -> np.ndarray:
    lens = seqs.str.len().to_numpy()
    if not (lens == lens[0]).all():
        raise AssertionError("windows have unequal lengths; equal-length ref_sequence expected")
    return LUT[np.frombuffer("".join(seqs.tolist()).encode(), dtype=np.uint8).reshape(len(seqs), int(lens[0]))]


def distinct_kmers(bases: np.ndarray, canonical: bool) -> tuple[np.ndarray, np.ndarray]:
    """Per window, the distinct valid K-mers. Returns (flat k-mer codes, per-window counts); row-major order."""
    n, L = bases.shape
    m = L - K + 1
    fwd = np.zeros((n, m), dtype=np.uint64)
    rc = np.zeros((n, m), dtype=np.uint64)
    bad = np.zeros((n, m), dtype=bool)
    for j in range(K):
        b = bases[:, j:j + m]
        bad |= b == 4
        b64 = np.where(b == 4, 0, b).astype(np.uint64)
        fwd |= np.left_shift(b64, np.uint64(2 * (K - 1 - j)))
        rc |= np.left_shift(np.uint64(3) - b64, np.uint64(2 * j))
    codes = np.minimum(fwd, rc) if canonical else fwd
    codes[bad] = SENTINEL
    codes.sort(axis=1)
    keep = codes != SENTINEL
    keep[:, 1:] &= codes[:, 1:] != codes[:, :-1]
    return codes[keep], keep.sum(axis=1)


def similarity(test_bases, train_bases, canonical: bool) -> dict:
    tr_k, tr_c = distinct_kmers(train_bases, canonical)
    tr_w = np.repeat(np.arange(len(tr_c)), tr_c)
    order = np.argsort(tr_k, kind="stable")
    sk, sw = tr_k[order], tr_w[order]
    te_k, te_c = distinct_kmers(test_bases, canonical)
    starts = np.concatenate([[0], np.cumsum(te_c)])
    lo = np.searchsorted(sk, te_k, side="left")
    hi = np.searchsorted(sk, te_k, side="right")
    n = len(te_c)
    frac = np.full(n, np.nan)
    best = np.full(n, -1, dtype=np.int64)
    ties = np.zeros(n, dtype=np.int64)
    tied_lists: list[np.ndarray] = [np.empty(0, dtype=np.int64)] * n
    for i in range(n):
        a, b = starts[i], starts[i + 1]
        lens = hi[a:b] - lo[a:b]
        tot = int(lens.sum())
        if te_c[i] == 0:
            continue
        if tot == 0:
            frac[i] = 0.0
            continue
        idx = np.repeat(lo[a:b] - (np.cumsum(lens) - lens), lens) + np.arange(tot)
        uniq, cnt = np.unique(sw[idx], return_counts=True)
        mx = int(cnt.max())
        tied = uniq[cnt == mx]
        frac[i] = mx / te_c[i]
        best[i], ties[i] = int(tied[0]), len(tied)
        tied_lists[i] = tied
    return {"frac": frac, "best": best, "ties": ties, "tied": tied_lists, "n_kmers": te_c, "train_kmer_entries": int(len(sk))}


# ---- integrity check ----------------------------------------------------------------------------------------------
def integrity(before_sha: Path | None, before_ckpt: Path | None) -> list[list]:
    if before_sha is None or before_ckpt is None:
        return [["result", "integrity check not run (no before-manifest given)"]]

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
    ca = {p.relative_to(ROOT).as_posix(): (p.stat().st_size, p.stat().st_mtime) for p in (ROOT / "results" / "checkpoints").rglob("*") if p.is_file()}
    c_changed = sorted(k for k in cb if k in ca and (ca[k][0] != cb[k][0] or abs(ca[k][1] - cb[k][1]) > 1e-3))
    c_added, c_removed = sorted(set(ca) - set(cb)), sorted(set(cb) - set(ca))
    ok = not (modified or removed or c_changed or c_added or c_removed)
    return [["files_hashed_sha256_before", len(before)], ["files_hashed_sha256_now", len(now)],
            ["modified_preexisting_files", "none" if not modified else ";".join(modified)], ["removed_files", "none" if not removed else ";".join(removed)],
            ["added_files_count", len(added)], ["added_files", ";".join(added)], ["checkpoints_size_mtime_before", len(cb)], ["checkpoints_size_mtime_now", len(ca)],
            ["checkpoints_changed_added_removed", "none" if not (c_changed or c_added or c_removed) else f"{c_changed};{c_added};{c_removed}"],
            ["result", "PASS: no pre-existing file changed" if ok else "FAIL"]]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--before-sha", type=Path)
    ap.add_argument("--before-ckpt", type=Path)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    HANDOFF.mkdir(parents=True, exist_ok=True)

    # ---- Step 1: load, split, length distribution -------------------------------------------------------------------
    df = pd.read_parquet(DATA, columns=["chrom", "pos", "ref", "alt", "gene", "label", "strand"] + SEQ_COLS)
    manifest = load_manifest(MANIFEST)
    train, test = apply_split(df, manifest)
    train, test = train.reset_index(drop=True), test.reset_index(drop=True)
    assert set(train.gene).isdisjoint(set(test.gene)), "SPLIT-B is not gene-disjoint"
    len_rows, stat_rows = [], []
    for col in SEQ_COLS:
        ln = df[col].str.len()
        for length, cnt in ln.value_counts().sort_index().items():
            len_rows.append([col, int(length), int(cnt)])
        has_n = df[col].str.upper().str.contains("N").sum()
        stat_rows.append([col, len(df), int(ln.min()), float(ln.median()), float(ln.mean()), int(ln.max()), int((ln == 100).sum()), int(has_n)])
    pd.DataFrame(len_rows, columns=["column", "length", "n_windows"]).to_csv(OUT / "length_distribution.csv", index=False)
    stats = pd.DataFrame(stat_rows, columns=["column", "n", "min_len", "median_len", "mean_len", "max_len", "n_len_100", "n_with_N"])
    stats.to_csv(OUT / "length_summary.csv", index=False)
    comp = [["train", len(train), train.gene.nunique(), int(train.label.sum())], ["test", len(test), test.gene.nunique(), int(test.label.sum())]]

    # ---- Step 2: exact duplicates (forward and reverse complement) --------------------------------------------------
    dup_rows = []
    for col in ("ref_sequence", "mut_sequence"):
        tr_set = set(train[col])
        fwd = test[col].isin(tr_set).to_numpy()
        rc = test[col].map(revcomp).isin(tr_set).to_numpy()
        dup_rows.append([col, len(test), int(fwd.sum()), int(rc.sum()), int((fwd | rc).sum())])
    dups = pd.DataFrame(dup_rows, columns=["window_column", "n_test_windows", "n_exact_in_train_forward", "n_revcomp_in_train", "n_either"])
    dups.to_csv(OUT / "exact_duplicates.csv", index=False)

    # ---- Step 3: 25-mer shingle similarity --------------------------------------------------------------------------
    te_b, tr_b = encode(test[WINDOW_COL]), encode(train[WINDOW_COL])
    results = {}
    for name, canon in (("forward", False), ("canonical", True)):
        results[name] = similarity(te_b, tr_b, canon)
    train_genes = train.gene.to_numpy()
    train_vids = variant_ids(train).to_numpy()
    test_vids = variant_ids(test).to_numpy()
    per = pd.DataFrame({"vid": test_vids, "gene": test.gene.to_numpy(), "label": test.label.to_numpy()})
    dist_rows, pair_rows, gene_rows = [], [], []
    for name, r in results.items():
        f = r["frac"]
        per[f"n_distinct_kmers_{name}"] = r["n_kmers"]
        per[f"max_shared_frac_{name}"] = f
        per[f"best_train_vid_{name}"] = [train_vids[b] if b >= 0 else "" for b in r["best"]]
        per[f"best_train_gene_{name}"] = [train_genes[b] if b >= 0 else "" for b in r["best"]]
        per[f"n_tied_train_windows_{name}"] = r["ties"]
        ok = ~np.isnan(f)
        dist_rows.append([name, "n_test_windows_scored", int(ok.sum()), None])
        for q in QUANTILES:
            dist_rows.append([name, f"quantile_{q:g}", float(np.nanquantile(f, q)), None])
        dist_rows.append([name, "mean", float(np.nanmean(f)), None])
        for t in THRESHOLDS:
            k = int((f >= t).sum())
            dist_rows.append([name, f"n_ge_{t:g}", k, k / int(ok.sum())])
        for t in THRESHOLDS:
            hits = np.flatnonzero(f >= t)
            weight: dict[tuple, float] = defaultdict(float)
            per_gene: dict[str, int] = defaultdict(int)
            for i in hits:
                genes = sorted(set(train_genes[r["tied"][i]]))
                for g in genes:
                    weight[(test.gene.iloc[i], g)] += 1.0 / len(genes)
                per_gene[test.gene.iloc[i]] += 1
            total = float(len(hits))
            for (tg, trg), w in sorted(weight.items(), key=lambda x: -x[1])[:TOP_PAIRS]:
                pair_rows.append([name, t, tg, trg, w, w / total if total else None, int(total), len(weight)])
            for tg, c in sorted(per_gene.items(), key=lambda x: -x[1]):
                gene_rows.append([name, t, tg, c, c / total, int((test.gene == tg).sum())])
    dist = pd.DataFrame(dist_rows, columns=["variant", "statistic", "value", "fraction_of_scored_windows"])
    pairs = pd.DataFrame(pair_rows, columns=["variant", "threshold", "test_gene", "train_gene", "n_hits_weighted", "share_of_hits", "total_hits_at_threshold", "n_distinct_pairs_at_threshold"])
    gene_hits = pd.DataFrame(gene_rows, columns=["variant", "threshold", "test_gene", "n_test_windows_hit", "share_of_hits", "n_test_windows_in_gene"])
    top_windows = per.sort_values("max_shared_frac_forward", ascending=False, kind="stable").head(20).reset_index(drop=True)
    top_rows = [[i + 1, r.vid, r.gene, int(r.label), r.max_shared_frac_forward, r.best_train_vid_forward, r.best_train_gene_forward, r.n_tied_train_windows_forward, r.max_shared_frac_canonical,
                 r.best_train_gene_canonical] for i, r in top_windows.iterrows()]
    top_cols = ["rank", "test_vid", "test_gene", "label", "max_shared_frac_forward", "best_train_vid_forward", "best_train_gene_forward", "n_tied_train_windows_forward",
                "max_shared_frac_canonical", "best_train_gene_canonical"]
    pd.DataFrame(top_rows, columns=top_cols).to_csv(OUT / "top_similarity_windows.csv", index=False)
    dist.to_csv(OUT / "similarity_distribution.csv", index=False)
    pairs.to_csv(OUT / "gene_pairs_high_similarity.csv", index=False)
    gene_hits.to_csv(OUT / "test_gene_hit_counts.csv", index=False)
    per.to_csv(OUT / "similarity_per_test_window.csv", index=False)

    figure(results)

    params = [["k", K], ["window_column", WINDOW_COL], ["sequence_columns_found", ";".join(SEQ_COLS)], ["thresholds", ";".join(str(t) for t in THRESHOLDS)], ["top_pairs_listed", TOP_PAIRS],
              ["kmer_encoding", "exact 2-bit-per-base integers (50 bits for K=25): a collision-free hash; k-mers containing N are skipped; distinct k-mers per window (set semantics)"],
              ["variants", "forward = k-mers as written (primary); canonical = min(k-mer, reverse complement)"],
              ["similarity", "max over single train windows of |distinct 25-mers of the test window found in that train window| / |distinct valid 25-mers of the test window|"],
              ["index", "inverted index over train windows: sorted k-mer array + searchsorted posting lists (no all-pairs comparison)"],
              ["gene_pair_attribution", "each high-similarity test window is attributed to the distinct train genes among its tied best train windows, weight 1/(number of distinct genes)"],
              ["split", f"data/splits/split_manifest_b.json ({manifest['content_hash']})"], ["train_windows", len(train)], ["test_windows", len(test)],
              ["gene_disjoint_assert", "test genes and train genes are disjoint (asserted)"], ["python", platform.python_version()],
              ["generated", datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z")]]
    pd.DataFrame(params, columns=["key", "value"]).to_csv(OUT / "parameters.csv", index=False)
    integ = integrity(args.before_sha, args.before_ckpt)
    na = [["af_coverage_and_bucket_cutoffs", "not applicable", "D3 items in the generic reporting note; D6 uses no allele frequency"],
          ["models", "not applicable", "D6 scores no model; it compares sequences only"],
          ["bucketing_or_skipped_groups", "none", "no test window was excluded; windows with no valid 25-mer would be unscored (count in similarity_distribution: n_test_windows_scored)"]]
    secs = [block("1. parameters", "Fixed before any result was inspected.", ["key", "value"], params),
            block("2. split_composition", "SPLIT-B train and test windows (gene-disjoint by construction; asserted).", ["set", "n_windows", "n_genes", "n_pathogenic"], comp),
            block("3. length_summary", "Length statistics of the four sequence columns over all 168,927 variants. ref_sequence is the window used for similarity.",
                  list(stats.columns), stats.values.tolist()),
            block("4. length_distribution", "Full length distribution of each sequence column.", ["column", "length", "n_windows"], len_rows),
            block("5. exact_duplicates", "SPLIT-B test windows whose exact sequence (forward) or reverse complement also occurs among train windows.", list(dups.columns), dups.values.tolist()),
            block("6. similarity_distribution", "Distribution of the per-test-window maximum shared 25-mer fraction, and the count / fraction of test windows at or above each threshold.",
                  list(dist.columns), dist.values.tolist()),
            block("7. gene_pairs_high_similarity", f"Top {TOP_PAIRS} (test gene, train gene) pairs per variant and threshold by weighted hits; share_of_hits = n_hits_weighted / total_hits_at_threshold.",
                  list(pairs.columns), pairs.values.tolist()),
            block("8. test_gene_hit_counts", "Number of test windows at or above each threshold, per test gene.", list(gene_hits.columns), gene_hits.values.tolist()),
            block("9. top_similarity_windows", "Descriptive: the 20 test windows with the highest forward max shared 25-mer fraction, with their best train window and gene (not an additional threshold).", top_cols, top_rows),
            block("10. not_applicable_or_skipped", "Items in the generic reporting note that do not apply to D6.", ["item", "status", "note"], na),
            block("11. integrity_check", "sha256 of every file under results/, scripts/, docs/ before vs. after (results/checkpoints: size + mtime only). Added files are this task's outputs.", ["check", "value"], integ)]
    header = ("# D6 cross-gene sequence-similarity export\n\n" f"Generated {datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %z')}. Numbers only; no interpretation. Floats have 6 decimals; blank = not defined.\n\n"
              "Per-window results (33,780 rows) are in results/d6_cross_gene_similarity/similarity_per_test_window.csv and are not embedded here.\n\n")
    (HANDOFF / "d6_export.md").write_text(header + "\n".join(secs))

    g = lambda v, s: dist[(dist.variant == v) & (dist.statistic == s)].iloc[0]  # noqa: E731
    print(f"report: {(HANDOFF / 'd6_export.md').relative_to(ROOT)}  (csv/figure: {OUT.relative_to(ROOT)}/)")
    print(f"1 windows: train {len(train):,} / {train.gene.nunique()} genes, test {len(test):,} / {test.gene.nunique()} genes; ref_sequence length min {int(stats.iloc[0].min_len)} max {int(stats.iloc[0].max_len)}")
    for _, r in dups.iterrows():
        print(f"2 exact duplicates ({r.window_column}): forward {r.n_exact_in_train_forward}, reverse complement {r.n_revcomp_in_train}, either {r.n_either} of {r.n_test_windows}")
    for v in ("forward", "canonical"):
        print(f"3 {v}: median {g(v, 'quantile_0.5').value:.4f}, p95 {g(v, 'quantile_0.95').value:.4f}, p99 {g(v, 'quantile_0.99').value:.4f}, max {g(v, 'quantile_1').value:.4f} | "
              + ", ".join(f">={t:g}: {int(g(v, f'n_ge_{t:g}').value)} ({g(v, f'n_ge_{t:g}').fraction_of_scored_windows:.6f})" for t in THRESHOLDS))
    for v in ("forward",):
        top = pairs[(pairs.variant == v) & (pairs.threshold == 0.5)].head(3)
        print("4 top pairs (forward, >=0.5): " + "; ".join(f"{a}->{b} {w:.1f}" for a, b, w in zip(top.test_gene, top.train_gene, top.n_hits_weighted)) if len(top) else "4 no windows >= 0.5 (forward)")
    print(f"5 integrity: {[x for x in integ if x[0] == 'result'][0][1]}")


def figure(results: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.6), sharey=True)
    bins = np.linspace(0, 1, 51)
    for ax, (name, r) in zip(axes, results.items()):
        f = r["frac"][~np.isnan(r["frac"])]
        ax.hist(f, bins=bins, color=SPLIT_COLORS["B"], edgecolor="white", alpha=0.9)
        for t in THRESHOLDS:
            ax.axvline(t, color="dimgray", linestyle="--", linewidth=1)
            ax.text(t, 0.97, f"{t:g}", transform=ax.get_xaxis_transform(), ha="right", va="top", fontsize=8, color="dimgray", rotation=90)
        ax.set_yscale("log")
        ax.set_xlim(0, 1)
        ax.set_xlabel("Max shared 25-mer fraction with any one train window")
        ax.set_title(f"SPLIT-B test windows, {name} 25-mers (n={len(f):,})", fontsize=11, weight="bold")
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("Test windows (log scale)")
    fig.suptitle("Cross-Gene 25-mer Similarity: SPLIT-B Test vs. Train Windows", fontsize=14, weight="bold")
    foot = (f"Window = ref_sequence (100 bp). Distinct valid 25-mers per window; inverted index over train windows; similarity = shared distinct 25-mers of the test window in the single most similar "
            f"train window / distinct 25-mers of the test window. forward = as written; canonical = min(k-mer, reverse complement). Dashed lines: {', '.join(f'{t:g}' for t in THRESHOLDS)}.")
    fig.tight_layout(rect=(0, 0.09, 1, 0.94))
    fig.text(0.5, 0.04, "\n".join(textwrap.fill(p, 170) for p in foot.split("\n")), ha="center", va="center", fontsize=8.5, style="italic", color="dimgray")
    fig.savefig(OUT / "d6_similarity_histogram.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
