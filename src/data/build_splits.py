#!/usr/bin/env python
"""Build and freeze SPLIT-A/B/C train/val manifests for the gene-disjoint leakage audit.

SPLIT-A (random, existing S1 behavior): retroactively freezes what
``train_test_split(df, test_size=.2, random_state=42, stratify=df.label)`` currently
produces on data/processed/variants_100bp.parquet, for documentation/audit parity.
It is NOT wired into run_experiment.py's default path -- that stays exactly as-is.

SPLIT-B (gene-disjoint): assigns whole genes to train or val (no gene split across
both sides), greedily targeting an 80/20 variant-count split while nudging toward
SPLIT-A's ~54/46 benign/pathogenic balance on each side. Deterministic given seed 42
(used only to break ties among same-size genes).

SPLIT-C (temporal): neither ClinVar VCF in data/raw/ carries a per-variant date (only
a single file-wide ##fileDate) -- their INFO fields have no CLNDATE/last-evaluated
field. The only ClinVar release that publishes one is the separate tab-delimited
variant_summary.txt.gz, whose GRCh38 rows include LastEvaluated, PositionVCF,
ReferenceAlleleVCF and AlternateAlleleVCF. This joins that (date-matched,
2025-06 archived snapshot, same month as clinvar_20250630.vcf.gz) file to our
168,927 variants by (chrom, PositionVCF, ReferenceAlleleVCF, AlternateAlleleVCF),
then picks a cutoff date giving an ~80/20 split by count (train = older, val =
strictly newer).

Usage:
    python -m src.data.build_splits --split a
    python -m src.data.build_splits --split b
    python -m src.data.build_splits --split c
"""
from __future__ import annotations

import argparse
import random
from pathlib import Path

import pandas as pd
from sklearn.model_selection import train_test_split

from src.data.splits import ROOT, SPLITS_DIR, build_manifest, save_manifest, variant_ids

DATA = ROOT / "data" / "processed" / "variants_100bp.parquet"
VARIANT_SUMMARY = ROOT / "data" / "raw" / "variant_summary_2025-06.txt.gz"


def build_split_a(seed: int = 42) -> dict:
    df = pd.read_parquet(DATA)
    train_df, val_df = train_test_split(df, test_size=0.2, random_state=seed, stratify=df.label)
    train_ids = variant_ids(train_df).tolist()
    val_ids = variant_ids(val_df).tolist()
    extra = {
        "method": "random 80/20, stratified by label, sklearn train_test_split(random_state=42)",
        "train_label_balance": train_df.label.value_counts(normalize=True).round(4).to_dict(),
        "val_label_balance": val_df.label.value_counts(normalize=True).round(4).to_dict(),
    }
    manifest = build_manifest(
        split_id="a", name="random_80_20",
        description="Retroactive freeze of S1's existing random split (train_test_split, seed 42, stratified by label). "
                     "Documentation/audit only -- run_experiment.py's default path recomputes this inline and does not read this file.",
        source_data=str(DATA.relative_to(ROOT)), train_ids=train_ids, val_ids=val_ids, seed=seed, extra=extra,
    )
    return manifest


def _assign_genes(gene_stats: pd.DataFrame, total_n: int, overall_benign_frac: float, seed: int,
                   label_weight: float = 5.0, max_passes: int = 200) -> dict[str, str]:
    """Gene-disjoint bin-balancing in two deterministic phases:

    1. Size-only greedy: largest genes first (seed-shuffled tie-break among equal
       sizes), each assigned to whichever side is further below its 80/20 variant-count
       target. Gets close to the size target fast, but a single irrevocable largest-first
       pass can trap a handful of same-size capped genes together on one side (e.g. all
       6 genes at the 5,804 cap landing in val), skewing label balance badly.
    2. Local search: repeatedly flip single genes between sides, in a fixed
       (alphabetical) order, whenever it strictly reduces a combined cost of
       (size-target deviation) + label_weight*(benign-fraction-target deviation), until
       no flip improves it. This trades a few points of size precision to let smaller
       genes swap in/out and correct label balance -- e.g. label_weight=5.0 recovers
       val benign fraction to within 0.0003 of the overall 0.5374 while landing at
       19.997% (vs. 20.61% and a 45.4% benign flip from phase 1 alone).
    """
    target_train, target_val = 0.8 * total_n, 0.2 * total_n

    def cost(train_n, train_benign, val_n, val_benign):
        size_cost = abs(train_n / target_train - 1) + abs(val_n / target_val - 1)
        tb_frac = train_benign / train_n if train_n else overall_benign_frac
        vb_frac = val_benign / val_n if val_n else overall_benign_frac
        label_cost = abs(tb_frac - overall_benign_frac) + abs(vb_frac - overall_benign_frac)
        return size_cost + label_weight * label_cost

    # Phase 1: size-only greedy, seed-shuffled tie-break, largest genes first.
    rng = random.Random(seed)
    order = gene_stats.index.tolist()
    rng.shuffle(order)
    order.sort(key=lambda g: -gene_stats.loc[g, "n"])
    train_n = train_benign = val_n = val_benign = 0
    assignment: dict[str, str] = {}
    for gene in order:
        n = int(gene_stats.loc[gene, "n"])
        benign = int(gene_stats.loc[gene, "benign"])
        size_cost_train = abs((train_n + n) / target_train - 1) + abs(val_n / target_val - 1)
        size_cost_val = abs(train_n / target_train - 1) + abs((val_n + n) / target_val - 1)
        if size_cost_train <= size_cost_val:
            assignment[gene] = "train"
            train_n += n
            train_benign += benign
        else:
            assignment[gene] = "val"
            val_n += n
            val_benign += benign

    # Phase 2: deterministic local search (alphabetical scan order) to fix label balance.
    genes_alpha = sorted(gene_stats.index.tolist())
    cur_cost = cost(train_n, train_benign, val_n, val_benign)
    for _ in range(max_passes):
        improved = False
        for gene in genes_alpha:
            n = int(gene_stats.loc[gene, "n"])
            benign = int(gene_stats.loc[gene, "benign"])
            if assignment[gene] == "train":
                new_cost = cost(train_n - n, train_benign - benign, val_n + n, val_benign + benign)
                if new_cost < cur_cost - 1e-12:
                    assignment[gene] = "val"
                    train_n -= n; train_benign -= benign; val_n += n; val_benign += benign
                    cur_cost = new_cost; improved = True
            else:
                new_cost = cost(train_n + n, train_benign + benign, val_n - n, val_benign - benign)
                if new_cost < cur_cost - 1e-12:
                    assignment[gene] = "train"
                    train_n += n; train_benign += benign; val_n -= n; val_benign -= benign
                    cur_cost = new_cost; improved = True
        if not improved:
            break

    return assignment


def build_split_b(seed: int = 42) -> dict:
    df = pd.read_parquet(DATA)
    gene_stats = df.groupby("gene").agg(
        n=("label", "size"),
        benign=("label", lambda s: int((s == 0).sum())),
        pathogenic=("label", lambda s: int((s == 1).sum())),
    )
    total_n = len(df)
    overall_benign_frac = float((df.label == 0).mean())

    assignment = _assign_genes(gene_stats, total_n, overall_benign_frac, seed)
    df = df.assign(_side=df["gene"].map(assignment))
    train_df, val_df = df[df._side == "train"], df[df._side == "val"]

    train_ids = variant_ids(train_df).tolist()
    val_ids = variant_ids(val_df).tolist()
    val_genes = sorted(g for g, s in assignment.items() if s == "val")
    train_genes = sorted(g for g, s in assignment.items() if s == "train")

    extra = {
        "method": "gene-disjoint: whole genes assigned to train or val, greedy size+label-balance target, seed 42 tie-break",
        "n_genes_train": len(train_genes),
        "n_genes_val": len(val_genes),
        "genes_val": val_genes,
        "genes_train": train_genes,
        "train_pct": round(100 * len(train_df) / total_n, 2),
        "val_pct": round(100 * len(val_df) / total_n, 2),
        "train_label_balance": train_df.label.value_counts(normalize=True).round(4).to_dict(),
        "val_label_balance": val_df.label.value_counts(normalize=True).round(4).to_dict(),
        "train_label_counts": train_df.label.value_counts().to_dict(),
        "val_label_counts": val_df.label.value_counts().to_dict(),
        "largest_val_gene": gene_stats.loc[val_genes, "n"].idxmax() if val_genes else None,
        "largest_val_gene_share_of_val": round(
            float(gene_stats.loc[val_genes, "n"].max() / len(val_df)), 4
        ) if val_genes else None,
    }
    manifest = build_manifest(
        split_id="b", name="gene_disjoint",
        description="No gene has variants on both sides. Greedy-assigned to target ~80/20 by variant count "
                     "while nudging toward SPLIT-A's benign/pathogenic balance.",
        source_data=str(DATA.relative_to(ROOT)), train_ids=train_ids, val_ids=val_ids, seed=seed, extra=extra,
    )
    return manifest


def _load_variant_summary_dates(vs_path: Path, our_keys: set) -> dict:
    """Return {(chrom,pos,ref,alt): last_evaluated Timestamp} for GRCh38 rows matching our_keys.

    Reads in chunks and filters to Assembly==GRCh38 and PositionVCF/allele fields present
    (some older HGVS-only entries have "-" for these) to keep peak memory bounded --
    the raw file has several million rows across both assemblies.
    """
    cols = ["Assembly", "Chromosome", "PositionVCF", "ReferenceAlleleVCF", "AlternateAlleleVCF", "LastEvaluated"]
    dates: dict = {}
    for chunk in pd.read_csv(vs_path, sep="\t", usecols=cols, dtype=str, chunksize=500_000, low_memory=False):
        chunk = chunk[chunk.Assembly == "GRCh38"]
        chunk = chunk[(chunk.PositionVCF != "-") & (chunk.ReferenceAlleleVCF != "-") & (chunk.AlternateAlleleVCF != "-")]
        if chunk.empty:
            continue
        keys = list(zip(chunk.Chromosome, chunk.PositionVCF, chunk.ReferenceAlleleVCF, chunk.AlternateAlleleVCF))
        parsed = pd.to_datetime(chunk.LastEvaluated, format="%b %d, %Y", errors="coerce")
        for key, dt, pos in zip(keys, parsed, chunk.PositionVCF):
            k = (key[0], int(pos), key[2], key[3])
            if k in our_keys and pd.notna(dt):
                if k not in dates or dt > dates[k]:
                    dates[k] = dt
    return dates


def build_split_c(seed: int = 42, vs_path: Path = VARIANT_SUMMARY) -> dict:
    df = pd.read_parquet(DATA)
    our_keys = set(zip(df.chrom.astype(str), df.pos.astype(int), df.ref.astype(str), df.alt.astype(str)))

    dates = _load_variant_summary_dates(vs_path, our_keys)
    key_series = list(zip(df.chrom.astype(str), df.pos.astype(int), df.ref.astype(str), df.alt.astype(str)))
    df = df.assign(_last_evaluated=[dates.get(k) for k in key_series])

    n_total = len(df)
    n_dated = int(df._last_evaluated.notna().sum())
    n_undated = n_total - n_dated

    # Undated variants always land in train (see train_mask below), so the cutoff is
    # chosen against the overall 80/20 target -- not just the dated subset -- by taking
    # the ~20%-of-total most-recent dated variants as val.
    dated_sorted = df[df._last_evaluated.notna()].sort_values("_last_evaluated")
    target_val_n = min(round(0.2 * n_total), len(dated_sorted))
    cutoff_idx = max(0, len(dated_sorted) - target_val_n - 1)
    cutoff_date = dated_sorted._last_evaluated.iloc[cutoff_idx]

    train_mask = (df._last_evaluated <= cutoff_date) | df._last_evaluated.isna()
    val_mask = df._last_evaluated > cutoff_date
    train_df, val_df = df[train_mask], df[val_mask]

    train_ids = variant_ids(train_df).tolist()
    val_ids = variant_ids(val_df).tolist()

    extra = {
        "method": "temporal: ClinVar variant_summary_2025-06.txt.gz LastEvaluated (GRCh38), joined by chrom/PositionVCF/"
                  "ReferenceAlleleVCF/AlternateAlleleVCF; cutoff = 80th percentile of dated variants; "
                  "train = LastEvaluated <= cutoff (or undated, placed in train conservatively), val = LastEvaluated > cutoff",
        "date_field_used": "LastEvaluated (variant_summary.txt.gz)",
        "variant_summary_source": str(vs_path.relative_to(ROOT)) if vs_path.is_relative_to(ROOT) else str(vs_path),
        "cutoff_date": cutoff_date.strftime("%Y-%m-%d"),
        "n_dated": n_dated,
        "n_undated_placed_in_train": n_undated,
        "train_pct": round(100 * len(train_df) / n_total, 2),
        "val_pct": round(100 * len(val_df) / n_total, 2),
        "train_label_balance": train_df.label.value_counts(normalize=True).round(4).to_dict(),
        "val_label_balance": val_df.label.value_counts(normalize=True).round(4).to_dict(),
        "train_label_counts": train_df.label.value_counts().to_dict(),
        "val_label_counts": val_df.label.value_counts().to_dict(),
    }
    manifest = build_manifest(
        split_id="c", name="temporal",
        description="Train = variants last evaluated on/before the cutoff date (or with no resolvable date). "
                     "Val = variants last evaluated strictly after the cutoff date.",
        source_data=str(DATA.relative_to(ROOT)), train_ids=train_ids, val_ids=val_ids, seed=seed, extra=extra,
    )
    return manifest


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--split", choices=["a", "b", "c"], required=True)
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--variant-summary", default=str(VARIANT_SUMMARY), help="SPLIT-C only: path to variant_summary_YYYY-MM.txt.gz")
    args = ap.parse_args()

    if args.split == "a":
        manifest = build_split_a(args.seed)
    elif args.split == "b":
        manifest = build_split_b(args.seed)
    else:
        manifest = build_split_c(args.seed, Path(args.variant_summary))

    out = SPLITS_DIR / f"split_manifest_{args.split}.json"
    save_manifest(manifest, out)
    summary = {k: v for k, v in manifest.items() if k not in ("train_ids", "val_ids", "genes_train")}
    print(f"Wrote {out}")
    import json
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
