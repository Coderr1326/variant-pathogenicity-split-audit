#!/usr/bin/env python
"""D3: allele-frequency (AF) matched re-scoring of already-trained models.

The models never receive AF. AF is joined locally from the raw ClinVar VCF that built the dataset
(data/raw/clinvar_20250630.vcf.gz, INFO AF_ESP / AF_EXAC / AF_TGP) by the variant key (CHROM, POS, REF, ALT); one AF per variant
= the MAX over the non-missing of the three fields; a variant with none of them is in the "missing" bucket (not the same as AF = 0).
No downloads, no network, no model inference, CPU only, add-only. Saved probabilities only (pred_loader): M3, NT, DNABERT-2 seeds
42/1/2/3 (SPLIT-B warmup 0.1; A/C single run, fixed LR). CNN/BiLSTM/CNN+BiLSTM/Ensemble have none and are skipped. M4 = AF-only baseline.

Matching procedure (per split test set, per resample r = 0..99): rng = default_rng(r); for each eligible bucket in BUCKET_LABELS order,
positives and negatives are the bucket's rows in manifest order; the majority class is downsampled to the minority count with
rng.choice(rows, k, replace=False) (positives drawn first, then negatives); the matched subset is the union over eligible buckets.
Eligible = at least MIN_PER_CLASS_BUCKET of BOTH classes. Two variants: all eligible buckets (incl. "missing") and AF-known only (excl. "missing").

    python scripts/run_d3_af_matched.py --before-sha <sha256 manifest> --before-ckpt <checkpoint size+mtime manifest> [--rebuild-af]
"""
from __future__ import annotations

import os

os.environ["CUDA_VISIBLE_DEVICES"] = ""

import argparse
import csv
import gzip
import hashlib
import io
import json
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
from matplotlib.patches import Patch
from scipy.stats import rankdata

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
import pred_loader as pl  # noqa: E402
import src.evaluation.generate_results as gr  # noqa: E402
from src.evaluation.generate_split_comparison import SPLIT_COLORS, WARMUP_MARK, WARMUP_NOTE  # noqa: E402
from src.training.run_experiment import metrics as full_metrics  # noqa: E402

# ---- parameters, all fixed before any result was looked at --------------------------------------------------------
BUCKET_LABELS = ["missing", "0", "(0,1e-4)", "[1e-4,1e-3)", "[1e-3,1e-2)", "[1e-2,5e-2)", ">=5e-2"]
CUTOFFS = (1e-4, 1e-3, 1e-2, 5e-2)       # bucket edges (see BUCKET_LABELS); AF == 0 is its own bucket; NaN -> "missing"
AF_FIELDS = ["AF_ESP", "AF_EXAC", "AF_TGP"]
MIN_PER_CLASS_BUCKET = 50                # a bucket is used for matching / within-bucket AUROC only with >= this many of BOTH classes
N_RESAMPLES = 100                        # resample seeds 0..99
CI_PCT = (2.5, 97.5)
ALPHA = 1.0                              # M4 smoothing toward the global train rate
VCF = ROOT / "data" / "raw" / "clinvar_20250630.vcf.gz"
CAVEAT = ("ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected "
          "and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models.")
SKIPPED_MODELS = ["CNN", "BiLSTM", "CNN+BiLSTM", "Ensemble"]
SKIP_NOTE = "Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities)."
OUT = ROOT / "results" / "d3_af_matched"
HANDOFF = ROOT / "results" / "handoff"
gr.TABLE_DIR = OUT
f4 = lambda v: "-" if v is None or pd.isna(v) else f"{v:.4f}"  # noqa: E731
RUNS = [("M3", "m3", None, "ABC"), ("NT", "nt", None, "ABC"), ("DNABERT-2 seed 42", "dnabert2", 42, "ABC"),
        ("DNABERT-2 seed 1", "dnabert2", 1, "B"), ("DNABERT-2 seed 2", "dnabert2", 2, "B"), ("DNABERT-2 seed 3", "dnabert2", 3, "B")]


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


def auroc(y, s) -> float:
    y = np.asarray(y)
    n1 = int(y.sum())
    n0 = len(y) - n1
    if n1 == 0 or n0 == 0:
        return float("nan")
    r = rankdata(s)
    return float((r[y == 1].sum() - n1 * (n1 + 1) / 2) / (n1 * n0))


def bucket_index(af: np.ndarray) -> np.ndarray:
    ok = ~np.isnan(af)
    b = np.select([af == 0, af < CUTOFFS[0], af < CUTOFFS[1], af < CUTOFFS[2], af < CUTOFFS[3]], [1, 2, 3, 4, 5], default=6)
    return np.where(ok, b, 0)


# ---- AF join from the raw ClinVar VCF -----------------------------------------------------------------------------
def build_af_table(rebuild: bool) -> pd.DataFrame:
    cache = OUT / "af_table.csv"
    data = pl._data()
    if cache.exists() and not rebuild:
        t = pd.read_csv(cache, dtype={"vid": str})
        if len(t) == len(data) and (t.vid.to_numpy() == data.vid.to_numpy()).all():
            return t
    keys = set(data.vid)
    found: dict[str, list] = {}
    with gzip.open(VCF, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            p = line.split("\t", 8)
            if len(p) < 8:
                continue
            chrom, pos, _id, ref, alts, _q, _f, info = p[:8]
            alt_list = alts.split(",")
            hit = [i for i, a in enumerate(alt_list) if f"{chrom}_{pos}_{ref}_{a}" in keys]
            if not hit:
                continue
            vals = {}
            for item in info.split(";"):
                k, _, v = item.partition("=")
                if k in AF_FIELDS:
                    vals[k] = v
            for i in hit:
                row = []
                for f in AF_FIELDS:
                    v = vals.get(f)
                    parts = v.split(",") if v is not None else []
                    pick = parts[i] if len(parts) == len(alt_list) else (parts[0] if len(alt_list) == 1 and parts else None)
                    try:
                        row.append(float(pick) if pick not in (None, "", ".") else np.nan)
                    except ValueError:
                        row.append(np.nan)
                found[f"{chrom}_{pos}_{ref}_{alt_list[i]}"] = row
    missing_keys = keys - set(found)
    if missing_keys:
        raise AssertionError(f"{len(missing_keys)} dataset variants were not found in the VCF by (CHROM,POS,REF,ALT)")
    arr = np.array([found[v] for v in data.vid], dtype=float)
    t = pd.DataFrame({"vid": data.vid.to_numpy(), "af_esp": arr[:, 0], "af_exac": arr[:, 1], "af_tgp": arr[:, 2]})
    t["af_max"] = np.where(np.isnan(arr).all(axis=1), np.nan, np.nanmax(np.where(np.isnan(arr), -1.0, arr), axis=1))
    OUT.mkdir(parents=True, exist_ok=True)
    t.to_csv(cache, index=False)
    return t


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
    ap.add_argument("--rebuild-af", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    HANDOFF.mkdir(parents=True, exist_ok=True)

    # ---- Steps 1/3: AF per variant, coverage, buckets ---------------------------------------------------------------
    data = pl._data()
    af = build_af_table(args.rebuild_af)
    af_by_vid = af.set_index("vid")
    bucket_by_vid = pd.Series(bucket_index(af.af_max.to_numpy()), index=af.vid.to_numpy())
    label_by_vid = pd.Series(data.label.to_numpy(), index=data.vid.to_numpy())
    manifests = {s: pl._manifest_split(s) for s in "ABC"}
    sets = {}  # name -> vid array (manifest order)
    for s in "ABC":
        _, train, test = manifests[s]
        sets[f"train_{s}"], sets[f"test_{s}"] = train.vid.to_numpy(), test.vid.to_numpy()

    cov_rows = []
    for name, vids in [("all_variants", data.vid.to_numpy())] + list(sets.items()):
        lab = label_by_vid.loc[vids].to_numpy()
        has = ~np.isnan(af_by_vid.loc[vids, "af_max"].to_numpy())
        for lname, m in (("all", np.ones(len(vids), bool)), ("pathogenic", lab == 1), ("benign", lab == 0)):
            cov_rows.append([name, lname, int(m.sum()), int((has & m).sum()), float((has & m).sum() / m.sum())])
    src_rows = []
    for f, col in zip(AF_FIELDS, ["af_esp", "af_exac", "af_tgp"]):
        has = ~af[col].isna().to_numpy()
        src_rows.append([f, int(has.sum()), float(has.mean()), int((has & (data.label.to_numpy() == 1)).sum()), int((has & (data.label.to_numpy() == 0)).sum())])
    cov = pd.DataFrame(cov_rows, columns=["set", "label", "n", "n_with_af", "frac_with_af"])
    cov.to_csv(OUT / "af_coverage.csv", index=False)
    pd.DataFrame(src_rows, columns=["source_field", "n_with_value", "frac_with_value", "n_pathogenic_with_value", "n_benign_with_value"]).to_csv(OUT / "af_coverage_by_source.csv", index=False)
    gr.write_table(pd.DataFrame({"Set": cov.set, "Label": cov.label, "n": cov.n, "With AF": cov.n_with_af, "Fraction with AF": cov.frac_with_af.map(f4)}),
                   "view_af_coverage", "AF Coverage (non-missing AF_ESP / AF_EXAC / AF_TGP) by Set and Label",
                   footnote="AF = max of the non-missing of AF_ESP, AF_EXAC, AF_TGP from the ClinVar VCF INFO; 'missing' means none of the three is present (not AF = 0). " + CAVEAT)

    # bucket x label counts (train sets and test sets)
    cnt_rows = []
    for name, vids in sets.items():
        b, lab = bucket_by_vid.loc[vids].to_numpy(), label_by_vid.loc[vids].to_numpy()
        for bi, bl in enumerate(BUCKET_LABELS):
            cnt_rows.append([name, bl, int(((b == bi) & (lab == 1)).sum()), int(((b == bi) & (lab == 0)).sum())])
    cnt = pd.DataFrame(cnt_rows, columns=["set", "bucket", "n_pathogenic", "n_benign"])
    cnt.to_csv(OUT / "bucket_label_counts.csv", index=False)
    view = pd.DataFrame({"Bucket": BUCKET_LABELS})
    for name in sets:
        t = cnt[cnt.set == name].set_index("bucket").loc[BUCKET_LABELS]
        view[name.replace("_", " ").replace("train", "Train").replace("test", "Test")] = [f"{p:,} / {b:,}" for p, b in zip(t.n_pathogenic, t.n_benign)]
    gr.write_table(view, "view_bucket_label_counts", "AF Bucket x Label Counts (pathogenic / benign)",
                   footnote="Buckets are fixed: missing; AF = 0; (0,1e-4); [1e-4,1e-3); [1e-3,1e-2); [1e-2,5e-2); >=5e-2. 'missing' is its own category, never dropped. " + CAVEAT)

    # ---- Step 4: M4 (AF-only baseline) -------------------------------------------------------------------------------
    frames = {}
    for s in "ABC":
        _, train, test = manifests[s]
        fr = pl.test_frame(s, "manifest").copy()
        fr["bucket"] = bucket_by_vid.loc[fr.vid].to_numpy()
        frames[s] = fr
    m4_scores, m4_rows = {}, []
    for s in "ABC":
        tr_vids = sets[f"train_{s}"]
        tb, tl = bucket_by_vid.loc[tr_vids].to_numpy(), label_by_vid.loc[tr_vids].to_numpy()
        g = float(tl.mean())
        rate = {bi: (float(tl[tb == bi].sum()) + ALPHA * g) / (float((tb == bi).sum()) + ALPHA) for bi in range(len(BUCKET_LABELS))}
        fr = frames[s]
        score = fr.bucket.map(rate).to_numpy(dtype=float)
        pred = (score >= 0.5).astype(int)
        met = full_metrics(fr.label.to_numpy(dtype=int).tolist(), pred.tolist(), score.tolist())
        d = OUT / "m4_af_only" / f"split_{s}"
        d.mkdir(parents=True, exist_ok=True)
        (d / "final_result.json").write_text(json.dumps({"model": "m4_af_only", "window": 100, "metrics": met}, indent=2))
        pd.DataFrame({"label": fr.label.to_numpy(dtype=int), "prediction": pred, "probability_pathogenic": score}).to_parquet(d / "predictions.parquet", index=False)
        m4_scores[s] = score
        m3 = json.loads((ROOT / "results" / "m3_gene_only" / f"split_{s}" / "final_result.json").read_text())["metrics"]
        m4_rows.append([s, met["accuracy"], met["weighted_f1"], met["roc_auc"], met["pathogenic_recall"], m3["roc_auc"], g, *[rate[i] for i in range(len(BUCKET_LABELS))]])
    m4 = pd.DataFrame(m4_rows, columns=["split", "accuracy", "weighted_f1", "roc_auc", "pathogenic_recall", "m3_roc_auc_reference", "global_train_rate"] + [f"rate_{b}" for b in BUCKET_LABELS])
    m4.to_csv(OUT / "m4_af_only_metrics.csv", index=False)
    gr.write_table(pd.DataFrame({"Split": m4.split, "Accuracy": m4.accuracy.map(f4), "Weighted F1": m4.weighted_f1.map(f4), "ROC-AUC": m4.roc_auc.map(f4),
                                 "Pathogenic recall": m4.pathogenic_recall.map(f4), "M3 ROC-AUC (ref.)": m4.m3_roc_auc_reference.map(f4)}),
                   "view_m4_af_only", "M4: AF-Only Baseline (pathogenic rate per AF bucket, train only, alpha = 1)",
                   footnote="Threshold 0.5. M4 uses no sequence and no gene. " + CAVEAT)

    # ---- Step 5: matched subsets ------------------------------------------------------------------------------------
    elig_rows, resamples = [], {}
    for s in "ABC":
        fr = frames[s]
        y, b = fr.label.to_numpy(dtype=int), fr.bucket.to_numpy()
        per = {}
        for bi, bl in enumerate(BUCKET_LABELS):
            pos, neg = np.flatnonzero((b == bi) & (y == 1)), np.flatnonzero((b == bi) & (y == 0))
            ok = min(len(pos), len(neg)) >= MIN_PER_CLASS_BUCKET
            per[bi] = (pos, neg, ok)
            elig_rows.append([s, bl, len(pos), len(neg), ok, "" if ok else f"smaller class has {min(len(pos), len(neg))} (< {MIN_PER_CLASS_BUCKET})"])
        for variant in ("all_eligible", "af_known_only"):
            buckets = [bi for bi in range(len(BUCKET_LABELS)) if per[bi][2] and (variant == "all_eligible" or bi != 0)]
            sets_r = []
            for r in range(N_RESAMPLES):
                rng = np.random.default_rng(r)
                parts = []
                for bi in buckets:
                    pos, neg, _ = per[bi]
                    k = min(len(pos), len(neg))
                    parts.append(pos if len(pos) == k else rng.choice(pos, k, replace=False))
                    parts.append(neg if len(neg) == k else rng.choice(neg, k, replace=False))
                sets_r.append(np.sort(np.concatenate(parts)) if parts else np.array([], dtype=int))
            resamples[(s, variant)] = (buckets, sets_r)
    elig = pd.DataFrame(elig_rows, columns=["split", "bucket", "n_pathogenic", "n_benign", "eligible", "reason_if_skipped"])
    elig.to_csv(OUT / "matching_eligibility.csv", index=False)
    gr.write_table(pd.DataFrame({"Split": elig.split, "Bucket": elig.bucket, "Pathogenic": elig.n_pathogenic, "Benign": elig.n_benign,
                                 "Used": elig.eligible.map({True: "yes", False: "skipped"}), "Reason": elig.reason_if_skipped.replace("", "-")}),
                   "view_matching_eligibility", f"Bucket Eligibility for Matching (>= {MIN_PER_CLASS_BUCKET} of both classes)",
                   footnote=CAVEAT)

    # run data
    run_data = {}
    for name, model, seed, splits in RUNS:
        for s in splits:
            d = pl.load_predictions(model, s, seed)
            fr_run = pl.test_frame(s, pl._order_kind(model, s))
            pos_in_run = pd.Series(np.arange(len(fr_run)), index=fr_run.vid.to_numpy())
            run_data[(name, s)] = (d.label.to_numpy(), d.score.to_numpy(), pos_in_run.loc[frames[s].vid.to_numpy()].to_numpy())
    for s in "ABC":
        y = frames[s].label.to_numpy(dtype=int)
        run_data[("M4", s)] = (y, m4_scores[s], np.arange(len(y)))
    run_names = ["M3", "M4", "NT", "DNABERT-2 seed 42", "DNABERT-2 seed 1", "DNABERT-2 seed 2", "DNABERT-2 seed 3"]
    res_rows, per_resample, wb_rows = [], [], []
    for s in "ABC":
        for name in run_names:
            if (name, s) not in run_data:
                continue
            y, sc, mp = run_data[(name, s)]
            full = auroc(y, sc)
            for variant in ("all_eligible", "af_known_only"):
                buckets, sets_r = resamples[(s, variant)]
                if not buckets:
                    res_rows.append([name, s, variant, full, None, None, None, None, 0])
                    continue
                vals = np.array([auroc(y[mp[idx]], sc[mp[idx]]) for idx in sets_r])
                for r, v in enumerate(vals):
                    per_resample.append([name, s, variant, r, len(sets_r[r]), float(v)])
                lo, hi = np.percentile(vals, CI_PCT)
                res_rows.append([name, s, variant, full, float(vals.mean()), float(lo), float(hi), float(vals.mean() - full), len(sets_r[0])])
            b = frames[s].bucket.to_numpy()
            yy = frames[s].label.to_numpy(dtype=int)
            for bi, bl in enumerate(BUCKET_LABELS):
                ok = elig[(elig.split == s) & (elig.bucket == bl)].eligible.iloc[0]
                if ok:
                    rows_b = np.flatnonzero(b == bi)
                    wb_rows.append([name, s, bl, len(rows_b), int(yy[rows_b].sum()), int(len(rows_b) - yy[rows_b].sum()), auroc(y[mp[rows_b]], sc[mp[rows_b]])])
    res = pd.DataFrame(res_rows, columns=["run", "split", "variant", "full_auroc", "matched_mean", "matched_p2_5", "matched_p97_5", "delta_matched_minus_full", "matched_n"])
    d2 = res[res.run.str.startswith("DNABERT-2 seed") & (res.split == "B")]
    agg = []
    for variant in ("all_eligible", "af_known_only"):
        t = d2[d2.variant == variant]
        for stat in ("mean", "min", "max"):
            agg.append([f"DNABERT-2 {stat} over seeds", "B", variant] + [None if t[c].isna().any() else float(getattr(t[c], stat)()) for c in ["full_auroc", "matched_mean", "matched_p2_5", "matched_p97_5", "delta_matched_minus_full"]] + [int(t.matched_n.iloc[0])])
    res_all = pd.concat([res, pd.DataFrame(agg, columns=res.columns)], ignore_index=True)
    res_all.to_csv(OUT / "matched_auroc.csv", index=False)
    pd.DataFrame(per_resample, columns=["run", "split", "variant", "resample", "matched_n", "auroc"]).to_csv(OUT / "matched_resample_aurocs.csv", index=False)
    wb = pd.DataFrame(wb_rows, columns=["run", "split", "bucket", "n", "n_pathogenic", "n_benign", "auroc"])
    wb.to_csv(OUT / "within_bucket_auroc.csv", index=False)

    def lab(r: str) -> str:
        return r.replace(" over seeds", "").replace("DNABERT-2", "DNABERT-2" + WARMUP_MARK) if r.startswith("DNABERT-2") else r

    for s in "ABC":
        t = res_all[(res_all.split == s) & (res_all.variant == "all_eligible")]
        gr.write_table(pd.DataFrame({"Run": t.run.map(lab), "Full test AUROC": t.full_auroc.map(f4), "Matched mean": t.matched_mean.map(f4), "Matched 2.5%": t.matched_p2_5.map(f4),
                                     "Matched 97.5%": t.matched_p97_5.map(f4), "Delta (matched - full)": t.delta_matched_minus_full.map(f4), "Matched n": t.matched_n}),
                       f"view_matched_auroc_split{s}", f"AF-Matched AUROC vs. Full Test Set: SPLIT-{s}",
                       footnote=f"Matched subset = per-bucket class downsampling over all eligible buckets incl. 'missing' ({N_RESAMPLES} resamples, seeds 0-{N_RESAMPLES - 1}); "
                                f"mean and {CI_PCT[0]}-{CI_PCT[1]} percentile of the resample AUROCs. DNABERT-2 SPLIT-A/C: single run (seed 42, fixed LR); SPLIT-B mean/min/max are over seeds 42/1/2/3. "
                                f"{WARMUP_NOTE} {CAVEAT} {SKIP_NOTE}")
        w = wb[wb.split == s]
        wide = w.pivot(index="run", columns="bucket", values="auroc").reindex([r for r in run_names if r in set(w.run)])
        cols = [b for b in BUCKET_LABELS if b in wide.columns]
        vv = pd.DataFrame({"Run": [lab(r) for r in wide.index]})
        for b_ in cols:
            vv[b_] = wide[b_].map(f4).to_numpy()
        gr.write_table(vv, f"view_within_bucket_auroc_split{s}", f"Plain Within-Bucket AUROC (no resampling): SPLIT-{s}",
                       footnote=f"Only buckets with >= {MIN_PER_CLASS_BUCKET} of both classes are shown. {WARMUP_NOTE} {CAVEAT} {SKIP_NOTE}")
    t = res_all[(res_all.variant == "af_known_only") & res_all.matched_mean.notna()]
    if len(t):
        gr.write_table(pd.DataFrame({"Split": t.split, "Run": t.run.map(lab), "Full test AUROC": t.full_auroc.map(f4), "Matched mean": t.matched_mean.map(f4), "Matched 2.5%": t.matched_p2_5.map(f4),
                                     "Matched 97.5%": t.matched_p97_5.map(f4), "Delta": t.delta_matched_minus_full.map(f4), "Matched n": t.matched_n}),
                       "view_matched_afknown_only", "AF-Matched AUROC, AF-Known Buckets Only (excluding 'missing')",
                       footnote=f"Same procedure as the main matched table but the 'missing' bucket is left out. {WARMUP_NOTE} {CAVEAT}")

    figure(res_all)

    # ---- export -----------------------------------------------------------------------------------------------------
    params = [["af_source", "data/raw/clinvar_20250630.vcf.gz INFO fields AF_ESP, AF_EXAC, AF_TGP (the release that built the dataset; no download)"],
              ["variant_key", "CHROM, POS, REF, ALT (dataset vid = chrom_pos_ref_alt); all 168,927 dataset variants found in the VCF"],
              ["af_combination", "max over the non-missing of AF_ESP, AF_EXAC, AF_TGP; none present = missing bucket (not AF = 0)"],
              ["bucket_cutoffs", "missing; ==0; (0,1e-4); [1e-4,1e-3); [1e-3,1e-2); [1e-2,5e-2); >=5e-2"], ["min_per_class_bucket", MIN_PER_CLASS_BUCKET],
              ["n_resamples", N_RESAMPLES], ["resample_seeds", "0..99, default_rng(r) per resample"], ["ci_percentiles", "2.5, 97.5 of the resample AUROCs"], ["m4_alpha", ALPHA],
              ["matching", "per bucket the majority class is downsampled to the minority count; union over eligible buckets; variants: all_eligible (incl. missing), af_known_only"],
              ["models_scored", "M3; M4; NT; DNABERT-2 seeds 42/1/2/3 (SPLIT-B warmup 0.1; SPLIT-A/C single run seed 42, fixed LR)"], ["python", platform.python_version()],
              ["generated", datetime.now().astimezone().strftime("%Y-%m-%d %H:%M:%S %z")]]
    pd.DataFrame(params, columns=["key", "value"]).to_csv(OUT / "parameters.csv", index=False)
    skipped = [[m, "model", "all", "no saved probabilities (weights and history only); predictions were not regenerated"] for m in SKIPPED_MODELS]
    skipped += [["DNABERT-2 seed 1 warmup 0", "run", "B", "collapsed run; never loaded"]]
    skipped += [["DNABERT-2 seeds 1, 2, 3", "run", "A;C", "no seed runs exist for SPLIT-A/C"]]
    skipped += [[r[1], "bucket", r[0], r[5]] for r in elig_rows if not r[4]]
    integ = integrity(args.before_sha, args.before_ckpt)
    secs = [block("1. af_definition_and_parameters", "Where AF comes from, how it is combined, and every fixed parameter.", ["key", "value"], params),
            block("2. af_coverage", "Fraction of variants with non-missing AF, by set and label. train_X/test_X = SPLIT-X manifest train/test set; all_variants = full 168,927-variant dataset.",
                  ["set", "label", "n", "n_with_af", "frac_with_af"], cov.values.tolist()),
            block("3. af_coverage_by_source", "Per-source coverage over all variants.", ["source_field", "n_with_value", "frac_with_value", "n_pathogenic_with_value", "n_benign_with_value"], src_rows),
            block("4. bucket_label_counts", "Pathogenic/benign counts per AF bucket for each split's train set and test set. 'missing' is its own category.",
                  ["set", "bucket", "n_pathogenic", "n_benign"], cnt.values.tolist()),
            block("5. m4_af_only", "M4 = AF-only baseline (pathogenic rate per bucket, train only, alpha=1, threshold 0.5). rate_* = fitted smoothed rates.", list(m4.columns), m4.values.tolist()),
            block("6. matching_eligibility", f"Buckets used for matching / within-bucket AUROC (>= {MIN_PER_CLASS_BUCKET} of both classes) per split test set.", list(elig.columns), elig.values.tolist()),
            block("7. matched_auroc", "Full-test AUROC vs. AF-matched AUROC (mean and 2.5-97.5 percentile over 100 resamples), delta = matched mean - full, matched_n = subset size. "
                  "variant all_eligible includes the 'missing' bucket; af_known_only excludes it. DNABERT-2 mean/min/max are over the four SPLIT-B seeds.", list(res_all.columns), res_all.values.tolist()),
            block("8. within_bucket_auroc", "Plain AUROC inside each eligible bucket (no resampling).", list(wb.columns), wb.values.tolist()),
            block("9. skipped", "Models, runs and buckets not used, with reasons.", ["item", "kind", "scope", "reason"], skipped),
            block("10. integrity_check", "sha256 of every file under results/, scripts/, docs/ before vs. after (results/checkpoints: size + mtime only). Added files are this task's outputs.", ["check", "value"], integ)]
    header = ("# D3 AF-matched export\n\n" f"Generated {datetime.now().astimezone().strftime('%Y-%m-%d %H:%M:%S %z')}. Numbers only; no mechanism is claimed. Floats have 6 decimals; blank = not defined.\n\n"
              f"Caveat: {CAVEAT}\n\nDNABERT-2 SPLIT-B results (seeds 42/1/2/3) used LR warmup 0.1 after the initial fixed-LR SPLIT-B attempt collapsed; DNABERT-2 SPLIT-A and SPLIT-C used the original fixed-LR schedule.\n\n")
    (HANDOFF / "d3_export.md").write_text(header + "\n".join(secs))

    g = lambda r, s, v="all_eligible": res_all[(res_all.run == r) & (res_all.split == s) & (res_all.variant == v)].iloc[0]  # noqa: E731
    print(f"report: {(HANDOFF / 'd3_export.md').relative_to(ROOT)}  (tables/figure: {OUT.relative_to(ROOT)}/)")
    c = cov[(cov.set == "all_variants")].set_index("label").frac_with_af
    print(f"1 AF coverage (non-missing): all {c['all']:.3f}, pathogenic {c['pathogenic']:.3f}, benign {c['benign']:.3f}; source ClinVar VCF AF_ESP/AF_EXAC/AF_TGP (max), no download")
    print("2 M4 AUROC A/B/C: " + " / ".join(f"{x:.4f}" for x in m4.roc_auc) + " | M4 accuracy: " + " / ".join(f"{x:.4f}" for x in m4.accuracy))
    for s in "ABC":
        used = [r[1] for r in elig_rows if r[0] == s and r[4]]
        print(f"3 SPLIT-{s}: buckets used {len(used)} of 7 ({', '.join(used)}); matched n = {g('M3', s).matched_n}")
    for r in ("NT", "DNABERT-2 seed 42", "M3", "M4"):
        print(f"4 {r}: full->matched(mean) A {g(r, 'A').full_auroc:.4f}->{g(r, 'A').matched_mean:.4f}, B {g(r, 'B').full_auroc:.4f}->{g(r, 'B').matched_mean:.4f}, C {g(r, 'C').full_auroc:.4f}->{g(r, 'C').matched_mean:.4f}")
    print(f"5 integrity: {[x for x in integ if x[0] == 'result'][0][1]}")


def figure(res: pd.DataFrame) -> None:
    models = [("M3 (gene-only)", "M3"), ("M4 (AF-only)", "M4"), ("NT", "NT"), (f"DNABERT-2{WARMUP_MARK}", "DNABERT-2 seed 42")]
    fig, axes = plt.subplots(1, 3, figsize=(17, 6.4), sharey=True)
    series = [("Full test set", "#6baed6", None), ("Matched (all eligible buckets)", "#08519c", None), ("Matched (AF-known buckets only)", "#fd8d3c", "//")]
    x, width = np.arange(len(models)), 0.26
    for ax, s in zip(axes, "ABC"):
        for i, (lab_, color, hatch) in enumerate(series):
            for j, (_, run) in enumerate(models):
                row = res[(res.run == run) & (res.split == s) & (res.variant == ("af_known_only" if i == 2 else "all_eligible"))]
                if row.empty:
                    continue
                r = row.iloc[0]
                v = r.full_auroc if i == 0 else r.matched_mean
                if pd.isna(v):
                    continue
                xp = x[j] + (i - 1) * width
                ax.bar(xp, v, width, color=color, hatch=hatch, edgecolor="black" if hatch else "none", linewidth=0.7)
                top = v
                if i > 0 and not pd.isna(r.matched_p2_5):
                    ax.errorbar(xp, v, yerr=[[max(0, v - r.matched_p2_5)], [max(0, r.matched_p97_5 - v)]], fmt="none", ecolor="black", capsize=2.5, elinewidth=1)
                    top = r.matched_p97_5
                ax.text(xp, top + 0.008, f"{v:.3f}", ha="center", va="bottom", fontsize=6.5, rotation=90)
        ax.axhline(0.5, color="dimgray", linestyle="--", linewidth=1)
        ax.set_xticks(x, [m[0] for m in models], rotation=15, ha="right")
        ax.set_title(f"SPLIT-{s}", fontsize=12, weight="bold", color=SPLIT_COLORS[s])
        ax.set_ylim(0, 1.12)
        ax.grid(axis="y", alpha=0.25)
    axes[0].set_ylabel("ROC-AUC")
    fig.legend(handles=[Patch(facecolor=c, hatch=h, edgecolor="black" if h else "none", label=l) for l, c, h in series], loc="upper center", bbox_to_anchor=(0.5, 0.935), ncol=3, fontsize=9.5, frameon=False)
    fig.suptitle("Full vs. AF-Matched ROC-AUC per Model and Split", fontsize=14, weight="bold")
    d2 = res[(res.run.str.startswith("DNABERT-2 seed")) & (res.split == "B") & (res.variant == "all_eligible")]
    foot = (f"{CAVEAT}\n{WARMUP_NOTE} DNABERT-2 SPLIT-B bar = seed 42 (matched-AUROC mean over seeds 42/1/2/3: {d2.matched_mean.mean():.3f}, range {d2.matched_mean.min():.3f}-{d2.matched_mean.max():.3f}). "
            f"Error bars: {CI_PCT[0]}-{CI_PCT[1]} percentile over {N_RESAMPLES} resamples. {SKIP_NOTE}")
    fig.tight_layout(rect=(0, 0.13, 1, 0.89))
    fig.text(0.5, 0.06, "\n".join(textwrap.fill(p, 200) for p in foot.split("\n")), ha="center", va="center", fontsize=8.5, style="italic", color="dimgray")
    fig.savefig(OUT / "d3_grouped_bar_full_vs_matched.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    main()
