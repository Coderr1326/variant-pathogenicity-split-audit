#!/usr/bin/env python
"""Summarize DNABERT-2 SPLIT-B seed runs and report which docs/seed_stopping_rule.md bucket they fall in.

Reads the existing seed-42 result plus every results/metrics/dnabert2_100bp_splitB_seed*_final.json.
Report-only: it never launches or queues anything, and it never modifies any file.
"""
from __future__ import annotations

import os

os.environ["CUDA_VISIBLE_DEVICES"] = ""  # CPU only

import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from src.evaluation.generate_results import checkpoint_best_epoch  # noqa: E402

METRICS_DIR = ROOT / "results" / "metrics"
PRED_DIR = ROOT / "results" / "predictions"

# Thresholds from docs/seed_stopping_rule.md (committed before any new run).
# "within +/-0.03 of each other" is implemented pairwise: max - min <= 0.03.
STABLE_RANGE, STABLE_MAX = 0.03, 0.75
INSTABILITY_RANGE, INSTABILITY_ANY = 0.10, 0.80
RUN_CAP, BASE_RUNS = 6, 4
COLLAPSE_LO, COLLAPSE_HI = 0.01, 0.99  # predicted-positive fraction ~0 or ~1
MAIN_WARMUP = 0.1


def load_run(final_path: Path, tag: str) -> dict:
    data = json.loads(final_path.read_text())
    m = data["metrics"]
    # The original seed-42 result predates the seed/warmup_ratio keys; it was trained with warmup 0.1
    # (its checkpoints carry a scheduler state, 84,470 steps ending at LR 0).
    seed = data.get("seed", 42)
    warmup = data.get("warmup_ratio", MAIN_WARMUP)
    pred = pd.read_parquet(PRED_DIR / f"dnabert2_100bp_splitB{tag}.parquet")
    pos_frac = float(pred["prediction"].astype(int).mean())
    epoch, _src = checkpoint_best_epoch("dnabert2", m["weighted_f1"], "_splitB" + tag)
    return {
        "seed": seed, "warmup": warmup, "accuracy": m["accuracy"], "weighted F1": m["weighted_f1"],
        "AUROC": m["roc_auc"], "path. recall": m["pathogenic_recall"],
        "best epoch": epoch if epoch is not None else "?", "pred-pos frac": pos_frac,
        "flag": "COLLAPSED" if (pos_frac < COLLAPSE_LO or pos_frac > COLLAPSE_HI) else "",
    }


def main() -> None:
    runs = []
    base = METRICS_DIR / "dnabert2_100bp_splitB_final.json"
    if base.exists():
        runs.append(load_run(base, ""))
    for path in sorted(METRICS_DIR.glob("dnabert2_100bp_splitB_seed*_final.json")):
        tag = path.name[len("dnabert2_100bp_splitB"):-len("_final.json")]
        runs.append(load_run(path, tag))
    if not runs:
        print("No DNABERT-2 SPLIT-B results found.")
        return

    df = pd.DataFrame(runs).sort_values(["warmup", "seed"], ascending=[False, True]).reset_index(drop=True)
    print(df.to_string(index=False, float_format=lambda v: f"{v:.4f}"))

    main_runs = df[df.warmup == MAIN_WARMUP]
    n = len(main_runs)
    print(f"\nAUROC over the warmup-{MAIN_WARMUP} runs only (n={n}, seeds {main_runs.seed.tolist()}):")
    if n == 0:
        return
    lo, hi, mean = main_runs.AUROC.min(), main_runs.AUROC.max(), main_runs.AUROC.mean()
    rng = hi - lo
    print(f"  min={lo:.4f}  max={hi:.4f}  mean={mean:.4f}  range={rng:.4f}")
    collapsed = main_runs[main_runs.flag == "COLLAPSED"].seed.tolist()
    if collapsed:
        print(f"  COLLAPSED warmup-{MAIN_WARMUP} runs (seeds): {collapsed}")

    print("\nRule bucket (docs/seed_stopping_rule.md):")
    if n < BASE_RUNS:
        print(f"  INCOMPLETE: {n} of {BASE_RUNS} required warmup-{MAIN_WARMUP} runs available; no bucket yet.")
    else:
        if rng > INSTABILITY_RANGE or hi >= INSTABILITY_ANY:
            bucket, action = "INSTABILITY", "rule satisfied; no additional seeds."
        elif rng <= STABLE_RANGE and hi <= STABLE_MAX:
            bucket, action = "STABLE EFFECT", "rule satisfied; no additional seeds."
        else:
            bucket = "IN BETWEEN"
            if n <= BASE_RUNS:
                action = f"rule says: add seeds 4 and 5 (hard cap of {RUN_CAP} total runs)."
            elif n < RUN_CAP:
                action = f"rule says: seeds 4 and 5 are in progress; {RUN_CAP - n} more run(s) allowed (hard cap {RUN_CAP})."
            else:
                action = f"hard cap of {RUN_CAP} runs reached; report as-is and state the cap in the writeup. No more seeds."
        print(f"  {bucket}  (range {rng:.4f}, max {hi:.4f}; stable needs range<={STABLE_RANGE} & max<={STABLE_MAX}; "
              f"instability is range>{INSTABILITY_RANGE} or any>={INSTABILITY_ANY})")
        print(f"  -> {action}")
        if n > RUN_CAP:
            print(f"  WARNING: {n} runs exceeds the hard cap of {RUN_CAP}.")

    sec = df[df.warmup != MAIN_WARMUP]
    if len(sec):
        print("\nSecondary check (not part of the bucket): runs with warmup != 0.1")
        for _, r in sec.iterrows():
            print(f"  seed {int(r.seed)} warmup {r.warmup}: AUROC {r.AUROC:.4f}, pred-pos frac {r['pred-pos frac']:.4f} {r.flag}")
    print("\n(report only -- this script never launches anything)")


if __name__ == "__main__":
    main()
