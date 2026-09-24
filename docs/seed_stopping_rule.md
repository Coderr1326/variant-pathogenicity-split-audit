# Seed stopping rule

Question: is DNABERT-2's SPLIT-B AUROC of 0.674 a stable effect or optimization fragility?

Runs: DNABERT-2, SPLIT-B, warmup 0.1, seeds 42 (existing), 1, 2, 3 => 4 runs.

Stable effect: all 4 AUROCs within +/-0.03 of each other AND none above 0.75.

Instability: range > 0.10 OR any seed >= 0.80.

In between: add seeds 4, 5 (hard cap of 6 total runs), then report whatever results and state the cap in the writeup. No stopping early on a preferred outcome, no adding seeds beyond the cap.

Report: min / max / mean / range only, no significance tests (n is too small).

Secondary check: seed 1 with --warmup-ratio 0 on SPLIT-B, to see whether the original collapse was specific to seed 42's batch order.

Thresholds were fixed before seeing any new results.
