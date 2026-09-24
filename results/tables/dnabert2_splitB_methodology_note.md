# DNABERT-2 SPLIT-B methodology note

† DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule.

- The first SPLIT-B attempt (fixed LR 2e-5, no warmup) collapsed: accuracy 0.5377, ROC-AUC 0.4999, tp=0 (every validation example predicted benign). It was discarded.
- The reported SPLIT-B result is the rerun with linear LR warmup (`--warmup-ratio 0.1`, warmup over the first 10% of total steps to 2e-5, then linear decay; `scripts/train_transformer_local.py`).
- NT needed no warmup: NT SPLIT-A/B/C all used the fixed-LR schedule (`--warmup-ratio 0`), so NT rows carry no marker.
- SPLIT-A and SPLIT-C DNABERT-2 runs used the original fixed-LR schedule (no scheduler). The DNABERT-2 SPLIT-A/C vs SPLIT-B comparison therefore also differs in LR schedule, not only in split.
- Reported SPLIT-B numbers (best epoch 7): accuracy 0.6188, weighted F1 0.5901, ROC-AUC 0.6740, pathogenic recall 0.3425.
- `results/metrics/dnabert2_100bp_splitB_epoch_times.csv` is append-only and contains BOTH runs: the first 10 rows (2026-09-18) are the discarded collapsed run, the last 10 rows (2026-09-19) are the reported rerun.

Where marked (†): confusion_matrix_dnabert2_splitB.png, local_results_summary_splitB, cross_split_comparison, leakage_delta, cross_split_comparison_chart.png. The .csv files carry no marker (kept machine-clean); this note applies to them.
