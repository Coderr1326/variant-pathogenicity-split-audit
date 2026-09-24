#!/usr/bin/env bash
set -euo pipefail
for w in 30 50 100 200; do
  python -m src.data.build_dataset --window "$w"
done
for w in 100 30 50 200; do
  for m in cnn bilstm cnn_bilstm ensemble dnabert2 nt; do
    python -m src.training.run_experiment --data "data/processed/variants_${w}bp.parquet" --model "$m" --window "$w"
  done
done
python -m src.evaluation.compare_results
python -m src.visualization.figures
