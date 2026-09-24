#!/usr/bin/env bash
set -euo pipefail
mkdir -p results/compute_gate
for m in cnn bilstm cnn_bilstm ensemble dnabert2 nt; do
  /usr/bin/time -f '%e' -o "results/compute_gate/${m}_100bp.seconds" \
    python -m src.training.run_experiment --data data/processed/variants_100bp.parquet --model "$m" --window 100 --fraction 0.05 --epochs 1 --batch-size 16 \
    > "results/compute_gate/${m}_100bp.jsonl"
done
python - <<'PY'
from pathlib import Path
import pandas as pd
rows=[]
for p in Path('results/compute_gate').glob('*.seconds'):
    sec=float(p.read_text()); model=p.stem.replace('_100bp',''); rows.append({'model':model,'window':100,'time_per_5pct_epoch_s':sec,'estimated_full_epoch_s':sec/0.05,'estimated_10_epoch_s':sec/0.05*10})
pd.DataFrame(rows).sort_values('estimated_10_epoch_s').to_csv('results/compute_gate/cost_table.csv',index=False)
print(pd.DataFrame(rows).sort_values('estimated_10_epoch_s').to_string(index=False))
PY
