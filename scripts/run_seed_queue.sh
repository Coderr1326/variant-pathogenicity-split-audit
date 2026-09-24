#!/usr/bin/env bash
# Sequential, resumable seed queue for DNABERT-2 SPLIT-B (see docs/seed_stopping_rule.md).
# One GPU job at a time. Stops with a non-zero exit at the first failed run.
# Does NOT queue seeds 4-5: the stopping rule decides that, after these four runs are reported.
#
#   nohup bash scripts/run_seed_queue.sh > logs/seeds/queue.log 2>&1 &
#
# A run is skipped when its final result JSON already exists. An interrupted run resumes from its own
# latest epoch checkpoint (train_transformer_local.py auto-resume, scheduler state included). Each run
# writes only to its own '_splitB_seed<N>...' paths, never to the seed-42 outputs. --restart is passed
# only if QUEUE_RESTART=1 is set explicitly (for a stale checkpoint dir you want to discard).
set -euo pipefail
cd "$(dirname "$0")/.."

PY="${PYTHON:-.venv/bin/python}"
LOG_DIR="logs/seeds"
mkdir -p "$LOG_DIR"

# "seed:warmup_ratio:tag" -- tag must contain seed<N> (enforced by the training script)
RUNS=(
  "1:0.1:_seed1"
  "2:0.1:_seed2"
  "3:0.1:_seed3"
  "1:0:_seed1_warmup0"
)

ts() { date '+%F %T'; }

for spec in "${RUNS[@]}"; do
  IFS=: read -r seed warmup tag <<< "$spec"
  name="dnabert2_splitB${tag}"
  final="results/metrics/dnabert2_100bp_splitB${tag}_final.json"
  log="$LOG_DIR/${name}.log"
  if [[ -f "$final" ]]; then
    echo "[$(ts)] SKIP  $name (final result exists: $final)"
    continue
  fi
  restart=()
  [[ "${QUEUE_RESTART:-0}" == "1" ]] && restart=(--restart)
  echo "[$(ts)] START $name (seed=$seed warmup_ratio=$warmup) log=$log"
  if PYTHONUNBUFFERED=1 "$PY" scripts/train_transformer_local.py \
        --model dnabert2 --data data/processed/variants_100bp.parquet \
        --split-manifest data/splits/split_manifest_b.json \
        --epochs 10 --batch-size 16 --warmup-ratio "$warmup" \
        --seed "$seed" --tag "$tag" "${restart[@]}" > "$log" 2>&1; then
    echo "[$(ts)] END   $name OK"
  else
    rc=$?
    echo "[$(ts)] END   $name FAILED (exit $rc) -- queue stopped; see $log" >&2
    exit "$rc"
  fi
done
echo "[$(ts)] QUEUE COMPLETE (seeds 4-5 are NOT auto-launched; apply docs/seed_stopping_rule.md first)"
