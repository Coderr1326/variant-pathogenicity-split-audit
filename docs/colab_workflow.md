# Colab workflow

## Local preparation

1. Run the local data pipeline through Part B. The current canonical file is `data/processed/variants_100bp.parquet`.
2. Upload the repository (including `src/`) and the processed parquet file to a Drive folder, for example `MyDrive/aicaps_reproduction/`.
3. Upload `results/compute_gate/cost_table.csv` if you want the notebook to use updated local timing values.

## Colab execution

1. Open `notebooks/train_transformers_colab.ipynb` in Colab and enable a T4 or A100 runtime.
2. Set `DRIVE_ROOT` to the Drive folder containing `data/processed/` and `results/`.
3. Set `MODEL` to `dnabert2` or `nt`; keep training batch size 16.
4. Run setup, data, model, and checkpoint cells in order. The notebook imports `src.training.run_experiment.PairSet`, `metrics`, `seed`, and the shared transformer loader, so tokenization and splitting are the same as local execution.
5. Check the Drive-space warning before training. Each epoch writes `results/checkpoints/<model>/epoch_NN.pt`; epoch wall times append to `results/metrics/<model>_100bp_epoch_times.csv`.

## Multi-session behavior

The default session budget is five hours. After an epoch checkpoint is saved, the notebook stops with:

`Stopping early at epoch N/10 to stay within session budget — resume in a new session.`

This is a clean, resumable stop. A Colab disconnect may not print that message; inspect the newest `epoch_NN.pt` and the last row of the epoch-time CSV. Reconnect by remounting Drive, rerunning setup/import/model cells, setting `AICAPS_RESUME_CHECKPOINT` to the newest checkpoint path, and rerunning the resume/training cell. NT is expected to require multiple sessions on free-tier Colab.

## Returning results locally

Download or sync the Drive `results/predictions/`, `results/metrics/`, and checkpoint folders into the local repository’s corresponding `results/` paths. Run `python -m src.evaluation.compare_results` and the visualization scripts locally to generate the final comparison tables and figures.

The only manual inputs are the Drive folder path, optional repository URL/path, and Colab GPU tier. No full training is started by the preparation work in this session.

The local capped-dataset classical compute gate is approximately 28–32 minutes per ten-epoch model. Colab T4/A100 transformer duration is intentionally not hardcoded: local RTX 3050 timings are not a reliable cross-GPU conversion. The notebook prints a local baseline, then records real per-epoch Colab timing; before NT it uses the observed DNABERT-2 Colab epoch time multiplied by the local NT/DNABERT-2 ratio (242.67/127.40 ≈ 1.90), or clearly labels the fallback estimate if no DNABERT-2 timing log exists yet.
