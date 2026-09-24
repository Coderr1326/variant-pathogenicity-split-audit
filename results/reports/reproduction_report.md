# Reproduction report

## A. Paper methodology

The paper’s explicit methodology is summarized in `reports/paper_specification.md`: ClinVar-labeled SNVs and small indels, centered windows, reference/mutant paired sequences, pretrained DNABERT-2 and Nucleotide Transformer encoders, classical sequence baselines, weighted cross-entropy, AdamW at 2e-5, batch size 16, 10 epochs, 80/20 validation split, and validation-F1 model selection.

## B. Implementation

The implementation uses sequence-level binary classification, as required by the paper’s Figure 1 and evaluation artifacts. Code is under `src/`; commands are in `README.md`.

## C. External reconstruction

Reference build: GRCh38. ClinVar release: `clinvar_20250706.vcf.gz`, the closest archived weekly release to 2025-07-01. The gene list is an explicit static substitute: 50 canonical cancer genes and 50 non-cancer Mendelian-disease genes, because the paper supplies no list.

## D. Assumptions

All assumptions use the required tags in `reports/paper_specification.md`. The key assumptions are seed 42, stratified split, AdamW defaults plus weight decay 0.01, linear schedule with no warmup, gradient clipping 1.0, dropout 0.1, and documented baseline widths/kernels.

## E. Dataset reproduction

The corrected canonical 100 bp dataset contains 168,927 deduplicated variants: 90,786 benign and 78,141 pathogenic. The uncapped builder scanned 3,512,796 ClinVar records, excluded 127,791 labels, dropped 1,405 unsupported variants, and recorded zero reference-allele failures. GENCODE v46 strand lookup yields 94,426 plus-strand and 74,501 minus-strand variants; reverse-complement examples were verified. A reproducible cap of 5,804 variants per gene reduced the uncapped 193,453 rows while preserving each gene’s observed label ratio. This remains about 4.1 times the approximate 40,800 source-size anchor, with 53.7%/46.3% class balance versus the paper’s approximate 59%/41% validation anchor.

## F. Model reproduction

DNABERT-2 uses `zhihan1996/DNABERT-2-117M`. Nucleotide Transformer attempts the 500M checkpoint and falls back to `InstaDeepAI/nucleotide-transformer-v2-100m-multi-species` because the available GPU has 6 GB VRAM. Classical models are implemented in `src/models/classical.py`.

## G. Training reproduction

AdamW, learning rate 2e-5, batch size 16, ten epochs, training-set-only class weights, validation F1 selection, mixed precision when CUDA is available, and clipping at 1.0. Full runtime and compute-gate estimates are recorded when the experiments run.

## H. Evaluation reproduction

Metrics are calculated from raw validation predictions: accuracy, ROC-AUC, weighted precision/recall/F1, and per-class precision/recall/F1. Confusion matrices are asserted to sum to the evaluated labels and are generated programmatically.

## I. Paper-vs-reproduction results

Generated automatically at `results/tables/paper_vs_reproduction.csv` by `src/evaluation/compare_results.py`. No paper values are used as model targets.

## J. Discrepancy analysis

The updated capped-dataset compute gate measured 100 bp one-epoch timings on 5% of the data: CNN 10.02 s, BiLSTM 8.46 s, CNN+BiLSTM 8.75 s, and ensemble 9.66 s. Extrapolated ten-epoch classical-baseline cost is approximately 28–32 minutes per model. Transformer training is prepared for Colab; the prior local estimates were approximately 7.1 hours for DNABERT-2 and 13.5 hours for NT on the 6 GB RTX 3050. Full training was not started.

## K. Fidelity classification

Methodology fidelity: MEDIUM — preprocessing now includes verified GENCODE reverse-strand handling, model/training/evaluation paths are implemented according to the stated procedure, but full 10-epoch results were not completed and the paper’s token-level/sequence-level conflict required an explicit resolution.

Dataset identity fidelity: LOW — the exact author gene list is unavailable, the capped reconstructed dataset remains 4.1× the approximate anchor and has different class balance, and the ClinVar release is reconstructed rather than author-pinned.

Overall classification: PARTIAL — the data pipeline and all model paths pass sanity checks, but full training/evaluation is deferred and dataset identity is low-confidence. Exact reproduction is not supportable because the paper provides no data-availability statement, gene list, release identifier, or checkpoint identifier.

## Colab preparation status

`notebooks/train_transformers_colab.ipynb` is the transformer execution artifact. It mounts Drive, imports the shared local `src/` modules, keeps training batch size 16, saves epoch checkpoints, supports resume, logs epoch wall time, checks Drive space using loaded-model parameter size, estimates NT multi-session cost, and stops after an epoch when the five-hour budget is reached. A local dry-run loaded the capped parquet and completed one transformer training step with CUDA. Full training was not run.
