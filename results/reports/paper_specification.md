# Paper specification

Paper: “Enhancing Genetic Variant Pathogenicity Classification using Mutation-Aware Transformer Architectures” (Anjana M. P. and Arun K. S., AICAPS 2026, DOI 10.1109/AICAPS68631.2026.11452691).

## Explicit methodology

- **Explicitly stated by paper:** The dataset uses ClinVar labels and paired reference and mutation-altered DNA sequences from fixed genomic windows. The reference FASTA is from NCBI and the variant data are from ClinVar.
- **Explicitly stated by paper:** The intended scope is SNVs plus small insertions and deletions; variants with conflicting interpretations or uncertain significance are excluded.
- **Explicitly stated by paper:** Windows are centered at the mutation position, reverse-complemented for reverse-strand genes, and validated against the reference allele.
- **Explicitly stated by paper:** The reported selected window is 100 bp; the window sweep is 30, 50, 100, and 200 bp.
- **Explicitly stated by paper:** Reference and mutated sequences are concatenated and tokenized with each model’s tokenizer. Maximum input length is twice the window size, with padding/truncation.
- **Explicitly stated by paper:** The data are randomly split 80%/20% into training and validation sets using a fixed seed. The paper does not state the seed or whether stratification was used.
- **Explicitly stated by paper:** DNABERT-2 and Nucleotide Transformer use pretrained encoders and a linear task head. The paper’s Figure 1 and Section IV describe one sequence-level pathogenicity prediction per variant.
- **Explicitly stated by paper:** AdamW, learning rate 2e-5, batch size 16, 10 epochs, class-weighted cross-entropy with weights computed from the training labels, validation F1 model selection, and mixed precision when GPU resources are available.
- **Explicitly stated by paper:** Baselines are CNN, BiLSTM, CNN+BiLSTM, and ensemble CNN-RNN. Their descriptions specify embeddings, convolutional motif extraction, recurrent global modeling, and an LSTM/BiLSTM/GRU ensemble, but not exact layer widths.

## Paper-internal inconsistency

- **Explicitly stated by paper:** Section III-B describes token-level mutation tagging and masked token loss.
- **Explicitly stated by paper:** Section IV, Figure 1, the confusion matrices, and the reported metrics describe sequence-level binary classification.
- **Implementation assumption (paper underspecified):** This reproduction implements sequence-level binary classification: pretrained encoder, pooled representation, linear head, and one Benign/Pathogenic prediction per variant. The inconsistency is retained in `reports/consistency_audit.md`.

## External reconstruction

- **Implementation assumption (paper underspecified):** Reference build is GRCh38/hg38, as required by the reproduction brief.
- **Recovered from external source:** GENCODE v46 human GRCh38 GTF, released/modified 2024-05-13, supplies gene-symbol strand lookup. The final 100 bp dataset has 94,426 `+` and 74,501 `-` variants. Reverse-strand reference and mutated windows are both reverse-complemented.
- **Recovered from external source:** The COSMIC Census page was reachable in the one permitted attempt, but no directly usable public gene-list download was exposed. The final substitute remains the static list of 50 canonical cancer genes and 50 non-cancer genes in `src/data/build_dataset.py`.
- **Implementation assumption (paper underspecified):** A per-gene cap of 5,804 variants was applied, equal to 3% of the uncapped 193,453-row dataset rounded up. Sampling is fixed-seed and stratified within each gene by its observed label ratio. The capped dataset contains 168,927 variants: 90,786 benign and 78,141 pathogenic.

### Strand verification examples

For each example, the first sequence is the genomic-oriented window and the second is the stored gene-oriented window. The stored sequence equals the exact reverse complement.

| Gene | Variant | Genomic window | Stored minus-strand window |
|---|---|---|---|
| BRCA1 | chr17:43039471 | `GGGAGAGATTCTGGCCTAGAAAACCTGGAGAAGGCTCTGTAGGGGTGAGCGGGAAGGATTGGGGCCTGTAGAGATGACAGATCAGGACTTGCCTAGTAGA` | `TCTACTAGGCAAGTCCTGATCTGTCATCTCTACAGAGGCCCCAATCCTTCCCGCTCACCCCTACAGAGCCTTCTCCAGGTTTTCTAGGCCAGAATCTCTCCC` |
| DMD | chrX:31119356 | `GTCTGATATGTTGTGAAAATGCAGTAAAACTGAAGTTTAAAAAAATAATTCGTAAATGTTACAGTGTTGGTGTTAAAACACAATATATTATGATACTCAA` | `TTGAGTATCATAATATATTGTGTTTTAACACCAACACTGTAACATTTACGAATTATTTTTTTAAACTTCAGTTTTACTGCATTTTCACAACATATCAGAC` |
| FBN1 | chr15:48408316 | `TCCCTAACATTTATGGGTATATAACTTTTAAGGATACTAATGCAGCATCAACCCAATTGTCCTTTATTTTGGCTCATCTAATTTACAGACATGATTTCTG` | `CAGAAATCATGTCTGTAAATTAGATGAGCCAAAATAAAGGACAATTGGGTTGATGCTGCATTAGTATCCTTAAAAGTTATATACCCATAAATGTTAGGGA` |
- **Recovered from external source:** The ClinVar GRCh38 archive contains weekly 2025 releases; `clinvar_20250630.vcf.gz` is the nearest available release to 2025-07-01 (one day before the target) and is recorded in the manifest.
- **Implementation assumption (paper underspecified):** DNABERT-2 is `zhihan1996/DNABERT-2-117M`. Nucleotide Transformer selection tries `InstaDeepAI/nucleotide-transformer-v2-500m-multi-species`, then falls back to `InstaDeepAI/nucleotide-transformer-v2-100m-multi-species` if the 6 GB GPU cannot support it.

## Reproduction controls

- **Implementation assumption (paper underspecified):** Random seed is 42, with stratified 80/20 splitting.
- **Implementation assumption (paper underspecified):** AdamW uses standard PyTorch/Transformers defaults for betas and epsilon, weight decay 0.01, linear decay with no warmup, gradient clipping at 1.0, and dropout 0.1 where a model exposes dropout.
- **Implementation assumption (paper underspecified):** Classical baseline widths and kernels are small reproducible defaults documented in the model configuration and manifest.

## Results reported by the paper

The body text is the comparison target. At 100 bp, the paper explicitly states that the CNN baseline achieved F1-score 0.702 and ROC-AUC 0.823 (Section IV-C). No concrete accuracy is given for this CNN baseline. At 100 bp, the paper reports DNABERT-2 overall accuracy 0.89, weighted F1 0.89, AUROC 0.95, benign precision/recall/F1 0.88/0.95/0.91, and pathogenic precision/recall/F1 0.91/0.81/0.86. Nucleotide Transformer is reported at accuracy 0.93, weighted F1 0.93, AUROC 0.97, with benign precision/recall/F1 0.95/0.93/0.94 and pathogenic precision/recall/F1 0.90/0.92/0.91; it remains pending in this reproduction.

The paper’s figures imply 8,161 evaluated variants (4,818 benign and 3,343 pathogenic); this is a validation-size sanity anchor, not a sampling target. If interpreted as 20% of the source data, the approximate total is 40,800.
