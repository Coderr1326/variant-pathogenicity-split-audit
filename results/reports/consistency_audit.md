# Consistency audit

## Abstract versus body

The abstract reports DNABERT-2 accuracy 0.89, F1 0.85, and AUROC 0.94, while the body reports accuracy 0.89, F1 0.89, and AUROC 0.95. The abstract reports Nucleotide Transformer accuracy 0.92, F1 0.91, and AUROC 0.967, while the body reports accuracy 0.93, weighted F1 0.93, and AUROC 0.97.

**Implementation assumption (paper underspecified):** The body-text values and confusion-matrix-derived values are the comparison target, as instructed by the reproduction brief. Abstract values are retained here as a documented discrepancy and are not silently substituted.

## Task formulation

**Explicitly stated by paper:** Section III-B describes token-level mutation classification with masked non-mutated positions.

**Explicitly stated by paper:** Figure 1, Section IV, confusion matrices, and all headline metrics describe sequence-level binary classification.

**Implementation assumption (paper underspecified):** The code uses sequence-level binary classification, one prediction per variant, because it is the only formulation consistent with the evaluation artifacts.

## Other audit notes

The paper does not provide a data-availability statement, exact gene list, ClinVar release, transformer checkpoint identifiers, split seed, or complete baseline architecture hyperparameters. These are recorded as reconstruction assumptions rather than published facts.

