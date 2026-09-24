# Gene-Composition Control: Macro Within-Gene AUROC, Locus-Shared vs. Locus-Novel

Run,Split,Genes used,Macro shared,Macro novel,Shared-novel
M3,A,40,0.5000,0.5000,0.0000
M3,C,35,0.5000,0.5000,0.0000
NT,A,40,0.9360,0.9177,0.0184
NT,C,35,0.8565,0.8368,0.0196
DNABERT-2 seed 42,A,40,0.9200,0.8969,0.0231
DNABERT-2 seed 42,C,35,0.7939,0.7707,0.0232


Genes with >= 30 pathogenic and >= 30 benign test variants in BOTH subsets. M3 is constant within a gene, so its within-gene AUROC is 0.5. DNABERT-2 SPLIT-A/C: single run, seed 42, fixed LR. Subsets are not random samples, so differences are descriptive. Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities).
