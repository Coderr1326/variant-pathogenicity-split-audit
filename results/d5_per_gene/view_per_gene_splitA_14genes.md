# Per-Gene AUROC on SPLIT-A Test Set (the same 14 genes as SPLIT-B)

Gene,n (SPLIT-A test),n pathogenic,n benign,M3 AUROC,NT AUROC,DNABERT-2 AUROC
ATM,1151,518,633,0.5000,0.9457,0.9122
BRCA2,1202,827,375,0.5000,0.9454,0.9276
DMD,1183,477,706,0.5000,0.9349,0.8822
NF1,1131,681,450,0.5000,0.9599,0.9426
CDH1,503,183,320,0.5000,0.9033,0.8848
ALK,472,57,415,0.5000,0.7892,0.6764
CACNA1S,355,86,269,0.5000,0.8324,0.7748
ARID1A,186,49,137,0.5000,0.8417,0.7599
BRAF,138,30,108,0.5000,0.8846,0.8086
F8,133,112,21,insufficient,insufficient,insufficient
CTNNB1,105,44,61,0.5000,0.9892,0.9769
ABL1,105,10,95,insufficient,insufficient,insufficient
FBXW7,20,4,16,insufficient,insufficient,insufficient
BCR,8,0,8,insufficient,insufficient,insufficient


insufficient = fewer than 30 pathogenic or benign variants in the gene. SPLIT-A has no per-seed DNABERT-2 runs; CIs are in per_gene_metrics_splitA_14genes.csv. Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities).
