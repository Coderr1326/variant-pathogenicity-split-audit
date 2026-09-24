# ROC-AUC by Locus-Overlap Subset

Run,A full,A shared,A novel,A shared-novel,B full,C full,C shared,C novel,C shared-novel
M3,0.6903,0.7056,0.6775,0.0281,0.5000,0.6044,0.6164,0.5961,0.0203
NT,0.9343,0.9506,0.9228,0.0278,0.9016,0.8519,0.8807,0.8349,0.0458
DNABERT-2† seed 42,0.9137,0.9345,0.8985,0.0360,0.6740,0.7772,0.8163,0.7570,0.0593
DNABERT-2† seed 1,-,-,-,-,0.6796,-,-,-,-
DNABERT-2† seed 2,-,-,-,-,0.6785,-,-,-,-
DNABERT-2† seed 3,-,-,-,-,0.6640,-,-,-,-
DNABERT-2† mean,-,-,-,-,0.6740,-,-,-,-
DNABERT-2† min,-,-,-,-,0.6640,-,-,-,-
DNABERT-2† max,-,-,-,-,0.6796,-,-,-,-


† DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. DNABERT-2 SPLIT-A and SPLIT-C have one run each (seed 42, fixed LR); seeds 1-3 exist for SPLIT-B only, so mean/min/max are over the four SPLIT-B seeds. 'shared-novel' is the point difference. Subsets are not random samples, so differences are descriptive. Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities).
