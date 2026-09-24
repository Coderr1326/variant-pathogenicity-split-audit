# M3 Gene-Only Baseline vs. Transformers (ROC-AUC)

Split,M3 AUROC,DNABERT-2 AUROC,NT AUROC,DNABERT-2 minus M3,NT minus M3
A,0.6903,0.9137,0.9343,0.2233,0.2440
B†,0.5000,0.6740,0.9016,0.1740,0.4016
C,0.6044,0.7772,0.8519,0.1729,0.2475


† DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. M3 SPLIT-B AUROC is 0.5 by construction (constant global prior; no test gene seen in train).
