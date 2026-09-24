# M3 Gene-Only Baseline Results at 100 bp (no sequence input)

Model,Split,Accuracy,Weighted Precision,Weighted Recall,Weighted F1,ROC-AUC,Benign P/R/F1,Pathogenic P/R/F1
M3 (gene-only),A,0.6417,0.6406,0.6417,0.6405,0.6903,0.658/0.695/0.676,0.621/0.580/0.599
M3 (gene-only),B,0.5377,0.2891,0.5377,0.3761,0.5000,0.538/1.000/0.699,0.000/0.000/0.000
M3 (gene-only),C,0.5833,0.5797,0.5833,0.5767,0.6044,0.598/0.697/0.644,0.558/0.449/0.498


SPLIT-B: chance floor: constant global prior (0 of 14 test genes seen in train); AUROC = 0.5 exactly and tp = 0 by construction (expected, not a bug).
