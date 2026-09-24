# Leakage Delta (D1): SPLIT-A minus SPLIT-B/C -- M3 Gene-Only Baseline

Model,D1 A-B Acc,D1 A-B F1,D1 A-B AUC,D1 A-B PathRecall,D1 A-C Acc,D1 A-C F1,D1 A-C AUC,D1 A-C PathRecall
M3 (gene-only),0.1040,0.2645,0.1903,0.5796,0.0584,0.0638,0.0860,0.1306


M3 SPLIT-B is a chance floor by construction (constant global prior), so its A-B deltas measure the loss of gene-identity signal, not model degradation.
