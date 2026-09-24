# Split-Shift Delta (SPLIT-A minus SPLIT-B/C)

Model,A->B Acc,A->B F1,A->B AUC,A->B PathRecall,A->C Acc,A->C F1,A->C AUC,A->C PathRecall
M3 (gene-only)‡,0.1040,0.2645,0.1903,0.5796,0.0584,0.0638,0.0860,0.1306
CNN,0.0198,0.0252,0.0355,0.0604,0.0735,0.0838,0.0959,0.1365
BiLSTM,-0.0084,0.0004,0.0015,0.0316,0.0620,0.1277,0.0777,0.2340
CNN+BiLSTM,-0.0073,-0.0005,0.0142,0.0783,0.0512,0.0503,0.0808,0.0354
Ensemble,0.0031,0.0095,0.0110,0.0641,0.0518,0.0568,0.0758,0.0970
DNABERT-2†,0.2424,0.2706,0.2397,0.4687,0.1456,0.1490,0.1364,0.2106
NT,0.0368,0.0373,0.0327,0.0640,0.0794,0.0828,0.0825,0.1739


† DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. DNABERT-2 SPLIT-B AUROC over seeds 42/1/2/3 (warmup 0.1): mean 0.6740, min-max 0.6640-0.6796 (A->B AUROC delta over these seeds: mean 0.2397, min-max 0.2341-0.2497); the table row uses seed 42. ‡ M3 SPLIT-B is a chance floor by construction (constant global prior; 0 of 14 test genes seen in train), so its A->B deltas measure loss of gene signal, not model degradation.
