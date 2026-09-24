# AF-Matched AUROC, AF-Known Buckets Only (excluding 'missing')

Split,Run,Full test AUROC,Matched mean,Matched 2.5%,Matched 97.5%,Delta,Matched n
A,M3,0.6903,0.5947,0.5884,0.6010,-0.0956,4222
A,M4,0.5762,0.5000,0.5000,0.5000,-0.0762,4222
A,NT,0.9343,0.8167,0.8111,0.8212,-0.1176,4222
A,DNABERT-2† seed 42,0.9137,0.7671,0.7631,0.7714,-0.1466,4222
B,M3,0.5000,0.5000,0.5000,0.5000,0.0000,3614
B,M4,0.5796,0.5000,0.5000,0.5000,-0.0796,3614
B,NT,0.9016,0.7695,0.7631,0.7749,-0.1320,3614
B,DNABERT-2† seed 42,0.6740,0.6355,0.6272,0.6431,-0.0384,3614
B,DNABERT-2† seed 1,0.6796,0.6515,0.6428,0.6598,-0.0280,3614
B,DNABERT-2† seed 2,0.6785,0.6349,0.6260,0.6414,-0.0436,3614
B,DNABERT-2† seed 3,0.6640,0.6310,0.6216,0.6394,-0.0330,3614
C,M3,0.6044,0.5784,0.5760,0.5810,-0.0260,11988
C,M4,0.5337,0.5000,0.5000,0.5000,-0.0337,11988
C,NT,0.8519,0.7926,0.7908,0.7945,-0.0593,11988
C,DNABERT-2† seed 42,0.7772,0.7312,0.7293,0.7329,-0.0460,11988
B,DNABERT-2† mean,0.6740,0.6382,0.6294,0.6459,-0.0358,3614
B,DNABERT-2† min,0.6640,0.6310,0.6216,0.6394,-0.0436,3614
B,DNABERT-2† max,0.6796,0.6515,0.6428,0.6598,-0.0280,3614


Same procedure as the main matched table but the 'missing' bucket is left out. † DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models.
