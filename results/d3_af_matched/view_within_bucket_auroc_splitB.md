# Plain Within-Bucket AUROC (no resampling): SPLIT-B

Run,missing,"(0,1e-4)","[1e-4,1e-3)"
M3,0.5000,0.5000,0.5000
M4,0.5000,0.5000,0.5000
NT,0.9210,0.7874,0.7256
DNABERT-2† seed 42,0.6718,0.6327,0.6417
DNABERT-2† seed 1,0.6727,0.6492,0.6578
DNABERT-2† seed 2,0.6756,0.6349,0.6346
DNABERT-2† seed 3,0.6589,0.6262,0.6428


Only buckets with >= 50 of both classes are shown. † DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models. Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities).
