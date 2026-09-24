# Plain Within-Bucket AUROC (no resampling): SPLIT-A

Run,missing,"(0,1e-4)","[1e-4,1e-3)","[1e-3,1e-2)"
M3,0.7067,0.6177,0.5473,0.5771
M4,0.5000,0.5000,0.5000,0.5000
NT,0.9529,0.8468,0.7665,0.6660
DNABERT-2† seed 42,0.9364,0.7924,0.7221,0.6673


Only buckets with >= 50 of both classes are shown. † DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models. Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities).
