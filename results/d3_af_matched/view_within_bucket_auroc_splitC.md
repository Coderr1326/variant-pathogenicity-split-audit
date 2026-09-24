# Plain Within-Bucket AUROC (no resampling): SPLIT-C

Run,missing,"(0,1e-4)","[1e-4,1e-3)","[1e-3,1e-2)"
M3,0.6140,0.5794,0.5730,0.6005
M4,0.5000,0.5000,0.5000,0.5000
NT,0.8953,0.8261,0.7552,0.6758
DNABERT-2† seed 42,0.8081,0.7548,0.7058,0.6629


Only buckets with >= 50 of both classes are shown. † DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models. Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities).
