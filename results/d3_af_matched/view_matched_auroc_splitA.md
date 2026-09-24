# AF-Matched AUROC vs. Full Test Set: SPLIT-A

Run,Full test AUROC,Matched mean,Matched 2.5%,Matched 97.5%,Delta (matched - full),Matched n
M3,0.6903,0.6893,0.6881,0.6903,-0.0010,30480
M4,0.5762,0.5000,0.5000,0.5000,-0.0762,30480
NT,0.9343,0.9356,0.9351,0.9359,0.0013,30480
DNABERT-2† seed 42,0.9137,0.9145,0.9139,0.9150,0.0008,30480


Matched subset = per-bucket class downsampling over all eligible buckets incl. 'missing' (100 resamples, seeds 0-99); mean and 2.5-97.5 percentile of the resample AUROCs. DNABERT-2 SPLIT-A/C: single run (seed 42, fixed LR); SPLIT-B mean/min/max are over seeds 42/1/2/3. † DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models. Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities).
