# AF-Matched AUROC vs. Full Test Set: SPLIT-C

Run,Full test AUROC,Matched mean,Matched 2.5%,Matched 97.5%,Delta (matched - full),Matched n
M3,0.6044,0.6017,0.5999,0.6030,-0.0026,30240
M4,0.5337,0.5000,0.5000,0.5000,-0.0337,30240
NT,0.8519,0.8567,0.8560,0.8576,0.0048,30240
DNABERT-2† seed 42,0.7772,0.7777,0.7766,0.7789,0.0004,30240


Matched subset = per-bucket class downsampling over all eligible buckets incl. 'missing' (100 resamples, seeds 0-99); mean and 2.5-97.5 percentile of the resample AUROCs. DNABERT-2 SPLIT-A/C: single run (seed 42, fixed LR); SPLIT-B mean/min/max are over seeds 42/1/2/3. † DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models. Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities).
