# AF-Matched AUROC vs. Full Test Set: SPLIT-B

Run,Full test AUROC,Matched mean,Matched 2.5%,Matched 97.5%,Delta (matched - full),Matched n
M3,0.5000,0.5000,0.5000,0.5000,0.0000,30204
M4,0.5796,0.5000,0.5000,0.5000,-0.0796,30204
NT,0.9016,0.9057,0.9051,0.9062,0.0041,30204
DNABERT-2† seed 42,0.6740,0.6676,0.6666,0.6689,-0.0064,30204
DNABERT-2† seed 1,0.6796,0.6703,0.6691,0.6714,-0.0093,30204
DNABERT-2† seed 2,0.6785,0.6708,0.6696,0.6720,-0.0077,30204
DNABERT-2† seed 3,0.6640,0.6557,0.6544,0.6568,-0.0083,30204
DNABERT-2† mean,0.6740,0.6661,0.6649,0.6673,-0.0079,30204
DNABERT-2† min,0.6640,0.6557,0.6544,0.6568,-0.0093,30204
DNABERT-2† max,0.6796,0.6708,0.6696,0.6720,-0.0064,30204


Matched subset = per-bucket class downsampling over all eligible buckets incl. 'missing' (100 resamples, seeds 0-99); mean and 2.5-97.5 percentile of the resample AUROCs. DNABERT-2 SPLIT-A/C: single run (seed 42, fixed LR); SPLIT-B mean/min/max are over seeds 42/1/2/3. † DNABERT-2 SPLIT-B required an LR warmup fix after an initial run experienced representation collapse (see methodology); SPLIT-A and SPLIT-C used the original fixed-LR schedule. ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models. Not shown: CNN, BiLSTM, CNN+BiLSTM, Ensemble (no saved probabilities).
