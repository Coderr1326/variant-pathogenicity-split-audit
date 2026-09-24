# D3 AF-matched export

Generated 2026-09-21 22:02:23 +0530. Numbers only; no mechanism is claimed. Floats have 6 decimals; blank = not defined.

Caveat: ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models.

DNABERT-2 SPLIT-B results (seeds 42/1/2/3) used LR warmup 0.1 after the initial fixed-LR SPLIT-B attempt collapsed; DNABERT-2 SPLIT-A and SPLIT-C used the original fixed-LR schedule.

## 1. af_definition_and_parameters

Where AF comes from, how it is combined, and every fixed parameter.

```csv
key,value
af_source,"data/raw/clinvar_20250630.vcf.gz INFO fields AF_ESP, AF_EXAC, AF_TGP (the release that built the dataset; no download)"
variant_key,"CHROM, POS, REF, ALT (dataset vid = chrom_pos_ref_alt); all 168,927 dataset variants found in the VCF"
af_combination,"max over the non-missing of AF_ESP, AF_EXAC, AF_TGP; none present = missing bucket (not AF = 0)"
bucket_cutoffs,"missing; ==0; (0,1e-4); [1e-4,1e-3); [1e-3,1e-2); [1e-2,5e-2); >=5e-2"
min_per_class_bucket,50
n_resamples,100
resample_seeds,"0..99, default_rng(r) per resample"
ci_percentiles,"2.5, 97.5 of the resample AUROCs"
m4_alpha,1.000000
matching,"per bucket the majority class is downsampled to the minority count; union over eligible buckets; variants: all_eligible (incl. missing), af_known_only"
models_scored,"M3; M4; NT; DNABERT-2 seeds 42/1/2/3 (SPLIT-B warmup 0.1; SPLIT-A/C single run seed 42, fixed LR)"
python,3.14.7
generated,2026-09-21 22:02:23 +0530
```

## 2. af_coverage

Fraction of variants with non-missing AF, by set and label. train_X/test_X = SPLIT-X manifest train/test set; all_variants = full 168,927-variant dataset.

```csv
set,label,n,n_with_af,frac_with_af
all_variants,all,168927,35687,0.211257
all_variants,pathogenic,78141,10714,0.137111
all_variants,benign,90786,24973,0.275075
train_A,all,135141,28525,0.211076
train_A,pathogenic,62513,8581,0.137267
train_A,benign,72628,19944,0.274605
test_A,all,33786,7162,0.211981
test_A,pathogenic,15628,2133,0.136486
test_A,benign,18158,5029,0.276958
train_B,all,135147,28948,0.214196
train_B,pathogenic,62525,8844,0.141447
train_B,benign,72622,20104,0.276831
test_B,all,33780,6739,0.199497
test_B,pathogenic,15616,1870,0.119749
test_B,benign,18164,4869,0.268058
train_C,all,135391,21518,0.158932
train_C,pathogenic,62728,4427,0.070575
train_C,benign,72663,17091,0.235209
test_C,all,33536,14169,0.422501
test_C,pathogenic,15413,6287,0.407902
test_C,benign,18123,7882,0.434917
```

## 3. af_coverage_by_source

Per-source coverage over all variants.

```csv
source_field,n_with_value,frac_with_value,n_pathogenic_with_value,n_benign_with_value
AF_ESP,9955,0.058931,3154,6801
AF_EXAC,27424,0.162342,9813,17611
AF_TGP,14495,0.085806,2462,12033
```

## 4. bucket_label_counts

Pathogenic/benign counts per AF bucket for each split's train set and test set. 'missing' is its own category.

```csv
set,bucket,n_pathogenic,n_benign
train_A,missing,53932,52684
train_A,0,60,79
train_A,"(0,1e-4)",5628,8961
train_A,"[1e-4,1e-3)",2447,4074
train_A,"[1e-3,1e-2)",414,3160
train_A,"[1e-2,5e-2)",28,1471
train_A,>=5e-2,4,2199
test_A,missing,13495,13129
test_A,0,14,16
test_A,"(0,1e-4)",1386,2209
test_A,"[1e-4,1e-3)",636,1034
test_A,"[1e-3,1e-2)",89,815
test_A,"[1e-2,5e-2)",6,389
test_A,>=5e-2,2,566
train_B,missing,53681,52518
train_B,0,62,74
train_B,"(0,1e-4)",5726,8996
train_B,"[1e-4,1e-3)",2564,4104
train_B,"[1e-3,1e-2)",454,3132
train_B,"[1e-2,5e-2)",33,1530
train_B,>=5e-2,5,2268
test_B,missing,13746,13295
test_B,0,12,21
test_B,"(0,1e-4)",1288,2174
test_B,"[1e-4,1e-3)",519,1004
test_B,"[1e-3,1e-2)",49,843
test_B,"[1e-2,5e-2)",1,330
test_B,>=5e-2,1,497
train_C,missing,58301,55572
train_C,0,30,65
train_C,"(0,1e-4)",3302,7683
train_C,"[1e-4,1e-3)",953,2799
train_C,"[1e-3,1e-2)",126,2711
train_C,"[1e-2,5e-2)",14,1488
train_C,>=5e-2,2,2345
test_C,missing,9126,10241
test_C,0,44,30
test_C,"(0,1e-4)",3712,3487
test_C,"[1e-4,1e-3)",2130,2309
test_C,"[1e-3,1e-2)",377,1264
test_C,"[1e-2,5e-2)",20,372
test_C,>=5e-2,4,420
```

## 5. m4_af_only

M4 = AF-only baseline (pathogenic rate per bucket, train only, alpha=1, threshold 0.5). rate_* = fitted smoothed rates.

```csv
split,accuracy,weighted_f1,roc_auc,pathogenic_recall,m3_roc_auc_reference,global_train_rate,rate_missing,rate_0,"rate_(0,1e-4)","rate_[1e-4,1e-3)","rate_[1e-3,1e-2)","rate_[1e-2,5e-2)",rate_>=5e-2
A,0.548274,0.508967,0.576244,0.863514,0.690329,0.462576,0.505852,0.431876,0.385775,0.375263,0.115934,0.018975,0.002025
B,0.551066,0.508205,0.579581,0.880251,0.500000,0.462644,0.505475,0.455932,0.388947,0.384535,0.126697,0.021396,0.002402
C,0.507156,0.504998,0.533670,0.592098,0.604367,0.463310,0.511982,0.317326,0.300607,0.254054,0.044561,0.009623,0.001049
```

## 6. matching_eligibility

Buckets used for matching / within-bucket AUROC (>= 50 of both classes) per split test set.

```csv
split,bucket,n_pathogenic,n_benign,eligible,reason_if_skipped
A,missing,13495,13129,true,
A,0,14,16,false,smaller class has 14 (< 50)
A,"(0,1e-4)",1386,2209,true,
A,"[1e-4,1e-3)",636,1034,true,
A,"[1e-3,1e-2)",89,815,true,
A,"[1e-2,5e-2)",6,389,false,smaller class has 6 (< 50)
A,>=5e-2,2,566,false,smaller class has 2 (< 50)
B,missing,13746,13295,true,
B,0,12,21,false,smaller class has 12 (< 50)
B,"(0,1e-4)",1288,2174,true,
B,"[1e-4,1e-3)",519,1004,true,
B,"[1e-3,1e-2)",49,843,false,smaller class has 49 (< 50)
B,"[1e-2,5e-2)",1,330,false,smaller class has 1 (< 50)
B,>=5e-2,1,497,false,smaller class has 1 (< 50)
C,missing,9126,10241,true,
C,0,44,30,false,smaller class has 30 (< 50)
C,"(0,1e-4)",3712,3487,true,
C,"[1e-4,1e-3)",2130,2309,true,
C,"[1e-3,1e-2)",377,1264,true,
C,"[1e-2,5e-2)",20,372,false,smaller class has 20 (< 50)
C,>=5e-2,4,420,false,smaller class has 4 (< 50)
```

## 7. matched_auroc

Full-test AUROC vs. AF-matched AUROC (mean and 2.5-97.5 percentile over 100 resamples), delta = matched mean - full, matched_n = subset size. variant all_eligible includes the 'missing' bucket; af_known_only excludes it. DNABERT-2 mean/min/max are over the four SPLIT-B seeds.

```csv
run,split,variant,full_auroc,matched_mean,matched_p2_5,matched_p97_5,delta_matched_minus_full,matched_n
M3,A,all_eligible,0.690329,0.689308,0.688118,0.690321,-0.001022,30480
M3,A,af_known_only,0.690329,0.594729,0.588440,0.601005,-0.095601,4222
M4,A,all_eligible,0.576244,0.500000,0.500000,0.500000,-0.076244,30480
M4,A,af_known_only,0.576244,0.500000,0.500000,0.500000,-0.076244,4222
NT,A,all_eligible,0.934311,0.935569,0.935124,0.935918,0.001259,30480
NT,A,af_known_only,0.934311,0.816689,0.811105,0.821235,-0.117621,4222
DNABERT-2 seed 42,A,all_eligible,0.913665,0.914454,0.913908,0.915045,0.000790,30480
DNABERT-2 seed 42,A,af_known_only,0.913665,0.767091,0.763135,0.771434,-0.146574,4222
M3,B,all_eligible,0.500000,0.500000,0.500000,0.500000,0.000000,30204
M3,B,af_known_only,0.500000,0.500000,0.500000,0.500000,0.000000,3614
M4,B,all_eligible,0.579581,0.500000,0.500000,0.500000,-0.079581,30204
M4,B,af_known_only,0.579581,0.500000,0.500000,0.500000,-0.079581,3614
NT,B,all_eligible,0.901574,0.905655,0.905074,0.906219,0.004081,30204
NT,B,af_known_only,0.901574,0.769525,0.763050,0.774937,-0.132050,3614
DNABERT-2 seed 42,B,all_eligible,0.673967,0.667583,0.666567,0.668866,-0.006383,30204
DNABERT-2 seed 42,B,af_known_only,0.673967,0.635543,0.627201,0.643110,-0.038423,3614
DNABERT-2 seed 1,B,all_eligible,0.679560,0.670255,0.669125,0.671387,-0.009305,30204
DNABERT-2 seed 1,B,af_known_only,0.679560,0.651515,0.642804,0.659797,-0.028045,3614
DNABERT-2 seed 2,B,all_eligible,0.678490,0.670831,0.669631,0.672003,-0.007658,30204
DNABERT-2 seed 2,B,af_known_only,0.678490,0.634878,0.625982,0.641442,-0.043612,3614
DNABERT-2 seed 3,B,all_eligible,0.664000,0.655651,0.654389,0.656832,-0.008350,30204
DNABERT-2 seed 3,B,af_known_only,0.664000,0.631038,0.621624,0.639421,-0.032962,3614
M3,C,all_eligible,0.604367,0.601721,0.599938,0.602979,-0.002645,30240
M3,C,af_known_only,0.604367,0.578415,0.575968,0.580955,-0.025952,11988
M4,C,all_eligible,0.533670,0.500000,0.500000,0.500000,-0.033670,30240
M4,C,af_known_only,0.533670,0.500000,0.500000,0.500000,-0.033670,11988
NT,C,all_eligible,0.851860,0.856687,0.855962,0.857621,0.004827,30240
NT,C,af_known_only,0.851860,0.792584,0.790764,0.794452,-0.059276,11988
DNABERT-2 seed 42,C,all_eligible,0.777240,0.777662,0.776593,0.778916,0.000423,30240
DNABERT-2 seed 42,C,af_known_only,0.777240,0.731236,0.729278,0.732939,-0.046004,11988
DNABERT-2 mean over seeds,B,all_eligible,0.674004,0.666080,0.664928,0.667272,-0.007924,30204
DNABERT-2 min over seeds,B,all_eligible,0.664000,0.655651,0.654389,0.656832,-0.009305,30204
DNABERT-2 max over seeds,B,all_eligible,0.679560,0.670831,0.669631,0.672003,-0.006383,30204
DNABERT-2 mean over seeds,B,af_known_only,0.674004,0.638244,0.629403,0.645943,-0.035761,3614
DNABERT-2 min over seeds,B,af_known_only,0.664000,0.631038,0.621624,0.639421,-0.043612,3614
DNABERT-2 max over seeds,B,af_known_only,0.679560,0.651515,0.642804,0.659797,-0.028045,3614
```

## 8. within_bucket_auroc

Plain AUROC inside each eligible bucket (no resampling).

```csv
run,split,bucket,n,n_pathogenic,n_benign,auroc
M3,A,missing,26624,13495,13129,0.706681
M3,A,"(0,1e-4)",3595,1386,2209,0.617672
M3,A,"[1e-4,1e-3)",1670,636,1034,0.547253
M3,A,"[1e-3,1e-2)",904,89,815,0.577108
M4,A,missing,26624,13495,13129,0.500000
M4,A,"(0,1e-4)",3595,1386,2209,0.500000
M4,A,"[1e-4,1e-3)",1670,636,1034,0.500000
M4,A,"[1e-3,1e-2)",904,89,815,0.500000
NT,A,missing,26624,13495,13129,0.952888
NT,A,"(0,1e-4)",3595,1386,2209,0.846801
NT,A,"[1e-4,1e-3)",1670,636,1034,0.766517
NT,A,"[1e-3,1e-2)",904,89,815,0.666010
DNABERT-2 seed 42,A,missing,26624,13495,13129,0.936423
DNABERT-2 seed 42,A,"(0,1e-4)",3595,1386,2209,0.792434
DNABERT-2 seed 42,A,"[1e-4,1e-3)",1670,636,1034,0.722054
DNABERT-2 seed 42,A,"[1e-3,1e-2)",904,89,815,0.667250
M3,B,missing,27041,13746,13295,0.500000
M3,B,"(0,1e-4)",3462,1288,2174,0.500000
M3,B,"[1e-4,1e-3)",1523,519,1004,0.500000
M4,B,missing,27041,13746,13295,0.500000
M4,B,"(0,1e-4)",3462,1288,2174,0.500000
M4,B,"[1e-4,1e-3)",1523,519,1004,0.500000
NT,B,missing,27041,13746,13295,0.920992
NT,B,"(0,1e-4)",3462,1288,2174,0.787356
NT,B,"[1e-4,1e-3)",1523,519,1004,0.725647
DNABERT-2 seed 42,B,missing,27041,13746,13295,0.671766
DNABERT-2 seed 42,B,"(0,1e-4)",3462,1288,2174,0.632706
DNABERT-2 seed 42,B,"[1e-4,1e-3)",1523,519,1004,0.641682
DNABERT-2 seed 1,B,missing,27041,13746,13295,0.672717
DNABERT-2 seed 1,B,"(0,1e-4)",3462,1288,2174,0.649222
DNABERT-2 seed 1,B,"[1e-4,1e-3)",1523,519,1004,0.657781
DNABERT-2 seed 2,B,missing,27041,13746,13295,0.675558
DNABERT-2 seed 2,B,"(0,1e-4)",3462,1288,2174,0.634931
DNABERT-2 seed 2,B,"[1e-4,1e-3)",1523,519,1004,0.634598
DNABERT-2 seed 3,B,missing,27041,13746,13295,0.658915
DNABERT-2 seed 3,B,"(0,1e-4)",3462,1288,2174,0.626229
DNABERT-2 seed 3,B,"[1e-4,1e-3)",1523,519,1004,0.642793
M3,C,missing,19367,9126,10241,0.613969
M3,C,"(0,1e-4)",7199,3712,3487,0.579361
M3,C,"[1e-4,1e-3)",4439,2130,2309,0.573029
M3,C,"[1e-3,1e-2)",1641,377,1264,0.600489
M4,C,missing,19367,9126,10241,0.500000
M4,C,"(0,1e-4)",7199,3712,3487,0.500000
M4,C,"[1e-4,1e-3)",4439,2130,2309,0.500000
M4,C,"[1e-3,1e-2)",1641,377,1264,0.500000
NT,C,missing,19367,9126,10241,0.895345
NT,C,"(0,1e-4)",7199,3712,3487,0.826057
NT,C,"[1e-4,1e-3)",4439,2130,2309,0.755242
NT,C,"[1e-3,1e-2)",1641,377,1264,0.675841
DNABERT-2 seed 42,C,missing,19367,9126,10241,0.808086
DNABERT-2 seed 42,C,"(0,1e-4)",7199,3712,3487,0.754820
DNABERT-2 seed 42,C,"[1e-4,1e-3)",4439,2130,2309,0.705837
DNABERT-2 seed 42,C,"[1e-3,1e-2)",1641,377,1264,0.662947
```

## 9. skipped

Models, runs and buckets not used, with reasons.

```csv
item,kind,scope,reason
CNN,model,all,no saved probabilities (weights and history only); predictions were not regenerated
BiLSTM,model,all,no saved probabilities (weights and history only); predictions were not regenerated
CNN+BiLSTM,model,all,no saved probabilities (weights and history only); predictions were not regenerated
Ensemble,model,all,no saved probabilities (weights and history only); predictions were not regenerated
DNABERT-2 seed 1 warmup 0,run,B,collapsed run; never loaded
"DNABERT-2 seeds 1, 2, 3",run,A;C,no seed runs exist for SPLIT-A/C
0,bucket,A,smaller class has 14 (< 50)
"[1e-2,5e-2)",bucket,A,smaller class has 6 (< 50)
>=5e-2,bucket,A,smaller class has 2 (< 50)
0,bucket,B,smaller class has 12 (< 50)
"[1e-3,1e-2)",bucket,B,smaller class has 49 (< 50)
"[1e-2,5e-2)",bucket,B,smaller class has 1 (< 50)
>=5e-2,bucket,B,smaller class has 1 (< 50)
0,bucket,C,smaller class has 30 (< 50)
"[1e-2,5e-2)",bucket,C,smaller class has 20 (< 50)
>=5e-2,bucket,C,smaller class has 4 (< 50)
```

## 10. integrity_check

sha256 of every file under results/, scripts/, docs/ before vs. after (results/checkpoints: size + mtime only). Added files are this task's outputs.

```csv
check,value
files_hashed_sha256_before,257
files_hashed_sha256_now,309
modified_preexisting_files,none
removed_files,none
added_files_count,52
added_files,results/d3_af_matched/af_coverage.csv;results/d3_af_matched/af_coverage_by_source.csv;results/d3_af_matched/af_table.csv;results/d3_af_matched/bucket_label_counts.csv;results/d3_af_matched/d3_grouped_bar_full_vs_matched.png;results/d3_af_matched/m4_af_only/split_A/final_result.json;results/d3_af_matched/m4_af_only/split_A/predictions.parquet;results/d3_af_matched/m4_af_only/split_B/final_result.json;results/d3_af_matched/m4_af_only/split_B/predictions.parquet;results/d3_af_matched/m4_af_only/split_C/final_result.json;results/d3_af_matched/m4_af_only/split_C/predictions.parquet;results/d3_af_matched/m4_af_only_metrics.csv;results/d3_af_matched/matched_auroc.csv;results/d3_af_matched/matched_resample_aurocs.csv;results/d3_af_matched/matching_eligibility.csv;results/d3_af_matched/parameters.csv;results/d3_af_matched/view_af_coverage.csv;results/d3_af_matched/view_af_coverage.md;results/d3_af_matched/view_af_coverage.png;results/d3_af_matched/view_bucket_label_counts.csv;results/d3_af_matched/view_bucket_label_counts.md;results/d3_af_matched/view_bucket_label_counts.png;results/d3_af_matched/view_m4_af_only.csv;results/d3_af_matched/view_m4_af_only.md;results/d3_af_matched/view_m4_af_only.png;results/d3_af_matched/view_matched_afknown_only.csv;results/d3_af_matched/view_matched_afknown_only.md;results/d3_af_matched/view_matched_afknown_only.png;results/d3_af_matched/view_matched_auroc_splitA.csv;results/d3_af_matched/view_matched_auroc_splitA.md;results/d3_af_matched/view_matched_auroc_splitA.png;results/d3_af_matched/view_matched_auroc_splitB.csv;results/d3_af_matched/view_matched_auroc_splitB.md;results/d3_af_matched/view_matched_auroc_splitB.png;results/d3_af_matched/view_matched_auroc_splitC.csv;results/d3_af_matched/view_matched_auroc_splitC.md;results/d3_af_matched/view_matched_auroc_splitC.png;results/d3_af_matched/view_matching_eligibility.csv;results/d3_af_matched/view_matching_eligibility.md;results/d3_af_matched/view_matching_eligibility.png;results/d3_af_matched/view_within_bucket_auroc_splitA.csv;results/d3_af_matched/view_within_bucket_auroc_splitA.md;results/d3_af_matched/view_within_bucket_auroc_splitA.png;results/d3_af_matched/view_within_bucket_auroc_splitB.csv;results/d3_af_matched/view_within_bucket_auroc_splitB.md;results/d3_af_matched/view_within_bucket_auroc_splitB.png;results/d3_af_matched/view_within_bucket_auroc_splitC.csv;results/d3_af_matched/view_within_bucket_auroc_splitC.md;results/d3_af_matched/view_within_bucket_auroc_splitC.png;results/d3_af_matched/within_bucket_auroc.csv;results/handoff/d3_export.md;scripts/run_d3_af_matched.py
checkpoints_size_mtime_before,92
checkpoints_size_mtime_now,92
checkpoints_changed_added_removed,none
result,PASS: no pre-existing file changed
```
