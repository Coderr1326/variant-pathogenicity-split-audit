# D7 locus-overlap export

Generated 2026-09-21 21:15:07 +0530. Numbers only; no mechanism is claimed. Floats have 6 decimals; blank = not defined. Subsets are not random samples, so differences are descriptive.

DNABERT-2 SPLIT-B results (seeds 42/1/2/3) used LR warmup 0.1 after the initial fixed-LR SPLIT-B attempt collapsed; DNABERT-2 SPLIT-A and SPLIT-C used the original fixed-LR schedule.

## 1. locus_reproduction_check

Locus-shared test variants per split; expected values are those recorded when 36.1% / 34.3% were first computed (docs/terminology.md). The script asserts every match.

```csv
split,n_test,n_locus_shared,n_locus_novel,frac_shared,pct_shared,expected_pct,expected_n_shared,matches_expected
A,33786,12186,21600,0.360682,36.068194,36.100000,12186,true
B,33780,0,33780,0.000000,0.000000,0.000000,0,true
C,33536,11495,22041,0.342766,34.276598,34.300000,11495,true
```

## 2. skipped_models

Models without saved probabilities; no predictions were regenerated (no-inference rule). The DNABERT-2 warmup-0 run is never loaded.

```csv
model,reason
CNN,no saved probabilities (only weights and per-epoch history)
BiLSTM,no saved probabilities (only weights and per-epoch history)
CNN+BiLSTM,no saved probabilities (only weights and per-epoch history)
Ensemble,no saved probabilities (only weights and per-epoch history)
```

## 3. subset_composition

Test-set subsets. top5_genes_count = gene:variant count, ordered by count.

```csv
split,subset,n,n_pathogenic,pathogenic_fraction,n_genes,top5_genes_count
A,full,33786,15628,0.462558,100,BRCA2:1202;DMD:1183;RYR1:1160;BRCA1:1159;TSC2:1158
A,shared,12186,6504,0.533727,93,MSH2:736;MSH6:697;TSC2:512;BRCA1:503;APC:452
A,novel,21600,9124,0.422407,100,RYR1:935;DMD:903;BRCA2:855;ATM:799;NF1:732
B,full,33780,15616,0.462285,14,ATM:5804;BRCA2:5804;DMD:5804;NF1:5804;CDH1:2673
C,full,33536,15413,0.459596,96,TSC2:1735;MSH2:1706;ATM:1287;MSH6:1271;MLH1:1189
C,shared,11495,5498,0.478295,91,MSH2:937;MSH6:715;TSC2:640;MLH1:539;APC:417
C,novel,22041,9915,0.449843,96,TSC2:1095;MET:917;ATM:916;MSH2:769;NF1:732
```

## 4. auroc_by_subset

AUROC within each subset; ci = 95% percentile bootstrap (1000 resamples of the subset's variants, seed 0). DNABERT-2 SPLIT-A/C: one run (seed 42, fixed LR); seeds 1/2/3 exist for SPLIT-B only (warmup 0.1). mean/min/max rows are over the four SPLIT-B seeds (no CI).

```csv
run,split,subset,n,n_pathogenic,auroc,ci_low,ci_high
M3,A,full,33786,15628,0.690329,0.684472,0.695922
M3,A,shared,12186,6504,0.705638,0.696176,0.714637
M3,A,novel,21600,9124,0.677519,0.670721,0.684655
M3,B,full,33780,15616,0.500000,0.500000,0.500000
M3,C,full,33536,15413,0.604367,0.598383,0.610219
M3,C,shared,11495,5498,0.616387,0.606266,0.626221
M3,C,novel,22041,9915,0.596068,0.588784,0.603500
NT,A,full,33786,15628,0.934311,0.931804,0.937132
NT,A,shared,12186,6504,0.950605,0.946772,0.953921
NT,A,novel,21600,9124,0.922844,0.918717,0.926579
NT,B,full,33780,15616,0.901574,0.898289,0.904907
NT,C,full,33536,15413,0.851860,0.847712,0.856286
NT,C,shared,11495,5498,0.880686,0.874359,0.887443
NT,C,novel,22041,9915,0.834865,0.829315,0.840156
DNABERT-2 seed 42,A,full,33786,15628,0.913665,0.910345,0.916939
DNABERT-2 seed 42,A,shared,12186,6504,0.934504,0.929726,0.938807
DNABERT-2 seed 42,A,novel,21600,9124,0.898473,0.893795,0.902852
DNABERT-2 seed 42,B,full,33780,15616,0.673967,0.668416,0.679613
DNABERT-2 seed 42,C,full,33536,15413,0.777240,0.772371,0.782064
DNABERT-2 seed 42,C,shared,11495,5498,0.816318,0.809201,0.823953
DNABERT-2 seed 42,C,novel,22041,9915,0.757026,0.751297,0.763532
DNABERT-2 seed 1,B,full,33780,15616,0.679560,0.674322,0.685136
DNABERT-2 seed 2,B,full,33780,15616,0.678490,0.672793,0.684408
DNABERT-2 seed 3,B,full,33780,15616,0.664000,0.658017,0.669705
DNABERT-2 mean over seeds,B,full,33780,15616,0.674004,,
DNABERT-2 min over seeds,B,full,33780,15616,0.664000,,
DNABERT-2 max over seeds,B,full,33780,15616,0.679560,,
```

## 5. shared_minus_novel

Point difference in AUROC, locus-shared minus locus-novel (descriptive; subsets are not random samples).

```csv
run,split,auroc_shared,auroc_novel,difference
M3,A,0.705638,0.677519,0.028119
M3,C,0.616387,0.596068,0.020319
NT,A,0.950605,0.922844,0.027761
NT,C,0.880686,0.834865,0.045821
DNABERT-2 seed 42,A,0.934504,0.898473,0.036031
DNABERT-2 seed 42,C,0.816318,0.757026,0.059292
```

## 6. within_gene_control

Macro-mean within-gene AUROC over genes with >= 30 pathogenic and >= 30 benign test variants in BOTH subsets (gene lists in section 7).

```csv
run,split,n_genes_used,macro_within_gene_auroc_shared,macro_within_gene_auroc_novel,difference
M3,A,40,0.500000,0.500000,0.000000
M3,C,35,0.500000,0.500000,0.000000
NT,A,40,0.936029,0.917670,0.018359
NT,C,35,0.856461,0.836811,0.019650
DNABERT-2 seed 42,A,40,0.920002,0.896896,0.023106
DNABERT-2 seed 42,C,35,0.793899,0.770715,0.023185
```

## 7. within_gene_control_genes

Genes used in section 6, with per-subset counts.

```csv
split,gene,n_shared,n_shared_pathogenic,n_shared_benign,n_novel,n_novel_pathogenic,n_novel_benign
A,APC,452,211,241,682,327,355
A,ATM,352,186,166,799,332,467
A,ATP7B,106,51,55,294,141,153
A,BRCA1,503,382,121,656,356,300
A,BRCA2,347,264,83,855,563,292
A,CDH1,240,81,159,263,102,161
A,CFTR,215,133,82,378,187,191
A,CHEK2,254,138,116,199,94,105
A,COL1A1,143,99,44,355,202,153
A,DMD,280,127,153,903,350,553
A,FBN1,395,309,86,715,420,295
A,HBB,120,49,71,153,40,113
A,KCNH2,146,74,72,285,114,171
A,KCNQ1,90,54,36,316,72,244
A,LDLR,366,297,69,303,175,128
A,MECP2,133,77,56,135,65,70
A,MEN1,161,63,98,196,103,93
A,MLH1,438,280,158,308,182,126
A,MSH2,736,452,284,341,193,148
A,MSH6,697,406,291,451,256,195
A,NF1,399,261,138,732,420,312
A,PALB2,342,205,137,357,187,170
A,PMS2,307,119,188,280,136,144
A,PTCH1,247,87,160,421,157,264
A,PTEN,269,196,73,183,92,91
A,RAD51C,109,41,68,165,62,103
A,RAD51D,88,35,53,118,54,64
A,RB1,188,92,96,320,131,189
A,RET,139,30,109,329,85,244
A,RYR1,225,53,172,935,241,694
A,SCN1A,211,175,36,459,315,144
A,SCN5A,101,39,62,383,143,240
A,SDHB,77,35,42,94,33,61
A,SLC26A4,78,41,37,183,90,93
A,SMAD4,104,30,74,186,57,129
A,STK11,176,65,111,150,49,101
A,TP53,317,224,93,188,85,103
A,TSC1,191,83,108,364,146,218
A,TSC2,512,195,317,646,290,356
A,VHL,114,56,58,102,42,60
C,APC,417,190,227,721,371,350
C,ATM,371,206,165,916,510,406
C,ATP7B,84,49,35,225,134,91
C,BRCA1,333,262,71,350,253,97
C,BRCA2,239,192,47,629,474,155
C,CDH1,219,84,135,301,156,145
C,CFTR,145,95,50,257,147,110
C,CHEK2,237,149,88,227,140,87
C,COL1A1,128,80,48,314,158,156
C,DMD,197,91,106,653,282,371
C,FBN1,263,160,103,566,309,257
C,KCNH2,151,70,81,223,85,138
C,KCNQ1,117,77,40,453,113,340
C,LDLR,225,168,57,171,97,74
C,MEN1,183,74,109,234,93,141
C,MLH1,539,191,348,650,185,465
C,MSH2,937,400,537,769,288,481
C,MSH6,715,435,280,556,373,183
C,NF1,403,244,159,732,396,336
C,PALB2,344,217,127,504,334,170
C,PMS2,331,175,156,392,201,191
C,PTCH1,290,95,195,583,245,338
C,PTEN,168,92,76,109,59,50
C,RAD51C,91,50,41,219,125,94
C,RAD51D,94,43,51,161,77,84
C,RB1,151,56,95,264,98,166
C,RET,172,58,114,428,185,243
C,RYR1,179,56,123,588,231,357
C,SCN5A,107,44,63,444,238,206
C,SDHB,81,36,45,125,47,78
C,STK11,207,74,133,265,135,130
C,TP53,322,198,124,197,72,125
C,TSC1,239,96,143,447,236,211
C,TSC2,640,220,420,1095,524,571
C,VHL,113,54,59,155,71,84
```

## 8. within_gene_control_per_gene_auroc

Per-gene AUROC behind section 6.

```csv
run,split,gene,subset,n,auroc
M3,A,APC,shared,452,0.500000
M3,A,APC,novel,682,0.500000
M3,A,ATM,shared,352,0.500000
M3,A,ATM,novel,799,0.500000
M3,A,ATP7B,shared,106,0.500000
M3,A,ATP7B,novel,294,0.500000
M3,A,BRCA1,shared,503,0.500000
M3,A,BRCA1,novel,656,0.500000
M3,A,BRCA2,shared,347,0.500000
M3,A,BRCA2,novel,855,0.500000
M3,A,CDH1,shared,240,0.500000
M3,A,CDH1,novel,263,0.500000
M3,A,CFTR,shared,215,0.500000
M3,A,CFTR,novel,378,0.500000
M3,A,CHEK2,shared,254,0.500000
M3,A,CHEK2,novel,199,0.500000
M3,A,COL1A1,shared,143,0.500000
M3,A,COL1A1,novel,355,0.500000
M3,A,DMD,shared,280,0.500000
M3,A,DMD,novel,903,0.500000
M3,A,FBN1,shared,395,0.500000
M3,A,FBN1,novel,715,0.500000
M3,A,HBB,shared,120,0.500000
M3,A,HBB,novel,153,0.500000
M3,A,KCNH2,shared,146,0.500000
M3,A,KCNH2,novel,285,0.500000
M3,A,KCNQ1,shared,90,0.500000
M3,A,KCNQ1,novel,316,0.500000
M3,A,LDLR,shared,366,0.500000
M3,A,LDLR,novel,303,0.500000
M3,A,MECP2,shared,133,0.500000
M3,A,MECP2,novel,135,0.500000
M3,A,MEN1,shared,161,0.500000
M3,A,MEN1,novel,196,0.500000
M3,A,MLH1,shared,438,0.500000
M3,A,MLH1,novel,308,0.500000
M3,A,MSH2,shared,736,0.500000
M3,A,MSH2,novel,341,0.500000
M3,A,MSH6,shared,697,0.500000
M3,A,MSH6,novel,451,0.500000
M3,A,NF1,shared,399,0.500000
M3,A,NF1,novel,732,0.500000
M3,A,PALB2,shared,342,0.500000
M3,A,PALB2,novel,357,0.500000
M3,A,PMS2,shared,307,0.500000
M3,A,PMS2,novel,280,0.500000
M3,A,PTCH1,shared,247,0.500000
M3,A,PTCH1,novel,421,0.500000
M3,A,PTEN,shared,269,0.500000
M3,A,PTEN,novel,183,0.500000
M3,A,RAD51C,shared,109,0.500000
M3,A,RAD51C,novel,165,0.500000
M3,A,RAD51D,shared,88,0.500000
M3,A,RAD51D,novel,118,0.500000
M3,A,RB1,shared,188,0.500000
M3,A,RB1,novel,320,0.500000
M3,A,RET,shared,139,0.500000
M3,A,RET,novel,329,0.500000
M3,A,RYR1,shared,225,0.500000
M3,A,RYR1,novel,935,0.500000
M3,A,SCN1A,shared,211,0.500000
M3,A,SCN1A,novel,459,0.500000
M3,A,SCN5A,shared,101,0.500000
M3,A,SCN5A,novel,383,0.500000
M3,A,SDHB,shared,77,0.500000
M3,A,SDHB,novel,94,0.500000
M3,A,SLC26A4,shared,78,0.500000
M3,A,SLC26A4,novel,183,0.500000
M3,A,SMAD4,shared,104,0.500000
M3,A,SMAD4,novel,186,0.500000
M3,A,STK11,shared,176,0.500000
M3,A,STK11,novel,150,0.500000
M3,A,TP53,shared,317,0.500000
M3,A,TP53,novel,188,0.500000
M3,A,TSC1,shared,191,0.500000
M3,A,TSC1,novel,364,0.500000
M3,A,TSC2,shared,512,0.500000
M3,A,TSC2,novel,646,0.500000
M3,A,VHL,shared,114,0.500000
M3,A,VHL,novel,102,0.500000
M3,C,APC,shared,417,0.500000
M3,C,APC,novel,721,0.500000
M3,C,ATM,shared,371,0.500000
M3,C,ATM,novel,916,0.500000
M3,C,ATP7B,shared,84,0.500000
M3,C,ATP7B,novel,225,0.500000
M3,C,BRCA1,shared,333,0.500000
M3,C,BRCA1,novel,350,0.500000
M3,C,BRCA2,shared,239,0.500000
M3,C,BRCA2,novel,629,0.500000
M3,C,CDH1,shared,219,0.500000
M3,C,CDH1,novel,301,0.500000
M3,C,CFTR,shared,145,0.500000
M3,C,CFTR,novel,257,0.500000
M3,C,CHEK2,shared,237,0.500000
M3,C,CHEK2,novel,227,0.500000
M3,C,COL1A1,shared,128,0.500000
M3,C,COL1A1,novel,314,0.500000
M3,C,DMD,shared,197,0.500000
M3,C,DMD,novel,653,0.500000
M3,C,FBN1,shared,263,0.500000
M3,C,FBN1,novel,566,0.500000
M3,C,KCNH2,shared,151,0.500000
M3,C,KCNH2,novel,223,0.500000
M3,C,KCNQ1,shared,117,0.500000
M3,C,KCNQ1,novel,453,0.500000
M3,C,LDLR,shared,225,0.500000
M3,C,LDLR,novel,171,0.500000
M3,C,MEN1,shared,183,0.500000
M3,C,MEN1,novel,234,0.500000
M3,C,MLH1,shared,539,0.500000
M3,C,MLH1,novel,650,0.500000
M3,C,MSH2,shared,937,0.500000
M3,C,MSH2,novel,769,0.500000
M3,C,MSH6,shared,715,0.500000
M3,C,MSH6,novel,556,0.500000
M3,C,NF1,shared,403,0.500000
M3,C,NF1,novel,732,0.500000
M3,C,PALB2,shared,344,0.500000
M3,C,PALB2,novel,504,0.500000
M3,C,PMS2,shared,331,0.500000
M3,C,PMS2,novel,392,0.500000
M3,C,PTCH1,shared,290,0.500000
M3,C,PTCH1,novel,583,0.500000
M3,C,PTEN,shared,168,0.500000
M3,C,PTEN,novel,109,0.500000
M3,C,RAD51C,shared,91,0.500000
M3,C,RAD51C,novel,219,0.500000
M3,C,RAD51D,shared,94,0.500000
M3,C,RAD51D,novel,161,0.500000
M3,C,RB1,shared,151,0.500000
M3,C,RB1,novel,264,0.500000
M3,C,RET,shared,172,0.500000
M3,C,RET,novel,428,0.500000
M3,C,RYR1,shared,179,0.500000
M3,C,RYR1,novel,588,0.500000
M3,C,SCN5A,shared,107,0.500000
M3,C,SCN5A,novel,444,0.500000
M3,C,SDHB,shared,81,0.500000
M3,C,SDHB,novel,125,0.500000
M3,C,STK11,shared,207,0.500000
M3,C,STK11,novel,265,0.500000
M3,C,TP53,shared,322,0.500000
M3,C,TP53,novel,197,0.500000
M3,C,TSC1,shared,239,0.500000
M3,C,TSC1,novel,447,0.500000
M3,C,TSC2,shared,640,0.500000
M3,C,TSC2,novel,1095,0.500000
M3,C,VHL,shared,113,0.500000
M3,C,VHL,novel,155,0.500000
NT,A,APC,shared,452,0.966726
NT,A,APC,novel,682,0.973020
NT,A,ATM,shared,352,0.950576
NT,A,ATM,novel,799,0.942784
NT,A,ATP7B,shared,106,0.904813
NT,A,ATP7B,novel,294,0.916099
NT,A,BRCA1,shared,503,0.939466
NT,A,BRCA1,novel,656,0.936573
NT,A,BRCA2,shared,347,0.951990
NT,A,BRCA2,novel,855,0.943180
NT,A,CDH1,shared,240,0.925382
NT,A,CDH1,novel,263,0.883205
NT,A,CFTR,shared,215,0.958005
NT,A,CFTR,novel,378,0.895736
NT,A,CHEK2,shared,254,0.980635
NT,A,CHEK2,novel,199,0.936069
NT,A,COL1A1,shared,143,0.955464
NT,A,COL1A1,novel,355,0.942212
NT,A,DMD,shared,280,0.954557
NT,A,DMD,novel,903,0.928520
NT,A,FBN1,shared,395,0.958870
NT,A,FBN1,novel,715,0.949556
NT,A,HBB,shared,120,0.929577
NT,A,HBB,novel,153,0.899336
NT,A,KCNH2,shared,146,0.939189
NT,A,KCNH2,novel,285,0.934852
NT,A,KCNQ1,shared,90,0.902778
NT,A,KCNQ1,novel,316,0.933459
NT,A,LDLR,shared,366,0.939150
NT,A,LDLR,novel,303,0.939955
NT,A,MECP2,shared,133,0.947820
NT,A,MECP2,novel,135,0.905495
NT,A,MEN1,shared,161,0.943149
NT,A,MEN1,novel,196,0.930160
NT,A,MLH1,shared,438,0.949458
NT,A,MLH1,novel,308,0.957091
NT,A,MSH2,shared,736,0.937368
NT,A,MSH2,novel,341,0.914683
NT,A,MSH6,shared,697,0.959567
NT,A,MSH6,novel,451,0.943329
NT,A,NF1,shared,399,0.968821
NT,A,NF1,novel,732,0.953289
NT,A,PALB2,shared,342,0.946306
NT,A,PALB2,novel,357,0.931771
NT,A,PMS2,shared,307,0.948462
NT,A,PMS2,novel,280,0.928462
NT,A,PTCH1,shared,247,0.914152
NT,A,PTCH1,novel,421,0.892757
NT,A,PTEN,shared,269,0.959883
NT,A,PTEN,novel,183,0.928213
NT,A,RAD51C,shared,109,0.922166
NT,A,RAD51C,novel,165,0.934075
NT,A,RAD51D,shared,88,0.861456
NT,A,RAD51D,novel,118,0.861111
NT,A,RB1,shared,188,0.945199
NT,A,RB1,novel,320,0.945555
NT,A,RET,shared,139,0.893272
NT,A,RET,novel,329,0.856557
NT,A,RYR1,shared,225,0.809237
NT,A,RYR1,novel,935,0.880643
NT,A,SCN1A,shared,211,0.934444
NT,A,SCN1A,novel,459,0.883003
NT,A,SCN5A,shared,101,0.929694
NT,A,SCN5A,novel,383,0.899913
NT,A,SDHB,shared,77,0.882313
NT,A,SDHB,novel,94,0.998013
NT,A,SLC26A4,shared,78,0.982861
NT,A,SLC26A4,novel,183,0.909916
NT,A,SMAD4,shared,104,0.977027
NT,A,SMAD4,novel,186,0.916769
NT,A,STK11,shared,176,0.916840
NT,A,STK11,novel,150,0.890281
NT,A,TP53,shared,317,0.923531
NT,A,TP53,novel,188,0.911479
NT,A,TSC1,shared,191,0.947345
NT,A,TSC1,novel,364,0.867789
NT,A,TSC2,shared,512,0.920877
NT,A,TSC2,novel,646,0.871029
NT,A,VHL,shared,114,0.962746
NT,A,VHL,novel,102,0.840873
NT,C,APC,shared,417,0.927104
NT,C,APC,novel,721,0.863596
NT,C,ATM,shared,371,0.920506
NT,C,ATM,novel,916,0.851082
NT,C,ATP7B,shared,84,0.819825
NT,C,ATP7B,novel,225,0.812695
NT,C,BRCA1,shared,333,0.832007
NT,C,BRCA1,novel,350,0.820464
NT,C,BRCA2,shared,239,0.853059
NT,C,BRCA2,novel,629,0.892950
NT,C,CDH1,shared,219,0.747972
NT,C,CDH1,novel,301,0.795889
NT,C,CFTR,shared,145,0.834947
NT,C,CFTR,novel,257,0.833271
NT,C,CHEK2,shared,237,0.911074
NT,C,CHEK2,novel,227,0.878900
NT,C,COL1A1,shared,128,0.896875
NT,C,COL1A1,novel,314,0.825584
NT,C,DMD,shared,197,0.812254
NT,C,DMD,novel,653,0.729971
NT,C,FBN1,shared,263,0.918022
NT,C,FBN1,novel,566,0.862176
NT,C,KCNH2,shared,151,0.785009
NT,C,KCNH2,novel,223,0.824723
NT,C,KCNQ1,shared,117,0.924026
NT,C,KCNQ1,novel,453,0.865669
NT,C,LDLR,shared,225,0.927527
NT,C,LDLR,novel,171,0.899276
NT,C,MEN1,shared,183,0.890156
NT,C,MEN1,novel,234,0.754823
NT,C,MLH1,shared,539,0.869426
NT,C,MLH1,novel,650,0.811555
NT,C,MSH2,shared,937,0.888259
NT,C,MSH2,novel,769,0.844263
NT,C,MSH6,shared,715,0.888957
NT,C,MSH6,novel,556,0.927907
NT,C,NF1,shared,403,0.921641
NT,C,NF1,novel,732,0.927572
NT,C,PALB2,shared,344,0.915871
NT,C,PALB2,novel,504,0.886844
NT,C,PMS2,shared,331,0.845311
NT,C,PMS2,novel,392,0.860775
NT,C,PTCH1,shared,290,0.812848
NT,C,PTCH1,novel,583,0.762891
NT,C,PTEN,shared,168,0.939788
NT,C,PTEN,novel,109,0.859322
NT,C,RAD51C,shared,91,0.850732
NT,C,RAD51C,novel,219,0.872851
NT,C,RAD51D,shared,94,0.836297
NT,C,RAD51D,novel,161,0.819264
NT,C,RB1,shared,151,0.841917
NT,C,RB1,novel,264,0.861999
NT,C,RET,shared,172,0.720811
NT,C,RET,novel,428,0.803381
NT,C,RYR1,shared,179,0.743031
NT,C,RYR1,novel,588,0.744273
NT,C,SCN5A,shared,107,0.809885
NT,C,SCN5A,novel,444,0.818145
NT,C,SDHB,shared,81,0.865432
NT,C,SDHB,novel,125,0.912984
NT,C,STK11,shared,207,0.758890
NT,C,STK11,novel,265,0.868604
NT,C,TP53,shared,322,0.897687
NT,C,TP53,novel,197,0.823222
NT,C,TSC1,shared,239,0.830201
NT,C,TSC1,novel,447,0.789079
NT,C,TSC2,shared,640,0.818496
NT,C,TSC2,novel,1095,0.794816
NT,C,VHL,shared,113,0.920276
NT,C,VHL,novel,155,0.787559
DNABERT-2 seed 42,A,APC,shared,452,0.949382
DNABERT-2 seed 42,A,APC,novel,682,0.943955
DNABERT-2 seed 42,A,ATM,shared,352,0.908958
DNABERT-2 seed 42,A,ATM,novel,799,0.913850
DNABERT-2 seed 42,A,ATP7B,shared,106,0.900178
DNABERT-2 seed 42,A,ATP7B,novel,294,0.895471
DNABERT-2 seed 42,A,BRCA1,shared,503,0.908983
DNABERT-2 seed 42,A,BRCA1,novel,656,0.924471
DNABERT-2 seed 42,A,BRCA2,shared,347,0.928715
DNABERT-2 seed 42,A,BRCA2,novel,855,0.926802
DNABERT-2 seed 42,A,CDH1,shared,240,0.900846
DNABERT-2 seed 42,A,CDH1,novel,263,0.868560
DNABERT-2 seed 42,A,CFTR,shared,215,0.931872
DNABERT-2 seed 42,A,CFTR,novel,378,0.888176
DNABERT-2 seed 42,A,CHEK2,shared,254,0.963581
DNABERT-2 seed 42,A,CHEK2,novel,199,0.916109
DNABERT-2 seed 42,A,COL1A1,shared,143,0.962121
DNABERT-2 seed 42,A,COL1A1,novel,355,0.918948
DNABERT-2 seed 42,A,DMD,shared,280,0.903916
DNABERT-2 seed 42,A,DMD,novel,903,0.874193
DNABERT-2 seed 42,A,FBN1,shared,395,0.946903
DNABERT-2 seed 42,A,FBN1,novel,715,0.906045
DNABERT-2 seed 42,A,HBB,shared,120,0.954585
DNABERT-2 seed 42,A,HBB,novel,153,0.936726
DNABERT-2 seed 42,A,KCNH2,shared,146,0.957395
DNABERT-2 seed 42,A,KCNH2,novel,285,0.912435
DNABERT-2 seed 42,A,KCNQ1,shared,90,0.853395
DNABERT-2 seed 42,A,KCNQ1,novel,316,0.954975
DNABERT-2 seed 42,A,LDLR,shared,366,0.934758
DNABERT-2 seed 42,A,LDLR,novel,303,0.927009
DNABERT-2 seed 42,A,MECP2,shared,133,0.957096
DNABERT-2 seed 42,A,MECP2,novel,135,0.903516
DNABERT-2 seed 42,A,MEN1,shared,161,0.936184
DNABERT-2 seed 42,A,MEN1,novel,196,0.897902
DNABERT-2 seed 42,A,MLH1,shared,438,0.940077
DNABERT-2 seed 42,A,MLH1,novel,308,0.939081
DNABERT-2 seed 42,A,MSH2,shared,736,0.912899
DNABERT-2 seed 42,A,MSH2,novel,341,0.861697
DNABERT-2 seed 42,A,MSH6,shared,697,0.937696
DNABERT-2 seed 42,A,MSH6,novel,451,0.930148
DNABERT-2 seed 42,A,NF1,shared,399,0.963435
DNABERT-2 seed 42,A,NF1,novel,732,0.928556
DNABERT-2 seed 42,A,PALB2,shared,342,0.933381
DNABERT-2 seed 42,A,PALB2,novel,357,0.904404
DNABERT-2 seed 42,A,PMS2,shared,307,0.931924
DNABERT-2 seed 42,A,PMS2,novel,280,0.895629
DNABERT-2 seed 42,A,PTCH1,shared,247,0.900503
DNABERT-2 seed 42,A,PTCH1,novel,421,0.870826
DNABERT-2 seed 42,A,PTEN,shared,269,0.963866
DNABERT-2 seed 42,A,PTEN,novel,183,0.917702
DNABERT-2 seed 42,A,RAD51C,shared,109,0.864419
DNABERT-2 seed 42,A,RAD51C,novel,165,0.893360
DNABERT-2 seed 42,A,RAD51D,shared,88,0.895957
DNABERT-2 seed 42,A,RAD51D,novel,118,0.900463
DNABERT-2 seed 42,A,RB1,shared,188,0.936707
DNABERT-2 seed 42,A,RB1,novel,320,0.912719
DNABERT-2 seed 42,A,RET,shared,139,0.870031
DNABERT-2 seed 42,A,RET,novel,329,0.747830
DNABERT-2 seed 42,A,RYR1,shared,225,0.785103
DNABERT-2 seed 42,A,RYR1,novel,935,0.829965
DNABERT-2 seed 42,A,SCN1A,shared,211,0.920476
DNABERT-2 seed 42,A,SCN1A,novel,459,0.873280
DNABERT-2 seed 42,A,SCN5A,shared,101,0.924318
DNABERT-2 seed 42,A,SCN5A,novel,383,0.874097
DNABERT-2 seed 42,A,SDHB,shared,77,0.895238
DNABERT-2 seed 42,A,SDHB,novel,94,0.954297
DNABERT-2 seed 42,A,SLC26A4,shared,78,0.974291
DNABERT-2 seed 42,A,SLC26A4,novel,183,0.932975
DNABERT-2 seed 42,A,SMAD4,shared,104,0.973874
DNABERT-2 seed 42,A,SMAD4,novel,186,0.917449
DNABERT-2 seed 42,A,STK11,shared,176,0.870963
DNABERT-2 seed 42,A,STK11,novel,150,0.865427
DNABERT-2 seed 42,A,TP53,shared,317,0.883113
DNABERT-2 seed 42,A,TP53,novel,188,0.850828
DNABERT-2 seed 42,A,TSC1,shared,191,0.895694
DNABERT-2 seed 42,A,TSC1,novel,364,0.817613
DNABERT-2 seed 42,A,TSC2,shared,512,0.878055
DNABERT-2 seed 42,A,TSC2,novel,646,0.852334
DNABERT-2 seed 42,A,VHL,shared,114,0.949200
DNABERT-2 seed 42,A,VHL,novel,102,0.896032
DNABERT-2 seed 42,C,APC,shared,417,0.808648
DNABERT-2 seed 42,C,APC,novel,721,0.772645
DNABERT-2 seed 42,C,ATM,shared,371,0.795822
DNABERT-2 seed 42,C,ATM,novel,916,0.807780
DNABERT-2 seed 42,C,ATP7B,shared,84,0.791837
DNABERT-2 seed 42,C,ATP7B,novel,225,0.701657
DNABERT-2 seed 42,C,BRCA1,shared,333,0.819751
DNABERT-2 seed 42,C,BRCA1,novel,350,0.837333
DNABERT-2 seed 42,C,BRCA2,shared,239,0.780696
DNABERT-2 seed 42,C,BRCA2,novel,629,0.825902
DNABERT-2 seed 42,C,CDH1,shared,219,0.766402
DNABERT-2 seed 42,C,CDH1,novel,301,0.758576
DNABERT-2 seed 42,C,CFTR,shared,145,0.822947
DNABERT-2 seed 42,C,CFTR,novel,257,0.777489
DNABERT-2 seed 42,C,CHEK2,shared,237,0.869738
DNABERT-2 seed 42,C,CHEK2,novel,227,0.851067
DNABERT-2 seed 42,C,COL1A1,shared,128,0.852865
DNABERT-2 seed 42,C,COL1A1,novel,314,0.750609
DNABERT-2 seed 42,C,DMD,shared,197,0.750363
DNABERT-2 seed 42,C,DMD,novel,653,0.686079
DNABERT-2 seed 42,C,FBN1,shared,263,0.868871
DNABERT-2 seed 42,C,FBN1,novel,566,0.797388
DNABERT-2 seed 42,C,KCNH2,shared,151,0.789418
DNABERT-2 seed 42,C,KCNH2,novel,223,0.774595
DNABERT-2 seed 42,C,KCNQ1,shared,117,0.893182
DNABERT-2 seed 42,C,KCNQ1,novel,453,0.853826
DNABERT-2 seed 42,C,LDLR,shared,225,0.808271
DNABERT-2 seed 42,C,LDLR,novel,171,0.771176
DNABERT-2 seed 42,C,MEN1,shared,183,0.733449
DNABERT-2 seed 42,C,MEN1,novel,234,0.736674
DNABERT-2 seed 42,C,MLH1,shared,539,0.770957
DNABERT-2 seed 42,C,MLH1,novel,650,0.719337
DNABERT-2 seed 42,C,MSH2,shared,937,0.777519
DNABERT-2 seed 42,C,MSH2,novel,769,0.759738
DNABERT-2 seed 42,C,MSH6,shared,715,0.817020
DNABERT-2 seed 42,C,MSH6,novel,556,0.858553
DNABERT-2 seed 42,C,NF1,shared,403,0.818074
DNABERT-2 seed 42,C,NF1,novel,732,0.814086
DNABERT-2 seed 42,C,PALB2,shared,344,0.827788
DNABERT-2 seed 42,C,PALB2,novel,504,0.803311
DNABERT-2 seed 42,C,PMS2,shared,331,0.811612
DNABERT-2 seed 42,C,PMS2,novel,392,0.840275
DNABERT-2 seed 42,C,PTCH1,shared,290,0.738246
DNABERT-2 seed 42,C,PTCH1,novel,583,0.694010
DNABERT-2 seed 42,C,PTEN,shared,168,0.838816
DNABERT-2 seed 42,C,PTEN,novel,109,0.873559
DNABERT-2 seed 42,C,RAD51C,shared,91,0.857073
DNABERT-2 seed 42,C,RAD51C,novel,219,0.798468
DNABERT-2 seed 42,C,RAD51D,shared,94,0.743730
DNABERT-2 seed 42,C,RAD51D,novel,161,0.756339
DNABERT-2 seed 42,C,RB1,shared,151,0.785526
DNABERT-2 seed 42,C,RB1,novel,264,0.719388
DNABERT-2 seed 42,C,RET,shared,172,0.700998
DNABERT-2 seed 42,C,RET,novel,428,0.662974
DNABERT-2 seed 42,C,RYR1,shared,179,0.662602
DNABERT-2 seed 42,C,RYR1,novel,588,0.637916
DNABERT-2 seed 42,C,SCN5A,shared,107,0.801948
DNABERT-2 seed 42,C,SCN5A,novel,444,0.667822
DNABERT-2 seed 42,C,SDHB,shared,81,0.808642
DNABERT-2 seed 42,C,SDHB,novel,125,0.815057
DNABERT-2 seed 42,C,STK11,shared,207,0.758687
DNABERT-2 seed 42,C,STK11,novel,265,0.814701
DNABERT-2 seed 42,C,TP53,shared,322,0.824332
DNABERT-2 seed 42,C,TP53,novel,197,0.825111
DNABERT-2 seed 42,C,TSC1,shared,239,0.736233
DNABERT-2 seed 42,C,TSC1,novel,447,0.681079
DNABERT-2 seed 42,C,TSC2,shared,640,0.758745
DNABERT-2 seed 42,C,TSC2,novel,1095,0.707554
DNABERT-2 seed 42,C,VHL,shared,113,0.795669
DNABERT-2 seed 42,C,VHL,novel,155,0.822938
```

## 9. parameters

Thresholds and resample counts, fixed in the script before any result was inspected.

```csv
key,value
min_class_count,30
n_bootstrap,1000
bootstrap_seed,0
ci,"95% percentile (2.5, 97.5)"
top_genes_listed,5
locus_definition,locus = chrom + '_' + pos; shared = locus occurs among the split's training variants (any alt allele); novel = otherwise
models,"M3; NT; DNABERT-2 seeds 42/1/2/3 (SPLIT-B warmup 0.1); DNABERT-2 SPLIT-A/C single run seed 42, fixed LR"
python,3.14.7
generated,2026-09-21 21:15:07 +0530
```

## 10. integrity_check

sha256 of every file under results/, scripts/, docs/ before vs. after (results/checkpoints: size + mtime only). Added files are this task's outputs.

```csv
check,value
files_hashed_sha256_before,225
files_hashed_sha256_now,257
modified_preexisting_files,none
removed_files,none
added_files,results/d7_locus_overlap/auroc_by_subset_long.csv;results/d7_locus_overlap/d7_grouped_bar_auroc_subsets.png;results/d7_locus_overlap/locus_overlap_fraction.csv;results/d7_locus_overlap/parameters.csv;results/d7_locus_overlap/shared_minus_novel.csv;results/d7_locus_overlap/subset_composition.csv;results/d7_locus_overlap/view_auroc_by_subset.csv;results/d7_locus_overlap/view_auroc_by_subset.md;results/d7_locus_overlap/view_auroc_by_subset.png;results/d7_locus_overlap/view_auroc_ci_by_subset.csv;results/d7_locus_overlap/view_auroc_ci_by_subset.md;results/d7_locus_overlap/view_auroc_ci_by_subset.png;results/d7_locus_overlap/view_locus_overlap_fraction.csv;results/d7_locus_overlap/view_locus_overlap_fraction.md;results/d7_locus_overlap/view_locus_overlap_fraction.png;results/d7_locus_overlap/view_subset_composition.csv;results/d7_locus_overlap/view_subset_composition.md;results/d7_locus_overlap/view_subset_composition.png;results/d7_locus_overlap/view_within_gene_control.csv;results/d7_locus_overlap/view_within_gene_control.md;results/d7_locus_overlap/view_within_gene_control.png;results/d7_locus_overlap/view_within_gene_control_genes_splitA.csv;results/d7_locus_overlap/view_within_gene_control_genes_splitA.md;results/d7_locus_overlap/view_within_gene_control_genes_splitA.png;results/d7_locus_overlap/view_within_gene_control_genes_splitC.csv;results/d7_locus_overlap/view_within_gene_control_genes_splitC.md;results/d7_locus_overlap/view_within_gene_control_genes_splitC.png;results/d7_locus_overlap/within_gene_control.csv;results/d7_locus_overlap/within_gene_control_genes.csv;results/d7_locus_overlap/within_gene_control_per_gene_auroc.csv;results/handoff/d7_export.md;scripts/run_d7_locus_overlap.py
checkpoints_size_mtime_before,92
checkpoints_size_mtime_now,92
checkpoints_changed_added_removed,none
result,PASS: no pre-existing file changed
```
