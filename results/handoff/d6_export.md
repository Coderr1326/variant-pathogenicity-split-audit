# D6 cross-gene sequence-similarity export

Generated 2026-09-21 22:19:29 +0530. Numbers only; no interpretation. Floats have 6 decimals; blank = not defined.

Per-window results (33,780 rows) are in results/d6_cross_gene_similarity/similarity_per_test_window.csv and are not embedded here.

## 1. parameters

Fixed before any result was inspected.

```csv
key,value
k,25
window_column,ref_sequence
sequence_columns_found,ref_sequence;mut_sequence;genomic_ref_sequence;genomic_mut_sequence
thresholds,0.5;0.8;0.95
top_pairs_listed,15
kmer_encoding,exact 2-bit-per-base integers (50 bits for K=25): a collision-free hash; k-mers containing N are skipped; distinct k-mers per window (set semantics)
variants,"forward = k-mers as written (primary); canonical = min(k-mer, reverse complement)"
similarity,max over single train windows of |distinct 25-mers of the test window found in that train window| / |distinct valid 25-mers of the test window|
index,inverted index over train windows: sorted k-mer array + searchsorted posting lists (no all-pairs comparison)
gene_pair_attribution,"each high-similarity test window is attributed to the distinct train genes among its tied best train windows, weight 1/(number of distinct genes)"
split,data/splits/split_manifest_b.json (sha256:7490bb4a5b7732905f583fa432268ed334ccc96246de1cb2e50d22f0f9de5924)
train_windows,135147
test_windows,33780
gene_disjoint_assert,test genes and train genes are disjoint (asserted)
python,3.14.7
generated,2026-09-21 22:19:28 +0530
```

## 2. split_composition

SPLIT-B train and test windows (gene-disjoint by construction; asserted).

```csv
set,n_windows,n_genes,n_pathogenic
train,135147,86,62525
test,33780,14,15616
```

## 3. length_summary

Length statistics of the four sequence columns over all 168,927 variants. ref_sequence is the window used for similarity.

```csv
column,n,min_len,median_len,mean_len,max_len,n_len_100,n_with_N
ref_sequence,168927,100,100.000000,100.000000,100,168927,0
mut_sequence,168927,51,100.000000,99.710780,149,135105,9
genomic_ref_sequence,168927,100,100.000000,100.000000,100,168927,0
genomic_mut_sequence,168927,51,100.000000,99.710780,149,135105,9
```

## 4. length_distribution

Full length distribution of each sequence column.

```csv
column,length,n_windows
ref_sequence,100,168927
mut_sequence,51,4
mut_sequence,52,4
mut_sequence,53,4
mut_sequence,54,12
mut_sequence,55,10
mut_sequence,56,24
mut_sequence,57,9
mut_sequence,58,9
mut_sequence,59,18
mut_sequence,60,6
mut_sequence,61,8
mut_sequence,62,19
mut_sequence,63,16
mut_sequence,64,12
mut_sequence,65,25
mut_sequence,66,17
mut_sequence,67,7
mut_sequence,68,20
mut_sequence,69,19
mut_sequence,70,16
mut_sequence,71,39
mut_sequence,72,40
mut_sequence,73,30
mut_sequence,74,46
mut_sequence,75,55
mut_sequence,76,34
mut_sequence,77,59
mut_sequence,78,70
mut_sequence,79,49
mut_sequence,80,85
mut_sequence,81,103
mut_sequence,82,77
mut_sequence,83,124
mut_sequence,84,143
mut_sequence,85,97
mut_sequence,86,211
mut_sequence,87,231
mut_sequence,88,134
mut_sequence,89,314
mut_sequence,90,327
mut_sequence,91,178
mut_sequence,92,410
mut_sequence,93,486
mut_sequence,94,261
mut_sequence,95,850
mut_sequence,96,1564
mut_sequence,97,838
mut_sequence,98,3609
mut_sequence,99,11972
mut_sequence,100,135105
mut_sequence,101,6626
mut_sequence,102,1466
mut_sequence,103,270
mut_sequence,104,930
mut_sequence,105,357
mut_sequence,106,114
mut_sequence,107,187
mut_sequence,108,195
mut_sequence,109,54
mut_sequence,110,89
mut_sequence,111,69
mut_sequence,112,48
mut_sequence,113,67
mut_sequence,114,59
mut_sequence,115,24
mut_sequence,116,51
mut_sequence,117,59
mut_sequence,118,29
mut_sequence,119,49
mut_sequence,120,48
mut_sequence,121,22
mut_sequence,122,29
mut_sequence,123,34
mut_sequence,124,10
mut_sequence,125,24
mut_sequence,126,21
mut_sequence,127,21
mut_sequence,128,22
mut_sequence,129,17
mut_sequence,130,5
mut_sequence,131,13
mut_sequence,132,15
mut_sequence,133,3
mut_sequence,134,11
mut_sequence,135,10
mut_sequence,136,2
mut_sequence,137,8
mut_sequence,138,6
mut_sequence,139,5
mut_sequence,140,3
mut_sequence,141,4
mut_sequence,142,5
mut_sequence,143,10
mut_sequence,144,13
mut_sequence,145,4
mut_sequence,146,9
mut_sequence,147,3
mut_sequence,148,3
mut_sequence,149,4
genomic_ref_sequence,100,168927
genomic_mut_sequence,51,4
genomic_mut_sequence,52,4
genomic_mut_sequence,53,4
genomic_mut_sequence,54,12
genomic_mut_sequence,55,10
genomic_mut_sequence,56,24
genomic_mut_sequence,57,9
genomic_mut_sequence,58,9
genomic_mut_sequence,59,18
genomic_mut_sequence,60,6
genomic_mut_sequence,61,8
genomic_mut_sequence,62,19
genomic_mut_sequence,63,16
genomic_mut_sequence,64,12
genomic_mut_sequence,65,25
genomic_mut_sequence,66,17
genomic_mut_sequence,67,7
genomic_mut_sequence,68,20
genomic_mut_sequence,69,19
genomic_mut_sequence,70,16
genomic_mut_sequence,71,39
genomic_mut_sequence,72,40
genomic_mut_sequence,73,30
genomic_mut_sequence,74,46
genomic_mut_sequence,75,55
genomic_mut_sequence,76,34
genomic_mut_sequence,77,59
genomic_mut_sequence,78,70
genomic_mut_sequence,79,49
genomic_mut_sequence,80,85
genomic_mut_sequence,81,103
genomic_mut_sequence,82,77
genomic_mut_sequence,83,124
genomic_mut_sequence,84,143
genomic_mut_sequence,85,97
genomic_mut_sequence,86,211
genomic_mut_sequence,87,231
genomic_mut_sequence,88,134
genomic_mut_sequence,89,314
genomic_mut_sequence,90,327
genomic_mut_sequence,91,178
genomic_mut_sequence,92,410
genomic_mut_sequence,93,486
genomic_mut_sequence,94,261
genomic_mut_sequence,95,850
genomic_mut_sequence,96,1564
genomic_mut_sequence,97,838
genomic_mut_sequence,98,3609
genomic_mut_sequence,99,11972
genomic_mut_sequence,100,135105
genomic_mut_sequence,101,6626
genomic_mut_sequence,102,1466
genomic_mut_sequence,103,270
genomic_mut_sequence,104,930
genomic_mut_sequence,105,357
genomic_mut_sequence,106,114
genomic_mut_sequence,107,187
genomic_mut_sequence,108,195
genomic_mut_sequence,109,54
genomic_mut_sequence,110,89
genomic_mut_sequence,111,69
genomic_mut_sequence,112,48
genomic_mut_sequence,113,67
genomic_mut_sequence,114,59
genomic_mut_sequence,115,24
genomic_mut_sequence,116,51
genomic_mut_sequence,117,59
genomic_mut_sequence,118,29
genomic_mut_sequence,119,49
genomic_mut_sequence,120,48
genomic_mut_sequence,121,22
genomic_mut_sequence,122,29
genomic_mut_sequence,123,34
genomic_mut_sequence,124,10
genomic_mut_sequence,125,24
genomic_mut_sequence,126,21
genomic_mut_sequence,127,21
genomic_mut_sequence,128,22
genomic_mut_sequence,129,17
genomic_mut_sequence,130,5
genomic_mut_sequence,131,13
genomic_mut_sequence,132,15
genomic_mut_sequence,133,3
genomic_mut_sequence,134,11
genomic_mut_sequence,135,10
genomic_mut_sequence,136,2
genomic_mut_sequence,137,8
genomic_mut_sequence,138,6
genomic_mut_sequence,139,5
genomic_mut_sequence,140,3
genomic_mut_sequence,141,4
genomic_mut_sequence,142,5
genomic_mut_sequence,143,10
genomic_mut_sequence,144,13
genomic_mut_sequence,145,4
genomic_mut_sequence,146,9
genomic_mut_sequence,147,3
genomic_mut_sequence,148,3
genomic_mut_sequence,149,4
```

## 5. exact_duplicates

SPLIT-B test windows whose exact sequence (forward) or reverse complement also occurs among train windows.

```csv
window_column,n_test_windows,n_exact_in_train_forward,n_revcomp_in_train,n_either
ref_sequence,33780,0,0,0
mut_sequence,33780,0,0,0
```

## 6. similarity_distribution

Distribution of the per-test-window maximum shared 25-mer fraction, and the count / fraction of test windows at or above each threshold.

```csv
variant,statistic,value,fraction_of_scored_windows
forward,n_test_windows_scored,33780.000000,
forward,quantile_0,0.000000,
forward,quantile_0.01,0.000000,
forward,quantile_0.05,0.000000,
forward,quantile_0.25,0.000000,
forward,quantile_0.5,0.000000,
forward,quantile_0.75,0.000000,
forward,quantile_0.9,0.000000,
forward,quantile_0.95,0.000000,
forward,quantile_0.99,0.000000,
forward,quantile_1,0.342105,
forward,mean,0.000743,
forward,n_ge_0.5,0.000000,0.000000
forward,n_ge_0.8,0.000000,0.000000
forward,n_ge_0.95,0.000000,0.000000
canonical,n_test_windows_scored,33780.000000,
canonical,quantile_0,0.000000,
canonical,quantile_0.01,0.000000,
canonical,quantile_0.05,0.000000,
canonical,quantile_0.25,0.000000,
canonical,quantile_0.5,0.000000,
canonical,quantile_0.75,0.000000,
canonical,quantile_0.9,0.000000,
canonical,quantile_0.95,0.000000,
canonical,quantile_0.99,0.000000,
canonical,quantile_1,0.342105,
canonical,mean,0.000829,
canonical,n_ge_0.5,0.000000,0.000000
canonical,n_ge_0.8,0.000000,0.000000
canonical,n_ge_0.95,0.000000,0.000000
```

## 7. gene_pairs_high_similarity

Top 15 (test gene, train gene) pairs per variant and threshold by weighted hits; share_of_hits = n_hits_weighted / total_hits_at_threshold.

```csv
variant,threshold,test_gene,train_gene,n_hits_weighted,share_of_hits,total_hits_at_threshold,n_distinct_pairs_at_threshold
```

## 8. test_gene_hit_counts

Number of test windows at or above each threshold, per test gene.

```csv
variant,threshold,test_gene,n_test_windows_hit,share_of_hits,n_test_windows_in_gene
```

## 9. top_similarity_windows

Descriptive: the 20 test windows with the highest forward max shared 25-mer fraction, with their best train window and gene (not an additional threshold).

```csv
rank,test_vid,test_gene,label,max_shared_frac_forward,best_train_vid_forward,best_train_gene_forward,n_tied_train_windows_forward,max_shared_frac_canonical,best_train_gene_canonical
1,13_32325701_C_T,BRCA2,0,0.342105,19_1219703_G_A,STK11,1,0.342105,STK11
2,13_32325741_C_T,BRCA2,0,0.342105,19_1219703_G_A,STK11,1,0.342105,STK11
3,13_32325750_T_C,BRCA2,0,0.315789,19_1219703_G_A,STK11,1,0.315789,STK11
4,13_32325696_A_G,BRCA2,0,0.315789,19_1219703_G_A,STK11,1,0.315789,STK11
5,7_140784730_C_T,BRAF,0,0.302632,9_132907554_A_G,TSC1,3,0.302632,TSC1
6,13_32380686_A_C,BRCA2,0,0.302632,3_10146298_C_T,VHL,2,0.302632,LDLR
7,7_140787300_CAAAAAAAA_C,BRAF,0,0.276316,9_84888982_A_ATTTTTTTTTTTTT,NTRK2,1,0.276316,NTRK2
8,7_140787279_C_T,BRAF,0,0.276316,9_84888982_A_ATTTTTTTTTTTTT,NTRK2,1,0.276316,NTRK2
9,7_140787300_C_CAA,BRAF,0,0.276316,9_84888982_A_ATTTTTTTTTTTTT,NTRK2,1,0.276316,NTRK2
10,7_140787300_C_CA,BRAF,0,0.276316,9_84888982_A_ATTTTTTTTTTTTT,NTRK2,1,0.276316,NTRK2
11,7_140787300_CA_C,BRAF,0,0.276316,9_84888982_A_ATTTTTTTTTTTTT,NTRK2,1,0.276316,NTRK2
12,7_140787300_C_CAAA,BRAF,0,0.276316,9_84888982_A_ATTTTTTTTTTTTT,NTRK2,1,0.276316,NTRK2
13,11_108367112_C_T,ATM,0,0.263158,22_29696915_C_T,NF2,1,0.263158,NF2
14,X_33283879_A_G,DMD,0,0.250000,17_43057356_A_G,BRCA1,1,0.250000,BRCA1
15,11_108317650_T_C,ATM,0,0.250000,9_84870126_T_TATATATAC,NTRK2,1,0.238095,NTRK2
16,11_108317650_TATACACAC_T,ATM,0,0.250000,9_84870126_T_TATATATAC,NTRK2,1,0.238095,NTRK2
17,11_108317654_C_T,ATM,0,0.250000,9_84870126_T_TATATATAC,NTRK2,1,0.238095,NTRK2
18,X_33283874_T_C,DMD,0,0.250000,17_43057356_A_G,BRCA1,1,0.250000,BRCA1
19,22_23290924_G_A,BCR,0,0.250000,17_7674639_C_A,TP53,1,0.250000,TP53
20,X_33283876_G_A,DMD,0,0.250000,17_43057356_A_G,BRCA1,1,0.250000,BRCA1
```

## 10. not_applicable_or_skipped

Items in the generic reporting note that do not apply to D6.

```csv
item,status,note
af_coverage_and_bucket_cutoffs,not applicable,D3 items in the generic reporting note; D6 uses no allele frequency
models,not applicable,D6 scores no model; it compares sequences only
bucketing_or_skipped_groups,none,no test window was excluded; windows with no valid 25-mer would be unscored (count in similarity_distribution: n_test_windows_scored)
```

## 11. integrity_check

sha256 of every file under results/, scripts/, docs/ before vs. after (results/checkpoints: size + mtime only). Added files are this task's outputs.

```csv
check,value
files_hashed_sha256_before,309
files_hashed_sha256_now,321
modified_preexisting_files,none
removed_files,none
added_files_count,12
added_files,results/d6_cross_gene_similarity/d6_similarity_histogram.png;results/d6_cross_gene_similarity/exact_duplicates.csv;results/d6_cross_gene_similarity/gene_pairs_high_similarity.csv;results/d6_cross_gene_similarity/length_distribution.csv;results/d6_cross_gene_similarity/length_summary.csv;results/d6_cross_gene_similarity/parameters.csv;results/d6_cross_gene_similarity/similarity_distribution.csv;results/d6_cross_gene_similarity/similarity_per_test_window.csv;results/d6_cross_gene_similarity/test_gene_hit_counts.csv;results/d6_cross_gene_similarity/top_similarity_windows.csv;results/handoff/d6_export.md;scripts/run_d6_cross_gene_similarity.py
checkpoints_size_mtime_before,92
checkpoints_size_mtime_now,92
checkpoints_changed_added_removed,none
result,PASS: no pre-existing file changed
```
