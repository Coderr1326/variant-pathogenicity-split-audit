# Locus Overlap: Test Variants Sharing a Locus with a Training Variant

Split,Test variants,Locus-shared,Locus-novel,Fraction shared,Recorded in docs (%)
A,33786,12186,21600,0.3607,36.1
B,33780,0,33780,0.0000,0.0
C,33536,11495,22041,0.3428,34.3


locus = chrom_pos; shared = the same locus occurs among the split's training variants (any alt allele). Reproduces the 36.1% / 34.3% figures in docs/terminology.md.
