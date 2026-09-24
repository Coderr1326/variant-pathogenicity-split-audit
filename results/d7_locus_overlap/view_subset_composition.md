# Test-Set Subset Composition (locus-shared vs. locus-novel)

Split,Subset,n,n pathogenic,Pathogenic frac.,Genes,Top 1,Top 2,Top 3,Top 4,Top 5
A,full,33786,15628,0.4626,100,"BRCA2 (1,202)","DMD (1,183)","RYR1 (1,160)","BRCA1 (1,159)","TSC2 (1,158)"
A,shared,12186,6504,0.5337,93,MSH2 (736),MSH6 (697),TSC2 (512),BRCA1 (503),APC (452)
A,novel,21600,9124,0.4224,100,RYR1 (935),DMD (903),BRCA2 (855),ATM (799),NF1 (732)
B,full,33780,15616,0.4623,14,"ATM (5,804)","BRCA2 (5,804)","DMD (5,804)","NF1 (5,804)","CDH1 (2,673)"
C,full,33536,15413,0.4596,96,"TSC2 (1,735)","MSH2 (1,706)","ATM (1,287)","MSH6 (1,271)","MLH1 (1,189)"
C,shared,11495,5498,0.4783,91,MSH2 (937),MSH6 (715),TSC2 (640),MLH1 (539),APC (417)
C,novel,22041,9915,0.4498,96,"TSC2 (1,095)",MET (917),ATM (916),MSH2 (769),NF1 (732)


SPLIT-B is gene-disjoint: every test variant is locus-novel, so only the full set is shown. Subsets are not random samples, so differences are descriptive.
