# Bucket Eligibility for Matching (>= 50 of both classes)

Split,Bucket,Pathogenic,Benign,Used,Reason
A,missing,13495,13129,yes,-
A,0,14,16,skipped,smaller class has 14 (< 50)
A,"(0,1e-4)",1386,2209,yes,-
A,"[1e-4,1e-3)",636,1034,yes,-
A,"[1e-3,1e-2)",89,815,yes,-
A,"[1e-2,5e-2)",6,389,skipped,smaller class has 6 (< 50)
A,>=5e-2,2,566,skipped,smaller class has 2 (< 50)
B,missing,13746,13295,yes,-
B,0,12,21,skipped,smaller class has 12 (< 50)
B,"(0,1e-4)",1288,2174,yes,-
B,"[1e-4,1e-3)",519,1004,yes,-
B,"[1e-3,1e-2)",49,843,skipped,smaller class has 49 (< 50)
B,"[1e-2,5e-2)",1,330,skipped,smaller class has 1 (< 50)
B,>=5e-2,1,497,skipped,smaller class has 1 (< 50)
C,missing,9126,10241,yes,-
C,0,44,30,skipped,smaller class has 30 (< 50)
C,"(0,1e-4)",3712,3487,yes,-
C,"[1e-4,1e-3)",2130,2309,yes,-
C,"[1e-3,1e-2)",377,1264,yes,-
C,"[1e-2,5e-2)",20,372,skipped,smaller class has 20 (< 50)
C,>=5e-2,4,420,skipped,smaller class has 4 (< 50)


ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models.
