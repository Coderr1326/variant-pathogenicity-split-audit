# M4: AF-Only Baseline (pathogenic rate per AF bucket, train only, alpha = 1)

Split,Accuracy,Weighted F1,ROC-AUC,Pathogenic recall,M3 ROC-AUC (ref.)
A,0.5483,0.5090,0.5762,0.8635,0.6903
B,0.5511,0.5082,0.5796,0.8803,0.5000
C,0.5072,0.5050,0.5337,0.5921,0.6044


Threshold 0.5. M4 uses no sequence and no gene. ClinVar benign assertions often use high allele frequency as evidence, so AF is partly how the labels were defined; a drop under matching is expected and does not by itself indicate a shortcut. The informative quantity is how the drop differs between models.
