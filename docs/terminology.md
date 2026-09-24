# Terminology: leakage, shortcut dependence, distribution shift

These three ideas are kept apart in the writeup. None of the mechanisms below has been established by the results so far; each section says what a given experiment probes, not what it proves.

## (a) Strict data leakage

**Definition.** The same variant (chrom, pos, ref, alt), or a duplicated input sequence, appears in both train and test.

**Status: not evidenced.** In the 168,927-variant dataset every variant ID is unique. In SPLIT-A, SPLIT-B and SPLIT-C there are 0 variant IDs shared between train and validation, 0 identical (ref, mut) sequence pairs, and 0 identical mutant sequences across the two sides.

**Related, but not strict leakage by this definition: locus overlap.** Two variants at the same position with different alt alleles share an identical reference window. 36.1% of SPLIT-A validation variants and 34.3% of SPLIT-C validation variants have a training variant at the same locus; SPLIT-B has 0% (it is gene-disjoint). Whether this overlap affects any result has not been tested.

**Addressed by:** the frozen split manifests and the checks above. No experiment isolates locus overlap on its own (that would need a locus-disjoint split).

## (b) Shortcut / memorization dependence

**Definition.** Performance that comes from information that does not transfer to unseen genes, for example per-gene label priors or gene-specific patterns memorized from same-gene training variants, rather than variant-level signal.

**Addressed by:**
- SPLIT-A vs SPLIT-B (gene-disjoint): the A->B split-shift delta.
- M3 gene-identity-only baseline (D2): how much AUROC gene identity alone provides (SPLIT-A and SPLIT-C; SPLIT-B is a constant prior, AUROC 0.5 by construction).
- D4 gene-label permutation test: whether M3's SPLIT-A/C AUROC exceeds what shuffled gene labels give.

**Not established.** A drop from A to B is consistent with dependence on gene identity, but it is also confounded: the SPLIT-B validation set covers 14 genes against 100, gene sizes and label mix differ, and DNABERT-2's SPLIT-B run needed LR warmup. The DNABERT-2 multi-seed runs ask whether the SPLIT-B number is stable; they do not test shortcuts.

## (c) Distribution shift

**Definition.** Test variants come from a different distribution than training variants, here along time.

**Addressed by:** SPLIT-C (temporal): train = ClinVar LastEvaluated on or before the cutoff date (or undated, placed in train), validation = strictly later. Compared through the A->C split-shift delta, and through M3 on SPLIT-C.

**Caveat.** LastEvaluated is the date a record was last reviewed, not the date a variant first appeared. Later-dated variants therefore mix newly submitted variants with older variants that were re-evaluated (and possibly reclassified), so SPLIT-C conflates novelty, re-review and label change. Locus overlap (see above) also remains under SPLIT-C.

**Not established.** Which of these factors drives any A->C drop.

## Summary

| Concept | Experiments that address it | Status |
|---|---|---|
| (a) Strict data leakage | Manifest checks (ID / pair / sequence overlap) | Not evidenced; locus overlap present in A and C, effect untested |
| (b) Shortcut / memorization dependence | A vs B, M3 (D2), D4 permutation | Probed; mechanism not established |
| (c) Distribution shift | A vs C, M3 on C | Probed; conflated with re-evaluation (LastEvaluated caveat) |
