# Results

Every number in this repository is generated from raw validation predictions — no paper values are ever copied in as results. Metrics JSONs are the source of truth; tables and figures are derived from them.

## Layout

| Directory | Contents |
|---|---|
| `metrics/` | One JSON per experiment: 6 models × 3 splits, transformer finals, 4-seed DNABERT-2 SPLIT-B runs (incl. warmup-0 control). Schema: `{model, window, metrics: {accuracy, roc_auc, weighted_*, per-class P/R/F1, confusion counts}}` |
| `tables/` | CSV (full precision) + MD + PNG renderings: cross-split comparison, split-shift deltas, local summaries, M3-vs-transformers, paper-vs-reproduction, gene counts |
| `figures/` | 29 publication-style figures: confusion matrices (all models × splits), model comparisons, training curves, permutation null |
| `reports/` | Narrative reports: paper specification extraction, reproduction report (fidelity self-assessment), discrepancy report, consistency audit |
| `m3_gene_only/` | Gene-identity-only baseline (smoothed per-gene rate) per split + D4 permutation null (200 gene-scramble permutations, seeds 0–199) |
| `d3_af_matched/` | Allele-frequency-matched evaluation: model AUROC on AF-balanced resampled subsets (100 resamples) + M4 AF-only baseline |
| `d5_per_gene/` | Per-gene decomposition: within-gene / between-gene / pooled AUROC, per-gene metrics with bootstrap CIs, heatmap + A-vs-B scatter |
| `d6_cross_gene_similarity/` | Exact + reverse-complement duplicate scan and 25-mer shingle similarity between SPLIT-B test genes and training genes |
| `d7_locus_overlap/` | Locus-shared vs locus-novel AUROC per model (bootstrap CIs), share of A→B drop attributable to locus overlap |
| `handoff/` | Full-precision markdown exports of every diagnostic (6-decimal tables) — the canonical numbers |
| `checkpoints/` | *(gitignored — 151 GB)* model weights, kept locally |

## Provenance rules

1. **Metrics JSONs are generated programmatically** from validation predictions; confusion matrices are asserted to sum to the evaluated label counts.
2. **Sanity-check runs are labelled as such** — early one-epoch/tiny metrics must not be read as final results (see reproduction report §J).
3. **The seed stopping rule was pre-registered** (commit `ff55707`, `docs/seed_stopping_rule.md`) before any seed run; the summarizer reports the rule bucket mechanically.
4. **No inference was re-run for diagnostics** — D3–D7 reuse saved prediction files; runs without saved probabilities are listed and skipped, not silently imputed.
5. **Split manifests are frozen and hash-verified** (`data/splits/*.json`); SPLIT-A's manifest is documentation-only (the training script recomputes it inline and the scripts assert equality).
