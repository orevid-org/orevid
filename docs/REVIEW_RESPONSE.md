# Response to the external review of the preprint (13-10-2026)

An AI agent reviewed the preprint draft as a journal referee (summary below). This file records, point
by point, what the data answer, what was done on 13-10 (work PC, needs only the release prediction
files) and **what is left for the home PC**, which has `data/`. Decision log entry: P34 in
`docs/DECISIONS.md`.

## Points raised and status

| # | Criticism (reviewer's severity) | Verdict | Status |
|---|---|---|---|
| 2.1 | External set not independent: used for the specificity weight and the CTD removal (high) | Mostly answered with existing numbers: the version frozen before the external set was consulted (CTD kept) gives 0.352 vs TxGNN 0.121 (P20); the specificity correction is not applied to benchmark rankings; no evidence rule was chosen on it; the prospective snapshot is already registered | **Done** (sentence in Evidence rules) |
| 2.2 | «Approval absent from PrimeKG» ≠ new discovery (high) | Correct, cannot be rebutted (median first trial 2008) | **Done**: wording narrowed in abstract and Discussion |
| 2.3 | Superiority only over TxGNN and popularity; TxGNN wins AUPRC in 5/9 areas (medium) | Already stated | **Done**: Discussion limited to the metrics and methods tested |
| 2.4 | CIs too narrow (clustering by disease), few splits, no multiplicity correction (medium-high) | Correct; fixable | Recall@100: **done** (`benchmark/CLUSTER_CI.md`, Table S4). Hit rates: **pending T3**. External exam on 5 seeds: **pending T2** |
| 2.5 | 3-gene / p < 1e-6 threshold: sensitivity and permutation null (medium-high) | Partial evidence exists (P21, P28, P29) | **Pending T4** |
| 2.6 | Open Targets supplies features and external labels (medium) | Design excludes OT drug and literature evidence from features, but no ablation | Stated in Limitations; **pending T1** |
| 3 / C | Hit rates are not the website's (no exclusions, specificity over test diseases) | Already stated in Methods | **Pending T5** |
| A | Tie handling in Recall@100 | Answered: ≤ 0.014 over 20 seeds, OREVID 0 | **Done** |
| B | Leakage through the masks | Reviewed before (A9, A11) | No action |
| 5 | Embedding baselines, strict temporal split, blinded expert review, experimental validation | Out of scope for now; stated | No action |

### Results obtained on 13-10 (`benchmark/cluster_ci.py`)

Disease-level bootstrap (10,000 resamples, within seeds / areas), paired by disease, Wilcoxon + Holm over
28 comparisons. OREVID − TxGNN Recall@100: zero-shot +0.399 (0.352–0.443), areas +0.380 (0.314–0.447),
external +0.238 (0.206–0.270), non-obvious +0.094 (0.064–0.123), recent +0.236 (0.141–0.332). OREVID ahead on
63 % of zero-shot diseases, behind on 2 %. TxGNN vs popularity: not significant after Holm except all area
diseases (p = 0.047). Full tables: `benchmark/CLUSTER_CI.md`.

Rerun: unzip `orevid-bench-predictions-v1.0.zip` (release) into `benchmark/`, then
`cd benchmark && py -3 cluster_ci.py` (~7 min; needs pandas + pyarrow).

## Preprint changes already made

In the Claude Docs document (source of truth, https://claude.ai/artifact/F5Rjejhgi21krq44CopKmf) and
re-exported to `docs/preprint/`: abstract Results (paired CI) and Conclusions («new uses» narrowed);
Evidence rules (0.352 before the external set was consulted); Benchmark (tie-breaking); Statistical
analysis (bootstrap, Wilcoxon + Holm, primary vs exploratory); Results (paired CIs, TxGNN vs popularity);
Discussion (metrics/methods scope, what the external exam measures); Limitations (Recall@100 CIs,
Open Targets overlap); new Supplementary Table S4.

## Pending on the home PC

Before starting: in `C:\dev\orevid` run `git pull`, then copy to `C:\dev\farmagraf` the files changed in
the 13-10 commit (`benchmark/cluster_ci.py`, `benchmark/CLUSTER_CI.md`, `benchmark/cluster_ci.json`,
`benchmark/README.md`, `README.md`, `docs/*.md`, `docs/preprint/*`), or switch to working in the public
clone only. Check that `data/eval/external_scores.npz` and `data/eval/external_positives.parquet` exist
(else `py -3 -m evaluation.external_exam evaluate`, a few minutes).

Order by value for the preprint: **T1 → T3 → T5 → T4 → T2**. Every decision rule below is fixed now,
before running; record the outcome whatever it is (P34 follow-up in DECISIONS.md).

### T1. External exam without Open Targets features (answers 2.6) — ~30 min

Open Targets enters the features through three metapaths of `PATHS_EXTRA` in `evaluation/txgnn.py`:
`target+ all` (ChEMBL mechanisms `moa` + OT gene sets), `interaction STRING all` and `similar all (OT)`
(both use `GENES_EXTRA`, which includes the `ot_assoc:*` gene sets).

1. In `evaluation/txgnn.py` add two metapaths that use PrimeKG genes only, and a feature set:
   ```python
   METAPATHS.update({
       "interaction STRING pkg": ("ppi", "drug_protein", "disease_protein",
                                  "target ↔ (STRING ≥ 700) PrimeKG disease gene", "string"),
   })
   FEATURE_SETS["extra_no_ot"] = PATHS + NEW_PATHS + ["interaction STRING pkg"]
   ```
   (`similar` and `target (all)` already cover the PrimeKG-only versions of the other two.) The name must
   start with `extra` so that `external_exam.score_ours` loads STRING; add `"extra_no_ot"` to the tuple in
   `txgnn.evaluate` (`extra=(feature_set in (...))`) too.
2. Run `py -3 -m evaluation.external_exam evaluate extra_no_ot` (ablation: prints, does not save) and,
   as a second reference, `py -3 -m evaluation.external_exam evaluate new` (no extra source at all).
3. Optional, same day: `py -3 -m evaluation.txgnn --paths extra_no_ot --rounds 5 --rank` (zero-shot) and
   `--split areas --seeds 1` to see the cost on the benchmark.

**Decision rule.** No material dependence if OREVID's external Recall@100 without Open Targets drops by
less than 0.05 (from 0.359) and stays above TxGNN by more than 0.10. Otherwise the preprint must say how
much of the external result depends on Open Targets features.

**Preprint.** Replace in Limitations «but a run without any Open Targets feature has not yet been done»
with the result; add a line to Table 1 (or a sentence in Results after the logistic regression).

### T3. Hit rates with a disease-level bootstrap (answers 2.4 for Table S3 / Figure 4) — ~1 h

`evaluation/evidence.py` sums hits over candidates (`res[k][0] += ...`). Keep the per-disease counts
instead: in `evaluate()` and `external()`, store `per_disease[(subset, section)][d] = [hits, n]`, then
bootstrap diseases (10,000 resamples; within seeds for the benchmark and within areas for the areas, as in
`benchmark/cluster_ci.py`) and report the ratio sum(hits)/sum(n) with its 2.5–97.5 percentiles. Write
the result to `docs/EVIDENCE.md` next to the Wilson intervals.

**Preprint.** Table S3: replace the Wilson intervals with the cluster intervals (or add them); abstract
«12.6% (95% CI 11.5–13.8%)» and «4.2% (3.6–4.9%)»; Figure 4 caption («Wilson intervals ignore the
clustering…»); Statistical analysis and Limitations sentences on Wilson intervals; the 4.8 % / 2.7 %
specific / generic genes in Results.

### T5. Website rules applied exactly on the external exam (answers section 3 and check C) — ~1–2 h

The external hit rates differ from the website in two ways (Methods, last paragraph of Benchmark and
evaluation): specificity over each split's test diseases instead of the 2,000 reference diseases, and no
exclusions. In `evidence.external()`:

1. Specificity: subtract 0.5 × each drug's mean standardised score over the 2,000 reference diseases of
   the zero-shot seed-1 training graph (same procedure as `predictor.py reference`), instead of
   `_transform` over the test diseases.
2. Exclusions as `web/build.py` `DiseasePages.disease().excluded`: same active moiety approved
   (`evaluation.twins.drug_family`), contraindicated for the disease (by name family too), diagnostic
   (`predictor.diagnostic_drugs`), veterinary only (`predictor.drug_status`).
3. Report the three sections and, separately, model-only (top 25 without rules) and rules-only numbers,
   so model, rules and full system are measured apart (check C).

**Decision rule.** None; this replaces an approximation by the exact figure. If the website numbers differ
from 12.6 % / 4.2 % by more than their CIs, the website's Method page and the abstract use the new ones.

**Preprint.** Remove the «with two differences» sentence from Benchmark and evaluation; update the
external row of Table S3, the abstract and Discussion («one in eight / one in 25»).

### T4. Sensitivity of the shared-genes thresholds and a permutation null (answers 2.5) — ~2–3 h

1. **Grid.** In `evidence.evaluate()` / `external()` make the rule's two parameters arguments
   (minimum shared genes, currently `len(com) >= 3`; p-value, currently `< 1e-6`) and run the grid
   genes ∈ {2, 3, 4, 5} × p ∈ {1e-4, 1e-6, 1e-8} on the benchmark (seeds 1–5), the areas and the
   external exam. Report hit rate and number of candidates shown per cell.
2. **Permutation null.** Keep the ranking and each drug's indications; replace the gene set of every
   disease with a random set of the same size, drawing genes with probability proportional to their
   degree (preserves disease size and gene popularity). 100 permutations; recompute the shared-genes
   section and its hit rate each time.

**Decision rule.** The 3 / 1e-6 rule stays if (a) its hit rate is within the range of the neighbouring
cells (no sharp optimum at the chosen point) and (b) the real hit rate is above the 99th percentile of the
permutation null on the external exam. If (a) fails, say that the threshold is tuned; if (b) fails, the
shared-genes line cannot be presented as biological evidence.

**Preprint.** New Supplementary Table S5 (grid) and one sentence in Results (Hit rates by evidence
line) with the null percentile; remove «0–2.5%» only if superseded.

### T2. External exam on five seeds (answers 2.4, «one split only») — OREVID ~15 min, TxGNN hours

`evaluation/external_exam.py` hard-codes `SPLIT, SEED = "complex_disease", 1`. Make the seed an
argument, save to `external_scores_s{seed}.npz` (keep seed 1 as the default file names so nothing else
breaks) and run seeds 2–5 for OREVID and popularity. TxGNN needs inference with the trained models of
seeds 2–5 in WSL (`txgnn_wsl.py --diseases data/eval/external_diseases.csv ...`, see the docstring); only
if those checkpoints were kept, otherwise skip TxGNN and report OREVID's spread alone.

**Preprint.** Table 1 external columns as mean ± SD over seeds (or the seed-1 value plus the range);
drop «it also uses one training split» from Limitations.

## After the pending tasks

1. Update the Claude Docs preprint (link above); keep `docs/preprint/` in sync: export markdown and docx
   from the doc, replace the four `&#91;embedded content: …\]` lines of the markdown export with
   `![](figures/figureN.png)` in order, keep the first comment line.
2. DECISIONS.md: close P34 (or add P35) with the results and decisions.
3. If T5 changes the website figures: update the Method page and redeploy (`web/deploy.sh`, see STATUS).
