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

Order: **T6 → T1 → T7 → T3 → T5 → T4 + T8 → T2** (T6–T8, website problems found on 13-10, are
described after T2). Every decision rule below is fixed now,
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

## Problems found on the website (13-10): disease subtypes and broad terms

Found by looking at the breast cancer page (https://orevid.org/#/d/36691). The graph has the subtypes
(MONDO: triple-negative #37655, HER2 positive #84130, oestrogen-receptor negative #37773, HER2 negative
#37775, invasive lobular #94666, inflammatory #94693, male #94716…), each with its own page, but:

1. **No navigation.** Nothing on a page links to its more specific forms or to its broader term; the
   ontology is used by the model and the direct extensions, never shown.
2. **Subtypes defined by a biomarker inherit drugs that cannot work.** The graph does not encode receptor
   status, and these pages have few genes (29–35, from GWAS), so the model copies the parent's drugs:
   tamoxifen is candidate #1 for triple-negative breast cancer; tamoxifen, formestane and testolactone
   (endocrine therapies) are candidates for ER-negative; for HER2-negative, four corticosteroids arrive
   through «non-Hodgkin lymphoma» and «allergic rhinitis» with 3–4 shared genes.
3. **Broad terms pass the shared-genes test with any cancer.** Breast cancer has 1,253 associated genes;
   it shares 1,085 with hereditary breast–ovarian cancer syndrome and 257 with non-small cell lung cancer,
   so p < 1e-6 is trivial. Its candidates (platinum drugs, alkylating agents, PARP inhibitors,
   pembrolizumab) are plausible only for subtypes (pembrolizumab is approved for triple-negative breast
   cancer since 2021; PARP inhibitors for BRCA-mutated). Also, «Approved (DrugBank)» lists only 2 drugs
   because PrimeKG spreads the indications over near-synonym nodes (breast cancer / breast carcinoma /
   breast neoplasm), and epirubicin and cyclophosphamide, approved for «breast carcinoma», are shown as
   direct extensions, i.e. as hypotheses.

Stated in the preprint's Limitations on 13-10 (two sentences after «comorbidities that read as shared
biology»). Remove or update them when T7 and T8 are done.

### T6. Ontology navigation and a «broad term» notice (presentation only) — ~2 h

No validation needed (display; DECISIONS rule 3 is met).

1. **Hierarchy with direction.** PrimeKG's `disease_disease` edges are undirected. Take parent → child
   from the Open Targets disease table (`parents`, already read in `sources/ids.py`, MONDO ids) mapped to
   graph nodes through each node's MONDO ids; keep only pairs where both ends are nodes. If PrimeKG
   groups several MONDO ids in one node, a parent inside the same node is ignored. Save it once in
   `data/prep/` (e.g. `disease_hierarchy.parquet`, columns parent, child) from `prepare.py`, or compute it
   in `web/build.py` if it is cheap.
2. **Page data.** In `web/build.py` `DiseasePages.disease()` add `broader` (direct parents) and
   `narrower` (direct children), each as `[node, name, has_page]`, children sorted by number of associated
   genes, plus `n_descendants`. Same output in `explore.py disease` (it uses the same class).
3. **Frontend** (`web/public/app.js`): under the header, «Broader term: …» and «More specific forms: …»
   as links (collapsed beyond ~12 children). A notice when `n_descendants ≥ 5` or `n_genes ≥ 500`:
   «This is a broad term. Candidates may apply only to one of its forms; see the more specific forms
   below.» All strings in `web/public/i18n/{en,es,ca}.json`; bump `?v=N` of `app.js` in `index.html`.
4. **Near-synonyms.** On the extensions table, mark a drug as «approved for an equivalent term» when it
   is indicated for a parent or child that has (almost) the same name tokens or a gene-set Jaccard
   ≥ 0.8 with this disease (breast cancer / breast carcinoma / breast neoplasm). Display only: it does not
   change the ranking or the evaluation.
5. **Tests** (`tests/`): breast cancer (36691) lists 37655 and 84130 among its descendants; every child
   lists its parent as broader; no cycles.

### T7. Biomarker-defined subtypes: marker rule (needs validation) — ~3 h

1. **Curated table** `resources/biomarker_subtypes.csv`, short and reviewed by hand: `node, name, marker,
   status (absent|present), drug_rule, source`. Start with breast cancer:
   - triple-negative (37655): ESR1, PGR, ERBB2 absent
   - oestrogen-receptor negative (37773): ESR1 absent; progesterone-receptor negative (37774): PGR absent
   - HER2 negative (37775): ERBB2 absent; HER2 positive (84130): ERBB2 present
   
   Then search the index for other marker-defined subtypes (names with «positive», «negative»,
   «-mutant», «ALK», «EGFR», «BRAF», «KRAS», «BCR-ABL1», «receptor») and add the clear ones.
2. **Drug rule per absent marker.** ESR1 absent → endocrine therapy: the drug targets ESR1, ESR2 or
   CYP19A1 (aromatase), or its ATC class starts with L02 (hormone antagonists and related agents; covers
   testolactone and androgens used in breast cancer). PGR absent → targets PGR. ERBB2 absent → targets
   ERBB2. For «present» markers no exclusion, only an annotation.
3. **Apply it** in `DiseasePages.disease()` as an exclusion with a reason, like the contraindicated
   drugs (P24): the drug goes to `excluded` with «the marker this drug acts on is absent in this subtype».
4. **Validation, fixed before running** (too few diseases for a hit-rate test):
   (a) list every flagged (drug, disease) pair over all 10,042 pages and check them by hand;
   (b) none of the flagged pairs may be an approved indication (DrugBank or confirmed Open Targets) for
   that subtype; if any is, the rule or the table is wrong. Flagged pairs with a phase 3 trial in that
   subtype are listed for the manual check;
   (c) report how many candidates and extensions it removes. Record as a new DECISIONS entry.
5. **Tests:** tamoxifen, formestane and testolactone are not candidates for 37655 or 37773; trastuzumab
   stays approved on 84130.

**Preprint.** Add the rule to Methods (Evidence rules, the paragraph on exclusions) and replace the
tamoxifen sentence in Limitations by what remains (markers not in the table; subtypes defined by
mutations the graph does not encode).

### T8. Shared genes in broad diseases (with T4) — ~1 h on top of T4

1. In the T4 grid, add an **effect-size** variant: fold enrichment (observed shared genes / expected
   under the hypergeometric) ≥ 2, besides p < 1e-6; and the existing «≥ 3 specific genes» annotation
   (P29) as a filter variant.
2. Split the hit rates of the shared-genes section by the disease's gene-set size
   (< 100, 100–499, ≥ 500 genes) on the benchmark, the areas and the external exam, with the T3
   disease-level bootstrap.

**Decision rule.** If, for diseases with ≥ 500 genes, the shared-genes hit rate on the external exam is
not above the hidden section (1.4 %) by more than its CI, broad pages stop showing shared-genes
candidates and point to the more specific forms (T6 notice) instead. A filter variant (fold enrichment or
specific genes) is adopted only under the P29 rule: the subset it removes hits < 5 % on both the benchmark
and the areas, and the external hit rate of what is kept does not fall.

**Preprint.** Methods (Evidence rules) if a rule changes; Results (Hit rates by evidence line) with the
split by gene-set size; the «broad terms» sentence in Limitations.

## After the pending tasks

1. Update the Claude Docs preprint (link above); keep `docs/preprint/` in sync: export markdown and docx
   from the doc, replace the four `&#91;embedded content: …\]` lines of the markdown export with
   `![](figures/figureN.png)` in order, keep the first comment line.
2. DECISIONS.md: close P34 (or add P35) with the results and decisions.
3. If T5 changes the website figures: update the Method page and redeploy (`web/deploy.sh`, see STATUS).

## Results of the pending tasks (home PC, 13–14 October 2026)

All decision rules above were applied as written. Raw outputs: `data/eval/` (not in git); generated reports:
`docs/EVIDENCE_CI.md` (T3, T4, T5, T8), `docs/BIOMARKERS.md` (T7), `docs/EXTERNAL_SEEDS.md` (T2).

### T6. Navigation and broad-term notice — done

`sources/hierarchy.py` builds `data/prep/disease_hierarchy.parquet` from the Open Targets disease table
(26,998 parent → child pairs between 17,615 disease nodes; edges inside cycles, created when several terms map
to one node, are dropped). Pages carry `broader`, `narrower` (children by number of genes) and
`n_descendants`; the frontend shows them as links, collapses beyond 12 children and shows the broad-term notice
when there are ≥ 5 descendants or ≥ 500 genes. Direct extensions indicated for a near-synonym (same name without
generic words, or gene-set Jaccard ≥ 0.8) are marked «approved for an equivalent term» (breast cancer:
epirubicin and cyclophosphamide, approved for «breast neoplasm»). Tests: `tests/test_hierarchy.py`.

### T1. External exam without Open Targets features — no material dependence

| Feature set | External Recall@100 | Zero-shot (5 seeds) | Zero-shot, hard | Areas |
|---|---|---|---|---|
| Production (`extra_no_ctd`) | 0.359 | 0.811 ± 0.042 | 0.539 | 0.617 |
| No Open Targets (`extra_no_ot`: PrimeKG + STRING) | 0.348 | 0.802 ± 0.048 | 0.522 | 0.628 |
| No external source (`new`: PrimeKG only) | 0.351 | — | — | — |
| TxGNN | 0.121 | 0.412 | 0.258 | 0.237 |

Rule: drop < 0.05 and margin over TxGNN > 0.10 → met (drop 0.011, margin 0.227).

### T3. Hit rates with disease-level intervals — done

The production rule reproduces the published counts exactly. Intervals from resampling whole diseases are two
to three times wider than Wilson's: external exam, direct extensions 12.6 % (9.8–15.6), shared genes 4.2 %
(3.1–5.5), hidden 1.4 % (1.1–1.8); benchmark shared genes 22.7 % (17.2–28.4); areas 25.2 % (20.8–29.8).
Full table: `docs/EVIDENCE_CI.md`.

### T5. Website rules applied exactly on the external exam — figures confirmed

With the website's specificity (2,000 reference diseases of the seed-1 training graph, model checked to
reproduce the published scores) and its exclusions (same active moiety, contraindicated, diagnostic and
veterinary drugs; Open Targets approvals are not excluded here because they are the answer): direct extensions
12.4 % (9.7–15.6), shared genes 4.3 % (3.1–5.6), hidden 1.5 % (1.1–1.9), all within the intervals of the
approximation, so the published figures stand. Check C, measured apart on the external exam: the model's top 25
without any rule hit 4.2 % (3.5–5.0); the rules without the model (every candidate with a line, any rank)
0.8 % (0.6–1.0); every candidate 0.1 %. Most of the enrichment comes from the model's ranking.

### T7. Biomarker rule — adopted after the manual check

`resources/biomarker_subtypes.csv` (16 rows: breast ER/PR/HER2 subtypes, CML BCR-ABL1, ALCL ALK); an absent
marker also applies to the descendants. Drug rule for an absent oestrogen receptor: ATC «endocrine therapy»
(L02) except anti-androgens, or targets ESR1, ESR2 or CYP19A1 unless it is a cytotoxic antineoplastic (L01).
The two exceptions were found by the pre-specified manual check (criterion a): enzalutamide (an anti-androgen,
studied in androgen-receptor-positive triple-negative breast cancer) and estramustine (a cytotoxic that acts on
microtubules) were flagged by the first version and are not endocrine therapy in this sense. Criterion (b): 0 of
the flagged pairs is an approved indication; the three with a phase 3 trial (lapatinib in triple-negative,
tamoxifen and letrozole in ER-negative breast cancer) are correct exclusions (no benefit without the receptor).
Criterion (c): 9 candidates and 0 direct extensions removed (tamoxifen from triple-negative; formestane,
tamoxifen and testolactone from ER-negative breast cancer; five BCR-ABL inhibitors from atypical CML,
BCR-ABL1 negative). Every pair: `docs/BIOMARKERS.md`. Tests: `tests/test_biomarkers.py`.

### T4. Shared-genes thresholds and permutation null — (a) met, (b) failed

(a) No sharp optimum: on the external exam the 12 cells of the grid (2–5 genes × p < 1e-4, 1e-6, 1e-8) give
3.9–4.6 % against 4.2 % for the production rule; on the benchmark 20.7–24.3 % (22.7 %), on the areas 22.5–27.8 %
(25.2 %). (b) With degree-preserving random gene sets (every disease keeps its size; genes drawn in proportion to
their number of diseases; ranking and indications unchanged; 100 permutations) the shared-genes section hits
4.5 % on average (99th percentile 4.8 %) against 4.2 % observed. **By the rule fixed in advance, the shared-genes
line cannot be presented as biological evidence.** Exploratory analyses added afterwards (not pre-specified,
`py -3 -m evaluation.evidence_ci null`) confirm it within every gene-set-size bin on the external exam (2.7 vs
2.2 %, 5.7 vs 5.5 %, 4.6 vs 4.5 %, all below the null's 99th percentile) and with ≥ 3 specific genes (4.8 vs
5.1 %); on the benchmark only diseases with < 100 genes are above the null (14.8 % vs 99th percentile 14.0 %).
The line selects drugs already used in diseases of similar size and research intensity; the specific genes it
lists add nothing measurable. Consequence (P38): candidates and rule unchanged (the hit rate is measured), but the
website and the preprint no longer describe the shared genes as evidence of a shared mechanism.

### T8. Shared genes in broad diseases — broad pages keep the section

External exam, shared-genes section by the disease's gene-set size: < 100 genes 2.7 % (1.7–3.9), 100–499 5.7 %
(3.3–8.4), ≥ 500 4.6 % (2.4–7.4). The ≥ 500 bin is above the hidden section's 1.4 % by more than its interval, so
by the rule broad pages keep shared-genes candidates (with the T6 notice). Fold enrichment ≥ 2 removes nothing
(p < 1e-6 already implies it). «≥ 3 specific genes» as a filter fails the P29 rule (the subset it removes hits
9.0 % on the benchmark and 9.8 % on the areas) and stays a display-only annotation.

### T2. External exam over five training seeds — done, OREVID and TxGNN

`py -3 -m evaluation.external_exam evaluate --seed N` (OREVID, seeds 2–5) and `evaluation/txgnn_external_seeds.sh`
(TxGNN inference with the models already trained on each zero-shot split, ~10 min per seed in WSL);
summary with `py -3 -m evaluation.external_exam seeds` → `docs/EXTERNAL_SEEDS.md`. Recall@100, mean ± SD over
the five training splits:

| Model | All | Non-obvious | Recent |
|---|---|---|---|
| Popularity | 0.130 ± 0.010 | 0.076 ± 0.006 | 0.098 ± 0.012 |
| TxGNN | 0.126 ± 0.004 | 0.086 ± 0.004 | 0.138 ± 0.024 |
| OREVID | 0.367 ± 0.008 | 0.184 ± 0.007 | 0.365 ± 0.021 |
| Logistic regression without popularity | 0.363 ± 0.010 | 0.183 ± 0.013 | 0.339 ± 0.016 |

Seed 1 (the published exam) is representative. Preprint: Table 1 external columns as mean ± SD, abstract and
Figure 2 caption with the five-split values, «one training split» removed from Statistical analysis and
Limitations.

### Preprint and website after T1–T8 (14-10)

Preprint (Claude Doc, re-exported to `docs/preprint/`): abstract (Open Targets ablation, disease-level CIs, the
null), Conclusions, Methods (biomarker rule, exact website rules, permutation test), Results (Open Targets
ablation, null and threshold grid, specific genes), new Supplementary Table S5, Table S3 and Figure 4 with
disease-level intervals, Discussion and Limitations. Website: navigation, broad-term notice, biomarker
exclusions, and the shared genes described as context rather than evidence (candidates help and Method page,
three languages), deployed with `web.build --patch`.
