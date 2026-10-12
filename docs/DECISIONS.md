# Decision log and validation status

To keep track of what is solid and what is not before putting anything in front of professionals. Statuses:

- **Validated**: checked with out-of-sample data or reproduced independently.
- **Partial**: there is evidence, but with a known (written-down) gap.
- **Not validated**: reasonable but without a test. Not presented as fact.
- **Error fixed**: it was done wrong and has been fixed (kept for the record).

Rules (since 11-10):
1. No rule or production change is adopted without validating it on all diseases.
2. Parameter choices are made with a criterion fixed beforehand and confirmed out of sample.
3. Nothing is published with a «not validated» item visible to the end user (or it is explained there as a limitation).
4. Independent review (an agent with no context that looks for errors) before every publication.
5. Automated tests (`py -3 -m pytest tests`) for the critical invariants.

## Evaluation (comparison with TxGNN and OREVID-bench)

| # | Decision | Evidence | Status |
|---|---|---|---|
| A1 | Reproduce TxGNN's zero-shot split | TxGNN's official code generates the same split; reviewer: identical splits | Validated |
| A2 | Disease areas with their own split (fixed seed), the same for both | TxGNN sets no seed | Partial (1 seed; diabetes 2 and adrenal 6 diseases) |
| A3 | Metrics from TxGNN's notebook | TxGNN here ≈ published; independent OREVID-bench reproduces the figures; reviewer recomputed them | Validated |
| A4 | Popularity baseline | Training graph only; reviewer recomputed it | Validated |
| A5 | «Hard» subset | Explicit definition (ontology ≤ 2 steps, ≥ 50 % of drugs) | Partial (arbitrary definition) |
| A6 | Improvements designed after seeing the areas | `base` was already winning; independent external exam | Partial (risk of overfitting to the areas) |
| A7 | LambdaRank | Decided with seeds 1-2, confirmed with 3-5 and areas | Validated |
| A8 | Popularity «minus its own» in training rows | Reviewer: not active at test time, cannot inflate | Validated |
| A9 | Extra sources without drug-disease relations | Reviewer: no EXTRA matrix has drug-disease edges | Validated |
| A10 | External exam (Open Targets not in PrimeKG) | Mostly old (median 2008) | Partial (not prospective) |
| A11 | Metapath «ontology 2-step» | It included the same disease (d → neighbour → d). It did not inflate the test but did inflate training | **Error fixed** (review); everything re-evaluated |
| A12 | TxGNN gives zeros to diseases with no edge at all in the training graph (1-4 per split) | Reviewer; ~0.02 of the difference in Recall@100 | Partial: **must be stated** in the benchmark README |
| A13 | Our model uses the validation labels for training; TxGNN only for early stopping | Legitimate, but asymmetric | Partial: **must be stated** |
| A14 | CTD drug-gene comes from the literature (like the OT literature excluded for targets) | No ablation test | **Not validated**: run an ablation in the external exam |

## Production (`explore.py disease`, website)

| # | Decision | Evidence | Status |
|---|---|---|---|
| P1 | Production model = LambdaRank with the whole graph | Same procedure as the evaluated one | Validated (by procedure) |
| P2 | Score shown | It is not a probability; the website says so | Validated (presentation) |
| P3 | Penalties (failed trials, literature, genetic direction) | Calibrated with the old model | **Removed from the ranking**; shown as warnings |
| P4 | Drugs indicated for forms of the disease, listed separately | Presentation decision | Validated (presentation; they are not hidden) |
| P5 | «No prior use» list | Test C: Recall@100 ~0.15-0.18 | Partial: **the website states it with figures** |
| P6 | Labels from the direction model | Reviewer: without labelled neighbours, AUROC 0.447 | **Error fixed**: removed |
| P7 | Checkpoint-blocker rule | First version flagged lymphomas; now only autoimmune diseases (123) | Error fixed; validated against the data |
| P8 | «Approved or investigational» filter | An escaping error left it empty | Error fixed; automated test |
| P9 | Targets: full model | Recall@20 0.37 strict (0.59 before, with drug leaks into the candidates) | Validated with the strict figure |
| P10 | Targets: essentiality/tissue/ChEMBL as annotations | They do not improve the ranking | Validated as annotation |
| P11 | Licences of the website data | CTD removed (P20). But PrimeKG includes **DrugBank (CC BY-NC 4.0)** data: indications, targets, names, status. Compatible with a free non-profit website, with attribution | **Decided (12-10)**: non-commercial use accepted by the project owner. Results under **CC BY-NC-SA 4.0** because PrimeKG also carries SIDER and DisGeNET (BY-NC-SA); ChEMBL comes via Open Targets (CC0). Code MIT. Legal notice on the website (`#/legal`); owner and contact still missing |
| P12 | Quality panel (60 random diseases + 24 test diseases) | `docs/PANEL.md`: found generic candidates and candidates without explanation | Done; repeated at every change |
| P13 | Ranking by specificity (z − 0.5 × the drug's mean over 2,000 diseases) | Weight chosen with seeds 1-2 (criterion: max. reduction of generic drugs losing ≤ 0.01); confirmed: seeds 3-5 Recall@100 0.828→0.818, generic drugs in the top 10 13 %→0 %; areas 0.612→0.607, 32 %→17 %; external 0.352→0.338 | Validated (cost documented) |
| P14 | Abstention rule: of the first 25, only those with an explainable path | Benchmark: explained ones hit 18.9 %, without a path 3.1 % | Validated (the benchmark's definition of «path» and the explainer's roughly coincide) |
| P15 | Risk warning = contraindicated > indicated in neighbouring/similar diseases, + checkpoints | 97.4 % are real contraindications (baseline 77 %) | Validated |
| P16 | «No prior use» = no indication, Open Targets approval or off-label use; no diagnostic agents | Reviewer: 26 % of the list were approved drugs | **Error fixed** |
| P17 | Wrong Open Targets «approved» entries (e.g. lincomycin for MS) | Cause: Open Targets reads FDA drug labels and turns other substances mentioned there into «approved». 241 of 11,175 approvals (2.2 %) only have labels with more than one drug; 84 have no approval report at all | Partial: limitation quantified; **fix at the source (sources/opentargets.py) at the next rebuild** |
| P18 | Shared phenotypes do not count as evidence for showing a candidate | Benchmark: candidates with phenotypes only hit 2.1 % (without a path 3.0 %; the rest 19.5 %) | Validated |
| P19 | Diagnostic agents removed from the candidates; experimental and withdrawn drugs labelled (not excluded) | Excluding experimental drugs would lose 5 % of real positives | Validated (presentation) |
| P20 | CTD removed from production (model, explanations, benchmark) | Ablation: zero-shot 0.813→0.811, areas 0.631→0.634, external 0.352→0.359 (noise) | Validated |
| P21 | Evidence for showing a candidate: target, neighbouring disease, nearby (2 steps) or similar with ≥ 3 genes in common | Criterion fixed (≥ 10 % hits on seeds 1-2), confirmed on 3-5: target 24 %, neighbour 47 %, nearby 18 %, similar≥3 10 %; discarded: similar 1-2 genes 2 %, interaction 0 %, similar drug 0 %, phenotypes 0 % | Validated (2nd review) |
| P22 | Twin nodes (same drug/disease in two nodes) re-linked by name | 2nd review: 15 % of «no use» were approved; IPF did not show 130 trials | **Error fixed**; tests |
| P23 | «vet_approved» is not approved; «approved and withdrawn» = «withdrawn in some market»; veterinary-only removed | 2nd review (an insecticide for schizophrenia) | **Error fixed**; tests |
| P24 | Drugs contraindicated in the same disease removed from the candidates | 2nd review: 17 cases (succinylcholine in Duchenne…) | **Error fixed**; test |
| P25 | A single implementation of the disease page (`web.build.DiseasePages`) for the website and `explore.py` | 2nd review: mismatched criteria | Done |
| P26 | Method figures = OREVID-bench, with a note that specificity costs ~0.01 | 2nd review: they said «exactly the configuration» and they were not | **Error fixed** |
| P27 | Evidence by profile (not by hierarchy) and confirmation on the 9 areas (`evaluation/evidence.py`). Replaces P14/P21 | P21 counted «target» first and absorbed the combinations: **target alone hits 1.8-3.4 %** (like no path, 1.5-2.5 %). «2 steps in the ontology» alone: 10.6 % on the benchmark but 3.5 %/1.1 % on the areas. Levels A>B not confirmed (they invert on the hard areas). New rule: only «indicated in a neighbouring disease» (→ «direct extensions», 50-59 %; 13-25 % hard) or «indicated in a disease with ≥ 3 genes in common» (→ table). Phase 3 no longer counts as an indication. No level letters | Validated (benchmark + areas) |
| P28 | Shared genes: additionally, significant (hypergeometric p < 1e-6) | Expert review: «shares 37 genes with breast cancer» is chance. S candidates with p < 1e-6: 22.7 % (benchmark) / 25.2 % (areas); with p ≥ 1e-6: 0-2.5 %. Final table: 21-27 % hits; hidden 2-4 % | Validated (benchmark + areas) |
| P29 | Mechanistic coherence and drug class (expert review) | Criterion fixed beforehand: remove a subset if it hits < 5 % on both sets. Drug target on the shared genes: 23 % vs 22 % (no effect); anti-infective in a non-infectious disease: 33 % vs 0 % (n=45, unstable; acne, ulcer…). No filter adopted. ≥ 3 **specific** shared genes (< 100 diseases): 32-33 % vs 9-11 % → shown (strong/weak chip), does not filter | Validated as annotation |
| P30 | Expert review of the disease page: table format (Open Targets style), ATC class, targets with action, current indications, shared genes with p, intermediate disease at 2 steps, CSV | Errors found: Fc receptors (FCGR*) and C1q as «targets» of antibodies (sarilumab did not show IL6R) → removed if there are others; salts and esters (prednisolone proposed for IPF where prednisolone phosphate is approved; estradiol benzoate/dienanthate duplicated) → `active_moiety`; status «withdrawn in some market» untranslated. Adverse effects coinciding with phenotypes (SIDER): they do not separate contraindications from indications (71 % vs baseline 78 %) → information, not a warning | **Errors fixed**; tests |

| P31 | 3rd independent review (agent in the role of a clinical pharmacologist) on 84 disease pages | **Critical**: the production S line counted Open Targets approvals (the validated one, PrimeKG only) → now only `indication`. «Approved» showed false Open Targets approvals → `sources/ot_approvals.py` classifies them: confirmed only by a label or agency for a single drug with the disease name written out (5,390), the rest «unconfirmed» (4,825; includes lincomycin-MS from a 3-molecule label, atomoxetine-Alzheimer's from «ad» = ADHD at the PMDA, leflunomide-MS from TTD; agreement with DrugBank: label/agency 26 %, TTD only 14 %, ATC only 4 %). Direct extensions now within the first 25 (as in the evaluation). External exam in `external()` and `docs/EVIDENCE.md` (table 4.2 %). **Others**: prodrugs/enantiomers (`drug_family`: teriflunomide, dexibuprofen, fosphenytoin), skin-test agents removed, ecotoxicology removed from target risks, targets without direct evidence collapsed with the validation figure, MS4A2 as an Fc receptor, whole-word translation, MONDO terms grouped, searches without «(disease)», extensions grouped by class, total trials, gene sources, legal notice (not a medical device, known errors, version). Risk-warning text with the baseline (77 %) and what it does not detect. **Limitations stated, not filtered**: gene-disease errors at the source (endocarditis with SLE genes → antibiotics), comorbidities (RA-Graves'), opposite mechanism and class toxicity without a warning | **Errors fixed**; tests (10) |

## Errors detected in this session (and how)

1. Leak in targets (subtracting the disease's own drugs): the result improved when information was removed.
2. Checkpoint rule too broad: by validating it against all diseases.
3. Risk labels without a reason: by looking at the distribution of the labels.
4. Empty filter and line breaks inside the code: escaping errors (heredocs). Measure: edit tool and tests.
5. Generic candidates (norfloxacin in 42 % of diseases) and candidates without explanation: systematic panel.
6. Independent reviews (1st): direction model = copy of neighbours; «biology only» with drug
   information; approved drugs in the «no use» list; ontology metapath including the same disease.
7. Independent review (2nd): twins, weak evidence accepted, veterinary drugs, contraindicated drugs in the list,
   Method figures misattributed. Fixed (P21-P26).
