# Project status and how to continue (updated 14-10-2026, afternoon)

Hand-over document between sessions and computers. Read this first, then `README.md` and
`docs/DECISIONS.md` (decision log with the validation status of every rule).

## Goal

Build the best open predictor of which existing drugs may work for each disease, with an honest
evaluation against the alternatives (TxGNN; Every Cure's MATRIX). Hypothesis generator, not a
clinical tool. Non-commercial use of some sources (data under CC BY-NC-SA 4.0; code MIT).

## Where everything is

| What | Where |
|---|---|
| Public code and docs | https://github.com/orevid-org/orevid (this repository) |
| Website | https://orevid.org (GitHub Pages, branch `gh-pages`; deploy: `DOMAIN=orevid.org REMOTE=https://github.com/orevid-org/orevid.git bash web/deploy.sh`) |
| Code archive | Zenodo concept DOI 10.5281/zenodo.23323226 (GitHub release integration) |
| Prospective snapshot (12-10-2026) | Zenodo DOI 10.5281/zenodo.23326593; plan in `docs/PROSPECTIVE_PLAN.md` |
| Preprint draft | `docs/preprint/` (Word, Markdown, figures). The editable source is a Claude Docs document in the maintainer's account |
| Contact | contact@orevid.org (forwarded by the registrar) |
| Data (`data/`, ~5 GB) | Not in git. Rebuild with `build_all.sh` (see README); several hours |
| Benchmark prediction files | GitHub release assets (`orevid-bench-predictions-v1.1.zip`) |

On the home PC the working copy is `C:\dev\farmagraf` (local history) and `sync_public.sh` copies it
to `C:\dev\orevid` (public clone). **On any other computer, work directly in a clone of
orevid-org/orevid**; when back on the home PC, pull there and copy the changes back to
`C:\dev\farmagraf` (or switch to working only in the public clone).

TxGNN runs only on the home PC (WSL: `~/txgnn-env`, `~/TxGNN`, `~/txgnn_data`; RTX 3060, ~95 min
per split). Its predictions are in `data/eval/txgnn_predictions/` and in the benchmark release.

## Current state (12-10-2026)

- Graph: 157,822 nodes, 9,457,076 relations (PrimeKG + Open Targets 26.06 + STRING 12.5;
  PubTator3 only as displayed context).
- Production model: LightGBM LambdaRank on 22 metapaths + popularity + drug degrees, masked
  training (`predictor.py`; evaluation in `evaluation/txgnn.py`, feature set `extra_no_ctd`, 5 rounds).
- Website shows only candidates with an evidence line (direct extension; ≥ 3 shared genes with
  hypergeometric p < 1e-5 over the genes with a disease, P40) among the top 25; hit rates in `docs/EVIDENCE.md`.
- OREVID-bench: `benchmark/` (leaderboard, fair comparison with TxGNN in `benchmark/FAIR_COMPARISON.md`).
- Key numbers (Recall@100): zero-shot 0.811 vs TxGNN 0.412; areas 0.617 vs 0.237; external exam
  0.359 vs 0.121. Hit rate on unseen approvals: direct extensions 12.6 %, shared genes 4.3 %.

## Done on 12-10-2026

1. **Fairer comparison with TxGNN** (P32): training-labels-only OREVID and exclusion of diseases
   where TxGNN outputs constant scores; conclusion unchanged.
2. **About page**, feedback links (GitHub issue templates + e-mail), Spanish added, English default.
3. **Prospective snapshot**: `evaluation/snapshot.py` froze the top 100 per disease with its
   website section; `docs/PROSPECTIVE_PLAN.md` (pre-registered analysis, SHA-256 of the files);
   archived on Zenodo (10.5281/zenodo.23326593). First outcome check: October 2027 (12 months).
4. **Preprint** (bioRxiv, English): written, two rounds of review by an AI agent (Fable) with all
   points applied; statistics (SD, paired differences, Wilson CIs), per-area table, ablations,
   related work, examples, AI-use statement. Draft exported to `docs/preprint/`.
5. **XGBoost, CatBoost and ensembles vs LightGBM** (`evaluation/xgboost_compare.py`,
   `docs/XGBOOST.md`, P33): protocol fixed before running, production reproduced exactly. None
   beats the production LightGBM LambdaRank (Recall@100 zero-shot 3–5 / areas / external: production
   0.825 / 0.617 / 0.359; XGBoost 0.810 / 0.612 / 0.348; CatBoost 0.812 / 0.560 / 0.336). The
   LightGBM + logistic-regression ensemble is better on the external exam (0.374; non-obvious 0.196
   vs 0.178) but worse elsewhere: production kept. Training rows are cached in `data/eval/xgb/`
   (~8 GB, delete when no longer needed).
6. Bug fixed: a train-only run of `evaluation/txgnn.py` overwrote the standard results CSV.
7. **Open Targets replied** to #4579 (David Ochoa): lincomycin–MS comes from ChEMBL's manual
   curation; the acronyms are the OnToma bug of #4573. He asked for a script or dataset to find more
   cases: `sources/ot_suspect_approvals.py` and `resources/opentargets/ot_suspect_approvals_26.06.csv`
   (14 acronym pairs, 13 wrong; 729 multi-molecule pairs, 355 high priority). Analysis, proposal and a
   reply posted by the maintainer on 13-10 (the analysis note is kept locally, not in the public repository).

## Done on 13-10-2026 (from the work PC, clone at `C:\dev-quimqu\orevid`)

1. **Reply to an external (AI) review of the preprint** (P34): `benchmark/cluster_ci.py` → `benchmark/CLUSTER_CI.md`,
   disease-level bootstrap intervals, paired differences, Wilcoxon + Holm, tie-breaking sensitivity. Needs only the
   release prediction files. The Claude Docs preprint was updated (Statistical analysis, Results, Discussion,
   Limitations, abstract, new Supplementary Table S4) and re-exported to `docs/preprint/`. Point-by-point
   answer and what is left: `docs/REVIEW_RESPONSE.md`.

## Done on 13–14 October 2026 (home PC): the eight review tasks T1–T8

Results and decisions in `docs/REVIEW_RESPONSE.md` (section «Results of the pending tasks») and DECISIONS P36–P39.

1. **T1** — no material dependence on Open Targets: without any Open Targets feature, external Recall@100 0.348
   (production 0.359; PrimeKG only 0.351), zero-shot 0.802, areas 0.628.
2. **T2** — external exam over five training seeds, OREVID and TxGNN: `docs/EXTERNAL_SEEDS.md`.
3. **T3, T5** — hit rates with disease-level intervals; the website's exact rules give the same figures
   (`evaluation/evidence_ci.py` → `docs/EVIDENCE_CI.md`). The model's ranking does most of the work (top 25 without
   rules 4.2 %; rules without the model 1.0 % after P40).
4. **T4, T8** — threshold grid: no sharp optimum; **permutation null failed**: random gene sets of the same size
   pick candidates that are right as often, also within gene-set-size bins and with specific genes (exploratory:
   `py -3 -m evaluation.evidence_ci null`). The shared genes are no longer presented as evidence of a shared
   mechanism (website wording and preprint changed; candidates unchanged, P38). Broad diseases keep the section.
5. **T6** — ontology navigation (`sources/hierarchy.py`), broad-term notice, near-synonym marker.
6. **T7** — biomarker rule (`resources/biomarker_subtypes.csv`), validated by hand (`docs/BIOMARKERS.md`):
   tamoxifen and other endocrine drugs no longer appear for triple-negative or ER-negative breast cancer.
7. Website pages updated with `py -3 -m web.build --patch` (navigation fields added, 10 biomarker pages rebuilt;
   checked to equal a full rebuild on 40 random pages) and deployed. Preprint updated in the Claude Doc and
   re-exported to `docs/preprint/`. Open Targets #4579: reply posted by the maintainer on 13-10.

## Done on 14-10-2026, night (home PC)

- **P40**: the shared-genes hypergeometric test now uses the genes with a disease as population (19,735 in the evaluation, 19,766 on the website; it counted diseases too)
  with the threshold re-matched to p < 1e-5 (protocol and results: `docs/HYPERGEOMETRIC_N.md`). All evidence
  figures recomputed (`docs/EVIDENCE.md`, `docs/EVIDENCE_CI.md`); the website rebuilt in full and deployed (Method
  page and help texts in three languages); preprint, README and decision log updated. P38 still holds.

## Done on 14-10-2026 (work PC)

- Preprint read end to end by the author: correct. Author name is **Joaquim Quadrada Jané** everywhere
  (preprint, `CITATION.cff`, website legal notice and About page in three languages), with ORCID
  [0009-0000-9624-332X](https://orcid.org/0009-0000-9624-332X) on the preprint and in `CITATION.cff`.
  PDF for bioRxiv exported from the Claude Doc (22 pages, A4).
- Account tasks done: ORCID created; Gmail «send as» for contact@orevid.org and the GitHub Pages TXT record verified.

## Next steps

**On a computer without `data/`** (the work PC): posting the preprint (step 1), the GitHub release (5), the biomarker table (3: editing the CSV; checking it needs `data/`) and
anything that only reads the benchmark prediction files (`benchmark/`, release zip). **Needs the home PC**
(`data/`, ~5 GB, and WSL for TxGNN): new evaluations, `evaluation/evidence_ci.py`, `web.build` (also `--patch`)
and the website deploy. The Claude Doc of the preprint is the source of truth for its text; re-export to
`docs/preprint/` after editing it (`docs/REVIEW_RESPONSE.md`, «After the pending tasks»).


1. Post the preprint to bioRxiv (submit.biorxiv.org, log in with ORCID; New Results, Bioinformatics, CC BY;
   single PDF exported from the Claude Doc). Then add the bioRxiv DOI to the README, the website, `CITATION.cff`
   (preferred-citation) and the ORCID record (Works).
2. Shared genes (P38): the line works as a selection of drugs used in comparable diseases, not as biology. Worth
   exploring a replacement that does beat the permutation null (e.g. shared genes weighted by specificity, or
   pathway-level overlap), pre-registered as before. Small diseases (< 100 genes) were the only ones above the
   null on the benchmark.
3. Biomarker table: extend `resources/biomarker_subtypes.csv` beyond breast cancer, CML and ALCL (e.g. EGFR- or
   ALK-mutant lung cancer if the graph has those nodes), with the same manual check (`evaluation/biomarker_check.py`).
4. Model: production stays (P33). Optional: a pre-specified test of the LightGBM + logistic ensemble on new uses.
5. GitHub release v1.1.0 with `orevid-bench-predictions-v1.1.zip` (the prospective plan cites v1.1.0).
6. ORCID record (optional): add the two Zenodo DOIs as Works, website and keywords.
7. Collect expert feedback from the website (GitHub issues) as new validation labels.
8. Open problems: drugs with no prior indication (Precision@K 0–2 %); wrong gene–disease
   associations and comorbidities that read as shared biology (lupus–endocarditis, RA–Graves).

## Rebuild order (from scratch)

```
bash build_all.sh          # downloads, sources, prepare.py, approvals, targets, predictor
py -3 -m pytest tests      # invariants
```

Details and individual steps: README. Any change in `sources/*` requires `prepare.py` again.

## Results to remember

- **Popularity baseline** (rank drugs by number of indications, ignoring the disease) is strong
  (zero-shot Recall@100 0.358, close to TxGNN's 0.412): every claim must be compared with it.
- Mechanistic evidence alone (drug acts on a disease gene) hits 1.8–3.4 %, like no evidence.
- Direct extensions (indication in an ontology neighbour) are largely rediscovery: 50–59 % when
  hiding treatments, 13–25 % on hard diseases, 12.6 % on unseen approvals.
- Validation labels change OREVID by −0.005 to +0.027; degenerate TxGNN diseases are few (8/335, 3/345).

### Things tried that do NOT work (do not repeat)

- CTD direction comes out inverted; PubTator and expression direction are null; IMPC always says
  «activate». CTD as features adds nothing (removed, P20).
- Direction as a metapath for finding candidates does not discriminate.
- Open Targets «affected_pathway»: identical 704-gene blocks for different diseases. Out.
- Open Targets «known_drug» and literature «treats» are the answer: never as features.
- GWAS+QTL direction: suggestive but not conclusive; in the graph as `dir_gwas`, not scored.
- Chemical structure features and drug combinations do not help.
- Filters suggested by reviewers (target must touch the shared genes; no anti-infectives in
  non-infectious diseases) failed validation; a three-level grading did not replicate.

### Data pitfalls

- Open Targets 26.09 merges multiple sclerosis (MONDO_0005301) into a measurement term
  (EFO_0803536): use 26.06. Always check `obsoleteTerms` before saying a term is missing.
- Open Targets approvals: some come from multi-drug labels (lincomycin–MS) or abbreviations
  («ad» = ADHD read as Alzheimer's); `sources/ot_approvals.py` classifies them.
- In 26.06, `clinical_report` has no `origin` (the `type` field gives the source).
- CTD's InChIKey column is empty; chemicals are linked by name.
- Dataverse, CTD and NCBI return 403 without a User-Agent.
- Use `py -3` (not `python`) on Windows. Browsers cache `app.js`: bump `?v=N` in `index.html`.

## External context

- Every Cure / MATRIX (github.com/everycure-org/matrix): temporal validation, off-label, failed
  trials as negatives. No popularity baseline found there.
- TxGNN (Nature Medicine 2024): GNN over PrimeKG; the direct comparison (same graph).
- Issues reported to Open Targets: #4573 (multiple sclerosis merged in 26.09) and approvals
  derived from multi-drug labels and abbreviations.
