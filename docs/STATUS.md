# Project status and how to continue (updated 13-10-2026)

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
  p < 1e-6) among the top 25; hit rates in `docs/EVIDENCE.md`.
- OREVID-bench: `benchmark/` (leaderboard, fair comparison with TxGNN in `benchmark/FAIR_COMPARISON.md`).
- Key numbers (Recall@100): zero-shot 0.811 vs TxGNN 0.412; areas 0.617 vs 0.237; external exam
  0.359 vs 0.121. Hit rate on unseen approvals: direct extensions 12.6 %, shared genes 4.2 %.

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
   draft reply (not posted): `docs/OPEN_TARGETS_FEEDBACK.md`.

## Done on 13-10-2026 (from the work PC, clone at `C:\dev-quimqu\orevid`)

1. **Reply to an external (AI) review of the preprint** (P34): `benchmark/cluster_ci.py` → `benchmark/CLUSTER_CI.md`,
   disease-level bootstrap intervals, paired differences, Wilcoxon + Holm, tie-breaking sensitivity. Needs only the
   release prediction files. The Claude Docs preprint was updated (Statistical analysis, Results, Discussion,
   Limitations, abstract, new Supplementary Table S4) and re-exported to `docs/preprint/`. Point-by-point
   answer and what is left: `docs/REVIEW_RESPONSE.md`.

## Next steps

00. **Home PC, first: the eight pending tasks (T1–T8) in `docs/REVIEW_RESPONSE.md`**, each with the files
    to change, the commands, a decision rule fixed in advance and the preprint text to update. T1–T5 answer
    the external review; T6–T8 fix website problems found on 13-10 (breast cancer page: no navigation to
    subtypes, tamoxifen #1 for triple-negative, broad terms with pan-cancer candidates). Order: T6 (ontology
    navigation, presentation only) → T1 (external exam without Open Targets) → T7 (biomarker rule for
    subtypes) → T3 (hit-rate CIs) → T5 (exact website rules) → T4 + T8 (thresholds, broad diseases) → T2
    (5 seeds). All need `data/`; T6–T8 need a website rebuild and deploy.
0. Post the reply to opentargets/issues#4579 (draft in `docs/OPEN_TARGETS_FEEDBACK.md`), after a look at the CSV.
1. Read the preprint draft end to end; add ORCID; post to bioRxiv (category Bioinformatics or
   Pharmacology and Toxicology). Then add the bioRxiv DOI to the README and the website.
2. Model: production stays (P33). Optional: a pre-specified test of the LightGBM + logistic
   ensemble on new uses (it gained on non-obvious external approvals); one sentence on P33 could be
   added to the preprint's Methods as a robustness check.
3. GitHub release v1.1.0 with `orevid-bench-predictions-v1.1.zip` (the prospective plan cites v1.1.0).
4. Pending account tasks: ORCID; Resend domain verification and Gmail «send as» for contact@orevid.org;
   verify the domain in the GitHub organisation's Pages settings (TXT record).
5. Collect expert feedback from the website (GitHub issues) as new validation labels.
6. Open problems: drugs with no prior indication (Precision@K 0–2 %); wrong gene–disease
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
