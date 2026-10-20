# Project status and how to continue (updated 20-10-2026)

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
| Preprint | Preprints.org, DOI 10.20944/preprints202610.1918.v1 (https://www.preprints.org/manuscript/202610.1918/v1). Source and exports in `docs/preprint/` (Word, Markdown, figures); the editable source is a Claude Docs document in the maintainer's account |
| Contact | contact@orevid.org (forwarded by the registrar) |
| Data (`data/`, ~5 GB) | Not in git. Rebuild with `build_all.sh` (see README); several hours |
| Benchmark prediction files | GitHub release assets (`orevid-bench-predictions-v1.1.zip`) |

On the home PC the working copy is `C:\dev\farmagraf` (local history) and `sync_public.sh` copies it
to `C:\dev\orevid` (public clone). **On any other computer, work directly in a clone of
orevid-org/orevid**; when back on the home PC, pull there and copy the changes back to
`C:\dev\farmagraf` (or switch to working only in the public clone).

TxGNN runs only on the home PC (WSL: `~/txgnn-env`, `~/TxGNN`, `~/txgnn_data`; RTX 3060, ~95 min
per split). Its predictions are in `data/eval/txgnn_predictions/` and in the benchmark release.

## Current state (18-10-2026)

- Graph: 157,822 nodes, 9,457,076 relations (PrimeKG + Open Targets 26.06 + STRING 12.5;
  PubTator3 only as displayed context).
- Production model: LightGBM LambdaRank on 22 metapaths + popularity + drug degrees, masked
  training (`predictor.py`; evaluation in `evaluation/txgnn.py`, feature set `extra_no_ctd`, 5 rounds). Unchanged
  since 12-10: every alternative tested since (P33, P43, P45–P47) failed its pre-registered rule.
- Website shows only candidates with an evidence line (direct extension; ≥ 3 shared genes with
  hypergeometric p < 1e-5 over the genes with a disease, P40) among the top 25; hit rates in `docs/EVIDENCE.md`.
  Each candidate also shows generic / biosimilar status and how much it was studied (P44). Logo, favicon and a
  GoatCounter visit counter (no cookies; dashboard https://orevid.goatcounter.com, account contact@orevid.org).
- OREVID-bench: `benchmark/` (leaderboard with popularity, TxGNN, RotatE and OREVID; fair comparison with TxGNN in
  `benchmark/FAIR_COMPARISON.md`; disease-level intervals in `benchmark/CLUSTER_CI.md`).
- Key numbers (Recall@100): zero-shot 0.811 vs TxGNN 0.412 (RotatE 0.289); areas 0.617 vs 0.237; external exam
  0.359 vs 0.121 (RotatE 0.156); temporal test at 2015 0.322 vs popularity 0.174 (`docs/TEMPORAL.md`). Hit rate on
  unseen approvals: direct extensions 12.6 %, shared genes 4.3 %.
- Preprint: v1 online on Preprints.org since 20-10 (preprints-238326, DOI 10.20944/preprints202610.1918.v1,
  registered with Crossref). The Claude Doc and `docs/preprint/` already hold the v2 text (see «Done on 16–18 October»).

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
  The website change (legal notice, About page) was published by committing only those files on top of `gh-pages`
  (`main`'s `web/public` matches the deployed site again; the next `web/deploy.sh` replaces that commit as usual).
- Account tasks done: ORCID created; Gmail «send as» for contact@orevid.org and the GitHub Pages TXT record verified.
- **GitHub release v1.1.0** (https://github.com/orevid-org/orevid/releases/tag/v1.1.0): the tag points to `7965ccf`,
  the commit that produced the prospective snapshot (as `docs/PROSPECTIVE_PLAN.md` requires), not to the current code.
  Zenodo version DOI **10.5281/zenodo.23362377** (concept 10.5281/zenodo.23323226). Asset
  `orevid-bench-predictions-v1.1.zip` = the v1.0 predictions of the four leaderboard models (the ranking model did
  not change) + the `results.json` of all five; the prediction files of the train-only variant are only on the home
  PC (add them in the next release). Zenodo took the author name from `CITATION.cff` at that commit; the creators
  were then corrected by hand to «Quadrada Jané, Joaquim» with ORCID in the three records (23323227 v1.0.0,
  23362377 v1.1.0, 23326593 snapshot; checked through the Zenodo API). From the work PC, zenodo.org fails in Chrome
  because the ESET SSL filter issues a `*.zenodo.org` certificate that does not cover the bare domain: do not bypass
  the warning (use another network or device).
- Biomarker table extended (8 rows, `resources/biomarker_subtypes.csv`); the graph has no molecular subtypes of lung,
  melanoma or colorectal cancer (EGFR, ALK, BRAF, KRAS…). New exclusion: acute myeloid leukaemia with mutated NPM1
  (99930) has no BCR::ABL1 (separate entity in WHO 2022 / ICC 2022) — its page proposed imatinib, ponatinib,
  nilotinib, dasatinib and radotinib through genes shared with CML. New annotations (marker present): Philadelphia-
  positive myelogenous leukaemia, PDGFRA- and JAK2-rearranged myeloid/lymphoid neoplasms, primary cutaneous CD30+
  T-cell lymphoproliferative disease, RELA fusion-positive ependymoma, NPM1. Not yet applied to the website.
- **bioRxiv rejected the preprint on 14-10-2026** (BIORXIV/2026/759044): «bioRxiv requires authors to have an
  organizational affiliation» that provides oversight of research; «Independent researcher» does not qualify (medRxiv
  has the same owner and rules). The submitted file is kept as `docs/preprint/OREVID_preprint_biorxiv_v1.pdf`.
  bioRxiv account notes, in case of a future submission with an affiliation: primary e-mail contact@orevid.org,
  ORCID linked; its fields reject accents, en dashes and ±; if «Approve Manuscript» stays inactive after changing the
  account e-mail, delete the author and add it again with «I AM THIS AUTHOR».
- **Going to Preprints.org (MDPI) instead**: no affiliation required (screening checks «genuine scholars», English,
  ethics, conflicts, AI disclosure; ~24 h; always CC BY 4.0; Word or LaTeX, not PDF). The Claude Doc was adapted:
  keywords on the first page, AI use also in Methods (end of «Statistical analysis and software»), «Conflicts of
  Interest» as its own paragraph just before the references. Word file and form values prepared (not in git).
  Preprints.org advises against posting the same paper on several servers: do not post it on bioRxiv later.
- **Submitted to Preprints.org on 14-10-2026: manuscript ID preprints-238326** (MDPI account; subject Medicine and
  Pharmacology › Pharmacology and Toxicology; no MDPI Topic and no «friendly journal» chosen; author marked as
  independent researcher, Spain, ORCID; Word + PDF + graphical abstract (Figure 1); related websites: GitHub, both
  Zenodo concept DOIs, orevid.org). The files Preprints.org generated were checked: identical to the ones uploaded.

## Done on 14-10-2026, evening (home PC)

- **Biomarker table applied (P41)**: `evaluation.biomarker_check` passes (0 approved conflicts, 14 candidates
  removed; the five BCR-ABL1 inhibitors leave the AML-with-mutated-NPM1 page). Gram-stain rule tried with a new
  `drugs:` rule type and **rejected** (criterion b: 14 approved conflicts through the ontology); antifungal analysis
  done (amphotericin B comes from a PrimeKG error; a class rule would remove TMP-SMX for pneumocystosis). Both in the
  hand-written section of `docs/BIOMARKERS.md`, which the script now keeps.
- **I1 temporal validation (P42)**: `evaluation/temporal.py` → `docs/TEMPORAL.md`. The pre-registered rule is met:
  at 2015 OREVID 0.322 vs popularity 0.174 Recall@100, paired +0.148 (0.123–0.173); non-obvious +0.043. Method page
  and preprint now claim it, with the undated-edge leakage (80 % of PrimeKG indications) as the limit.
- **I2 calibration (P43)**: `evaluation/calibration.py` → `docs/CALIBRATION.md`. Not adopted (Brier gain 2.5 % < 5 %).
- **I3 practical information (P44)**: `sources/practical.py` (Orange Book, Purple Book, EMA; `download`, `build`) →
  `data/ext/practical.parquet`; columns «Generic» and «Studied», filters and CSV fields on the website;
  `docs/PRACTICAL.md` (match rates, manual check of 50 drugs, two matching errors fixed before release).
  Website rebuilt in full and deployed (15-10, early morning).
- **I4 expression-signature reversal (P45)**: `evaluation/lincs_reversal.py` → `docs/LINCS.md` (Enrichr libraries
  `LINCS_L1000_Chem_Pert_Consensus_Sigs`, `Disease_Perturbations_from_GEO_up/down` in `data/raw/lincs/`). 2,410 drugs,
  239 diseases; at chance for drugs with no prior use (temporal 2015: 0.055 vs random 0.055; OREVID 0.210).
  Negative result; the LambdaRank-feature test (b) was not run, by the rule written before running.
- **I6 shared-biology line (P46)**: `evaluation/biology_lines.py` → `docs/BIOLOGY_LINES.md`. Pathway overlap chosen
  on the benchmark (29.8 % vs 23.0 %, half the candidates) but 5.0 % on the external exam, not above the null's
  99th percentile (5.32 %): not adopted. Pathways are the best lead for a future, separately pre-registered test.
- **I5 embedding baselines (P47)**: `evaluation/kge_bench.py` (+ `.sh`, WSL env `~/pykeen-env`) → `docs/KGE.md`,
  `benchmark/predictions/kge_rotate/`. Grid on zero_shot_1, then RotatE only on every task (16–17 October; TransE and
  ComplEx dropped as dominated, documented). RotatE is below popularity on zero-shot (0.289 vs 0.358), level with
  TxGNN on the areas and above TxGNN on the external exam (0.156 vs 0.121; non-obvious 0.127 vs 0.085). The
  OREVID + RotatE ensemble fails the P33 rule (`evaluation/kge_ensemble.py`). Leaderboard and `CLUSTER_CI.md`
  include RotatE (intervals only; the Holm family is unchanged). The prediction files go in the next release.

## Done on 16–18 October 2026 (home PC)

- **I5 finished (P47)**: RotatE run on seeds 2–5 and the nine areas (two WSL processes; RotatE at dimension 256
  saturates the GPU, so parallel runs gain nothing: ~45 min per run). Leaderboard, `CLUSTER_CI.md`, README, Method
  page (three languages) and preprint updated; the OREVID + RotatE ensemble was tested and not adopted.
- **Website**: OREVID logo in the header (letters as web text and the E as three bars the height of the capitals,
  so it aligns with any font and works in dark mode; the original SVG had the bars too high, a corrected copy is in
  `web/public/img/orevid_logo.svg`), «Open REpurposing EVIDence» next to it (hidden on narrow screens), favicon,
  and the GoatCounter counter: every `#/` route counts as a page, sent once the page has its title; privacy section
  of the legal notice updated in three languages.
- **Preprint v2 text ready, not yet posted** (changes against the submitted v1, commit `6d36551`): temporal
  validation (Methods paragraph, Results subsection, one sentence in the abstract, Limitations and Discussion);
  RotatE row in Table 1, a Results sentence and the Discussion paragraph that said no embedding baseline was run;
  the LINCS and shared-biology negative results in Limitations; practical information in «Website». No earlier
  number or figure changes.
- **Open Targets #4573**: David Ochoa answered (18-10) that they will keep multiple sclerosis as a disease and the
  severity scores as measurements, fix due in 26.12; see «Data pitfalls» for what to check before upgrading.

## Done on 20-10-2026 (work PC)

- **Preprint v1 online** on Preprints.org: DOI **10.20944/preprints202610.1918.v1**, registered with Crossref.
  The DOI is now in the README (badge and «Citing»), `CITATION.cff` (`preferred-citation`, so GitHub's «Cite this
  repository» gives the preprint) and the website: About page (how to cite) and a line at the top of the Method page,
  in three languages. The website change was published by committing only the three `i18n/*.json` on top of
  `gh-pages` (as on 14-10; `main`'s `web/public` matches the deployed site).
- At the next GitHub release, check that Zenodo still takes the software's title and type from the top-level fields of
  `CITATION.cff` and not from `preferred-citation`; if it does not, fix the record by hand.

## Next steps

**On a computer without `data/`** (the work PC): posting the preprint (step 1), the GitHub release (5), the biomarker table (3: editing the CSV; checking it needs `data/`) and
anything that only reads the benchmark prediction files (`benchmark/`, release zip). **Needs the home PC**
(`data/`, ~5 GB, and WSL for TxGNN): new evaluations, `evaluation/evidence_ci.py`, `web.build` (also `--patch`)
and the website deploy. The Claude Doc of the preprint is the source of truth for its text; re-export to
`docs/preprint/` after editing it (`docs/REVIEW_RESPONSE.md`, «After the pending tasks»).


1. **Preprints.org** (preprints-238326, online 20-10, DOI 10.20944/preprints202610.1918.v1): DOI added to README,
   website and `CITATION.cff`, and to ORCID (20-10; the Zenodo records were already there through DataCite). New versions go through the
   User Dashboard («My Preprints»; keep the same server; each version gets its own DOI, `.v2`).
   **Post v2 now that v1 is online**: export Word from the Claude Doc (`docs/preprint/OREVID_preprint_draft.docx` is
   the current export) with the changes listed in «Done on 16–18 October».
   Then choose a journal (independent authors are accepted by journals).
   **Journal plan (14-10):** first PCI Mathematical and Computational Biology (free, open reviews; then Peer
   Community Journal, free), Briefings in Bioinformatics as the fallback; JOSS for the software after April 2027
   (repository public for 6 months). In parallel, write to the TxGNN authors (Zitnik lab), Every Cure and David Ochoa.
   **Before PCI:** PCI does not accept AI-generated illustrations, and figures 1–4 were drawn by the AI assistant in
   the Claude Doc. `evaluation/figures.py` (written with AI assistance, to be read and run by the author) draws the
   four figures from the result files (`benchmark/predictions/*/results.json`, `docs/EVIDENCE_CI.md`): `py -3 -m pip
   install matplotlib`, `py -3 -m evaluation.figures` → `docs/preprint/figures/figure{1..4}.png|.pdf`. Then replace the
   four embedded charts of the Claude Doc by these images, add to Methods that the figures were produced by that
   script, number the lines (Word: Layout → Line Numbers → Continuous) and use PCI's conflict sentence («The authors
   declare that they comply with the PCI rule of having no financial conflicts of interest…»). Post the new version
   on Preprints.org (same server) before submitting to PCI.
1b. **Improvement plan**: done (I1–I6, P42–P47). Left over: the next OREVID-bench release must include
   `benchmark/predictions/kge_rotate/` (and the train-only files, item 5); optional TxGNN on the temporal cutoffs;
   pathway overlap with its own null (P46) as a new pre-registered test.
2. Shared genes (P38): the line works as a selection of drugs used in comparable diseases, not as biology. The
   three pre-registered replacements failed (P46); pathway overlap came closest and would need its own null
   (permuting pathway sets) in a new pre-registered test.
3. Biomarker table: applied on 14-10 (P41); Gram-stain and class rules rejected (`docs/BIOMARKERS.md`). New rows
   go through `py -3 -m evaluation.biomarker_check` and `py -3 -m web.build --patch` as before.
4. Model: production stays (P33). Optional: a pre-specified test of the LightGBM + logistic ensemble on new uses.
5. Next release (v1.2.0, current code): include the train-only and RotatE (`benchmark/predictions/kge_rotate/`) prediction files and update `CITATION.cff`
   `version`/`date-released` first (the name and ORCID are already right there).
6. ORCID record: done (preprint added by hand; both Zenodo concept DOIs added automatically by DataCite, with their versions; website links and keywords filled in).
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
  (EFO_0803536): use 26.06. Always check `obsoleteTerms` before saying a term is missing. On 18-10-2026 David
  Ochoa answered #4573 that they will keep the disease and the measurement apart, prototype the fix and roll it
  out in **26.12**. Another user (madaraviv) listed in that issue 92 MONDO ids merged into an EFO id in 26.09
  (e.g. MASLD into «non-alcoholic fatty liver», aortic valve stenosis, influenza, Chagas disease, peripheral
  neuropathy, Castleman disease). Before moving to 26.12: check that MONDO_0005301 is back as a disease, and that
  those 92 ids map to the same graph diseases as in 26.06 (`sources.ids.diseases()` resolves merges through
  `obsoleteTerms`; compare the mapping and the number of positives of the external exam).
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
