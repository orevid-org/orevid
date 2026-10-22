# Project status and how to continue (updated 22-10-2026)

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
| Code archive | Zenodo concept DOI 10.5281/zenodo.23323226 (GitHub release integration); v1.2.0 = 10.5281/zenodo.23545853 |
| Prospective snapshot (12-10-2026) | Zenodo DOI 10.5281/zenodo.23326593; plan in `docs/PROSPECTIVE_PLAN.md` |
| Preprint | Preprints.org, DOI 10.20944/preprints202610.1918.v2 (https://www.preprints.org/manuscript/202610.1918/v2; v1: `….v1`). Source and exports in `docs/preprint/` (Word, Markdown, figures); the editable source is a Claude Docs document in the maintainer's account |
| Contact | contact@orevid.org (forwarded by the registrar) |
| Data (`data/`, ~5 GB) | Not in git. Rebuild with `build_all.sh` (see README); several hours |
| Benchmark prediction files | GitHub release assets (`orevid-bench-predictions-v1.2.zip`, release v1.2.0) |
| Benchmark on Hugging Face | https://huggingface.co/datasets/orevid/orevid-bench (organization `orevid`, admin account `quimqu`; tasks, predictions of the six models, `bench.py`, dataset card; viewer: one config per model, one split per task) |
| Tool registry | bio.tools entry https://bio.tools/orevid (account `joaquim0`, e-mail contact@orevid.org; editable only by that account) |
| Bluesky | @joaquim.orevid.org (personal account, e-mail contact@orevid.org; handle = DNS TXT `_atproto.joaquim` at Namecheap, never delete) |

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
- Preprint: v1 online on Preprints.org since 20-10 (preprints-238326, DOI 10.20944/preprints202610.1918.v1), v2
  online since 22-10 (preprints-239415, DOI 10.20944/preprints202610.1918.v2); README, `CITATION.cff` and website
  point to v2. v2 has the old figure 3; the Claude Doc has the fixed one (see 21-10).

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
  `CITATION.cff` and not from `preferred-citation`; if it does not, fix the record by hand. (Checked for v1.2.0 on
  21-10: correct.)

## Done on 20-10-2026, afternoon (work PC)

- **Google Search Console**: domain property orevid.org verified with a TXT record at Namecheap (Advanced DNS, host
  `@`, next to the SPF record of the e-mail forwarding: never delete either). Indexing of the home page requested.
  Decision: be findable for researchers (home, About, Method), but **do not make the disease pages indexable**: today
  they live under `#/` routes, which search engines treat as one page; if static per-disease pages are ever built,
  mark them `noindex`. The risk is a patient landing on a page and self-medicating, not only a legal one.
- **First-visit notice** (`index.html` `#notice`, `app.js` `showNotice`, strings `notice_*` in the three locales): a
  modal on the first visit, whatever page it starts on, saying the candidates are hypotheses for researchers and must
  not be used to start, stop or change a treatment. Only «I understand» dismisses it (remembered in localStorage,
  key `notice_ok_v1`; bump the key to show it again after a wording change); Escape does not; the legal-notice link
  closes it for that visit. Deployed by committing the changed files on top of `gh-pages`.

## Done on 21-10-2026 (work PC)

- **Figure 3 fixed**: the x axis starts at -0.02 (`evaluation/figures.py`), so the popularity point of «Diabetes»
  (0.0) is no longer cut in half; `docs/preprint/figures/figure3.png|.pdf` regenerated with the Spyder venv's Python
  (the other three figures are unchanged). The v2 submitted to Preprints.org still has the old figure 3: replace the
  image in the Claude Doc and use the new one in the PCI submission (and in any v3).
- Preprints.org shows v2 as checked, but its DOI (`…1918.v2`) was not registered yet at doi.org / Crossref.
- **Release v1.2.0 postponed to the home PC**: the prediction files of RotatE and of the train-only variant
  (`benchmark/predictions/kge_rotate/*.parquet`, `orevid_lambdarank_train_only/*.parquet`) exist only there (this
  clone has their `results.json` only). Steps in «Next steps», item 5.

## Done on 21-10-2026, evening (home PC)

- Work-PC changes of 20–21 October brought into `C:\dev\farmagraf` (DOI, v2 files, figures, first-visit notice).
- **Figure 3 replaced in the Claude Doc** (rev 122): the doc still had the old image (identical to the one in v2);
  it now has the fixed `docs/preprint/figures/figure3.png`. Use the Claude Doc / that image for PCI.
- **Release v1.2.0 published** by the author (21-10, 22:39), tag `v1.2.0` → `6c9b734` (the commit with
  `CITATION.cff` 1.2.0). Before it: the 15 + 14 prediction files checked, RotatE re-evaluated (identical
  `results.json`), `CITATION.cff` and `benchmark/README.md` updated; zip and notes in `data/release/`
  (`orevid-bench-predictions-v1.2.zip`, 247 MB, 95 files; `notes_v1.2.0.md`).
- **Zenodo checked**: version record 10.5281/zenodo.23545853 has title «OREVID: Open Repurposing Evidence», type
  software, version v1.2.0, creator «Quadrada Jané, Joaquim» with ORCID: `preferred-citation` did not leak into it.
- **Zip attached** (21-10, with the GitHub CLI, now installed on this PC and logged in as quimqu):
  `orevid-bench-predictions-v1.2.zip`, 247,270,717 bytes, same size as the local file; release notes checked against
  `data/release/notes_v1.2.0.md`. Release v1.2.0 complete. `gh` lets future releases and uploads be done from here
  (`gh release create|upload … -R orevid-org/orevid`).

## Done on 22-10-2026 (work PC)

- **Preprint v2 online** on Preprints.org (dashboard «Online», 22-10): DOI **10.20944/preprints202610.1918.v2**.
  README (badge, documentation line, «Citing»), `CITATION.cff` (`preferred-citation` doi and url) and the website
  (Method and About pages, three languages) now point to v2; the `.v1` DOI keeps working. Pushed only once the v2
  DOI resolved at doi.org; the website change was published by committing only the three `i18n/*.json` on top of
  `gh-pages` (as on 14-10 and 20-10; `main`'s `web/public` matches the deployed site).
- ORCID: Crossref asked (e-mail of 21-10) for permission to add works automatically; once granted, v2 should appear
  by itself (the v1 entry was added by hand; ORCID groups entries that share a DOI).
- **Outreach channels that do not depend on followers** (the website is meant for researchers; no general press,
  disease pages stay non-indexable):
  - **bio.tools** entry `orevid` (summary, EDAM function and topics, links, downloads, documentation, v2 preprint,
    credit with ORCID), filled from a JSON prepared from this repository.
  - **Hugging Face**: OREVID-bench published as a dataset, `orevid/orevid-bench` (the v1.2.0 zip + `tasks/`,
    `disease_areas/`, `bench.py`, `LEADERBOARD.md` and a dataset card; CC BY-NC-SA 4.0). Uploaded with a
    fine-grained token limited to the `orevid` organization, logged in on the work PC (`~/.venvs/hf`). On the
    external exam the best value is OREVID without popularity (0.364), not production (0.359): the card says so.
  - **Bluesky** account created; Bluesky suspended it the same day for following 123 accounts at once with no
    posts (anti-spam), restored after an appeal. Follow at most 5–10 accounts a day; a short introduction post
    first; the announcement thread (6 posts, drafted) on launch day, pinned to the profile.
  - **Awesome lists**: pull requests opened from `quimqu` forks (branch `add-orevid`), each saying the author
    is the contributor: [everycure-org/awesome-ai-drug-repurposing#11](https://github.com/everycure-org/awesome-ai-drug-repurposing/pull/11)
    (Graph ML Methods), [AstraZeneca/awesome-drug-discovery-knowledge-graphs#12](https://github.com/AstraZeneca/awesome-drug-discovery-knowledge-graphs/pull/12)
    (new «Drug Repurposing» section in `chapters/applications.md`; so far only the maintainer's own PRs were merged
    there) and [yboulaamane/awesome-drug-discovery#28](https://github.com/yboulaamane/awesome-drug-discovery/pull/28)
    (Network Pharmacology). Answer any review comment there. Set aside:
    `robert-haas/awesome-biomedical-knowledge-graphs` (downloadable graphs only, no outside contributions); its
    survey gave the three candidates N1–N3 of `docs/IMPROVEMENT_PLAN.md`.
  - **E-mails sent** (22-10, from contact@orevid.org) to ReDO (info@anticancerfund.org, for Pan Pantziarka),
    Cures Within Reach (research@), CURE ID (CDRC@c-path.org), EURORDIS / REMEDi4ALL (Claudia Fuchs) and Patrick
    Aloy (Bioteque, IRB; in Catalan). **Aloy answered the same day** (his address: paloy@irbbarcelona.org): the
    embeddings may be used if the article is cited and they are not modified; he suggested BQsupports; the IRB
    websites were offline (UB/PCB security). Terms and caveats in `docs/IMPROVEMENT_PLAN.md`, N3.
  - **Open Targets community forum**: first post (announcement) held for moderation, the account «temporarily on
    hold» (automatic: new account with many links); a short note was sent to the moderators. Do not repost.
  - **Bluesky**: introduction post published (22-10).
  - Planned: launch day with the r/bioinformatics post (to be drafted) and the Bluesky thread together;
    **BOSC 2027** (Open Bioinformatics Foundation, within ISMB; abstracts usually due around April).

## Night of 22–23 October 2026 (home PC) — in progress, updated during the night

Registered first (`docs/IMPROVEMENT_PLAN.md`, «Registered 22-10-2026, before any run», commit `389458a`), then run:
- **R1b adopted (P48)**: hub genes (≥ 500 diseases) no longer count for the shared-genes line; R1a rejected
  (`docs/R1_R2.md`). Code changed in `web/build.py` and the evaluation; website not rebuilt yet.
- **R2 rejected (P49)**: removes an approved pair (sulfanilamide, vulvovaginal candidiasis).
- **N1 positive (P50)**: shared genes contain curated mechanism proteins (DrugMechDB) in 25.1 % of pairs vs 16.9 %
  for random genes of the same disease (`docs/DRUGMECHDB.md`).
- **N2 (a) above random** on drugs with no prior use (`docs/PROXIMITY.md`); (b) the LambdaRank-feature test is running.
- **N3 not run**: Bioteque still offline. **N4** (drug status audit): FDA labels downloading, code in
  `sources/drug_status.py`.
- Fixed: trials of a candidate counted over all the nodes of the same drug (belimumab in systemic sclerosis).

## Next steps

**On a computer without `data/`** (the work PC): posting the preprint (step 1), the biomarker table (3: editing the CSV; checking it needs `data/`) and
anything that only reads the benchmark prediction files (`benchmark/`, release zip). **Needs the home PC**
(`data/`, ~5 GB, and WSL for TxGNN): the GitHub release (5: the RotatE and train-only prediction files are only there), new evaluations, `evaluation/evidence_ci.py`, `web.build` (also `--patch`)
and the website deploy. The Claude Doc of the preprint is the source of truth for its text; re-export to
`docs/preprint/` after editing it (`docs/REVIEW_RESPONSE.md`, «After the pending tasks»).


1. **Preprints.org** (preprints-238326, online 20-10, DOI 10.20944/preprints202610.1918.v1): DOI added to README,
   website and `CITATION.cff`, and to ORCID (20-10; the Zenodo records were already there through DataCite). New versions go through the
   User Dashboard («My Preprints»; keep the same server; each version gets its own DOI, `.v2`).
   **v2 submitted on 20-10 at 12:24 (follow-up version, new ID preprints-239415), «Pending Check»**: files
   `docs/preprint/OREVID_preprint_v2.docx` (A4, continuous line numbers) and `OREVID_preprint_v2.pdf` (made from it
   with Word), graphical abstract = the new `figures/figure1.png`; abstract updated (temporal-validation sentence);
   the `supplementary.docx` Preprints.org had carried over from v1 was removed (the Word holds Tables S1–S5).
   **v2 online on 22-10** (DOI `…1918.v2`); README, website and `CITATION.cff` point to it (see 22-10). Left: check
   that v2 appears in ORCID (Crossref permission), or add it by hand.
   Then choose a journal (independent authors are accepted by journals).
   **Journal plan (14-10):** first PCI Mathematical and Computational Biology (free, open reviews; then Peer
   Community Journal, free), Briefings in Bioinformatics as the fallback; BOSC 2027 abstract (talk or poster);
   JOSS for the software after April 2027
   (repository public for 6 months). **Outreach, 20-10:** e-mails sent to the TxGNN authors (Marinka Zitnik,
   marinka@hms.harvard.edu, cc Kexin Huang) and to Every Cure (info@everycure.org), with the v1 link; both say OREVID
   is non-profit and open. David Ochoa: a reply on opentargets/issues#4573 (preprint link + the 26.12 check) was
   drafted. Answer any reply mentioning that v2 adds the temporal validation and RotatE.
   **PCI preparation: done on 20-10.** PCI does not accept AI-generated illustrations («Authors must not generate
   illustrations using AI», PCI code of conduct): the author read and ran `evaluation/figures.py` on the work PC
   (with the Spyder venv's Python, `~/.venvs/spyder`, which has matplotlib) → `docs/preprint/figures/figure{1..4}.png|.pdf`;
   the four charts of the Claude Doc were replaced by these images; Methods says the figures come from that script;
   Conflicts of Interest uses PCI's wording (no financial conflict of interest with the content of the article, and
   no non-financial one); the Word export has continuous line numbers. Post v2 on Preprints.org (same server)
   before submitting to PCI. Figure 3 was fixed on 21-10 (see that day): put the new image in the Claude Doc.
1b. **Improvement plan**: done (I1–I6, P42–P47). Left over: optional TxGNN on the temporal cutoffs;
   pathway overlap with its own null (P46) as a new pre-registered test. **Next candidates, not registered yet**
   (`docs/IMPROVEMENT_PLAN.md`, «Next candidates»): N4 drug status audit against current FDA / EMA / DrugCentral
   data (display only; fixes the approval errors of item 9; its checks are already fixed in the plan),
   N1 DrugMechDB as an evaluation set for the evidence shown,
   N2 network proximity with NeDRex as an external reference, N3 Bioteque embeddings without drug–disease edges.
   Each needs its own section with fixed criteria, committed before running (home PC).
2. Shared genes (P38): the line works as a selection of drugs used in comparable diseases, not as biology. The
   three pre-registered replacements failed (P46); pathway overlap came closest and would need its own null
   (permuting pathway sets) in a new pre-registered test.
3. Biomarker table: applied on 14-10 (P41); Gram-stain and class rules rejected (`docs/BIOMARKERS.md`). New rows
   go through `py -3 -m evaluation.biomarker_check` and `py -3 -m web.build --patch` as before.
4. Model: production stays (P33). Optional: a pre-specified test of the LightGBM + logistic ensemble on new uses.
5. Release v1.2.0: done 21-10 (zip attached, Zenodo checked). Next release: same steps, with `gh` from this PC.
6. ORCID record: done (preprint added by hand; both Zenodo concept DOIs added automatically by DataCite, with their versions; website links and keywords filled in).
7. Collect expert feedback from the website (GitHub issues) as new validation labels.
8. Open problems: drugs with no prior indication (Precision@K 0–2 %); wrong gene–disease
   associations and comorbidities that read as shared biology (lupus–endocarditis, RA–Graves).
9. **Corrections from a review of the visited pages (22-10, work PC; nothing changed on the website).** The eight
   disease pages visitors opened (GoatCounter) were read from the published `gh-pages` data and checked by hand:
   candidiasis (35873, also oral 83958 and vulvovaginal 83986), multiple sclerosis (37722; GoatCounter shows it
   with the home-page title), angiosarcoma (33713), systemic sclerosis (33225), Angelman syndrome (29997),
   Castleman disease (32631), Rosai-Dorfman disease (33528), POEMS syndrome (33287). Mostly sensible (angiosarcoma
   is the best example: propranolol first, with trials and papers). To do from the home PC; rule changes follow the
   usual validation (DECISIONS «Rules»), data fixes and display-only changes are labelled as such:
   - **Approval classifier (`sources/ot_approvals.py`), false «unconfirmed»** (→ plan N4): ublituximab for MS (FDA-approved
     12-2022), butoconazole and terconazole for candidiasis (FDA-approved for vulvovaginal candidiasis). The page
     says «do not treat them as treatments»: wrong for these. Check why a single-drug label was not found.
   - **False «approved»** (→ plan N4): griseofulvin for candidiasis (Open Targets; its label says it is *not* effective against
     Candida: a negation read as an indication); benzoic acid for systemic sclerosis (DrugBank/PrimeKG, probably a
     combination product); olaratumab shown as approved (withdrawn 2019; appears as an angiosarcoma extension);
     zinc oxide as a cutaneous-candidiasis drug (comes from a miconazole/zinc oxide combination).
   - **Inconsistency:** systemic sclerosis lists belimumab among phase 3 trials, but its candidate row says
     «trials: 0». Check how `n_trials` and the trials section are computed (exact disease vs descendants, sources).
   - **Shared-genes candidates from unrelated diseases:** MS gets insulin glargine / lispro (via type 1 diabetes),
     T-cell lymphoma drugs (forodesine, pralatrexate, romidepsin, mogamulizumab) and mesalazine; systemic sclerosis
     gets theophylline and pirbuterol (asthma); Castleman (only 3 genes: IL6, JAK1, NOTCH1) gets paclitaxel via
     ovarian cancers. Candidate rules to pre-register and test: a minimum disease gene-set size for the line;
     requiring the source disease to share more than the generic immune genes. Until then the page's caveat stands.
   - **Infectious diseases:** extensions through a parent term mix pathogens (vulvovaginal candidiasis gets
     metronidazole, tinidazole, secnidazole and sulfonamides via «vulvovaginitis»; antibiotics even favour
     candidiasis); the «targets» section shows human targets (glucocorticoid receptor, tubulins) that make no
     sense for a pathogen. Options: restrict parent-term extensions to the same pathogen class (antifungal /
     antibacterial / antiviral / antiparasitic, by ATC), and hide or relabel the targets section for infections.
     The earlier filter «no anti-infectives in non-infectious diseases» failed validation; this is a different
     rule and needs its own test.
   - **Safety knowledge missing:** in Angelman syndrome carbamazepine, phenytoin and vigabatrin are reported to
     aggravate seizures, yet they appear without a warning (vigabatrin is the second candidate). Eculizumab reaches
     MS through neuromyelitis optica, a distinct disease where drugs do not transfer. Consider a curated
     «known to worsen» list (with sources) for the risk column; display only, labelled.
   - **Gaps:** POEMS syndrome shows nothing (no genes, no link to plasma-cell neoplasms, where its real treatments
     come from); chronic mucocutaneous candidiasis shows nothing although azoles are used.

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
