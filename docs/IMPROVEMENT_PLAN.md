# Improvement plan: more useful and more credible predictions (written 14-10-2026, before any result)

This plan fixes, before running anything, what will be built, how it will be evaluated and which result
decides each change. Its git commit date is the registration date; any deviation will be reported as such
(same practice as `docs/PROSPECTIVE_PLAN.md`). Everything here needs `data/` and runs on the home PC.

Starting point (see STATUS «Results to remember», DECISIONS P27–P40): the ranking model does most of the
work (top 25 without rules 4.2 % on the external exam; rules without the model 1.0 %); direct extensions are
largely rediscovery; shared genes do not beat a permutation null (P38); mechanistic evidence alone is like no
evidence; drugs with no prior indication are not found (Precision@K 0–2 %); and the external exam measures
mostly old approvals missing from PrimeKG, not discovery.

| # | Work | Why | Effort | Order |
|---|---|---|---|---|
| I1 | Temporal validation | answers «does it find new uses?», the main open question of reviewers | medium | 1 |
| I2 | Calibrated probability per candidate | turns a score into «X in 100 similar candidates became real uses» | low–medium | 2 (needs I1) |
| I3 | Practical repurposing information | generic / patent status, route, «already studied?», «underexplored» flag | low | 1 (independent) |
| I4 | Expression-signature reversal (LINCS L1000) | an independent signal for drugs with no prior indication | high, uncertain | 4 |
| I5 | Knowledge-graph embedding baselines | completes the comparison asked for by the review | medium (GPU) | 3 |
| I6 | A replacement for shared genes that beats the null | a «biological» evidence line that really is one | medium, uncertain | 4 |

General rules (as in DECISIONS «Rules»): criteria fixed here before running; every result reported, including
null ones; intervals by disease-level bootstrap (`benchmark/cluster_ci.py`, `evaluation/evidence_ci.py`);
nothing reaches the website without passing its rule, except items marked «display only».

---

## I1. Temporal validation

**Question.** Trained only on what was known before a cutoff year, does OREVID rank the drug–disease pairs
that reached phase 3 or approval after that year above the popularity baseline and above the other drugs?

**Dates available.** `data/prep/attributes.parquet`, relation `clinical`: highest phase per pair
(`display_relation`: approved / phase 3 / 2 / 1) and `year` = first trial start (`sources/opentargets.py`
`clinical()`); relation `ot_assoc`: first year with evidence per gene–disease association. PrimeKG edges
(indications from DrugBank, DisGeNET associations…) are undated.

**Cutoffs.** Development cutoff **C1 = 2012** (anything tuned is tuned here only); confirmatory cutoff
**C2 = 2015** (reported once, never used to choose). Outcome window: from the cutoff to the Open Targets
26.06 data.

**Training graph at cutoff C.**
1. Remove every drug–disease edge (indication, contraindication, off-label, clinical) whose pair has a dated
   Open Targets clinical record with first year ≥ C and no dated record before C.
2. Remove Open Targets gene–disease associations whose first evidence year is ≥ C.
3. PrimeKG edges without a date stay (they cannot be dated). Report the share of PrimeKG indications that
   have no Open Targets date, as the size of the residual leakage.

**Positives (primary endpoint).** Pairs whose highest phase is 3 or approved, with first trial year ≥ C, no
dated clinical record before C, and no PrimeKG drug–disease edge for the pair. **Secondary:** confirmed
approvals only (`sources/ot_approvals.py`, single-drug label or agency report); and phase ≥ 2.
**Candidates:** all drugs except those with any drug–disease relation to the disease dated before C or in
PrimeKG. Pairs without any year are excluded from positives and counted.

**Models.** OREVID production configuration (`extra_no_ctd`, 5 masked rounds, LambdaRank) retrained on the
cutoff graph; popularity baseline computed on the cutoff graph; logistic regression without popularity.
Sensitivity: the same with `extra_no_ot` (no Open Targets features), which removes the dated-but-imperfect
literature component. TxGNN is optional (one run per cutoff if time allows).

**Metrics.** Recall@100 per disease and the hit rate of the top 25 by website section (website rules applied
on the cutoff graph), each with disease-level 95 % intervals; non-obvious subset (no indication within two
ontology steps before C).

**Decision rules.**
- The preprint may say «OREVID ranks later phase 3 / approved uses above the popularity baseline» only if
  the paired per-disease difference has a 95 % interval excluding 0 at C2. Otherwise it says the opposite.
- Any model change aimed at new uses (features, objective, ensemble — e.g. the LightGBM + logistic ensemble of
  P33) is chosen on C1 and adopted only if, at C2, its paired improvement in Recall@100 has a 95 % interval
  excluding 0 and zero-shot Recall@100 does not fall by more than 0.02.

**Implementation.** New `evaluation/temporal.py` (masking function analogous to `split_edges` returning the
0/1/2 label array; reuse `masked_rows`, `fit_ranker`, `score_paths`); outputs `docs/TEMPORAL.md` and
`data/eval/temporal_C{year}.npz`; disease-level intervals with the code of `benchmark/cluster_ci.py`.

**Preprint.** New Results subsection «Temporal validation»; replaces «the external exam is not prospective»
in Limitations by what the temporal test shows and what it still cannot (undated PrimeKG edges).

---

## I2. Calibrated probability per candidate

**Goal.** Next to each candidate, an estimate such as «about 6 in 100 candidates like this reached phase 3
or approval within ~10 years», with an interval, instead of a section-wide hit rate.

**Unit and features.** Every pair in the top 25 of a disease (shown or hidden). Features fixed now: website
section (extension / shared genes / hidden), rank bucket (1–5, 6–10, 11–25), disease gene-set size bin
(< 100, 100–499, ≥ 500), hard flag (near-identical neighbour), ≥ 3 specific shared genes.

**Fit and test.** Logistic regression on those features (no interactions), fitted on the C1 temporal set
(outcome: primary endpoint of I1); evaluated on C2: Brier score, log loss, calibration intercept and slope,
reliability diagram, all with disease-level bootstrap. **Baseline:** the section-only rates the website shows
today.

**Decision rule.** Shown on the website only if, at C2, the Brier score is at least 5 % (relative) better than
the section-only baseline and the calibration slope is within 0.8–1.25. The text states the horizon (≈ 11
years for C2) and that approval is a slow, lower-bound outcome. Otherwise the website keeps section rates.

**Implementation.** `evaluation/calibration.py` → `docs/CALIBRATION.md` and a small coefficients file used by
`web/build.py`.

---

## I3. Practical repurposing information (display only)

Does not change the ranking; validated by spot checks, not by hit rates (DECISIONS rule 3 applies: anything
uncertain is labelled as such on the page).

1. **Generic and exclusivity status.** Sources: FDA Orange Book (products, patents, exclusivity files:
   ANDA product present = generic in the US; latest patent / exclusivity expiry), FDA Purple Book
   (biologics, biosimilars), EMA medicines table (generic / biosimilar flags). Matched to graph drugs by
   active ingredient through `evaluation.twins.active_moiety` / `drug_family` and DrugBank synonyms.
   Fields: `generic_us`, `generic_eu`, `biosimilar`, `exclusivity_until`, `routes`. Check: match rate per
   source and a manual check of 50 random drugs against the source pages, written to `docs/PRACTICAL.md`.
   Check the reuse terms of each source before publishing (attribution on the legal page).
2. **«Already studied?».** Per shown pair: number of trials for the pair (Open Targets clinical, any phase,
   with the latest phase) and number of PubTator3 papers co-mentioning drug and disease (already in
   `attributes.parquet` as `pubtator_drug_disease`).
3. **«Underexplored» flag**, thresholds fixed now: shown candidate or extension with **0 trials** for the pair
   and **≤ 2 co-mention papers**. Exploratory analysis (labelled as such, outside the pre-registered
   prospective plan): outcome rates of flagged vs unflagged shown pairs at the 12- and 24-month checks.
4. Website: columns «Generic», «Studied» and a filter «underexplored / generic available»; texts in three
   languages; CSV export includes the fields.

---

## I4. Expression-signature reversal (LINCS L1000)

**Idea.** Drugs whose induced expression profile is the reverse of the disease's profile (Connectivity Map
logic). Independent of known indications, so it can reach drugs with no prior use.

**Data.** LINCS L1000 consensus signatures (GEO GSE92742 and GSE70138, level 5, or the clue.io consensus
signatures); disease signatures from CREEDS (manual and automatic, mapped to MONDO via DOID/UMLS) and, if
needed, GEO disease-vs-control series. Drugs matched to the graph by InChIKey / name. Report coverage
(drugs × diseases with both signatures) before any evaluation.

**Score.** Weighted connectivity (signed overlap of the top/bottom 150 genes), the reversal score per
drug–disease pair.

**Evaluation and rule (fixed now).** (a) Standalone, on drugs with no prior indication (novelty scenario C of
OREVID-bench, and the temporal I1 set): Recall@100 and Precision@K against random and against the current
exploratory list. (b) As an added feature in the LambdaRank model: adopted only if Recall@100 on the
no-prior-use subset improves with a 95 % interval excluding 0 on **both** the benchmark and the areas, and
overall zero-shot Recall@100 does not fall by more than 0.01. Otherwise documented as a negative result.

---

## I5. Knowledge-graph embedding baselines

**Models.** TransE, ComplEx and RotatE (PyKEEN), trained on each split's PrimeKG training edges (the same
graph TxGNN sees), scoring (drug, indication, disease). Hyperparameters: PyKEEN defaults plus a small grid
fixed now (embedding dim 128/256, learning rate 1e-3/1e-2, 100 epochs, negative sampler «basic»), chosen on
zero-shot seed 1 only, then run unchanged on seeds 2–5, the nine areas and the external exam.

**Use.** Predictions exported in OREVID-bench format (`benchmark/predictions/<model>/`), evaluated with
`bench.py` and `cluster_ci.py`; added to the leaderboard and to the preprint's Table 1. An OREVID + embedding
ensemble follows the P33 rule (selection on seeds 1–2, confirmation on seeds 3–5, areas and external).

---

## I6. A shared-biology line that beats the null

Candidates, fixed now: (a) **pathway overlap**: indicated for a disease whose genes are enriched in the same
Reactome pathways (hypergeometric over pathways, not genes); (b) **specificity-weighted gene overlap**
(inverse-document-frequency weights, genes in many diseases count less); (c) **genetic-only genes** (Open
Targets genetic association and somatic mutation evidence only).

**Rule.** A candidate line replaces «shared genes» only if, on the external exam, (1) its hit rate is above
the 99th percentile of the degree-preserving permutation null of `evaluation/evidence_ci.py`, (2) its hit
rate is not below the current 4.3 % beyond the interval, and (3) on the benchmark and areas it is not lower
than the current line. Selection among (a)–(c) on the benchmark only; confirmation on the external exam and
on the temporal C2 set. If none passes, the website keeps today's wording (P38).

---

## Next candidates (noted 22-10-2026, NOT registered)

Ideas from the survey of biomedical knowledge graphs by R. Haas
(https://github.com/robert-haas/awesome-biomedical-knowledge-graphs), checked against the sources OREVID already
uses and the rules in STATUS. **Nothing here is registered:** each item gets its own section with fixed criteria,
committed before any run (as I1–I6), and only then is evaluated. Most graphs in that survey were set aside:
Hetionet, DRKG, RTX-KG2, ROBOKOP, SPOKE, the Translator, SemMedDB, INDRA and Implicitome carry literature
«treats» edges (the answer, never a feature); Reactome is already tested (I6, P46); Monarch's curated
gene–disease links (OMIM, Orphanet) already reach the graph through Open Targets' genetic associations.

| # | Work | Why |
|---|---|---|
| N1 | **DrugMechDB** as an evaluation set for the evidence shown | first direct measure of whether the evidence paths the website shows resemble real mechanisms (the weak point P38 left) |
| N2 | **Network proximity** (drug targets ↔ disease module), with **NeDRex** as an external reference | a signal independent of known indications, aimed at drugs with no prior use |
| N3 | **Bioteque** embeddings restricted to sources without drug–disease edges | the same target, from a resource built on 150 sources |

**N1. DrugMechDB** (https://github.com/SuLab/DrugMechDB, CC0; curated drug → target → … → disease paths for
indications sampled from DrugCentral, written from DrugBank, Wikipedia and the literature). Use it **only for
evaluation**: it encodes indications, so it never enters the features or the training labels. Question: for the
pairs it covers, do the paths OREVID shows (target genes, shared genes, pathways) contain the genes / proteins
of the curated mechanism more often than a matched null (random paths of the same length and degree)? To fix at
registration: identifier mapping to the graph and its coverage (report before any result), the overlap measure,
the null, and what the Method page may say depending on the outcome.

**N2. Network proximity and NeDRex** (https://github.com/repotrial/nedrexapi_v2d). Proximity as in network
medicine: mean shortest distance from a drug's targets to the disease's genes in the STRING / PrimeKG
protein network, z-scored against degree-matched random gene sets. Two uses: (a) a standalone score and an
added LambdaRank feature on the no-prior-use subset (novelty scenario C and the temporal I1 set), with an
adoption rule like I4; (b) NeDRex's own ranking (closeness / TrustRank over its disease modules) run on
OREVID-bench as an external baseline, if its API can be restricted to each split's training data. To fix at
registration: which disease genes (all, genetic only), the network and threshold, the number of random sets,
how NeDRex's data are kept free of the test drug–disease edges (if they cannot be, it is reported as not
comparable and not run).

**N3. Bioteque** (IRB Barcelona, https://bioteque.irbbarcelona.org; Fernández-Torras et al., Nat Commun 2022).
Precomputed node embeddings from 150 sources. Only metapaths without any drug–disease or indication edge
(e.g. compound–target, gene–gene, gene–pathway) may be used; the list is fixed at registration and checked
against Bioteque's metapath descriptions. Use: drug–drug and drug–disease-gene similarity as features, with
the I4 adoption rule. To fix at registration: the metapaths and how drugs and genes are matched.

*Terms, from Patrick Aloy (e-mail of 22-10-2026):* the embeddings may be used provided the Bioteque article is cited
and the embeddings are not modified; no other restriction. Computing similarities from them leaves them unmodified;
OREVID-bench will not redistribute them (users download them from Bioteque; the code documents how). He also
pointed to **BQsupports** (https://bqsupports.irbbarcelona.org, Bioinformatics 2023,
doi:10.1093/bioinformatics/btad581), which, given a set of pairs (e.g. drug–gene interactions), finds the
metapaths / embeddings that best explain them. It may be used to choose the metapaths **only with pairs that are
not drug–disease** (e.g. drug–target interactions) or with training pairs of the development split, never with
test indications; the choice is fixed before the evaluation. On 22-10 the IRB websites (Bioteque, BQsupports) were
offline (access cut by UB/PCB for security reasons); the download waits until they are back.

All three need `data/` and run on the home PC.

---

## What to update after each item

DECISIONS (new entries P41…, one per item), STATUS, the Claude Doc of the preprint (re-export to
`docs/preprint/`; new version on Preprints.org), the website (Method page and help texts in three languages)
and, for I1/I5, OREVID-bench (leaderboard and a release with the new prediction files).
