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
| N4 | **Drug status audit** against current regulatory sources (data quality, display only) | the review of 22-10 (STATUS, next steps 9) found wrong «approved» and «unconfirmed» labels and a withdrawn drug shown as approved |

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

**N4. Drug status audit** (added 22-10-2026). Not a model change: it fixes what the website says about each drug
and about each «approved» pair, so it is validated by manual checks, not by hit rates (like I3).

*Where the errors found on 22-10 come from* (STATUS, next steps 9):
- **Extraction errors upstream:** griseofulvin «approved» for candidiasis comes from Open Targets / ChEMBL label
  mining, which read «not effective against … candidiasis» as an indication (negation); zinc oxide for candidiasis
  comes from a PrimeKG indication of the miconazole / zinc oxide combination given to each ingredient (as lincomycin
  for MS); benzoic acid for systemic sclerosis is probably a name-mapping error for potassium *aminobenzoate*.
- **Our classifier too strict:** ublituximab (MS), butoconazole and terconazole (candidiasis) are approved and
  Open Targets has them right; `sources/ot_approvals.py` marks them «unconfirmed», probably because the labels use
  acronyms (RMS, VVC) and the rule rejects acronyms of ≤ 3 letters (written for the PMDA «ad» = ADHD error).
- **Outdated status:** olaratumab shown as approved (withdrawn 2019); PrimeKG's drug status comes from DrugBank
  (~2021) and does not record withdrawals well; newer approvals are missing.
- Check also the website label «Approved (DrugBank)» for PrimeKG indications: PrimeKG takes drug–disease edges
  mainly from DrugCentral; correct the label if so.

*Sources.* Already downloaded for I3 (`data/raw/practical/`): FDA Orange Book (marketing status, discontinued
products), FDA Purple Book (licensed / discontinued biologics), EMA medicines table (authorised / withdrawn /
refused / suspended, and the therapeutic indication text). To add: **DrugCentral** (current dump; FDA / EMA / PMDA
approvals with dates, label indications mapped to disease codes; CC BY-SA 4.0: check the licence terms before
publishing derived fields) and **current FDA labels** (DailyMed / openFDA «indications and usage» section).
Optional: AEMPS CIMA (Spain).

*What to build.*
1. Per drug shown on the website: current US and EU status (marketed, discontinued, withdrawn with year and,
   when given, reason; authorised / withdrawn / refused / suspended in the EU), with the source and download date.
2. Per «approved» pair: confirmed if a current FDA label or the EMA indication text names the disease (or a
   synonym, via the graph's disease names and the MONDO / UMLS cross-references), **and** the sentence is not
   negated («not effective», «not indicated», «has not been established»). Otherwise flagged «not found in a
   current label».
3. With 2 in place, revisit `ot_approvals.py`: accept a short acronym when a current label confirms the full
   disease name; reject label mentions that are negated.
4. Combination products: an indication that only a combination product carries is shown on the combination,
   not on each ingredient.

*Checks before release (fixed now).* Report the match rate of each source (all drugs, approved drugs). Manual
check of 50 random «approved» pairs and of 50 random pairs that change category (confirmed ↔ flagged) against the
source pages, written to a new `docs/DRUG_STATUS.md`; the release goes ahead only if ≤ 5 % of the checked pairs are
wrong after the change and none of the cases of 22-10 remains wrong. The validated definitions used in the
evaluations (PrimeKG indications) do not change; if the audit removes PrimeKG indications, the effect on the
benchmark is reported, not silently applied.

All four need `data/` and run on the home PC.

---

## Registered 22-10-2026, before any run: R1, R2, N1, N2 (and the state of N3, N4)

Written and pushed before running any of them (the commit time is the proof). Unless said otherwise the test unit
is the one of `evaluation/evidence_ci.py`: the 25 highest-ranked candidates per test disease, sections «direct
extensions», «shared genes» and «hidden»; benchmark = zero-shot seeds 1–5, areas = nine areas (seed 1), external =
Open Targets approvals absent from PrimeKG (training split zero-shot seed 1); 95 % intervals resample diseases.
The **P29 rule** reused below: a filter that removes a subset from a shown section is adopted only if the removed
candidates hit < 5 % on both the benchmark and the areas.

**R1. Shared-genes candidates from unrelated diseases** (STATUS, next steps 9: insulin for MS through type 1
diabetes, paclitaxel for Castleman disease through ovarian cancers). Two variants, fixed now:
- **R1a, small gene sets:** the shared-genes line does not apply to a test disease with fewer than 10 associated
  genes (production gene sets, `GENES_EXTRA`); its shared-genes candidates become hidden.
- **R1b, hub genes:** genes associated with ≥ 500 diseases in the training graph do not count; the line needs ≥ 3
  shared non-hub genes and the hypergeometric p < 10⁻⁵ recomputed over the non-hub genes with at least one disease.
Each variant is adopted if (i) it meets the P29 rule and (ii) on the external exam the removed candidates hit no more
than the candidates the section keeps. If both pass, both are applied together only if the combination also passes;
otherwise the one that removes more candidates. Reported in `docs/R1_R2.md`.

**R2. Infectious diseases: anti-infectives of the wrong pathogen class.** Disease class from the ontology ancestors
(`sources/hierarchy.py`), by the MONDO labels «bacterial infectious disease», «fungal infectious disease», «viral
infectious disease» and «parasitic infection» (a label not found in the graph is reported); a disease with more than
one class is left alone. Drug classes from the DrugBank ATC texts (`atc_1`–`atc_4`, lower case): antibacterial
(«antibacterial», «antibiotic», «antimycobacterial»), antifungal («antimycotic», «antifungal»), antiviral («antiviral»),
antiparasitic («antiprotozoal», «anthelmintic», «antimalarial», «ectoparasiticide»); a drug belongs to every class
its texts name (so fluconazole, listed with an antibacterial combination, stays in the antifungal class). Rule: in
a disease of one class, a direct-extension or shared-genes candidate that belongs to some anti-infective class but
not to the disease's class becomes hidden; other drugs are untouched. Adopted only if (i) it meets the P29 rule (a set
with fewer than 20 removed candidates counts as uninformative, and then the external exam must also show the
removed candidates below 5 %), and (ii) as the biomarker rule (P37), **no removed pair is an approved indication**
of that disease in the full graph (PrimeKG indication or confirmed Open Targets approval); the removed pairs are
listed for a manual check. Expected conflict, stated now: trimethoprim–sulfamethoxazole for pneumocystosis
(*Pneumocystis* is a fungus); if it or any other approved pair is removed, the rule is rejected, as the Gram-stain
rule was. Display only, independent of the rule: for diseases of one infection class, the «targets» section says
that it lists human proteins.

**N1. DrugMechDB as an evaluation of the shared-genes evidence** (`indication_paths.json`, CC0, from
github.com/SuLab/DrugMechDB). Mapping: drugs by DrugBank id, diseases MESH → MONDO (SSSOM exact matches) → graph
node, proteins UniProt → gene (HGNC); coverage reported before any result. For each covered pair: the mechanism
genes *M* = proteins of the curated path except the drug's direct targets (the nodes linked to the drug in the
path); the evidence genes *S* = the disease's genes shared with the diseases where the drug is indicated (the disease
itself and its twins excluded) under the production rule (≥ 3 genes, p < 10⁻⁵, `GENES_EXTRA`), i.e. what the
shared-genes line would show. Pairs with empty *M* or *S* are counted and left out. Observed: the share of pairs with
*S* ∩ *M* ≠ ∅. Null: for each pair, |*S*| genes drawn at random from the same disease's genes, 1,000 draws; the
pooled share per draw gives the null distribution. Outcome: if the observed share is above the null's 99th
percentile, the Method page may say that the shared genes overlap curated mechanisms more than random genes of the
same disease; otherwise the P38 wording stays. Reported in `docs/DRUGMECHDB.md`.

**N2. Network proximity** (Guney et al. 2016 closest distance). Network: PrimeKG `protein_protein` ∪ STRING (score ≥
700) among genes. Disease genes: `GENES_EXTRA` of the training graph in the network (diseases with < 5 such genes get
no score). Drug targets: `drug_protein:target` in the network. Score: mean over the drug's targets of the shortest
distance to any disease gene, z-scored against **50** random disease gene sets of the same size and degree
distribution (genes binned by degree, ~100 genes per bin; 50 instead of 1,000 for compute); proximity = −z.
Drugs without a target in the network get the lowest score. (a) Standalone, as I4: benchmark scenario C (positives
without a training indication; the other positives removed from the candidates) and the temporal 2015 set (primary
positives of drugs with no indication at the cutoff), candidates restricted to drugs with a target in the network;
Recall@100 and Precision@K against random (100 / candidates) and against OREVID's ranking of the same candidates.
(b) As a LambdaRank feature: run only if (a) is above random on both the benchmark and the temporal set; adopted only
if Recall@100 on the no-prior-use subset improves with a 95 % interval excluding 0 on both the benchmark and the
areas, and overall zero-shot Recall@100 does not fall by more than 0.01. **NeDRex:** its knowledge base includes
drug–indication edges (DrugBank, DrugCentral); if its API cannot exclude them per split, it is reported as not
comparable and not run. Reported in `docs/PROXIMITY.md`.

**N3. Bioteque:** the IRB websites were still offline on 22-10 at night (no download possible): not run tonight.

**N4. Drug status audit:** registered above (sources, what to build, checks before release). Tonight: as much as
the sources allow; the release gate (≤ 5 % wrong in the two manual checks of 50 pairs, none of the 22-10 cases
still wrong) applies before anything reaches the website.

*Added 23-10-2026, after the first run (P53, not released), before any second run.* The second version follows
`docs/DRUG_STATUS.md` («What a second version needs»): disease matching with MONDO / UMLS synonyms, plurals and all
descendants of the label term; only indication sentences, with qualifiers such as the diabetes type; drugs matched by
UNII / RxNorm ingredient. Release in parts, each with its own check and the same gate (≤ 5 % wrong, none of the
22-10 cases wrong): first (a) «negated» pairs (all 62 reviewed), (b) withdrawn / revoked / refused EU medicines and
discontinued-only US drugs, (c) «unconfirmed» Open Targets approvals that a current label confirms (ublituximab,
butoconazole, terconazole); the full «not found in a current label» flag only after its own check. **The check is
done against the source pages**, by a person or by an agent that opens them and quotes the sentence it relies on,
with the URL kept for every verdict; a check from memory, without the pages, does not count (the first run's did not).

---

## Registered 23-10-2026, before any run: G1 (curated gene–disease validity, roadmap item 1)

Written on the work PC and pushed before any GenCC file was downloaded or any number computed (the commit time is the
proof). To run on the home PC. Same test unit, sets, intervals and P29 rule as «Registered 22-10-2026» above; the
production rule today is the one after R1b (P48: hub genes, ≥ 500 diseases, do not count).

**Question.** Does restricting the shared-genes line to gene–disease links that curators have classified as valid
remove wrong candidates without losing right ones?

**Source.** GenCC (https://thegencc.org), the submissions file (`gencc-submissions.tsv`), downloaded once into
`data/raw/gencc/` with a `VERSION` file (download date, file size, licence text as published on the site). GenCC
already gathers ClinGen, Genomics England PanelApp, Orphanet, Gene2Phenotype and other submitters, so ClinGen and
PanelApp are not downloaded separately. If the licence does not allow publishing derived fields, the rule may still
be evaluated and applied, but no GenCC classification is shown on the website.

**Mapping, fixed now.**
- Genes: GenCC `gene_curie` (HGNC) → NCBI Gene through the HGNC file already in `data/raw/hgnc/` (graph genes are NCBI
  keys, `sources/ids.py`).
- Diseases: GenCC `disease_curie` when it is MONDO → the graph disease node whose MONDO group contains it
  (`sources.ids.primekg_keys`); OMIM / Orphanet curies → MONDO through `sources.ids.xref_to_mondo()`. No propagation to
  ancestors or descendants.
- Classification per (gene, disease) pair, over all its submissions:
  **supported** = at least one Definitive, Strong or Moderate, and no ClinGen submission Disputed or Refuted;
  **contradicted** = no Definitive / Strong / Moderate and at least one Disputed, Refuted or «No Known Disease
  Relationship»; **limited** = everything else (Limited, Animal model only, Supportive).
- A graph disease is **covered** if it has ≥ 3 supported genes after mapping.
- **Report before any result:** GenCC pairs read and mapped (genes, diseases), covered graph diseases, covered test
  diseases in the benchmark, areas and external exam, and how many of the graph's gene–disease edges are supported /
  limited / contradicted, by edge source (`disease_protein`, Open Targets genetic, animal, expression).

**Variants, fixed now** (they change only the shared-genes line, in the evaluation and on the website, as R1b; the
model's features and its ranking are not touched):
- **G1a, contradicted links out:** a gene–disease edge classified contradicted does not count for the line, for any
  disease.
- **G1b, curated genes for the test disease:** for a covered test disease, only its supported genes count on its side
  of the comparison (the other disease keeps its gene set); the hypergeometric test is recomputed with that gene set
  (same population N of non-hub genes with at least one disease, ≥ 3 shared genes, p < 10⁻⁵). Diseases not covered are
  unchanged.
- **G1c, curated genes on both sides:** as G1b, and the disease the candidate drug is indicated for also uses only
  its supported genes when that disease is covered.

**Decision rules, fixed now.** A variant changes the section in two ways: candidates **removed** (shared-genes
candidates that no longer qualify become hidden) and candidates **added** (hidden candidates that now qualify). It is
adopted only if all of these hold:
1. P29 rule on the removed candidates: they hit < 5 % on both the benchmark and the areas;
2. on the external exam the removed candidates hit no more than the candidates the section keeps;
3. if a set has ≥ 20 added candidates, their hit rate on that set is not below the current section's (production,
   after R1b) beyond the 95 % interval;
4. the section's hit rate after the change is not lower than today's on the external exam beyond the interval.

A set with fewer than 20 removed candidates is uninformative for rule 1; if both the benchmark and the areas are
uninformative, the external exam must show the removed candidates below 5 %, otherwise the variant is not adopted
(«not enough data» is reported as such). **Choice:** among the variants that pass, G1c if it passes, else G1b, else
G1a (the more curated, the better, if the data allow it); a combination of G1a with G1b or G1c is adopted only if
the combination itself passes. If none passes, the website keeps today's rule.

**Reported, not deciding:** the degree-preserving permutation null of `evaluation/evidence_ci.py` for the adopted
variant (does curation make the line beat the null, which P38 says it does not today?); the temporal C2 set (I1) with
the same table; and the pages reviewed on 22-10 (STATUS, next steps 9). Expected now, before running: GenCC covers
mostly Mendelian and rare diseases, so multiple sclerosis (insulin through type 1 diabetes) and Castleman disease
(paclitaxel through IL6, JAK1, NOTCH1) will probably not change; for those the roadmap's item 4 (drug-target MR) is
the better route.

**Display only, independent of the rule** (if the licence allows): on each disease page, the GenCC classification
next to the genes that have one, with a link to GenCC; a gene without one says nothing.

**Implementation.** `sources/gencc.py` (download, mapping, classification → `data/ext/gencc.parquet`) and
`evaluation/g1_gencc.py` (the variants, reusing `evaluation.evidence_ci`) → `docs/G1_GENCC.md`; decision in DECISIONS
P54; if adopted, the same rule in `web/build.py`, the Method page (three languages) and the README.

---

## Registered 23-10-2026, before any run: D1 (drugs that may worsen a disease; roadmap item 2)

Written on the home PC and pushed before OnSIDES is downloaded or any number computed.

**Question.** Can OREVID show, per disease, approved drugs that may **worsen** it, with a reason and a source, and are
those warnings right more often than a naive list?

**Sources.** (i) PrimeKG contraindications (`contraindication`, DrugBank) and the validated risk rule (P15:
contraindicated in more neighbouring or similar diseases than indicated, plus checkpoint inhibitors in immune
diseases); (ii) **OnSIDES** v3.2.1 (tatonetti-lab, adverse events extracted from current US / EU / UK / JP labels;
data CC BY 4.0, attribution by citing the paper), downloaded once into `data/raw/onsides/` with a `VERSION` file.
Curated drugs-to-avoid lists (porphyria, Brugada, long QT, G6PD, myasthenia gravis) are **not** part of D1: each has
its own terms; they get their own item after their licences are checked.

**Mapping, fixed now.** OnSIDES ingredients → graph drugs by RxNorm / name through the active moiety (as N4); OnSIDES
adverse-event terms (MedDRA) → graph diseases through MONDO cross-references to MedDRA or UMLS in
`data/raw/mondo/mondo.sssom.tsv`, exact matches only; if that file has none, by exact normalised name of the disease
or of its Open Targets exact synonyms. Report before any result: OnSIDES pairs read, drugs and terms mapped, graph
diseases reached.

**Signals, fixed now.** For a disease *d* and an approved drug *f*: **S1** the P15 risk rule (already on the website);
**S2** an adverse event of *f* in a current label (OnSIDES) maps to *d* or to a direct broader term of *d*; **S3** both.

**Test.** Diseases with ≥ 3 PrimeKG contraindications; for each, its own contraindication edges are hidden from the
signals (S1 uses other diseases' contraindications only; S2 uses labels, not PrimeKG). Positives: the hidden
contraindicated drugs. Candidates: approved drugs not indicated for *d*. Baseline: «contraindicated somewhere», drugs
ranked by the number of diseases they are contraindicated for, flagging as many drugs per disease as the signal does.
Diseases are split by node id parity into a development half (even) and a confirmation half (odd), the latter
reported once. Precision = share of flagged pairs that are hidden contraindications; disease-level bootstrap.

**Decision rule.** A signal is shown in the new section «May worsen this disease» only if, on the confirmation half,
its precision is above the baseline's with a 95 % interval of the difference excluding 0, and it flags at least 200
pairs. Then, before release, **a manual check of 50 random flagged pairs (not contraindications already in PrimeKG)
against the source pages** (the label section on DailyMed or the EMA product page, opened, URL and sentence kept):
released only if ≤ 5 % are wrong. Reported, not deciding: the Angelman syndrome page of the 22-10 review
(carbamazepine, phenytoin, vigabatrin), the share of flags explained by confounding by indication (the drug is
indicated for a neighbour of *d*), and recall of all PrimeKG contraindications.

**Display (if released).** A section per disease page listing the flagged drugs with the reason («contraindicated in
related diseases», «adverse event in the current label: …») and the source; never a score; research use only.
Implementation: `sources/onsides.py`, `evaluation/d1_worsen.py` → `docs/D1_WORSEN.md`; DECISIONS P56.

---

## What to update after each item

DECISIONS (new entries P41…, one per item), STATUS, the Claude Doc of the preprint (re-export to
`docs/preprint/`; new version on Preprints.org), the website (Method page and help texts in three languages)
and, for I1/I5, OREVID-bench (leaderboard and a release with the new prediction files).
