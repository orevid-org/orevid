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

## Registered 23-10-2026, before any run: R3 (infections by pathogen class; roadmap item 8)

Written on the home PC and pushed before MONDO or ChEMBL is downloaded or any number computed. R2 (P49) failed for
two reasons this item removes: PrimeKG's ontology merges «bacterial» and «viral infectious disease», and the DrugBank
ATC *texts* mix classes (corticosteroids read as antibacterial through combination products).

**Disease class.** From the current MONDO release (`mondo.obo`, version recorded in `data/raw/mondo/VERSION`): the
MONDO ids each graph disease node carries (PrimeKG groups several), their `is_a` ancestors, and the classes
«bacterial infectious disease», «viral infectious disease», «fungal infectious disease» and «parasitic infection»
(ids reported; a label not found is reported). A node with more than one class, or none, is left alone.
Pneumocystosis (any disease whose name contains «pneumocyst») is left alone, stated now: *Pneumocystis* is a fungus
treated with antibacterial and antiparasitic drugs.

**Drug class.** From ATC **codes** in ChEMBL (`molecule.atc_classifications`, ChEMBL release recorded), mapped to
graph drugs through Open Targets' ChEMBL → DrugBank cross-references (`sources/ids.py`, salts and parents
included). Classes by code prefix, nothing else: antibacterial J01, J04; antifungal J02, D01; antiviral J05, D06BB,
S01AD; antiparasitic P01, P02, P03. Combination and mixed groups (A07A, G01, D06A, S01A except S01AD, D07C, S01C …)
give no class. A drug belongs to every class its codes give.

**Which classes count as wrong**, fixed now from standard practice (cross-class uses that are standard are allowed):
in a viral disease, antibacterial, antifungal and antiparasitic drugs; in a bacterial disease, antiviral and
antifungal drugs; in a fungal disease, antiviral and antibacterial drugs; in a parasitic disease, antiviral drugs
only (antibacterials such as doxycycline, clindamycin or metronidazole and antifungals such as amphotericin B are
standard there).

**Exemption.** A drug indicated for any disease of the same class in the graph the model sees (the training graph in
the evaluation, the full graph on the website) is never out of class (e.g. a sulfonamide with a fungal indication).

**Rule.** In a disease with one class, a direct-extension or shared-genes candidate that has an anti-infective class,
none of them the disease's, whose classes are all «wrong» for that disease class, and that is not exempt, becomes
hidden. Other drugs are untouched.

**Test and decision rule, the same as R2.** Same top 25, sections and sets as `evaluation/evidence_ci.py` (benchmark
seeds 1–5, nine areas, external exam). Adopted only if (i) the P29 rule holds (removed candidates hit < 5 % on the
benchmark and the areas; a set with fewer than 20 removed candidates counts as uninformative, and then the external
exam must also show the removed candidates below 5 %), and (ii) **no removed pair is an approved indication** of that
disease in the full graph (PrimeKG indication or confirmed Open Targets approval). Reported first: classes found,
drugs per class, test diseases per class; then every removed pair for a manual check. Expected conflicts, stated now:
sulfanilamide for vulvovaginal candidiasis (P49; exempt only if it has another fungal indication in the training
graph); any approved pair removed rejects the rule. If adopted, the website applies it with the full graph and the
Method page says so. Implementation: `sources/infections.py`, `evaluation/r3_infections.py` → `docs/R3_INFECTIONS.md`;
DECISIONS P57.

---

## Registered 23-10-2026, before any run: M1 (direction from protein QTL colocalisation; roadmap item 4, first step)

Written on the home PC and pushed before the colocalisation files are downloaded (one 100 MB part was opened to read
the columns: `betaRatioSignAverage` is ±1, the left side is always a GWAS, the method is COLOC with eCAVIAR).

**Question.** When a protein's genetically higher level goes with higher risk of a disease, do drugs that lower that
protein treat the disease more often than drugs that raise it (and the reverse)? This is the sign of a
single-instrument Mendelian randomisation, as Open Targets computes it; a full MR (several instruments, pleiotropy
tests) is a later step only if this one passes.

**Data, fixed now.** Open Targets 26.06 `colocalisation` (all 200 parts, streamed; only rows whose right side is a
`pqtl` or `eqtl` kept), `credible_set` light table (study, type, `isTransQtl`) and `study` (QTL gene, GWAS
`diseaseIds`). Kept: h4 ≥ 0.8, `betaRatioSignAverage` = ±1, **cis** QTL only. GWAS disease ids → graph disease nodes
by the Open Targets disease id of the graph node or its MONDO cross-reference (`sources/ids.py`), exact ids only, no
propagation. Per (gene, disease) the sign is kept only if all its colocalisations agree. Drug mechanisms: Open
Targets `drug_mechanism_of_action` (ChEMBL), drugs → graph by `sources/ids.py`; action types that **lower** the
target (inhibitor, antagonist, blocker, negative allosteric modulator, inverse agonist, degrader, antisense,
RNAi inhibitor, disrupting agent) or **raise** it (agonist, activator, positive allosteric modulator, opener,
stabiliser, partial agonist); others ignored. A drug with targets of both kinds for one (gene, disease) is dropped.

**Signal.** For (drug, disease): **agrees** if the drug lowers a target whose higher level raises the disease risk, or
raises one whose higher level lowers it; **opposes** otherwise. Primary: pQTL. Reported, not deciding: eQTL, and both.

**Test (a), deciding.** Positives: the drug's indications for that disease in the full graph (PrimeKG indication or
confirmed Open Targets approval). Statistic: indication rate among «agrees» pairs ÷ rate among «opposes» pairs,
95 % interval by resampling diseases. **Direction validated** only if the lower bound is > 1 and there are at least
30 indication pairs among the directed pairs. Reported: coverage first (genes, drugs, diseases, pairs, positives),
the same ratio on the Open Targets approvals absent from PrimeKG (the external exam's out-of-graph positives), and
per therapeutic area if ≥ 10 positives.

**Test (b), deciding for a risk warning.** Positives: PrimeKG contraindications. Statistic: contraindication rate
among «opposes» ÷ rate among «agrees», same interval and thresholds.

**What follows.** Nothing changes on the website from this item. If (a) passes, a display line («human genetics
agrees with the direction of this drug») and a model feature are registered separately (I4 rule); if (b) passes, a
risk warning is registered separately with the D1 release check. Implementation: `sources/coloc.py`,
`evaluation/m1_direction.py` → `docs/M1_DIRECTION.md`; DECISIONS P58.

---

## Registered 23-10-2026, before any run: P1 (measured off-targets; roadmap item 5)

Written on the home PC and pushed before any activity is downloaded (one ChEMBL API call was made to read the fields).

**Question.** Does a drug that measurably hits a disease's genes outside its known mechanism rank its new uses
better, above all for drugs with no prior use?

**Sources, fixed now.** (i) **ChEMBL 37** activities (web services) of every ChEMBL molecule that maps to a graph drug
(`sources/ids.py`, salts and parents included): `pchembl_value` ≥ 6, standard type IC50, Ki, Kd or EC50 (no HTS
«Potency»), no data-validity comment, target organism *Homo sapiens*, target type SINGLE PROTEIN, target → UniProt →
HGNC → graph gene. (ii) **IUPHAR / BPS Guide to Pharmacology** `interactions.csv` (version recorded): human targets,
ligands mapped through its ChEMBL ids, affinity (median of the reported range) ≥ 6. A drug's **off-targets** are its
measured targets minus the targets it already has in the graph (`drug_protein`). The Broad Repurposing Hub is not
used now (licence to check). Coverage reported first: drugs with an off-target, off-targets per drug, pairs.

**Feature.** For (drug, disease): log(1 + number of the disease's genes, in the split's training graph
(`GENES_EXTRA`), that are off-targets of the drug). It uses no drug–disease edge.

**(a) Gate.** Standalone, on the benchmark's no-prior-use positives (zero-shot seeds 1–5, scenario C of
`docs/NOVELTY.md`): Recall@100 of the feature (ties broken at random) against random ranking, difference with a
disease-level 95 % interval. (b) is run only if the interval excludes 0.

**(b) Decision rule, the same as N2 (P51).** The production rows and model with and without the column
(`evaluation/proximity_feature.py` pattern, paired by disease): adopted only if Recall@100 on the no-prior-use subset
improves with a 95 % interval excluding 0 on both the benchmark (zero-shot 1–5) and the areas, and overall zero-shot
Recall@100 does not fall by more than 0.01. If adopted: the temporal test (2015 cutoff) is reported before the
website model changes, and the licences (ChEMBL CC BY-SA 3.0, GtoPdb CC BY-SA 4.0) are added to the legal page;
only derived features are published. Implementation: `sources/offtargets.py`, `evaluation/p1_offtargets.py` →
`docs/P1_OFFTARGETS.md`; DECISIONS P59.

---

## Registered 23-10-2026, after M1 and before any new data: M1b (prospective confirmation of the genetic direction)

M1 (P58) did not validate the pQTL direction (21 positives, 30 required); the eQTL and pQTL + eQTL ratios (4.5 and
5.1, intervals above 1) were reported, not deciding, and were seen before this registration, so they need data the
analysis has not seen. **Test:** exactly the M1 code and filters (`evaluation/m1_direction.py`, colocalisations of
Open Targets 26.06, unchanged), with positives replaced by the drug–disease approvals that are confirmed (by
`sources/ot_approvals.py`) in a later Open Targets release (26.09 is skipped: its multiple sclerosis merge, see
README) and absent from both PrimeKG and Open Targets 26.06. Primary: pQTL + eQTL; ratio agree ÷ oppose with the
95 % interval resampling diseases. **Validated** if the lower bound is > 1 with at least 15 new positives among the
directed pairs; if there are fewer, the test waits for the next release (positives accumulate, nothing else
changes). Then, and only then, a display line («human genetics agrees with the direction of this drug») is
registered with its own release check.

---

## Registered 23-10-2026, after M1, by the author's decision: M2 (genetic direction as a model feature)

**Deviation, decided by the author after seeing M1 (P58):** the 30-positive minimum of M1 was set by the assistant
without the author; the author judges that it should not stop a signal of this size, and takes **pQTL + eQTL** as
the M1 result (indication rate 5.14 times higher when the drug's action agrees with the genetic direction, 95 %
interval 2.49–12.34, 117 indication pairs, 1,351 diseases), since pQTL alone had too few cases. Recorded here and in
DECISIONS; M1b (prospective) stays as an extra confirmation and no longer blocks the next steps.

**M2 test, fixed now before any run.** One column added to the production model: for (drug, disease), the number of
the drug's targets whose action agrees with the colocalisation direction for that disease minus the number that
oppose it (pQTL + eQTL, exactly M1's filters and mappings, expanded to the disease's twin nodes; missing when the
drug has no directed target for the disease). It uses no drug–disease edge (GWAS, QTL and mechanisms only), so the
same value is valid in every split. Same procedure as N2 (b) (`evaluation/proximity_feature.py` pattern): production
rows and model with and without the column, paired by disease, zero-shot 1–5 and the nine areas. **Adopted** only if
Recall@100 on all positives improves with a 95 % interval excluding 0 on the benchmark **or** the areas while not
falling (point estimate ≥ −0.005) on the other, and the no-prior-use subset does not fall by more than 0.01.
Reported: the same restricted to test diseases with at least one directed pair (where the column can act).

**If adopted:** the temporal test (2015) is reported, then the website model is retrained and a display line
(«human genetics agrees / disagrees with the direction of this drug») is registered with a release check of 50
shown pairs against the Open Targets pages. Implementation: `evaluation/m2_direction_feature.py` →
`docs/M2_DIRECTION_FEATURE.md`; DECISIONS P61.

---

## Registered 23-10-2026, before any change to the website: G2 (genetic direction shown as information)

M1 (with the author's decision, P58) found that drugs used for a disease agree with the QTL colocalisation direction
5.1 times more often than they oppose it; M2 (P61) found too little coverage on the benchmark to change the ranking.
G2 shows the direction as **information**, never as a score and without changing the order.

**What is shown.** On a disease page, for a listed candidate whose target has a cis-QTL colocalising with a GWAS of
that disease (M1's data and filters, pQTL + eQTL, twins included): one line per target, «Human genetics: higher
*GENE* goes with higher / lower risk; this drug lowers / raises it → agrees / disagrees», the QTL type, and a link to
the Open Targets target–disease page. With the reading of M1 on the Method page (5.1 times; based on known uses;
association, not proof). Coverage reported before release (pages and candidates with a line).

**Release check, fixed now.** 50 random shown lines, checked against Open Targets (its GraphQL API, the data behind
the platform pages, release 26.06): (1) a GWAS credible set of a study mapped to the disease colocalises (h4 ≥ 0.8)
with a cis QTL of the gene, with the stated direction (`betaRatioSignAverage`); (2) the drug's mechanism on that gene
in Open Targets / ChEMBL is of the stated kind (lowers or raises). A line is wrong if (1) or (2) fails. **Released
only if ≤ 2 of 50 are wrong**; the queries and answers are kept. Implementation: `web/build.py` (`genetics` field),
`web/public/app.js`, `evaluation/g2_check.py` → `docs/G2_GENETIC_DIRECTION.md`; DECISIONS P62.

---

## Registered 23-10-2026, evening, before any download: O1, O2, P1b, R4, X1 (roadmap items 9, 5, 12, 13)

Written on the home PC before the data are downloaded (only file listings and licence pages were read). Every
source also enters the internal graph as a layer (`layers.py`, P60) whatever its result. Following the author's view
(M1, P58), no minimum counts are set: a test passes on its interval alone, and the counts are reported beside it.

**O1. CIViC (CC0), oncology exam.** CIViC clinical evidence items, status accepted, type *Predictive*, significance
sensitivity / response, levels A–B (validated or clinical), single-drug therapies only; drugs → graph by name and
synonyms (`sources/practical.py` moieties), diseases by DOID → MONDO (SSSOM, exact) → graph node. Positives: those
drug–cancer pairs that are neither PrimeKG indications nor confirmed Open Targets approvals. For each such cancer,
the production model's ranking (full graph, `predictor.py`) of all approved drugs minus the known indications,
against popularity (number of indications). **Claim** («OREVID ranks CIViC-supported uses above popularity») if the
paired Recall@100 difference, resampling diseases, has a 95 % interval above 0. Reported: the subset without any Open
Targets clinical record for the pair.

**O2. PRISM Repurposing 19Q4 (CC BY 4.0), drug sensitivity by cancer lineage.** Secondary screen dose–response
(AUC; lower = more sensitive), per drug and cell-line lineage the median AUC; drugs → graph by name / Broad id
synonyms. Lineage → graph cancers: a table written from the cell-line metadata before any AUC is read
(`resources/prism_lineages.csv`). Test: for each lineage, among the screened drugs, do the drugs approved for a cancer
of that lineage (PrimeKG or confirmed Open Targets) have a lower median AUC than the other screened drugs? AUROC per
lineage, mean over lineages against 0.5 with a 95 % interval resampling lineages. **Passes** if the interval is above
0.5; then (only then) a feature test is registered.

**P1b. Measured and annotated targets again, with the temporal set (the deviation of P1).** The P1 gate on the
benchmark and, as the roadmap asked, on the temporal 2015 no-prior-use positives, with the targets of ChEMBL +
GtoPdb + the **Broad Drug Repurposing Hub** (CC BY 4.0, research use; its `target` annotations). **Passes** if the
difference with random has a 95 % interval above 0 on either set; then the P1 model test (b) is run.

**R4. Orphadata prevalence (CC BY 4.0), display.** Point prevalence class per ORPHA code (product 9), ORPHA → MONDO
(SSSOM exact) → graph node; shown on rare-disease pages as information («prevalence: < 1 / 1,000,000, Orphanet»).
Release check: 50 random shown values against orpha.net, ≤ 2 wrong.

**X1. CPIC (CC0), display.** CPIC gene–drug pairs with a guideline (levels A and B); shown on candidates as «a
pharmacogenomic guideline exists for this drug (CPIC: gene)». Release check: 50 random shown pairs against
cpicpgx.org, ≤ 2 wrong. PharmGKB (CC BY-SA 4.0) when its download service answers.

Not obtainable openly (recorded, not downloaded): CredibleMeds (long QT; needs an account), the porphyria drug
databases and the G6PD / myasthenia lists (no open licence). Implementation: `sources/civic.py`, `sources/prism.py`,
`sources/repurposing_hub.py`, `sources/orphadata.py`, `sources/cpic.py`, evaluations `evaluation/o1_civic.py`,
`evaluation/o2_prism.py`, `evaluation/p1_offtargets.py` (P1b); DECISIONS P63–P67.

---

## Registered 23-10-2026, evening, before any download: T1 (trial outcomes; roadmap item 6)

AACT snapshots now need an account (not created); the same registry is read from the **ClinicalTrials.gov API v2**
(public domain): interventional studies with posted results. Kept: phase 2, 2/3 or 3; exactly one drug or biological
intervention that is not a placebo; conditions and interventions by their MeSH terms (`derivedSection`), MeSH → MONDO
(SSSOM exact) → graph disease, intervention names → graph drugs (moieties, synonyms). Outcome per trial: **met** if
any analysis of a primary outcome reports a p-value < 0.05, **not met** if primary analyses exist and all are ≥ 0.05,
otherwise unknown. The trials enter the internal graph as a layer (drug–disease, phase, outcome, NCT id).

**Test (validation of OREVID's lists, no model change).** For drug–disease pairs with a met / not-met trial, from the
website's lists as built: is the share of «met» higher among the **shown candidates** (extensions and shared-genes
candidates) than among the drugs ranked 1–25 and not shown, and than among ranks 26–100? Difference of shares with a
95 % interval resampling diseases. **Claim** («shown candidates with a posted trial met their primary endpoint more
often») only if the interval of shown − not shown is above 0. Reported: the same by phase, and approved indications
left out. Implementation: `sources/ctgov.py`, `evaluation/t1_trials.py` → `docs/T1_TRIALS.md`; DECISIONS P68.

---

## Registered 23-10-2026, evening, before any download: P2 (Perturb-seq knockdown reversal; roadmap item 14)

**Data.** Replogle et al. 2022 genome-wide Perturb-seq, K562, normalised pseudo-bulk (`K562_gwps_normalized_bulk_01.h5ad`,
figshare 20029387, CC BY 4.0): per knocked-down gene, the expression z-scores of every measured gene. Per knockdown,
its signature = the 100 most increased and the 100 most decreased genes. Disease signatures: CREEDS (as P45).
**Score** for (drug, disease): over the drug's targets that it **lowers** (M1's action types) and that were knocked
down, the mean reversal of the disease signature (P45's formula: disease-up ∩ knockdown-down + disease-down ∩
knockdown-up − same-direction overlaps); missing otherwise. **Gate**, as P1b: Recall@100 against random on the
benchmark no-prior-use positives and on the temporal 2015 no-prior-use positives, difference with a 95 % interval
resampling diseases; **passes** if the interval is above 0 on either set, and then the model test of P1 (b) is run
with this column. Implementation: `sources/perturbseq.py`, `evaluation/p2_perturbseq.py` → `docs/P2_PERTURBSEQ.md`;
DECISIONS P69. Not feasible on this PC and not downloaded: the single-cell files (66 GB), Tahoe-100M and the CELLxGENE
census (terabytes).

---

## Registered 23-10-2026, night, after the R4 check: R4b (prevalence through single-MONDO nodes only)

R4's check (50 shown values against orpha.net, read on 23-10): every value matched its Orphanet page, but **3 of 50
were on the wrong disease** (gate ≤ 2): a PrimeKG *grouped* node (several MONDO ids, named after the group) received the
prevalence of one subtype (dystrophic epidermolysis bullosa ← nails-only form; «intellectual developmental disorder,
autosomal recessive» ← a cataract–aciduria syndrome; «SEMD with joint laxity» ← Isidor type). An earlier sample drawn
before «Unknown» classes were hidden also showed a broader node with a narrower Orphanet entity (proximal spinal
muscular atrophy ← its autosomal dominant adult-onset form). Also decided: only known classes are shown.
**R4b:** an Orphanet prevalence is attached only to graph disease nodes that carry exactly **one** MONDO id, mapped
exactly (SSSOM) to that ORPHA code. New check: 50 new random shown values (other seed), the same criteria, ≤ 2 wrong.

---

## Registered 23-10-2026, night, before running: O2b (PRISM as a model column) and P1b (b)

O2's gate passed (P63) and P1b's gate passed on the temporal set (P65), so their model tests are run now, both
against the current production model (`extra_qtl`). **O2b column:** 1.5 − the drug's median PRISM AUC in the lineage
of the disease (`resources/prism_lineages.csv`, descendants), 0 when not screened or no lineage. **Rule: M2's**
(Recall@100 on all positives up with a 95 % interval above 0 on the benchmark or the areas, the other not below
−0.005, no-prior-use not down by more than 0.01). **P1b (b):** P1's registered rule (N2's), with the ChEMBL + GtoPdb +
Repurposing Hub off-targets. Implementation: `evaluation/o2_feature.py`, `evaluation/p1_offtargets.py`.

---

## Registered 23-10-2026, night, before running: B1 (NeDRex-style TrustRank baseline in OREVID-bench; roadmap item 16)

NeDRex's web service would need ≥ 1,300 remote jobs (one per test disease and task); instead its published
drug-ranking method is implemented locally, on exactly the benchmark's data, and labelled as such («TrustRank, as in
NeDRex», not NeDRex itself). **Method, fixed now:** for each task and test disease, personalised PageRank (damping
0.85, 50 iterations) on the PrimeKG protein–protein network of the task's training edges, restarting at the disease's
genes (training `disease_protein` edges); a drug's score = the sum of the TrustRank of its targets
(`drug_protein`, training edges), so drugs with no target score 0. It goes to `benchmark/predictions/trustrank/`
and the leaderboard as a baseline; no decision depends on it. Implementation: `benchmark/trustrank.py`; DECISIONS P70.

---

## Registered 24-10-2026, after G2 and before any new query: G2b (only lines stable between Open Targets releases)

G2 failed (P62): 9 of 50 lines were not supported by Open Targets 26.09, all single-study eQTL colocalisations. **G2b
shows a genetic line only if it is stable:** besides OREVID's 26.06 data, release 26.09 (the API) has a colocalisation
(h4 ≥ 0.8) of a cis QTL of the same type (pQTL or eQTL) of the gene with a GWAS mapped to the disease, with the stated
direction, and **no** such colocalisation with the opposite direction. Every one of the website's genetic lines is
queried once; the supporting 26.09 GWAS study is kept and linked on the card.

**Release check, fixed now.** Since the stability filter uses the same API, the check looks at what the filter does
not: 50 random lines that pass it, each read by a person (the assistant, answers kept) for (1) the supporting GWAS
study's reported trait is the page's disease or a form of it (not another disease or a measurement), (2) the QTL is of
the stated gene, (3) the drug's mechanism on that gene is of the stated kind (Open Targets / ChEMBL). Released only
if ≤ 2 of 50 are wrong. Implementation: `evaluation/release_checks.py` (`g2b`), `web/annotate.py`; DECISIONS P71.

---

## Registered 24-10-2026, after the G2b check and before any new query: G2c (case–control GWAS of the page's disease only)

G2b failed (P71): stable lines, but 7 of 50 supported by a GWAS that is not the page's disease (quantitative traits,
questionnaires, combined phenotypes, Open Targets mapping). **G2c** keeps a line only if, besides G2b's stability, at
least one supporting 26.09 GWAS study (i) has cases (`nCases` > 0: a case–control study, not a quantitative trait or a
questionnaire item) and (ii) has **every** one of its disease ids among the page's disease ids (no combined or
extra phenotypes). **Check:** 50 new random lines that pass (other seed), read by the assistant with the same three
criteria as G2b (trait = the page's disease or a form of it; QTL of the gene; drug mechanism), ≤ 2 wrong. If it fails,
no further version tonight: the genetic line stays off the pages. DECISIONS P72.

---

## Registered 24-10-2026, before any download: L1–L6 (more open sources, author's request)

Every source enters the internal graph as a layer (`layers.py`, P60) whatever its result; no minimum counts (author's
view, P58); a claim passes on its 95 % interval alone.

**L1. DrugCentral (CC BY-SA 4.0, dump of 02-09-2026) and MED-RT (VA / NLM, public domain).** Curated indications,
contraindications and off-label uses (DrugCentral `omop_relationship`; MED-RT `may_treat`, `CI_with`), DrugCentral
first-approval dates. Drugs by DrugBank / UNII / name, diseases by UMLS / SNOMED / MeSH → MONDO (SSSOM exact). **Exam:**
indications in either source that are neither PrimeKG indications nor confirmed Open Targets approvals; production
model vs popularity, Recall@100 per disease, paired difference resampling diseases; **claim** if the interval is above
0. Reported: how many PrimeKG indications get an approval date (for a future temporal test).

**L2. EMA orphan designations (medicines-output report).** Drug–rare-disease designations by name; **exam** as L1 on
the designations not already indications, and on those dated after 2023 (later than PrimeKG); claim on the first.

**L3. Orphadata genes (product 6) and HPO phenotypes with frequency (product 4), CC BY 4.0.** Layers; reported: rare
diseases that gain genes they had none of. No model test now (it would change the model's relations; registered
separately if wanted).

**L4. Human Protein Atlas (CC BY-SA 4.0) tissue specificity and Pharos / TCRD target development level.** Gene
attributes for target prioritisation (roadmap 11); no test now.

**L5. DepMap CRISPR gene effect (CC BY 4.0).** Per gene and cancer lineage the median effect. **Gate:** for each PRISM
lineage (`resources/prism_lineages.csv`), do the targets of drugs approved for its cancers have a stronger (more
negative) dependency in that lineage than in the others? AUROC per lineage of −(effect in the lineage − the gene's
median over lineages), mean over lineages above 0.5 with a 95 % interval → then a feature test.

**L6. OFFSIDES (nSIDES, FAERS-derived adverse events).** D1's S2 test with OFFSIDES instead of OnSIDES (confirmation
half, the same rule as D1). Implementation: `sources/drugcentral.py`, `sources/medrt.py`, `sources/orphan_ema.py`,
`sources/orphadata.py` (genes, phenotypes), `sources/hpa.py`, `sources/pharos.py`, `sources/depmap.py`,
`sources/offsides.py`; evaluations `evaluation/l_exams.py`, `evaluation/l5_depmap.py`, `evaluation/d1_worsen.py`;
DECISIONS P73–P78.

---

## Registered 24-10-2026, before running: L3b (Orphanet genes as a model column)

**Column:** for (drug, disease), log(1 + the drug's targets (`drug_protein|moa`, the split's training graph) that are
Orphanet genes of the disease (`data/ext/orphadata_genes.parquet`, single-MONDO nodes)). It uses no drug–disease edge.
**Rule: M2's** (P61), against the current production model `extra_qtl`: Recall@100 on all positives up with a 95 %
interval above 0 on the benchmark or the areas, the other not below −0.005, no-prior-use not down by more than 0.01.
Reported: the diseases that have Orphanet genes and none in PrimeKG. Implementation: `evaluation/l3b_feature.py`;
DECISIONS P79.

---

## Registered 24-10-2026, before any model run: N5–N7 (three more sources as model columns)

Coverage was counted before registering (no result looked at): current HPO annotations add 34,063 disease–phenotype
pairs to 3,061 graph diseases (248 had none); MGI mouse models add 1,104 disease–gene pairs (477 diseases); B3DB
gives blood–brain-barrier permeability for 1,996 graph drugs. Each column is tested alone against the production model
of the time (`extra_qtl_orpha`, P79) with **M2's rule** (P61). All three enter the internal graph whatever the result.

- **N5. Current HPO phenotypes** (`phenotype.hpoa` 2026-09-02; OMIM / ORPHA → MONDO, SSSOM exact): column
  log(1 + Σ over diseases *m* where the drug is indicated in the split's graph of the number of phenotypes *m* shares
  with the disease), phenotypes = PrimeKG's ∪ current HPO ∪ Orphanet product 4; the disease itself and its twins excluded.
- **N6. MGI mouse models** (`MGI_DO.rpt`; DOID → MONDO, mouse symbol → human symbol): column log(1 + the drug's targets
  among the disease's mouse-model genes).
- **N7. B3DB** (CC0): for diseases under MONDO «nervous system disorder» or «psychiatric disorder» (mondo.obo), +1 if
  the drug crosses the blood–brain barrier (BBB+), −1 if not (BBB−), 0 unknown; 0 for every other disease.

Implementation: `sources/hpo_current.py`, `sources/mgi.py`, `sources/b3db.py`, `evaluation/n_columns.py`; DECISIONS
P80–P82.

---

## Registered 24-10-2026, before running: T2 (temporal test with dated drugs)

The temporal test (I1, P42) leaves the 80 % of PrimeKG indications that Open Targets does not date in the training
graph. **T2** dates drugs with DrugCentral's first approval (any agency; `data/ext/drugcentral_approvals.parquet`, P73):
in the cutoff graph at year C, every PrimeKG drug–disease edge (indication, contraindication, off-label use) of a drug
first approved in C or later is also removed, as such a drug could not have them at C. Everything else as I1
(production model `extra_qtl_orpha`, C1 = 2012 development, C2 = 2015 confirmation, the same positives, candidates
and metrics). **Decision rule (as P42):** the claim «OREVID ranks later phase 3 / approved uses above popularity»
holds under T2 if, at C2, the paired Recall@100 difference OREVID − popularity has a 95 % interval above 0.
Reported: the residual share of PrimeKG indications with neither an Open Targets year nor a dated drug; the same
analysis for the drugs first approved at C or later («new drugs», which have no use at C) against random ranking;
the I1 numbers of the same model for comparison. Caveat stated now: the QTL and Orphanet columns and the undated
gene–disease data use 2026 knowledge. Implementation: `evaluation/temporal.py` (`dated`); `docs/TEMPORAL.md`
(section T2); DECISIONS P83.

---

## Registered 24-10-2026, before running: T3 (TxGNN on the dated temporal test)

**Question.** With the leakage-corrected temporal test (T2, P83), would TxGNN have ranked the uses that later reached
phase 3 or approval better or worse than OREVID?

**Setup, fixed now.** TxGNN's official code and the benchmark's hyperparameters (`evaluation/txgnn_wsl.py`), trained
on exactly T2's cutoff graph: every PrimeKG edge that T2 removes at year C (first trial at C or later, or a drug first
approved at C or later) is marked «test» in TxGNN's own split column, so TxGNN never trains on it (the same mechanism
as the disease-area splits). TxGNN then scores all drugs for T2's test diseases, and the scores are evaluated with
T2's candidates and positives (`data/eval/temporal_C<C>_<PRODUCTION>_dated.npz`). C1 = 2012 development, C2 = 2015
confirmation. TxGNN uses PrimeKG only, as in every comparison; OREVID's `extra_no_ot` (no Open Targets feature) is
also run on T2 so that a same-data comparison is reported.

**Decision rule.** At C2, primary positives: the claim «OREVID ranks later phase 3 / approved uses above TxGNN» if
the paired Recall@100 difference OREVID − TxGNN has a 95 % interval above 0 (and «below» if under 0). Reported:
TxGNN − popularity, the confirmed approvals, non-obvious uses and the new drugs (first approved at C or later), and
the same with `extra_no_ot`. Implementation: `evaluation/t3_txgnn.py`, `evaluation/txgnn_wsl.py` (`--label`,
`--tag`); `docs/T3_TXGNN_TEMPORAL.md`; DECISIONS P84.

---

## What to update after each item

DECISIONS (new entries P41…, one per item), STATUS, the Claude Doc of the preprint (re-export to
`docs/preprint/`; new version on Preprints.org), the website (Method page and help texts in three languages)
and, for I1/I5, OREVID-bench (leaderboard and a release with the new prediction files).
