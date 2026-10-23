# Roadmap: new sources and new directions, ranked by potential (written 23-10-2026)

Written on the work PC, from a review of OREVID's known weak points against the public resources available in
2026. **Nothing here is registered.** Each item gets its own section in `docs/IMPROVEMENT_PLAN.md`, with fixed
criteria and decision rules, committed before any run (as R1–R2 and N1–N4); results go to `docs/DECISIONS.md`.

**Potential** = how much the item can improve what OREVID gives researchers (better candidates, fewer wrong ones,
safer warnings, more credibility), weighed by how likely it is to work and how much it costs. Items are ranked by
potential; the effort column says what each costs.

## The weak points the items answer

| # | Weak point | Where it shows |
|---|---|---|
| W1 | Drugs with no prior indication are not found | Precision@K 0–2 % for every model tested (STATUS, open problems) |
| W2 | No direction of effect: will the drug improve or worsen the disease? | genetic direction for 1 % of the pairs; the direction model copied neighbours' labels (`docs/DIRECTION.md`) |
| W3 | Noisy gene–disease associations | shared-genes outliers: insulin for MS, paclitaxel for Castleman (STATUS item 9; R1b did not remove them) |
| W4 | Infections: the graph is human, the pathogen is not in it | R2 failed because PrimeKG has no pathogen classes (P49) |
| W5 | Safety knowledge missing from the graph | antiepileptics that aggravate seizures in Angelman syndrome, shown without a warning |
| W6 | Little evidence from real patients | evaluation labels are approvals; trial outcomes and real-world use are unused |
| W7 | Outdated or wrongly extracted drug status | N4 (`docs/DRUG_STATUS.md`); second version pending |

## Shared rules for every item

- **No leakage.** Any source that encodes drug–disease therapeutic relations (trial results, case reports,
  literature «treats», real-world use) is used only as displayed evidence or as evaluation labels, never as a
  feature or a training label.
- **Licences** checked before anything is published (DrugBank, OncoKB, CredibleMeds and UK Biobank individual data
  have restrictive terms; GenCC, CIViC, Open Targets are open).
- **Research use only.** Nothing that invites patients (no supplements, no symptom-based diagnosis).
- Every claim compared with the popularity baseline and with disease-level intervals, as now.

## Ranking

| Rank | Item | Answers | Potential | Effort | Depends on |
|---|---|---|---|---|---|
| 1 | Curated gene–disease validity (GenCC, ClinGen, PanelApp) | W3 | high | low | — |
| 2 | «Drugs that may worsen this disease» + drugs-to-avoid lists + OnSIDES | W5, W2 | high | medium | — |
| 3 | OREVID-bench as an open leaderboard on Hugging Face | credibility | high | low | — |
| 4 | Drug-target Mendelian randomisation with pQTL (UKB-PPP, deCODE) + FinnGen / UKB GWAS | W2, targets | high | high | — |
| 5 | Full polypharmacology: ChEMBL bioactivities, IUPHAR / BPS, Broad Repurposing Hub | W1 | medium–high | medium | — |
| 6 | Trial outcomes (AACT results, CTIS / EudraCT) | W6 | medium–high | medium | — |
| 7 | N4 second version (drug status against current labels) | W7 | medium–high | medium | in progress |
| 8 | Disease → infectious agent (MONDO / DO, NCBI Taxonomy) + drug spectrum | W4 | medium | medium | — |
| 9 | Oncology: PRISM / DepMap repurposing screen + CIViC | W1 (cancer) | medium–high for cancer pages | medium | — |
| 10 | Gene burden of loss-of-function (Genebass) | W2 | medium | low–medium | 4 helps |
| 11 | Target prioritisation as a product of its own | targets | medium | medium | 4, 10 |
| 12 | Trial-feasibility support (Orphadata prevalence / natural history + I3) | usefulness | medium | low–medium | — |
| 13 | Stratified use: CPIC / PharmGKB, CIViC in the biomarker table | W5, precision | medium | medium | 9 for CIViC |
| 14 | Perturbation and single-cell atlases (Replogle Perturb-seq, Tahoe-100M, CELLxGENE) | W1 | uncertain | high | — |
| 15 | Real-world evidence networks (OHDSI / EHDEN, Nordic registries, UK Biobank prescriptions) | W6 | very high, long term | very high (partners, applications) | — |
| 16 | Already open: N3 Bioteque, NeDRex baseline, pathway overlap with its own null | W1, credibility | low–medium | low–medium | IRB websites (N3) |
| 17 | Binding-site similarity on AlphaFold / PDBe-KB structures | W1 | low–medium | high | 5 first |

---

## 1. Curated gene–disease validity (GenCC, ClinGen, PanelApp)

**Registered 23-10-2026 as G1** (`docs/IMPROVEMENT_PLAN.md`, «Registered 23-10-2026»): GenCC only (it already gathers
ClinGen and PanelApp), variants G1a–G1c and their decision rules. To run on the home PC.

**What.** GenCC aggregates gene–disease validity classifications from ClinGen, Orphanet, PanelApp, Genomics England
and others, graded definitive / strong / moderate / limited / disputed / refuted. ClinGen gives the curated
evidence; PanelApp (Genomics England) gives expert-reviewed gene panels per condition.

**Why first.** The shared-genes line and many meta-paths rest on gene–disease associations, part of which are weak
(text-mined, GWAS proximity, comorbidity). A validity grade lets a rule keep only well-supported genes, or weight
them, for the diseases GenCC covers (mostly Mendelian and rare). Cheap: a few files, a mapping to MONDO (GenCC uses
MONDO already).

**What to test.** (a) The shared-genes line restricted to genes with ≥ moderate validity where the disease has
such genes; (b) the effect on the insulin–MS and paclitaxel–Castleman cases is reported but does not decide;
(c) adoption rule as I6 / R1: hit rates on benchmark, areas and external exam not lower beyond the interval, and the
removed candidates hit less than the kept ones.

**Limits.** Coverage is good for rare Mendelian diseases, poor for common complex ones (MS, diabetes): for those,
item 4 is the better route. Licence: GenCC data are openly downloadable (check the exact terms at download).

## 2. «Drugs that may worsen this disease»

**What.** A new section per disease: approved drugs that may worsen it, each with its reason and source. Inputs:
contraindications in neighbouring diseases (the rule already validated at 97.4 %, `docs/DIRECTION.md`), the
checkpoint-inhibitor rule, **curated drugs-to-avoid lists** (acute porphyrias: the NAPOS drug database; Brugada
syndrome: brugadadrugs.org; long QT: CredibleMeds; G6PD deficiency; myasthenia gravis; malignant hyperthermia),
**OnSIDES** (adverse effects extracted from current FDA labels, replacing SIDER, which dates from 2015) and, later,
the direction from item 4.

**Why high.** It is the other side of repurposing and arguably more useful to clinicians: a warning, not a
prescription, so less risk of misuse. Most of the machinery exists and it can be validated the same way (hide
known contraindications and recover them).

**What to test.** Recovery of held-out PrimeKG contraindications and of drugs on the curated lists, by disease,
against the popularity of «contraindicated anywhere»; a manual check against sources before release (with the
lesson of N4: against the source pages, by a person or an agent that opens them).

**Limits.** Each curated list has its own terms (CredibleMeds needs registration); some lists cover few diseases.
Confounding by indication in adverse-effect data (`docs/DIRECTION.md`) applies: OnSIDES is a reason to show, not a
score.

## 3. OREVID-bench as an open leaderboard

**What.** A Hugging Face space (or a results folder plus a submission template) where anyone submits prediction
files for the 15 tasks and gets the same evaluation (`bench.py`), the popularity baseline beside it, the external
exam and disease-level intervals.

**Why high.** Nobody offers an honest, reproducible repurposing benchmark with a popularity baseline and an
out-of-graph exam; it makes OREVID the reference evaluator of the field and invites other groups (TxGNN, Every
Cure, Bioteque) to compare. Low effort: the dataset, `bench.py` and six models are already published.

**What to decide.** Automatic evaluation (a space running `bench.py`) or reviewed pull requests; rules against
using test labels (the submitter declares the training data; the external exam stays hidden-label if possible).

**State 23-10-2026 (home PC): submission route ready, nothing published.** Chosen as the first step: pull requests
evaluated automatically (no hosting, every result auditable in the repository). `benchmark/SUBMITTING.md` fixes the
rules (training edges only, no outside drug–disease relations, settings on validation or `zero_shot_1`, open code and
DOI for the predictions; the external exam is honour-based because its labels are public), `bench.py validate`
checks a prediction folder, `submission_check.py` downloads the zip, checks its MD5, extracts only task files and
evaluates with the `main` code, and `.github/workflows/bench-submission.yml` runs it on any pull request that adds one
`benchmark/submissions/<model>/submission.json` (nothing from the pull request is executed). `leaderboard.py` lists
accepted submissions. Tested locally with the popularity predictions (same numbers as the leaderboard). **Author's
decisions:** announce it (Hugging Face dataset card, README of the dataset, preprint), and whether to add a Hugging
Face Space later; a hidden-label track would be a prospective round like the 12-10 snapshot.

## 4. Drug-target Mendelian randomisation with protein QTLs

**What.** For each drug target, genetic variants near its gene that change the protein's level (cis-pQTL: UK
Biobank Pharma Proteomics Project, ~2,900 proteins; deCODE, ~4,900 aptamers; eQTLGen / GTEx for expression) used as
instruments, with disease GWAS summary statistics (FinnGen, UK Biobank, consortia via the GWAS Catalog): does
lowering the protein lower or raise the risk of the disease? With colocalisation to avoid confounding by linkage.

**Why high.** It gives **direction** (W2), which OREVID lacks almost everywhere, and it is the method the industry
uses to prioritise targets (Open Targets' own genetics points the same way). It also gives on-target safety: a
protein whose lowering raises the risk of another disease is a warning.

**What to test.** (a) Do drugs whose target MR agrees with the indication direction rank above the rest among
held-out approvals and in the temporal test? (b) As a risk signal: do contraindications show the opposite
direction more often than chance? Adoption rules as I4.

**Limits.** Heavy (many GWAS, careful MR: instrument strength, pleiotropy, colocalisation); coverage limited to
proteins with cis-pQTL and diseases with good GWAS (common diseases). Summary statistics are public; individual
UK Biobank data are not needed.

## 5. Full polypharmacology (ChEMBL bioactivities, IUPHAR / BPS, Broad Repurposing Hub)

**What.** Today OREVID uses each drug's mechanism-of-action targets (DrugBank, ChEMBL). Add every measured
off-target with potency (ChEMBL activities with pChEMBL ≥ 6 or 7), the IUPHAR / BPS Guide to Pharmacology
(curated, quantitative) and the Broad Drug Repurposing Hub annotations.

**Why.** Many repurposings come from secondary targets (W1): a drug with no indication near the disease may still
hit one of its genes.

**What to test.** Meta-paths through measured off-targets as added features, with the I4 rule on the no-prior-use
subset and the temporal set; report coverage first.

**Limits.** Assay heterogeneity; potency thresholds fixed before running. ChEMBL is CC BY-SA 3.0: check what derived
fields may be published.

**State 23-10-2026:** registered and run as P1 (P59, `docs/P1_OFFTARGETS.md`): ChEMBL + GtoPdb give 14,722 drug–gene
pairs; standalone on no-prior-use positives 0.042 vs random 0.013 with an interval including 0, so the gate failed and
the model test was not run. Repurposing Hub not used (licence).

## 6. Trial outcomes (AACT results, CTIS / EudraCT)

**What.** OREVID counts trials per pair; it does not know whether they worked. ClinicalTrials.gov results (AACT:
primary outcome, statistical analyses, «terminated for lack of efficacy») and EU trials (CTIS, EudraCT results).

**Why.** «Tried in phase 2 and failed» and «met its primary endpoint» are very different for a researcher, and the
second is the closest open evidence to clinical efficacy (W6). It can also add evaluation labels (positive phase 2
results as a softer outcome than approval).

**Use.** Displayed evidence and evaluation labels only, never features (leakage). A classifier of «positive /
negative / unclear» trial outcome, checked by hand against the trial pages before release.

## 7. N4 second version

Already described: `docs/DRUG_STATUS.md` and `docs/IMPROVEMENT_PLAN.md` («Added 23-10-2026»). Release in parts,
checked against the source pages. Listed here for completeness: it fixes errors visitors see today.

## 8. Infections: disease → infectious agent, drug → spectrum

**What.** MONDO and the Disease Ontology link many infectious diseases to their agent (NCBI Taxonomy). Drugs get a
spectrum from ATC (J01 antibacterials, J02 antifungals, J04 antimycobacterials, J05 antivirals, P antiparasitics)
and from ChEMBL target organisms.

**Why.** It makes the failed R2 possible (P49): extensions and candidates of the wrong pathogen class can be flagged
or removed, and the «targets» section can be hidden or relabelled for infections, where human targets make no
sense.

**What to test.** R2 again, registered anew, with the approved pairs as the check that nothing right is removed.

**State 23-10-2026:** registered and run as R3 (P57, `docs/R3_INFECTIONS.md`): the MONDO classes and ATC codes work, but
the rule removes one approved pair (sulfanilamide, vulvovaginal candidiasis) and too few candidates to show a gain.
Rejected. ChEMBL target organisms were not used.

## 9. Oncology: PRISM / DepMap repurposing screen and CIViC

**What.** The PRISM repurposing screen measured the effect of ~4,500 drugs (most non-oncology) on ~580 cancer cell
lines; DepMap gives the lines' genetics. CIViC (CC0) curates clinical evidence linking tumour variants to drugs.

**Why.** Direct experimental evidence for W1 in cancer, the largest group of diseases in the graph; CIViC brings
the biomarker table up to date for cancer (subtypes the graph lacks: EGFR, ALK, BRAF, KRAS…, STATUS 14-10).

**What to test.** Do drugs that kill the cell lines of a cancer type rank its held-out approvals and later phase 3
pairs above the rest? Feature or displayed evidence decided by the I4 rule.

**Limits.** Cell-line killing is far from clinical benefit; mapping cell lines to graph diseases needs care.

## 10. Gene burden of loss-of-function (Genebass)

UK Biobank exome burden tests per gene and phenotype: if losing a gene protects from a disease, an inhibitor of
its product might help. Complements item 4 for genes without pQTL; same tests and rules.

## 11. Target prioritisation as a product of its own

The website already has a «targets» section (`docs/TARGETS.md`). With items 4 and 10 it can say, per disease,
which protein, in which direction and with which modality (Open Targets tractability), with the evidence behind
it. Interest from industry and rare-disease foundations. Validated as now (hidden approved-drug targets) plus the
temporal test.

## 12. Trial-feasibility support

Per candidate pair: prevalence and natural history (Orphadata, CC BY 4.0), earlier trials and their outcomes
(item 6), generic and cheap or not (I3 already). What a foundation needs before funding a repurposing trial.
Display only.

## 13. Stratified use: pharmacogenomics and tumour variants

CPIC guidelines and PharmGKB (CC BY-SA) for patient groups at risk; CIViC (item 9) for variant-defined subtypes in
the biomarker table. Display only, checked by hand.

## 14. Perturbation and single-cell atlases

Genome-wide Perturb-seq (Replogle et al., CRISPRi knock-downs), Tahoe-100M (single-cell responses to ~1,100 drugs
in 50 cancer lines), CELLxGENE Census (disease vs normal by cell type). The next generation of the signature
reversal that LINCS failed (I4, negative): cell-type resolution might change the result, but it is uncertain and
heavy. Only after items 1–9.

## 15. Real-world evidence networks

OHDSI / EHDEN network studies, Nordic registries, UK Biobank prescriptions: are people who take drug X less likely
to develop disease Y? The validation clinicians value most (W6), but it needs partners, applications and months.
Long term; the outreach contacts (CURE ID, Every Cure, REMEDi4ALL) are the way in.

## 16. Already open

N3 Bioteque (terms agreed with Patrick Aloy; waiting for the IRB websites), NeDRex as a benchmark baseline
(comparable, ≥ 1,300 remote jobs), pathway overlap with its own null (P46).

## 17. Binding-site similarity

Pocket similarity between drug targets and disease proteins on AlphaFold / PDBe-KB structures, for off-target
hypotheses. High effort, overlaps with item 5, which should come first.

## Not proposed

Supplements and natural products (they draw patients), symptom-based diagnosis (another product, another risk),
anything that turns a candidate into advice for an individual patient.
