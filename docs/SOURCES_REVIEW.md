# Sources reviewed (23–24 October 2026)

Every open source considered for OREVID, what was done with it and why. Sources that were obtained are layers of the
internal graph (`docs/GRAPH_LAYERS.md`) whatever their result; decisions are in `docs/DECISIONS.md`.

## In the model (ranking)

| Source | How | Decision |
|---|---|---|
| PrimeKG, Open Targets 26.06, STRING 12.5 | core graph and metapaths | production since 12-10 |
| Open Targets QTL colocalisations (pQTL / eQTL) | column «qtl direction» | M2 not adopted by the rule; adopted by the author (P61) |
| Orphadata genes (product 6) | column «orphanet genes» | L3b borderline; adopted by the author (P79) |

## Tested as model columns, not adopted

ChEMBL 37 + GtoPdb + Broad Repurposing Hub off-targets (P59, P65), network proximity (P51), LINCS / CREEDS reversal (P45),
Perturb-seq knockdown reversal (P69), PRISM lineage sensitivity (P63), RotatE embeddings (P47). Also current HPO phenotypes (N5,
P80; it harmed new drugs), MGI mouse models (N6, P81) and B3DB blood–brain barrier (N7, P82).

## Used as exams (independent validation of OREVID)

| Source | Result |
|---|---|
| CIViC oncology evidence (P64) | 0.55 vs popularity 0.21 |
| DrugCentral + MED-RT indications (P73) | 0.39 vs 0.13 on 863 diseases |
| EMA orphan designations (P74) | +0.15; from 2023 on +0.08 |
| ClinicalTrials.gov posted results (P68) | shown candidates do not meet endpoints more often (claim not supported) |

## Shown on the website as information

Orphanet prevalence (P66), CPIC pharmacogenomic guidelines (P67), genetic direction lines checked against Open Targets
(P72), current US / EU status (P55), generics (P44).

## In the internal graph only

GenCC (P54), OnSIDES (P56), OFFSIDES (P78), MONDO infection classes and ATC (P57), DrugMechDB (P50), DepMap CRISPR
(P77), Human Protein Atlas tissue specificity and TCRD target levels (P76), PharmGKB, Orphadata phenotypes (P75),
label checks (P55), CTD, PubTator3, Perturb-seq, PRISM, CIViC, Repurposing Hub, trials.

## Not obtainable or not usable

| Source | Reason |
|---|---|
| Bioteque | server does not answer from the home PC (IRB); terms agreed with P. Aloy |
| AACT snapshots, CredibleMeds | need an account (the ClinicalTrials.gov API was used instead of AACT) |
| Genebass | requester-pays cloud bucket |
| DepMap portal (current release) | bot check; Public 24Q2 from figshare used |
| FDA orphan designations (OOPD) | the site moved; not found |
| WHO Essential Medicines List | no machine-readable download |
| SemMedDB | UMLS licence |
| DisGeNET (current), COSMIC, OncoKB, KEGG, TTD, GBD | commercial or restrictive licences |
| Porphyria, myasthenia drug lists | no open licence (G6PD covered by CPIC) |
| Tahoe-100M, CELLxGENE census, single-cell Perturb-seq | tens of GB to terabytes |
| NIH RePORTER | drug and disease only in free text; too noisy without text mining |
| Hetionet, DRKG, Monarch KG, BindingDB, STITCH, WikiPathways | mostly repeat sources already in the graph |
