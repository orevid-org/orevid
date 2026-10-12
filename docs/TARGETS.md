# Targets: which gene or protein to hit, in which direction and with what

Command: `py -3 explore.py targets "multiple sclerosis" [--top 15] [--with-evidence]`.
Code: `evaluation/targets.py` (model and evaluation), `sources/targets.py` (annotations for each gene).
Model: `py -3 -m evaluation.targets train` (saves `data/prep/targets_model.pkl`).

## What it answers

For a disease, a list of candidate target genes/proteins, separating:

- **already exploited**: targets of drugs that are already approved for the disease;
- **new targets**, and for each one:
  - the **evidence** linking it to the disease (genetics, somatic mutations, animal models,
    expression, PrimeKG) or whether it only gets there through the network or through similar diseases;
  - the **direction** required by human genetics (activate or inhibit) when known, with a warning if
    existing drugs go against it;
  - **what it can be hit with** (Open Targets tractability): small molecule, antibody, PROTAC,
    oligonucleotide, with the level (clinical / evidence / predicted);
  - **known risks** of the target and whether the gene tolerates loss of function poorly (LOEUF < 0.35);
  - the **drugs already acting on it** (ChEMBL) and in which direction: a starting point for repurposing.

With `--with-evidence` only genes with direct evidence in the disease are shown: this is what matters
when looking for targets that modify the disease, not just the symptoms.

## Sources (all were already in the project)

- Open Targets 26.06, `target` table (CC0): tractability, safety, genetic constraint (LOEUF),
  chemical probes, class and location (`sources/targets.py` -> `data/ext/targets.parquet`, 31,261 genes).
- Open Targets gene-disease associations by evidence type (without the literature, see below),
  PrimeKG, genetic direction and GWAS+QTL, STRING, ChEMBL mechanisms.

## Evaluation

For each disease, the positives are the targets (ChEMBL mechanism) of the drugs already approved for
it (PrimeKG + Open Targets): 1,199 diseases, ~10 targets per disease. Cross-validation on
two halves of the diseases. Recall@20: of the real targets, how many appear among the first 20 genes.

| Method | Recall@20 | Recall@100 | Precision@K |
|---|---|---|---|
| Genetics only (Open Targets) | 0.06 | 0.17 | 0.02 |
| Total evidence (sum) | 0.09 | 0.22 | 0.04 |
| Target popularity (how many drugs act on it) | 0.23 | 0.55 | 0.12 |
| «Biology only» model (⚠ included drug information; strict figure: 0.24, see the end) | 0.33 | 0.67 | 0.19 |
| + target popularity | 0.35 | 0.70 | 0.21 |
| + targets of drugs in diseases with genes in common | 0.44 | 0.77 | 0.30 |
| **Full model** (+ targets in neighbouring diseases in the ontology) | **0.59** | **0.82** | **0.41** |

90 % of the real targets are among the candidates (genes with evidence, their STRING neighbours and
targets of drugs for similar or neighbouring diseases).

Reading:

- **Genetics alone recovers few targets of approved drugs** (6 % among the first 20): many
  approved drugs act on targets without direct genetic evidence (they often treat symptoms).
- The biological model recovers 5-6 times more (strict version, at the end: Recall@20 0.24). The jump to the full model comes from copying targets
  from neighbouring diseases (forms of the same disease): useful, but not very new.

## Leaks found and avoided (important)

1. **The Open Targets literature** (Europe PMC) is left out: it often mentions the target because there
   already is a drug.
2. **«Clinical» tractability** (there is already an approved drug on the target): it may come from the drug for
   the same disease. It is limited to «evidence / predicted».
3. **«Target popularity» minus the drugs of the same disease: NO.** It looked correct,
   but it creates a leak: the model recognises each gene by its fixed features (STRING degree,
   pleiotropy…) and learns «this gene comes out lower than usual ⇒ it is positive». Recall@20 rose
   from 0.32 to 0.62, falsely. How many drugs act on a gene is a fact about the gene that does not depend on the
   disease, and it is used in full.

## Example: multiple sclerosis (`--with-evidence`)

Already exploited: potassium channels (targets of dalfampridine). New targets with direct evidence,
among the first: SCN3A, **VDR** (vitamin D), BCHE, **PDE4A** (ibudilast, a PDE4 inhibitor,
has been trialled in progressive MS; the GWAS direction says «activate», and the output warns about it), **IL12A/IL12B**
(ustekinumab was tried and did not work), **CNR1** (cannabinoids; one is approved for
MS spasticity in Europe), TNFSF14.

## Limitations and possible improvements

- The exam rewards rediscovering targets already exploited in similar diseases.
- Sources that would help are missing: cellular essentiality (DepMap, CC BY 4.0), expression by tissue
  (GTEx/HPA), protein structure and pockets (AlphaFold/PDB) and compound bioactivity
  (ChEMBL activities: molecules that bind even if they are not drugs).
- Genetic direction is incomplete (many targets «unknown»).

## Correction after the independent review (11-10)

The «biology only» model in the evaluation table **was not free of drug information**: 38 %
of the positives entered the candidates only because they were targets of drugs for similar or
neighbouring diseases, and there were variables that give away «it is already a target of some drug» (chemical probes, safety
risks, ligands, tractability). Strict version (`py -3 -m evaluation.targets --strict`):
candidates only through gene-disease evidence or STRING neighbours, and without those variables.

| Method (strict) | Recall@20 | Recall@100 | Precision@K |
|---|---|---|---|
| Genetics only (Open Targets) | 0.043 | 0.079 | 0.021 |
| Target popularity | 0.147 | 0.318 | 0.068 |
| **Model: biology only** | **0.241** | **0.408** | **0.130** |
| Full model | 0.372 | 0.481 | 0.276 |

Coverage: only 52 % of the real targets have direct biological evidence or are STRING neighbours.
**The honest figure for the biological model is Recall@20 0.24** (not 0.33), 5-6 times genetics alone.
