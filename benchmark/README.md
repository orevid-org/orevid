# OREVID-bench: an honest benchmark for drug repurposing on PrimeKG

*Which existing drugs could treat this disease?* Many models answer by ranking drugs for a disease.
OREVID-bench checks how well they do it, on the same knowledge graph (PrimeKG) and the same splits as
TxGNN (Huang et al., *Nature Medicine* 2024), with three things that published evaluations usually miss:

1. **A popularity baseline.** Ranking drugs by how many diseases they already treat, ignoring the
   disease, is hard to beat: on TxGNN's zero-shot split it reaches Recall@100 0.36 (TxGNN: 0.41).
   Every model should be compared with it.
2. **Ranking over all drugs, not pairs.** TxGNN's AUPRC compares each true drug with one random
   *approved* drug. What a researcher needs is the top of the full list: Recall@100 and
   Precision@K over all 7,957 drugs.
3. **Hard diseases and novelty.** ~60 % of the «zero-shot» test diseases have a near-identical
   sibling in the ontology with the same drugs (e.g. *atopic eczema* / *dermatitis, atopic*); copying
   it is enough. OREVID-bench reports a *hard* subset without them, and novelty scenarios without
   blockbuster drugs and for drugs with no prior use. An *external* task uses approvals that are
   not in PrimeKG at all.

## Results

See [LEADERBOARD.md](LEADERBOARD.md). Summary (Recall@100, mean over diseases):

| Model | Zero-shot | Zero-shot, hard | Disease areas | External (not in PrimeKG) |
|---|---|---|---|---|
| Popularity (ignores the disease) | 0.358 | 0.203 | 0.185 | 0.130 |
| TxGNN (official code, trained on the same splits) | 0.412 | 0.258 | 0.237 | 0.121 |
| OREVID (LambdaRank over meta-paths) | **0.811** | **0.539** | **0.617** | **0.359** |

TxGNN was retrained here with its official code and hyper-parameters; on the zero-shot split it
reproduces its published numbers (Recall@100 0.412 vs 0.405 published; AUPRC 0.910 vs 0.904).

## Tasks

| Task | What is hidden | Test diseases |
|---|---|---|
| `zero_shot_{1..5}` | TxGNN's `complex_disease` split, **reproduced exactly**: 5 % of the diseases lose all their drug–disease edges; 5 % / 11.875 % of every other relation go to test / validation. | 64–71 |
| `area_{9 areas}` | All drug–disease edges of a branch of the disease ontology (cancer, mental health, cardiovascular, anaemia, adrenal, autoimmune, metabolic, diabetes, neurodegenerative), plus random edges of their 2-hop neighbourhood up to 5 % of the graph. TxGNN draws these without a seed; OREVID-bench fixes it (seed 0). | 2–180 |
| `external_ot2606` | Nothing extra: positives are drug–disease pairs **approved according to Open Targets 26.06 (ChEMBL) that are not in PrimeKG**. Train on `zero_shot_1`. Candidates exclude drugs PrimeKG already relates to the disease. Subsets: `recent` (first trial ≥ 2018), `non_obvious` (the drug is not indicated for any disease ≤ 2 ontology steps away). | 645 (3,718 pairs) |

Subsets of the split tasks: `all`; `hard` (no ontology neighbour ≤ 2 steps with ≥ 50 % of the
disease's drugs); novelty `A_approved_only` (candidates and positives restricted to drugs with a
training indication), `B_no_blockbusters` (the 10 % of drugs with most indications removed),
`C_no_prior_use` (only positives with no training indication count; the other positives are removed).

## Metrics

- **Recall@100**, **Precision@K** (K = number of the disease's drugs), **Avg Rank** (mean normalised
  rank of the positives, 1 = top) over all candidate drugs, averaged over diseases.
- **AUPRC / AUROC** with TxGNN's protocol (`reproduce/gather result-txgnn.ipynb`): positives vs the
  same number of random negatives (seed 42, with replacement) among drugs with a training indication,
  pooled over diseases.
- Ties are broken at random (seed 0).

## How to use it

Requirements: Python ≥ 3.10, numpy, pandas, scipy, scikit-learn, pyarrow.

```bash
# 1. PrimeKG as used by TxGNN (Harvard Dataverse doi:10.7910/DVN/IXA7BM)
curl -L -o kg.csv   https://dataverse.harvard.edu/api/access/datafile/7144484
curl -L -o node.csv https://dataverse.harvard.edu/api/access/datafile/7144482

# 2. (optional) rebuild the tasks; they are already in tasks/. Checks the files' MD5.
python bench.py build --kg kg.csv --nodes node.csv

# 3. Train your model on each task's training edges (see below) and write predictions
# 4. Evaluate
python bench.py evaluate --pred predictions/my_model
python bench.py baseline          # popularity baseline -> predictions/popularity
```

**Training edges.** `tasks/<task>.split.npz` holds the row numbers (in `load_kg(kg.csv)` order) of
the validation and test edges; everything else is training. In Python:

```python
from bench import load_kg, train_edges
kg = load_kg("kg.csv")
mask = train_edges(kg, "tasks", "zero_shot_1")   # boolean, True = you may use this edge
```

Do not use anything that encodes drug–disease therapeutic relations for the test diseases (e.g.
clinical-trial databases, drug labels, literature saying «drug treats disease»).

**Predictions.** One file per task, `predictions/<model>/<task>.parquet` (or `.csv`) with columns
`disease` and `drug` (PrimeKG `node_index`) and `score` (higher = more likely to treat). Score every
drug in `tasks/<task>.json` → `drugs` for every test disease; missing pairs get the lowest score.

The prediction files of the four models in the leaderboard (~180 MB, not kept in git) are attached to
the GitHub release as `orevid-bench-predictions-v1.0.zip`: unzip it here to re-evaluate them with
`python bench.py evaluate --pred predictions/<model>`.

## Caveats

- `area_*`: TxGNN's own disease-area splits cannot be reproduced (no seed); OREVID-bench's are an
  approximation with a fixed seed. One seed only; diabetes (2 diseases) and adrenal (6) are noisy.
- `external_ot2606` is out-of-sample but not prospective: most of those approvals are older than
  PrimeKG (median first trial 2008); they are simply missing from it.
- Labels come from PrimeKG/DrugBank and Open Targets; they reward rediscovering known treatments.
- TxGNN outputs the same score for every drug in test diseases with no edge left in the training graph
  (8 of 335 zero-shot diseases, 3 of 345 in the areas); its ranking there is random.
- OREVID normally uses the validation labels of each split as training examples; TxGNN uses them only
  for early stopping. Both are allowed by the protocol, but it is an asymmetry.
- Both differences are **measured** in [`FAIR_COMPARISON.md`](FAIR_COMPARISON.md) (`fair_comparison.py`):
  with OREVID trained on the training labels only and the degenerate diseases removed, zero-shot
  Recall@100 is 0.809 vs TxGNN 0.421 (hard subset 0.528 vs 0.275; areas 0.627 vs 0.238). The validation
  labels add 0.01–0.02 to OREVID; removing the degenerate diseases helps both models similarly.
- Uncertainty: [`CLUSTER_CI.md`](CLUSTER_CI.md) (`cluster_ci.py`) gives Recall@100 intervals from a bootstrap
  that resamples whole diseases, paired differences per disease (Wilcoxon, Holm-adjusted) and the effect of
  the random tie-breaking (at most 0.014 over 20 seeds; none for OREVID).

## Data and licences

- PrimeKG (Chandak, Huang & Zitnik, *Sci Data* 2023), CC0 1.0. Open Targets Platform 26.06, CC0 1.0.
- Disease-area lists (`disease_areas/`) from the TxGNN repository (MIT licence).
- OREVID-bench code: MIT licence; task files and predictions: CC BY-NC-SA 4.0 (see LICENSE and LICENSE-DATA.md in the main repository).
