# Comparison with TxGNN (09-10/10/2026)

> **Historical document.** It records how TxGNN was reproduced and the first comparison (gradient boosting,
> `--paths extra`). The current model (LambdaRank, `--paths extra_no_ctd --rounds 5 --rank`) and its figures
> are in [`benchmark/LEADERBOARD.md`](../benchmark/LEADERBOARD.md) (zero-shot Recall@100 0.811 vs TxGNN 0.412).

Summary: with the same graph, the same splits and the same metric code, our model
(paths through the graph + gradient boosting) **puts the good drugs at the top of the list far
better than TxGNN on all 14 splits** (Recall@100 0.79 against 0.41 on zero-shot; double or more
on almost all areas), also without using popularity. On AUPRC/AUROC with sampled
negatives, we win on zero-shot and TxGNN wins narrowly on most areas.

Split-by-split progress table: `docs/TXGNN_PROGRESS.md`.

## How to reproduce it

```
py -3 -m evaluation.txgnn --paths extra --rounds 5            # our model, zero-shot seeds 1-5
py -3 -m evaluation.txgnn --split areas --seeds 1 --paths extra --rounds 5
powershell -ExecutionPolicy Bypass -File evaluation\txgnn_queue.ps1   # TxGNN on WSL (~23 h of GPU)
py -3 -m evaluation.txgnn_progress                            # table with both
```

- Data: `data/raw/txgnn/kg.csv` (Dataverse, file 7144484, 980 MB) and `node.csv` (7144482). It is
  PrimeKG as used by TxGNN; the first 129,375 nodes match our
  `data/prep/nodes.pkl` exactly.
- TxGNN's published figures: `resources/txgnn/published_results.csv` (TxGNN repository, MIT).
  Disease lists per area: `resources/txgnn/*.csv`.
- **TxGNN run here** (`evaluation/txgnn_wsl.py`, `train_txgnn.sh`, `txgnn_queue.ps1`): its
  official code (github.com/mims-harvard/TxGNN) and the hyperparameters of `reproduce/train.py`
  (pre-training 1 epoch, fine-tuning 500), on WSL Ubuntu with Python 3.10, PyTorch 1.13.1+cu116,
  DGL 0.9.1 (cu116, with the `nvidia-*-cu11` libraries on the `LD_LIBRARY_PATH`), pandas 1.5.3, on
  an RTX 3060. ~95 min per split. The per-disease predictions are in
  `data/eval/txgnn_predictions/` and are evaluated with exactly our own code.
  - Pitfalls: WSL has 7 GB; TxGNN's evaluation re-reads the whole `kg.csv` and crashes,
    so after training a minimal `kg.csv` (names only) is put in place and the model is saved before
    evaluating. Do not run anything else on WSL at the same time. `evaluation/keep_awake.ps1` prevents Windows
    from going to sleep while training.

## Protocol

- **Zero-shot (`complex_disease`)**: 5 % of the diseases with drugs lose all their drug-disease
  edges; from the remaining relations 5 % and 11.875 % are removed. Reproduced **exactly**:
  TxGNN's code generates the same split as ours (71 diseases and 444 indications for
  seed 1; same number of training, validation and test edges).
- **Disease areas** (9 DOID branches): all diseases in the area lose their drug-disease
  edges and random edges from the neighbourhood are removed up to 5 % of the graph. TxGNN sets no
  seed; we use our (approximate) split **for both models**. One seed.
- **Metrics** (TxGNN's `gather result-txgnn.ipynb` notebook), «indication» relation:
  - AUPRC and AUROC: for each disease, the positives and the same number of random negatives
    among the drugs with some indication in training (seed 42).
  - Recall@100, Precision@K (K = number of positives), mean rank (1 = at the top): over
    **all** 7,957 drugs. These are the ones that correspond to real use.
- **Hard diseases**: test diseases with no disease ≤ 2 steps away in the ontology that already has
  ≥ 50 % of their drugs. ~40 % of the zero-shot set; the rest have a near-identical sibling
  (e.g. «atopic eczema» ↔ «dermatitis, atopic») and copying its drugs is enough.

## Our model

Gradient boosting over ~24 metapaths (DWPC) + popularity + drug degrees, trained on
diseases that are not test diseases (`evaluation/txgnn.py`):

- `base`: 9 PrimeKG paths (target, interaction, «indicated in a disease with genes in common», ontology…).
- `new`: + adverse effects = phenotypes, target → pathway/process/phenotype, 2-step ontology and
  «similar drugs» (by targets, effects or interactions) to those that score.
- `extra`: + OREVID biological sources (Open Targets genes, ChEMBL, CTD drug-gene, STRING).
  **No drug-disease relation at all** (no trials, no CTD therapeutic, no PubTator, no FAERS).
- `--rounds 5`: **masked training**: in each round the drug-disease edges of whole
  neighbourhoods (siblings in the ontology and diseases with almost identical genes) are hidden and those diseases
  are added as examples, so that the model learns under the same conditions as the test.

## Final results

> **Update (11-10)**: after the independent review, the current model does not use CTD and has
> the ontology metapath fixed. Current figures in `benchmark/LEADERBOARD.md` (Recall@100: zero-shot
> 0.811, hard 0.539, areas 0.617, external 0.359). The website additionally applies the specificity ranking
> (0.80 / 0.53 / 0.60 / 0.35). The tables in this section are from the previous version (with CTD).

**Zero-shot, mean of 5 seeds** (± standard deviation across seeds):

| | AUPRC | AUROC | Recall@100 | Precision@K | Mean rank |
|---|---|---|---|---|---|
| TxGNN published | 0.904 | 0.910 | 0.405 | 0.072 | 0.782 |
| **TxGNN run here** | 0.910 ± 0.04 | 0.909 ± 0.04 | 0.412 ± 0.04 | 0.079 ± 0.02 | 0.783 ± 0.04 |
| Popularity | 0.807 | 0.756 | 0.358 | 0.097 | 0.905 |
| Ours, without popularity | 0.928 | 0.896 | 0.771 | 0.501 | 0.902 |
| **Ours** (extra + masked) | **0.955** ± 0.01 | **0.938** ± 0.02 | **0.789** ± 0.04 | **0.486** ± 0.05 | **0.970** ± 0.003 |
| *Hard:* popularity | 0.692 | 0.615 | 0.203 | 0.039 | 0.835 |
| *Hard:* TxGNN | 0.837 | 0.830 | 0.258 | 0.046 | 0.693 |
| *Hard:* **ours** | **0.866** | **0.831** | **0.509** | **0.183** | **0.928** |

TxGNN run here reproduces almost exactly what was published: the comparison is fair.

**Disease areas** (one seed). Recall@100 and AUPRC, TxGNN → ours:

| Area (test diseases) | Recall@100 | Recall@100 hard | AUPRC | Mean rank |
|---|---|---|---|---|
| Cancer (180) | 0.35 → **0.73** | 0.23 → **0.57** | 0.95 → 0.95 | 0.81 → **0.95** |
| Mental health (39) | 0.21 → **0.43** | 0.20 → **0.35** | **0.92** → 0.86 | 0.68 → **0.89** |
| Cardiovascular (53) | 0.11 → **0.34** | 0.09 → **0.29** | **0.77** → 0.74 | 0.49 → **0.89** |
| Metabolic (28) | 0.08 → **0.53** | 0.01 → **0.32** | **0.93** → 0.90 | 0.64 → **0.91** |
| Autoimmune (13) | 0.33 → **0.72** | 0.29 → **0.47** | **0.91** → 0.90 | 0.74 → **0.93** |
| Anaemia (12) | 0.08 → **0.67** | 0.00 → **0.44** | 0.65 → **0.98** | 0.53 → **0.97** |
| Neurodegenerative (12) | 0.48 → **0.81** | 0.19 → **0.62** | **0.98** → 0.91 | 0.90 → **0.96** |
| Adrenal (6) | 0.47 → **0.91** | 0.23 → **0.89** | **0.98** → 0.97 | 0.58 → **0.99** |
| Diabetes (2) | 0.03 → **0.23** | 0.03 → **0.23** | 0.69 → **0.70** | 0.69 → **0.86** |

Without popularity («ours, graph only»), Recall@100 stays above TxGNN on all 14
splits (e.g. mental health 0.42, metabolic 0.67, cardiovascular 0.35).

## Reading

1. **For producing a list of candidates** (Recall@100, Precision@K, mean rank), our model
   beats TxGNN on all 14 splits, often by double or more, also on the hard diseases and
   also without popularity. The advantage comes from the biological paths, not from popularity.
2. **TxGNN has a good signal for separating a good drug from any approved drug** (AUPRC/AUROC):
   it wins narrowly on 6 of 9 areas, and clearly on mental health and neurodegenerative. But in the
   full list it puts many compounds with no use at all ahead (mean rank 0.49 on cardiovascular).
3. **TxGNN barely beats popularity** in ranking (zero-shot: Recall@100 0.41 against 0.36;
   mean rank 0.78 against 0.91). They did not publish it as a baseline.
4. **Zero-shot has a catch**: ~60 % of the test diseases have a near-identical sibling.
5. **Our weak spot**: when the good drugs have no other use (diabetes: insulins,
   gliptins, gliflozins…), without popularity we get 0.08 Recall@100. Masked training
   and the extra sources raise it to 0.23, but this is where we are weakest.

## Limitations

- Areas: a single seed and an approximate split (the same for both); diabetes (2 diseases)
  and adrenal (6) are very noisy.
- The improvements (`new`, `extra`, `--rounds`) were designed after having seen results on the areas
  (in good faith, and all are reported). The `base` version was already winning on Recall@100 before any improvement.
  The definitive test is a future exam (indications later than PrimeKG).
- Ties in our model (drugs with score 0) are broken at random.
- The exam rewards rediscovery: well-studied drugs have more recorded indications.

## Next steps

1. **Novelty**: ranking restricted to drugs with some use, and removing the «wildcards» (those that
   treat many diseases). How much of our advantage is biology?
2. **Drugs with no prior use**: who is better at finding the good ones that have no other indication?
3. **Future exam**: indications approved after PrimeKG (Open Targets), which no model has seen.
4. Combine: our ranking with the TxGNN signal as one more variable.

## Additional improvements and tests (night of 10 to 11/10)

### Model: LambdaRank (now the production model)

LightGBM with a ranking objective (`--rank`), which learns to order the drugs within each
disease. Decided with zero-shot seeds 1-2 and **confirmed with data not used for
the decision** (seeds 3-5 and the 9 areas):

| | GBM | **LambdaRank** | TxGNN |
|---|---|---|---|
| Zero-shot seeds 3-5: Recall@100 | 0.801 | **0.828** | 0.440 |
| Zero-shot seeds 3-5: Precision@K | 0.484 | **0.584** | 0.095 |
| Hard (seeds 3-5): Precision@K | 0.200 | **0.331** | 0.062 |
| Areas where Recall@100 improves | — | 8 of 9 (not diabetes) | — |

With LambdaRank we also beat TxGNN on AUPRC in cardiovascular (0.775 vs 0.772) and adrenal.

### Tests that did not help (do not repeat)

- **Chemical structure** (`sources/structure.py`, Morgan/Tanimoto of ChEMBL SMILES; set
  `--paths structure`): does not improve the overall ranking (0.755 vs 0.773 with GBM, seeds 1-2) nor
  drugs with no prior use. It stays in the graph (`structure` relation) in case it helps to look for derivatives.
- **Combining** the ranking with the model without popularity: +5 points on drugs with no prior use but
  −1 overall. Instead of mixing, the `disease` command gives two lists.

### Novelty (`docs/NOVELTY.md`)

Recall@100 TxGNN → ours (ranking): approved drugs only 0.60 → 0.85 (zero-shot);
without «wildcard» drugs 0.29 → 0.74. **The advantage does not come from popular drugs.**
Drugs with no prior use: 0.15 → 0.15 (zero-shot) and 0.08 → 0.16 (areas): the frontier for both.

### External exam (`docs/EXTERNAL_EXAM.md`)

3,718 approved indications in Open Targets that are not in PrimeKG (no model has seen them), same
training graph for both. Recall@100: **TxGNN 0.12, ours 0.35**; recent (first trial
≥ 2018) 0.12 → 0.33; non-obvious 0.08 → 0.17. TxGNN does not beat popularity (0.13) and its mean
rank (0.56) is almost that of chance.
