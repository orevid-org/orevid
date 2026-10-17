# OREVID — Open Repurposing Evidence

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23323226.svg)](https://doi.org/10.5281/zenodo.23323226)

OREVID ranks existing drugs for each disease using a biomedical knowledge graph (drugs, genes,
diseases, phenotypes, trials, adverse effects) and shows **the evidence behind every candidate**.
It only shows candidates backed by evidence that has been checked against held-out data, and it says
how often each kind of evidence was right.

- Website: **https://orevid.org** (English, Spanish, Catalan)
- Benchmark: [`benchmark/`](benchmark/) — OREVID-bench, a public and reproducible drug-repurposing
  benchmark on PrimeKG, with TxGNN run with its official code on the same splits
- Licences: code MIT, data and results CC BY-NC-SA 4.0 (see [Licences](#licences))

> **Research use only.** OREVID generates hypotheses for researchers. It is not a medical device,
> not medical advice, and must not be used to make decisions about individual patients. The data
> contain known errors, and suggestions may be wrong or harmful.

## What it does, in two sentences

For each disease, OREVID orders the ~8,000 drugs in PrimeKG by how likely they are to be useful, using
paths in a knowledge graph, and keeps only those supported by checked evidence: mainly, the drug
already treats a neighbouring disease in the ontology, or a disease that significantly shares
associated genes with this one. Validated against TxGNN (Harvard, *Nature Medicine* 2024) on the same
data, it recovers about twice as many known treatments; for genuinely new indications about 1 in
20–25 of its table candidates proves right, so it is a tool to prioritise research hypotheses.

## Results

Recall@100 = fraction of a disease's correct drugs found among the top 100 of 7,957 candidates
(mean over diseases). Full tables: [`benchmark/LEADERBOARD.md`](benchmark/LEADERBOARD.md).

| Test | Popularity baseline | TxGNN | RotatE (embedding) | OREVID |
|---|---|---|---|---|
| Zero-shot diseases (5 seeds) | 0.358 | 0.412 | 0.289 | **0.811** |
| Same, without near-identical diseases («hard») | 0.203 | 0.258 | 0.162 | **0.539** |
| Whole disease areas hidden (9) | 0.185 | 0.237 | 0.219 | **0.617** |
| External: Open Targets approvals not in PrimeKG | 0.130 | 0.121 | 0.156 | **0.359** |

RotatE: knowledge-graph embedding (PyKEEN), settings chosen on zero-shot seed 1 ([`docs/KGE.md`](docs/KGE.md)).

What the website shows, and how often it was right among each disease's top 25
([`docs/EVIDENCE.md`](docs/EVIDENCE.md)):

| Section | Hidden-treatment tests | New approvals (external exam) |
|---|---|---|
| Direct extensions (indicated for a neighbouring disease) | 13–59 % | 12.6 % |
| Candidates table (≥ 3 shared genes, hypergeometric p < 10⁻⁵ over the 19,766 genes with a disease) | 22–27 % | 4.3 % |
| Hidden (none of the above) | 2–4 % | 1.5 % |

Chance would be below 0.1 %. Random gene sets of the same size select candidates that are right as often as the real shared genes, so the shared genes are not evidence of a shared mechanism: the table works as a selection of drugs already used in comparable diseases ([`docs/EVIDENCE_CI.md`](docs/EVIDENCE_CI.md), P38). The hidden-treatment test is the usual one (also TxGNN's) but rewards
rediscovering what is known; the external exam is the more realistic one. The two protocol differences
with TxGNN are measured, not only stated, in [`benchmark/FAIR_COMPARISON.md`](benchmark/FAIR_COMPARISON.md):
with OREVID trained on the training labels only (as TxGNN) and the diseases where TxGNN outputs a
constant score removed, zero-shot Recall@100 is **0.809 vs 0.421** (hard 0.528 vs 0.275; areas 0.627 vs 0.238).
Intervals that resample whole diseases are in [`benchmark/CLUSTER_CI.md`](benchmark/CLUSTER_CI.md): the
paired difference with TxGNN is +0.40 (95% CI 0.35–0.44) on zero-shot and +0.24 (0.21–0.27) on the external exam.

**Test against the future** ([`docs/TEMPORAL.md`](docs/TEMPORAL.md), pre-registered in
[`docs/IMPROVEMENT_PLAN.md`](docs/IMPROVEMENT_PLAN.md)): with the graph and the model as known in 2015, the uses
that later reached phase 3 or approval are ranked well above popularity (Recall@100 0.322 vs 0.174, paired
+0.148, 95% CI 0.123–0.173; development cutoff 2012 agrees). The advantage on non-obvious uses is small (+0.043),
and 80 % of PrimeKG indications are undated, so some leakage remains. A calibrated probability per candidate was
tested and not adopted ([`docs/CALIBRATION.md`](docs/CALIBRATION.md)). Each candidate on the website also shows
whether a generic or biosimilar exists and how much the pair has been studied, as information only
([`docs/PRACTICAL.md`](docs/PRACTICAL.md)).
Two pre-registered attempts gave negative results: expression-signature reversal (LINCS L1000) ranks drugs
with no prior use at chance ([`docs/LINCS.md`](docs/LINCS.md)), and no replacement for the shared-genes line beat
the permutation null on the external exam ([`docs/BIOLOGY_LINES.md`](docs/BIOLOGY_LINES.md)).

## How it works

1. **Graph.** PrimeKG plus Open Targets 26.06 (mechanisms with direction, clinical stages, stopped
   trials, gene–disease evidence), STRING, PubTator3 and others: 157,822 nodes and 9,457,076
   relations (`prepare.py`, `sources/`).
2. **Features.** For each disease–drug pair, degree-weighted path counts (DWPC, as in Hetionet) over
   «metapaths» such as *drug targets a disease gene* or *drug treats a disease sharing genes*
   (`graph.py`, `evaluation/txgnn.py`).
3. **Ranking.** A LightGBM LambdaRank model with masked training (the disease's own treatments are
   hidden while learning), then reordered by specificity so drugs that score high everywhere do not
   dominate (`predictor.py`; decisions P13, P27–P29 in [`docs/DECISIONS.md`](docs/DECISIONS.md)).
4. **Evidence filter.** Of the top 25, only candidates with checked evidence are shown; every rule was
   fixed in advance and confirmed on an independent split ([`evaluation/evidence.py`](evaluation/evidence.py)).
   «Acts on a disease gene» on its own is right 2–3 % of the time, the same as no evidence, so it does
   not count.
5. **Risk warnings** are a transparent rule, not a model: contraindicated in more related diseases
   than indicated, or an immune-checkpoint blocker in an autoimmune disease
   ([`docs/DIRECTION.md`](docs/DIRECTION.md)). No warning does not mean safe.

Other documents: [`docs/TARGETS.md`](docs/TARGETS.md) (which gene to target),
[`docs/NOVELTY.md`](docs/NOVELTY.md) and [`docs/EXTERNAL_EXAM.md`](docs/EXTERNAL_EXAM.md) (novelty tests),
[`docs/TXGNN_EVALUATION.md`](docs/TXGNN_EVALUATION.md) (how TxGNN was reproduced),
[`docs/PANEL.md`](docs/PANEL.md) (automatic quality panel over all diseases),
[`docs/XGBOOST.md`](docs/XGBOOST.md) (XGBoost, CatBoost and ensembles against the production LightGBM),
[`docs/PROSPECTIVE_PLAN.md`](docs/PROSPECTIVE_PLAN.md) (prospective validation, snapshot of 12 October 2026),
[`resources/opentargets/`](resources/opentargets/) and `sources/ot_suspect_approvals.py` (suspect Open Targets approvals, opentargets/issues#4579),
[`docs/REVIEW_RESPONSE.md`](docs/REVIEW_RESPONSE.md) (answer to an external review, with [`docs/EVIDENCE_CI.md`](docs/EVIDENCE_CI.md), [`docs/BIOMARKERS.md`](docs/BIOMARKERS.md) and [`docs/EXTERNAL_SEEDS.md`](docs/EXTERNAL_SEEDS.md)),
[`docs/preprint/`](docs/preprint/) (preprint draft) and [`docs/STATUS.md`](docs/STATUS.md) (status and next steps).

## Data sources

| Source | Licence | What it contributes |
|---|---|---|
| **PrimeKG** (Chandak et al., *Sci Data* 2023) | distributed CC0 by Harvard; contains DrugBank (CC BY-NC 4.0), SIDER and DisGeNET (CC BY-NC-SA 4.0) and other sources | The base graph: 129k nodes, 4 M relations from 20 sources |
| **Open Targets Platform 26.06** | CC0 (includes ChEMBL mechanisms and molecules) | Mechanisms with direction, clinical stages and year, stopped trials, FAERS signals, gene–disease evidence by type |
| **STRING 12.5** | CC BY 4.0 | Human protein interactions (confidence ≥ 700) |
| **PubTator3** (NCBI) | public domain | Relations mined from all of PubMed, shown only as context (never used for scoring or validation) |
| **HGNC, MONDO, HPO** | CC0 / CC BY 4.0 / HPO licence | Identifier mapping and ontologies |
| **CTD** | non-commercial | Used in analyses only (decision P20); nothing from CTD is in the published model or data |
| **FDA Orange Book and Purple Book** | US public domain | Generic / biosimilar status, US patent or exclusivity date and routes, display only (P44) |
| **EMA medicines table** | © European Medicines Agency, reuse with acknowledgement | EU generic / biosimilar status (centralised procedure only), display only (P44) |

Pitfalls found in the data and how they are handled are recorded in [`docs/DECISIONS.md`](docs/DECISIONS.md).
Two were reported upstream to Open Targets: multiple sclerosis merged into a measurement term in 26.09
(that is why 26.06 is used), and approvals derived from multi-drug labels or abbreviations (handled by
[`sources/ot_approvals.py`](sources/ot_approvals.py)).

## Quick start

Requirements: Python 3.10+ and the packages in `requirements.txt` (`pip install -r requirements.txt`);
about 7 GB of disk for the raw data. Commands below use `python`; on Windows use `py -3`.

```bash
python download.py --all              # PrimeKG and TxGNN's kg.csv/node.csv (Harvard Dataverse, ~1.5 GB)
python -m sources.download            # Open Targets, HGNC, MONDO, STRING, CTD (~4 GB)
# PubTator3: download https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/relation2pubtator3.gz to data/raw/pubtator/
bash build_all.sh                     # graph, models, evaluations, benchmark, website, tests (2–3 h)
```

`build_all.sh` lists every step in order. The TxGNN predictions used in the comparison come from a
separate GPU run (`evaluation/train_txgnn.sh`, Linux/WSL with CUDA; about 22 hours for all splits); the
benchmark already includes TxGNN's scores (`benchmark/predictions/txgnn/results.json`).

Exploring from the command line:

```bash
python explore.py disease "multiple sclerosis" --top 20   # the same page as the website
python explore.py targets "multiple sclerosis"            # which gene or protein to target
python explore.py drug "ocrelizumab"
python explore.py path "vitamin d" "multiple sclerosis"
```

Website locally: `python -m web.build --test` (24 diseases) or `python -m web.build` (all, ~1 h), then
`python web/serve.py` and open http://localhost:8796. Tests: `python -m pytest tests -q` (most need the
built data and are skipped without it). Publishing the website: [`docs/DEPLOY.md`](docs/DEPLOY.md).

## Repository layout

| Path | Content |
|---|---|
| `download.py`, `sources/` | Downloads and converts each source to `data/ext/` |
| `prepare.py`, `graph.py` | Merges everything into the graph (`data/prep/`) and computes metapath features |
| `predictor.py` | Production model and drug status/filters |
| `evaluation/` | Comparison with TxGNN, novelty tests, external exam, evidence rules, targets, direction |
| `benchmark/` | OREVID-bench: tasks, evaluation code, predictions' metrics, leaderboard |
| `web/` | Website builder (`build.py`), quality panel, local server, static site in `web/public/` |
| `resources/` | Small inputs: TxGNN disease-area lists and published results, a PubMed extraction pilot |
| `tests/` | Integrity tests |

## Licences

- **Code**: MIT ([`LICENSE`](LICENSE)).
- **Data and results** (website data, benchmark tasks and predictions, any ranking produced from the
  sources): **CC BY-NC-SA 4.0**, because they derive from sources with non-commercial and share-alike
  terms. Details and third-party notices: [`LICENSE-DATA.md`](LICENSE-DATA.md),
  [`THIRD_PARTY_NOTICES.md`](THIRD_PARTY_NOTICES.md).

## Citing

If you use OREVID or OREVID-bench, please cite it (DOI [10.5281/zenodo.23323226](https://doi.org/10.5281/zenodo.23323226),
all versions; see [`CITATION.cff`](CITATION.cff)) and the original sources, in particular PrimeKG and TxGNN for the
benchmark.

## Reporting errors

Wrong or dangerous suggestions, data errors and questions: open an issue at
https://github.com/orevid-org/orevid/issues or write to contact@orevid.org.

## Known limitations

- The «truth» in every test is recorded approvals: a candidate counted as wrong may be right but not
  yet approved, and known treatments are favoured.
- Effect direction cannot be predicted for new drugs; the risk warning relies on recorded
  contraindications and does not detect opposite mechanisms or class toxicities.
- Source errors propagate (e.g. wrong gene–disease associations, comorbidities read as shared
  biology); the website shows where every candidate comes from so they can be spotted.
- Drugs with no prior use remain hard for every model tested (Precision@K 0–2 %).
