# Correcting the population size of the shared-genes test (P40)

## The problem

The shared-genes line uses a hypergeometric test: how likely is it that disease *d* (with *n_d* associated
genes) and disease *m* (with *n_m*) share *k* genes by chance, drawing from a population of *N* genes. Since
P28, both the website (`web/build.py`, `Explainer.N_genes`) and the evaluation (`evaluation/evidence.py`)
computed *N* as the number of nodes with at least one gene–disease link. Because the gene–disease matrix is
symmetric, that count includes the diseases as well as the genes (36,769 instead of the number of genes), so
every p-value was smaller than it should be, and the p-values shown on the website were wrong. Found during
T4 (P38), where «N = genes only» was one of the variants.

## Protocol (written before running)

**Fix.** *N* = number of gene nodes with at least one associated disease, the population the gene sets are
drawn from. The statistic is corrected everywhere: website, evaluation scripts and the p-values displayed.

**Threshold.** The value 10⁻⁶ was chosen for the wrong *N*, and P38 showed that the threshold carries no
biological meaning (no sharp optimum; random gene sets do as well). The aim is therefore to correct the
statistic while keeping the validated selection as far as possible, not to tune a new rule. The threshold is
chosen among {10⁻⁶, 3·10⁻⁶, 10⁻⁵, 3·10⁻⁵, 10⁻⁴} as the one whose shown shared-genes candidates (top 25, not
direct extensions) agree best with production (Jaccard of the (disease, drug) pairs), pooled over the
benchmark (zero-shot seeds 1–5) and the nine areas. The minimum of 3 shared genes is kept.

**Acceptance (external exam, not used for the choice).** The corrected rule is adopted if, on the external
exam, (a) the shared-genes hit rate stays within production's disease-level 95% interval (3.1–5.5 %) and (b) the
number of shared-genes candidates changes by less than 10 %. The pairs it removes and adds are reported with
their hit rates. If it is not accepted, production keeps the old rule and the issue stays documented.

**If adopted.** Website rebuilt in full (the shown candidates change slightly) and deployed; evaluation scripts,
`docs/EVIDENCE.md`, `docs/EVIDENCE_CI.md`, the website's Method page and help texts, the README, the preprint
(abstract, Methods, Results, Table S3, Table S5, Figure 4) and the decision log updated with the new figures.

## Results (14 October 2026, `py -3 -m evaluation.hypergeom_n`)

*N*: 36,769 linked nodes before, **19,735 genes** now (training graph of the evaluation; 19,766 in the website's
full graph).

Agreement (Jaccard of the shown shared-genes pairs) with production, benchmark and areas pooled:

| Threshold with corrected N | 10⁻⁶ | 3·10⁻⁶ | **10⁻⁵** | 3·10⁻⁵ | 10⁻⁴ |
|---|---|---|---|---|---|
| Jaccard | 0.912 | 0.942 | **0.968** | 0.966 | 0.955 |

Chosen: **p < 10⁻⁵**. Shared-genes section, hit rate (disease-level 95% interval; candidates):

| Test | Before (N = linked nodes, p < 10⁻⁶) | After (N = genes, p < 10⁻⁵) | Pairs removed | Pairs added |
|---|---|---|---|---|
| Benchmark | 22.7 % (17.3–28.3; 1,992) | 23.0 % (17.6–28.8; 1,970) | 23, none correct | 1, not correct |
| Areas | 25.2 % (20.9–29.9; 2,725) | 26.3 % (21.9–30.8; 2,608) | 123, none correct | 6, none correct |
| External exam | 4.2 % (3.1–5.5; 4,119) | 4.3 % (3.1–5.6; 3,936) | 191, 5 correct (2.6 %) | 8, none correct |

Acceptance: external hit rate within production's interval (yes) and candidates changed by 4.4 % (< 10 %):
**adopted**. The pairs removed are below the 5 % of the P29 rule on every test, so the change is also consistent
with it. Production code: `web/build.py` (`MAX_P_SHARED = 1e-5`, `Explainer.N_genes` counts gene nodes only),
`evaluation/evidence.py` (`MAX_P`, `n_genes`), `evaluation/evidence_ci.py` (base rule; the old rule kept as the
variant «before P40»). The website was rebuilt in full; all evidence figures below and in the other documents are
recomputed with the corrected rule.
