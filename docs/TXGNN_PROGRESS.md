# Progress: TxGNN against our model

Updated: 10/10 23:01. Regenerated automatically (`evaluation/txgnn_progress.py`).

Same split and same metric code for both (see `docs/TXGNN_EVALUATION.md`). «Ours» = gradient boosting with extra sources and masked training. «Hard» = test diseases with no near-identical sibling with known drugs.

## Recall@100: of every 100 good drugs, how many appear among the first 100

| Split | Diseases (hard) | Popularity | TxGNN → ours | Hard: TxGNN → ours | TxGNN published |
|---|---|---|---|---|---|
| Zero-shot (seed 1) | 71 (29) | 0.33 | 0.37 → **0.82** | 0.20 → **0.58** | 0.40 |
| Mental health (seed 1) | 39 (33) | 0.02 | 0.21 → **0.43** | 0.20 → **0.35** | 0.18 |
| Diabetes (seed 1) | 2 (2) | 0.00 | 0.03 → **0.23** | 0.03 → **0.23** | 0.24 |
| Neurodegenerative (seed 1) | 12 (6) | 0.06 | 0.48 → **0.81** | 0.19 → **0.62** | 0.54 |
| Cardiovascular (seed 1) | 53 (49) | 0.12 | 0.11 → **0.34** | 0.09 → **0.29** | 0.09 |
| Anaemia (seed 1) | 12 (6) | 0.24 | 0.08 → **0.67** | 0.00 → **0.44** | 0.29 |
| Autoimmune (seed 1) | 13 (6) | 0.36 | 0.33 → **0.72** | 0.29 → **0.47** | 0.34 |
| Adrenal (seed 1) | 6 (4) | 0.52 | 0.47 → **0.91** | 0.23 → **0.89** | 0.50 |
| Metabolic (seed 1) | 28 (19) | 0.01 | 0.08 → **0.53** | 0.01 → **0.32** | 0.21 |
| Cancer (seed 1) | 180 (108) | 0.32 | 0.35 → **0.73** | 0.23 → **0.57** | 0.36 |
| Zero-shot (seed 2) | 64 (27) | 0.31 | 0.37 → **0.72** | 0.18 → **0.36** | 0.40 |
| Zero-shot (seed 3) | 66 (25) | 0.38 | 0.43 → **0.80** | 0.30 → **0.49** | 0.40 |
| Zero-shot (seed 4) | 66 (22) | 0.41 | 0.47 → **0.82** | 0.27 → **0.53** | 0.40 |
| Zero-shot (seed 5) | 68 (26) | 0.36 | 0.42 → **0.79** | 0.34 → **0.58** | 0.40 |

## Other metrics (TxGNN → ours; best in bold)

AUPRC and AUROC: positives against the same number of random negatives among drugs with some indication. Precision@K: hits among the first K (K = good drugs). Mean rank: 1 = at the top.

| Split | AUPRC | AUROC | Precision@K | Mean rank | AUPRC hard |
|---|---|---|---|---|---|
| Zero-shot (seed 1) | 0.93 → **0.97** | 0.93 → **0.95** | 0.06 → **0.53** | 0.77 → **0.97** | 0.86 → **0.90** |
| Mental health (seed 1) | **0.92** → 0.86 | **0.92** → 0.80 | 0.09 → **0.17** | 0.68 → **0.89** | **0.94** → 0.84 |
| Diabetes (seed 1) | 0.69 → **0.70** | **0.80** → 0.64 | 0.00 → **0.18** | 0.69 → **0.86** | 0.69 → **0.70** |
| Neurodegenerative (seed 1) | **0.98** → 0.91 | **0.98** → 0.91 | 0.21 → **0.61** | 0.90 → **0.96** | **0.97** → 0.88 |
| Cardiovascular (seed 1) | **0.77** → 0.74 | **0.77** → 0.70 | 0.01 → **0.14** | 0.49 → **0.89** | **0.77** → 0.73 |
| Anaemia (seed 1) | 0.65 → **0.98** | 0.74 → **0.98** | 0.00 → **0.36** | 0.53 → **0.97** | 0.64 → **0.94** |
| Autoimmune (seed 1) | **0.91** → 0.90 | **0.91** → 0.84 | 0.11 → **0.40** | 0.74 → **0.93** | 0.79 → **0.80** |
| Adrenal (seed 1) | **0.98** → 0.97 | **0.98** → 0.94 | 0.07 → **0.50** | 0.58 → **0.99** | **0.97** → 0.97 |
| Metabolic (seed 1) | **0.93** → 0.90 | **0.92** → 0.84 | 0.02 → **0.21** | 0.64 → **0.91** | **0.90** → 0.84 |
| Cancer (seed 1) | **0.95** → 0.95 | **0.96** → 0.92 | 0.02 → **0.48** | 0.81 → **0.95** | **0.95** → 0.90 |
| Zero-shot (seed 2) | 0.85 → **0.94** | 0.84 → **0.91** | 0.05 → **0.44** | 0.74 → **0.97** | 0.70 → **0.77** |
| Zero-shot (seed 3) | 0.92 → **0.96** | 0.93 → **0.95** | 0.10 → **0.50** | 0.82 → **0.97** | 0.86 → **0.88** |
| Zero-shot (seed 4) | 0.93 → **0.96** | 0.93 → **0.94** | 0.09 → **0.53** | 0.83 → **0.97** | 0.87 → **0.89** |
| Zero-shot (seed 5) | 0.92 → **0.95** | 0.92 → **0.93** | 0.09 → **0.42** | 0.76 → **0.97** | **0.89** → 0.89 |

## TxGNN queue

| Split | Status |
|---|---|
| Zero-shot (seed 1) | ✅ finished (02:19:58) |
| Mental health (seed 1) | ✅ finished (03:53:44) |
| Diabetes (seed 1) | ✅ finished (05:27:01) |
| Neurodegenerative (seed 1) | ✅ finished (06:59:16) |
| Cardiovascular (seed 1) | ✅ finished (08:33:36) |
| Anaemia (seed 1) | ✅ finished (10:07:44) |
| Autoimmune (seed 1) | ✅ finished (11:46:30) |
| Adrenal (seed 1) | ✅ finished (13:18:47) |
| Metabolic (seed 1) | ✅ finished (14:53:03) |
| Cancer (seed 1) | ✅ finished (16:31:35) |
| Zero-shot (seed 2) | ✅ finished (18:08:36) |
| Zero-shot (seed 3) | ✅ finished (19:43:08) |
| Zero-shot (seed 4) | ✅ finished (21:20:18) |
| Zero-shot (seed 5) | ✅ finished (22:56:09) |
