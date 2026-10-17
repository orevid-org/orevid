# Knowledge-graph embedding baselines (I5)

Protocol: `docs/IMPROVEMENT_PLAN.md` (I5) and `evaluation/kge_bench.py` (details fixed before running: batch size
4,096, Adam, one «basic» negative per positive, each model's default PyKEEN loss, seed 0). Training edges: each
task's PrimeKG training edges in TxGNN's directed form (3,367,655 for zero_shot_1). Score: (drug, indication,
disease). Trained in WSL on an RTX 3060 (TransE ≈ 18–27 min per run, ComplEx 30–50 min, RotatE 25–43 min).

## Grid on zero_shot_1 (selection; 15-10-2026)

Recall@100 over the 71 test diseases of zero_shot_1 (OREVID-bench metric). These numbers chose the settings, so
they are optimistic for the embeddings; the confirmation is seeds 2–5, the nine areas and the external exam.

| Model | dim 128, lr 1e-3 | dim 128, lr 1e-2 | dim 256, lr 1e-3 | dim 256, lr 1e-2 | Chosen |
|---|---|---|---|---|---|
| TransE | **0.160** | 0.095 | 0.152 | 0.091 | 128, 1e-3 |
| ComplEx | 0.003 | 0.192 | 0.050 | **0.202** | 256, 1e-2 |
| RotatE | 0.249 | 0.161 | **0.272** | 0.133 | 256, 1e-3 |

Same seed, for reference: popularity 0.334, TxGNN 0.373, OREVID 0.843.

External exam (trained on zero_shot_1 with the chosen settings, not used for the choice): TransE 0.080, ComplEx
0.117, RotatE 0.156 (popularity 0.130, TxGNN 0.121, OREVID 0.359 on the leaderboard).

## Deviation from the plan (decided 16-10-2026, after the grid, before any confirmation run)

The plan said to run the three models unchanged on seeds 2–5, the areas and the external exam. Only **RotatE**
was run: TransE and ComplEx were below RotatE on seed 1 and on the external exam with their chosen settings, and
the final runs cost about 45 min of GPU each (the GPU, not the CPU, is the bottleneck for RotatE at dimension
256). TransE and ComplEx are reported from the grid only.

## RotatE on every task (dimension 256, learning rate 1e-3)

Recall@100 with 95 % intervals resampling diseases (`benchmark/CLUSTER_CI.md`; full metrics in
`benchmark/LEADERBOARD.md`). Seed 1 chose the settings.

| Test | Popularity | TxGNN | RotatE | OREVID |
|---|---|---|---|---|
| Zero-shot (5 seeds) | 0.358 (0.314–0.402) | 0.412 (0.369–0.457) | 0.289 (0.250–0.328) | 0.811 (0.772–0.847) |
| Zero-shot, hard | 0.203 (0.150–0.259) | 0.258 (0.197–0.324) | 0.162 (0.114–0.214) | 0.539 (0.465–0.613) |
| Areas (9) | 0.185 (0.132–0.239) | 0.237 (0.181–0.295) | 0.219 (0.170–0.270) | 0.617 (0.566–0.665) |
| Areas, hard | 0.141 (0.078–0.213) | 0.141 (0.086–0.204) | 0.154 (0.099–0.216) | 0.490 (0.414–0.565) |
| External | 0.130 (0.109–0.152) | 0.121 (0.100–0.143) | 0.156 (0.134–0.180) | 0.359 (0.327–0.390) |
| External, non-obvious | 0.075 (0.058–0.094) | 0.085 (0.066–0.105) | 0.127 (0.105–0.150) | 0.178 (0.152–0.206) |
| External, recent | 0.097 (0.045–0.154) | 0.119 (0.061–0.184) | 0.149 (0.084–0.219) | 0.355 (0.266–0.448) |

Per seed: 0.272, 0.209, 0.368, 0.307, 0.288. Novelty scenario C (drugs with no prior use): 0.046 zero-shot,
0.012 areas.

Paired differences per disease (exploratory, not part of the Holm-adjusted family of `CLUSTER_CI.md`; Wilcoxon p
unadjusted):

| Test | Comparison | Difference (95% CI) | Ahead | Behind | p |
|---|---|---|---|---|---|
| zero-shot | RotatE − TxGNN | -0.123 (-0.169 to -0.080) | 17 % | 39 % | 0.000 |
| zero-shot | RotatE − popularity | -0.069 (-0.119 to -0.017) | 25 % | 34 % | 0.010 |
| zero-shot, hard | RotatE − TxGNN | -0.096 (-0.158 to -0.037) | 15 % | 28 % | 0.004 |
| zero-shot, hard | RotatE − popularity | -0.041 (-0.110 to +0.028) | 23 % | 29 % | 0.221 |
| areas | RotatE − TxGNN | -0.018 (-0.064 to +0.027) | 27 % | 23 % | 0.459 |
| areas | RotatE − popularity | +0.035 (-0.012 to +0.080) | 34 % | 16 % | 0.001 |
| areas, hard | RotatE − TxGNN | +0.013 (-0.036 to +0.059) | 24 % | 19 % | 0.169 |
| areas, hard | RotatE − popularity | +0.013 (-0.051 to +0.072) | 30 % | 12 % | 0.018 |
| external, all | RotatE − TxGNN | +0.036 (+0.007 to +0.065) | 24 % | 13 % | 0.004 |
| external, all | RotatE − popularity | +0.026 (-0.001 to +0.054) | 23 % | 16 % | 0.040 |
| external, non_obvious | RotatE − TxGNN | +0.042 (+0.014 to +0.070) | 21 % | 10 % | 0.001 |
| external, non_obvious | RotatE − popularity | +0.051 (+0.024 to +0.079) | 21 % | 10 % | 0.000 |
| external, recent | RotatE − TxGNN | +0.030 (-0.056 to +0.116) | 14 % | 9 % | 0.519 |
| external, recent | RotatE − popularity | +0.052 (-0.027 to +0.134) | 12 % | 9 % | 0.186 |

**Reading.** RotatE is below popularity and TxGNN on the zero-shot split, where test diseases have lost every
drug edge and an embedding has little to place them with; it is level with TxGNN on the areas; and it is above
both on the external exam, especially on non-obvious approvals (+0.042 over TxGNN, interval +0.014 to +0.070).
OREVID is far above all three everywhere.

## OREVID + RotatE ensemble (P33 rule)

`py -3 -m evaluation.kge_ensemble`: z(OREVID) + w · z(RotatE) per disease, w chosen on seeds 1–2.

| Group | OREVID | + 0.1 · RotatE | + 0.25 · RotatE | + 0.5 · RotatE | + 1.0 · RotatE |
|---|---|---|---|---|---|
| selection (seeds 1-2), all | 0.790 | 0.789 | 0.789 | 0.792 | 0.818 |
| selection (seeds 1-2), hard | 0.511 | 0.507 | 0.507 | 0.514 | 0.581 |
| seeds 3-5, all | 0.825 | 0.823 | 0.824 | 0.821 | 0.816 |
| seeds 3-5, hard | 0.558 | 0.554 | 0.557 | 0.559 | 0.558 |
| areas, all | 0.617 | 0.617 | 0.610 | 0.606 | 0.604 |
| areas, hard | 0.490 | 0.491 | 0.484 | 0.481 | 0.484 |
| external, all | 0.359 | 0.366 | 0.372 | 0.372 | 0.370 |
| external, non_obvious | 0.178 | 0.187 | 0.194 | 0.190 | 0.195 |
| external, recent | 0.355 | 0.355 | 0.338 | 0.328 | 0.365 |

Chosen on seeds 1-2: w = 1.0. Gains in the confirmation groups: seeds 3-5 -0.008, areas -0.012, external +0.012 (mean -0.003); hard / non-obvious changes: seeds 3-5 hard -0.000, areas hard -0.007, external non_obvious +0.016.

**Decision (P33 rule):** not adopted; production unchanged.

The weight that looked best on seeds 1–2 (+0.028) lost on seeds 3–5 and the areas and helped only on the external
exam (+0.012, non-obvious +0.016): the rule is not met.
