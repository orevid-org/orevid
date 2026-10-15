# Knowledge-graph embedding baselines (I5) — in progress

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

## Still to run

`bash evaluation/kge_bench.sh final` in WSL (seeds 2–5 and the nine areas for each model, 39 runs), then
`py -3 benchmark/bench.py evaluate --pred benchmark/predictions/kge_<model>`, `benchmark/cluster_ci.py`, the
leaderboard and this report; DECISIONS P47; the preprint paragraph that says no embedding baseline was run.
