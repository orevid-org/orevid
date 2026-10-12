# XGBoost against LightGBM for the ranking model

Question (12 October 2026): would XGBoost rank drugs better than the production LightGBM LambdaRank
model? On many tabular problems XGBoost does better than LightGBM, so it is worth measuring here
instead of assuming either way.

Script: `evaluation/xgboost_compare.py`. Raw results: `data/eval/xgb/results.csv` (not in git;
tables below).

## Protocol (written before running)

**What stays the same.** Same training rows as `evaluate()` in `evaluation/txgnn.py`: 1,500
training diseases (all positives, 10% of negatives), the five masked-neighbourhood rounds, the
same 22 metapath features plus popularity and drug degrees (`extra_no_ctd`), the same test
diseases and the same metrics (TxGNN notebook: Recall@100, Precision@K, AUPRC with sampled
negatives; external exam as in `evaluation/external_exam.py`). Only the learner changes. As a
check, the production configuration retrained by the script must reproduce the published OREVID
scores.

**Candidates.** 17 configurations, the same kind of small grid for both libraries so that neither
gets more tuning than the other. All use learning rate 0.05 and L2 regularisation 1.

| Library | Objective | Grid |
|---|---|---|
| LightGBM | LambdaRank | 300 trees × 31 leaves (production), 600 × 31, 300 × 15, 300 × 63, 600 × 63 |
| LightGBM | binary | 300 × 31 |
| XGBoost | rank:ndcg | 300 or 600 trees × depth 4, 6, 8; 600 × depth 6 with row and column subsampling 0.8 |
| XGBoost | rank:pairwise | 300 × depth 6, 600 × depth 8 |
| XGBoost | binary:logistic | 300 × depth 6, 600 × depth 8 |

LightGBM keeps `min_child_samples=50`; XGBoost keeps its defaults (`min_child_weight=1`, NDCG pairs
by top-k). XGBoost runs on the GPU (RTX 3060), LightGBM on the CPU.

**Selection.** On zero-shot seeds 1–2 only: the configuration with the highest mean Recall@100
(ties broken by Recall@100 on the hard subset), separately for LightGBM and for XGBoost.

**Confirmation.** The production model, the best LightGBM and the best XGBoost on data not used
for the selection: zero-shot seeds 3–5, the nine disease areas and the external exam (trained on
zero-shot seed 1, like the published one).

**Decision rule.** XGBoost replaces LightGBM in production only if its best configuration beats the
production model on mean Recall@100 in all three confirmation groups (seeds 3–5, areas, external),
by at least 0.01 on average over the three, and loses no more than 0.01 on any hard or non-obvious
subset. Otherwise production stays as it is. The same rule applies to a tuned LightGBM. Any change
of production model would also require re-running the evidence hit rates (`evaluation/evidence.py`),
the website build and the benchmark export, and is left to the maintainer's decision.

**Second round (added after the first tuning, before any confirmation run).** The first tuning
showed no XGBoost configuration above the production LightGBM, so two other families that often
win at learning to rank were added, with the same selection on seeds 1–2 and the same decision
rule: CatBoost with ranking losses (YetiRank: 600 trees × depth 6 or 8, 1,200 × depth 6, 600 × depth 6
with learning rate 0.1; PairLogit: 600 × depth 6; all on the GPU, L2 1) and rank-averaged
ensembles (mean of per-disease rank percentiles) of the production LightGBM with the best XGBoost
of the first round, with CatBoost YetiRank 600 × 6, with both, and with the logistic regression
without popularity (the published «graph only» model). The best of each family joins the
confirmation.

## Results (run on 12–13 October 2026)

**Check.** The production configuration retrained by the script reproduced the published OREVID
scores exactly (maximum absolute difference 0 on all 14 splits and on the external exam), so the
comparison runs on exactly the published setting.

**Selection (seeds 1–2).** Best LightGBM: LambdaRank with 15 leaves (0.792, production 0.790).
Best XGBoost: rank:pairwise, 300 trees, depth 6 (0.783). Best CatBoost: PairLogit, 600 × depth 6
(0.776). Best ensemble: production LightGBM + logistic regression without popularity (0.785). XGBoost
with the NDCG objective (0.759–0.779) and binary objectives (0.732–0.758) were weaker, and longer or
deeper XGBoost models tended to overfit (binary 600 × depth 8: 0.651 on seed 2).

**Confirmation (seeds 3–5, areas, external exam), mean Recall@100.**

| Group | Production LightGBM | LightGBM 15 leaves | XGBoost pairwise | CatBoost PairLogit | LightGBM + logistic |
|---|---|---|---|---|---|
| Zero-shot seeds 3–5 | **0.825** | 0.817 | 0.810 | 0.812 | 0.817 |
| Zero-shot seeds 3–5, hard | **0.558** | 0.551 | 0.533 | 0.544 | 0.552 |
| Areas (mean of 9) | **0.617** | 0.610 | 0.612 | 0.560 | **0.617** |
| Areas, hard | 0.490 | 0.485 | 0.479 | 0.411 | **0.496** |
| External exam | 0.359 | 0.360 | 0.348 | 0.336 | **0.374** |
| External, non-obvious | 0.178 | 0.175 | 0.156 | 0.152 | **0.196** |
| External, recent | 0.355 | **0.365** | 0.323 | 0.319 | 0.326 |
| Precision@K, zero-shot 3–5 | **0.602** | 0.594 | 0.505 | 0.529 | 0.569 |
| AUPRC, zero-shot 3–5 | 0.953 | 0.951 | 0.950 | **0.956** | 0.947 |

**Decision: production stays as it is (LightGBM LambdaRank, 300 trees, 31 leaves).** No candidate
met the pre-specified rule (better in all three confirmation groups, +0.01 on average, no hard or
non-obvious subset worse by more than 0.01):

- XGBoost was worse in all three groups (−0.015 zero-shot, −0.005 areas, −0.011 external) and
  clearly worse in Precision@K (0.505 against 0.602). On this problem LightGBM's LambdaRank is the
  better learner; XGBoost's pairwise objective comes closest.
- CatBoost was worse, especially on the areas (0.560).
- The tuned LightGBM (15 leaves) was within noise of production, slightly lower on two groups.
- The ensemble of the production model with the logistic regression without popularity tied on
  the areas and was better on the external exam (+0.015; non-obvious +0.018), but worse on the
  zero-shot seeds (−0.008), on recent approvals (−0.029) and in Precision@K. It does not meet the
  rule. Its gain on non-obvious external approvals, the closest proxy for new uses, is worth a
  dedicated, pre-specified test later (for example on the prospective outcomes).

Full tables (every configuration and task): below, generated by
`py -3 -m evaluation.xgboost_compare report`.

### Tuning: zero-shot seeds 1-2 (mean)

| Model | Recall@100 | Recall@100 hard | Precision@K | AUPRC | Fit (s) |
|---|---|---|---|---|---|
| lgb_rank_n300_l15 | 0.792 | 0.523 | 0.537 | 0.950 | 19 |
| lgb_rank_n300_l31 | 0.790 | 0.511 | 0.529 | 0.952 | 21 |
| lgb_rank_n600_l31 | 0.788 | 0.511 | 0.530 | 0.952 | 38 |
| lgb_rank_n600_l63 | 0.787 | 0.512 | 0.517 | 0.953 | 44 |
| ens_lgb_logit | 0.785 | 0.520 | 0.526 | 0.945 | 6 |
| lgb_rank_n300_l63 | 0.783 | 0.493 | 0.519 | 0.955 | 23 |
| xgb_pairwise_n300_d6 | 0.783 | 0.499 | 0.474 | 0.952 | 20 |
| xgb_pairwise_n600_d8 | 0.780 | 0.494 | 0.465 | 0.952 | 42 |
| ens_lgb_xgb | 0.779 | 0.484 | 0.485 | 0.959 | 20 |
| xgb_ndcg_n600_d8 | 0.779 | 0.501 | 0.472 | 0.951 | 44 |
| lgb_binary_n300_l31 | 0.778 | 0.481 | 0.511 | 0.957 | 10 |
| cb_pairlogit_n600_d6 | 0.776 | 0.492 | 0.476 | 0.950 | 60 |
| xgb_ndcg_n300_d8 | 0.774 | 0.488 | 0.456 | 0.953 | 23 |
| cb_yetirank_n1200_d6 | 0.774 | 0.479 | 0.490 | 0.948 | 66 |
| xgb_ndcg_n600_d4 | 0.772 | 0.475 | 0.488 | 0.947 | 39 |
| ens_lgb_xgb_cb | 0.772 | 0.473 | 0.500 | 0.957 | 0 |
| xgb_ndcg_n300_d6 | 0.770 | 0.472 | 0.465 | 0.953 | 21 |
| xgb_ndcg_n300_d4 | 0.768 | 0.463 | 0.495 | 0.948 | 20 |
| cb_yetirank_n600_d6_lr10 | 0.768 | 0.463 | 0.491 | 0.947 | 34 |
| cb_yetirank_n600_d8 | 0.767 | 0.461 | 0.399 | 0.948 | 36 |
| xgb_ndcg_n600_d6_sub | 0.766 | 0.465 | 0.463 | 0.952 | 42 |
| ens_lgb_cb | 0.763 | 0.452 | 0.490 | 0.956 | 20 |
| xgb_ndcg_n600_d6 | 0.759 | 0.449 | 0.460 | 0.951 | 41 |
| xgb_binary_n300_d6 | 0.758 | 0.444 | 0.463 | 0.956 | 5 |
| logit_nopop | 0.758 | 0.474 | 0.478 | 0.921 | 0 |
| cb_yetirank_n600_d6 | 0.747 | 0.415 | 0.460 | 0.947 | 34 |
| xgb_binary_n600_d8 | 0.732 | 0.412 | 0.448 | 0.952 | 10 |

### Confirmation: Recall@100 per task

| Task | Subset | cb_pairlogit_n600_d6 | ens_lgb_logit | lgb_rank_n300_l15 | lgb_rank_n300_l31 | xgb_pairwise_n300_d6 |
|---|---|---|---|---|---|---|
| zero_shot_1 | all | 0.819 | 0.834 | 0.836 | 0.843 | 0.820 |
| zero_shot_1 | hard | 0.570 | 0.641 | 0.613 | 0.630 | 0.572 |
| external | all | 0.336 | 0.374 | 0.360 | 0.359 | 0.348 |
| external | non_obvious | 0.152 | 0.196 | 0.175 | 0.178 | 0.156 |
| external | recent | 0.319 | 0.326 | 0.365 | 0.355 | 0.323 |
| zero_shot_3 | all | 0.789 | 0.803 | 0.813 | 0.817 | 0.799 |
| zero_shot_3 | hard | 0.462 | 0.500 | 0.526 | 0.536 | 0.488 |
| zero_shot_4 | all | 0.834 | 0.818 | 0.816 | 0.826 | 0.816 |
| zero_shot_4 | hard | 0.590 | 0.551 | 0.539 | 0.547 | 0.530 |
| zero_shot_5 | all | 0.813 | 0.828 | 0.820 | 0.831 | 0.815 |
| zero_shot_5 | hard | 0.579 | 0.605 | 0.589 | 0.589 | 0.581 |
| area_cell_proliferation | all | 0.708 | 0.791 | 0.747 | 0.737 | 0.753 |
| area_cell_proliferation | hard | 0.527 | 0.666 | 0.593 | 0.576 | 0.603 |
| area_mental_health | all | 0.344 | 0.460 | 0.439 | 0.427 | 0.399 |
| area_mental_health | hard | 0.284 | 0.394 | 0.365 | 0.349 | 0.321 |
| area_cardiovascular | all | 0.359 | 0.368 | 0.377 | 0.378 | 0.368 |
| area_cardiovascular | hard | 0.325 | 0.324 | 0.334 | 0.335 | 0.324 |
| area_anemia | all | 0.706 | 0.692 | 0.722 | 0.704 | 0.722 |
| area_anemia | hard | 0.486 | 0.412 | 0.486 | 0.468 | 0.486 |
| area_adrenal_gland | all | 0.578 | 0.911 | 0.911 | 0.911 | 0.944 |
| area_adrenal_gland | hard | 0.392 | 0.892 | 0.892 | 0.892 | 0.917 |
| area_autoimmune | all | 0.733 | 0.732 | 0.726 | 0.726 | 0.726 |
| area_autoimmune | hard | 0.449 | 0.462 | 0.449 | 0.449 | 0.449 |
| area_metabolic_disorder | all | 0.552 | 0.656 | 0.570 | 0.582 | 0.540 |
| area_metabolic_disorder | hard | 0.365 | 0.518 | 0.391 | 0.409 | 0.347 |
| area_diabetes | all | 0.244 | 0.087 | 0.145 | 0.236 | 0.247 |
| area_diabetes | hard | 0.244 | 0.087 | 0.145 | 0.236 | 0.247 |
| area_neurodigenerative | all | 0.814 | 0.856 | 0.856 | 0.849 | 0.808 |
| area_neurodigenerative | hard | 0.628 | 0.712 | 0.712 | 0.699 | 0.615 |

### Confirmation: summary (mean Recall@100)

| Group | cb_pairlogit_n600_d6 | ens_lgb_logit | lgb_rank_n300_l15 | lgb_rank_n300_l31 | xgb_pairwise_n300_d6 |
|---|---|---|---|---|---|
| zero-shot seeds 3-5, all | 0.812 | 0.817 | 0.817 | 0.825 | 0.810 |
| zero-shot seeds 3-5, hard | 0.544 | 0.552 | 0.551 | 0.558 | 0.533 |
| areas (mean of 9), all | 0.560 | 0.617 | 0.610 | 0.617 | 0.612 |
| areas (mean of 9), hard | 0.411 | 0.496 | 0.485 | 0.490 | 0.479 |
| external, all | 0.336 | 0.374 | 0.360 | 0.359 | 0.348 |
| external, non_obvious | 0.152 | 0.196 | 0.175 | 0.178 | 0.156 |
| external, recent | 0.319 | 0.326 | 0.365 | 0.355 | 0.323 |

### Confirmation: mean Precision@K

| Group | cb_pairlogit_n600_d6 | ens_lgb_logit | lgb_rank_n300_l15 | lgb_rank_n300_l31 | xgb_pairwise_n300_d6 |
|---|---|---|---|---|---|
| zero-shot seeds 3-5 | 0.529 | 0.569 | 0.594 | 0.602 | 0.505 |
| areas (mean of 9) | 0.371 | 0.359 | 0.394 | 0.408 | 0.379 |
| external | 0.081 | 0.096 | 0.094 | 0.090 | 0.081 |

### Confirmation: mean AUPRC

| Group | cb_pairlogit_n600_d6 | ens_lgb_logit | lgb_rank_n300_l15 | lgb_rank_n300_l31 | xgb_pairwise_n300_d6 |
|---|---|---|---|---|---|
| zero-shot seeds 3-5 | 0.956 | 0.947 | 0.951 | 0.953 | 0.950 |
| areas (mean of 9) | 0.884 | 0.870 | 0.892 | 0.897 | 0.881 |
