"""Does XGBoost rank better than LightGBM? (docs/XGBOOST.md)

Same training rows, features, masked training rounds and metrics as `evaluate()` in
evaluation/txgnn.py and `score_ours()` in evaluation/external_exam.py; only the learner changes.
The protocol was written in docs/XGBOOST.md before running:

    tune     zero-shot seeds 1-2: every candidate configuration (LightGBM and XGBoost)
    confirm  zero-shot seeds 3-5, the 9 disease areas and the external exam: the production
             LightGBM, the best LightGBM and the best XGBoost of `tune`
    report   writes docs/XGBOOST.md tables from data/eval/xgb/results.csv

    py -3 -m evaluation.xgboost_compare tune
    py -3 -m evaluation.xgboost_compare confirm
    py -3 -m evaluation.xgboost_compare report

Training rows are cached in data/eval/xgb/ (several GB; delete when done). XGBoost runs on the GPU
(device="cuda") when there is one.
"""
import json
import sys
import time

import numpy as np
import pandas as pd

from evaluation.txgnn import (AREAS, CACHE, FEATURE_SETS, INDICATION, PRODUCTION, TxGraph, area_diseases,
                              drug_degrees, load_kg, masked_rows, metrics, near_duplicates, node_types, score_paths,
                              split_edges)

OUT = CACHE / "xgb"
RESULTS = OUT / "results.csv"
ROUNDS, N_TRAIN = 5, 1500

# --- candidate configurations (fixed before running, docs/XGBOOST.md) --------------------
PRODUCTION_MODEL = "lgb_rank_n300_l31"
CONFIGS = {
    # LightGBM: the production model and the same small grid around it as for XGBoost
    "lgb_rank_n300_l31": ("lgb", dict(objective="lambdarank", n_estimators=300, num_leaves=31)),
    "lgb_rank_n600_l31": ("lgb", dict(objective="lambdarank", n_estimators=600, num_leaves=31)),
    "lgb_rank_n300_l15": ("lgb", dict(objective="lambdarank", n_estimators=300, num_leaves=15)),
    "lgb_rank_n300_l63": ("lgb", dict(objective="lambdarank", n_estimators=300, num_leaves=63)),
    "lgb_rank_n600_l63": ("lgb", dict(objective="lambdarank", n_estimators=600, num_leaves=63)),
    "lgb_binary_n300_l31": ("lgb", dict(objective="binary", n_estimators=300, num_leaves=31)),
    # XGBoost: ranking objectives (NDCG, pairwise) and binary classification
    "xgb_ndcg_n300_d4": ("xgb", dict(objective="rank:ndcg", n_estimators=300, max_depth=4)),
    "xgb_ndcg_n300_d6": ("xgb", dict(objective="rank:ndcg", n_estimators=300, max_depth=6)),
    "xgb_ndcg_n300_d8": ("xgb", dict(objective="rank:ndcg", n_estimators=300, max_depth=8)),
    "xgb_ndcg_n600_d4": ("xgb", dict(objective="rank:ndcg", n_estimators=600, max_depth=4)),
    "xgb_ndcg_n600_d6": ("xgb", dict(objective="rank:ndcg", n_estimators=600, max_depth=6)),
    "xgb_ndcg_n600_d8": ("xgb", dict(objective="rank:ndcg", n_estimators=600, max_depth=8)),
    "xgb_ndcg_n600_d6_sub": ("xgb", dict(objective="rank:ndcg", n_estimators=600, max_depth=6,
                                         subsample=0.8, colsample_bytree=0.8)),
    "xgb_pairwise_n300_d6": ("xgb", dict(objective="rank:pairwise", n_estimators=300, max_depth=6)),
    "xgb_pairwise_n600_d8": ("xgb", dict(objective="rank:pairwise", n_estimators=600, max_depth=8)),
    "xgb_binary_n300_d6": ("xgb", dict(objective="binary:logistic", n_estimators=300, max_depth=6)),
    "xgb_binary_n600_d8": ("xgb", dict(objective="binary:logistic", n_estimators=600, max_depth=8)),
}
# Second round (added after the first tuning, before any confirmation): CatBoost ranking losses and
# rank-averaged ensembles. Same selection on seeds 1-2 (docs/XGBOOST.md).
CONFIGS.update({
    "cb_yetirank_n600_d6": ("cb", dict(loss_function="YetiRank", iterations=600, depth=6, learning_rate=0.05)),
    "cb_yetirank_n600_d6_lr10": ("cb", dict(loss_function="YetiRank", iterations=600, depth=6, learning_rate=0.1)),
    "cb_yetirank_n600_d8": ("cb", dict(loss_function="YetiRank", iterations=600, depth=8, learning_rate=0.05)),
    "cb_yetirank_n1200_d6": ("cb", dict(loss_function="YetiRank", iterations=1200, depth=6, learning_rate=0.05)),
    "cb_pairlogit_n600_d6": ("cb", dict(loss_function="PairLogit", iterations=600, depth=6, learning_rate=0.05)),
    "logit_nopop": ("logit", {}),  # the published model without popularity (evaluation/txgnn.py fit_model)
    "ens_lgb_xgb": ("ens", ["lgb_rank_n300_l31", "xgb_pairwise_n300_d6"]),
    "ens_lgb_cb": ("ens", ["lgb_rank_n300_l31", "cb_yetirank_n600_d6"]),
    "ens_lgb_xgb_cb": ("ens", ["lgb_rank_n300_l31", "xgb_pairwise_n300_d6", "cb_yetirank_n600_d6"]),
    "ens_lgb_logit": ("ens", ["lgb_rank_n300_l31", "logit_nopop"]),
})
FITTED = {}  # models already trained on the current split (ensembles reuse them)


def fit(name, X, y, groups, seed):
    """Trains configuration `name` and returns its scoring function (cached per split in FITTED)."""
    if name not in FITTED:
        FITTED[name] = _fit(name, X, y, groups, seed)
    return FITTED[name]


def _fit(name, X, y, groups, seed):
    lib, p = CONFIGS[name]
    if lib == "ens":  # mean of the per-disease rank percentiles of the members
        from scipy.stats import rankdata
        parts = [fit(m, X, y, groups, seed) for m in p]
        return lambda Z: np.mean([rankdata(f(Z)) / len(Z) for f in parts], axis=0)
    if lib == "logit":
        from evaluation.txgnn import fit_model
        nc = len(FEATURE_SETS[PRODUCTION])
        model, _ = fit_model(X[:, :nc], y)
        return lambda Z: model(Z[:, :nc])
    if lib == "cb":
        import catboost
        clf = catboost.CatBoostRanker(**p, metric_period=p["iterations"], allow_writing_files=False, l2_leaf_reg=1.0, random_seed=seed, task_type="GPU" if DEVICE == "cuda" else "CPU",
                                      verbose=0)
        clf.fit(X, y.astype(int), group_id=np.repeat(np.arange(len(groups)), groups))
        return clf.predict
    p = dict(p)
    if lib == "lgb":
        import lightgbm as lgb
        common = dict(learning_rate=0.05, min_child_samples=50, reg_lambda=1.0, random_state=seed, verbose=-1)
        if p["objective"] == "lambdarank":
            clf = lgb.LGBMRanker(**p, **common)
            clf.fit(X, y.astype(int), group=groups)
            return clf.predict
        clf = lgb.LGBMClassifier(**p, **common)
        clf.fit(X, y.astype(int))
        return lambda Z: clf.predict_proba(Z)[:, 1]
    import xgboost as xgb
    common = dict(learning_rate=0.05, reg_lambda=1.0, tree_method="hist", device=DEVICE, random_state=seed)
    if p["objective"].startswith("rank:"):
        clf = xgb.XGBRanker(**p, **common)
        clf.fit(X, y.astype(int), qid=np.repeat(np.arange(len(groups)), groups))
        return clf.predict
    clf = xgb.XGBClassifier(**p, **common)
    clf.fit(X, y.astype(int))
    return lambda Z: clf.predict_proba(Z)[:, 1]


def _device():
    try:
        import xgboost as xgb
        xgb.XGBRegressor(n_estimators=1, device="cuda").fit(np.zeros((4, 1)), np.zeros(4))
        return "cuda"
    except Exception:
        return "cpu"


DEVICE = None


# --- training and test data, exactly as evaluate() builds them ---------------------------

def build(kg, types, split, seed):
    """Training rows (with masked rounds), test features and labels for one split; cached."""
    f = OUT / f"{split}_{seed}_{PRODUCTION}.npz"  # keyed by the feature set: a stale cache must not be reused
    if f.exists():
        z = np.load(f, allow_pickle=True)
        return {k: z[k] for k in z.files}
    t0 = time.time()
    paths = FEATURE_SETS[PRODUCTION]
    label = split_edges(kg, split, seed)
    g = TxGraph(kg, label == 0, types, extra=True)
    drug_pos = np.full(g.N, -1)
    drug_pos[g.drugs] = np.arange(len(g.drugs))
    ind = kg.relation.to_numpy() == "indication"
    x, y = kg.x_index.to_numpy(), kg.y_index.to_numpy()
    test_ind = ind & (label == 2)
    if split in AREAS:
        test_ind &= np.isin(y, area_diseases(split))
    test_diseases = pd.unique(y[test_ind])
    test_pos = [drug_pos[x[test_ind & (y == d)]] for d in test_diseases]
    degree = np.diff(g.M(INDICATION).indptr)[g.drugs].astype(np.float64)
    known = ind & (label < 2) & ~np.isin(y, test_diseases)
    train_diseases = pd.unique(y[known])
    train_diseases = np.random.default_rng(seed).choice(train_diseases, min(N_TRAIN, len(train_diseases)), replace=False)
    train_pos = {d: drug_pos[x[known & (y == d)]] for d in train_diseases}
    is_train = ind & (label == 0)
    F_test = score_paths(g, test_diseases, paths)
    rows_X, rows_y = [], []
    rng = np.random.default_rng(seed)
    gf = drug_degrees(g)
    for i in range(0, len(train_diseases), 256):
        block = train_diseases[i:i + 256]
        F_train = score_paths(g, block, paths)
        for c, d in enumerate(block):
            lab = np.zeros(len(g.drugs), bool)
            lab[train_pos[d]] = True
            pop = degree.copy()
            pop[drug_pos[x[is_train & (y == d)]]] -= 1
            sel = lab | (rng.random(len(lab)) < 0.1)
            rows_X.append(np.column_stack([F_train[sel, c, :], pop[sel], gf[sel]]))
            rows_y.append(lab[sel])
    X, Y = np.concatenate(rows_X), np.concatenate(rows_y)
    Xm, Ym, Gm = masked_rows(g, test_diseases, known, x, y, drug_pos, paths, ROUNDS, seed)
    nd = near_duplicates(g, test_diseases, test_pos, drug_pos)
    data = dict(X=np.concatenate([X, Xm]), Y=np.concatenate([Y, Ym]), groups=np.r_[[len(v) for v in rows_y], Gm],
                F_test=F_test.astype(np.float64), degree=degree, gf=gf, pool_neg=degree > 0,
                test_pos=np.array([np.asarray(p) for p in test_pos], dtype=object), hard=np.flatnonzero(nd < 0.5),
                diseases=test_diseases, drugs=g.drugs)
    if split == "complex_disease" and seed == 1:  # the external exam is trained on this split
        data.update(external_features(kg, g, paths))
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(f, **data)
    print(f"  data {split} seed {seed}: {len(data['Y']):,} rows, {len(test_diseases)} test diseases "
          f"({time.time() - t0:.0f} s)", flush=True)
    return data


def external_features(kg, g, paths):
    """Features of the external exam diseases (evaluation/external_exam.py) on the seed-1 graph."""
    z = np.load(CACHE / "external_scores.npz", allow_pickle=True)
    diseases = z["diseases"]
    assert np.array_equal(z["drugs"], g.drugs)
    F = np.concatenate([score_paths(g, diseases[i:i + 256], paths) for i in range(0, len(diseases), 256)], axis=1)
    return {"F_ext": F.astype(np.float32)}


def scores_for(predict, F, degree, gf):
    return np.column_stack([predict(np.column_stack([F[:, c, :], degree, gf])) for c in range(F.shape[1])])


def external_metrics(S):
    """Recall@100, Precision@K and mean rank on the external exam, as evaluation/external_exam.py."""
    from evaluation.novelty import measure
    z = np.load(CACHE / "external_scores.npz", allow_pickle=True)
    p = pd.read_parquet(CACHE / "external_positives.parquet")
    drugs = z["drugs"]
    pos_f = {int(f): j for j, f in enumerate(drugs)}
    rng = np.random.default_rng(0)
    rows = []
    for c, (d, known) in enumerate(zip(z["diseases"], z["known"])):
        cand = np.ones(len(drugs), bool)
        cand[[pos_f[int(f)] for f in known if int(f) in pos_f]] = False
        pp = p[p.disease == d]
        for name, s in {"all": pp, "recent": pp[pp["year"] >= 2018], "non_obvious": pp[~pp.obvious]}.items():
            pos = np.array([pos_f[int(f)] for f in s.drug if int(f) in pos_f], int)
            if len(pos):
                m = measure(S[:, c].astype(np.float64), pos, cand, rng=rng)
                if m:
                    rows.append({"subset": name, **m})
    d = pd.DataFrame(rows)
    return {sub: e[["Recall@100", "Precision@K", "Avg Rank"]].mean().to_dict() | {"n_diseases": len(e)}
            for sub, e in d.groupby("subset")}


# --- runs -------------------------------------------------------------------------------

def done():
    if not RESULTS.exists():
        return set()
    r = pd.read_csv(RESULTS)
    return set(zip(r.stage, r.task, r.model, r.subset))


def record(rows):
    df = pd.DataFrame(rows)
    df.to_csv(RESULTS, mode="a", header=not RESULTS.exists(), index=False)


def run(stage, tasks, models, kg, types):
    seen = done()
    for split, seed in tasks:
        task = f"zero_shot_{seed}" if split == "complex_disease" else f"area_{split}"
        todo = [m for m in models if (stage, task, m, "all") not in seen]
        external = stage == "confirm" and split == "complex_disease" and seed == 1
        todo_ext = [m for m in models if ("confirm", "external", m, "all") not in seen] if external else []
        if not todo and not todo_ext:
            continue
        d = build(kg, types, split, seed)
        FITTED.clear()
        for m in sorted(set(todo) | set(todo_ext)):
            t0 = time.time()
            predict = fit(m, d["X"], d["Y"], d["groups"], seed)
            t_fit = time.time() - t0
            rows = []
            if m in todo:
                S = scores_for(predict, d["F_test"], d["degree"], d["gf"])
                pos = list(d["test_pos"])
                for sub, cols in (("all", np.arange(S.shape[1])), ("hard", d["hard"])):
                    r = metrics(S[:, cols], [pos[i] for i in cols], d["pool_neg"])
                    rows.append({"stage": stage, "task": task, "model": m, "subset": sub, "n_diseases": len(cols),
                                 "fit_seconds": round(t_fit, 1), **{k: float(v) for k, v in r.items()}})
                if m == PRODUCTION_MODEL:  # must reproduce the published OREVID scores exactly
                    check_production(split, seed, S)
            if m in todo_ext:
                Se = scores_for(predict, d["F_ext"].astype(np.float64), d["degree"], d["gf"])
                if m == PRODUCTION_MODEL:
                    ref = np.load(CACHE / "external_scores.npz", allow_pickle=True)["Ours_ranking"]
                    print(f"  check against the published external scores: max |difference| = "
                          f"{float(np.abs(ref - Se.astype(np.float32)).max()):.2e}", flush=True)
                for sub, r in external_metrics(Se).items():
                    rows.append({"stage": stage, "task": "external", "model": m, "subset": sub,
                                 "n_diseases": r.pop("n_diseases"), "fit_seconds": round(t_fit, 1), **r})
            record(rows)
            a = next(r for r in rows if r["subset"] == "all")
            print(f"  {stage} {task} {m}: Recall@100 {a['Recall@100']:.3f} ({t_fit:.0f} s)", flush=True)


def check_production(split, seed, S):
    f = CACHE / "scores" / f"{split}_{seed}_{PRODUCTION}_r5_rank.npz"
    if f.exists():
        ref = np.load(f, allow_pickle=True)["LambdaRank_+_masked"]
        diff = float(np.abs(ref - S.astype(np.float32)).max())
        print(f"  check against the published scores ({split} {seed}): max |difference| = {diff:.2e}", flush=True)


def best(stage_rows, prefix):
    r = stage_rows[(stage_rows.subset == "all") & stage_rows.model.str.startswith(prefix)]
    m = r.groupby("model")[["Recall@100"]].mean()
    h = stage_rows[(stage_rows.subset == "hard") & stage_rows.model.str.startswith(prefix)].groupby("model")["Recall@100"].mean()
    m["hard"] = h
    return m.sort_values(["Recall@100", "hard"], ascending=False).index[0]


def main():
    global DEVICE
    sys.stdout.reconfigure(encoding="utf-8")
    stage = sys.argv[1] if len(sys.argv) > 1 else "tune"
    if stage == "report":
        return report()
    DEVICE = _device()
    print(f"XGBoost device: {DEVICE}", flush=True)
    kg = load_kg()
    types = node_types(kg)
    if stage == "tune":
        run("tune", [("complex_disease", 1), ("complex_disease", 2)], list(CONFIGS), kg, types)
    elif stage == "confirm":
        r = pd.read_csv(RESULTS)
        r = r[r.stage == "tune"]
        chosen = list(dict.fromkeys([PRODUCTION_MODEL] + [best(r, f) for f in ("lgb_", "xgb_", "cb_", "ens_")
                                                          if r.model.str.startswith(f).any()]))
        (OUT / "chosen.json").write_text(json.dumps(chosen))
        print("confirming:", chosen, flush=True)
        tasks = [("complex_disease", s) for s in (1, 3, 4, 5)] + [(a, 1) for a in AREAS]
        run("confirm", tasks, chosen, kg, types)


# --- report -----------------------------------------------------------------------------

def report():
    r = pd.read_csv(RESULTS)
    L = []
    tune = r[r.stage == "tune"]
    if len(tune):
        L += ["### Tuning: zero-shot seeds 1-2 (mean)\n",
              "| Model | Recall@100 | Recall@100 hard | Precision@K | AUPRC | Fit (s) |", "|---|---|---|---|---|---|"]
        a = tune[tune.subset == "all"].groupby("model")[["Recall@100", "Precision@K", "AUPRC", "fit_seconds"]].mean()
        a["hard"] = tune[tune.subset == "hard"].groupby("model")["Recall@100"].mean()
        for m, v in a.sort_values("Recall@100", ascending=False).iterrows():
            L.append(f"| {m} | {v['Recall@100']:.3f} | {v['hard']:.3f} | {v['Precision@K']:.3f} | {v['AUPRC']:.3f} | {v['fit_seconds']:.0f} |")
        L.append("")
    conf = r[r.stage == "confirm"]
    if len(conf):
        models = list(dict.fromkeys(conf.model))
        L += ["### Confirmation: Recall@100 per task\n",
              "| Task | Subset | " + " | ".join(models) + " |", "|---|---|" + "---|" * len(models)]
        for (task, sub), e in conf.groupby(["task", "subset"], sort=False):
            v = e.set_index("model")["Recall@100"]
            L.append(f"| {task} | {sub} | " + " | ".join(f"{v.get(m, np.nan):.3f}" for m in models) + " |")
        L += ["", "### Confirmation: summary (mean Recall@100)\n",
              "| Group | " + " | ".join(models) + " |", "|---|" + "---|" * len(models)]
        groups = {"zero-shot seeds 3-5": conf.task.isin([f"zero_shot_{s}" for s in (3, 4, 5)]),
                  "areas (mean of 9)": conf.task.str.startswith("area_"),
                  "external": conf.task == "external"}
        for gname, sel in groups.items():
            for sub in ("all", "hard", "non_obvious", "recent"):
                e = conf[sel & (conf.subset == sub)]
                if len(e):
                    v = e.groupby("model")["Recall@100"].mean()
                    L.append(f"| {gname}, {sub} | " + " | ".join(f"{v.get(m, np.nan):.3f}" for m in models) + " |")
        for metric in ("Precision@K", "AUPRC"):
            L += ["", f"### Confirmation: mean {metric}\n", "| Group | " + " | ".join(models) + " |", "|---|" + "---|" * len(models)]
            for gname, sel in groups.items():
                e = conf[sel & (conf.subset == "all")]
                if len(e) and metric in e and e[metric].notna().any():
                    v = e.groupby("model")[metric].mean()
                    L.append(f"| {gname} | " + " | ".join(f"{v.get(m, np.nan):.3f}" for m in models) + " |")
    text = "\n".join(L) + "\n"
    (OUT / "report.md").write_text(text, "utf-8")
    print(text)


if __name__ == "__main__":
    main()
