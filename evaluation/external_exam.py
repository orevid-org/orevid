"""External exam: approved indications according to Open Targets 26.06 that are not in PrimeKG.

Neither model has seen them: neither TxGNN (trained on PrimeKG) nor ours (which does not use
any drug-disease relation from Open Targets). Both use the same training graph: the zero-shot
split of seed 1, where we already have TxGNN trained.

- Positives: (drug, disease) pairs «approved» in Open Targets without any relation in PrimeKG.
- Candidates: all PrimeKG drugs except those PrimeKG already relates to the disease.
- Subsets: all; recent (first trial ≥ 2018); non-obvious (the drug is not indicated in
  PrimeKG for any disease ≤ 2 steps away in the ontology).

    py -3 -m evaluation.external_exam prepare     # list of diseases for TxGNN
    (in WSL) txgnn_wsl.py --diseases data/eval/external_diseases.csv ...
    py -3 -m evaluation.external_exam evaluate
"""
import sys
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp

from evaluation.txgnn import (PRODUCTION, ROOT, CACHE, RAW, FEATURE_SETS, INDICATION, DD, TxGraph, _id, fit_gbm, fit_model,
                              fit_ranker, load_kg, masked_rows, drug_degrees, score_paths, split_edges, node_types)
from evaluation.novelty import measure

SPLIT, SEED = "complex_disease", 1
N_PKG = 129375


def new_positives(kg):
    """DataFrame (drug, disease, year) with the Open Targets approvals that are not in PrimeKG."""
    a = pd.read_parquet(ROOT / "data" / "prep" / "attributes.parquet",
                        columns=["relation", "display_relation", "x_index", "y_index", "year"])
    ap = a[(a.relation == "clinical") & (a.display_relation == "approved") & (a.x_index < N_PKG) & (a.y_index < N_PKG)]
    dd = kg[kg.relation.isin(DD)]
    known = set(zip(dd.x_index, dd.y_index))
    typ = node_types(kg)
    ap = ap[(typ[ap.x_index] == "drug") & (typ[ap.y_index] == "disease")]
    ap = ap[[(f, m) not in known for f, m in zip(ap.x_index, ap.y_index)]]
    return ap.rename(columns={"x_index": "drug", "y_index": "disease"}).groupby(
        ["drug", "disease"], as_index=False)["year"].min()


def prepare():
    kg = load_kg()
    p = new_positives(kg)
    nodes = pd.read_csv(RAW / "node.csv", sep="\t", dtype=str).set_index(pd.Index(range(N_PKG)))
    m = pd.DataFrame({"node_index": sorted(p.disease.unique())})
    m["node_id"] = nodes.node_id.to_numpy()[m.node_index]
    m.to_csv(CACHE / "external_diseases.csv", index=False)
    print(f"{len(p)} new indications in {len(m)} diseases -> {CACHE / 'external_diseases.csv'}")


def score_ours(kg, types, diseases, feature_set=None, rounds=5, n_train=1500, seed=SEED):
    feature_set = feature_set or PRODUCTION
    """Trains our model as evaluate() does (same split) and scores the given diseases."""
    t0 = time.time()
    paths = FEATURE_SETS[feature_set]
    label = split_edges(kg, SPLIT, seed)
    g = TxGraph(kg, label == 0, types, extra=feature_set.startswith("extra"))
    drug_pos = np.full(g.N, -1)
    drug_pos[g.drugs] = np.arange(len(g.drugs))
    ind = kg.relation.to_numpy() == "indication"
    x, y = kg.x_index.to_numpy(), kg.y_index.to_numpy()
    test_diseases = pd.unique(y[ind & (label == 2)])
    degree = np.diff(g.M(INDICATION).indptr)[g.drugs].astype(np.float64)
    known = ind & (label < 2) & ~np.isin(y, test_diseases)
    train_diseases = pd.unique(y[known])
    train_diseases = np.random.default_rng(seed).choice(train_diseases, min(n_train, len(train_diseases)), replace=False)
    is_train = ind & (label == 0)
    gf = drug_degrees(g)
    rng = np.random.default_rng(seed)
    rows_X, rows_y = [], []
    for i in range(0, len(train_diseases), 256):
        block = train_diseases[i:i + 256]
        F = score_paths(g, block, paths)
        for c, d in enumerate(block):
            lab = np.zeros(len(g.drugs), bool)
            lab[drug_pos[x[known & (y == d)]]] = True
            pop = degree.copy()
            pop[drug_pos[x[is_train & (y == d)]]] -= 1
            sel = lab | (rng.random(len(lab)) < 0.1)
            rows_X.append(np.column_stack([F[sel, c, :], pop[sel], gf[sel]]))
            rows_y.append(lab[sel])
    X, Y = np.concatenate(rows_X), np.concatenate(rows_y)
    nc = len(paths)
    graph_only, _ = fit_model(X[:, :nc], Y)
    Xm, Ym, Gm = masked_rows(g, test_diseases, known, x, y, drug_pos, paths, rounds, seed)
    gbm = fit_gbm(np.concatenate([X, Xm]), np.concatenate([Y, Ym]), seed)
    rank = fit_ranker(np.concatenate([X, Xm]), np.concatenate([Y, Ym]), np.r_[[len(v) for v in rows_y], Gm], seed)
    res = {"Popularity": np.repeat(degree[:, None], len(diseases), axis=1),
           "Ours without popularity": np.zeros((len(g.drugs), len(diseases)), np.float32),
           "Ours": np.zeros((len(g.drugs), len(diseases)), np.float32),
           "Ours (ranking)": np.zeros((len(g.drugs), len(diseases)), np.float32)}
    for i in range(0, len(diseases), 256):
        block = diseases[i:i + 256]
        F = score_paths(g, block, paths)
        for c in range(len(block)):
            Z = np.column_stack([F[:, c, :], degree, gf])
            res["Ours without popularity"][:, i + c] = graph_only(Z[:, :nc])
            res["Ours"][:, i + c] = gbm(Z)
            res["Ours (ranking)"][:, i + c] = rank(Z)
    print(f"our model: {len(diseases)} diseases scored ({time.time() - t0:.0f} s)")
    return g, res


def txgnn(diseases, drugs, seed=SEED):
    f = CACHE / "txgnn_predictions" / f"TxGNN_{seed}_{SPLIT}_external.parquet"
    if not f.exists():
        print("(there are no TxGNN predictions for the external exam)")
        return None
    d = pd.read_parquet(f)
    d = d[d.relation == "indication"]
    nodes = pd.read_csv(RAW / "node.csv", sep="\t", dtype=str)
    idx = {(t, _id(i)): int(n) for n, i, t in zip(nodes.node_index, nodes.node_id, nodes.node_type)}
    mi = np.array([idx.get(("disease", _id(i)), -1) for i in d.disease_id])
    fi = np.array([idx.get(("drug", _id(i)), -1) for i in d.drug_id])
    pos_m = {m: j for j, m in enumerate(diseases)}
    pos_f = np.full(N_PKG, -1)
    pos_f[drugs] = np.arange(len(drugs))
    res = np.full((len(drugs), len(diseases)), np.nan, np.float32)
    ok = np.isin(mi, list(pos_m)) & (fi >= 0)
    res[pos_f[fi[ok]], [pos_m[m] for m in mi[ok]]] = d.score.to_numpy()[ok]
    missing = np.isnan(res).all(axis=0).sum()
    if missing:
        print(f"⚠ TxGNN has no predictions for {missing} diseases")
    return res


def evaluate():
    sys.stdout.reconfigure(encoding="utf-8")
    kg = load_kg()
    types = node_types(kg)
    p = new_positives(kg)
    diseases = np.array(sorted(p.disease.unique()))
    args = sys.argv[2:]
    seed = int(args[args.index("--seed") + 1]) if "--seed" in args else SEED
    args = [x for i, x in enumerate(args) if x != "--seed" and (i == 0 or args[i - 1] != "--seed")]
    feature_set = args[0] if args else PRODUCTION
    g, scores = score_ours(kg, types, diseases, feature_set=feature_set, seed=seed)
    tx = txgnn(diseases, g.drugs, seed)
    if tx is not None:
        scores["TxGNN"] = tx
    pos_f = np.full(N_PKG, -1)
    pos_f[g.drugs] = np.arange(len(g.drugs))
    # Known in PrimeKG (the whole split, also validation and test): they are not candidates.
    dd = kg[kg.relation.isin(DD)]
    known = dd.groupby("y_index").x_index.apply(lambda s: pos_f[s.to_numpy()]).to_dict()
    # «Obvious»: the drug is indicated (training graph) for a disease ≤ 2 steps away in the ontology.
    O = g.M("disease_disease")
    O2 = ((O + O @ O) > 0).astype(np.float64).tocsr()
    I = g.M(INDICATION)
    rows = []
    rng = np.random.default_rng(0)
    obvious = []
    for c, d in enumerate(diseases):
        pp = p[p.disease == d]
        neighbors = O2.indices[O2.indptr[d]:O2.indptr[d + 1]]
        neighbor_drugs = set(I[neighbors].indices) if len(neighbors) else set()
        obvious += [(f, d, f in neighbor_drugs) for f in pp.drug]
        cand = np.ones(len(g.drugs), bool)
        cand[known.get(d, [])] = False
        subsets = {
            "all": pp,
            "recent (≥ 2018)": pp[pp["year"] >= 2018],
            "non-obvious": pp[[f not in neighbor_drugs for f in pp.drug]],
        }
        for name, s in subsets.items():
            pos = pos_f[s.drug.to_numpy()]
            if not len(pos):
                continue
            for model, S in scores.items():
                col = S[:, c]
                if np.isnan(col).all():
                    continue
                m = measure(np.nan_to_num(col, nan=np.nanmin(col) - 1), pos, cand, rng=rng)
                if m:
                    rows.append({"subset": name, "model": model, **m})
    # For OREVID-bench (evaluation/export_bench.py): scores, positives and known pairs.
    if feature_set != PRODUCTION:  # ablations: only show, do not save
        d = pd.DataFrame(rows)
        print(d[d.subset == "all"].groupby("model")[["Recall@100", "Precision@K", "Avg Rank"]].mean().round(3))
        return
    ob = pd.DataFrame(obvious, columns=["drug", "disease", "obvious"])
    suffix = "" if seed == SEED else f"_s{seed}"  # other training seeds (REVIEW_RESPONSE.md T2): separate files
    p.merge(ob, on=["drug", "disease"]).to_parquet(CACHE / f"external_positives{suffix}.parquet", index=False)
    np.savez_compressed(CACHE / f"external_scores{suffix}.npz", diseases=diseases, drugs=g.drugs,
                        degree=np.diff(g.M(INDICATION).indptr)[g.drugs],
                        known=np.array([g.drugs[known.get(d, np.array([], int))] for d in diseases], dtype=object),
                        **{k.replace(" ", "_").replace("(", "").replace(")", ""): v.astype(np.float32) for k, v in scores.items()})
    d = pd.DataFrame(rows)
    if seed != SEED:
        print(d.groupby(["subset", "model"])["Recall@100"].mean().round(3).to_string())
        return
    order = [m for m in ["Popularity", "TxGNN", "Ours without popularity", "Ours", "Ours (ranking)"]
             if m in d.model.unique()]
    L = ["# External exam: Open Targets indications that are not in PrimeKG\n",
         "Generated by `evaluation/external_exam.py`. Same training graph for both models (zero-shot, seed 1). "
         "Candidates: all PrimeKG drugs except those PrimeKG already relates to the disease. "
         "Mean per disease.\n",
         f"Positives: {len(p)} approved indications in {len(diseases)} diseases "
         f"({int((p['year'] >= 2018).sum())} with first trial ≥ 2018).\n",
         "| Subset | Metric | " + " | ".join(order) + " |", "|---|---|" + "---|" * len(order)]
    for sub in ["all", "recent (≥ 2018)", "non-obvious"]:
        e = d[d.subset == sub]
        if e.empty:
            continue
        n = e[e.model == "Ours"].shape[0]
        for met in ["Recall@100", "Precision@K", "Avg Rank"]:
            v = e.groupby("model")[met].mean().reindex(order)
            best = v.drop("Popularity").idxmax()
            L.append(f"| {sub} ({n} diseases) | {met} | " + " | ".join(
                (f"**{x:.2f}**" if k == best else f"{x:.2f}") for k, x in v.items()) + " |")
    (ROOT / "docs" / "EXTERNAL_EXAM.md").write_text("\n".join(L) + "\n", "utf-8")
    print("\n".join(L))


def seeds():
    """External exam over the five zero-shot training seeds (REVIEW_RESPONSE.md T2): mean ± SD of the
    per-seed Recall@100 (mean over diseases), from external_scores[_s<seed>].npz."""
    sys.stdout.reconfigure(encoding="utf-8")
    rows = []
    for seed in (1, 2, 3, 4, 5):
        suffix = "" if seed == SEED else f"_s{seed}"
        f = CACHE / f"external_scores{suffix}.npz"
        if not f.exists():
            print(f"(seed {seed} missing: py -3 -m evaluation.external_exam evaluate --seed {seed})")
            continue
        z = dict(np.load(f, allow_pickle=True))
        if "TxGNN" not in z:  # TxGNN's predictions may have been made after OREVID's run
            tx = txgnn(z["diseases"], z["drugs"], seed)
            if tx is not None:
                z["TxGNN"] = tx
        p = pd.read_parquet(CACHE / f"external_positives{suffix}.parquet")
        drugs = z["drugs"]
        pos_f = {int(x): j for j, x in enumerate(drugs)}
        rng = np.random.default_rng(0)
        models = {"Popularity": "Popularity", "TxGNN": "TxGNN", "Ours_without_popularity": "Logistic regression without popularity",
                  "Ours_ranking": "OREVID"}
        for c, (d, known) in enumerate(zip(z["diseases"], z["known"])):
            cand = np.ones(len(drugs), bool)
            cand[[pos_f[int(x)] for x in known if int(x) in pos_f]] = False
            e = p[p.disease == d]
            for sub, s_ in {"all": e, "recent": e[e["year"] >= 2018], "non-obvious": e[~e.obvious]}.items():
                pos = np.array([pos_f[int(x)] for x in s_.drug if int(x) in pos_f], int)
                if not len(pos):
                    continue
                for key, name in models.items():
                    if key not in z:
                        continue
                    col = z[key][:, c].astype(np.float64)
                    if np.isnan(col).all():
                        continue
                    m = measure(np.nan_to_num(col, nan=np.nanmin(col) - 1), pos, cand, rng=rng)
                    if m:
                        rows.append({"seed": seed, "subset": sub, "model": name, **m})
    d = pd.DataFrame(rows)
    per_seed = d.groupby(["subset", "model", "seed"])["Recall@100"].mean().unstack("seed")
    L = ["# External exam over five training seeds\n",
         "Generated by `py -3 -m evaluation.external_exam seeds` (REVIEW_RESPONSE.md T2). Recall@100 (mean over "
         "diseases) per training seed (the zero-shot splits 1–5), and mean ± SD over seeds. Seed 1 is the published "
         "exam. TxGNN: inference with the models trained on each zero-shot split.\n",
         "| Subset | Model | " + " | ".join(f"seed {s}" for s in per_seed.columns) + " | mean ± SD |",
         "|---|---|" + "---|" * (len(per_seed.columns) + 1)]
    for (sub, model), r in per_seed.iterrows():
        v = r.dropna()
        L.append(f"| {sub} | {model} | " + " | ".join("—" if x != x else f"{x:.3f}" for x in r)
                 + f" | {v.mean():.3f} ± {v.std(ddof=1) if len(v) > 1 else 0:.3f} |")
    (ROOT / "docs" / "EXTERNAL_SEEDS.md").write_text("\n".join(L) + "\n", "utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    {"prepare": prepare, "evaluate": evaluate, "seeds": seeds}[sys.argv[1]]()
