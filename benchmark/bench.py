"""OREVID-bench: an honest, reproducible benchmark for drug repurposing on PrimeKG.

Task: for a disease, rank all 7,957 PrimeKG drugs so that the drugs that treat it come first.

    python bench.py build    --kg kg.csv --nodes node.csv          # tasks/ (deterministic)
    python bench.py baseline --tasks tasks                          # predictions/popularity/
    python bench.py evaluate --tasks tasks --pred predictions/<model>

See README.md. Only numpy, pandas, scipy and scikit-learn are needed.
"""
import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.metrics import average_precision_score, roc_auc_score

HERE = Path(__file__).resolve().parent
KG_MD5 = "aac8191d4fbc5bf09cdf8c3c78b4e75f"     # PrimeKG kg.csv, Harvard Dataverse file 7144484
NODES_MD5 = "4924de04fb3deefa1e0a8dada424538e"  # PrimeKG node.csv, Harvard Dataverse file 7144482
DD = ["contraindication", "indication", "off-label use"]
AREAS = ["cell_proliferation", "mental_health", "cardiovascular", "anemia", "adrenal_gland",
         "autoimmune", "metabolic_disorder", "diabetes", "neurodigenerative"]
SEEDS = [1, 2, 3, 4, 5]


# --- knowledge graph ----------------------------------------------------------------------

def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        while b := f.read(1 << 20):
            h.update(b)
    return h.hexdigest()


def load_kg(kg_csv):
    """One row per edge, exactly as TxGNN's kg_directed.csv: for relations between nodes of the same
    type, one direction per pair (first occurrence); otherwise the direction of the relation's first
    row. Original row order is kept (the splits depend on it)."""
    kg = pd.read_csv(kg_csv, usecols=["relation", "display_relation", "x_index", "x_id", "x_type",
                                      "y_index", "y_id", "y_type"], dtype={"x_id": str, "y_id": str},
                     low_memory=False)
    keep = []
    for rel in kg.relation.unique():
        d = kg[kg.relation == rel]
        if "_" in rel and rel.split("_")[0] == rel.split("_")[1]:
            pair = np.minimum(d.x_index, d.y_index).astype(np.int64) * 10**6 + np.maximum(d.x_index, d.y_index)
            keep.append(d.index[~pair.duplicated().to_numpy()])
        else:
            keep.append(d.index[(d.x_type == d.x_type.iloc[0]).to_numpy()])
    return kg.loc[np.sort(np.concatenate(keep))].reset_index(drop=True)


def _sample(n, frac, seed):
    """Positions chosen by pandas.DataFrame.sample(frac=frac, random_state=seed)."""
    return np.random.RandomState(seed).choice(n, size=int(round(frac * n)), replace=False)


def split_zero_shot(kg, seed):
    """TxGNN's «complex_disease» split, reproduced exactly: 0 train, 1 valid, 2 test (per row)."""
    lab = np.zeros(len(kg), np.int8)
    is_dd = kg.relation.isin(DD).to_numpy()
    rel = kg.relation.to_numpy()
    diseases = pd.unique(kg.y_index.to_numpy()[is_dd])
    np.random.seed(seed)
    np.random.shuffle(diseases)
    n = len(diseases)
    cut = [int(0.83125 * n), int((0.83125 + 0.11875) * n)]
    y = kg.y_index.to_numpy()
    lab[is_dd & np.isin(y, diseases[cut[0]:cut[1]])] = 1
    lab[is_dd & np.isin(y, diseases[cut[1]:])] = 2
    for r in pd.unique(rel[~is_dd]):
        pos = np.flatnonzero(rel == r)
        p = _sample(len(pos), 0.05, seed)
        lab[pos[p]] = 2
        rest = np.setdiff1d(np.arange(len(pos)), p)
        lab[pos[rest[_sample(len(rest), 0.11875 / 0.95, 1)]]] = 1
    return lab


def _norm_id(x):
    """TxGNN's convert2str for node ids («8840» -> «8840.0»)."""
    x = str(x)
    try:
        return x if "_" in x else str(float(x))
    except ValueError:
        return x


def area_diseases(area, nodes_csv):
    """Disease nodes of a TxGNN disease area (lists in disease_areas/, matched by MONDO id)."""
    nodes = pd.read_csv(nodes_csv, sep="\t", dtype=str)
    nodes = nodes[nodes.node_type == "disease"]
    id2idx = {}
    for i, idx in zip(nodes.node_id, nodes.node_index.astype(int)):
        id2idx[_norm_id(i)] = idx
        if "_" in i:
            for part in i.split("_"):
                id2idx.setdefault(_norm_id(part), idx)
    d = pd.read_csv(HERE / "disease_areas" / f"{area}.csv", dtype=str)
    return np.unique([id2idx[_norm_id(i)] for i in d.node_id if _norm_id(i) in id2idx])


def split_area(kg, area, nodes_csv, frac=0.05):
    """Disease-area split: all drug-disease edges of the area's diseases are hidden, plus random edges
    from their 2-hop neighbourhood up to 5 % of the graph (fixed seed; TxGNN's own version has no
    seed), and a 12.5 % validation sample of each relation (as TxGNN)."""
    nodes = area_diseases(area, nodes_csv)
    is_dd = kg.relation.isin(DD).to_numpy()
    x, y, rel = kg.x_index.to_numpy(), kg.y_index.to_numpy(), kg.relation.to_numpy()
    test = is_dd & (np.isin(x, nodes) | np.isin(y, nodes))
    N = int(max(x.max(), y.max())) + 1
    A = sp.csr_matrix((np.ones(len(kg)), (x, y)), shape=(N, N))
    A = (A + A.T).tocsr()
    inside = np.zeros(N, bool)
    inside[nodes] = True
    for _ in range(2):
        inside |= np.asarray(A[inside].sum(axis=0)).ravel() > 0
    cand = np.flatnonzero(inside[x] & inside[y] & ~test)
    n_rand = int(round(2 * len(kg) * frac)) - 2 * int(test.sum())
    pick = np.random.default_rng(0).choice(cand, min(len(cand), max(n_rand, 0) // 2), replace=False)
    test[pick] = True
    lab = np.where(test, 2, 0).astype(np.int8)
    for r in pd.unique(rel):
        pos = np.flatnonzero((rel == r) & (lab == 0))
        lab[pos[_sample(len(pos), 0.125, 1)]] = 1
    return lab, nodes


def _sym(kg, mask, rel, N):
    s = kg[mask & (kg.relation == rel).to_numpy()]
    x, y = s.x_index.to_numpy(), s.y_index.to_numpy()
    m = sp.coo_matrix((np.ones(2 * len(s)), (np.r_[x, y], np.r_[y, x])), shape=(N, N)).tocsr()
    m.data[:] = 1
    return m


def build_task(kg, lab, name, area_nodes=None):
    """Everything needed to evaluate a split, without the graph."""
    N = int(max(kg.x_index.max(), kg.y_index.max())) + 1
    types = np.empty(N, dtype=object)
    types[kg.x_index.to_numpy()] = kg.x_type.to_numpy()
    types[kg.y_index.to_numpy()] = kg.y_type.to_numpy()
    drugs = np.flatnonzero(types == "drug")
    train = lab == 0
    ind = kg.relation.to_numpy() == "indication"
    x, y = kg.x_index.to_numpy(), kg.y_index.to_numpy()
    test_ind = ind & (lab == 2)
    if area_nodes is not None:
        test_ind &= np.isin(y, area_nodes)
    diseases = pd.unique(y[test_ind])
    I = _sym(kg, train, "indication", N)
    popularity = np.diff(I.indptr)[drugs]
    O = _sym(kg, train, "disease_disease", N)
    O2 = (O + O @ O).tocsr()
    rows = []
    for d in diseases:
        pos = x[test_ind & (y == d)]
        nb = O2.indices[O2.indptr[d]:O2.indptr[d + 1]]
        nb = nb[nb != d]
        s = set(pos)
        best = max((len(s.intersection(I.indices[I.indptr[v]:I.indptr[v + 1]])) / len(s) for v in nb), default=0.0)
        rows.append({"disease": int(d), "positives": [int(p) for p in sorted(set(pos))], "near_duplicate_overlap": best})
    return {"name": name, "drugs": drugs.tolist(), "train_popularity": popularity.tolist(),
            "diseases": rows, "n_test_diseases": len(rows)}


def build(kg_csv, nodes_csv, out, ot_external=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for f, ref in [(kg_csv, KG_MD5), (nodes_csv, NODES_MD5)]:
        if md5(f) != ref:
            sys.exit(f"{f}: MD5 does not match the PrimeKG release this benchmark was built on ({ref}).")
    kg = load_kg(kg_csv)
    print(f"directed KG: {len(kg):,} edges", flush=True)
    for seed in SEEDS:
        lab = split_zero_shot(kg, seed)
        np.savez_compressed(out / f"zero_shot_{seed}.split.npz", valid=np.flatnonzero(lab == 1), test=np.flatnonzero(lab == 2))
        (out / f"zero_shot_{seed}.json").write_text(json.dumps(build_task(kg, lab, f"zero_shot_{seed}")))
        print(f"  zero_shot_{seed}", flush=True)
    for area in AREAS:
        lab, nodes = split_area(kg, area, nodes_csv)
        np.savez_compressed(out / f"area_{area}.split.npz", valid=np.flatnonzero(lab == 1), test=np.flatnonzero(lab == 2))
        (out / f"area_{area}.json").write_text(json.dumps(build_task(kg, lab, f"area_{area}", nodes)))
        print(f"  area_{area}", flush=True)


def train_edges(kg, task_dir, name):
    """Boolean mask of the training edges of a task (for model builders)."""
    z = np.load(Path(task_dir) / f"{name}.split.npz")
    m = np.ones(len(kg), bool)
    m[z["valid"]] = False
    m[z["test"]] = False
    return m


# --- metrics ------------------------------------------------------------------------------

def _jitter(s, rng):
    return s + rng.random(len(s)) * 1e-9 * (np.abs(s).max() + 1e-9)


def txgnn_metrics(S, positives, pool_neg):
    """TxGNN's metrics (reproduce/gather result-txgnn.ipynb): AUPRC/AUROC of positives vs as many
    random negatives (seed 42, with replacement) among drugs with a training indication; Recall@100,
    Precision@K and average normalised rank over all drugs. Ties are broken at random."""
    rng = np.random.default_rng(0)
    pp, pn, rec, prk, rank = [], [], [], [], []
    n = S.shape[0]
    for c, pos in enumerate(positives):
        s = _jitter(S[:, c], rng)
        np.random.seed(42)
        opts = np.intersect1d(np.setdiff1d(np.arange(n), pos), np.flatnonzero(pool_neg))
        neg = np.random.choice(opts, len(pos))
        pp += s[pos].tolist()
        pn += s[neg].tolist()
        order = np.argsort(s)[::-1]
        rec.append(len(np.intersect1d(order[:100], pos)) / len(pos))
        prk.append(len(np.intersect1d(order[:len(pos)], pos)) / len(pos))
        rank.append(np.mean(np.argsort(np.argsort(s))[pos]) / n)
    y = [1] * len(pp) + [0] * len(pn)
    return {"AUPRC": average_precision_score(y, pp + pn), "AUROC": roc_auc_score(y, pp + pn),
            "Recall@100": float(np.mean(rec)), "Precision@K": float(np.mean(prk)), "Avg Rank": float(np.mean(rank))}


def ranking(s, pos, cand, rng, k=100):
    """Recall@k, Precision@K and normalised rank of the positives among the candidates."""
    pos = np.intersect1d(pos, np.flatnonzero(cand))
    if not len(pos):
        return None
    idx = np.flatnonzero(cand)
    v = _jitter(s[idx], rng)
    r = np.empty(len(idx))
    r[np.argsort(-v)] = np.arange(len(idx))
    pr = r[np.searchsorted(idx, pos)]
    return {"Recall@100": float(np.mean(pr < k)), "Precision@K": float(np.mean(pr < len(pos))),
            "Avg Rank": float(1 - np.mean(pr) / len(idx))}


def load_predictions(path, task):
    """Scores matrix (drugs × test diseases). Missing pairs get the lowest score."""
    f = Path(path)
    d = pd.read_parquet(f) if f.suffix == ".parquet" else pd.read_csv(f)
    drugs = np.asarray(task["drugs"])
    dis = [r["disease"] for r in task["diseases"]]
    pd_ = {v: i for i, v in enumerate(drugs)}
    pdis = {v: i for i, v in enumerate(dis)}
    S = np.full((len(drugs), len(dis)), np.nan)
    d = d[d.disease.isin(pdis) & d.drug.isin(pd_)]
    S[d.drug.map(pd_).to_numpy(), d.disease.map(pdis).to_numpy()] = d.score.to_numpy()
    lo = np.nanmin(S) if np.isfinite(S).any() else 0
    return np.nan_to_num(S, nan=lo - 1), int(np.isnan(S).all(axis=0).sum())


def evaluate_external(task, S):
    """External task: positives never seen by any model; known pairs are not candidates."""
    drugs = np.asarray(task["drugs"])
    where = {v: i for i, v in enumerate(drugs)}
    rng = np.random.default_rng(0)
    out = {}
    for sub in ["positives", "recent", "non_obvious"]:
        v = []
        for c, r in enumerate(task["diseases"]):
            pos = np.array([where[p] for p in r[sub]], int)
            if not len(pos):
                continue
            cand = np.ones(len(drugs), bool)
            cand[[where[p] for p in r["exclude"] if p in where]] = False
            m = ranking(S[:, c], pos, cand, rng)
            if m:
                v.append(m)
        if v:
            out[{"positives": "all"}.get(sub, sub)] = {m: float(np.mean([x[m] for x in v])) for m in v[0]} | {"n_diseases": len(v)}
    out["n_test_diseases"] = len(task["diseases"])
    return out


def evaluate_task(task, S):
    if task.get("type") == "external":
        return evaluate_external(task, S)
    drugs = np.asarray(task["drugs"])
    where = {v: i for i, v in enumerate(drugs)}
    positives = [np.array([where[p] for p in r["positives"]]) for r in task["diseases"]]
    pop = np.asarray(task["train_popularity"])
    out = {"all": txgnn_metrics(S, positives, pop > 0)}
    hard = [i for i, r in enumerate(task["diseases"]) if r["near_duplicate_overlap"] < 0.5]
    if hard:
        out["hard"] = txgnn_metrics(S[:, hard], [positives[i] for i in hard], pop > 0)
    # Novelty: A approved drugs only; B without the top-10 % «blockbuster» drugs; C drugs with no
    # training indication (the other positives are removed from the list).
    rng = np.random.default_rng(0)
    approved = pop > 0
    blockbuster = pop >= np.quantile(pop[approved], 0.9)
    nov = {"A_approved_only": [], "B_no_blockbusters": [], "C_no_prior_use": []}
    for c, pos in enumerate(positives):
        others = np.zeros(len(drugs), bool)
        others[pos[pop[pos] > 0]] = True
        for key, (p, cand) in {"A_approved_only": (pos, approved), "B_no_blockbusters": (pos, ~blockbuster),
                               "C_no_prior_use": (pos[pop[pos] == 0], ~others)}.items():
            m = ranking(S[:, c], p, cand, rng)
            if m:
                nov[key].append(m)
    for key, v in nov.items():
        if v:
            out[key] = {m: float(np.mean([x[m] for x in v])) for m in v[0]} | {"n_diseases": len(v)}
    out["n_test_diseases"] = len(positives)
    out["n_hard"] = len(hard)
    return out


def validate(task_dir, pred_dir):
    """Checks a prediction folder before evaluation: one file per task, columns, ids, finite scores, no duplicates,
    every drug scored for every test disease. Returns the list of problems (empty = valid)."""
    task_dir, pred_dir = Path(task_dir), Path(pred_dir)
    problems = []
    for tf in sorted(task_dir.glob("*.json")):
        task = json.loads(tf.read_text())
        name = task["name"]
        pf = next((p for p in [pred_dir / f"{name}.parquet", pred_dir / f"{name}.csv"] if p.exists()), None)
        if pf is None:
            problems.append(f"{name}: no file")
            continue
        d = pd.read_parquet(pf) if pf.suffix == ".parquet" else pd.read_csv(pf)
        if missing := {"disease", "drug", "score"} - set(d.columns):
            problems.append(f"{name}: missing columns {sorted(missing)}")
            continue
        score = pd.to_numeric(d.score, errors="coerce")
        if not np.isfinite(score).all():
            problems.append(f"{name}: {int((~np.isfinite(score)).sum())} scores not finite")
        if n := int(d.duplicated(["disease", "drug"]).sum()):
            problems.append(f"{name}: {n} duplicated (disease, drug) pairs")
        drugs, dis = set(task["drugs"]), {r["disease"] for r in task["diseases"]}
        if n := int((~d.drug.isin(drugs)).sum()):
            problems.append(f"{name}: {n} rows with a drug not in the task's drugs")
        if n := int((~d.disease.isin(dis)).sum()):
            problems.append(f"{name}: {n} rows with a disease not among the test diseases")
        ok = d[d.drug.isin(drugs) & d.disease.isin(dis)].drop_duplicates(["disease", "drug"])
        if (n := len(drugs) * len(dis) - len(ok)) > 0:
            problems.append(f"{name}: {n:,} of {len(drugs) * len(dis):,} (disease, drug) pairs not scored "
                            f"({len(dis - set(ok.disease))} test diseases with no score at all)")
    print("\n".join(problems) if problems else f"valid: every task in {task_dir} has a complete prediction file")
    return problems


def evaluate(task_dir, pred_dir):
    task_dir, pred_dir = Path(task_dir), Path(pred_dir)
    res = {}
    for tf in sorted(task_dir.glob("*.json")):
        task = json.loads(tf.read_text())
        pf = next((p for p in [pred_dir / f"{task['name']}.parquet", pred_dir / f"{task['name']}.csv"] if p.exists()), None)
        if pf is None:
            continue
        S, missing = load_predictions(pf, task)
        res[task["name"]] = evaluate_task(task, S) | {"diseases_without_predictions": missing}
        print(f"  {task['name']}: Recall@100 {res[task['name']]['all']['Recall@100']:.3f}", flush=True)
    (pred_dir / "results.json").write_text(json.dumps(res, indent=1))
    return res


def summary(res):
    """Mean over the 5 zero-shot seeds and over the 9 areas."""
    rows = []
    for group, keys in [("zero-shot (5 seeds)", [k for k in res if k.startswith("zero_shot")]),
                        ("disease areas (9)", [k for k in res if k.startswith("area_")]),
                        ("external (Open Targets 26.06)", [k for k in res if k.startswith("external")])]:
        if not keys:
            continue
        for sub in ["all", "hard", "A_approved_only", "B_no_blockbusters", "C_no_prior_use", "recent", "non_obvious"]:
            v = [res[k][sub] for k in keys if sub in res[k]]
            if v:
                rows.append({"group": group, "subset": sub,
                             **{m: float(np.mean([x[m] for x in v])) for m in v[0] if m != "n_diseases"}})
    return pd.DataFrame(rows)


def baseline_popularity(task_dir, out):
    """Popularity: rank drugs by their number of training indications (ignores the disease)."""
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    for tf in sorted(Path(task_dir).glob("*.json")):
        t = json.loads(tf.read_text())
        drugs, pop = np.asarray(t["drugs"]), np.asarray(t["train_popularity"], float)
        dis = [r["disease"] for r in t["diseases"]]
        pd.DataFrame({"disease": np.repeat(dis, len(drugs)), "drug": np.tile(drugs, len(dis)),
                      "score": np.tile(pop, len(dis))}).to_parquet(out / f"{t['name']}.parquet", index=False)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    s = p.add_subparsers(dest="cmd", required=True)
    b = s.add_parser("build"); b.add_argument("--kg", required=True); b.add_argument("--nodes", required=True)
    b.add_argument("--out", default=str(HERE / "tasks"))
    b = s.add_parser("baseline"); b.add_argument("--tasks", default=str(HERE / "tasks"))
    b.add_argument("--out", default=str(HERE / "predictions" / "popularity"))
    b = s.add_parser("evaluate"); b.add_argument("--tasks", default=str(HERE / "tasks")); b.add_argument("--pred", required=True)
    b = s.add_parser("validate"); b.add_argument("--tasks", default=str(HERE / "tasks")); b.add_argument("--pred", required=True)
    a = p.parse_args()
    if a.cmd == "build":
        build(a.kg, a.nodes, a.out)
    elif a.cmd == "baseline":
        baseline_popularity(a.tasks, a.out)
    elif a.cmd == "validate":
        sys.exit(1 if validate(a.tasks, a.pred) else 0)
    else:
        res = evaluate(a.tasks, a.pred)
        print(summary(res).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
