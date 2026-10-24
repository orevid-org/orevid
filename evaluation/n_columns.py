"""N5–N7: three sources as model columns (docs/IMPROVEMENT_PLAN.md, «Registered 24-10-2026: N5–N7»), each against the
production model with M2's rule (evaluation/m2_direction_feature.py procedure).

    py -3 -m evaluation.n_columns n5|n6|n7          # -> data/eval/<item>_feature/*.json, docs/<ITEM>_FEATURE.md
"""
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

import evaluation.proximity_feature as pf
import evaluation.txgnn as tx
import evaluation.xgboost_compare as xc

ROOT = Path(__file__).resolve().parent.parent


class HPOSimilar:
    """N5: log(1 + Σ over diseases where the drug is indicated (graph h) of the phenotypes they share with the disease)."""
    ph = None

    def __init__(self, g):
        import scipy.sparse as sp
        from evaluation.twins import twins
        from graph import Graph
        if HPOSimilar.ph is None:
            frames = [pd.read_parquet(ROOT / "data" / "ext" / f, columns=["disease", "phenotype"])
                      for f in ("hpo_current.parquet", "orphadata_phenotypes.parquet")]
            o = pd.concat(frames)
            P = g.M("disease_phenotype_positive").tocsr()
            x, y = P.nonzero()
            extra = sp.csr_matrix((np.ones(len(o), np.float32), (o.disease, o.phenotype)), shape=(g.N, g.N))
            HPOSimilar.ph = ((sp.csr_matrix((np.ones(len(x), np.float32), (x, y)), shape=(g.N, g.N)) + extra) > 0).astype(np.float32).tocsr()
            HPOSimilar.tw = twins(Graph(), "disease")
        self.g = g

    def column(self, h, drugs, diseases):
        import scipy.sparse as sp
        from graph import INDICATION
        D = np.asarray(diseases, int)
        S = (HPOSimilar.ph[D] @ HPOSimilar.ph.T).tolil()  # diseases × nodes: shared phenotypes
        for c, d in enumerate(D):
            for t in HPOSimilar.tw.get(int(d), {int(d)}):
                S[c, t] = 0
        I = (h.M(INDICATION).tocsr()[np.asarray(drugs, int)] > 0).astype(np.float32)
        return np.log1p((I @ sp.csr_matrix(S).T).toarray()).astype(np.float32)


class TargetsIn:
    """N6 (and L3b): log(1 + the drug's targets among the disease's genes of a table)."""

    def __init__(self, g, table):
        import scipy.sparse as sp
        o = pd.read_parquet(ROOT / "data" / "ext" / table)
        o = o[(o.disease < g.N) & (o.gene < g.N)]
        self.G = sp.csr_matrix((np.ones(len(o), np.float32), (o.disease, o.gene)), shape=(g.N, g.N))
        self.T = (g.M("drug_protein|moa").tocsr()[g.drugs] > 0).astype(np.float32)

    def column(self, h, drugs, diseases):
        return np.log1p((self.T @ self.G[np.asarray(diseases, int)].T).toarray()).astype(np.float32)


class BBB:
    """N7: ±1 for BBB+ / BBB− drugs in nervous-system or psychiatric diseases, else 0."""

    def __init__(self, g):
        from graph import Graph
        from sources.infections import mondo, mondo_ids
        parents, _, _ = mondo()
        roots = {5071, 2025}  # MONDO nervous system disorder, psychiatric disorder

        def under(m, seen=None):
            seen = seen or set()
            if m in roots:
                return True
            return any(under(p, seen | {m}) for p in parents.get(m, ()) if p not in seen)
        full = Graph()
        self.cns = {int(k) for k, i, t in zip(full.nodes.node_index, full.nodes.node_id, full.nodes.node_type)
                    if t == "disease" and any(under(m) for m in mondo_ids(i))}
        b = pd.read_parquet(ROOT / "data" / "ext" / "b3db.parquet")
        self.bbb = dict(zip(b.drug.astype(int), b.bbb.astype(float)))

    def column(self, h, drugs, diseases):
        v = np.array([self.bbb.get(int(f), 0.0) for f in drugs], np.float32)
        return np.stack([v if int(d) in self.cns else np.zeros_like(v) for d in diseases], axis=1)


ITEMS = {"n5": ("current HPO phenotypes", lambda g: HPOSimilar(g)),
         "n6": ("MGI mouse-model genes", lambda g: TargetsIn(g, "mgi_models.parquet")),
         "n7": ("B3DB blood–brain barrier", lambda g: BBB(g))}
ITEM = sys.argv[1] if len(sys.argv) > 1 and sys.argv[1] in ITEMS else "n5"
OUT = ROOT / "data" / "eval" / f"{ITEM}_feature"


STATE = {"G": None, "drugs": None}
_orig = tx.score_paths


def patched(h, diseases, paths):
    F = _orig(h, diseases, paths)
    if STATE["G"] is None:
        STATE["G"] = ITEMS[ITEM][1](h)  # the first call of build() is on the split's training graph
    return np.concatenate([F, STATE["G"].column(h, h.drugs, diseases)[:, :, None]], axis=2)


def run(split, seed, kg, types):
    task = f"zero_shot_{seed}" if split == "complex_disease" else f"area_{split}"
    f = OUT / f"{task}.json"
    if f.exists():
        return json.loads(f.read_text())
    t0 = time.time()
    STATE["G"] = None  # rebuilt on each split's training graph
    xc.OUT = OUT / "rows"
    d = xc.build(kg, types, split, seed)
    X, Y, G = d["X"], d["Y"], d["groups"]
    k = d["F_test"].shape[2] - 1
    directed = (d["F_test"][:, :, k] != 0).any(axis=0)
    res = {}
    for name, Xs, Fs in ((W, X, d["F_test"]), (P0, np.delete(X, k, axis=1), d["F_test"][:, :, :k])):
        predict = xc._fit(xc.PRODUCTION_MODEL, Xs, Y, G, seed)
        S = xc.scores_for(predict, Fs, d["degree"], d["gf"])
        a, n = pf.recall_rows(S, list(d["test_pos"]), d["degree"], np.random.default_rng(0))
        res[name] = {"all": a.tolist(), "no_prior_use": n.tolist()}
    out = {"task": task, "diseases": [int(x) for x in d["diseases"]], "directed": directed.tolist(), **res}
    OUT.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out), "utf-8")
    print(f"  {task}: done, {int(directed.sum())} of {len(directed)} test diseases with a directed pair "
          f"({time.time() - t0:.0f} s)", flush=True)
    return out


W, P0 = f"with {ITEMS[ITEM][0]}", f"production (same rows, no {ITEMS[ITEM][0]})"
GROUPS = {"benchmark (zero-shot 1-5)": [f"zero_shot_{s}" for s in range(1, 6)], "areas": [f"area_{a}" for a in tx.AREAS]}


def summary():
    rng = np.random.default_rng(0)
    res = {}
    for gname, tasks in GROUPS.items():
        rs = [json.loads((OUT / f"{t}.json").read_text()) for t in tasks if (OUT / f"{t}.json").exists()]
        if len(rs) < len(tasks):
            res[gname] = {"done": f"{len(rs)}/{len(tasks)}"}
            continue
        cat = lambda key, sub, only=False: np.concatenate(
            [np.array(r[key][sub], float)[np.array(r["directed"]) if only else slice(None)] for r in rs])
        mean_task = lambda key, sub: float(np.mean([np.nanmean(np.array(r[key][sub], float)) for r in rs]))
        res[gname] = {"all_diff": pf.boot_diff(cat(W, "all"), cat(P0, "all"), rng),
                      "novel_diff": pf.boot_diff(cat(W, "no_prior_use"), cat(P0, "no_prior_use"), rng),
                      "directed_all_diff": pf.boot_diff(cat(W, "all", True), cat(P0, "all", True), rng),
                      "all_with": mean_task(W, "all"), "all_prod": mean_task(P0, "all"),
                      "novel_with": mean_task(W, "no_prior_use"), "novel_prod": mean_task(P0, "no_prior_use"),
                      "directed_diseases": int(sum(sum(r["directed"]) for r in rs)),
                      "diseases": int(sum(len(r["directed"]) for r in rs))}
    b, a = res.get("benchmark (zero-shot 1-5)", {}), res.get("areas", {})
    ok = False
    if b.get("all_diff") and a.get("all_diff"):
        up = lambda x: x["all_diff"][1] > 0
        flat = lambda x: x["all_diff"][0] >= -0.005
        novel_ok = all(x["novel_diff"] is None or x["novel_diff"][0] >= -0.01 for x in (b, a))
        ok = bool(((up(b) and flat(a)) or (up(a) and flat(b))) and novel_ok)
    res["adopted"] = ok
    return res


def report():
    res = summary()
    fmt = lambda x: "—" if not x else f"{x[0]:+.4f} ({x[1]:+.4f} to {x[2]:+.4f}; {x[3]} diseases)"
    L = ["# O2 feature test: PRISM lineage sensitivity as a model column\n",
         f"Generated by `py -3 -m evaluation.n_columns {ITEM}`. Rule fixed before running (M2's): `docs/IMPROVEMENT_PLAN.md`, "
         "«Registered 24-10-2026: N5–N7». Production rows and model with and without the "
         "column, paired by disease; mean of the task means; differences pooled over the diseases of the group.\n",
         "| Group | Recall@100, all (production → with) | Difference, all | Difference, no prior use | "
         "Difference, diseases where the column is non-zero |", "|---|---|---|---|---|"]
    for gname in GROUPS:
        r = res[gname]
        if "done" in r:
            L.append(f"| {gname} | {r['done']} tasks done | | | |")
            continue
        L.append(f"| {gname} | {r['all_prod']:.4f} → {r['all_with']:.4f} | {fmt(r['all_diff'])} | {fmt(r['novel_diff'])} | "
                 f"{fmt(r['directed_all_diff'])} ({r['directed_diseases']} of {r['diseases']}) |")
    L += ["", f"**Decision (rule fixed before running):** {'adopted' if res['adopted'] else 'not adopted'}.\n"]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / f"{ITEM.upper()}_FEATURE.md").write_text(text, "utf-8")
    (OUT / "summary.json").write_text(json.dumps(res, indent=1), "utf-8") if OUT.exists() else None
    print(text)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    args = [a for a in sys.argv[1:] if a not in ITEMS]
    if args == ["report"]:
        return report()
    tx.score_paths = patched
    xc.score_paths = patched
    tasks = [(args[i], int(args[i + 1])) for i in range(0, len(args), 2)] if args else pf.TASKS
    kg = tx.load_kg()
    types = tx.node_types(kg)
    for split, seed in tasks:
        run(split, seed, kg, types)
    report()


if __name__ == "__main__":
    main()
