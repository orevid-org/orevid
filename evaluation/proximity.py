"""N2: network proximity between a drug's targets and a disease's genes (Guney et al., Nat Commun 2016).

Protocol and decision rule: docs/IMPROVEMENT_PLAN.md, «Registered 22-10-2026, before any run» (N2). Network:
PrimeKG protein_protein ∪ STRING ≥ 700 among genes; disease genes GENES_EXTRA of the training graph; drug targets
drug_protein:target. Closest distance = mean over the drug's targets of the shortest distance to any disease gene,
z-scored against 50 random gene sets of the same size and degree distribution; proximity = −z.

    py -3 -m evaluation.proximity        # -> data/eval/proximity.json, docs/PROXIMITY.md
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.sparse.csgraph import dijkstra

from evaluation.evidence_ci import node_types_all

ROOT = Path(__file__).resolve().parent.parent
N_RANDOM, MIN_GENES, BIN = 50, 5, 100
TASKS = ROOT / "benchmark" / "tasks"
PRED = ROOT / "benchmark" / "predictions" / "orevid_lambdarank"


class Proximity:
    def __init__(self, g, seed=0):
        from evaluation.txgnn import GENES_EXTRA
        is_gene = node_types_all()[: g.N] == "gene/protein"
        A = (g.M("protein_protein") + g.M("string")).tocsr()
        mask = sp.diags(is_gene.astype(np.float32))
        A = (mask @ A @ mask).tocsr()
        A.data[:] = 1
        A.eliminate_zeros()
        self.A = A
        deg = np.diff(A.indptr)
        self.net = np.flatnonzero(is_gene & (deg > 0))
        self.in_net = np.zeros(g.N, bool)
        self.in_net[self.net] = True
        # degree bins of ~BIN genes, for the random sets
        order = self.net[np.argsort(deg[self.net], kind="stable")]
        self.bin_of = np.full(g.N, -1)
        self.bins = [order[i:i + BIN] for i in range(0, len(order), BIN)]
        for b, genes in enumerate(self.bins):
            self.bin_of[genes] = b
        self.G = g.M(GENES_EXTRA).tocsr()
        T = g.M("drug_protein:target").tocsr()[g.drugs]
        T = T.multiply(self.in_net[None, :]).tocsr()
        T.data[:] = 1
        n_t = np.asarray(T.sum(axis=1)).ravel()
        self.has_target = n_t > 0
        self.Tn = sp.diags(1 / np.maximum(n_t, 1)) @ T  # row-normalised: mean over targets
        # distances computed on the gene network only (same result: other nodes are unreachable), much faster
        self.local = np.full(g.N, -1)
        self.local[self.net] = np.arange(len(self.net))
        self.A_net = A[self.net][:, self.net].tocsr()
        self.Tn_net = self.Tn.tocsc()[:, self.net].tocsr()
        self.rng = np.random.default_rng(seed)

    def _dc(self, genes):
        d = dijkstra(self.A_net, directed=False, indices=self.local[genes], unweighted=True, min_only=True)
        finite = np.isfinite(d)
        d[~finite] = (d[finite].max() + 1) if finite.any() else 50
        return np.asarray(self.Tn_net @ d).ravel()

    def score(self, disease):
        """−z per drug (np.nan if the disease has < MIN_GENES genes in the network); drugs without targets: nan."""
        genes = np.intersect1d(self.G[disease].indices, self.net)
        if len(genes) < MIN_GENES:
            return None
        obs = self._dc(genes)
        bins = self.bin_of[genes]
        rand = np.empty((N_RANDOM, len(obs)))
        for k in range(N_RANDOM):
            pick = np.array([self.rng.choice(self.bins[b]) for b in bins])
            rand[k] = self._dc(np.unique(pick))
        sd = rand.std(axis=0)
        z = np.where(sd > 1e-6, (obs - rand.mean(axis=0)) / np.maximum(sd, 1e-6), 0.0)  # no spread: no evidence
        s = -z
        s[~self.has_target] = np.nan
        return s


def recall_at(scores, pos, k=100, rng=None):
    order = np.argsort(-(scores + rng.random(len(scores)) * 1e-9))
    return len(set(order[:k].tolist()) & pos) / len(pos), len(set(order[:len(pos)].tolist()) & pos) / len(pos)


def boot(v, rng, B=2000):
    v = np.asarray(v, float)
    if not len(v):
        return None
    m = rng.choice(v, (B, len(v))).mean(axis=1)
    return float(v.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def evaluate_rows(rows, rng):
    out = {}
    for key in ("prox", "orevid", "random"):
        out[key] = boot([r[key] for r in rows], rng)
    out["prox_p"] = boot([r["prox_p"] for r in rows], rng)
    out["prox_minus_random"] = boot([r["prox"] - r["random"] for r in rows], rng)
    out["diseases"] = len(rows)
    return out


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.txgnn import TxGraph, load_kg, node_types, split_edges
    t0 = time.time()
    kg = load_kg()
    types = node_types(kg)
    rng = np.random.default_rng(0)
    res = {}
    # (a1) benchmark, zero-shot seeds 1-5, scenario C
    rowsC, rowsAll, n_dis = [], [], 0
    for s in range(1, 6):
        task = json.loads((TASKS / f"zero_shot_{s}.json").read_text())
        g = TxGraph(kg, split_edges(kg, "complex_disease", s) == 0, types, extra=True)
        P = Proximity(g, seed=s)
        drugs = np.asarray(task["drugs"])
        assert np.array_equal(np.sort(drugs), np.sort(g.drugs)), "task drugs differ from the graph's"
        pos_of = {int(f): j for j, f in enumerate(g.drugs)}
        idx = np.array([pos_of[int(f)] for f in drugs])  # graph order -> task order
        pop = np.asarray(task["train_popularity"])
        pred = pd.read_parquet(PRED / f"zero_shot_{s}.parquet")
        for r in task["diseases"]:
            d = int(r["disease"])
            sc = P.score(d)
            if sc is None:
                continue
            n_dis += 1
            sc = sc[idx]
            has = ~np.isnan(sc)
            pos_all = np.isin(drugs, r["positives"])
            ore = pred[pred.disease == d].set_index("drug").score.reindex(drugs).to_numpy()
            ore = np.nan_to_num(ore, nan=np.nanmin(ore) - 1 if np.isfinite(ore).any() else 0)
            for rows, novel in ((rowsC, True), (rowsAll, False)):
                cand = has & ~(pos_all & (pop > 0)) if novel else has.copy()
                pos = pos_all & cand & ((pop == 0) if novel else True)
                if not pos.any():
                    continue
                ci = np.flatnonzero(cand)
                Pset = set(np.flatnonzero(pos[ci]).tolist())
                a, b = recall_at(sc[ci], Pset, rng=rng)
                c, e = recall_at(ore[ci], Pset, rng=rng)
                rows.append({"prox": a, "prox_p": b, "orevid": c, "orevid_p": e, "random": min(1.0, 100 / len(ci))})
        print(f"  zero_shot_{s}: {n_dis} diseases with a score so far ({time.time() - t0:.0f} s)", flush=True)
    res["benchmark, scenario C"] = evaluate_rows(rowsC, rng)
    res["benchmark, all positives (exploratory)"] = evaluate_rows(rowsAll, rng)
    # (a2) temporal 2015, primary positives of drugs with no indication at the cutoff
    from evaluation.temporal import clinical, cutoff_graph, cutoff_mask
    z = dict(np.load(ROOT / "data" / "eval" / "temporal_C2015_extra_no_ctd.npz", allow_pickle=True))
    clin, attrs = clinical()
    keep, _, _ = cutoff_mask(kg, clin, 2015)
    g = cutoff_graph(kg, types, keep, attrs, 2015)
    assert np.array_equal(g.drugs, z["drugs"])
    P = Proximity(g, seed=2015)
    prior = z["degree"] > 0
    rowsT, rowsTall = [], []
    for c, d in enumerate(z["D"]):
        sc = P.score(int(d))
        if sc is None:
            continue
        has = ~np.isnan(sc)
        pos_j = np.asarray(z["pos|primary"][c], int)
        ore = z["S_spec"][:, c].astype(float)
        for rows, novel in ((rowsT, True), (rowsTall, False)):
            cand = z["cand"][:, c] & has
            pj = pos_j
            if novel:
                cand &= ~prior | ~np.isin(np.arange(len(cand)), pos_j)
                pj = pos_j[~prior[pos_j]]
            ci = np.flatnonzero(cand)
            Pset = set(np.flatnonzero(np.isin(ci, pj)).tolist())
            if not Pset:
                continue
            a, b = recall_at(sc[ci], Pset, rng=rng)
            e, f = recall_at(ore[ci], Pset, rng=rng)
            rows.append({"prox": a, "prox_p": b, "orevid": e, "orevid_p": f, "random": min(1.0, 100 / len(ci))})
        if (c + 1) % 100 == 0:
            print(f"  temporal {c + 1}/{len(z['D'])} ({time.time() - t0:.0f} s)", flush=True)
    res["temporal 2015, no prior indication"] = evaluate_rows(rowsT, rng)
    res["temporal 2015, all positives (exploratory)"] = evaluate_rows(rowsTall, rng)
    bc, tc = res["benchmark, scenario C"], res["temporal 2015, no prior indication"]
    above = bool(bc["prox"] and tc["prox"] and bc["prox"][0] > bc["random"][0] and tc["prox"][0] > tc["random"][0])
    out = {"results": res, "above_random_both": above}
    (ROOT / "data" / "eval" / "proximity.json").write_text(json.dumps(out, indent=1), "utf-8")
    report(out)


def report(out=None):
    out = out or json.loads((ROOT / "data" / "eval" / "proximity.json").read_text("utf-8"))
    f = lambda x: "—" if not x else f"{x[0]:.3f} ({x[1]:.3f}–{x[2]:.3f})"
    L = ["# N2. Network proximity (drug targets ↔ disease genes)\n",
         "Generated by `py -3 -m evaluation.proximity`. Protocol and rule: `docs/IMPROVEMENT_PLAN.md`, «Registered "
         "22-10-2026, before any run» (N2). Recall@100 among candidate drugs with a target in the protein network "
         "(mean over diseases, 95 % interval resampling diseases). Random = 100 / candidates. OREVID = the production "
         "ranking of the same candidates.\n",
         "| Test | Diseases | Proximity | Proximity − random | Random | OREVID | Proximity Precision@K |",
         "|---|---|---|---|---|---|---|"]
    for k, v in out["results"].items():
        L.append(f"| {k} | {v['diseases']} | {f(v['prox'])} | {f(v['prox_minus_random'])} | {f(v['random'])} | "
                 f"{f(v['orevid'])} | {f(v['prox_p'])} |")
    L += ["", f"**Decision:** proximity is {'above' if out['above_random_both'] else 'not above'} random on both the "
              f"benchmark (scenario C) and the temporal set, so (b), the LambdaRank-feature test, "
              f"{'is run' if out['above_random_both'] else 'is not run (rule fixed before running)'}.\n",
          "**NeDRex:** see the section below.\n"]
    text = "\n".join(L) + "\n"
    p = ROOT / "docs" / "PROXIMITY.md"
    old = p.read_text("utf-8") if p.exists() else ""
    tail = "## NeDRex"
    kept = old[old.index(tail):] if tail in old else ""
    p.write_text(text + ("\n" + kept if kept else ""), "utf-8")
    print(text)


if __name__ == "__main__":
    report() if sys.argv[1:] == ["report"] else main()
