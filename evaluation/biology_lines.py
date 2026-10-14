"""A shared-biology line that beats the permutation null (docs/IMPROVEMENT_PLAN.md, I6).

The «shared genes» section (indicated for a disease sharing ≥ 3 genes, hypergeometric p < 1e-5) does not beat a
degree-preserving permutation null (P38, P40). Three replacements were named in the plan; the details the plan
left open are fixed here, before running (this file is committed before its first run):

(a) **pathway overlap**: a disease is enriched in a Reactome pathway (PrimeKG `pathway_protein`) when ≥ 3 of its
    genes are in it and the hypergeometric p < 1e-3 (population: genes in at least one pathway; draws: the
    disease's genes in at least one pathway). The line: indicated for a disease sharing ≥ 3 enriched pathways,
    hypergeometric p < 1e-5 over the pathways enriched in at least one disease (the gene rule, on pathways).
(b) **specificity-weighted gene overlap**: each shared gene weighs log(diseases with genes / diseases of the
    gene) (inverse document frequency); the line: ≥ 3 shared genes and weighted overlap ≥ τ, with τ set on the
    benchmark so that the section has as many candidates as the production rule there (a re-weighting, not a
    larger section).
(c) **genetic-only genes**: the production rule (≥ 3 shared genes, p < 1e-5) on Open Targets genetic
    association and somatic mutation evidence only (`ot_assoc:genetic_association|ot_assoc:somatic_mutation`),
    population the genes with at least one such disease.

Same top 25 and sections as `evaluation/evidence_ci.py` (the line replaces «shared genes»; direct extensions are
unchanged). Selection among (a)-(c): the highest hit rate of the section on the benchmark (zero-shot seeds 1-5,
pooled). The chosen line replaces «shared genes» only if (rule of the plan) on the external exam (1) its hit
rate is above the 99th percentile of the P40 permutation null (`data/eval/evidence_ci.json`), (2) the upper end of
its 95 % interval is not below the production 4.3 %, and (3) its hit rates on the benchmark and the areas are not
lower than production's; and, as the plan's confirmation on the temporal C2 set, (4) its hit rate at the 2015
cutoff (`evaluation/temporal.py`, primary positives) is not lower than production's. If not, the website keeps
the P38 wording.

    py -3 -m evaluation.biology_lines          # -> data/eval/biology_lines.json, docs/BIOLOGY_LINES.md
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.stats import hypergeom

from evaluation.evidence_ci import Genes, bench_task, bootstrap, counts, external_task, fmt, node_types_all

ROOT = Path(__file__).resolve().parent.parent
LINES = ("production", "pathways", "idf", "genetic")
GENETIC = "ot_assoc:genetic_association|ot_assoc:somatic_mutation"
ENRICH_MIN, ENRICH_P = 3, 1e-3


def symmetric(rows, cols, n):
    A = sp.csr_matrix((np.ones(len(rows), np.float32), (rows, cols)), shape=(n, n))
    return ((A + A.T) > 0).astype(np.float32).tocsr()


def pathway_matrix(g, G):
    """Graph-sized symmetric matrix disease <-> enriched pathway."""
    from evaluation.txgnn import GENES_EXTRA  # noqa: F401  (same gene sets as production)
    types = node_types_all()
    is_gene = types == "gene/protein"
    P = g.M("pathway_protein").tocsr()
    in_path = is_gene & (np.diff(P.indptr) > 0)
    pathways = np.flatnonzero((types == "pathway") & (np.diff(P.indptr) > 0))
    Gc = G.tocsr()
    dis = np.flatnonzero((types == "disease") & (np.diff(Gc.indptr) > 0))
    Dg = Gc[dis][:, np.flatnonzero(in_path)]  # diseases × genes in a pathway
    Pg = P[pathways][:, np.flatnonzero(in_path)]  # pathways × genes
    K = (Dg @ Pg.T).tocoo()  # shared genes disease × pathway
    n_draw = np.asarray(Dg.sum(axis=1)).ravel()
    size = np.asarray(Pg.sum(axis=1)).ravel()
    N = int(in_path.sum())
    m = K.data >= ENRICH_MIN
    r, c, k = K.row[m], K.col[m], K.data[m]
    p = hypergeom.sf(k - 1, N, size[c], n_draw[r])
    ok = p < ENRICH_P
    return symmetric(dis[r[ok]], pathways[c[ok]], g.N), "pathway"


class Lines(Genes):
    """Genes with the population N of a different node type (pathways)."""

    def __init__(self, g, D, G, population_type="gene/protein"):
        super().__init__(g, D, G=G)
        deg = np.diff(G.tocsc().indptr)
        self.N = int(((deg > 0) & (node_types_all() == population_type)).sum())


def idf_weights(G):
    G = G.tocsc()
    types = node_types_all()
    deg = np.diff(G.indptr).astype(float)  # diseases per gene (the matrix is symmetric)
    n_dis = int(((types == "disease") & (np.diff(G.tocsr().indptr) > 0)).sum())
    w = np.zeros(G.shape[0])
    m = (types == "gene/protein") & (deg > 0)
    w[m] = np.log(n_dis / deg[m])
    return w


def weighted_overlap(gn, G, w):
    G = G.tocsr()
    W = (G[gn.D] @ sp.diags(w.astype(np.float32)) @ G[gn.M].T).toarray()
    posM = {m: j for j, m in enumerate(gn.M)}
    for i, d in enumerate(gn.D):
        if d in posM:
            W[i, posM[d]] = 0
    return W


def table_from(gn, sig):
    return np.asarray(gn.IndM @ sig.T.astype(np.float32)) > 0


def tables(t, tau=None):
    """drugs × D table of every line (and the weighted overlaps for choosing τ)."""
    from evaluation.txgnn import GENES_EXTRA
    g = t["g"]
    G = g.M(GENES_EXTRA)
    gn = Genes(g, t["D"])
    out = {"production": gn.table("base")}
    Gp, ptype = pathway_matrix(g, G)
    out["pathways"] = Lines(g, t["D"], Gp, ptype).table("base")
    W = weighted_overlap(gn, G, idf_weights(G))
    t["_idf"] = (gn, W)
    if tau is not None:
        out["idf"] = table_from(gn, (gn.K >= 3) & (W >= tau))
    out["genetic"] = Genes(g, t["D"], G=g.M(GENETIC)).table("base")
    return out


def section_counts(t, T, hits_key="hits"):
    return counts(t["top"], t["V"], T, t[hits_key])["table"]


def choose_tau(bench, prod_n):
    """τ such that the idf section has as many candidates on the benchmark as production."""
    allw = np.concatenate([W[gn.K >= 3] for gn, W in (t["_idf"] for t in bench)])
    grid = np.unique(np.quantile(allw, np.linspace(0, 0.999, 400)))
    best = None
    for tau in grid:
        n = sum(section_counts(t, table_from(t["_idf"][0], (t["_idf"][0].K >= 3) & (t["_idf"][1] >= tau)))[:, 1].sum()
                for t in bench)
        if best is None or abs(n - prod_n) < best[1]:
            best = (float(tau), abs(n - prod_n), int(n))
    return best[0], best[2]


def temporal_c2(kg, types, line, tau):
    """Hit rate of the line's section at the 2015 cutoff (primary positives), against production."""
    from evaluation.temporal import clinical, cutoff_graph, cutoff_mask, top25
    z = dict(np.load(ROOT / "data" / "eval" / "temporal_C2015_extra_no_ctd.npz", allow_pickle=True))
    clin, attrs = clinical()
    keep, _, _ = cutoff_mask(kg, clin, 2015)
    g = cutoff_graph(kg, types, keep, attrs, 2015)
    assert np.array_equal(g.drugs, z["drugs"])
    pos = list(z["pos|primary"])
    t = {"g": g, "D": z["D"], "V": z["V"], "top": top25(z),
         "hits": [np.isin(np.arange(len(z["drugs"])), p) for p in pos]}
    T = tables(t, tau)
    return {k: section_counts(t, T[k]) for k in ("production", line)}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.txgnn import AREAS, load_kg, node_types
    t0 = time.time()
    kg = load_kg()
    types = node_types(kg)
    rng = np.random.default_rng(0)
    groups = {"benchmark": [("complex_disease", s) for s in (1, 2, 3, 4, 5)], "areas": [(a, 1) for a in AREAS]}
    tasks = {k: [] for k in groups}
    for group, lst in groups.items():
        for split, seed in lst:
            t = bench_task(kg, types, split, seed)
            t["T"] = tables(t)
            t["g"] = None  # free memory; the tables are kept
            tasks[group].append(t)
            print(f"  {group} {split} {seed} ({time.time() - t0:.0f} s)", flush=True)
    prod_n = sum(section_counts(t, t["T"]["production"])[:, 1].sum() for t in tasks["benchmark"])
    tau, n_idf = choose_tau(tasks["benchmark"], prod_n)
    print(f"  tau {tau:.3f}: {n_idf} idf candidates against {prod_n} production on the benchmark", flush=True)
    for lst in tasks.values():
        for t in lst:
            gn, W = t["_idf"]
            t["T"]["idf"] = table_from(gn, (gn.K >= 3) & (W >= tau))
            t["_idf"] = None
    res = {}
    for group, lst in tasks.items():
        for line in LINES:
            res[f"{group}|{line}"] = bootstrap([tuple(section_counts(t, t["T"][line]).T) for t in lst], rng)
    chosen = max(("pathways", "idf", "genetic"), key=lambda k: res[f"benchmark|{k}"][0])
    ext = external_task(kg, types)
    T = tables(ext, tau)
    for line in LINES:
        res[f"external|{line}"] = bootstrap([tuple(section_counts(ext, T[line]).T)], rng)
        res[f"external, non-obvious|{line}"] = bootstrap([tuple(section_counts(ext, T[line], "hits_new").T)], rng)
    print(f"  external ({time.time() - t0:.0f} s)", flush=True)
    del ext, T
    c2 = temporal_c2(kg, types, chosen, tau)
    for k, v in c2.items():
        res[f"temporal 2015|{k}"] = bootstrap([tuple(v.T)], rng)
    null = json.loads((ROOT / "data" / "eval" / "evidence_ci.json").read_text("utf-8"))["null_table_external"]
    p99 = float(np.percentile([h / n if n else 0.0 for h, n in null], 99))
    e, ep = res[f"external|{chosen}"], res["external|production"]
    checks = {"(1) external above the null's 99th percentile": e[0] > p99,
              "(2) external interval reaches production's rate": e[2] >= ep[0],
              "(3) benchmark not lower": res[f"benchmark|{chosen}"][0] >= res["benchmark|production"][0],
              "(3) areas not lower": res[f"areas|{chosen}"][0] >= res["areas|production"][0],
              "(4) temporal 2015 not lower": res[f"temporal 2015|{chosen}"][0] >= res["temporal 2015|production"][0]}
    out = {"results": res, "chosen": chosen, "tau": tau, "null_p99": p99, "checks": checks,
           "adopted": bool(all(checks.values()))}
    (ROOT / "data" / "eval" / "biology_lines.json").write_text(json.dumps(out, indent=1), "utf-8")
    report(out)


def report(out=None):
    out = out or json.loads((ROOT / "data" / "eval" / "biology_lines.json").read_text("utf-8"))
    R = out["results"]
    names = {"production": "Shared genes (production)", "pathways": "(a) Shared enriched pathways",
             "idf": f"(b) Specificity-weighted genes (τ = {out['tau']:.2f})", "genetic": "(c) Genetic-only genes"}
    L = ["# A shared-biology line against the null (I6)\n",
         "Generated by `py -3 -m evaluation.biology_lines`. Protocol and rule: `docs/IMPROVEMENT_PLAN.md` (I6) and "
         "the module docstring (details fixed before running). Hit rate in % of the section that replaces «shared "
         "genes» among the 25 highest-ranked candidates per disease (95% interval resampling diseases; candidates).\n",
         "| Line | Benchmark | Areas | External | External, non-obvious | Temporal 2015 |", "|---|---|---|---|---|---|"]
    for k in LINES:
        L.append(f"| {names[k]} | " + " | ".join(fmt(R.get(f"{grp}|{k}")) for grp in
                                                 ("benchmark", "areas", "external", "external, non-obvious",
                                                  "temporal 2015")) + " |")
    L += ["", f"Chosen on the benchmark: **{names[out['chosen']]}**. Permutation null of the production line on the "
              f"external exam (P40): 99th percentile {100 * out['null_p99']:.2f} %.\n"]
    L += [f"- {k}: {'yes' if v else 'no'}" for k, v in out["checks"].items()]
    L += ["", f"**Decision:** {'the line replaces «shared genes»' if out['adopted'] else 'not adopted; the website keeps the P38 wording'}."]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "BIOLOGY_LINES.md").write_text(text, "utf-8")
    print(text)


if __name__ == "__main__":
    report() if sys.argv[1:] == ["report"] else main()
