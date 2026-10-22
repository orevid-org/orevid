"""Evidence levels of the candidates (the ones the web shows) and their hit rate in the validation.

Each candidate has a profile of four signals (the model paths an expert can check):
    V   indicated for a neighbouring disease in the ontology (1 step)
    v2  indicated for a disease 2 steps away in the ontology (without V)
    D   acts on a gene associated with the disease
    S   indicated for a disease with ≥ 3 associated genes in common

Levels (fixed with the complex_disease bench, seeds 1-5, before looking at the areas; DECISIONS.md P27):
    A  convergent evidence: S and, in addition, D, V or v2
    B  one solid line: S alone, or V (with or without D)
    C  weak: only v2
    -  insufficient (not shown): only D, v2 + D, or none. «Acts on a disease gene» with nothing
       else hits 2-3 %, almost the same as having no path at all (1.5-2 %).

Confirmation criterion, fixed before looking at them: on the 9 hidden areas (different split), the
order A > B > C > - holds and «-» stays below 5 %.

First pass (S = ≥ 3 genes, without significance), hits among the top 25:
                     bench all   bench hard   areas all   areas hard
    A                  41.7 %       26.0 %       32.1 %      15.2 %
    B                  20.3 %       17.0 %       19.5 %      20.9 %
    C                  10.6 %        5.2 %        3.5 %       1.1 %
    -                   1.6 %        2.5 %        1.5 %       1.8 %
«-» < 5 %: confirmed. A > B: NOT confirmed (it flips on the hard areas). C: does not reach 10 % out
of sample. Decision: the web shows no level letters; it shows only A and B (rule `has_evidence`:
indicated in a neighbouring disease or in one with ≥ 3 genes in common). The lines are shown as facts.

Final version (S, in addition, significant: hypergeometric p < 1e-5 over the genes with at least one disease;
P28, corrected in P40), the one this script reproduces and saves in data/eval/log_evidence.txt and
docs/EVIDENCE.md: A 45.3/29.2/37.4/19.5 %; per web section, direct extensions («neighbour»)
49.9/13.3/59.4/24.7 %, table 23.0/27.4/26.3/22.5 %, hidden 4.2/2.9/1.7/1.6 %.

External exam (new approvals the graph did not know, the most realistic test; `external()`), among the
top 25 not known: neighbour 12.6 %, table with specific genes 4.8 % and with generic genes 2.8 % (whole
table 4.3 %), hidden 1.5 % (random < 0.1 %). The test of hiding treatments overestimates the ability to
find new indications (DECISIONS.md P31). A permutation null shows that the shared genes themselves carry no
measurable signal (P38, evaluation/evidence_ci.py).

    py -3 -m evaluation.evidence            # table per level: complex_disease bench and areas
"""
import sys
from collections import defaultdict

import numpy as np
from scipy.stats import hypergeom

LEVELS = ["A", "B", "C", "-"]
# Hypergeometric test of the shared genes: population N = genes with at least one disease, p < 1e-5 (P40,
# docs/HYPERGEOMETRIC_N.md; before: N = every linked node, genes and diseases, p < 1e-6).
MAX_P = 1e-5


def n_genes(G):
    """Number of gene nodes with at least one associated disease (G is symmetric: diseases excluded)."""
    from evaluation.evidence_ci import node_types_all
    deg = np.diff(G.tocsc().indptr)
    return int(((deg > 0) & (node_types_all()[:len(deg)] == "gene/protein")).sum())
SCREEN = ["neighbour", "table", "hidden"]


def level(V, v2, D, S):
    if S and (D or V or v2):
        return "A"
    if S or V:
        return "B"
    if v2 and not D:
        return "C"
    return "-"


def _transform(S, w=0.5):
    """The web score: z per disease minus w × the drug's mean (specificity)."""
    Z = (S - S.mean(axis=0)) / (S.std(axis=0) + 1e-9)
    return Z - w * Z.mean(axis=1, keepdims=True)


def evaluate(split, seed, top=25, significance=True):
    """{(subset, level): [hits, candidates]} among the `top` first of each test disease."""
    from evaluation.txgnn import (CACHE, FEATURE_SETS, GENES_EXTRA, INDICATION, PRODUCTION, TxGraph, load_kg,
                                  score_paths, near_duplicates, split_edges, node_types)
    kg = load_kg()
    paths = FEATURE_SETS[PRODUCTION]
    ix = {c: paths.index(c) for c in paths}
    z = np.load(CACHE / "scores" / f"{split}_{seed}_{PRODUCTION}_r5_rank.npz", allow_pickle=True)
    S = _transform(z["LambdaRank_+_masked"].astype(float))
    diseases, positives = z["diseases"], [np.asarray(p, int) for p in z["positives"]]
    g = TxGraph(kg, split_edges(kg, split, seed) == 0, node_types(kg), extra=True)
    F = score_paths(g, diseases, paths)
    from evaluation.evidence_ci import drop_hubs
    G, I = drop_hubs(g.M(GENES_EXTRA)), g.M(INDICATION).tocsr()  # hub genes out (R1b, P48)
    drug_pos = np.full(g.N, -1)
    drug_pos[g.drugs] = np.arange(len(g.drugs))
    nd = near_duplicates(g, diseases, positives, drug_pos)
    ng = np.diff(G.indptr)
    N = n_genes(G)
    res = defaultdict(lambda: [0, 0])
    for c, d in enumerate(diseases):
        pos = set(positives[c])
        genes_d = set(G.indices[G.indptr[d]:G.indptr[d + 1]])
        for f in np.argsort(-S[:, c])[:top]:
            v = lambda n: F[f, c, ix[n]] > 0
            dr = g.drugs[f]
            dis = I.indices[I.indptr[dr]:I.indptr[dr + 1]]
            # S: ≥ 3 genes in common with a disease where it is indicated, and significant (p < MAX_P; P28, P40)
            shared = 0
            for m in dis:
                if m == d:
                    continue
                k = len(genes_d.intersection(G.indices[G.indptr[m]:G.indptr[m + 1]]))
                if k >= 3 and (not significance or hypergeom.sf(k - 1, N, ng[m], len(genes_d)) < MAX_P):
                    shared = max(shared, k)
            V = v("ontology")
            n = level(V, v("ontology 2-step") and not V, v("target") or v("target (all)") or v("target+ all"),
                      shared >= 3)
            # as the web shows it: V -> «direct extension», S without V -> candidate table, the rest hidden
            screen = "neighbour" if V else ("table" if shared >= 3 else "hidden")
            for sub in ("all",) + (("hard",) if nd[c] < 0.5 else ()):
                for k in (n, screen):
                    res[(sub, k)][0] += f in pos
                    res[(sub, k)][1] += 1
    return res


def table(res, title):
    print(f"\n{title}")
    for sub in ("all", "hard"):
        row = []
        for n in LEVELS + SCREEN:
            h, k = res.get((sub, n), (0, 0))
            row.append(f"{n}: {h / k:5.1%} ({k:5d})" if k else f"{n}:    —")
        print(f"  {sub:9s} " + "   ".join(row))


def external(top=25):
    """External exam (Open Targets approvals that are not in PrimeKG; evaluation/external_exam.py): hit
    rate among the `top` first candidates not known, per web section. It is the most realistic test for
    new candidates (P31). Needs data/eval/external_scores.npz (py -3 -m evaluation.external_exam evaluate)."""
    import pandas as pd
    from evaluation.external_exam import SEED, SPLIT
    from evaluation.txgnn import (CACHE, FEATURE_SETS, GENES_EXTRA, INDICATION, PRODUCTION, TxGraph, load_kg,
                                  score_paths, split_edges, node_types)
    kg = load_kg()
    z = np.load(CACHE / "external_scores.npz", allow_pickle=True)
    diseases, drugs = z["diseases"], z["drugs"]
    S = _transform(z["Ours_ranking"].astype(float))
    all_pos = pd.read_parquet(CACHE / "external_positives.parquet")
    g = TxGraph(kg, split_edges(kg, SPLIT, SEED) == 0, node_types(kg), extra=True)
    paths = FEATURE_SETS[PRODUCTION]
    ix = {c: paths.index(c) for c in paths}
    F = score_paths(g, diseases, paths)
    from evaluation.evidence_ci import drop_hubs
    G, I = drop_hubs(g.M(GENES_EXTRA)), g.M(INDICATION).tocsr()  # hub genes out (R1b, P48)
    ng, gene_degree = np.diff(G.indptr), np.diff(G.tocsc().indptr)
    N = n_genes(G)
    res = defaultdict(lambda: [0, 0])
    for c, d in enumerate(diseases):
        pos, known = set(all_pos[all_pos.disease == d].drug), set(z["known"][c])
        genes_d = set(G.indices[G.indptr[d]:G.indptr[d + 1]])
        for j in [j for j in np.argsort(-S[:, c]) if drugs[j] not in known][:top]:
            f = drugs[j]
            Sb = Se = False
            for m in I.indices[I.indptr[f]:I.indptr[f + 1]]:
                com = genes_d.intersection(G.indices[G.indptr[m]:G.indptr[m + 1]]) if m != d else set()
                if len(com) >= 3 and hypergeom.sf(len(com) - 1, N, ng[m], len(genes_d)) < MAX_P:
                    Sb = True
                    Se |= sum(gene_degree[x] < 100 for x in com) >= 3
            k = "neighbour" if F[j, c, ix["ontology"]] > 0 else ("table, specific genes" if Se else
                                                                 ("table, generic genes" if Sb else "hidden"))
            for kk in (k, "all") + (("whole table",) if k.startswith("table") else ()):
                res[kk][0] += f in pos
                res[kk][1] += 1
    print(f"\nExternal exam ({len(diseases)} diseases), hits among the top {top} not known:")
    for k in ("neighbour", "whole table", "table, specific genes", "table, generic genes", "hidden", "all"):
        h, n = res[k]
        print(f"  {k:24s} {h / max(n, 1):5.1%} ({n})")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.txgnn import AREAS
    total = defaultdict(lambda: [0, 0])
    for s in (1, 2, 3, 4, 5):
        for k, (h, n) in evaluate("complex_disease", s).items():
            total[k][0] += h
            total[k][1] += n
    table(total, "complex_disease bench, seeds 1-5 (where the levels were fixed)")
    areas = defaultdict(lambda: [0, 0])
    for a in AREAS:
        r = evaluate(a, 1)
        for k, (h, n) in r.items():
            areas[k][0] += h
            areas[k][1] += n
        print(f"  {a} done", flush=True)
    table(areas, "Confirmation: 9 hidden areas, seed 1")
    external()


if __name__ == "__main__":
    main()
