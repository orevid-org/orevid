"""G1: curated gene–disease validity (GenCC) for the shared-genes line.

Protocol and decision rules: docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, before any run: G1» (fixed and pushed
on the work PC before any GenCC file was downloaded). Production rule today: ≥ 3 shared genes, p < 1e-5, hub genes
(≥ 500 diseases) left out (R1b, P48). Variants: G1a contradicted links out; G1b only the supported genes of a covered
test disease; G1c as G1b on both sides. Same top 25, sections, sets and intervals as evaluation/evidence_ci.py.

    py -3 -m evaluation.g1_gencc        # -> data/eval/g1_gencc.json, docs/G1_GENCC.md
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import hypergeom

from evaluation.evidence_ci import (HUB_DISEASES, Genes, bench_task, bootstrap, drop_hubs, external_task, fmt,
                                    node_types_all)

ROOT = Path(__file__).resolve().parent.parent
MIN_SUPPORTED = 3
EDGE_SOURCES = ["disease_protein", "ot_assoc:genetic_association", "ot_assoc:genetic_literature",
                "ot_assoc:somatic_mutation", "ot_assoc:animal_model", "ot_assoc:rna_expression"]


def load_gencc():
    t = pd.read_parquet(ROOT / "data" / "ext" / "gencc.parquet")
    sup = t[t["class"] == "supported"].groupby("disease").gene.apply(set).to_dict()
    covered = {d: s for d, s in sup.items() if len(s) >= MIN_SUPPORTED}
    contra = t[t["class"] == "contradicted"]
    return t, covered, list(zip(contra.disease.astype(int), contra.gene.astype(int)))


def custom_table(gn, D, d_sets, m_sets=None):
    """Shared-genes table with the D side (and optionally the M side) given as gene sets; same rule as production."""
    n = gn.g.N
    def mat(sets):
        r, c = [], []
        for i, s in enumerate(sets):
            r += [i] * len(s)
            c += list(s)
        return sp.csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=(len(sets), n))
    DS = mat(d_sets)
    GM = gn._GM if m_sets is None else mat(m_sets)
    K = (DS @ GM.T).toarray().astype(np.int32)
    posM = {m: j for j, m in enumerate(gn.M)}
    for i, d in enumerate(D):
        if d in posM:
            K[i, posM[d]] = 0
    ngD = np.array([len(s) for s in d_sets])[:, None]
    ngM = (np.asarray(GM.sum(axis=1)).ravel() if m_sets is not None else gn.ngM.ravel())[None, :]
    p = np.ones(K.shape)
    m = K >= 3
    p[m] = hypergeom.sf(K[m] - 1, gn.N, np.broadcast_to(ngM, K.shape)[m], np.broadcast_to(ngD, K.shape)[m])
    sig = (K >= 3) & (p < 1e-5)
    return np.asarray(gn.IndM @ sig.T.astype(np.float32)) > 0


def tables(t, covered, contra):
    from evaluation.txgnn import GENES_EXTRA
    g, D = t["g"], t["D"]
    G = drop_hubs(g.M(GENES_EXTRA))
    gn = Genes(g, D)  # production (hub genes out)
    gn._GM = G[gn.M]
    prod = gn.table("base")
    deg = np.diff(g.M(GENES_EXTRA).tocsc().indptr)
    hub = (node_types_all()[:len(deg)] == "gene/protein") & (deg >= HUB_DISEASES)
    # G1a: contradicted edges out (both directions)
    cc = [(d, x) for d, x in contra if d < g.N and x < g.N]
    if cc:
        r = [a for a, b in cc] + [b for a, b in cc]
        c = [b for a, b in cc] + [a for a, b in cc]
        C = sp.csr_matrix((np.ones(len(r), np.float32), (r, c)), shape=G.shape)
        Ga = (G - G.multiply(C > 0)).tocsr()
        Ga.eliminate_zeros()
    else:
        Ga = G
    gna = Genes(g, D, G=Ga)
    gna._GM = Ga[gna.M]
    t1a = gna.table("base")
    # G1b / G1c: supported genes (non-hub) of covered diseases, on the production matrix or on G1a's
    sup = {d: {x for x in s if x < g.N and not hub[x]} for d, s in covered.items()}

    def curated(base, gbase):
        row = lambda d: set(base.indices[base.indptr[d]:base.indptr[d + 1]].tolist())
        d_sets = [sup[int(d)] if int(d) in sup else row(int(d)) for d in D]
        m_sets = [sup[int(m)] if int(m) in sup else row(int(m)) for m in gbase.M]
        return custom_table(gbase, D, d_sets), custom_table(gbase, D, d_sets, m_sets)
    t1b, t1c = curated(G, gn)
    t1ab, t1ac = curated(Ga, gna)
    n_cov = int(sum(int(d) in sup for d in D))
    return prod, {"G1a": t1a, "G1b": t1b, "G1c": t1c, "G1a+G1b": t1ab, "G1a+G1c": t1ac}, n_cov


def split_counts(top, V, T0, T1, hits):
    out = {k: np.zeros((len(top), 2), int) for k in ("removed", "added", "kept", "after", "before")}
    for c, idx in enumerate(top):
        v, a, b, h = V[idx, c], T0[idx, c], T1[idx, c], hits[c][idx]
        for k, msk in (("removed", ~v & a & ~b), ("added", ~v & ~a & b), ("kept", ~v & a & b), ("after", ~v & b),
                       ("before", ~v & a)):
            out[k][c] = (h[msk].sum(), msk.sum())
    return out


def edge_report(gencc):
    from graph import Graph
    g = Graph()
    cls = {(int(d), int(x)): c for d, x, c in zip(gencc.disease, gencc.gene, gencc["class"])}
    rows = []
    for rel in EDGE_SOURCES:
        M = g.M(rel).tocoo()
        types = g.types
        k = {"supported": 0, "limited": 0, "contradicted": 0, "not in GenCC": 0}
        for a, b in zip(M.row, M.col):
            if types[a] == "disease" and types[b] == "gene/protein":
                k[cls.get((int(a), int(b)), "not in GenCC")] += 1
        rows.append({"source": rel, **k})
    return rows


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.txgnn import AREAS, load_kg, node_types
    t0 = time.time()
    gencc, covered, contra = load_gencc()
    kg = load_kg()
    types = node_types(kg)
    rng = np.random.default_rng(0)
    store, cov = {}, {}
    groups = {"benchmark": [("complex_disease", s) for s in (1, 2, 3, 4, 5)], "areas": [(a, 1) for a in AREAS]}
    tasks = [(grp, sp_, se) for grp, lst in groups.items() for sp_, se in lst] + [("external", None, None)]
    for grp, split, seed in tasks:
        t = external_task(kg, types) if grp == "external" else bench_task(kg, types, split, seed)
        prod, variants, n_cov = tables(t, covered, contra)
        cov[grp] = cov.get(grp, 0) + n_cov
        cov[grp + " diseases"] = cov.get(grp + " diseases", 0) + len(t["D"])
        for name, T in variants.items():
            for k, v in split_counts(t["top"], t["V"], prod, T, t["hits"]).items():
                store.setdefault((grp, name, k), []).append(tuple(v.T))
        print(f"  {grp} {split} {seed} ({time.time() - t0:.0f} s)", flush=True)
    res = {"|".join(k): bootstrap(v, rng) for k, v in store.items()}
    out = {"results": res, "coverage": cov, "gencc_pairs": gencc["class"].value_counts().to_dict(),
           "covered_graph_diseases": len(covered), "edges": edge_report(gencc)}
    (ROOT / "data" / "eval" / "g1_gencc.json").write_text(json.dumps(out, indent=1), "utf-8")
    report(out)


def decide(R, name):
    get = lambda grp, k: R.get(f"{grp}|{name}|{k}")
    notes, ok = [], True
    inf = [grp for grp in ("benchmark", "areas") if (get(grp, "removed") or (0, 0, 0, 0))[3] >= 20]
    if inf:
        r1 = all(get(grp, "removed")[0] < 0.05 for grp in inf) and len(inf) == 2
        if len(inf) == 1:  # one informative set: it must pass, and the external exam too
            r1 = get(inf[0], "removed")[0] < 0.05 and (get("external", "removed") or (0, 0, 0, 0))[0] < 0.05
    else:
        e = get("external", "removed")
        r1 = bool(e and e[3] >= 20 and e[0] < 0.05)
        notes.append("benchmark and areas uninformative (< 20 removed)" + ("" if r1 else "; not enough data"))
    ok &= r1
    er, ek = get("external", "removed"), get("external", "kept")
    r2 = er is None or er[3] == 0 or (ek and er[0] <= ek[0])
    ok &= bool(r2)
    r3 = True
    for grp in ("benchmark", "areas", "external"):
        a, b = get(grp, "added"), get(grp, "before")
        if a and a[3] >= 20 and b and a[0] < b[1]:
            r3 = False
    ok &= r3
    af, bf = get("external", "after"), get("external", "before")
    r4 = not (af and bf and af[0] < bf[1])
    ok &= r4
    return bool(ok), {"rule 1 (P29 on removed)": bool(r1), "rule 2 (external removed ≤ kept)": bool(r2),
                      "rule 3 (added not below the section)": bool(r3), "rule 4 (external section not lower)": bool(r4)}, notes


def report(out=None):
    out = out or json.loads((ROOT / "data" / "eval" / "g1_gencc.json").read_text("utf-8"))
    R = out["results"]
    L = ["# G1. Curated gene–disease validity (GenCC) for the shared-genes line\n",
         "Generated by `py -3 -m evaluation.g1_gencc`. Protocol and decision rules: `docs/IMPROVEMENT_PLAN.md`, "
         "«Registered 23-10-2026, before any run: G1». GenCC submissions file downloaded 23-10-2026 (`data/raw/gencc/VERSION`; "
         "CC0 1.0, attribution to GenCC and its submitters requested).\n", "## Coverage (reported before any result)\n",
         f"- GenCC graph pairs: " + ", ".join(f"{k} {v:,}" for k, v in out["gencc_pairs"].items()),
         f"- Graph diseases covered (≥ {MIN_SUPPORTED} supported genes): {out['covered_graph_diseases']:,}",
         "- Covered test diseases: " + ", ".join(f"{g_} {out['coverage'][g_]:,} of {out['coverage'][g_ + ' diseases']:,}"
                                                  for g_ in ("benchmark", "areas", "external")), "",
         "Graph gene–disease edges by GenCC class:\n", "| Edge source | Supported | Limited | Contradicted | Not in GenCC |",
         "|---|---|---|---|---|"]
    L += [f"| {e['source']} | {e['supported']:,} | {e['limited']:,} | {e['contradicted']:,} | {e['not in GenCC']:,} |" for e in out["edges"]]
    L += ["", "## Variants\n", "Hit rate in % (95 % interval resampling diseases; candidates) of the shared-genes candidates each "
          "variant removes, adds and keeps, and of the section before and after.\n",
          "| Variant | Set | Removed | Added | Kept | Section before | Section after |", "|---|---|---|---|---|---|---|"]
    decisions = {}
    for name in ("G1a", "G1b", "G1c", "G1a+G1b", "G1a+G1c"):
        for grp in ("benchmark", "areas", "external"):
            L.append(f"| {name} | {grp} | " + " | ".join(fmt(R.get(f"{grp}|{name}|{k}")) for k in
                                                       ("removed", "added", "kept", "before", "after")) + " |")
        ok, rules, notes = decide(R, name)
        decisions[name] = ok
        L.append(f"| **{name}** | rules | " + "; ".join(f"{k}: {'yes' if v else 'no'}" for k, v in rules.items())
                 + (f" ({'; '.join(notes)})" if notes else "") + f" | → **{'passes' if ok else 'does not pass'}** | | | |")
    single = [k for k in ("G1c", "G1b", "G1a") if decisions[k]]
    chosen = single[0] if single else None
    if chosen and decisions["G1a"] and chosen != "G1a" and decisions.get(f"G1a+{chosen}"):
        chosen = f"G1a+{chosen}"
    L += ["", f"**Decision:** {(chosen + ' adopted') if chosen else 'no variant adopted; the website keeps the rule of P48'}."]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "G1_GENCC.md").write_text(text, "utf-8")
    print(text)
    return chosen


if __name__ == "__main__":
    report() if sys.argv[1:] == ["report"] else main()
