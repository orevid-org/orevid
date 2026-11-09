"""G3: disease genes without animal models, for the shared-genes line.

Protocol and decision rules: docs/IMPROVEMENT_PLAN.md, «Registered 09-11-2026, on the work PC, before any run: G3 …».
First the report fixed before any result (genes per source; diseases whose genes are ≥ 50 % animal-model only;
shown shared-genes candidates resting on an animal-only gene), then the variants:
- G3a: a gene–disease edge whose only source is `ot_assoc:animal_model` does not count (both sides);
- G3b: as G3a, and edges whose only source is `ot_assoc:rna_expression` do not count either.
Each variant is the production rule (≥ 3 shared non-hub genes, p < 1e-5, population = non-hub genes with a disease)
recomputed over the remaining edges, judged with the four rules of G1 (`evaluation/g1_gencc.py`). G3b if it passes,
else G3a. G3c (model column) is a separate run.

    py -3 -m evaluation.g3_animal_genes          # -> data/eval/g3_animal_genes.json, docs/G3_ANIMAL_GENES.md
    py -3 -m evaluation.g3_animal_genes report   # report only, from the json
"""
import json
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.stats import hypergeom

from evaluation.evidence_ci import (HUB_DISEASES, Genes, bench_task, bootstrap, counts, drop_hubs, external_task, fmt,
                                    node_types_all, permuted_genes)
from evaluation.g1_gencc import decide, split_counts

ROOT = Path(__file__).resolve().parent.parent
PRIMEKG = "disease_protein"
GENETIC = "ot_assoc:genetic_association|ot_assoc:genetic_literature|ot_assoc:somatic_mutation"
ANIMAL = "ot_assoc:animal_model"
EXPRESSION = "ot_assoc:rna_expression"
SOURCES = {"PrimeKG (disease_protein)": PRIMEKG, "Open Targets genetic": GENETIC, "Open Targets animal model": ANIMAL,
           "Open Targets expression": EXPRESSION}
VARIANT_RELS = {"G3a": f"{PRIMEKG}|{GENETIC}|{EXPRESSION}", "G3b": f"{PRIMEKG}|{GENETIC}"}
PERMUTATIONS = 100


def matrix(g, rel):
    """Union of the pieces of `rel` that exist in the graph, binary and symmetric (as Graph.M)."""
    pieces = [r for r in rel.split("|") if r in g.relations]
    m = sum(g.M(r) for r in pieces).tocsr()
    m.data[:] = 1
    return m


def source_stats(g, diseases):
    """Per-source gene counts for `diseases`; animal-only and expression-only edges; share of animal-only genes."""
    D = np.asarray(sorted({int(x) for x in diseases}))
    rows = {k: matrix(g, r)[D] for k, r in SOURCES.items()}
    allm = matrix(g, "|".join(SOURCES.values()))[D]
    non_animal = matrix(g, VARIANT_RELS["G3a"])[D]
    human = matrix(g, VARIANT_RELS["G3b"])[D]
    n_all = np.diff(allm.indptr)
    animal_only = n_all - np.diff(non_animal.indptr)
    expr_only = np.diff(non_animal.indptr) - np.diff(human.indptr)
    has = n_all > 0
    out = {"diseases": int(len(D)), "with genes": int(has.sum()),
           "edges": {k: int(m.nnz) for k, m in rows.items()},
           "diseases with the source": {k: int((np.diff(m.indptr) > 0).sum()) for k, m in rows.items()},
           "median genes (diseases with genes)": float(np.median(n_all[has])) if has.any() else 0.0,
           "animal-only edges": int(animal_only.sum()), "expression-only edges": int(expr_only.sum()),
           "all edges": int(n_all.sum()),
           "≥ 50 % animal-only": int((has & (animal_only >= 0.5 * n_all)).sum()),
           "≥ 50 % animal- or expression-only": int((has & (animal_only + expr_only >= 0.5 * n_all)).sum()),
           "left with no gene under G3a": int((has & (np.diff(non_animal.indptr) == 0)).sum()),
           "left with no gene under G3b": int((has & (np.diff(human.indptr) == 0)).sum())}
    return out


class WebRule:
    """The website's shared-genes test (web/build.py, Explainer.explain) on a given gene–disease matrix."""

    def __init__(self, G, types):
        self.G = G.tocsr()
        deg = np.diff(self.G.tocsc().indptr)
        self.hub = (deg >= HUB_DISEASES) & (types[:len(deg)] == "gene/protein")
        self.N = int(((deg > 0) & (types[:len(deg)] == "gene/protein") & ~self.hub).sum())

    def genes(self, x):
        s = self.G.indices[self.G.indptr[x]:self.G.indptr[x + 1]]
        return {int(y) for y in s if not self.hub[y]}

    def valid(self, d, m):
        a, b = self.genes(d), self.genes(m)
        k = len(a & b)
        return k >= 3 and hypergeom.sf(k - 1, self.N, len(b), len(a)) < 1e-5


def website(gfull):
    """Shown shared-genes candidates resting on at least one animal-only gene, and how many keep the line with the
    diseases listed on the page under G3a and G3b (an upper bound on the candidates that would lose it: only the
    three listed diseases are re-tested; a full rebuild re-tests all)."""
    from web.build import DEST, valid_similar
    types = gfull.types
    G0 = matrix(gfull, "|".join(SOURCES.values()))
    A = matrix(gfull, ANIMAL)
    E = matrix(gfull, EXPRESSION)
    Ga, Gb = matrix(gfull, VARIANT_RELS["G3a"]), matrix(gfull, VARIANT_RELS["G3b"])
    animal_only = (A - A.multiply(Ga)).tocsr()
    animal_only.eliminate_zeros()
    w0, wa, wb = WebRule(G0, types), WebRule(Ga, types), WebRule(Gb, types)
    by_name = {}
    for x in np.flatnonzero(types == "disease"):
        by_name.setdefault(gfull.name(int(x)), []).append(int(x))
    ao = lambda x: set(animal_only.indices[animal_only.indptr[x]:animal_only.indptr[x + 1]].tolist())
    tot = Counter()
    nf = None
    for fj in sorted((DEST / "d").glob("*.json")):
        for k, v in json.loads(fj.read_text("utf-8")).items():
            d = int(k)
            rec = []
            for c in v.get("candidates") or []:
                via = [m[0] for e in c.get("expl", []) if e.get("t") == "similar" for m in e.get("m", []) if valid_similar(m)]
                if not via:
                    continue
                tot["shown shared-genes candidates"] += 1
                nodes = [max(by_name.get(n, []), key=lambda m: len(w0.genes(d) & w0.genes(m)), default=None) for n in via]
                nodes = [m for m in nodes if m is not None]
                shared_ao = any((w0.genes(d) & w0.genes(m)) & (ao(d) | ao(m)) for m in nodes)
                keep_a = any(wa.valid(d, m) for m in nodes)
                keep_b = any(wb.valid(d, m) for m in nodes)
                tot["check: valid with the production rule recomputed"] += any(w0.valid(d, m) for m in nodes)
                tot["rest on ≥ 1 animal-only gene"] += shared_ao
                tot["keep the line under G3a (listed diseases)"] += keep_a
                tot["keep the line under G3b (listed diseases)"] += keep_b
                rec.append((c["name"], ", ".join(via), shared_ao, keep_a, keep_b))
            if v["name"] == "neurofibromatosis":
                nf = {"id": d, "genes": {s: int(len(gfull.neighbors(d, r))) for s, r in SOURCES.items()},
                      "animal-only": len(ao(d)), "G3a genes": len(wa.genes(d)), "G3b genes": len(wb.genes(d)),
                      "all (non-hub)": len(w0.genes(d)), "candidates": rec}
    return dict(tot), nf


def task_tables(t, mats):
    g, D = t["g"], t["D"]
    prod = Genes(g, D).table("base")
    out = {}
    for name, G in mats.items():
        out[name] = Genes(g, D, G=G).table("base")
    return prod, out


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from evaluation.txgnn import AREAS, load_kg, node_types
    from web.build import DEST
    t0 = time.time()
    gfull = Graph()
    pages = [int(k) for fj in sorted((DEST / "d").glob("*.json")) for k in json.loads(fj.read_text("utf-8"))]
    stats = {"website pages": source_stats(gfull, pages)}
    web, nf = website(gfull)
    print(f"  website ({time.time() - t0:.0f} s): {web}", flush=True)
    kg = load_kg()
    types = node_types(kg)
    rng = np.random.default_rng(0)
    store, test_d = {}, {}
    groups = {"benchmark": [("complex_disease", s) for s in (1, 2, 3, 4, 5)], "areas": [(a, 1) for a in AREAS]}
    tasks = [(grp, sp_, se) for grp, lst in groups.items() for sp_, se in lst] + [("external", None, None)]
    ext = None
    for grp, split, seed in tasks:
        t = external_task(kg, types) if grp == "external" else bench_task(kg, types, split, seed)
        g = t["g"]
        test_d.setdefault(grp, []).append((g, t["D"]))
        mats = {k: matrix(g, r) for k, r in VARIANT_RELS.items()}
        prod, variants = task_tables(t, mats)
        for name, T in variants.items():
            for k, v in split_counts(t["top"], t["V"], prod, T, t["hits"]).items():
                store.setdefault((grp, name, k), []).append(tuple(v.T))
        if grp == "external":
            ext = (t, mats)
        print(f"  {grp} {split} {seed} ({time.time() - t0:.0f} s)", flush=True)
    for grp, lst in test_d.items():  # per-source report on each set's own graph (first task's graph, all its diseases)
        g = lst[0][0]
        stats[f"{grp} test diseases"] = source_stats(g, np.concatenate([D for _, D in lst]))
    res = {"|".join(k): bootstrap(v, rng) for k, v in store.items()}
    chosen = next((n for n in ("G3b", "G3a") if decide(res, n)[0]), None)
    null = []
    if chosen:  # P38: permutation null of the adopted variant's section on the external exam
        t, mats = ext
        prng = np.random.default_rng(1)
        for k in range(PERMUTATIONS):
            gp = Genes(t["g"], t["D"], G=permuted_genes(t["g"], mats[chosen], prng))
            cnt = counts(t["top"], t["V"], gp.table("base"), t["hits"])
            null.append((int(cnt["table"][:, 0].sum()), int(cnt["table"][:, 1].sum())))
            if (k + 1) % 10 == 0:
                print(f"  permutation {k + 1}/{PERMUTATIONS} ({time.time() - t0:.0f} s)", flush=True)
        Tc = Genes(t["g"], t["D"], G=mats[chosen]).table("base")
        obs = counts(t["top"], t["V"], Tc, t["hits"])["table"].sum(axis=0)
        null = {"observed": [int(obs[0]), int(obs[1])], "null": null}
    out = {"stats": stats, "website": web, "neurofibromatosis": nf, "results": res, "chosen": chosen, "null": null}
    (ROOT / "data" / "eval" / "g3_animal_genes.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), "utf-8")
    report(out)


def report(out=None):
    out = out or json.loads((ROOT / "data" / "eval" / "g3_animal_genes.json").read_text("utf-8"))
    R, S = out["results"], out["stats"]
    sets = list(S)
    L = ["# G3. Disease genes without animal models (shared-genes line)\n",
         "Generated by `py -3 -m evaluation.g3_animal_genes`. Protocol and decision rules: `docs/IMPROVEMENT_PLAN.md`, "
         "«Registered 09-11-2026, on the work PC, before any run: G3 …». Variants change only the shared-genes line; "
         "the model is not touched (G3c is a separate run).\n", "## Genes per source (reported before any result)\n",
         "Disease–gene edges; a disease's genes are the union of the four sources. «Animal-only» = the edge's only "
         "source is `ot_assoc:animal_model` (IMPC/MGI mouse models); «expression-only» = its only source left after "
         "removing animal models is `ot_assoc:rna_expression`. Test diseases are counted on each set's own graph "
         "(benchmark: seed 1 graph, diseases of the five seeds).\n",
         "| | " + " | ".join(sets) + " |", "|---|" + "---|" * len(sets)]
    row = lambda lab, f: L.append(f"| {lab} | " + " | ".join(f(S[s]) for s in sets) + " |")
    row("Diseases (with genes)", lambda x: f"{x['diseases']:,} ({x['with genes']:,})")
    for k in S[sets[0]]["edges"]:
        row(f"{k}: edges (diseases)", lambda x, k=k: f"{x['edges'][k]:,} ({x['diseases with the source'][k]:,})")
    row("Median genes per disease", lambda x: f"{x['median genes (diseases with genes)']:.0f}")
    row("Animal-only edges", lambda x: f"{x['animal-only edges']:,} ({100 * x['animal-only edges'] / max(x['all edges'], 1):.1f} %)")
    row("Expression-only edges", lambda x: f"{x['expression-only edges']:,} ({100 * x['expression-only edges'] / max(x['all edges'], 1):.1f} %)")
    row("Diseases with ≥ 50 % animal-only genes", lambda x: f"{x['≥ 50 % animal-only']:,}")
    row("Diseases with ≥ 50 % animal- or expression-only genes", lambda x: f"{x['≥ 50 % animal- or expression-only']:,}")
    row("Diseases left with no gene, G3a", lambda x: f"{x['left with no gene under G3a']:,}")
    row("Diseases left with no gene, G3b", lambda x: f"{x['left with no gene under G3b']:,}")
    w = out["website"]
    n = max(w.get("shown shared-genes candidates", 0), 1)
    L += ["", "Website (pages as deployed): " + "; ".join(f"{k} {v:,} ({100 * v / n:.1f} %)" if k != "shown shared-genes candidates"
                                                        else f"{k} {v:,}" for k, v in w.items()) + ".\n"]
    L += ["## Variants\n", "Hit rate in % (95 % interval resampling diseases; candidates) of the shared-genes candidates each "
          "variant removes, adds and keeps, and of the section before and after (top 25; benchmark zero-shot seeds 1–5, "
          "the nine areas, the external exam).\n",
          "| Variant | Set | Removed | Added | Kept | Section before | Section after |", "|---|---|---|---|---|---|---|"]
    for name in ("G3a", "G3b"):
        for grp in ("benchmark", "areas", "external"):
            L.append(f"| {name} | {grp} | " + " | ".join(fmt(R.get(f"{grp}|{name}|{k}")) for k in
                                                       ("removed", "added", "kept", "before", "after")) + " |")
        ok, rules, notes = decide(R, name)
        L.append(f"| **{name}** | rules | " + "; ".join(f"{k}: {'yes' if v else 'no'}" for k, v in rules.items())
                 + (f" ({'; '.join(notes)})" if notes else "") + f" | → **{'passes' if ok else 'does not pass'}** | | | |")
    chosen = out["chosen"]
    L += ["", f"**Decision (G3b if it passes, else G3a):** {(chosen + ' adopted') if chosen else 'neither passes; the website keeps the rule of P48'}.\n"]
    if out.get("null"):
        obs, null = out["null"]["observed"], out["null"]["null"]
        rates = np.array([h / max(c, 1) for h, c in null])
        L += ["## Permutation null of the adopted variant (P38; reported, not deciding)\n",
              f"External exam, {len(null)} permutations (every disease keeps its number of {chosen} genes, genes redrawn "
              f"in proportion to their degree): observed {100 * obs[0] / max(obs[1], 1):.2f} % ({obs[0]} of {obs[1]}); "
              f"null mean {100 * rates.mean():.2f} %, 99th percentile {100 * np.percentile(rates, 99):.2f} %.\n"]
    nf = out.get("neurofibromatosis")
    if nf:
        L += ["## «neurofibromatosis» (the E1 reply's example)\n",
              f"Genes by source: " + ", ".join(f"{k} {v}" for k, v in nf["genes"].items())
              + f"; animal-only {nf['animal-only']}; non-hub genes counted by the line: all {nf['all (non-hub)']}, "
              f"G3a {nf['G3a genes']}, G3b {nf['G3b genes']}.\n",
              "| Shared-genes candidate | Through | Rests on an animal-only gene | Keeps the line, G3a | G3b |", "|---|---|---|---|---|"]
        L += [f"| {a} | {b} | {'yes' if c else 'no'} | {'yes' if x else 'no'} | {'yes' if y else 'no'} |" for a, b, c, x, y in nf["candidates"]]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "G3_ANIMAL_GENES.md").write_text(text, "utf-8")
    print(text)
    return chosen


if __name__ == "__main__":
    report() if sys.argv[1:] == ["report"] else main()
