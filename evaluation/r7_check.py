"""Validation of R7 (docs/IMPROVEMENT_PLAN.md, «Registered 08-11-2026: R7»): no direct extension from a *feature
parent* (a MONDO parent of a syndrome with ≥ 20 HPO phenotypes that shares no significant genes with it, p ≥ 1e-3).
Criteria fixed before running: the dropped extensions must hit below 5 % on the benchmark (zero-shot seeds 1-5) and on
the nine areas (top 25 per disease, as evaluation/evidence_ci.py); external exam reported; manual check of 20 random
dropped extensions on the website; pages that change.

    py -3 -m evaluation.r7_check     # -> docs/R7_FEATURE_PARENTS.md
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parent.parent
MIN_PHENOTYPES, MAX_P = 20, 1e-3
EXAMPLES = [("neurofibromatosis", "hypertension", "wrong"), ("Barth syndrome", "syndromic dyslipidemia", "wrong"),
            ("CHARGE syndrome", "congenital hypogonadotropic hypogonadism", "wrong"),
            ("22q11.2 deletion syndrome", "heart disease", "wrong"), ("anorexia nervosa", "anxiety disorder", "wrong"),
            ("Chediak-Higashi syndrome", "constitutional neutropenia", "wrong"), ("breast carcinoma", "carcinoma", "right"),
            ("Ewing sarcoma", "sarcoma", "right"), ("atopic eczema", "dermatitis", "right"),
            ("Bethlem myopathy", "muscular dystrophy", "right"), ("brain ischemia", "cerebrovascular disorder", "right"),
            ("acral lentiginous melanoma (disease)", "melanoma", "right"), ("hemoglobin C disease", "hemoglobinopathy", "right"),
            ("interstitial lung disease", "lung disease", "right")]


def feature_parents(G):
    """{(parent, child)} pairs of the MONDO hierarchy where the parent is only a feature of the child."""
    from graph import GENES_ALL
    from sources.hierarchy import load
    genes = G.M(GENES_ALL).tocsr()
    n_genes = int((np.diff(genes.tocsc().indptr) > 0).sum())
    phen = G.M("disease_phenotype_positive").tocsr()
    h = load()

    def stats(p, c):
        gc, gp = set(genes.indices[genes.indptr[c]:genes.indptr[c + 1]]), set(genes.indices[genes.indptr[p]:genes.indptr[p + 1]])
        sh = len(gc & gp)
        pv = float(hypergeom.sf(sh - 1, n_genes, len(gp), len(gc))) if gc and gp else 1.0
        return int(phen.indptr[c + 1] - phen.indptr[c]), sh, pv

    flagged = set()
    for p, c in zip(h.parent.astype(int), h.child.astype(int)):
        nph, _, pv = stats(p, c)
        if nph >= MIN_PHENOTYPES and pv >= MAX_P:
            flagged.add((p, c))
    return flagged, stats


def task_rates(t, flagged):
    """Hits and counts of direct extensions in the top 25, dropped by R7 or kept."""
    g = t["g"]
    ind, dd = g.M("indication").tocsr(), g.M("disease_disease").tocsr()
    out = {"dropped": [0, 0], "kept": [0, 0], "no_neighbour_found": 0}
    for c, d in enumerate(t["D"]):
        d = int(d)
        nb = set(dd.indices[dd.indptr[d]:dd.indptr[d + 1]])
        for j in t["top"][c]:
            if not t["V"][j, c]:
                continue
            f = int(g.drugs[j])
            via = nb & set(ind.indices[ind.indptr[f]:ind.indptr[f + 1]])
            if not via:
                out["no_neighbour_found"] += 1
                continue
            k = "dropped" if all((m, d) in flagged for m in via) else "kept"
            out[k][0] += int(t["hits"][c][j])
            out[k][1] += 1
    return out


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph, INDICATION
    from evaluation.evidence_ci import bench_task, external_task
    from evaluation.txgnn import AREAS, load_kg, node_types
    G = Graph()
    flagged, stats = feature_parents(G)
    print(f"feature-parent pairs: {len(flagged)}", flush=True)
    kg = load_kg()
    types = node_types(kg)
    groups = {"benchmark": [("complex_disease", s) for s in (1, 2, 3, 4, 5)], "areas": [(a, 1) for a in AREAS]}
    res = {}
    for name, tasks in groups.items():
        tot = {"dropped": [0, 0], "kept": [0, 0], "no_neighbour_found": 0}
        for split, seed in tasks:
            r = task_rates(bench_task(kg, types, split, seed), flagged)
            for k in ("dropped", "kept"):
                tot[k] = [tot[k][0] + r[k][0], tot[k][1] + r[k][1]]
            tot["no_neighbour_found"] += r["no_neighbour_found"]
            print(f"  {name} {split} {seed}: {r}", flush=True)
        res[name] = tot
    res["external"] = task_rates(external_task(kg, types), flagged)
    print(res, flush=True)

    # website: direct extensions on the built pages whose every neighbour with the drug is a feature parent
    dest = ROOT / "web" / "public" / "data" / "d"
    ind, dd = G.M(INDICATION).tocsr(), G.M("disease_disease").tocsr()
    by_name = {}
    for x in np.flatnonzero(G.types == "drug"):
        by_name.setdefault(G.name(int(x)), []).append(int(x))
    dropped_rows, pages_changed = [], set()
    for fj in sorted(dest.glob("*.json")):
        for k, v in json.loads(fj.read_text("utf-8")).items():
            d = int(k)
            nb = set(dd.indices[dd.indptr[d]:dd.indptr[d + 1]])
            for c in v.get("extensions") or []:
                via = {m for f in by_name.get(c["name"], []) for m in ind.indices[ind.indptr[f]:ind.indptr[f + 1]]} & nb
                if via and all((m, d) in flagged for m in via):
                    dropped_rows.append((v["name"], c["name"], ", ".join(sorted(G.name(m) for m in via))))
                    pages_changed.add(d)
    rng = np.random.default_rng(0)
    sample = [dropped_rows[i] for i in rng.choice(len(dropped_rows), min(20, len(dropped_rows)), replace=False)]

    rate = lambda x: f"{100 * x[0] / x[1]:.1f} % ({x[0]} of {x[1]})" if x[1] else "– (0)"
    names = {G.name(int(x)): int(x) for x in np.flatnonzero(G.types == "disease")}
    L = ["# R7. No direct extension from a feature parent\n",
         "Generated by `py -3 -m evaluation.r7_check`. Protocol: `docs/IMPROVEMENT_PLAN.md`, «Registered 08-11-2026: R7». "
         f"Feature parent: a MONDO parent of a disease with ≥ {MIN_PHENOTYPES} HPO phenotypes that shares no significant "
         f"genes with it (hypergeometric p ≥ {MAX_P:g}). Feature-parent pairs: {len(flagged):,}.\n",
         "## The 14 pairs the thresholds were chosen on (no hit rate looked at)\n",
         "| Child | Parent | Judged | Phenotypes of the child | Shared genes | p | Feature parent? |", "|---|---|---|---|---|---|---|"]
    for c, p, lab in EXAMPLES:
        if c in names and p in names:
            nph, sh, pv = stats(names[p], names[c])
            L.append(f"| {c} | {p} | {lab} | {nph} | {sh} | {pv:.1e} | {'yes' if (names[p], names[c]) in flagged else 'no'} |")
    L += ["", "## Hit rate of direct extensions in the top 25 (criterion: dropped < 5 % on benchmark and areas)\n",
          "| Test | Dropped by R7 | Kept |", "|---|---|---|"]
    for name in ("benchmark", "areas", "external"):
        L.append(f"| {name} | {rate(res[name]['dropped'])} | {rate(res[name]['kept'])} |")
    ok = all(res[n]["dropped"][1] == 0 or res[n]["dropped"][0] / res[n]["dropped"][1] < 0.05 for n in ("benchmark", "areas"))
    L += ["", f"**Hit-rate criterion: {'met' if ok else 'NOT met'}.** Extensions whose neighbour could not be found in "
          f"the training graph (left as kept): benchmark {res['benchmark']['no_neighbour_found']}, areas "
          f"{res['areas']['no_neighbour_found']}, external {res['external']['no_neighbour_found']}.\n",
          f"## Website\n", f"Direct extensions dropped: {len(dropped_rows)} on {len(pages_changed)} pages.\n",
          "Random sample of 20 for the manual check:\n", "| Disease | Drug | Feature parent(s) it came from |",
          "|---|---|---|"]
    L += [f"| {a} | {b} | {c} |" for a, b, c in sample]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "R7_FEATURE_PARENTS.md").write_text(text, "utf-8")
    print(text)


if __name__ == "__main__":
    main()
