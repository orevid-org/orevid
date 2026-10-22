"""N1: do the shared genes the website shows overlap curated drug mechanisms (DrugMechDB) more than chance?

Protocol and outcome rule: docs/IMPROVEMENT_PLAN.md, «Registered 22-10-2026, before any run» (N1), fixed before
R1b was adopted; the primary result uses the rule registered then (≥ 3 shared genes, p < 1e-5, hub genes counted),
the rule with R1b (hub genes left out, P48) is reported as secondary. DrugMechDB is used only for evaluation.

    py -3 -m evaluation.drugmechdb       # -> data/eval/drugmechdb.json, docs/DRUGMECHDB.md
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "drugmechdb" / "indication_paths.json"
DRAWS, MIN_SHARED, MAX_P = 1000, 3, 1e-5


def uniprot_to_node(g):
    h = pd.read_csv(ROOT / "data" / "raw" / "hgnc" / "hgnc_complete_set.txt", sep="\t", low_memory=False,
                    usecols=["entrez_id", "uniprot_ids"])
    gene_node = {str(i): int(n) for n, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type)
                 if t == "gene/protein"}
    out = {}
    for e, u in zip(h.entrez_id, h.uniprot_ids):
        if e == e and u == u:
            node = gene_node.get(str(int(e)))
            if node is not None:
                for x in str(u).split("|"):
                    out.setdefault(x.strip(), node)
    return out


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.evidence_ci import drop_hubs, node_types_all
    from evaluation.twins import twins
    from evaluation.txgnn import GENES_EXTRA
    from graph import INDICATION, Graph
    from sources.ids import xref_to_mondo
    g = Graph()
    paths = json.loads(RAW.read_text("utf-8"))
    drug_node = {f"DB:{i}": int(n) for n, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type) if t == "drug"}
    mondo_node = {}
    for n, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type):
        if t == "disease" and str(i).replace("_", "").isdigit():
            for p in str(i).split("_"):
                mondo_node.setdefault(f"MONDO:{int(p)}", int(n))
    to_mondo = xref_to_mondo()
    up = uniprot_to_node(g)
    tw = twins(g, "disease")
    I = g.M(INDICATION).tocsr()
    rules = {"registered (hub genes counted)": g.M(GENES_EXTRA).tocsr(), "with R1b (hub genes left out)": drop_hubs(g.M(GENES_EXTRA))}
    cov = defaultdict(int)
    pairs = []  # (drug node, disease node, mechanism gene nodes)
    for p in paths:
        cov["paths"] += 1
        f = drug_node.get(p["graph"].get("drugbank") or "")
        mesh = p["graph"].get("disease_mesh") or ""
        ds = {mondo_node[k] for k in to_mondo.get(mesh.replace("MESH:", "MESH:"), ()) if k in mondo_node}
        if f is None:
            cov["drug not in graph"] += 1
            continue
        if len(ds) != 1:
            cov["disease not mapped (or ambiguous)"] += 1
            continue
        d = ds.pop()
        drug_ids = {n["id"] for n in p["nodes"] if n["label"] == "Drug"}
        direct = {l["target"] for l in p["links"] if l["source"] in drug_ids}
        M = {up[n["id"].split(":", 1)[1]] for n in p["nodes"]
             if n["label"] == "Protein" and n["id"] not in direct and n["id"].split(":", 1)[1] in up}
        if not M:
            cov["no protein beyond the direct targets"] += 1
            continue
        pairs.append((f, d, M))
    cov["pairs with drug, disease and mechanism genes"] = len(pairs)
    rng = np.random.default_rng(0)
    res = {}
    for name, G in rules.items():
        deg = np.diff(G.tocsc().indptr)
        N = int(((deg > 0) & (node_types_all()[:len(deg)] == "gene/protein")).sum())
        ng = np.diff(G.indptr)
        rows = []
        for f, d, M in pairs:
            genes_d = set(G.indices[G.indptr[d]:G.indptr[d + 1]])
            same = tw.get(d, {d})
            S = set()
            for m in I.indices[I.indptr[f]:I.indptr[f + 1]]:
                if m in same:
                    continue
                sh = genes_d.intersection(G.indices[G.indptr[m]:G.indptr[m + 1]])
                if len(sh) >= MIN_SHARED and hypergeom.sf(len(sh) - 1, N, ng[m], len(genes_d)) < MAX_P:
                    S |= sh
            if S:
                rows.append((np.array(sorted(genes_d)), len(S), bool(S & M), M))
        obs = float(np.mean([r[2] for r in rows])) if rows else None
        null = np.zeros(DRAWS)
        for k in range(DRAWS):
            null[k] = np.mean([bool(set(rng.choice(gd, s, replace=False).tolist()) & M) for gd, s, _, M in rows]) if rows else 0
        res[name] = {"pairs": len(rows), "observed": obs, "null_mean": float(null.mean()),
                     "null_p99": float(np.percentile(null, 99)), "p": float((null >= obs).mean()) if rows else None,
                     "median_S": float(np.median([r[1] for r in rows])) if rows else None}
        print(name, res[name], flush=True)
    out = {"coverage": dict(cov), "results": res}
    (ROOT / "data" / "eval" / "drugmechdb.json").write_text(json.dumps(out, indent=1), "utf-8")
    report(out)


def report(out=None):
    out = out or json.loads((ROOT / "data" / "eval" / "drugmechdb.json").read_text("utf-8"))
    L = ["# N1. Shared genes against curated drug mechanisms (DrugMechDB)\n",
         "Generated by `py -3 -m evaluation.drugmechdb`. Protocol and outcome rule: `docs/IMPROVEMENT_PLAN.md`, «Registered "
         "22-10-2026, before any run» (N1). DrugMechDB (Gonzalez-Cavazos et al., Sci Data 2023; CC0) is used only to "
         "evaluate; nothing from it enters the model or the website.\n", "## Coverage\n"]
    L += [f"- {k}: {v:,}" for k, v in out["coverage"].items()]
    L += ["", "## Result\n",
          "Pairs where the shared-genes evidence (genes the disease shares with other diseases the drug treats, production "
          "rule) contains at least one protein of the curated mechanism beyond the drug's direct targets, against the same "
          "number of genes drawn at random from the same disease's genes (1,000 draws).\n",
          "| Rule | Pairs with shared-genes evidence | Observed | Null mean | Null 99th percentile | p | Median genes shown |",
          "|---|---|---|---|---|---|---|"]
    for k, v in out["results"].items():
        o = "—" if v["observed"] is None else f"{100 * v['observed']:.1f} %"
        L.append(f"| {k} | {v['pairs']:,} | {o} | {100 * v['null_mean']:.1f} % | {100 * v['null_p99']:.1f} % | "
                 f"{v['p'] if v['p'] is not None else '—'} | {v['median_S']} |")
    prim = next(iter(out["results"].values()))
    above = prim["observed"] is not None and prim["observed"] > prim["null_p99"]
    L += ["", f"**Outcome (registered rule, primary row):** the shared genes "
              f"{'overlap curated mechanisms more than random genes of the same disease' if above else 'do not overlap curated mechanisms more than random genes of the same disease'}"
              f" ({'above' if above else 'not above'} the null's 99th percentile). "
              f"{'The Method page may say so.' if above else 'The P38 wording stays.'}"]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "DRUGMECHDB.md").write_text(text, "utf-8")
    print(text)


if __name__ == "__main__":
    report() if sys.argv[1:] == ["report"] else main()
