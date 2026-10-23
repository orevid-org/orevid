"""M1: does the direction of a protein QTL colocalisation agree with the drugs that treat a disease?
(docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, before any run: M1»; data from sources/coloc.py)

    py -3 -m evaluation.m1_direction      # -> data/eval/m1_direction.json, docs/M1_DIRECTION.md
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
B = 5000
LOWER = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "INVERSE AGONIST", "DEGRADER",
         "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "DISRUPTING AGENT"}
RAISE = {"AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR", "OPENER", "STABILISER", "PARTIAL AGONIST"}


def gene_disease_signs(g, qtl_types):
    """(gene node, disease node) -> +1 (higher level, higher risk) / -1, only when every colocalisation agrees;
    also the Open Targets disease ids behind each pair."""
    from evaluation.twins import twins
    from sources import ids
    c = pd.read_parquet(ROOT / "data" / "ext" / "coloc_qtl.parquet")
    c = c[c.rightStudyType.isin(qtl_types) & (c.trans == False) & c.gene.notna() & (c.disease_ids != "")]  # noqa: E712
    n = g.nodes
    gene_node = {f"NCBI:{i}": int(k) for k, i, t in zip(n.node_index, n.node_id, n.node_type) if t == "gene/protein"}
    key_node = {}
    for k, i, t, src in zip(n.node_index, n.node_id, n.node_type, n.node_source):
        if t != "disease":
            continue
        keys = ids.primekg_keys("disease", i) if src in ("MONDO", "MONDO_grouped") else [str(i)]
        for key in keys:
            key_node.setdefault(key, int(k))
    tw = twins(g, "disease")
    canon = lambda d: min(tw.get(d, {d}))
    mapping, _, _, _ = ids.diseases()
    ens = ids.genes()
    signs, ot_ids = defaultdict(set), defaultdict(set)
    for gene, dids, s in zip(c.gene, c.disease_ids, c.betaRatioSignAverage):
        gk = ens.get(gene, (None,))[0]
        gn = gene_node.get(gk)
        if gn is None:
            continue
        for did in dids.split("|"):
            d = key_node.get(mapping.get(did, ""))
            if d is not None:
                signs[(gn, canon(d))].add(int(np.sign(s)))
                ot_ids[(gn, canon(d))].add(did)
    return {k: next(iter(v)) for k, v in signs.items() if len(v) == 1}, ot_ids, canon, len(signs)


def drug_actions(g):
    """drug node -> {gene node: -1 lowers / +1 raises} (drugs with both kinds on a gene: that gene dropped)."""
    from sources import ids
    m = pd.read_parquet(ROOT / "data" / "raw" / "opentargets" / "drug_mechanism_of_action",
                        columns=["actionType", "chemblIds", "targets"])
    m = m.explode("chemblIds").explode("targets").dropna(subset=["chemblIds", "targets"])
    n = g.nodes
    key_node = {(f"DB:{i}" if s == "DrugBank" else f"CHEMBL:{i}"): int(k)
                for k, i, t, s in zip(n.node_index, n.node_id, n.node_type, n.node_source) if t == "drug"}
    gene_node = {f"NCBI:{i}": int(k) for k, i, t in zip(n.node_index, n.node_id, n.node_type) if t == "gene/protein"}
    mol, ens = ids.drugs(), ids.genes()
    acts = defaultdict(lambda: defaultdict(set))
    for a, cid, tgt in zip(m.actionType, m.chemblIds, m.targets):
        s = -1 if a in LOWER else 1 if a in RAISE else 0
        f = key_node.get(mol.get(cid, (None,))[0])
        gn = gene_node.get(ens.get(tgt, (None,))[0])
        if s and f is not None and gn is not None:
            acts[f][gn].add(s)
    return {f: {gn: next(iter(v)) for gn, v in d.items() if len(v) == 1} for f, d in acts.items()}


def ratio(rows, pos_key, rng, swap=False):
    """rate(agree) / rate(oppose) (swap: the reverse), 95 % interval resampling diseases."""
    df = pd.DataFrame(rows)
    if df.empty:
        return None
    a, o = ("oppose", "agree") if swap else ("agree", "oppose")
    per = df.groupby("disease").apply(lambda x: pd.Series({
        "an": (x.dir == a).sum(), "ah": ((x.dir == a) & x[pos_key]).sum(),
        "on": (x.dir == o).sum(), "oh": ((x.dir == o) & x[pos_key]).sum()}), include_groups=False)
    v = per.to_numpy(float)
    idx = rng.integers(0, len(v), (B, len(v)))
    s = v[idx].sum(axis=1)
    r = (s[:, 1] / np.maximum(s[:, 0], 1)) / np.maximum(s[:, 3] / np.maximum(s[:, 2], 1), 1e-12)
    tot = v.sum(axis=0)
    point = (tot[1] / max(tot[0], 1)) / max(tot[3] / max(tot[2], 1), 1e-12)
    return {"ratio": float(point), "low": float(np.percentile(r, 2.5)), "high": float(np.percentile(r, 97.5)),
            "pairs_agree": int(tot[0]), "pos_agree": int(tot[1]), "pairs_oppose": int(tot[2]), "pos_oppose": int(tot[3]),
            "diseases": int(len(v))}


def run(g, qtl_types, acts, ind, contra, pkg_ind, rng, areas):
    signs, ot_ids, canon, n_pairs = gene_disease_signs(g, qtl_types)
    by_dis = defaultdict(dict)
    for (gn, d), s in signs.items():
        by_dis[d][gn] = s
    rows = []
    for f, tg in acts.items():
        for d, gs in by_dis.items():
            dirs = {gs[gn] * -tg[gn] for gn in tg.keys() & gs.keys()}  # +1: agrees
            if len(dirs) != 1:
                continue
            area = sorted({a for gn in tg.keys() & gs.keys() for did in ot_ids[(gn, d)] for a in areas.get(did, ())})
            rows.append({"drug": f, "disease": d, "dir": "agree" if dirs == {1} else "oppose",
                         "ind": (f, d) in ind, "contra": (f, d) in contra, "pkg": (f, d) in pkg_ind,
                         "ot_only": (f, d) in ind and (f, d) not in pkg_ind, "area": area[0] if area else ""})
    cover = {"gene_disease_pairs_with_coloc": n_pairs, "gene_disease_pairs_consistent_sign": len(signs),
             "genes": len({k[0] for k in signs}), "diseases": len(by_dis),
             "drugs_with_directed_pairs": len({r["drug"] for r in rows}), "directed_pairs": len(rows),
             "indication_pairs": sum(r["ind"] for r in rows), "contraindication_pairs": sum(r["contra"] for r in rows)}
    res = {"coverage": cover, "a": ratio(rows, "ind", rng), "b": ratio(rows, "contra", rng, swap=True),
           "a_ot_only": ratio([r for r in rows if not r["pkg"]], "ot_only", rng)}
    per_area = {}
    for a in sorted({r["area"] for r in rows} - {""}):
        sub = [r for r in rows if r["area"] == a]
        if sum(r["ind"] for r in sub) >= 10:
            per_area[a] = ratio(sub, "ind", rng)
    res["a_by_area"] = per_area
    return res


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.twins import twins
    from graph import CONTRA, INDICATION, Graph
    g = Graph()
    tw = twins(g, "disease")
    canon = lambda d: min(tw.get(d, {d}))
    I, C = g.M(INDICATION).tocsr(), g.M(CONTRA).tocsr()
    drugs = set(int(f) for f in g.drugs)
    edges = lambda M: {(int(f), canon(int(d))) for f in drugs for d in M.indices[M.indptr[f]:M.indptr[f + 1]]
                       if g.types[d] == "disease"}
    pkg_ind, contra = edges(I), edges(C)
    ot = pd.read_parquet(ROOT / "data" / "ext" / "ot_approvals.parquet")
    ot = ot[ot.status == "confirmed"]
    ind = pkg_ind | {(int(f), canon(int(d))) for f, d in zip(ot.x_index, ot.y_index)}
    dis = pd.read_parquet(ROOT / "data" / "raw" / "opentargets" / "disease", columns=["id", "therapeuticAreas", "name"])
    tname = dict(zip(dis.id, dis.name))
    areas = {i: [tname.get(a, a) for a in (t if t is not None else [])] for i, t in zip(dis.id, dis.therapeuticAreas)}
    acts = drug_actions(g)
    rng = np.random.default_rng(0)
    out = {name: run(g, types, acts, ind, contra, pkg_ind, rng, areas)
           for name, types in (("pQTL", ["pqtl"]), ("eQTL", ["eqtl"]), ("pQTL+eQTL", ["pqtl", "eqtl"]))}
    (ROOT / "data" / "eval" / "m1_direction.json").write_text(json.dumps(out, indent=1), "utf-8")
    report(out)


def fmt(r):
    if not r:
        return "–"
    return (f"{r['ratio']:.2f} ({r['low']:.2f}–{r['high']:.2f}); agree {r['pos_agree']}/{r['pairs_agree']:,}, "
            f"oppose {r['pos_oppose']}/{r['pairs_oppose']:,}; {r['diseases']} diseases")


def report(out=None):
    out = out or json.loads((ROOT / "data" / "eval" / "m1_direction.json").read_text("utf-8"))
    L = ["# M1. Direction from protein QTL colocalisation\n",
         "Generated by `py -3 -m evaluation.m1_direction`. Protocol and decision rules: `docs/IMPROVEMENT_PLAN.md`, "
         "«Registered 23-10-2026, before any run: M1». Data: Open Targets 26.06 colocalisations (h4 ≥ 0.8, consistent "
         "effect-ratio sign, cis QTL only; `sources/coloc.py`). «Agrees»: the drug lowers a protein whose higher level "
         "goes with higher risk, or raises one whose higher level goes with lower risk. Ratios with 95 % intervals "
         "resampling diseases; counts are positives / pairs.\n"]
    for name, r in out.items():
        cv = r["coverage"]
        L += [f"## {name}{' (deciding)' if name == 'pQTL' else ' (reported)'}\n",
              "Coverage: " + ", ".join(f"{k.replace('_', ' ')} {v:,}" for k, v in cv.items()) + ".\n",
              "| Test | Ratio (95 % CI); counts |", "|---|---|",
              f"| (a) indication rate, agree ÷ oppose | {fmt(r['a'])} |",
              f"| (b) contraindication rate, oppose ÷ agree | {fmt(r['b'])} |",
              f"| (a) on Open Targets approvals absent from PrimeKG | {fmt(r['a_ot_only'])} |"]
        for a, v in r["a_by_area"].items():
            L.append(f"| (a) {a} | {fmt(v)} |")
        L.append("")
    p = out["pQTL"]
    dec_a = bool(p["a"] and p["a"]["low"] > 1 and p["a"]["pos_agree"] + p["a"]["pos_oppose"] >= 30)
    dec_b = bool(p["b"] and p["b"]["low"] > 1 and p["b"]["pos_agree"] + p["b"]["pos_oppose"] >= 30)
    L += [f"**Decision (pQTL):** (a) direction {'validated' if dec_a else 'not validated'}; (b) risk signal "
          f"{'validated' if dec_b else 'not validated'}. Nothing changes on the website from M1; a validated test only "
          "allows a separate registration (display line or model feature; risk warning with the D1 release check).\n"]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "M1_DIRECTION.md").write_text(text, "utf-8")
    print(text)


if __name__ == "__main__":
    report() if sys.argv[1:] == ["report"] else main()
