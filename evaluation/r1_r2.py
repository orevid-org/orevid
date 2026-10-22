"""R1 (shared-genes candidates from unrelated diseases) and R2 (anti-infectives of the wrong pathogen class).

Protocol and decision rules: docs/IMPROVEMENT_PLAN.md, «Registered 22-10-2026, before any run» (fixed and pushed
before the first run). Same top 25 and sections as evaluation/evidence_ci.py.

    py -3 -m evaluation.r1_r2        # -> data/eval/r1_r2.json, docs/R1_R2.md
"""
import json
import sys
import time
from collections import deque
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from evaluation.evidence_ci import Genes, bench_task, bootstrap, external_task, fmt, node_types_all

ROOT = Path(__file__).resolve().parent.parent
MIN_GENES, HUB = 10, 500
CLASS_LABELS = {"bacterial": "bacterial infectious disease", "fungal": "fungal infectious disease",
                "viral": "viral infectious disease", "parasitic": "parasitic infection"}
DRUG_WORDS = {"bacterial": ("antibacterial", "antibiotic", "antimycobacterial"), "fungal": ("antimycotic", "antifungal"),
              "viral": ("antiviral",), "parasitic": ("antiprotozoal", "anthelmintic", "antimalarial", "ectoparasiticide")}


def split_counts(top, V, T_prod, T_var, hits):
    """Per disease (h, n) of the shared-genes candidates the variant removes and of those it keeps."""
    rem, kept = np.zeros((len(top), 2), int), np.zeros((len(top), 2), int)
    for c, idx in enumerate(top):
        v, tp, tv, h = V[idx, c], T_prod[idx, c], T_var[idx, c], hits[c][idx]
        r, k = ~v & tp & ~tv, ~v & tp & tv
        rem[c] = (h[r].sum(), r.sum())
        kept[c] = (h[k].sum(), k.sum())
    return rem, kept


def removed_from_shown(top, V, T, drop, hits):
    """(h, n) per disease of shown candidates (extensions or shared genes) that `drop` (drugs × D) hides."""
    out = np.zeros((len(top), 2), int)
    pairs = []
    for c, idx in enumerate(top):
        shown = V[idx, c] | T[idx, c]
        r = shown & drop[idx, c]
        out[c] = (hits[c][idx][r].sum(), r.sum())
        pairs += [(c, int(j)) for j in idx[r]]
    return out, pairs


# --- R1 -------------------------------------------------------------------------------------------

def r1_tables(t):
    from evaluation.txgnn import GENES_EXTRA
    g = t["g"]
    G = g.M(GENES_EXTRA).tocsr()
    gn = Genes(g, t["D"], hubs=False)  # the rule before R1b
    prod = gn.table("base")
    small = gn.ng[t["D"]] < MIN_GENES
    r1a = prod & ~small[None, :]
    deg = np.diff(G.tocsc().indptr)
    hub = (node_types_all() == "gene/protein") & (deg >= HUB)
    keep = sp.diags((~hub).astype(np.float32))
    G2 = (keep @ G @ keep).tocsr()
    G2.eliminate_zeros()
    r1b = Genes(g, t["D"], G=G2, hubs=False).table("base")
    return prod, {"R1a": r1a, "R1b": r1b, "R1a+R1b": r1a & r1b}, int(hub.sum())


# --- R2 -------------------------------------------------------------------------------------------

def disease_classes(gc):
    from sources.hierarchy import load
    h = load()
    parents = {}
    for p, c in zip(h.parent.astype(int), h.child.astype(int)):
        parents.setdefault(c, []).append(p)
    names = gc.nodes.node_name.str.lower()
    roots, missing = {}, []
    for cls, lab in CLASS_LABELS.items():
        hit = gc.nodes[(gc.nodes.node_type == "disease") & (names == lab)].node_index.tolist()
        if hit:
            roots[int(hit[0])] = cls
        else:
            missing.append(lab)

    def classes(d):
        seen, q, out = {d}, deque([d]), set()
        while q:
            x = q.popleft()
            if x in roots:
                out.add(roots[x])
            for p in parents.get(x, []):
                if p not in seen:
                    seen.add(p)
                    q.append(p)
        return out
    return classes, missing


def drug_classes(gc):
    d = gc.descriptions
    out = {}
    for i in d.index:
        txt = " ".join(str(d.loc[i, c]) for c in ("atc_1", "atc_2", "atc_3", "atc_4") if c in d and d.loc[i, c] == d.loc[i, c]).lower()
        cl = {k for k, words in DRUG_WORDS.items() if any(w in txt for w in words)}
        if cl:
            out[int(i)] = cl
    return out


def r2_drop(t, dis_cls, drug_cls):
    drugs = t["g"].drugs
    drop = np.zeros((len(drugs), len(t["D"])), bool)
    for c, d in enumerate(t["D"]):
        dc = dis_cls(int(d))
        if len(dc) != 1:
            continue
        (cls,) = dc
        for j, f in enumerate(drugs):
            k = drug_cls.get(int(f))
            if k and cls not in k:
                drop[j, c] = True
    return drop


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.txgnn import AREAS, load_kg, node_types
    from graph import INDICATION, Graph
    t0 = time.time()
    kg = load_kg()
    types = node_types(kg)
    gc = Graph()
    dis_cls, missing = disease_classes(gc)
    drug_cls = drug_classes(gc)
    import pandas as pd
    ot = pd.read_parquet(ROOT / "data" / "ext" / "ot_approvals.parquet")
    ot_conf = set(zip(ot[ot.status == "confirmed"].x_index.astype(int), ot[ot.status == "confirmed"].y_index.astype(int)))
    rng = np.random.default_rng(0)
    store = {}
    r2_pairs = []
    n_hubs = None
    groups = {"benchmark": [("complex_disease", s) for s in (1, 2, 3, 4, 5)], "areas": [(a, 1) for a in AREAS]}
    tasks = [(grp, split, seed) for grp, lst in groups.items() for split, seed in lst] + [("external", None, None)]
    infect = {}
    for grp, split, seed in tasks:
        t = external_task(kg, types) if grp == "external" else bench_task(kg, types, split, seed)
        hits = t["hits"]
        prod, variants, n_hubs = r1_tables(t)
        for name, T in variants.items():
            rem, kept = split_counts(t["top"], t["V"], prod, T, hits)
            store.setdefault((grp, name, "removed"), []).append(tuple(rem.T))
            store.setdefault((grp, name, "kept"), []).append(tuple(kept.T))
        drop = r2_drop(t, dis_cls, drug_cls)
        rem, pairs = removed_from_shown(t["top"], t["V"], prod, drop, hits)
        store.setdefault((grp, "R2", "removed"), []).append(tuple(rem.T))
        infect[grp] = infect.get(grp, 0) + int(sum(len(dis_cls(int(d))) == 1 for d in t["D"]))
        for c, j in pairs:
            f, d = int(t["g"].drugs[j]), int(t["D"][c])
            r2_pairs.append({"set": grp, "split": split, "disease": gc.name(d), "drug": gc.name(f),
                             "hit": bool(hits[c][j]), "class": sorted(dis_cls(d))[0],
                             "drug_classes": sorted(drug_cls.get(f, ())),
                             "approved": bool(d in set(int(x) for x in gc.neighbors(f, INDICATION)) or (f, d) in ot_conf)})
        print(f"  {grp} {split} {seed} ({time.time() - t0:.0f} s)", flush=True)
    res = {"|".join(k): bootstrap(v, rng) for k, v in store.items()}
    out = {"results": res, "hubs": n_hubs, "missing_labels": missing, "infectious_test_diseases": infect,
           "r2_pairs": r2_pairs}
    (ROOT / "data" / "eval" / "r1_r2.json").write_text(json.dumps(out, indent=1), "utf-8")
    report(out)


def p29(R, name):
    b, a = R.get(f"benchmark|{name}|removed"), R.get(f"areas|{name}|removed")
    return bool(b and a and b[0] < 0.05 and a[0] < 0.05), b, a


def report(out=None):
    out = out or json.loads((ROOT / "data" / "eval" / "r1_r2.json").read_text("utf-8"))
    R = out["results"]
    L = ["# R1 and R2: two filters on the candidates shown\n",
         "Generated by `py -3 -m evaluation.r1_r2`. Rules fixed before running: `docs/IMPROVEMENT_PLAN.md`, «Registered "
         "22-10-2026, before any run». Hit rate in % (95 % interval resampling diseases; candidates) among the 25 "
         "highest-ranked candidates per disease.\n",
         f"## R1. Shared-genes candidates (R1a: test diseases with < {MIN_GENES} genes; R1b: hub genes, ≥ {HUB} diseases, "
         f"do not count — {out['hubs']:,} hub genes)\n",
         "| Variant | Set | Removed from the section | Kept in the section |", "|---|---|---|---|"]
    decisions = {}
    for name in ("R1a", "R1b", "R1a+R1b"):
        for grp in ("benchmark", "areas", "external"):
            L.append(f"| {name} | {grp} | {fmt(R.get(f'{grp}|{name}|removed'))} | {fmt(R.get(f'{grp}|{name}|kept'))} |")
        ok29, _, _ = p29(R, name)
        er, ek = R.get(f"external|{name}|removed"), R.get(f"external|{name}|kept")
        ok_ext = bool(er is None or (ek and er[0] <= ek[0]))
        decisions[name] = ok29 and ok_ext
        L.append(f"| **{name}** | rule | P29 (removed < 5 % on benchmark and areas): {'yes' if ok29 else 'no'} | "
                 f"external: removed ≤ kept: {'yes' if ok_ext else 'no'} → **{'adopted' if decisions[name] else 'not adopted'}** |")
    if decisions["R1a"] and decisions["R1b"]:
        chosen = "R1a+R1b" if decisions["R1a+R1b"] else max(("R1a", "R1b"), key=lambda k: R[f"benchmark|{k}|removed"][3])
    else:
        chosen = next((k for k in ("R1a", "R1b") if decisions[k]), None)
    L += ["", f"**R1 decision:** {chosen + ' adopted' if chosen else 'neither variant adopted'}.\n",
          "## R2. Anti-infectives of the wrong pathogen class\n",
          "Infection classes found in the graph: " + ", ".join(k for k, lab in CLASS_LABELS.items() if lab not in out["missing_labels"])
          + "; not found (PrimeKG merges them into «infectious disease»): " + (", ".join(out["missing_labels"]) or "none") + ".",
          "Test diseases with one infection class: " + ", ".join(f"{k} {v}" for k, v in out["infectious_test_diseases"].items()) + ".\n",
          "| Set | Removed from the shown sections |", "|---|---|"]
    for grp in ("benchmark", "areas", "external"):
        L.append(f"| {grp} | {fmt(R.get(f'{grp}|R2|removed'))} |")
    pairs = out["r2_pairs"]
    approved = [p for p in pairs if p["approved"]]
    inf = lambda grp: (R.get(f"{grp}|R2|removed") or (0, 0, 0, 0))
    informative = [g for g in ("benchmark", "areas") if inf(g)[3] >= 20]
    below = all(inf(g)[0] < 0.05 for g in informative) and (len(informative) == 2 or inf("external")[0] < 0.05)
    ok = bool(below and not approved and any(inf(g)[3] for g in ("benchmark", "areas", "external")))
    L += ["", f"Informative sets (≥ 20 removed candidates): {', '.join(informative) or 'none'}. "
              f"Approved indications among the removed pairs: **{len(approved)}**.\n",
          f"**R2 decision:** {'adopted' if ok else 'not adopted'}.\n", "Removed pairs (manual check):\n",
          "| Set | Disease (class) | Drug (classes) | Correct | Approved |", "|---|---|---|---|---|"]
    seen = set()
    for p in pairs:
        k = (p["disease"], p["drug"])
        if k in seen:
            continue
        seen.add(k)
        L.append(f"| {p['set']} | {p['disease']} ({p['class']}) | {p['drug']} ({', '.join(p['drug_classes'])}) | "
                 f"{'yes' if p['hit'] else ''} | {'**yes**' if p['approved'] else ''} |")
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "R1_R2.md").write_text(text, "utf-8")
    print(text[:6000])


if __name__ == "__main__":
    report() if sys.argv[1:] == ["report"] else main()
