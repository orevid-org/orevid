"""D1: drugs that may worsen a disease (docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, before any run: D1»).

Signals for a disease d and an approved drug f, with d's own contraindications (and its twins') hidden:
  S1  the P15 rule: contraindicated in more ontology-neighbour (≤ 2 steps) and gene-similar (top 30) diseases than
      indicated in them, or an immune-checkpoint blocker in an autoimmune disease (as web/build.py);
  S2  an adverse event of f in a current label (OnSIDES, sources/onsides.py) maps to d or a direct broader term of d;
  S3  both.
Positives: the hidden PrimeKG contraindications of d. Candidates: approved drugs not indicated for d. Baseline:
«contraindicated somewhere», drugs ranked by the number of other diseases they are contraindicated for, as many flags
per disease as the signal. Diseases with ≥ 3 contraindications; development half = even node id, confirmation = odd.

    py -3 -m evaluation.d1_worsen      # -> data/eval/d1_worsen.json, docs/D1_WORSEN.md
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
B = 5000
ANGELMAN = 29997


def main(source="onsides"):
    """source: onsides (D1) or offsides (L6, the same test with FAERS-derived signals)."""
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.direction import Direction
    from evaluation.twins import twins
    from graph import CONTRA, INDICATION
    from predictor import drug_status
    from sources.hierarchy import load as load_hierarchy
    from evaluation.direction import checkpoint_rule
    D = Direction()
    g = D.g
    status = drug_status(g)
    approved = np.array(sorted(f for f in g.drugs if status.get(int(f)) == "approved"), dtype=int)
    pos_of = {int(f): j for j, f in enumerate(approved)}
    tw = twins(g, "disease")
    C = D.contra.tocsr()
    A = D.approved.tocsr()
    blocked, autoimmune = checkpoint_rule(g)
    blocked = {int(x) for x in blocked}
    # contraindicated-somewhere counts per drug (rows of C restricted to diseases)
    is_dis = g.types == "disease"
    contra_count = np.asarray(C[approved][:, np.flatnonzero(is_dis)].sum(axis=1)).ravel()
    h = load_hierarchy()
    parents = defaultdict(list)
    for a_, b_ in zip(h.parent.astype(int), h.child.astype(int)):
        parents[b_].append(a_)
    ons = pd.read_parquet(ROOT / "data" / "ext" / f"{source}.parquet")
    ade = defaultdict(set)  # disease -> drugs with a label adverse event naming it
    for f, d in zip(ons.drug.astype(int), ons.disease.astype(int)):
        ade[d].add(f)
    test = [int(d) for d in np.flatnonzero(is_dis & (np.diff(C.indptr) >= 3)) if g.origin[d] == "PrimeKG"]
    rng = np.random.default_rng(0)
    rows = []
    angel = {}
    for d in test:
        same = tw.get(d, {d})
        pos = {int(x) for t in same for x in C.indices[C.indptr[t]:C.indptr[t + 1]]} & set(pos_of)
        if len(pos) < 3:
            continue
        ind_d = {int(x) for t in same for x in A.indices[A.indptr[t]:A.indptr[t + 1]]}
        cand = np.array([f for f in approved if int(f) not in ind_d], dtype=int)
        if not len(cand):
            continue
        o = np.setdiff1d(D.O2[d].indices, list(same))
        sim = (D.G @ D.G[d].T).toarray().ravel()
        sim[list(same)] = 0
        s = np.argsort(-sim)[:30]
        s = s[sim[s] > 0]
        nc = np.asarray(C[cand][:, o].sum(axis=1)).ravel() + np.asarray(C[cand][:, s].sum(axis=1)).ravel()
        ni = np.asarray(A[cand][:, o].sum(axis=1)).ravel() + np.asarray(A[cand][:, s].sum(axis=1)).ravel()
        s1 = (nc > ni) | (np.isin(cand, list(blocked)) & (d in autoimmune))
        ade_d = set().union(ade.get(d, set()), *(ade.get(p, set()) for p in parents.get(d, ())))
        s2 = np.isin(cand, list(ade_d))
        hits = np.isin(cand, list(pos))
        # baseline score: contraindicated somewhere, without d and its twins
        own = np.asarray(C[cand][:, list(same)].sum(axis=1)).ravel()
        base = contra_count[[pos_of[int(f)] for f in cand]] - own + rng.random(len(cand)) * 1e-6
        order = np.argsort(-base)
        rec = {"disease": d, "half": "development" if d % 2 == 0 else "confirmation", "candidates": len(cand), "positives": int(hits.sum())}
        for name, flag in (("S1", s1), ("S2", s2), ("S3", s1 & s2)):
            k = int(flag.sum())
            rec[f"{name}_n"] = k
            rec[f"{name}_h"] = int(hits[flag].sum())
            rec[f"{name}_base_h"] = int(hits[order[:k]].sum()) if k else 0
        rows.append(rec)
        if d == ANGELMAN:
            angel = {n: sorted(g.name(int(f)) for f in cand[fl]) for n, fl in (("S1", s1), ("S2", s2))}
            angel["hidden contraindications"] = sorted(g.name(f) for f in pos)
    df = pd.DataFrame(rows)
    res = {}
    for half, sub in df.groupby("half"):
        for name in ("S1", "S2", "S3"):
            n, h, b = sub[f"{name}_n"].to_numpy(), sub[f"{name}_h"].to_numpy(), sub[f"{name}_base_h"].to_numpy()
            idx = rng.integers(0, len(sub), (B, len(sub)))
            ns = n[idx].sum(axis=1)
            diff = (h[idx].sum(axis=1) - b[idx].sum(axis=1)) / np.maximum(ns, 1)
            res[f"{half}|{name}"] = {"diseases": int((n > 0).sum()), "flags": int(n.sum()), "precision": float(h.sum() / max(n.sum(), 1)),
                                     "baseline": float(b.sum() / max(n.sum(), 1)),
                                     "diff": [float((h.sum() - b.sum()) / max(n.sum(), 1)), float(np.percentile(diff, 2.5)), float(np.percentile(diff, 97.5))]}
    out = {"results": res, "test_diseases": len(df), "angelman": angel,
           "onsides_pairs": int(len(ons)), "base_rate": float(df.positives.sum() / df.candidates.sum())}
    out["source"] = source
    (ROOT / "data" / "eval" / ("d1_worsen.json" if source == "onsides" else f"d1_worsen_{source}.json")).write_text(json.dumps(out, indent=1), "utf-8")
    report(out)


def report(out=None):
    out = out or json.loads((ROOT / "data" / "eval" / "d1_worsen.json").read_text("utf-8"))
    R = out["results"]
    L = ["# D1. Drugs that may worsen a disease\n",
         "Generated by `py -3 -m evaluation.d1_worsen`. Protocol and decision rule: `docs/IMPROVEMENT_PLAN.md`, «Registered "
         "23-10-2026, before any run: D1». Each disease's own contraindications (and its twins') are hidden; precision = "
         "share of flagged drugs that are hidden PrimeKG contraindications; baseline = «contraindicated somewhere» with "
         "the same number of flags; 95 % interval resampling diseases.\n",
         f"Test diseases (≥ 3 contraindications): {out['test_diseases']:,}; base rate (contraindications among candidates): "
         f"{100 * out['base_rate']:.2f} %; OnSIDES drug–disease pairs: {out['onsides_pairs']:,}.\n",
         "| Half | Signal | Diseases flagged | Flags | Precision | Baseline | Difference (95 % CI) |", "|---|---|---|---|---|---|---|"]
    names = {"S1": "S1 P15 rule", "S2": "S2 label adverse event (OnSIDES)", "S3": "S3 both"}
    dec = {}
    for half in ("development", "confirmation"):
        for s in ("S1", "S2", "S3"):
            v = R.get(f"{half}|{s}")
            if not v:
                continue
            d = v["diff"]
            L.append(f"| {half} | {names[s]} | {v['diseases']:,} | {v['flags']:,} | {100 * v['precision']:.1f} % | "
                     f"{100 * v['baseline']:.1f} % | {100 * d[0]:+.1f} ({100 * d[1]:+.1f} to {100 * d[2]:+.1f}) |")
            if half == "confirmation":
                dec[s] = bool(d[1] > 0 and v["flags"] >= 200)
    L += ["", "**Decision (confirmation half):** " + "; ".join(f"{names[s]}: {'passes' if ok else 'does not pass'}" for s, ok in dec.items())
          + ". A signal that passes still needs the manual check of 50 flagged pairs against the source pages before release.\n",
          "Angelman syndrome (review of 22-10; reported, not deciding):\n"]
    for k, v in out.get("angelman", {}).items():
        L.append(f"- {k}: {', '.join(v[:25]) if v else 'none'}{' …' if len(v) > 25 else ''}")
    text = "\n".join(L) + "\n"
    src = out.get("source", "onsides")
    if src != "onsides":
        text = text.replace("OnSIDES", src.upper()).replace("# D1. Drugs that may worsen a disease", f"# L6. D1 with {src.upper()} (FAERS-derived)")
    (ROOT / "docs" / ("D1_WORSEN.md" if src == "onsides" else f"L6_{src.upper()}.md")).write_text(text, "utf-8")
    print(text)


if __name__ == "__main__":
    report() if sys.argv[1:] == ["report"] else main(sys.argv[1] if sys.argv[1:] else "onsides")
