"""Temporal validation (docs/IMPROVEMENT_PLAN.md, I1; decision rules fixed there before running).

Trained only on what was known before a cutoff year C, does OREVID rank the drug–disease pairs that reached
phase 3 or approval after C above the popularity baseline and above the other drugs?

Cutoff graph at C:
  1. PrimeKG drug–disease edges (indication, contraindication, off-label use) whose pair has an Open Targets
     clinical record with first trial year ≥ C are removed (`clinical` rows hold the first year of the pair).
  2. Open Targets gene–disease associations whose first evidence year is ≥ C are removed (undated ones stay).
  3. Undated PrimeKG edges stay; the share of PrimeKG indications without an Open Targets date is reported as
     the residual leakage.
Positives (primary): pairs with highest phase 3 or approved, first trial year ≥ C, between a PrimeKG drug and
a PrimeKG disease, that are not a drug–disease edge of the cutoff graph (interpretation of «no PrimeKG edge»: an
edge left in the cutoff graph; the strict reading, no PrimeKG edge at all, is reported as a variant).
Secondary: confirmed approvals only (sources/ot_approvals.py) and phase ≥ 2.
Candidates: every PrimeKG drug except those with a drug–disease edge in the cutoff graph or an Open Targets
clinical record (any phase) before C or without a year.

    py -3 -m evaluation.temporal run 2012 [<PRODUCTION>|extra_no_ot]   # -> data/eval/temporal_C2012_<fs>.npz
    py -3 -m evaluation.temporal report                                # -> docs/TEMPORAL.md
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from evaluation.txgnn import PRODUCTION

ROOT = Path(__file__).resolve().parent.parent
CUTOFFS = {"C1 (development)": 2012, "C2 (confirmation)": 2015}
TOP, B = 25, 10_000
PHASE = {"approved": 4, "phase 3": 3, "phase 2": 2, "phase 1": 1}


def clinical():
    a = pd.read_parquet(ROOT / "data" / "prep" / "attributes.parquet", columns=["relation", "display_relation", "x_index",
                                                                                 "y_index", "year"])
    c = a[a.relation == "clinical"].rename(columns={"x_index": "drug", "y_index": "disease"})
    c["phase"] = c.display_relation.map(PHASE).fillna(0).astype(int)
    return c[["drug", "disease", "phase", "year"]], a


def cutoff_mask(kg, clin, C):
    """kg rows kept in the cutoff training graph (PrimeKG drug–disease edges first tried at C or later removed)."""
    from evaluation.txgnn import DD
    first = {(int(d), int(m)): y for d, m, y in zip(clin.drug, clin.disease, clin.year) if y == y}
    dd = kg.relation.isin(DD).to_numpy()
    years = np.array([first.get((int(x), int(y)), first.get((int(y), int(x)), np.nan)) if r else np.nan
                      for x, y, r in zip(kg.x_index.to_numpy(), kg.y_index.to_numpy(), dd)])
    keep = ~(dd & (years >= C))
    ind = (kg.relation == "indication").to_numpy()
    leak = float(np.isnan(years[ind]).mean())
    return keep, leak, int((dd & (years >= C)).sum())


def cutoff_graph(kg, types, keep, attrs, C, extra=True):
    """TxGraph on the kept edges, with the Open Targets associations first seen at C or later removed."""
    from evaluation.txgnn import TxGraph
    g = TxGraph(kg, keep, types, extra=extra)
    if extra:
        ot = attrs[(attrs.relation == "ot_assoc") & (attrs.year >= C)]
        for rel in list(g._m):
            if rel.startswith("ot_assoc:"):
                s = ot[ot.display_relation == rel.split(":", 1)[1]]
                x, y = s.x_index.to_numpy(), s.y_index.to_numpy()
                R = sp.coo_matrix((np.ones(2 * len(s)), (np.r_[x, y], np.r_[y, x])), shape=(g.N, g.N)).tocsr()
                R.data[:] = 1
                m = g._m[rel]
                g._m[rel] = (m - m.multiply(R)).tocsr()
                g._m[rel].eliminate_zeros()
    return g


def run(C, feature_set=None, rounds=5, n_train=1500, seed=1):
    feature_set = feature_set or PRODUCTION
    from evaluation.evidence import _transform
    from evaluation.evidence_ci import Genes
    from evaluation.external_exam import N_PKG
    from evaluation.txgnn import (DD, FEATURE_SETS, INDICATION, drug_degrees, fit_model, fit_ranker, load_kg,
                                  masked_rows, near_duplicates, node_types, score_paths)
    t0 = time.time()
    kg = load_kg()
    types = node_types(kg)
    clin, attrs = clinical()
    keep, leak, removed = cutoff_mask(kg, clin, C)
    g = cutoff_graph(kg, types, keep, attrs, C)
    paths = FEATURE_SETS[feature_set]
    drugs = g.drugs
    drug_pos = np.full(g.N, -1)
    drug_pos[drugs] = np.arange(len(drugs))
    x, y = kg.x_index.to_numpy(), kg.y_index.to_numpy()
    ind = (kg.relation == "indication").to_numpy()
    known = ind & keep
    # --- training (as evaluate(), with every disease known at C; no disease is hidden) ---
    rng = np.random.default_rng(seed)
    train_diseases = pd.unique(y[known])
    train_diseases = np.random.default_rng(seed).choice(train_diseases, min(n_train, len(train_diseases)), replace=False)
    degree = np.diff(g.M(INDICATION).indptr)[drugs].astype(np.float64)
    gf = drug_degrees(g)
    rows_X, rows_y = [], []
    for i in range(0, len(train_diseases), 256):
        block = train_diseases[i:i + 256]
        F = score_paths(g, block, paths)
        for c, d in enumerate(block):
            lab = np.zeros(len(drugs), bool)
            lab[drug_pos[x[known & (y == d)]]] = True
            pop = degree.copy()
            pop[drug_pos[x[known & (y == d)]]] -= 1
            sel = lab | (rng.random(len(lab)) < 0.1)
            rows_X.append(np.column_stack([F[sel, c, :], pop[sel], gf[sel]]))
            rows_y.append(lab[sel])
    X, Y = np.concatenate(rows_X), np.concatenate(rows_y)
    Xm, Ym, Gm = masked_rows(g, np.array([], int), known, x, y, drug_pos, paths, rounds, seed)
    rank = fit_ranker(np.concatenate([X, Xm]), np.concatenate([Y, Ym]), np.r_[[len(v) for v in rows_y], Gm], seed)
    nc = len(paths)
    logit, _ = fit_model(X[:, :nc], Y)
    print(f"  trained at C={C} ({time.time() - t0:.0f} s)", flush=True)
    # --- positives and candidates ---
    is_drug = np.zeros(g.N, bool)
    is_drug[drugs] = True
    pkg_disease = (np.arange(g.N) < N_PKG) & (np.r_[types, np.full(g.N - len(types), "")] == "disease")
    c_ = clin[is_drug[clin.drug.to_numpy()] & pkg_disease[clin.disease.to_numpy()]]
    dd = kg.relation.isin(DD).to_numpy()
    edge_cut = set(zip(x[dd & keep], y[dd & keep])) | set(zip(y[dd & keep], x[dd & keep]))
    edge_all = set(zip(x[dd], y[dd])) | set(zip(y[dd], x[dd]))
    before = c_[(c_.year < C) | c_.year.isna()]
    after = c_[(c_.year >= C).to_numpy() & np.array([(d, m) not in edge_cut for d, m in zip(c_.drug, c_.disease)], bool)]
    conf = pd.read_parquet(ROOT / "data" / "ext" / "ot_approvals.parquet")
    confirmed = set(zip(conf.x_index[conf.status == "confirmed"], conf.y_index[conf.status == "confirmed"]))
    D = np.array(sorted(set(after[after.phase >= 3].disease)))
    col = {d: i for i, d in enumerate(D)}
    cand = np.ones((len(drugs), len(D)), bool)
    for d, m in edge_cut:
        if m in col and drug_pos[d] >= 0:
            cand[drug_pos[d], col[m]] = False
    for d, m in zip(before.drug, before.disease):
        if m in col:
            cand[drug_pos[d], col[m]] = False

    def positives(frame):
        out = [[] for _ in D]
        for d, m in zip(frame.drug, frame.disease):
            if m in col and cand[drug_pos[d], col[m]]:
                out[col[m]].append(drug_pos[d])
        return [np.array(sorted(set(p)), int) for p in out]

    pos = {"primary": positives(after[after.phase >= 3]),
           "strict (no PrimeKG edge at all)": positives(after[(after.phase >= 3).to_numpy() & np.array(
               [(d, m) not in edge_all for d, m in zip(after.drug, after.disease)], bool)]),
           "confirmed approvals": positives(after[(after.phase == 4).to_numpy() & np.array(
               [(d, m) in confirmed for d, m in zip(after.drug, after.disease)], bool)]),
           "phase 2 or later": positives(after[after.phase >= 2])}
    # test diseases: at least one primary positive left among the candidates
    keep_d = np.array([len(v) > 0 for v in pos["primary"]])
    D, cand = D[keep_d], cand[:, keep_d]
    pos = {k: [v[i] for i in np.flatnonzero(keep_d)] for k, v in pos.items()}
    # non-obvious: the drug is not indicated (cutoff graph) for a disease up to two ontology steps away
    O = g.M("disease_disease")
    O2 = ((O + O @ O) > 0).tocsr()
    I = g.M(INDICATION).tocsr()
    nonobv = []
    for c, d in enumerate(D):
        nb = O2.indices[O2.indptr[d]:O2.indptr[d + 1]]
        near = set(I[nb].indices) if len(nb) else set()
        nonobv.append(np.array([p for p in pos["primary"][c] if drugs[p] not in near], int))
    pos["non-obvious"] = nonobv
    # --- scores ---
    S_rank = np.zeros((len(drugs), len(D)), np.float32)
    S_logit = np.zeros_like(S_rank)
    V = np.zeros((len(drugs), len(D)), bool)
    io = paths.index("ontology")
    for i in range(0, len(D), 200):
        F = score_paths(g, D[i:i + 200], paths)
        for c in range(F.shape[1]):
            Z = np.column_stack([F[:, c, :], degree, gf])
            S_rank[:, i + c] = rank(Z)
            S_logit[:, i + c] = logit(Z[:, :nc])
            V[:, i + c] = F[:, c, io] > 0
    gn = Genes(g, D)
    T = gn.table("base")
    Ts = gn.table("base+specific3")
    hard = near_duplicates(g, D, pos["primary"], drug_pos) < 0.5
    out = ROOT / "data" / "eval" / f"temporal_C{C}_{feature_set}.npz"
    np.savez_compressed(out, D=D, drugs=drugs, degree=degree, cand=cand, S_rank=S_rank, S_logit=S_logit,
                        S_spec=_transform(S_rank.astype(float)).astype(np.float32), V=V, T=T, Ts=Ts,
                        sizes=gn.ng[D], hard=hard, leak=leak, removed=removed,
                        **{f"pos|{k}": as_objects(v) for k, v in pos.items()})
    print(f"  C={C} {feature_set}: {len(D)} test diseases, {sum(map(len, pos['primary']))} primary positives, "
          f"{removed} PrimeKG edges removed, {100 * leak:.0f} % of PrimeKG indications undated "
          f"({time.time() - t0:.0f} s) -> {out}", flush=True)


def as_objects(arrays):
    out = np.empty(len(arrays), dtype=object)
    out[:] = arrays
    return out


# --- analysis -----------------------------------------------------------------------------------

def recall_per_disease(S, cand, pos, k=100):
    from evaluation.novelty import measure
    rng = np.random.default_rng(0)
    r = np.full(S.shape[1], np.nan)
    for c in range(S.shape[1]):
        if len(pos[c]):
            m = measure(S[:, c].astype(np.float64), pos[c], cand[:, c], k=k, rng=rng)
            if m:
                r[c] = m["Recall@100"]
    return r


def boot_mean(v, rng):
    v = v[~np.isnan(v)]
    idx = rng.integers(0, len(v), (B, len(v)))
    m = v[idx].mean(axis=1)
    return float(v.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5)), len(v)


def boot_ratio(h, n, rng):
    idx = rng.integers(0, len(h), (B, len(h)))
    r = h[idx].sum(axis=1) / np.maximum(n[idx].sum(axis=1), 1)
    return float(h.sum() / max(n.sum(), 1)), float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5)), int(n.sum())


def top25(z):
    S, cand = z["S_spec"], z["cand"]
    return [np.flatnonzero(cand[:, c])[np.argsort(-S[cand[:, c], c])][:TOP] for c in range(S.shape[1])]


def pairs(z, key="primary"):
    """Every pair in the top 25 of each disease with the features of I2 and its outcome (for calibration)."""
    pos = z[f"pos|{key}"]
    rows = []
    for c, idx in enumerate(top25(z)):
        p = set(pos[c].tolist())
        for r, j in enumerate(idx, 1):
            sec = "extension" if z["V"][j, c] else ("shared genes" if z["T"][j, c] else "hidden")
            rows.append({"disease": int(z["D"][c]), "drug": int(z["drugs"][j]), "section": sec,
                         "rank": "1-5" if r <= 5 else ("6-10" if r <= 10 else "11-25"),
                         "size": "< 100" if z["sizes"][c] < 100 else ("100-499" if z["sizes"][c] < 500 else "≥ 500"),
                         "hard": bool(z["hard"][c]), "specific": bool(z["Ts"][j, c]), "y": int(j in p)})
    return pd.DataFrame(rows)


def analyse(C, feature_set=None):
    feature_set = feature_set or PRODUCTION
    z = dict(np.load(ROOT / "data" / "eval" / f"temporal_C{C}_{feature_set}.npz", allow_pickle=True))
    rng = np.random.default_rng(0)
    res = {"C": C, "feature_set": feature_set, "diseases": len(z["D"]), "leak": float(z["leak"]),
           "removed": int(z["removed"]), "recall": {}, "paired": {}, "sections": {}}
    deg = np.repeat(z["degree"][:, None], len(z["D"]), axis=1).astype(np.float32)
    for key in [k[4:] for k in z if k.startswith("pos|")]:
        pos = list(z[f"pos|{key}"])
        n_pos = sum(map(len, pos))
        r = {m: recall_per_disease(S, z["cand"], pos) for m, S in
             (("OREVID", z["S_rank"]), ("Popularity", deg), ("Logistic regression without popularity", z["S_logit"]))}
        res["recall"][key] = {m: boot_mean(v, rng) for m, v in r.items()} | {"positives": n_pos}
        d = r["OREVID"] - r["Popularity"]
        res["paired"][key] = boot_mean(d, rng)
    pos = list(z["pos|primary"])
    for c_name, sub in (("all", None), ("non-obvious", "non-obvious")):
        pp = list(z[f"pos|{sub}"]) if sub else pos
        h = {s: [] for s in ("extension", "shared genes", "hidden", "top 25")}
        n = {s: [] for s in h}
        for c, idx in enumerate(top25(z)):
            v, t = z["V"][idx, c], z["T"][idx, c]
            hit = np.isin(idx, pp[c])
            for s, m in (("extension", v), ("shared genes", ~v & t), ("hidden", ~v & ~t), ("top 25", np.ones_like(v))):
                h[s].append(int(hit[m].sum()))
                n[s].append(int(m.sum()))
        res["sections"][c_name] = {s: boot_ratio(np.array(h[s]), np.array(n[s]), rng) for s in h}
    return res


def report():
    sys.stdout.reconfigure(encoding="utf-8")
    out = {}
    for label, C in CUTOFFS.items():
        for fs in (PRODUCTION, "extra_no_ot"):
            f = ROOT / "data" / "eval" / f"temporal_C{C}_{fs}.npz"
            if f.exists():
                out[f"{label}|{fs}"] = analyse(C, fs)
    (ROOT / "data" / "eval" / "temporal.json").write_text(json.dumps(out, indent=1), "utf-8")
    fmt = lambda x: f"{x[0]:.3f} ({x[1]:.3f}–{x[2]:.3f})"
    pct = lambda x: f"{100 * x[0]:.1f} % ({100 * x[1]:.1f}–{100 * x[2]:.1f}; {x[3]:,})"
    L = ["# Temporal validation (I1)\n",
         "Generated by `py -3 -m evaluation.temporal report`. Protocol and decision rules: `docs/IMPROVEMENT_PLAN.md` (I1), "
         "fixed before running. The model, popularity and the logistic regression are trained on the graph known at each "
         "cutoff; positives are drug–disease pairs whose first trial started at or after the cutoff and that reached phase 3 "
         "or approval (Open Targets 26.06). Means over test diseases with 95 % intervals from resampling diseases. "
         "C1 = 2012 was the development cutoff, C2 = 2015 the confirmatory one (reported once).\n"]
    for key, r in out.items():
        label, fs = key.split("|")
        L += [f"## {label}, cutoff {r['C']}, features `{fs}`\n",
              f"{r['diseases']:,} test diseases; {r['removed']:,} PrimeKG drug–disease edges removed (first trial at or after "
              f"the cutoff); {100 * r['leak']:.0f} % of PrimeKG indications have no Open Targets date and stay in the training "
              f"graph (residual leakage).\n",
              "| Positives | n | OREVID | Popularity | Logistic regression without popularity | OREVID − popularity (paired) |",
              "|---|---|---|---|---|---|"]
        for k, v in r["recall"].items():
            L.append(f"| {k} | {v['positives']:,} | {fmt(v['OREVID'])} | {fmt(v['Popularity'])} | "
                     f"{fmt(v['Logistic regression without popularity'])} | {fmt(r['paired'][k])} |")
        L += ["", "Hit rate among the top 25 per disease, by website section (primary positives):\n",
              "| Positives | Direct extensions | Shared genes | Hidden | Whole top 25 |", "|---|---|---|---|---|"]
        for k, v in r["sections"].items():
            L.append(f"| {k} | " + " | ".join(pct(v[s]) for s in ("extension", "shared genes", "hidden", "top 25")) + " |")
        L.append("")
    c2 = out.get(f"C2 (confirmation)|{PRODUCTION}")
    if c2:
        lo, hi = c2["paired"]["primary"][1:3]
        verdict = ("above" if lo > 0 else ("below" if hi < 0 else "not distinguishable from"))
        L += [f"**Decision rule at C2 (primary positives):** the paired OREVID − popularity difference in Recall@100 has a "
              f"95 % interval of {lo:.3f} to {hi:.3f}, so OREVID ranks later phase 3 / approved uses **{verdict}** the "
              "popularity baseline.\n"]
    (ROOT / "docs" / "TEMPORAL.md").write_text("\n".join(L) + "\n", "utf-8")
    print("\n".join(L))


if __name__ == "__main__":
    if sys.argv[1] == "run":
        sys.stdout.reconfigure(encoding="utf-8")
        run(int(sys.argv[2]), sys.argv[3] if len(sys.argv) > 3 else PRODUCTION)
    elif sys.argv[1] == "report":
        report()
