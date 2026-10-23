"""P1: measured off-targets (ChEMBL, GtoPdb) as a drug–disease feature (docs/IMPROVEMENT_PLAN.md, «Registered
23-10-2026, before any run: P1»; data from sources/offtargets.py).

Feature: log(1 + number of the disease's genes, in the split's training graph, that are off-targets of the drug).
(a) gate: standalone on the benchmark's no-prior-use positives against random; (b) only if (a) passes: added LambdaRank
column, the N2 rule (evaluation/proximity_feature.py pattern).

    py -3 -m evaluation.p1_offtargets            # (a), then (b) if (a) passes -> docs/P1_OFFTARGETS.md
    py -3 -m evaluation.p1_offtargets report
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

import evaluation.proximity_feature as pf
import evaluation.txgnn as tx
import evaluation.xgboost_compare as xc
from evaluation.proximity import TASKS as TASK_DIR, boot, recall_at

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "eval" / "p1_feature"
GATE = ROOT / "data" / "eval" / "p1_gate.json"
GATE_B = ROOT / "data" / "eval" / "p1b_gate.json"


class OffTargets:
    """drugs (graph order) × diseases: log1p of the disease genes among the drug's off-targets."""

    def __init__(self, g, hub=False):
        m = pd.read_parquet(ROOT / "data" / "ext" / "measured_targets.parquet")
        if hub:  # P1b: plus the Repurposing Hub's annotated targets
            m = pd.concat([m[["drug", "gene"]], pd.read_parquet(ROOT / "data" / "ext" / "repurposing_hub.parquet")[["drug", "gene"]]])
        m = m[(m.drug < g.N) & (m.gene < g.N)]
        Mt = sp.csr_matrix((np.ones(len(m)), (m.drug, m.gene)), shape=(g.N, g.N))
        known = sum((g.M(r) for r in g.relations if r.startswith("drug_protein")), sp.csr_matrix((g.N, g.N)))
        off = (Mt > 0).astype(np.float32) - (Mt > 0).multiply(known > 0).astype(np.float32)
        off.eliminate_zeros()
        self.off = off.tocsr()[g.drugs]
        self.G = (g.M(tx.GENES_EXTRA) > 0).astype(np.float32).tocsr()
        self.n_drugs = int((np.diff(self.off.indptr) > 0).sum())
        self.n_pairs = int(self.off.nnz)

    def column(self, diseases):
        return np.log1p((self.off @ self.G[np.asarray(diseases, int)].T).toarray()).astype(np.float32)


def gate():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.txgnn import TxGraph, load_kg, node_types, split_edges
    kg = load_kg()
    types = node_types(kg)
    rng = np.random.default_rng(0)
    rows, cover = [], {}
    for s in range(1, 6):
        task = json.loads((TASK_DIR / f"zero_shot_{s}.json").read_text())
        g = TxGraph(kg, split_edges(kg, "complex_disease", s) == 0, types, extra=True)
        O = OffTargets(g)
        cover[f"zero_shot_{s}"] = {"drugs_with_off_targets": O.n_drugs, "drug_off_target_pairs": O.n_pairs}
        drugs = np.asarray(task["drugs"])
        pos_of = {int(f): j for j, f in enumerate(g.drugs)}
        idx = np.array([pos_of[int(f)] for f in drugs])
        pop = np.asarray(task["train_popularity"])
        D = [int(r["disease"]) for r in task["diseases"]]
        S = O.column(D)[idx]
        for c, r in enumerate(task["diseases"]):
            pos_all = np.isin(drugs, r["positives"])
            cand = ~(pos_all & (pop > 0))
            pos = pos_all & cand & (pop == 0)
            if not pos.any():
                continue
            ci = np.flatnonzero(cand)
            P = set(np.flatnonzero(pos[ci]).tolist())
            a, _ = recall_at(S[ci, c], P, rng=rng)
            rows.append({"feature": a, "random": min(1.0, 100 / len(ci)), "nonzero": bool((S[ci, c] > 0).any())})
        print(f"  zero_shot_{s}: {len(rows)} diseases so far", flush=True)
    res = {"feature": boot([r["feature"] for r in rows], rng), "random": boot([r["random"] for r in rows], rng),
           "difference": boot([r["feature"] - r["random"] for r in rows], rng), "diseases": len(rows),
           "diseases_with_a_nonzero_score": sum(r["nonzero"] for r in rows), "coverage": cover}
    res["passes"] = bool(res["difference"] and res["difference"][1] > 0)
    GATE.write_text(json.dumps(res, indent=1), "utf-8")
    return res


def gate_b():
    """P1b: the gate with ChEMBL + GtoPdb + Repurposing Hub targets, on the benchmark and the temporal 2015 set."""
    res = gate_sets(lambda g: OffTargets(g, hub=True))
    GATE_B.write_text(json.dumps(res, indent=1), "utf-8")
    return res


def gate_sets(make_column):
    """Standalone gate of a column (make_column(graph).column(diseases) -> drugs × diseases, graph drug order) on the
    benchmark and the temporal 2015 no-prior-use positives, against random (P1b, P2)."""
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.temporal import clinical, cutoff_graph, cutoff_mask
    from evaluation.txgnn import TxGraph, load_kg, node_types, split_edges
    kg = load_kg()
    types = node_types(kg)
    rng = np.random.default_rng(0)
    rows = []
    for s in range(1, 6):
        task = json.loads((TASK_DIR / f"zero_shot_{s}.json").read_text())
        g = TxGraph(kg, split_edges(kg, "complex_disease", s) == 0, types, extra=True)
        O = make_column(g)
        drugs = np.asarray(task["drugs"])
        pos_of = {int(f): j for j, f in enumerate(g.drugs)}
        idx = np.array([pos_of[int(f)] for f in drugs])
        pop = np.asarray(task["train_popularity"])
        S = O.column([int(r["disease"]) for r in task["diseases"]])[idx]
        for c, r in enumerate(task["diseases"]):
            pos_all = np.isin(drugs, r["positives"])
            cand = ~(pos_all & (pop > 0))
            pos = pos_all & cand & (pop == 0)
            if pos.any():
                ci = np.flatnonzero(cand)
                a, _ = recall_at(S[ci, c], set(np.flatnonzero(pos[ci]).tolist()), rng=rng)
                rows.append({"feature": a, "random": min(1.0, 100 / len(ci))})
    bench = {"feature": boot([r["feature"] for r in rows], rng), "random": boot([r["random"] for r in rows], rng),
             "difference": boot([r["feature"] - r["random"] for r in rows], rng), "diseases": len(rows)}
    z = dict(np.load(ROOT / "data" / "eval" / "temporal_C2015_extra_no_ctd.npz", allow_pickle=True))
    clin, attrs = clinical()
    keep, _, _ = cutoff_mask(kg, clin, 2015)
    g = cutoff_graph(kg, types, keep, attrs, 2015)
    assert np.array_equal(g.drugs, z["drugs"])
    O = make_column(g)
    S = O.column([int(d) for d in z["D"]])
    prior = z["degree"] > 0
    rows = []
    for c in range(len(z["D"])):
        pj = np.asarray(z["pos|primary"][c], int)
        cand = z["cand"][:, c] & (~prior | ~np.isin(np.arange(len(prior)), pj))
        pj = pj[~prior[pj]]
        ci = np.flatnonzero(cand)
        P = set(np.flatnonzero(np.isin(ci, pj)).tolist())
        if P:
            a, _ = recall_at(S[ci, c], P, rng=rng)
            rows.append({"feature": a, "random": min(1.0, 100 / len(ci))})
    temporal = {"feature": boot([r["feature"] for r in rows], rng), "random": boot([r["random"] for r in rows], rng),
                "difference": boot([r["feature"] - r["random"] for r in rows], rng), "diseases": len(rows)}
    res = {"benchmark": bench, "temporal": temporal, "drugs_with_off_targets": getattr(O, "n_drugs", None),
           "pairs": getattr(O, "n_pairs", None)}
    res["passes"] = bool(any(x["difference"] and x["difference"][1] > 0 for x in (bench, temporal)))
    return res


STATE = {"O": None}
_orig = tx.score_paths


def patched(h, diseases, paths):
    F = _orig(h, diseases, paths)
    if STATE["O"] is None:  # the first call of build() is on the split's training graph
        STATE["O"] = OffTargets(h, hub=GATE_B.exists())  # P1b passed: with the Repurposing Hub targets
    return np.concatenate([F, STATE["O"].column(diseases)[:, :, None]], axis=2)


def run(split, seed, kg, types):
    task = f"zero_shot_{seed}" if split == "complex_disease" else f"area_{split}"
    f = OUT / f"{task}.json"
    if f.exists():
        return json.loads(f.read_text())
    t0 = time.time()
    STATE["O"] = None
    xc.OUT = OUT / "rows"
    d = xc.build(kg, types, split, seed)
    X, Y, G = d["X"], d["Y"], d["groups"]
    k = d["F_test"].shape[2] - 1
    res = {}
    for name, Xs, Fs in (("with off-targets", X, d["F_test"]),
                         ("production (same rows, no off-targets)", np.delete(X, k, axis=1), d["F_test"][:, :, :k])):
        predict = xc._fit(xc.PRODUCTION_MODEL, Xs, Y, G, seed)
        S = xc.scores_for(predict, Fs, d["degree"], d["gf"])
        a, n = pf.recall_rows(S, list(d["test_pos"]), d["degree"], np.random.default_rng(0))
        res[name] = {"all": a.tolist(), "no_prior_use": n.tolist()}
    out = {"task": task, "diseases": [int(x) for x in d["diseases"]], **res}
    OUT.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out), "utf-8")
    print(f"  {task}: done ({time.time() - t0:.0f} s)", flush=True)
    return out


def feature_summary():
    rng = np.random.default_rng(0)
    W, P0 = "with off-targets", "production (same rows, no off-targets)"
    groups = {"benchmark (zero-shot 1-5)": [f"zero_shot_{s}" for s in range(1, 6)], "areas": [f"area_{a}" for a in tx.AREAS]}
    res = {}
    for gname, tasks in groups.items():
        rs = [json.loads((OUT / f"{t}.json").read_text()) for t in tasks if (OUT / f"{t}.json").exists()]
        if len(rs) < len(tasks):
            res[gname] = {"done": f"{len(rs)}/{len(tasks)}"}
            continue
        cat = lambda key, sub: np.concatenate([np.array(r[key][sub], float) for r in rs])
        mean_task = lambda key, sub: float(np.mean([np.nanmean(np.array(r[key][sub], float)) for r in rs]))
        res[gname] = {"no_prior_use_diff": pf.boot_diff(cat(W, "no_prior_use"), cat(P0, "no_prior_use"), rng),
                      "all_diff": pf.boot_diff(cat(W, "all"), cat(P0, "all"), rng),
                      "all_with": mean_task(W, "all"), "all_prod": mean_task(P0, "all"),
                      "novel_with": mean_task(W, "no_prior_use"), "novel_prod": mean_task(P0, "no_prior_use")}
    b, a = res.get("benchmark (zero-shot 1-5)", {}), res.get("areas", {})
    res["complete"] = "done" not in b and "done" not in a
    res["adopted"] = bool(b.get("no_prior_use_diff") and a.get("no_prior_use_diff") and b["no_prior_use_diff"][1] > 0
                          and a["no_prior_use_diff"][1] > 0 and b["all_with"] - b["all_prod"] >= -0.01)
    return res


def report():
    meta = json.loads((ROOT / "data" / "ext" / "measured_targets.json").read_text("utf-8"))
    gt = json.loads(GATE.read_text("utf-8")) if GATE.exists() else None
    f3 = lambda x: "—" if not x else f"{x[0]:.3f} ({x[1]:.3f}–{x[2]:.3f})"
    L = ["# P1. Measured off-targets (ChEMBL, GtoPdb)\n",
         "Generated by `py -3 -m evaluation.p1_offtargets`. Protocol and decision rules: `docs/IMPROVEMENT_PLAN.md`, "
         "«Registered 23-10-2026, before any run: P1». Feature: log(1 + the disease's genes, in the split's training "
         "graph, that are measured targets of the drug outside its known targets).\n",
         f"Sources: {meta['chembl']} activities (pChEMBL ≥ 6, IC50 / Ki / Kd / EC50, human single proteins; "
         f"{meta['chembl_molecules_queried']:,} molecules queried, {meta['chembl_activity_rows']:,} activity rows) and "
         f"{meta['gtopdb']} (human, affinity ≥ 6). Drug–gene pairs: {meta['pairs']:,} ({meta['drugs']:,} drugs, "
         f"{meta['genes']:,} genes; " + ", ".join(f"{k} {v:,}" for k, v in meta["pairs_by_source"].items()) + ").\n"]
    fd = lambda x: "—" if not x else f"{x[0]:+.3f} ({x[1]:+.3f} to {x[2]:+.3f})"
    if gt:
        cv = gt["coverage"]["zero_shot_1"]
        L += ["## (a) Gate: standalone on the no-prior-use positives (zero-shot 1–5)\n",
              f"Off-targets outside the graph's targets (zero-shot 1 training graph): {cv['drugs_with_off_targets']:,} "
              f"drugs, {cv['drug_off_target_pairs']:,} pairs. Diseases with no-prior-use positives: {gt['diseases']} "
              f"({gt['diseases_with_a_nonzero_score']} with a non-zero score for some candidate).\n",
              "| Ranking | Recall@100 (95 % CI) |", "|---|---|",
              f"| off-target feature | {f3(gt['feature'])} |", f"| random | {f3(gt['random'])} |",
              f"| difference (feature − random) | {fd(gt['difference'])} |",
              "", f"**Gate:** {'passes → (b) run' if gt['passes'] else 'does not pass → (b) not run; not adopted'}.\n"]
    gb = json.loads(GATE_B.read_text("utf-8")) if GATE_B.exists() else None
    if gb:
        L += ["## P1b: with the Repurposing Hub's targets, on the benchmark and the temporal 2015 set\n",
              f"Off-targets outside the graph's targets (2015 graph): {gb['drugs_with_off_targets']:,} drugs, {gb['pairs']:,} pairs.\n",
              "| Set | Diseases | Feature | Random | Difference |", "|---|---|---|---|---|"]
        for k in ("benchmark", "temporal"):
            x = gb[k]
            L.append(f"| {k} | {x['diseases']} | {f3(x['feature'])} | {f3(x['random'])} | {fd(x['difference'])} |")
        L += ["", f"**P1b gate:** {'passes → (b) run' if gb['passes'] else 'does not pass on either set; not adopted'}.\n"]
    if (gt and gt["passes"]) or (gb and gb["passes"]):
        fs = feature_summary()
        fmt = lambda x: "—" if not x else f"{x[0]:+.3f} ({x[1]:+.3f} to {x[2]:+.3f}; {x[3]} diseases)"
        L += ["## (b) As an added LambdaRank feature\n",
              "| Group | Recall@100, all (production → with) | Difference | Recall@100, no prior use | Difference |",
              "|---|---|---|---|---|"]
        for gname in ("benchmark (zero-shot 1-5)", "areas"):
            r = fs[gname]
            if "done" in r:
                L.append(f"| {gname} | {r['done']} tasks done | | | |")
                continue
            L.append(f"| {gname} | {r['all_prod']:.3f} → {r['all_with']:.3f} | {fmt(r['all_diff'])} | "
                     f"{r['novel_prod']:.3f} → {r['novel_with']:.3f} | {fmt(r['no_prior_use_diff'])} |")
        L += ["", "**Decision (rule fixed before running):** " + ("pending (tasks still running)" if not fs["complete"]
              else "adopted" if fs["adopted"] else "not adopted") + ".\n"]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "P1_OFFTARGETS.md").write_text(text, "utf-8")
    print(text)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    if sys.argv[1:] == ["report"]:
        return report()
    if sys.argv[1:] == ["b"]:
        gate_b()
        return report()
    gt = json.loads(GATE.read_text("utf-8")) if GATE.exists() else gate()
    gb = json.loads(GATE_B.read_text("utf-8")) if GATE_B.exists() else {"passes": False}
    if gt["passes"] or gb["passes"]:
        tx.score_paths = patched
        xc.score_paths = patched
        kg = tx.load_kg()
        types = tx.node_types(kg)
        for split, seed in pf.TASKS:
            run(split, seed, kg, types)
    report()


if __name__ == "__main__":
    main()
