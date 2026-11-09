"""G3c: metapath 22 («similar all (OT)», Table S2) recomputed without animal-model-only gene–disease edges (the G3a gene
sets), production model retrained (docs/IMPROVEMENT_PLAN.md, «Registered 09-11-2026, on the work PC, before any run:
G3 …»). Same rows and model as production (evaluation/xgboost_compare.py, as O2): the column is computed both ways on
the same graph, and the model is fitted once with the production column and once with the G3a column in its place.

Rule (fixed before running): adopted only if the paired difference in Recall@100 (bootstrap over diseases) has a 95 %
interval whose lower bound is above −0.005 on both the benchmark (zero-shot seeds 1–5) and the areas, and the
external-exam Recall@100 is not lower than the lower bound of production's 95 % interval.

    py -3 -m evaluation.g3c_feature [split seed ...]   # -> data/eval/g3c_feature/*.json, docs/G3C_FEATURE.md
    py -3 -m evaluation.g3c_feature report
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import evaluation.proximity_feature as pf
import evaluation.txgnn as tx
import evaluation.xgboost_compare as xc
from graph import METAPATHS

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "eval" / "g3c_feature"
COL = "similar all (OT)"
G3A_REL = "genes_g3a"
G3A_PATH = "similar all (G3a)"
G3A_PIECES = ("disease_protein", "ot_assoc:genetic_association", "ot_assoc:genetic_literature", "ot_assoc:somatic_mutation",
              "ot_assoc:rna_expression")
METAPATHS[G3A_PATH] = ("similar", tx.INDICATION, G3A_REL, "as «similar all (OT)» without animal-model-only gene edges")
_orig = tx.score_paths


def patched(h, diseases, paths):
    """Production features plus the G3a version of «similar all (OT)» as a last column."""
    if G3A_REL not in h._m:
        m = sum(h.M(r) for r in G3A_PIECES if r in h.relations).tocsr()
        m.data[:] = 1
        h.register(G3A_REL, m)
    F = _orig(h, diseases, paths)
    return np.concatenate([F, _orig(h, diseases, [G3A_PATH])], axis=2)


def external_recall(S):
    """Per external-exam disease: Recall@100 on all positives (as evaluation/external_exam.py)."""
    from evaluation.novelty import measure
    z = np.load(tx.CACHE / "external_scores.npz", allow_pickle=True)
    p = pd.read_parquet(tx.CACHE / "external_positives.parquet")
    drugs = z["drugs"]
    pos_f = {int(f): j for j, f in enumerate(drugs)}
    rng = np.random.default_rng(0)
    out = []
    for c, (d, known) in enumerate(zip(z["diseases"], z["known"])):
        cand = np.ones(len(drugs), bool)
        cand[[pos_f[int(f)] for f in known if int(f) in pos_f]] = False
        pos = np.array([pos_f[int(f)] for f in p[p.disease == d].drug if int(f) in pos_f], int)
        m = measure(S[:, c].astype(np.float64), pos, cand, rng=rng) if len(pos) else None
        out.append(m["Recall@100"] if m else np.nan)
    return np.array(out)


W, P0 = "G3a column", "production column"


def run(split, seed, kg, types):
    task = f"zero_shot_{seed}" if split == "complex_disease" else f"area_{split}"
    f = OUT / f"{task}.json"
    if f.exists():
        return json.loads(f.read_text())
    t0 = time.time()
    xc.OUT = OUT / "rows"
    d = xc.build(kg, types, split, seed)
    X, Y, G = d["X"], d["Y"], d["groups"]
    k = d["F_test"].shape[2] - 1  # the G3a column (last path feature)
    j = tx.FEATURE_SETS[tx.PRODUCTION].index(COL)  # the production column it replaces
    Xw = np.delete(X, k, axis=1)
    Xw[:, j] = X[:, k]
    Fw = np.delete(d["F_test"], k, axis=2)
    Fw[:, :, j] = d["F_test"][:, :, k]
    variants = {W: (Xw, Fw), P0: (np.delete(X, k, axis=1), d["F_test"][:, :, :k])}
    res = {}
    for name, (Xs, Fs) in variants.items():
        predict = xc._fit(xc.PRODUCTION_MODEL, Xs, Y, G, seed)
        S = xc.scores_for(predict, Fs, d["degree"], d["gf"])
        a, n = pf.recall_rows(S, list(d["test_pos"]), d["degree"], np.random.default_rng(0))
        res[name] = {"all": a.tolist(), "no_prior_use": n.tolist()}
        if "F_ext" in d:
            E = d["F_ext"]
            Ee = np.delete(E, k, axis=2)
            if name == W:
                Ee[:, :, j] = E[:, :, k]
            res[name]["external"] = external_recall(xc.scores_for(predict, Ee, d["degree"], d["gf"])).tolist()
    changed = float((d["F_test"][:, :, k] != d["F_test"][:, :, j]).mean())
    out = {"task": task, "diseases": [int(x) for x in d["diseases"]], "changed_cells": changed, **res}
    OUT.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out), "utf-8")
    print(f"  {task}: done; {100 * changed:.2f} % of test cells change ({time.time() - t0:.0f} s)", flush=True)
    return out


GROUPS = {"benchmark (zero-shot 1-5)": [f"zero_shot_{s}" for s in range(1, 6)], "areas": [f"area_{a}" for a in tx.AREAS]}


def summary():
    rng = np.random.default_rng(0)
    res = {}
    for gname, tasks in GROUPS.items():
        rs = [json.loads((OUT / f"{t}.json").read_text()) for t in tasks if (OUT / f"{t}.json").exists()]
        if len(rs) < len(tasks):
            res[gname] = {"done": f"{len(rs)}/{len(tasks)}"}
            continue
        cat = lambda key, sub: np.concatenate([np.array(r[key][sub], float) for r in rs])
        mean_task = lambda key, sub: float(np.mean([np.nanmean(np.array(r[key][sub], float)) for r in rs]))
        res[gname] = {"all_diff": pf.boot_diff(cat(W, "all"), cat(P0, "all"), rng),
                      "novel_diff": pf.boot_diff(cat(W, "no_prior_use"), cat(P0, "no_prior_use"), rng),
                      "all_with": mean_task(W, "all"), "all_prod": mean_task(P0, "all"),
                      "novel_with": mean_task(W, "no_prior_use"), "novel_prod": mean_task(P0, "no_prior_use"),
                      "changed_cells": float(np.mean([r["changed_cells"] for r in rs]))}
    z1 = OUT / "zero_shot_1.json"
    if z1.exists():
        r = json.loads(z1.read_text())
        a, b = np.array(r[W]["external"], float), np.array(r[P0]["external"], float)
        m = ~np.isnan(b)
        bs = rng.choice(b[m], (5000, int(m.sum()))).mean(axis=1)
        res["external"] = {"with": float(np.nanmean(a)), "prod": float(np.nanmean(b)),
                           "prod_ci": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                           "diff": pf.boot_diff(a, b, rng)}
    b, a, e = res.get("benchmark (zero-shot 1-5)", {}), res.get("areas", {}), res.get("external")
    ok = bool(b.get("all_diff") and a.get("all_diff") and e and b["all_diff"][1] > -0.005 and a["all_diff"][1] > -0.005
              and e["with"] >= e["prod_ci"][0])
    res["adopted"] = ok
    return res


def report():
    res = summary()
    fmt = lambda x: "—" if not x else f"{x[0]:+.4f} ({x[1]:+.4f} to {x[2]:+.4f}; {x[3]} diseases)"
    L = ["# G3c. «similar all (OT)» without animal-model-only gene edges, model retrained\n",
         "Generated by `py -3 -m evaluation.g3c_feature`. Rule fixed before running: `docs/IMPROVEMENT_PLAN.md`, "
         "«Registered 09-11-2026, on the work PC, before any run: G3 …». Production rows and model (LambdaRank), fitted "
         "with the production column and with the G3a column in its place, paired by disease; mean of the task means; "
         "differences pooled over the diseases of the group.\n",
         "| Group | Recall@100, all (production → G3a) | Difference, all | Difference, no prior use | Test cells that change |",
         "|---|---|---|---|---|"]
    for gname in GROUPS:
        r = res[gname]
        if "done" in r:
            L.append(f"| {gname} | {r['done']} tasks done | | | |")
            continue
        L.append(f"| {gname} | {r['all_prod']:.4f} → {r['all_with']:.4f} | {fmt(r['all_diff'])} | {fmt(r['novel_diff'])} | "
                 f"{100 * r['changed_cells']:.2f} % |")
    e = res.get("external")
    if e:
        L.append(f"| external exam (seed 1) | {e['prod']:.4f} (95 % {e['prod_ci'][0]:.4f}–{e['prod_ci'][1]:.4f}) → "
                 f"{e['with']:.4f} | {fmt(e['diff'])} | | |")
    L += ["", f"**Decision (rule fixed before running):** {'adopted' if res['adopted'] else 'not adopted'}.\n"]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "G3C_FEATURE.md").write_text(text, "utf-8")
    if OUT.exists():
        (OUT / "summary.json").write_text(json.dumps(res, indent=1), "utf-8")
    print(text)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    if sys.argv[1:] == ["report"]:
        return report()
    tx.score_paths = patched
    xc.score_paths = patched
    args = sys.argv[1:]
    tasks = [(args[i], int(args[i + 1])) for i in range(0, len(args), 2)] if args else pf.TASKS
    kg = tx.load_kg()
    types = tx.node_types(kg)
    for split, seed in tasks:
        run(split, seed, kg, types)
    report()


if __name__ == "__main__":
    main()
