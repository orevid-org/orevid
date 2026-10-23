"""M2: the genetic direction of M1 as a LambdaRank column (docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, after M1,
by the author's decision: M2»).

Column for (drug, disease): targets whose action agrees with the QTL colocalisation direction minus those that oppose
it (pQTL + eQTL, M1's filters; evaluation/m1_direction.py), NaN when the drug has no directed target for the disease.
Production rows and model with and without the column (evaluation/proximity_feature.py pattern).

    py -3 -m evaluation.m2_direction_feature [split seed ...]   # -> data/eval/m2_feature/*.json, docs/M2_DIRECTION_FEATURE.md
    py -3 -m evaluation.m2_direction_feature report
"""
import json
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np

import evaluation.proximity_feature as pf
import evaluation.txgnn as tx
import evaluation.xgboost_compare as xc

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "eval" / "m2_feature"


class GeneticDirection:
    def __init__(self):
        from evaluation.m1_direction import drug_actions, gene_disease_signs
        from evaluation.twins import twins
        from graph import Graph
        g = Graph()
        signs, _, _, _ = gene_disease_signs(g, ["pqtl", "eqtl"])
        tw = twins(g, "disease")
        self.by_dis = defaultdict(dict)
        for (gene, d), s in signs.items():
            for d2 in tw.get(d, {d}):  # signs are keyed by the twin group's smallest node
                self.by_dis[d2][gene] = s
        self.acts = drug_actions(g)

    def column(self, drugs, diseases):
        out = np.full((len(drugs), len(diseases)), np.nan, np.float32)
        pos = {int(f): j for j, f in enumerate(drugs)}
        for c, d in enumerate(diseases):
            gs = self.by_dis.get(int(d))
            if not gs:
                continue
            for f, tg in self.acts.items():
                j = pos.get(f)
                common = tg.keys() & gs.keys()
                if j is not None and common:
                    out[j, c] = sum(gs[gn] * -tg[gn] for gn in common)  # +1 per agreeing target
        return out


STATE = {"G": None, "drugs": None}
_orig = tx.score_paths


def patched(h, diseases, paths):
    F = _orig(h, diseases, paths)
    if STATE["G"] is None:
        STATE["G"] = GeneticDirection()
    return np.concatenate([F, STATE["G"].column(h.drugs, diseases)[:, :, None]], axis=2)


def run(split, seed, kg, types):
    task = f"zero_shot_{seed}" if split == "complex_disease" else f"area_{split}"
    f = OUT / f"{task}.json"
    if f.exists():
        return json.loads(f.read_text())
    t0 = time.time()
    xc.OUT = OUT / "rows"
    d = xc.build(kg, types, split, seed)
    X, Y, G = d["X"], d["Y"], d["groups"]
    k = d["F_test"].shape[2] - 1
    directed = np.isfinite(d["F_test"][:, :, k]).any(axis=0)
    res = {}
    for name, Xs, Fs in (("with direction", X, d["F_test"]),
                         ("production (same rows, no direction)", np.delete(X, k, axis=1), d["F_test"][:, :, :k])):
        predict = xc._fit(xc.PRODUCTION_MODEL, Xs, Y, G, seed)
        S = xc.scores_for(predict, Fs, d["degree"], d["gf"])
        a, n = pf.recall_rows(S, list(d["test_pos"]), d["degree"], np.random.default_rng(0))
        res[name] = {"all": a.tolist(), "no_prior_use": n.tolist()}
    out = {"task": task, "diseases": [int(x) for x in d["diseases"]], "directed": directed.tolist(), **res}
    OUT.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out), "utf-8")
    print(f"  {task}: done, {int(directed.sum())} of {len(directed)} test diseases with a directed pair "
          f"({time.time() - t0:.0f} s)", flush=True)
    return out


W, P0 = "with direction", "production (same rows, no direction)"
GROUPS = {"benchmark (zero-shot 1-5)": [f"zero_shot_{s}" for s in range(1, 6)], "areas": [f"area_{a}" for a in tx.AREAS]}


def summary():
    rng = np.random.default_rng(0)
    res = {}
    for gname, tasks in GROUPS.items():
        rs = [json.loads((OUT / f"{t}.json").read_text()) for t in tasks if (OUT / f"{t}.json").exists()]
        if len(rs) < len(tasks):
            res[gname] = {"done": f"{len(rs)}/{len(tasks)}"}
            continue
        cat = lambda key, sub, only=False: np.concatenate(
            [np.array(r[key][sub], float)[np.array(r["directed"]) if only else slice(None)] for r in rs])
        mean_task = lambda key, sub: float(np.mean([np.nanmean(np.array(r[key][sub], float)) for r in rs]))
        res[gname] = {"all_diff": pf.boot_diff(cat(W, "all"), cat(P0, "all"), rng),
                      "novel_diff": pf.boot_diff(cat(W, "no_prior_use"), cat(P0, "no_prior_use"), rng),
                      "directed_all_diff": pf.boot_diff(cat(W, "all", True), cat(P0, "all", True), rng),
                      "all_with": mean_task(W, "all"), "all_prod": mean_task(P0, "all"),
                      "novel_with": mean_task(W, "no_prior_use"), "novel_prod": mean_task(P0, "no_prior_use"),
                      "directed_diseases": int(sum(sum(r["directed"]) for r in rs)),
                      "diseases": int(sum(len(r["directed"]) for r in rs))}
    b, a = res.get("benchmark (zero-shot 1-5)", {}), res.get("areas", {})
    ok = False
    if b.get("all_diff") and a.get("all_diff"):
        up = lambda x: x["all_diff"][1] > 0
        flat = lambda x: x["all_diff"][0] >= -0.005
        novel_ok = all(x["novel_diff"] is None or x["novel_diff"][0] >= -0.01 for x in (b, a))
        ok = bool(((up(b) and flat(a)) or (up(a) and flat(b))) and novel_ok)
    res["adopted"] = ok
    return res


def report():
    res = summary()
    fmt = lambda x: "—" if not x else f"{x[0]:+.4f} ({x[1]:+.4f} to {x[2]:+.4f}; {x[3]} diseases)"
    L = ["# M2. Genetic direction (QTL colocalisation) as a model feature\n",
         "Generated by `py -3 -m evaluation.m2_direction_feature`. Rule fixed before running: `docs/IMPROVEMENT_PLAN.md`, "
         "«Registered 23-10-2026, after M1, by the author's decision: M2». Production rows and model with and without the "
         "column, paired by disease; mean of the task means; differences pooled over the diseases of the group.\n",
         "| Group | Recall@100, all (production → with) | Difference, all | Difference, no prior use | "
         "Difference, diseases with a directed pair |", "|---|---|---|---|---|"]
    for gname in GROUPS:
        r = res[gname]
        if "done" in r:
            L.append(f"| {gname} | {r['done']} tasks done | | | |")
            continue
        L.append(f"| {gname} | {r['all_prod']:.4f} → {r['all_with']:.4f} | {fmt(r['all_diff'])} | {fmt(r['novel_diff'])} | "
                 f"{fmt(r['directed_all_diff'])} ({r['directed_diseases']} of {r['diseases']}) |")
    L += ["", f"**Decision (rule fixed before running):** {'adopted' if res['adopted'] else 'not adopted'}.\n"]
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "M2_DIRECTION_FEATURE.md").write_text(text, "utf-8")
    (OUT / "summary.json").write_text(json.dumps(res, indent=1), "utf-8") if OUT.exists() else None
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
