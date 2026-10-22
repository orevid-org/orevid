"""N2 (b): network proximity as an added LambdaRank feature (docs/IMPROVEMENT_PLAN.md, «Registered 22-10-2026», N2).

Run because (a) was above random on both the benchmark and the temporal set (docs/PROXIMITY.md). Same training rows,
masked rounds and test features as production (evaluation/xgboost_compare.build), with one more column: the proximity
of `evaluation/proximity.py` (−z of the closest distance, 50 random gene sets), computed once per disease on each
split's training graph (it depends only on genes, targets and the protein network, not on drug–disease edges, so the
masked graphs give the same value). The baseline is the production model fitted on exactly the same rows without that
column. Rule: adopted only if Recall@100 on the no-prior-use subset (positives with no training indication; the other
positives removed from the candidates) improves with a 95 % interval excluding 0 on both the benchmark (zero-shot seeds
1-5) and the areas, and overall zero-shot Recall@100 does not fall by more than 0.01.

    py -3 -m evaluation.proximity_feature [split seed ...]   # -> data/eval/prox_feature/*.json, docs/PROXIMITY.md
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.sparse.csgraph import dijkstra

import evaluation.txgnn as tx
import evaluation.xgboost_compare as xc
from evaluation.proximity import N_RANDOM, MIN_GENES, Proximity

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "eval" / "prox_feature"
TASKS = [("complex_disease", s) for s in (1, 2, 3, 4, 5)] + [(a, 1) for a in tx.AREAS]


class FastProximity(Proximity):
    """Same score as Proximity, with all gene-to-gene distances of the network computed once."""

    def __init__(self, g, seed=0):
        super().__init__(g, seed)
        n = len(self.net)
        self.D = np.empty((n, n), np.uint8)
        for i in range(0, n, 1000):
            d = dijkstra(self.A_net, directed=False, indices=np.arange(i, min(i + 1000, n)), unweighted=True)
            d[~np.isfinite(d)] = 255
            self.D[i:i + 1000] = np.minimum(d, 255).astype(np.uint8)
        self.cache = {}

    def _dc(self, genes):
        d = self.D[self.local[genes]].min(axis=0).astype(np.float64)
        finite = d < 255
        d[~finite] = (d[finite].max() + 1) if finite.any() else 50
        return np.asarray(self.Tn_net @ d).ravel()

    def column(self, diseases):
        out = np.full((len(self.has_target), len(diseases)), np.nan, np.float32)
        for c, d in enumerate(diseases):
            d = int(d)
            if d not in self.cache:
                s = self.score(d)
                self.cache[d] = None if s is None else s.astype(np.float32)
            if self.cache[d] is not None:
                out[:, c] = self.cache[d]
        return out


STATE = {"P": None}
_orig = tx.score_paths


def patched(h, diseases, paths):
    F = _orig(h, diseases, paths)
    if STATE["P"] is None:  # the first call of build() is on the split's training graph
        t0 = time.time()
        STATE["P"] = FastProximity(h, seed=0)
        print(f"    distances of {len(STATE['P'].net):,} network genes ({time.time() - t0:.0f} s)", flush=True)
    return np.concatenate([F, STATE["P"].column(diseases)[:, :, None]], axis=2)


def recall_rows(S, test_pos, degree, rng):
    """Per test disease: Recall@100 on all positives and on positives with no training indication."""
    from evaluation.novelty import measure
    n = S.shape[0]
    allr, novel = [], []
    for c, pos in enumerate(test_pos):
        pos = np.asarray(pos, int)
        if not len(pos):
            allr.append(np.nan), novel.append(np.nan)
            continue
        m = measure(S[:, c].astype(np.float64), pos, np.ones(n, bool), rng=rng)
        allr.append(m["Recall@100"] if m else np.nan)
        pn = pos[degree[pos] == 0]
        cand = np.ones(n, bool)
        cand[pos[degree[pos] > 0]] = False
        m = measure(S[:, c].astype(np.float64), pn, cand, rng=rng) if len(pn) else None
        novel.append(m["Recall@100"] if m else np.nan)
    return np.array(allr), np.array(novel)


def run(split, seed, kg, types):
    task = f"zero_shot_{seed}" if split == "complex_disease" else f"area_{split}"
    f = OUT / f"{task}.json"
    if f.exists():
        return json.loads(f.read_text())
    t0 = time.time()
    STATE["P"] = None
    xc.OUT = OUT / "rows"
    d = xc.build(kg, types, split, seed)
    X, Y, G = d["X"], d["Y"], d["groups"]
    k = d["F_test"].shape[2] - 1  # index of the proximity column among the path features
    res = {}
    for name, Xs, Fs in (("with proximity", X, d["F_test"]),
                         ("production (same rows, no proximity)", np.delete(X, k, axis=1), d["F_test"][:, :, :k])):
        predict = xc._fit(xc.PRODUCTION_MODEL, Xs, Y, G, seed)
        S = xc.scores_for(predict, Fs, d["degree"], d["gf"])
        a, n = recall_rows(S, list(d["test_pos"]), d["degree"], np.random.default_rng(0))
        res[name] = {"all": a.tolist(), "no_prior_use": n.tolist()}
    out = {"task": task, "diseases": [int(x) for x in d["diseases"]], **res,
           "with_proximity_score": int(sum(v is not None for v in STATE["P"].cache.values())),
           "diseases_scored": len(STATE["P"].cache)}
    OUT.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(out), "utf-8")
    print(f"  {task}: done ({time.time() - t0:.0f} s)", flush=True)
    return out


def boot_diff(a, b, rng, B=5000):
    m = ~np.isnan(a) & ~np.isnan(b)
    d = (a - b)[m]
    if not len(d):
        return None
    bs = rng.choice(d, (B, len(d))).mean(axis=1)
    return float(d.mean()), float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5)), int(len(d))


def report():
    rng = np.random.default_rng(0)
    W, P0 = "with proximity", "production (same rows, no proximity)"
    groups = {"benchmark (zero-shot 1-5)": [f"zero_shot_{s}" for s in range(1, 6)], "areas": [f"area_{a}" for a in tx.AREAS]}
    rows, res = [], {}
    for gname, tasks in groups.items():
        rs = [json.loads((OUT / f"{t}.json").read_text()) for t in tasks if (OUT / f"{t}.json").exists()]
        if len(rs) < len(tasks):
            rows.append(f"| {gname} | {len(rs)}/{len(tasks)} tasks done | | | | |")
            continue
        cat = lambda key, sub: np.concatenate([np.array(r[key][sub], float) for r in rs])
        mean_task = lambda key, sub: float(np.mean([np.nanmean(np.array(r[key][sub], float)) for r in rs]))
        dn = boot_diff(cat(W, "no_prior_use"), cat(P0, "no_prior_use"), rng)
        da = boot_diff(cat(W, "all"), cat(P0, "all"), rng)
        res[gname] = {"no_prior_use_diff": dn, "all_diff": da, "all_with": mean_task(W, "all"), "all_prod": mean_task(P0, "all"),
                      "novel_with": mean_task(W, "no_prior_use"), "novel_prod": mean_task(P0, "no_prior_use")}
        f = lambda x: "—" if not x else f"{x[0]:+.3f} ({x[1]:+.3f} to {x[2]:+.3f}; {x[3]} diseases)"
        rows.append(f"| {gname} | {res[gname]['all_prod']:.3f} → {res[gname]['all_with']:.3f} | {f(da)} | "
                    f"{res[gname]['novel_prod']:.3f} → {res[gname]['novel_with']:.3f} | {f(dn)} |")
    b, a = res.get("benchmark (zero-shot 1-5)"), res.get("areas")
    adopted = bool(b and a and b["no_prior_use_diff"] and a["no_prior_use_diff"] and b["no_prior_use_diff"][1] > 0
                   and a["no_prior_use_diff"][1] > 0 and b["all_with"] - b["all_prod"] >= -0.01)
    L = ["## (b) Proximity as an added LambdaRank feature\n",
         "`py -3 -m evaluation.proximity_feature`: production rows and model, with and without the proximity column "
         "(paired by disease; mean of the task means; differences pooled over the diseases of the group).\n",
         "| Group | Recall@100, all (production → with proximity) | Difference | Recall@100, no prior use | Difference |",
         "|---|---|---|---|---|"] + rows
    L += ["", f"**Decision (rule fixed before running):** {'adopted' if adopted else 'not adopted'}."]
    text = "\n".join(L) + "\n"
    p = ROOT / "docs" / "PROXIMITY.md"
    old = p.read_text("utf-8") if p.exists() else ""
    head = old.split("## (b) Proximity as an added LambdaRank feature")[0].rstrip()
    p.write_text(head + "\n\n" + text, "utf-8")
    (OUT / "summary.json").write_text(json.dumps({"results": res, "adopted": adopted}, indent=1), "utf-8")
    print(text)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    tx.score_paths = patched
    xc.score_paths = patched
    if sys.argv[1:] == ["report"]:
        return report()
    args = sys.argv[1:]
    tasks = [(args[i], int(args[i + 1])) for i in range(0, len(args), 2)] if args else TASKS
    kg = tx.load_kg()
    types = tx.node_types(kg)
    for split, seed in tasks:
        run(split, seed, kg, types)
    report()


if __name__ == "__main__":
    main()
