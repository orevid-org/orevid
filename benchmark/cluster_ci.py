"""Uncertainty of the Recall@100 comparisons, with diseases as the unit of resampling.

The benchmark reports means over diseases (and over seeds or areas). Positives of the same disease
are not independent, so the confidence intervals here resample whole diseases (cluster bootstrap):

- zero-shot: diseases resampled within each seed, then the mean of the 5 seed means (as reported);
- disease areas: diseases resampled within each area, then the mean of the 9 area means (as reported)
  and the mean weighted by disease;
- external exam: the 645 diseases resampled (one training split).

Model comparisons are paired by disease (same resample for both models), with the share of diseases
where each model is ahead and a Wilcoxon signed-rank test; the p-values are Holm-adjusted over all the
comparisons in the table. Tie-breaking: the benchmark breaks ties at random with seed 0; the means are
recomputed with 20 other seeds to show how much that choice moves them.

    python cluster_ci.py          # -> CLUSTER_CI.md, cluster_ci.json (needs the prediction files, see README)
"""
import json
from pathlib import Path

import numpy as np
from scipy.stats import wilcoxon

from bench import _jitter, load_predictions

HERE = Path(__file__).resolve().parent
MODELS = {"popularity": "Popularity", "txgnn": "TxGNN", "orevid_lambdarank": "OREVID",
          "orevid_no_popularity": "Logistic regression without popularity"}
PAIRS = [("orevid_lambdarank", "txgnn"), ("orevid_lambdarank", "popularity"), ("txgnn", "popularity"),
         ("orevid_no_popularity", "txgnn")]
B = 10_000
TIE_SEEDS = range(1, 21)


def split_recall(S, positives, seed=0):
    """Recall@100 per disease exactly as bench.txgnn_metrics (same random tie-breaking sequence)."""
    rng = np.random.default_rng(seed)
    out = []
    for c, pos in enumerate(positives):
        s = _jitter(S[:, c], rng)
        order = np.argsort(s)[::-1]
        out.append(len(np.intersect1d(order[:100], pos)) / len(pos))
    return np.array(out)


def external_recall(task, S, seed=0):
    """Recall@100 per disease and subset exactly as bench.evaluate_external."""
    drugs = np.asarray(task["drugs"])
    where = {v: i for i, v in enumerate(drugs)}
    rng = np.random.default_rng(seed)
    out = {}
    for sub in ["positives", "recent", "non_obvious"]:
        v = {}
        for c, r in enumerate(task["diseases"]):
            pos = np.array([where[p] for p in r[sub]], int)
            if not len(pos):
                continue
            cand = np.ones(len(drugs), bool)
            cand[[where[p] for p in r["exclude"] if p in where]] = False
            pos = np.intersect1d(pos, np.flatnonzero(cand))
            if not len(pos):
                continue
            idx = np.flatnonzero(cand)
            x = _jitter(S[idx, c], rng)
            rk = np.empty(len(idx))
            rk[np.argsort(-x)] = np.arange(len(idx))
            v[r["disease"]] = float(np.mean(rk[np.searchsorted(idx, pos)] < 100))
        out[{"positives": "all"}.get(sub, sub)] = v
    return out


def collect(seed=0):
    """{test: {model: {stratum: {disease: recall}}}} for every test and model."""
    res = {}
    for tf in sorted((HERE / "tasks").glob("*.json")):
        task = json.loads(tf.read_text())
        S = {m: load_predictions(f, task)[0] for m in MODELS
             if (f := HERE / "predictions" / m / f"{task['name']}.parquet").exists()}
        if task.get("type") == "external":
            for m, s in S.items():
                for sub, v in external_recall(task, s, seed).items():
                    res.setdefault(f"external, {sub}", {}).setdefault(m, {})["external"] = v
            continue
        drugs = np.asarray(task["drugs"])
        where = {v: i for i, v in enumerate(drugs)}
        dis = [r["disease"] for r in task["diseases"]]
        positives = [np.array([where[p] for p in r["positives"]]) for r in task["diseases"]]
        hard = [i for i, r in enumerate(task["diseases"]) if r["near_duplicate_overlap"] < 0.5]
        group = "zero-shot" if task["name"].startswith("zero_shot") else "areas"
        for m, s in S.items():
            r_all = split_recall(s, positives, seed)
            r_hard = split_recall(s[:, hard], [positives[i] for i in hard], seed)
            res.setdefault(group, {}).setdefault(m, {})[task["name"]] = dict(zip(dis, r_all.tolist()))
            res.setdefault(f"{group}, hard", {}).setdefault(m, {})[task["name"]] = dict(zip([dis[i] for i in hard], r_hard.tolist()))
    return res


def strata(test):
    """Arrays per stratum (seed or area) in a fixed disease order."""
    return {m: {k: np.array(list(v.values())) for k, v in by.items()} for m, by in test.items()}


def stat(arrs, weighted=False):
    if weighted:
        return float(np.mean(np.concatenate(list(arrs.values()))))
    return float(np.mean([a.mean() for a in arrs.values()]))


def bootstrap(test, models, weighted=False, seed=0):
    """Point estimate and 95% percentile intervals for each model and for paired differences."""
    rng = np.random.default_rng(seed)
    st = strata(test)
    keys = list(next(iter(st.values())).keys())
    draws = {k: rng.integers(0, len(st[models[0]][k]), size=(B, len(st[models[0]][k]))) for k in keys}
    boot = {}
    for m in models:
        per = [st[m][k][draws[k]] for k in keys]                    # B × n_k each
        if weighted:
            boot[m] = np.concatenate(per, axis=1).mean(axis=1)
        else:
            boot[m] = np.mean([p.mean(axis=1) for p in per], axis=0)
    return {m: (stat(st[m], weighted), *np.percentile(boot[m], [2.5, 97.5])) for m in models}, boot


def paired(test, a, b):
    """Per-disease differences a − b: share ahead/behind and Wilcoxon signed-rank p."""
    da = np.concatenate([np.array(list(v.values())) for v in test[a].values()])
    db = np.concatenate([np.array(list(v.values())) for v in test[b].values()])
    d = da - db
    p = wilcoxon(da, db, zero_method="wilcox").pvalue if np.any(d != 0) else 1.0
    return {"n": len(d), "ahead": float(np.mean(d > 0)), "behind": float(np.mean(d < 0)), "p": float(p)}


def holm(ps):
    order = np.argsort(ps)
    adj, run = np.empty(len(ps)), 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(ps) - r) * ps[i]))
        adj[i] = run
    return adj


def fmt(x):
    return f"{x[0]:.3f} ({x[1]:.3f}–{x[2]:.3f})"


def main():
    res = collect()
    tests = [("zero-shot", False), ("zero-shot, hard", False), ("areas", False), ("areas", True),
             ("areas, hard", False), ("external, all", False), ("external, non_obvious", False),
             ("external, recent", False)]
    rows, comps = [], []
    for test, w in tests:
        models = [m for m in MODELS if m in res[test]]
        est, boot = bootstrap(res[test], models, weighted=w)
        name = test + (" (weighted by disease)" if w else "")
        rows.append((name, {m: est[m] for m in models}))
        for a, b in PAIRS:
            if a in models and b in models:
                diff = boot[a] - boot[b]
                c = {"test": name, "a": a, "b": b, "diff": est[a][0] - est[b][0],
                     "lo": float(np.percentile(diff, 2.5)), "hi": float(np.percentile(diff, 97.5))}
                if not w:
                    c |= paired(res[test], a, b)
                comps.append(c)
    tested = [i for i, c in enumerate(comps) if "p" in c]
    for i, p in zip(tested, holm(np.array([comps[i]["p"] for i in tested]))):
        comps[i]["p_holm"] = float(p)

    # Tie-breaking: recompute the means with other seeds.
    ties = {}
    for s in TIE_SEEDS:
        r = collect(seed=s)
        for test, w in tests:
            for m in res[test]:
                ties.setdefault((test, w, m), []).append(stat(strata(r[test])[m], w))
    tie_rows = [(t, w, m, min(v), max(v)) for (t, w, m), v in ties.items()]

    L = ["# Recall@100 with disease-level uncertainty\n",
         "Generated by `cluster_ci.py`. 95% intervals from a bootstrap that resamples whole diseases "
         f"({B:,} resamples; within seeds for zero-shot and within areas for the disease areas), so the "
         "positives of a disease are never treated as independent. Differences are paired by disease. "
         "Wilcoxon signed-rank tests on the per-disease differences, Holm-adjusted over every comparison "
         "in the second table.\n",
         "## Recall@100 (95% cluster-bootstrap interval)\n",
         "| Test | " + " | ".join(MODELS.values()) + " |", "|---|" + "---|" * len(MODELS)]
    for name, est in rows:
        L.append(f"| {name} | " + " | ".join(fmt(est[m]) if m in est else "–" for m in MODELS) + " |")
    L += ["", "## Paired differences\n",
          "| Test | Comparison | Difference (95% CI) | Diseases | Ahead | Behind | Holm p |", "|---|---|---|---|---|---|---|"]
    for c in comps:
        extra = (f"{c['n']} | {c['ahead']:.0%} | {c['behind']:.0%} | {c['p_holm']:.1e}" if "p" in c else "– | – | – | –")
        L.append(f"| {c['test']} | {MODELS[c['a']]} − {MODELS[c['b']]} | {c['diff']:+.3f} ({c['lo']:+.3f} to {c['hi']:+.3f}) | {extra} |")
    L += ["", f"## Random tie-breaking\n",
          f"Range of the reported mean over {len(TIE_SEEDS)} other tie-breaking seeds (the benchmark uses seed 0).\n",
          "| Test | Model | Min | Max |", "|---|---|---|---|"]
    for t, w, m, lo, hi in tie_rows:
        L.append(f"| {t}{' (weighted)' if w else ''} | {MODELS[m]} | {lo:.3f} | {hi:.3f} |")
    (HERE / "CLUSTER_CI.md").write_text("\n".join(L) + "\n", encoding="utf-8")
    (HERE / "cluster_ci.json").write_text(json.dumps(
        {"estimates": {n: {m: list(v) for m, v in e.items()} for n, e in rows}, "comparisons": comps,
         "tie_breaking": [{"test": t, "weighted": w, "model": m, "min": lo, "max": hi} for t, w, m, lo, hi in tie_rows]},
        indent=1))
    print("-> CLUSTER_CI.md")


if __name__ == "__main__":
    main()
