"""Correcting N in the shared-genes hypergeometric test (docs/HYPERGEOMETRIC_N.md, P40).

Production used N = every node with a gene-disease link (genes and diseases). This script compares, with the
corrected N (gene nodes only), the thresholds {1e-6, 3e-6, 1e-5, 3e-5, 1e-4}: the one whose shown shared-genes
candidates agree best with production on the benchmark and the areas is chosen, and it is accepted if on the
external exam the hit rate stays within production's 95% interval and the number of candidates changes < 10 %.

    py -3 -m evaluation.hypergeom_n        # -> data/eval/hypergeom_n.json (printed summary)
"""
import json
import sys
import time
from pathlib import Path

import numpy as np

from evaluation.evidence_ci import Genes, bench_task, bootstrap, external_task

ROOT = Path(__file__).resolve().parent.parent
THRESHOLDS = (1e-6, 3e-6, 1e-5, 3e-5, 1e-4)
PRODUCTION_INTERVAL = (0.031, 0.055)  # external exam, shared genes, production (docs/EVIDENCE_CI.md)


def table(gn, N, t, k=3):
    if N not in gn.__dict__.setdefault("_pv", {}):
        gn._pv[N] = gn.pvalues(N)
    sig = (gn.K >= k) & (gn._pv[N] < t)
    return np.asarray(gn.IndM @ sig.T.astype(np.float32)) > 0


def shown(t_, T):
    """(disease column, drug index) pairs shown as shared-genes candidates, and per-disease hits/candidates."""
    pairs, h, n = set(), [], []
    for c, idx in enumerate(t_["top"]):
        sel = idx[~t_["V"][idx, c] & T[idx, c]]
        pairs |= {(c, int(j)) for j in sel}
        h.append(int(t_["hits"][c][sel].sum()))
        n.append(len(sel))
    return pairs, np.array(h), np.array(n)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.txgnn import AREAS, load_kg, node_types
    t0 = time.time()
    kg = load_kg()
    types = node_types(kg)
    rng = np.random.default_rng(0)
    tasks = {"benchmark": [("complex_disease", s) for s in (1, 2, 3, 4, 5)], "areas": [(a, 1) for a in AREAS]}
    res = {"N": {}, "groups": {}}
    inter = {t: 0 for t in THRESHOLDS}
    union = {t: 0 for t in THRESHOLDS}
    parts = {}  # (group, rule) -> list of (h, n); rule = "production" | t | ("removed"|"added", t)

    def run(group, t_, k):
        gn = Genes(t_["g"], t_["D"])
        res["N"] = {"production (all linked nodes)": gn.N_linked, "genes only": gn.N_genes}
        base, hb, nb = shown(t_, table(gn, gn.N_linked, 1e-6))  # production before P40
        parts.setdefault((group, "production"), []).append((hb, nb))
        for t in THRESHOLDS:
            new, h, n = shown(t_, table(gn, gn.N_genes, t))
            parts.setdefault((group, t), []).append((h, n))
            if group != "external":
                inter[t] += len(base & new)
                union[t] += len(base | new)
            for name, s in (("removed", base - new), ("added", new - base)):
                hh = sum(bool(t_["hits"][c][j]) for c, j in s)
                parts.setdefault((group, name, t), []).append((np.array([hh]), np.array([len(s)])))
        print(f"  {k} ({time.time() - t0:.0f} s)", flush=True)

    for group, ts in tasks.items():
        for split, seed in ts:
            run(group, bench_task(kg, types, split, seed), f"{split} {seed}")
    run("external", external_task(kg, types), "external")

    jac = {t: inter[t] / union[t] for t in THRESHOLDS}
    best = max(THRESHOLDS, key=lambda t: jac[t])
    res["jaccard_benchmark_areas"] = {f"{t:.0e}": round(v, 4) for t, v in jac.items()}
    res["chosen"] = best
    for (group, *rule), v in parts.items():
        key = f"{group}|{'|'.join(str(r) for r in rule)}"
        res["groups"][key] = bootstrap(v, rng)
    ext_prod = res["groups"]["external|production"]
    ext_new = res["groups"][f"external|{best}"]
    ok_rate = PRODUCTION_INTERVAL[0] <= ext_new[0] <= PRODUCTION_INTERVAL[1]
    ok_n = abs(ext_new[3] - ext_prod[3]) / ext_prod[3] < 0.10
    res["accepted"] = bool(ok_rate and ok_n)
    (ROOT / "data" / "eval" / "hypergeom_n.json").write_text(json.dumps(res, indent=1), "utf-8")
    print(json.dumps(res["N"]), "\nJaccard with production (benchmark + areas):", res["jaccard_benchmark_areas"])
    print(f"chosen threshold {best:.0e}")
    for group in ("benchmark", "areas", "external"):
        for rule in ("production", best, ("removed", best), ("added", best)):
            key = f"{group}|{rule if not isinstance(rule, tuple) else '|'.join(map(str, rule))}"
            r = res["groups"].get(key)
            if r:
                print(f"  {key:32s} {100 * r[0]:5.1f} % ({100 * r[1]:.1f}–{100 * r[2]:.1f}; {r[3]:,})")
    print(f"external: rate within production interval {ok_rate}; candidates within 10 % {ok_n} -> "
          f"{'ACCEPTED' if res['accepted'] else 'NOT accepted'}")


if __name__ == "__main__":
    main()
