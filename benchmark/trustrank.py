"""TrustRank baseline, as in NeDRex (Sadegh et al., Nat Commun 2021), implemented locally on OREVID-bench's data
(docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, night: B1»). Not the NeDRex service: the same idea on the task's
training edges. Personalised PageRank (damping 0.85, 50 iterations) on the protein–protein network, restarting at the
disease's genes; a drug's score = the sum of the TrustRank of its targets.

    python trustrank.py --kg kg.csv          # -> predictions/trustrank/<task>.parquet
"""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from bench import HERE, load_kg, train_edges

DAMPING, ITER = 0.85, 50


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--kg", required=True)
    p.add_argument("--tasks", default=str(HERE / "tasks"))
    p.add_argument("--out", default=str(HERE / "predictions" / "trustrank"))
    a = p.parse_args()
    kg = load_kg(a.kg)
    n = int(max(kg.x_index.max(), kg.y_index.max())) + 1
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    for tf in sorted(Path(a.tasks).glob("*.json")):
        task = json.loads(tf.read_text())
        mask = train_edges(kg, a.tasks, task.get("train_split") or task["name"])
        e = kg[mask]

        def sym(rel):
            d = e[e.relation == rel]
            x, y = d.x_index.to_numpy(), d.y_index.to_numpy()
            m = sp.coo_matrix((np.ones(2 * len(d)), (np.r_[x, y], np.r_[y, x])), shape=(n, n)).tocsr()
            m.data[:] = 1
            return m
        ppi, dp, tg = sym("protein_protein"), sym("disease_protein"), sym("drug_protein")
        deg = np.asarray(ppi.sum(axis=1)).ravel()
        W = (ppi @ sp.diags(1 / np.maximum(deg, 1))).tocsr()  # column-stochastic on proteins with edges
        drugs = np.asarray(task["drugs"])
        T = tg[drugs]
        rows = []
        for r in task["diseases"]:
            d = r["disease"]
            seeds = dp.indices[dp.indptr[d]:dp.indptr[d + 1]]
            if not len(seeds):
                score = np.zeros(len(drugs))
            else:
                s = np.zeros(n)
                s[seeds] = 1 / len(seeds)
                v = s.copy()
                for _ in range(ITER):
                    v = (1 - DAMPING) * s + DAMPING * (W @ v)
                score = np.asarray(T @ v).ravel()
            rows.append(pd.DataFrame({"disease": d, "drug": drugs, "score": score}))
        pd.concat(rows).to_parquet(out / f"{task['name']}.parquet", index=False)
        print(f"  {task['name']}: {len(task['diseases'])} diseases", flush=True)


if __name__ == "__main__":
    main()
