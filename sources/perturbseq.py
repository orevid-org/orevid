"""Knockdown signatures from the genome-wide Perturb-seq of Replogle et al. 2022 (K562, normalised pseudo-bulk;
figshare 20029387, CC BY 4.0) for P2 (docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, evening: P2»).
Per knocked-down gene: the 100 most increased and 100 most decreased measured genes (z-scores), genes -> graph by symbol.

    py -3 -m sources.perturbseq    # -> data/ext/perturbseq.parquet (knocked-down gene node, gene node, up / down)
"""
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "perturbseq" / "K562_gwps_normalized_bulk_01.h5ad"
DEST = ROOT / "data" / "ext" / "perturbseq.parquet"
TOP = 100


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    g = Graph()
    sym = {str(s).upper(): int(k) for k, s, t in zip(g.nodes.node_index, g.nodes.node_name, g.nodes.node_type)
           if t == "gene/protein"}
    with h5py.File(RAW, "r") as f:
        X = f["X"][:]
        pert = [x.decode() for x in f["obs"]["gene_transcript"][:]]
        cats = [x.decode() for x in f["var"]["__categories"]["gene_name"][:]]
        names = [cats[i] if i >= 0 else "" for i in f["var"]["gene_name"][:]]
    measured = np.array([sym.get(n.upper(), -1) for n in names])
    rows = []
    for i, p in enumerate(pert):
        gene = p.split("_")[1] if "_" in p else p
        kd = sym.get(gene.upper())
        if kd is None or gene.lower().startswith("non-targeting"):
            continue
        z = X[i]
        order = np.argsort(z)
        for direction, idx in (("down", order[:TOP]), ("up", order[::-1][:TOP])):
            rows += [(kd, int(measured[j]), direction, float(z[j])) for j in idx if measured[j] >= 0]
    out = pd.DataFrame(rows, columns=["knockdown", "gene", "direction", "z"])
    out.to_parquet(DEST, index=False)
    print(f"Perturb-seq: {len(pert):,} perturbations, {out.knockdown.nunique():,} knocked-down genes in the graph, "
          f"{len(out):,} signature rows -> {DEST}")
    return out


if __name__ == "__main__":
    build()
