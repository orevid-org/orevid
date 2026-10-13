"""Disease hierarchy with direction (parent -> child) between graph nodes, for the website's
navigation between a disease and its broader and more specific forms (docs/REVIEW_RESPONSE.md, T6).

PrimeKG's `disease_disease` edges are undirected. The direction comes from the Open Targets disease
table (`parents`, MONDO and EFO ids), mapped to graph nodes with the same keys as the rest of the
graph (`sources.ids.diseases`, `sources.ot_approvals.node_index`). Pairs whose ends fall in the same
node (PrimeKG groups several MONDO ids in one node) are dropped, and so is every edge inside a cycle
(cycles appear when several terms map to one node). Display only: it does not enter features, rules or evaluation.

    py -3 -m sources.hierarchy     # -> data/prep/disease_hierarchy.parquet (parent, child)
"""
import sys
from pathlib import Path

import pandas as pd

from sources import ids

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "data" / "prep" / "disease_hierarchy.parquet"


def build():
    from sources.ot_approvals import node_index
    mapping, _, parents, _ = ids.diseases()
    index = node_index()
    nodes = pd.read_pickle(ROOT / "data" / "prep" / "nodes.pkl")
    disease = set(nodes.node_index[nodes.node_type == "disease"].astype(int))

    def node(ot_id):
        k = mapping.get(ot_id)
        n = index.get(k) if k else None
        return n if n in disease else None

    pairs = set()
    for child, ps in parents.items():
        c = node(child)
        if c is None:
            continue
        for p in ps:
            q = node(p)
            if q is not None and q != c:
                pairs.add((q, c))
    df = pd.DataFrame(sorted(pairs), columns=["parent", "child"])
    # Mapping several Open Targets terms to one node can create cycles: drop every edge inside a
    # strongly connected component (both directions of a pair are the smallest case).
    import numpy as np
    import scipy.sparse as sp
    from scipy.sparse.csgraph import connected_components
    n = int(max(df.parent.max(), df.child.max())) + 1
    A = sp.csr_matrix((np.ones(len(df)), (df.parent, df.child)), shape=(n, n))
    _, comp = connected_components(A, directed=True, connection="strong")
    df = df[comp[df.parent.to_numpy()] != comp[df.child.to_numpy()]].reset_index(drop=True)
    DEST.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(DEST, index=False)
    return df


def load():
    """(parent, child) pairs; built on first use."""
    return pd.read_parquet(DEST) if DEST.exists() else build()


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    df = build()
    print(f"{len(df):,} parent -> child pairs between {len(set(df.parent) | set(df.child)):,} disease nodes -> {DEST}")


if __name__ == "__main__":
    main()
