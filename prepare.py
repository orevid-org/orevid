"""Converts PrimeKG (data/raw) into sparse matrices per relation (data/prep).

Each relation is saved as a symmetric N x N matrix (scipy .npz). If a relation
has subtypes (e.g. drug_protein -> target / enzyme / transporter / carrier),
each subtype is also saved as "relation:subtype".

The external sources of CORE_EXT (data/ext/<name>_nodes|_edges.parquet, sources/*.py) are merged by canonical
key (sources/ids.py): nodes that do not exist in PrimeKG are appended at the
end, and the attributes of the external edges (phase, year, score…) are saved
to data/prep/attributes.parquet.

    py -3 prepare.py
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from sources.ids import primekg_keys

ROOT = Path(__file__).parent
RAW = ROOT / "data" / "raw"
EXT = ROOT / "data" / "ext"
PREP = ROOT / "data" / "prep"
# A fixed list: other tables in data/ext go to the internal graph (layers.py), never here, so that a new source
# cannot add nodes and shift the indices every model and cache relies on.
CORE_EXT = ["ctd", "direction", "gwas_direction", "opentargets", "pubmed", "pubtator", "string", "structure"]


def file_name(rel):
    return "rel_" + re.sub(r"[^a-z0-9]+", "_", rel.lower()).strip("_") + ".npz"


def matrix(x, y, n, values=None):
    x, y = np.concatenate([x, y]), np.concatenate([y, x])
    d = np.ones(len(x), np.int8) if values is None else np.concatenate([values, values])
    m = sp.coo_matrix((d, (x, y)), shape=(n, n)).tocsr()
    if values is None:
        m.data[:] = 1  # duplicated edges must not count twice
    return m


def merge_external(nodes, edges):
    index = {}
    for i, t, nid in zip(nodes.node_index, nodes.node_type, nodes.node_id):
        for c in primekg_keys(t, nid):
            index.setdefault(c, i)

    new = [pd.read_parquet(f) for f in (EXT / f"{s}_nodes.parquet" for s in CORE_EXT) if f.exists()]
    if new:
        new = pd.concat(new).drop_duplicates("key")
        new = new[~new.key.isin(index)].reset_index(drop=True)
        new["node_index"] = np.arange(len(nodes), len(nodes) + len(new))
        index.update(zip(new.key, new.node_index))
        new = new.rename(columns={"key": "node_id"}).assign(origin=new.node_source)
        nodes = pd.concat([nodes, new[nodes.columns]], ignore_index=True)

    external, attributes = [edges], []
    for f in (f for f in (EXT / f"{s}_edges.parquet" for s in CORE_EXT) if f.exists()):  # pubmed: second pass
        a = pd.read_parquet(f)
        a["x_index"] = a.x.map(index)
        a["y_index"] = a.y.map(index)
        lost = a.x_index.isna() | a.y_index.isna()
        a = a[~lost].astype({"x_index": np.int32, "y_index": np.int32})
        print(f"{f.name}: {len(a):,} edges ({lost.sum():,} without a node)")
        external.append(a[["relation", "display_relation", "x_index", "y_index", "source"]])
        attributes.append(a.drop(columns=["x", "y"]))
    if attributes:
        pd.concat(attributes, ignore_index=True).to_parquet(PREP / "attributes.parquet", index=False)
    return nodes, pd.concat(external, ignore_index=True)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    PREP.mkdir(parents=True, exist_ok=True)

    nodes = pd.read_csv(RAW / "nodes.tab", sep="\t", dtype={"node_id": str}, low_memory=False)
    assert (nodes.node_index.to_numpy() == np.arange(len(nodes))).all(), "node_index is not 0..N-1"
    edges = pd.read_csv(
        RAW / "edges.csv",
        dtype={"relation": str, "display_relation": str, "x_index": np.int32, "y_index": np.int32},
    )
    nodes["origin"] = "PrimeKG"
    edges["source"] = "PrimeKG"
    print(f"PrimeKG: {len(nodes):,} nodes, {len(edges) // 2:,} edges")
    nodes, edges = merge_external(nodes, edges)
    edges["relation"] = edges.relation.astype("category")
    edges["display_relation"] = edges.display_relation.astype("category")

    n = len(nodes)
    nodes.to_pickle(PREP / "nodes.pkl")
    print(f"\nnodes: {n:,}")
    print(nodes.groupby(["node_type", "origin"]).size().to_string(), "\n")
    print(f"edges (rows): {len(edges):,}")

    relations = {}

    def save(name, df):
        m = matrix(df.x_index.to_numpy(), df.y_index.to_numpy(), n)
        sp.save_npz(PREP / file_name(name), m)
        relations[name] = {"file": file_name(name), "edges": int(m.nnz // 2)}
        print(f"  {name:45s} {m.nnz // 2:>10,}")

    for rel, df in edges.groupby("relation", observed=True):
        save(rel, df)
        subtypes = df.display_relation.unique()
        if len(subtypes) > 1:
            for sub, dfs in df.groupby("display_relation", observed=True):
                save(f"{rel}:{sub}", dfs)

    # Union graph (to search for paths): the value of each cell is the label code.
    label = (edges.relation.astype(str) + ":" + edges.display_relation.astype(str)).astype("category")
    codes = label.cat.codes.to_numpy().astype(np.int16) + 1
    u = pd.DataFrame({"x": edges.x_index, "y": edges.y_index, "c": codes})
    u = pd.concat([u, u.rename(columns={"x": "y", "y": "x"})]).drop_duplicates(["x", "y"])
    union = sp.coo_matrix((u.c.to_numpy(), (u.x.to_numpy(), u.y.to_numpy())), shape=(n, n)).tocsr()
    sp.save_npz(PREP / "all.npz", union)
    print(f"\nunion graph: {union.nnz // 2:,} undirected edges")

    meta = {"relations": relations, "labels": {str(i + 1): e for i, e in enumerate(label.cat.categories)}}
    (PREP / "relations.json").write_text(json.dumps(meta, ensure_ascii=False, indent=1), "utf-8")

    drug_features = RAW / "drug_features.tab"
    if drug_features.exists():
        pd.read_csv(drug_features, sep="\t", low_memory=False).to_pickle(PREP / "drugs.pkl")
        print("drug descriptions: ok")


if __name__ == "__main__":
    main()
