"""STRING protein interactions (human) -> data/ext/.

    py -3 -m sources.string_db

New relation: string (gene ↔ gene), with two subtypes according to the STRING combined
confidence: "≥900" (very high) and "700-899" (high). The ENSP identifiers
are mapped to Entrez genes by symbol (STRING preferred_name -> HGNC).
"""
import sys
from pathlib import Path

import pandas as pd

from sources import ids
from sources.download import STRING_VERSION

RAW = ids.RAW / "string"
EXT = Path(__file__).parent.parent / "data" / "ext"
THRESHOLD = 700


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    info = pd.read_csv(RAW / f"9606.protein.info.v{STRING_VERSION}.txt.gz", sep="\t",
                       usecols=[0, 1], names=["ensp", "symbol"], header=0)
    h = pd.read_csv(ids.RAW / "hgnc" / "hgnc_complete_set.txt", sep="\t", low_memory=False,
                    usecols=["symbol", "entrez_id", "prev_symbol", "alias_symbol"])
    h = h.dropna(subset=["entrez_id"])
    by_symbol = {}
    # Old symbols and aliases only if they do not clash with a current symbol.
    for col in ("alias_symbol", "prev_symbol"):
        for s, e in zip(h[col], h.entrez_id):
            if isinstance(s, str):
                for x in s.split("|"):
                    by_symbol[x] = f"NCBI:{int(e)}"
    by_symbol.update({s: f"NCBI:{int(e)}" for s, e in zip(h.symbol, h.entrez_id)})
    gene = {p: by_symbol.get(s) for p, s in zip(info.ensp, info.symbol)}
    name = dict(zip(info.ensp, info.symbol))
    print(f"STRING {STRING_VERSION}: {sum(v is not None for v in gene.values()):,} of {len(gene):,} proteins with a gene")

    rows = []
    for block in pd.read_csv(RAW / f"9606.protein.links.v{STRING_VERSION}.txt.gz", sep=" ", chunksize=2_000_000):
        block = block[(block.combined_score >= THRESHOLD) & (block.protein1 < block.protein2)]  # each pair once
        rows.append(block)
    l = pd.concat(rows, ignore_index=True)
    l["x"] = l.protein1.map(gene)
    l["y"] = l.protein2.map(gene)
    lost = l.x.isna() | l.y.isna()
    l = l[~lost & (l.x != l.y)]
    l["display_relation"] = (l.combined_score >= 900).map({True: "≥900", False: "700-899"})
    l = l.rename(columns={"combined_score": "score"}).sort_values("score", ascending=False).drop_duplicates(["x", "y"])
    edges = l.assign(relation="string", source=f"STRING {STRING_VERSION}")[
        ["x", "y", "relation", "display_relation", "score", "source"]]

    # Genes that are not in PrimeKG (with the symbol as name)
    primekg = ids.all_primekg_keys()
    new = {}
    for p, g in gene.items():
        if g and g not in primekg and g not in new:
            new[g] = name[p]
    used = set(edges.x) | set(edges.y)
    nodes = pd.DataFrame([(g, "gene/protein", n, "STRING") for g, n in new.items() if g in used],
                         columns=["key", "node_type", "node_name", "node_source"])

    EXT.mkdir(parents=True, exist_ok=True)
    edges.to_parquet(EXT / "string_edges.parquet", index=False)
    nodes.to_parquet(EXT / "string_nodes.parquet", index=False)
    print(f"edges: {len(edges):,} ({lost.sum():,} without a gene); new nodes: {len(nodes):,}")
    print(edges.display_relation.value_counts().to_string())


if __name__ == "__main__":
    main()
