"""DepMap CRISPR gene effect (Public 24Q2, figshare 25880521, CC BY 4.0) for L5 (docs/IMPROVEMENT_PLAN.md, «Registered
24-10-2026: L1–L6»). Per gene and Oncotree lineage the median Chronos effect over cell lines (more negative = the
lineage depends more on the gene); genes -> graph by Entrez id (the column header «SYMBOL (entrez)»).

    py -3 -m sources.depmap      # -> data/ext/depmap_lineage.parquet (gene node, lineage, median effect, cell lines)
"""
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "depmap"
DEST = ROOT / "data" / "ext" / "depmap_lineage.parquet"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    g = Graph()
    entrez = {str(i): int(k) for k, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type) if t == "gene/protein"}
    m = pd.read_csv(RAW / "Model.csv", usecols=["ModelID", "OncotreeLineage"])
    e = pd.read_csv(RAW / "CRISPRGeneEffect.csv", index_col=0)
    cols = {c: entrez.get(m_.group(1)) for c in e.columns if (m_ := re.search(r"\((\d+)\)", c))}
    e = e[[c for c, v in cols.items() if v is not None]]
    e.columns = [cols[c] for c in e.columns]
    lin = m.set_index("ModelID").OncotreeLineage.reindex(e.index)
    med = e.groupby(lin.values).median().T
    n = e.groupby(lin.values).size()
    out = med.stack().reset_index()
    out.columns = ["gene", "lineage", "effect"]
    out["cell_lines"] = out.lineage.map(n)
    out.to_parquet(DEST, index=False)
    print(f"DepMap: {e.shape[0]:,} cell lines, {e.shape[1]:,} graph genes, {out.lineage.nunique()} lineages -> {DEST}")


if __name__ == "__main__":
    build()
