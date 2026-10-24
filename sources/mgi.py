"""MGI mouse models of human disease (MGI_DO.rpt, CC BY 4.0) for N6 (docs/IMPROVEMENT_PLAN.md, «Registered 24-10-2026:
N5–N7»): DO disease -> MONDO (SSSOM exact) -> graph disease; mouse gene -> human gene by symbol (upper case).

    py -3 -m sources.mgi           # -> data/ext/mgi_models.parquet (disease, gene)
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.ids import xref_to_mondo
    g = Graph()
    n = g.nodes
    to_mondo = xref_to_mondo()
    mondo_node = {}
    for k, i, t in zip(n.node_index, n.node_id, n.node_type):
        if t == "disease" and str(i).replace("_", "").isdigit():
            for p in str(i).split("_"):
                mondo_node.setdefault(f"MONDO:{int(p)}", int(k))
    sym = {str(s).upper(): int(k) for k, s, t in zip(n.node_index, n.node_name, n.node_type) if t == "gene/protein"}
    m = pd.read_csv(ROOT / "data" / "raw" / "mgi" / "MGI_DO.rpt", sep="\t", dtype=str)
    m = m[m["Common Organism Name"].str.startswith("mouse")]
    rows = [{"disease": mondo_node[mo], "gene": sym[str(s).upper()], "mgi": mgi}
            for doid, s, mgi in zip(m["DO Disease ID"], m.Symbol, m["Mouse MGI ID"])
            for mo in to_mondo.get(doid, ()) if mo in mondo_node and str(s).upper() in sym]
    out = pd.DataFrame(rows).drop_duplicates(["disease", "gene"])
    out.to_parquet(ROOT / "data" / "ext" / "mgi_models.parquet", index=False)
    print(f"MGI: {len(out):,} disease–gene pairs from mouse models ({out.disease.nunique():,} diseases)")


if __name__ == "__main__":
    build()
