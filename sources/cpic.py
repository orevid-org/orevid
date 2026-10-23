"""CPIC gene–drug pairs (CC0; api.cpicpgx.org, files in data/raw/cpic/) for X1 (docs/IMPROVEMENT_PLAN.md,
«Registered 23-10-2026, evening»). Drugs -> graph by the DrugBank id CPIC gives; genes by symbol.

    py -3 -m sources.cpic        # -> data/ext/cpic.parquet (drug node, gene node, CPIC level, guideline)
"""
import json
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "cpic"
DEST = ROOT / "data" / "ext" / "cpic.parquet"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    g = Graph()
    n = g.nodes
    drug_node = {str(i): int(k) for k, i, t, s in zip(n.node_index, n.node_id, n.node_type, n.node_source)
                 if t == "drug" and s == "DrugBank"}
    gene_node = {str(s).upper(): int(k) for k, s, t in zip(n.node_index, n.node_name, n.node_type) if t == "gene/protein"}
    drugs = {d["drugid"]: d for d in json.loads((RAW / "drugs.json").read_text("utf-8"))}
    guide = {x["id"]: x for x in json.loads((RAW / "guidelines.json").read_text("utf-8"))}
    rows, lost = [], 0
    for p in json.loads((RAW / "pairs.json").read_text("utf-8")):
        if p.get("removed"):
            continue
        d = drugs.get(p["drugid"], {})
        f, gn = drug_node.get(d.get("drugbankid") or ""), gene_node.get(p["genesymbol"].upper())
        if f is None or gn is None:
            lost += 1
            continue
        gl = guide.get(p.get("guidelineid"), {})
        rows.append({"drug": f, "gene": gn, "drug_name": d.get("name"), "gene_symbol": p["genesymbol"],
                     "cpic_level": p["cpiclevel"], "guideline": gl.get("name", ""), "url": gl.get("url", "")})
    out = pd.DataFrame(rows)
    out.to_parquet(DEST, index=False)
    print(f"CPIC: {len(out):,} gene–drug pairs mapped ({lost} not mapped); levels {out.cpic_level.value_counts().to_dict()}; "
          f"with a guideline (A/B): {out.cpic_level.isin(['A', 'B', 'A/B']).sum()} -> {DEST}")
    return out


if __name__ == "__main__":
    build()
