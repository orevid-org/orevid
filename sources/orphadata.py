"""Orphanet point prevalence (Orphadata product 9, CC BY 4.0) for R4 (docs/IMPROVEMENT_PLAN.md, «Registered
23-10-2026, evening»). ORPHA code -> MONDO (SSSOM exact match) -> graph disease node; the worldwide point-prevalence
class (validated) if there is one, else the first point-prevalence class. R4b: only graph nodes with one MONDO id;
only known classes are kept.

    py -3 -m sources.orphadata   # -> data/ext/orphadata_prevalence.parquet (disease node, ORPHA, class, name)
"""
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "orphadata" / "en_product9_prev.xml"
DEST = ROOT / "data" / "ext" / "orphadata_prevalence.parquet"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.ids import xref_to_mondo
    g = Graph()
    to_mondo = xref_to_mondo()
    mondo_node = {}
    for k, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type):
        if t == "disease" and str(i).isdigit():  # R4b: single-MONDO nodes only (grouped nodes took a subtype's value)
            mondo_node.setdefault(f"MONDO:{int(i)}", int(k))
    rows, n_dis = [], 0
    for d in ET.parse(RAW).getroot().iter("Disorder"):
        code, name = d.findtext("OrphaCode"), d.findtext("Name")
        if not code:
            continue
        n_dis += 1
        best = None
        for p in d.iter("Prevalence"):
            if p.findtext("PrevalenceType/Name") != "Point prevalence" or p.find("PrevalenceClass/Name") is None:
                continue
            item = (p.findtext("PrevalenceGeographic/Name") == "Worldwide",
                    p.findtext("PrevalenceValidationStatus/Name") == "Validated", p.findtext("PrevalenceClass/Name"))
            if best is None or item[:2] > best[:2]:
                best = item
        if best is None or best[2] == "Unknown":
            continue
        nodes = {mondo_node[m] for m in to_mondo.get(f"ORPHANET:{code}", ()) if m in mondo_node}
        if len(nodes) == 1:
            rows.append({"disease": nodes.pop(), "orpha": code, "orpha_name": name, "prevalence": best[2],
                         "worldwide": best[0], "validated": best[1]})
    out = pd.DataFrame(rows).drop_duplicates("disease")
    out.to_parquet(DEST, index=False)
    print(f"Orphadata: {n_dis:,} disorders, {len(out):,} graph diseases with a point-prevalence class -> {DEST}")
    return out


if __name__ == "__main__":
    build()
