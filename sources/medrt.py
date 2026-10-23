"""MED-RT (VA / NLM, public domain; Core_MEDRT_XML.zip, release 2026.10.05) for L1 (docs/IMPROVEMENT_PLAN.md,
«Registered 24-10-2026: L1–L6»): RxNorm drug -> MeSH disease associations may_treat, may_prevent, CI_with, induces,
may_diagnose. Drugs -> graph by name (moieties, synonyms); MED-RT gives MeSH *concept* ids (M…), turned into descriptor ids
(D…) with the NLM descriptor file (data/raw/mesh/desc2026.gz, concept_to_descriptor.json), then MeSH -> MONDO (SSSOM
exact) -> graph node.

    py -3 -m sources.medrt       # -> data/ext/medrt.parquet
"""
import json
import re
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "medrt" / "Core_MEDRT_XML.zip"
DEST = ROOT / "data" / "ext" / "medrt.parquet"
KEEP = {"may_treat", "may_prevent", "CI_with", "induces", "may_diagnose"}
ASSN = re.compile(r"<association>\s*<namespace>[^<]*</namespace>\s*<name>([^<]+)</name>\s*<from_namespace>RxNorm</from_namespace>"
                  r"\s*<from_name>([^<]+)</from_name>\s*<from_code>[^<]*</from_code>\s*<to_namespace>MeSH</to_namespace>"
                  r"\s*<to_name>([^<]+)</to_name>\s*<to_code>([^<]+)</to_code>")


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.ids import xref_to_mondo
    from sources.practical import moiety
    from sources.prism import drug_index
    g = Graph()
    z = zipfile.ZipFile(RAW)
    xml = z.read(next(n for n in z.namelist() if n.endswith(".xml"))).decode("utf-8", "replace")
    by_name = drug_index(g)
    to_mondo = xref_to_mondo()
    mondo_node = {}
    for k, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type):
        if t == "disease" and str(i).replace("_", "").isdigit():
            for p in str(i).split("_"):
                mondo_node.setdefault(f"MONDO:{int(p)}", int(k))
    c2d = json.loads((ROOT / "data" / "raw" / "mesh" / "concept_to_descriptor.json").read_text())
    rows, n = [], 0
    for rel, drug, dname, mesh in ASSN.findall(xml):
        if rel not in KEEP:
            continue
        n += 1
        fs = by_name.get(moiety(drug), set())
        ds = {mondo_node[m] for m in to_mondo.get(f"MESH:{c2d.get(mesh, mesh)}", ()) if m in mondo_node}
        for f in fs:
            for d in ds:
                rows.append({"drug": f, "disease": d, "relation": rel, "rxnorm_name": drug, "mesh": mesh, "mesh_name": dname})
    out = pd.DataFrame(rows).drop_duplicates(["drug", "disease", "relation"])
    out.to_parquet(DEST, index=False)
    print(f"MED-RT: {n:,} drug–disease associations kept; {len(out):,} graph pairs {out.relation.value_counts().to_dict()} "
          f"({out.drug.nunique():,} drugs, {out.disease.nunique():,} diseases) -> {DEST}")


if __name__ == "__main__":
    build()
