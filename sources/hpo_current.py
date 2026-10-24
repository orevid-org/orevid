"""Current HPO disease annotations (phenotype.hpoa 2026-09-02; HPO licence, attribution) for N5 (docs/IMPROVEMENT_PLAN.md,
«Registered 24-10-2026: N5–N7»). OMIM / ORPHA ids -> MONDO (SSSOM exact) -> graph disease; HPO id -> graph phenotype.
Negated (NOT) annotations and non-phenotype aspects are left out.

    py -3 -m sources.hpo_current   # -> data/ext/hpo_current.parquet (disease, phenotype, frequency, source)
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
    hpo = {int(i): int(k) for k, i, t in zip(n.node_index, n.node_id, n.node_type) if t == "effect/phenotype" and str(i).isdigit()}
    h = pd.read_csv(ROOT / "data" / "raw" / "hpo" / "phenotype.hpoa", sep="\t", comment="#", low_memory=False)
    h = h[(h.qualifier != "NOT") & (h.aspect == "P")]
    rows = []
    for db, hp, freq in zip(h.database_id, h.hpo_id, h.frequency):
        pre, _, num = db.partition(":")
        key = {"OMIM": f"OMIM:{num}", "ORPHA": f"ORPHANET:{num}"}.get(pre)
        ph = hpo.get(int(hp.split(":")[1]))
        if key and ph is not None:
            for m in to_mondo.get(key, ()):
                if m in mondo_node:
                    rows.append({"disease": mondo_node[m], "phenotype": ph, "frequency": freq, "source": db})
    out = pd.DataFrame(rows).drop_duplicates(["disease", "phenotype"])
    out.to_parquet(ROOT / "data" / "ext" / "hpo_current.parquet", index=False)
    print(f"HPO current: {len(out):,} disease–phenotype pairs ({out.disease.nunique():,} diseases)")


if __name__ == "__main__":
    build()
