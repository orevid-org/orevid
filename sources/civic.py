"""CIViC clinical evidence (CC0; nightly ClinicalEvidenceSummaries.tsv) for O1 (docs/IMPROVEMENT_PLAN.md, «Registered
23-10-2026, evening»). All accepted evidence items go to the internal graph; O1 uses the predictive, supporting,
sensitivity / response items of levels A–B with a single drug. Drugs -> graph by name (moieties, synonyms); diseases
by DOID -> MONDO (SSSOM exact) -> graph node; genes by the molecular profile's first word (gene symbol).

    py -3 -m sources.civic       # -> data/ext/civic.parquet
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "civic" / "ClinicalEvidenceSummaries.tsv"
DEST = ROOT / "data" / "ext" / "civic.parquet"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.ids import xref_to_mondo
    from sources.practical import moiety
    from sources.prism import drug_index
    g = Graph()
    to_mondo = xref_to_mondo()
    mondo_node = {}
    for k, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type):
        if t == "disease" and str(i).replace("_", "").isdigit():
            for p in str(i).split("_"):
                mondo_node.setdefault(f"MONDO:{int(p)}", int(k))
    gene_node = {str(s).upper(): int(k) for k, s, t in zip(g.nodes.node_index, g.nodes.node_name, g.nodes.node_type)
                 if t == "gene/protein"}
    by_name = drug_index(g)
    c = pd.read_csv(RAW, sep="\t", low_memory=False)
    c = c[c.evidence_status == "accepted"]
    rows = []
    for r in c.itertuples():
        doid = f"DOID:{int(r.doid)}" if r.doid == r.doid else ""
        ds = {mondo_node[m] for m in to_mondo.get(doid, ()) if m in mondo_node}
        d = ds.pop() if len(ds) == 1 else None
        therapies = [x.strip() for x in str(r.therapies).split(",")] if r.therapies == r.therapies else []
        gene = gene_node.get(str(r.molecular_profile).split(" ")[0].upper())
        for th in therapies or [None]:
            fs = by_name.get(moiety(th), set()) if th else {None}
            for f in fs or {None}:
                rows.append({"evidence_id": int(r.evidence_id), "drug": f, "disease": d, "gene": gene,
                             "therapy": th, "n_therapies": len(therapies), "interaction": r.therapy_interaction_type,
                             "disease_name": r.disease, "molecular_profile": r.molecular_profile,
                             "evidence_type": r.evidence_type, "direction": r.evidence_direction,
                             "significance": r.significance, "level": r.evidence_level, "url": r.evidence_civic_url})
    out = pd.DataFrame(rows)
    out.to_parquet(DEST, index=False)
    o1 = out[(out.evidence_type == "Predictive") & (out.direction == "Supports") & (out.significance == "Sensitivity/Response")
             & out.level.isin(["A", "B"]) & (out.n_therapies == 1)]
    print(f"CIViC: {c.shape[0]:,} accepted evidence items, {len(out):,} rows; O1 items {o1.evidence_id.nunique():,}, "
          f"with drug and disease mapped {o1.dropna(subset=['drug', 'disease']).evidence_id.nunique():,} -> {DEST}")
    return out


if __name__ == "__main__":
    build()
