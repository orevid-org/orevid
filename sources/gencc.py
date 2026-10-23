"""GenCC gene–disease validity (G1, docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026»; CC0 1.0, attribution requested).

Mapping and classification fixed in the registration:
- genes: GenCC `gene_curie` (HGNC) -> NCBI Gene (HGNC file in data/raw/hgnc/) -> graph gene node;
- diseases: GenCC `disease_curie` (MONDO) -> the graph disease node whose MONDO group contains it; no propagation;
- per (gene, disease): «supported» = at least one Definitive, Strong or Moderate and no ClinGen Disputed or Refuted;
  «contradicted» = no Definitive / Strong / Moderate and at least one Disputed, Refuted or «No Known Disease
  Relationship»; «limited» = the rest (Limited, Supportive, Animal model only).

    py -3 -m sources.gencc        # -> data/ext/gencc.parquet (gene node, disease node, class)
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "gencc" / "gencc-submissions.tsv"
DEST = ROOT / "data" / "ext" / "gencc.parquet"
STRONG = {"Definitive", "Strong", "Moderate"}
AGAINST = {"Disputed Evidence", "Refuted Evidence", "No Known Disease Relationship"}


def classify(group):
    cls = set(group.classification_title)
    clingen_against = bool(((group.submitter_title == "ClinGen") & group.classification_title.isin(
        {"Disputed Evidence", "Refuted Evidence"})).any())
    if cls & STRONG and not clingen_against:
        return "supported"
    if not cls & STRONG and cls & AGAINST:
        return "contradicted"
    return "limited"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    g = Graph()
    t = pd.read_csv(RAW, sep="\t", dtype=str, low_memory=False,
                    usecols=["gene_curie", "disease_curie", "classification_title", "submitter_title"])
    h = pd.read_csv(ROOT / "data" / "raw" / "hgnc" / "hgnc_complete_set.txt", sep="\t", dtype=str,
                    usecols=["hgnc_id", "entrez_id"], low_memory=False)
    hgnc_entrez = dict(zip(h.hgnc_id, h.entrez_id))
    gene_node = {str(i): int(n) for n, i, ty in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type)
                 if ty == "gene/protein"}
    mondo_node = {}
    for n, i, ty in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type):
        if ty == "disease" and str(i).replace("_", "").isdigit():
            for p in str(i).split("_"):
                mondo_node.setdefault(int(p), int(n))
    t["gene"] = t.gene_curie.map(lambda c: gene_node.get(str(hgnc_entrez.get(c, "")).split(".")[0]))
    t["disease"] = t.disease_curie.map(lambda c: mondo_node.get(int(c.split(":")[1])) if str(c).startswith("MONDO:") else None)
    n_sub, n_pairs = len(t), t.groupby(["gene_curie", "disease_curie"]).ngroups
    m = t.dropna(subset=["gene", "disease"])
    rows = [{"gene": int(gn), "disease": int(dn), "class": classify(grp), "submissions": len(grp)}
            for (gn, dn), grp in m.groupby(["gene", "disease"])]
    out = pd.DataFrame(rows)
    DEST.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(DEST, index=False)
    print(f"{n_sub:,} submissions, {n_pairs:,} GenCC gene–disease pairs; genes mapped {t.gene.notna().mean():.1%}, "
          f"diseases mapped {t.disease.notna().mean():.1%}; {len(out):,} graph pairs "
          f"({out.gene.nunique():,} genes, {out.disease.nunique():,} diseases)")
    print(out["class"].value_counts().to_string())
    sup = out[out["class"] == "supported"].groupby("disease").size()
    print(f"covered graph diseases (≥ 3 supported genes): {(sup >= 3).sum():,}")
    return out


if __name__ == "__main__":
    build()
