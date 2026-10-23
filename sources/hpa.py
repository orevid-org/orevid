"""Gene attributes for target prioritisation (L4, docs/IMPROVEMENT_PLAN.md, «Registered 24-10-2026: L1–L6»): Human
Protein Atlas (CC BY-SA 4.0; proteinatlas.tsv) tissue specificity, protein class, subcellular and secretome location;
target development level (Tclin / Tchem / Tbio / Tdark) from TCRD / Pharos, as distributed in DrugCentral
(`target_component.tdl`, data/ext/drugcentral_tdl.parquet). Genes -> graph by symbol. Also the gene–tissue layer
(tissues where the RNA is enriched).

    py -3 -m sources.hpa         # -> data/ext/gene_attributes.parquet, data/ext/hpa_tissue.parquet
"""
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "hpa" / "proteinatlas.tsv.zip"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    g = Graph()
    sym = {str(s).upper(): int(k) for k, s, t in zip(g.nodes.node_index, g.nodes.node_name, g.nodes.node_type)
           if t == "gene/protein"}
    z = zipfile.ZipFile(RAW)
    h = pd.read_csv(z.open("proteinatlas.tsv"), sep="\t", low_memory=False,
                    usecols=["Gene", "Ensembl", "Protein class", "RNA tissue specificity", "RNA tissue distribution",
                             "RNA tissue specific nTPM", "Subcellular location", "Secretome location"])
    h["node"] = h.Gene.str.upper().map(sym)
    tdl = pd.read_parquet(ROOT / "data" / "ext" / "drugcentral_tdl.parquet")
    h = h.merge(tdl.rename(columns={"gene": "Gene"})[["Gene", "tdl"]], on="Gene", how="outer")
    h["node"] = h.Gene.str.upper().map(sym)
    attrs = h.dropna(subset=["node"]).astype({"node": int}).drop_duplicates("node")
    attrs.to_parquet(ROOT / "data" / "ext" / "gene_attributes.parquet", index=False)
    rows = []
    for node, spec in zip(attrs.node, attrs["RNA tissue specific nTPM"]):
        if isinstance(spec, str):
            for part in spec.split(";"):
                tissue, _, ntpm = part.partition(":")
                rows.append({"gene": node, "tissue": tissue.strip(), "ntpm": pd.to_numeric(ntpm.strip(), errors="coerce")})
    tissue = pd.DataFrame(rows)
    tissue.to_parquet(ROOT / "data" / "ext" / "hpa_tissue.parquet", index=False)
    print(f"Gene attributes: {len(attrs):,} graph genes (TDL {attrs.tdl.notna().sum():,}: {attrs.tdl.value_counts().to_dict()}; "
          f"tissue specificity {attrs['RNA tissue specificity'].notna().sum():,}); tissue-enriched gene–tissue pairs {len(tissue):,} "
          f"({tissue.tissue.nunique()} tissues)")


if __name__ == "__main__":
    build()
