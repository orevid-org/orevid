"""What has to be done with each gene in each disease: activate it or inhibit it -> data/ext/.

    py -3 -m sources.direction

New relations (disease -> gene, display: activate / inhibit):
    dir_genetic    human variants with direction (Open Targets: gene burden, ClinVar/EVA,
                   Gene2Phenotype, Orphanet, ClinGen, Genomics England, UniProt).
                   Variant that switches the gene off (LoF) and increases risk -> activate;
                   protective LoF -> inhibit; variant that enhances it (GoF) the other way round.
    dir_animal     knockout mice (IMPC): switching the gene off gives a phenotype similar
                   to the disease -> activate.
    dir_expression Expression Atlas: gene up in the disease -> inhibit; down -> activate.
                   It is the weakest: changed expression can be a consequence and not a cause.

When the evidence for a pair is contradictory, the majority wins (weight = Open Targets
score); if they tie, the pair is discarded.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

from sources import ids
from sources.opentargets import Nodes, disease, gene

OT = ids.OT
EXT = Path(__file__).parent.parent / "data" / "ext"
GENETIC = ["evidence_gene_burden", "evidence_eva", "evidence_gene2phenotype", "evidence_orphanet",
           "evidence_clingen", "evidence_genomics_england", "evidence_uniprot_variants"]


def sign(gene_effect, trait_effect):
    """+1 = the gene must be activated, -1 = it must be inhibited, 0 = unknown."""
    g = {"LoF": -1, "GoF": 1}.get(gene_effect, 0)
    t = {"risk": 1, "protect": -1}.get(trait_effect, 0)
    return -g * t


def read(table, columns):
    d = ds.dataset(sorted(str(f) for f in (OT / table).glob("*.parquet")))  # without half-downloaded .part files
    available = [c for c in columns if c in d.schema.names]
    return d.to_table(columns=available).to_pandas()


def genetic():
    parts = []
    for t in GENETIC:
        if not (OT / t).exists():
            print(f"  {t}: not downloaded")
            continue
        d = read(t, ["diseaseId", "targetId", "score", "directionOnTarget", "directionOnTrait"])
        if "directionOnTarget" not in d:
            print(f"  {t}: no direction fields")
            continue
        d["sign"] = [sign(a, b) for a, b in zip(d.directionOnTarget, d.directionOnTrait)]
        n = (d.sign != 0).sum()
        print(f"  {t}: {len(d):,} evidence rows, {n:,} with direction")
        parts.append(d[d.sign != 0].assign(source=t.removeprefix("evidence_")))
    return pd.concat(parts, ignore_index=True).assign(relation="dir_genetic") if parts else None


def animal():
    if not (OT / "evidence_impc").exists() or list((OT / "evidence_impc").glob("*.part")):
        print("  evidence_impc: not downloaded")
        return None
    d = read("evidence_impc", ["diseaseId", "targetId", "score"])
    print(f"  evidence_impc: {len(d):,} evidence rows (knockout -> activate)")
    return d.assign(sign=1, source="impc", relation="dir_animal")


def expression():
    d = read("evidence_expression_atlas", ["diseaseId", "targetId", "score", "log2FoldChangeValue"])
    d = d[d.log2FoldChangeValue.abs() >= 1]
    d["sign"] = -np.sign(d.log2FoldChangeValue).astype(int)
    print(f"  evidence_expression_atlas: {len(d):,} evidence rows with |log2FC| >= 1")
    return d.assign(source="expression_atlas", relation="dir_expression")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    nodes = Nodes()
    parts = [p for p in (genetic(), animal(), expression()) if p is not None]
    d = pd.concat(parts, ignore_index=True)
    dis = {x: disease(x, nodes) for x in d.diseaseId.unique()}
    gene_keys = {x: gene(x, nodes) for x in d.targetId.unique()}
    d["x"] = d.diseaseId.map(dis)
    d["y"] = d.targetId.map(gene_keys)
    d = d.dropna(subset=["x", "y"])
    d = d[d.x.str.startswith(("MONDO:", "OT:"))]  # only diseases (not HPO phenotypes)

    # Weighted vote per pair and evidence type
    d["vote"] = d.sign * d.score.fillna(0.1)
    v = d.groupby(["x", "y", "relation"]).agg(vote=("vote", "sum"), n=("vote", "size"),
                                               sources=("source", lambda s: ",".join(sorted(set(s))))).reset_index()
    v = v[v.vote != 0]
    v["display_relation"] = np.where(v.vote > 0, "activate", "inhibit")
    v["score"] = v.vote.abs()
    edges = v[["x", "y", "relation", "display_relation", "score", "n", "sources"]].assign(source="Open Targets (direction)")

    EXT.mkdir(parents=True, exist_ok=True)
    edges.to_parquet(EXT / "direction_edges.parquet", index=False)
    nodes.table().to_parquet(EXT / "direction_nodes.parquet", index=False)
    print(f"\nEdges: {len(edges):,}")
    print(edges.groupby(["relation", "display_relation"]).size().to_string())


if __name__ == "__main__":
    main()
