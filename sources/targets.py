"""How each protein can be attacked: Open Targets 26.06 "target" table (CC0).

For each gene in the graph (data/ext/targets.parquet):
- tractability per modality (SM small molecule, AB antibody, PR PROTAC, OC others:
  oligonucleotides…), summarised as 3 = clinical (there are already drugs in the clinic), 2 = evidence
  (structure with ligand, high-quality pocket, accessible location…), 1 = predicted
  (druggable family, medium pocket…), 0 = nothing;
- known safety risks (number and names), genetic constraint (LOEUF: the lower,
  the less it tolerates being switched off), chemical probes, protein class and whether it is in the membrane or
  secreted.

    py -3 -m sources.targets
"""
from pathlib import Path

import numpy as np
import pandas as pd

from sources.ids import OT, genes, primekg_keys

ROOT = Path(__file__).parent.parent
DEST = ROOT / "data" / "ext" / "targets.parquet"

CLINICAL = {"Approved Drug", "Advanced Clinical", "Phase 1 Clinical"}
EVIDENCE = {"Structure with Ligand", "High-Quality Ligand", "High-Quality Pocket", "UniProt loc high conf",
            "GO CC high conf", "Literature", "UniProt Ubiquitination", "Database Ubiquitination",
            "Small Molecule Binder"}
MODALITIES = {"SM": "small molecule", "AB": "antibody", "PR": "PROTAC", "OC": "other (oligonucleotide…)"}


def level(ids):
    if ids & CLINICAL:
        return 3
    if ids & EVIDENCE:
        return 2
    return 1 if ids else 0


# Open Targets prioritisation factors (target_prioritisation). Without maxClinicalStage: it comes from the
# drugs (leakage in the target evaluation).
PRIORITY = ["isInMembrane", "isSecreted", "hasPocket", "hasLigand", "hasSmallMoleculeBinder", "geneticConstraint",
            "paralogMaxIdentityPercentage", "mouseOrthologMaxIdentityPercentage", "isCancerDriverGene", "mouseKOScore",
            "hasHighQualityChemicalProbes", "tissueSpecificity", "tissueDistribution", "celltypeSpecificity",
            "celltypeDistribution"]


def essentiality():
    """DepMap via Open Targets (target_essentiality): per gene, whether it is common essential, the mean
    effect of knocking it out (Chronos; < 0 = the cells depend on it) and the fraction of cell
    lines that depend on it (effect < -0.5)."""
    import pyarrow.compute as pc
    import pyarrow.parquet as pq
    rows = []
    for f in sorted((OT / "target_essentiality").glob("*.parquet")):
        t = pq.read_table(f)
        ids = t.column("id").combine_chunks()
        ge = t.column("geneEssentiality").combine_chunks()
        g1 = pc.list_flatten(ge)
        p1 = pc.list_parent_indices(ge)
        essential = g1.field("isEssential")
        dm = g1.field("depMapEssentiality")
        g2 = pc.list_flatten(dm)
        p2 = pc.list_parent_indices(dm)
        sc = g2.field("screens")
        g3 = pc.list_flatten(sc)
        p3 = pc.list_parent_indices(sc)
        effect = g3.field("geneEffect").to_numpy(zero_copy_only=False)
        gene_of_screen = p1.to_numpy()[p2.to_numpy()[p3.to_numpy()]]
        d = pd.DataFrame({"gene": gene_of_screen, "effect": effect})
        r = d.groupby("gene").effect.agg(depmap_effect="mean", dependent_fraction=lambda v: np.mean(v < -0.5))
        ess = pd.Series(essential.to_numpy(zero_copy_only=False), index=p1.to_numpy()).groupby(level=0).max()
        r["essential"] = ess.reindex(r.index).fillna(False).astype(int)
        r.index = ids.to_numpy(zero_copy_only=False)[r.index]
        rows.append(r)
    return pd.concat(rows)


def main():
    t = pd.read_parquet(OT / "target", columns=["id", "approvedSymbol", "tractability", "safetyLiabilities",
                                                "constraint", "chemicalProbes", "targetClass", "subcellularLocations",
                                                "proteinIds"])
    ens = genes()
    nodes = pd.read_pickle(ROOT / "data" / "prep" / "nodes.pkl")
    index = {}
    for i, tp, nid in zip(nodes.node_index, nodes.node_type, nodes.node_id):
        if tp != "gene/protein":
            continue
        for c in (primekg_keys(tp, nid) if not str(nid).startswith("NCBI:") else [nid]):
            index.setdefault(c, i)
    rows = []
    for r in t.itertuples():
        key = ens.get(r.id, (None,))[0]
        if key not in index:
            continue
        f = {"node_index": index[key], "ensembl": r.id, "symbol": r.approvedSymbol}
        prot = r.proteinIds if r.proteinIds is not None else []
        f["uniprot"] = next((x["id"] for x in prot if x["source"] == "uniprot_swissprot"),
                            next((x["id"] for x in prot if "uniprot" in str(x["source"])), ""))
        tr = r.tractability if r.tractability is not None else []
        for mod in MODALITIES:
            ids = {x["id"] for x in tr if x["modality"] == mod and x["value"]}
            f[f"tract_{mod}"] = level(ids)
        liabilities = r.safetyLiabilities if r.safetyLiabilities is not None else []
        events = sorted({s["event"] for s in liabilities if s["event"]})
        f["n_safety"] = len(liabilities)
        f["safety"] = "; ".join(events[:5])
        c = [x for x in (r.constraint if r.constraint is not None else []) if x["constraintType"] == "lof"]
        f["loeuf"] = c[0]["oeUpper"] if c and c[0]["oeUpper"] is not None else np.nan
        f["n_probes"] = len(r.chemicalProbes) if r.chemicalProbes is not None else 0
        cl = r.targetClass if r.targetClass is not None else []
        f["target_class"] = next((x["label"] for x in cl if x["level"] == "l1"), "")
        loc = " ".join(str(x["location"]).lower() for x in (r.subcellularLocations if r.subcellularLocations is not None else []))
        f["membrane_or_secreted"] = int("membrane" in loc or "secreted" in loc or "extracellular" in loc)
        rows.append(f)
    d = pd.DataFrame(rows).drop_duplicates("node_index")
    pr = pd.read_parquet(OT / "target_prioritisation").set_index("targetId")[PRIORITY]
    d = d.join(pr, on="ensembl")
    e = essentiality()
    d = d.join(e, on="ensembl")
    print(f"  with DepMap data: {d.depmap_effect.notna().mean():.1%}; common essential: {d.essential.fillna(0).mean():.1%}")
    DEST.parent.mkdir(parents=True, exist_ok=True)
    d.to_parquet(DEST)
    print(f"{len(d)} genes with target data -> {DEST}")
    for mod, name in MODALITIES.items():
        print(f"  {name:24s} clinical {np.mean(d[f'tract_{mod}'] == 3):5.1%}  evidence {np.mean(d[f'tract_{mod}'] == 2):5.1%}")
    print(f"  with safety risks: {np.mean(d.n_safety > 0):.1%};  with LOEUF: {d.loeuf.notna().mean():.1%}")


if __name__ == "__main__":
    main()
