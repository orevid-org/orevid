"""Broad Drug Repurposing Hub (20200324; metadata CC BY 4.0 on the website, «non-commercial use only» in the file:
redistribution to check) — annotated targets for P1b (docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, evening»).
Compounds -> graph drugs by InChIKey (Open Targets molecules) and, failing that, by name (moieties, synonyms).

    py -3 -m sources.repurposing_hub   # -> data/ext/repurposing_hub.parquet (drug node, gene node, moa, phase)
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "repurposing_hub"
DEST = ROOT / "data" / "ext" / "repurposing_hub.parquet"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources import ids
    from sources.practical import moiety
    from sources.prism import drug_index
    g = Graph()
    n = g.nodes
    drugs = pd.read_csv(RAW / "repurposing_drugs_20200324.txt", sep="\t", comment="!", dtype=str)
    samples = pd.read_csv(RAW / "repurposing_samples_20200324.txt", sep="\t", comment="!", dtype=str,
                          usecols=["pert_iname", "InChIKey"])
    key_node = {(f"DB:{i}" if s == "DrugBank" else f"CHEMBL:{i}"): int(k)
                for k, i, t, s in zip(n.node_index, n.node_id, n.node_type, n.node_source) if t == "drug"}
    m = pd.read_parquet(ids.OT / "drug_molecule", columns=["id", "inchiKey"])
    mol = ids.drugs()
    inchi_node = {ik: key_node[mol[c][0]] for c, ik in zip(m.id, m.inchiKey) if isinstance(ik, str) and c in mol
                  and mol[c][0] in key_node}
    by_name = drug_index(g)
    gene_node = {str(s).upper(): int(k) for k, s, t in zip(n.node_index, n.node_name, n.node_type) if t == "gene/protein"}
    node_of = {}
    for name, ik in zip(samples.pert_iname, samples.InChIKey):
        if isinstance(ik, str) and ik in inchi_node:
            node_of.setdefault(name, {inchi_node[ik]})
    how = {"InChIKey": 0, "name": 0, "none": 0}
    rows = []
    for name, moa, tgt, phase in zip(drugs.pert_iname, drugs.moa, drugs.target, drugs.clinical_phase):
        fs = node_of.get(name)
        if fs:
            how["InChIKey"] += 1
        else:
            fs = by_name.get(moiety(name), set())
            how["name" if fs else "none"] += 1
        for t in (str(tgt).split("|") if isinstance(tgt, str) else []):
            gn = gene_node.get(t.upper())
            if gn is not None:
                for f in fs:
                    rows.append({"drug": f, "gene": gn, "hub_name": name, "moa": moa, "phase": phase})
    out = pd.DataFrame(rows).drop_duplicates(["drug", "gene"])
    out.to_parquet(DEST, index=False)
    print(f"Repurposing Hub: {len(drugs):,} compounds (mapped by InChIKey {how['InChIKey']:,}, by name {how['name']:,}, "
          f"not mapped {how['none']:,}); {len(out):,} drug–gene pairs ({out.drug.nunique():,} drugs) -> {DEST}")
    return out


if __name__ == "__main__":
    build()
