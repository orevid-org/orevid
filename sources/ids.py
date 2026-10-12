"""Identifier mappings between sources.

Canonical node keys (prepare.py uses them for merging):
    drug       DB:DB00001        (DrugBank) or CHEMBL:CHEMBL25 if it has none
    gene       NCBI:1234         (Entrez)
    disease    MONDO:5301        (MONDO number without leading zeros)
    phenotype  HPO:1250
    others     OT:EFO_0000180, MEDDRA:10019211…
"""
from collections import Counter, defaultdict
from functools import cache
from pathlib import Path

import pandas as pd

RAW = Path(__file__).parent.parent / "data" / "raw"
OT = RAW / "opentargets"


def primekg_keys(node_type, node_id):
    """Keys of a PrimeKG node (a grouped disease has several)."""
    if node_type == "drug":
        return [f"DB:{node_id}"]
    if node_type == "gene/protein":
        return [f"NCBI:{node_id}"]
    if node_type == "disease":
        return [f"MONDO:{int(p)}" for p in str(node_id).split("_")]
    if node_type == "effect/phenotype":
        return [f"HPO:{int(node_id)}"]
    return [f"{node_type}:{node_id}"]


def _xref(x):
    """Normalises "mesh:D1" / "MESH:D1" / "Orphanet:5" to an upper-case prefix."""
    pre, _, val = str(x).partition(":")
    return f"{pre.upper()}:{val}"


@cache
def genes():
    """Ensembl (ENSG…) -> (NCBI key, symbol), via HGNC."""
    h = pd.read_csv(RAW / "hgnc" / "hgnc_complete_set.txt", sep="\t", low_memory=False,
                    usecols=["symbol", "entrez_id", "ensembl_gene_id"])
    h = h.dropna(subset=["entrez_id", "ensembl_gene_id"])
    return {e: (f"NCBI:{int(n)}", s) for e, n, s in zip(h.ensembl_gene_id, h.entrez_id, h.symbol)}


@cache
def all_primekg_keys():
    n = pd.read_csv(RAW / "nodes.tab", sep="\t", dtype={"node_id": str}, low_memory=False)
    keys = {}
    for t, i, name in zip(n.node_type, n.node_id, n.node_name):
        for c in primekg_keys(t, i):
            keys[c] = (t, name)
    return keys


@cache
def drugs():
    """Open Targets molecules: ChEMBL id -> (key, name, type, maximum phase)."""
    m = pd.read_parquet(OT / "drug_molecule",
                        columns=["id", "name", "drugType", "parentId", "childChemblIds",
                                 "crossReferences", "maximumClinicalStage"])
    db = {}
    for cid, refs in zip(m.id, m.crossReferences):
        if refs is None:
            continue
        for r in refs:
            if r["source"] == "drugbank" and len(r["ids"]):
                db[cid] = sorted(r["ids"])[0]
    # Salts and forms: if the form has no DrugBank id but the parent (or a child) has one, it inherits it.
    for cid, parent, children in zip(m.id, m.parentId, m.childChemblIds):
        if cid in db:
            continue
        if isinstance(parent, str) and parent in db:
            db[cid] = db[parent]
        elif children is not None:
            for f in children:
                if f in db:
                    db[cid] = db[f]
                    break
    res = {}
    for cid, name, drug_type, phase in zip(m.id, m.name, m.drugType, m.maximumClinicalStage):
        key = f"DB:{db[cid]}" if cid in db else f"CHEMBL:{cid}"
        res[cid] = (key, str(name).capitalize() if isinstance(name, str) else cid, drug_type, phase)
    return res


@cache
def xref_to_mondo():
    """"MESH:D009103", "OMIM:126200", "DOID:2377"… -> set of MONDO keys (via SSSOM)."""
    s = pd.read_csv(RAW / "mondo" / "mondo.sssom.tsv", sep="\t", comment="#")
    s = s[s.predicate_id == "skos:exactMatch"]
    to_mondo = defaultdict(set)
    for mondo, obj in zip(s.subject_id, s.object_id):
        to_mondo[_xref(obj)].add(f"MONDO:{int(mondo.split(':')[1])}")
    return to_mondo


@cache
def diseases():
    """Open Targets disease id -> canonical key (MONDO, HPO or OT:…).

    Also returns the names and parents, to create the nodes that have no MONDO.
    """
    to_mondo = xref_to_mondo()
    d = pd.read_parquet(OT / "disease", columns=["id", "name", "dbXRefs", "parents"])
    mapping, via = {}, Counter()
    for did, refs in zip(d.id, d.dbXRefs):
        pre, _, num = did.partition("_")
        if pre == "MONDO":
            mapping[did] = f"MONDO:{int(num)}"; via["MONDO direct"] += 1
        elif pre == "HP":
            mapping[did] = f"HPO:{int(num)}"; via["HPO direct"] += 1
        elif (m := to_mondo.get(f"{pre.upper()}:{num}")) and len(m) == 1:
            mapping[did] = next(iter(m)); via["EFO->MONDO"] += 1
        else:
            votes = Counter(x for r in (refs if refs is not None else []) for x in to_mondo.get(_xref(r), ()))
            if votes:
                (best, n), *rest = votes.most_common(2) + [(None, 0)]
                if n > rest[0][1]:  # no tie
                    mapping[did] = best; via["xref->MONDO"] += 1
                    continue
            if pre in ("EFO", "Orphanet", "MONDO"):
                mapping[did] = f"OT:{did}"; via["new node"] += 1
            else:
                via[f"discarded ({pre})"] += 1
    names = dict(zip(d.id, d.name))
    parents = {i: list(p) for i, p in zip(d.id, d.parents) if p is not None}
    return mapping, names, parents, via
