"""Measured drug targets from ChEMBL activities and the IUPHAR / BPS Guide to Pharmacology (P1, docs/IMPROVEMENT_PLAN.md,
«Registered 23-10-2026: P1»).

ChEMBL 37 (web services): pChEMBL ≥ 6, IC50 / Ki / Kd / EC50, no data-validity comment, human single-protein targets.
GtoPdb (interactions.csv): human targets, ligands through their ChEMBL ids, affinity (median, else the middle of the
range) ≥ 6. Targets -> UniProt -> HGNC -> graph gene; drugs -> graph through sources/ids.py. Raw files are cached in
data/raw/chembl/ and data/raw/gtopdb/; only the derived table is used.

    py -3 -m sources.offtargets       # -> data/ext/measured_targets.parquet (drug node, gene node, source, best pX)
"""
import io
import json
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
CHEMBL = ROOT / "data" / "raw" / "chembl"
GTOP = ROOT / "data" / "raw" / "gtopdb"
DEST = ROOT / "data" / "ext" / "measured_targets.parquet"
API = "https://www.ebi.ac.uk/chembl/api/data"
TYPES = "IC50,Ki,Kd,EC50"


def get(url, tries=5):
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return json.load(r)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))


def pages(url, key):
    out = []
    while url:
        j = get(url)
        out += j[key]
        url = API.rsplit("/chembl", 1)[0] + j["page_meta"]["next"] if j["page_meta"]["next"] else None
    return out


def chembl_activities(mols):
    f = CHEMBL / "activities_p6.parquet"
    if f.exists():
        return pd.read_parquet(f)

    def batch(ids):
        q = urllib.parse.urlencode({"molecule_chembl_id__in": ",".join(ids), "pchembl_value__gte": 6,
                                    "standard_type__in": TYPES, "target_organism": "Homo sapiens",
                                    "data_validity_comment__isnull": "true", "limit": 1000,
                                    "only": "molecule_chembl_id,target_chembl_id,pchembl_value"})
        return pages(f"{API}/activity.json?{q}", "activities")

    mols = sorted(mols)
    chunks = [mols[i:i + 40] for i in range(0, len(mols), 40)]
    rows = []
    with ThreadPoolExecutor(4) as ex:
        for k, r in enumerate(ex.map(batch, chunks)):
            rows += r
            if k % 25 == 0:
                print(f"  ChEMBL activities: {k + 1}/{len(chunks)} batches, {len(rows):,} rows", flush=True)
    d = pd.DataFrame(rows, columns=["molecule_chembl_id", "target_chembl_id", "pchembl_value"])
    d["pchembl_value"] = d.pchembl_value.astype(float)
    CHEMBL.mkdir(parents=True, exist_ok=True)
    d.to_parquet(f, index=False)
    return d


def chembl_targets(tids):
    f = CHEMBL / "targets.json"
    if f.exists():
        return json.loads(f.read_text("utf-8"))
    tids = sorted(tids)
    out = {}
    for i in range(0, len(tids), 50):
        q = urllib.parse.urlencode({"target_chembl_id__in": ",".join(tids[i:i + 50]), "limit": 1000,
                                    "only": "target_chembl_id,target_type,organism,target_components"})
        for t in pages(f"{API}/target.json?{q}", "targets"):
            if t["target_type"] == "SINGLE PROTEIN" and t["organism"] == "Homo sapiens":
                acc = [c["accession"] for c in t["target_components"] if c.get("accession")]
                if len(acc) == 1:
                    out[t["target_chembl_id"]] = acc[0]
    f.write_text(json.dumps(out), "utf-8")
    return out


def gtopdb():
    GTOP.mkdir(parents=True, exist_ok=True)
    tabs = {}
    for name in ("interactions.csv", "ligand_id_mapping.csv"):
        f = GTOP / name
        if not f.exists():
            urllib.request.urlretrieve(f"https://www.guidetopharmacology.org/DATA/{name}", f)
        tabs[name] = pd.read_csv(f, skiprows=1, dtype=str, low_memory=False)
    version = (GTOP / "interactions.csv").open(encoding="utf-8").readline().strip().strip('"').lstrip('# ')
    return tabs["interactions.csv"], tabs["ligand_id_mapping.csv"], version


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources import ids
    g = Graph()
    n = g.nodes
    key_node = {(f"DB:{i}" if s == "DrugBank" else f"CHEMBL:{i}"): int(k)
                for k, i, t, s in zip(n.node_index, n.node_id, n.node_type, n.node_source) if t == "drug"}
    gene_node = {str(i): int(k) for k, i, t in zip(n.node_index, n.node_id, n.node_type) if t == "gene/protein"}
    h = pd.read_csv(ROOT / "data" / "raw" / "hgnc" / "hgnc_complete_set.txt", sep="\t", dtype=str, low_memory=False,
                    usecols=["entrez_id", "uniprot_ids"]).dropna()
    uni_gene = {}
    for e, us in zip(h.entrez_id, h.uniprot_ids):
        for u in us.split("|"):
            if e in gene_node:
                uni_gene.setdefault(u, gene_node[e])
    mol = {cid: key_node[k[0]] for cid, k in ids.drugs().items() if k[0] in key_node}
    act = chembl_activities(set(mol))
    tgt = chembl_targets(set(act.target_chembl_id))
    act["drug"] = act.molecule_chembl_id.map(mol)
    act["gene"] = act.target_chembl_id.map(tgt).map(uni_gene)
    a = act.dropna(subset=["drug", "gene"]).groupby(["drug", "gene"]).pchembl_value.max().reset_index()
    a["source"] = "ChEMBL"
    inter, lig, version = gtopdb()
    lig_chembl = dict(zip(lig["Ligand id"], lig["ChEMBl ID"]))
    i = inter[inter["Target Species"] == "Human"].copy()
    med = pd.to_numeric(i["Affinity Median"], errors="coerce")
    mid = (pd.to_numeric(i["Affinity High"], errors="coerce") + pd.to_numeric(i["Affinity Low"], errors="coerce")) / 2
    i["px"] = med.fillna(mid)
    i["drug"] = i["Ligand ID"].map(lig_chembl).map(mol)
    i["gene"] = i["Target UniProt ID"].map(uni_gene)
    i = i[i.px >= 6].dropna(subset=["drug", "gene"]).groupby(["drug", "gene"]).px.max().reset_index()
    i = i.rename(columns={"px": "pchembl_value"})
    i["source"] = "GtoPdb"
    out = pd.concat([a, i], ignore_index=True).astype({"drug": int, "gene": int})
    out.to_parquet(DEST, index=False)
    status = get(f"{API}/status.json")["chembl_db_version"]
    meta = {"chembl": status, "gtopdb": version, "chembl_molecules_queried": len(mol),
            "chembl_activity_rows": len(act), "pairs": len(out), "drugs": int(out.drug.nunique()),
            "genes": int(out.gene.nunique()), "pairs_by_source": out.source.value_counts().to_dict()}
    (ROOT / "data" / "ext" / "measured_targets.json").write_text(json.dumps(meta, indent=1), "utf-8")
    print(json.dumps(meta, indent=1))
    return out


if __name__ == "__main__":
    build()
