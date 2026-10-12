"""CTD (Comparative Toxicogenomics Database) -> data/ext/.

    py -3 -m sources.ctd

Non-commercial use. CTD's terms ask to cite them and to notify them if anything
made with their data is published (http://ctdbase.org/about/legal.jsp).

New relations (only for the chemicals that are drugs in the graph):
    ctd_gene     drug -> gene (display: increases / decreases / other)
                 Interactions curated from the literature: the drug increases or
                 decreases the expression or activity of the gene. They come from
                 any organism; the number of articles is kept.
    ctd_disease  drug -> disease (display: therapeutic / marker/mechanism)
                 Only the direct (curated) evidence. "therapeutic" is literature
                 saying that the drug treats the disease: it must NOT be used
                 as a signal in validation, because it would be the answer.

Mapping chemicals to drugs: by InChIKey (with the Open Targets molecules)
and, failing that, by name or an exact synonym.
"""
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from sources import ids

RAW = ids.RAW / "ctd"
EXT = Path(__file__).parent.parent / "data" / "ext"


def read(name, columns):
    return pd.read_csv(RAW / name, sep="\t", comment="#", header=None, names=columns,
                       dtype=str, low_memory=False)


def chemical_map():
    """CTD chemical id (MeSH, without prefix) -> graph drug key."""
    q = read("CTD_chemicals.tsv.gz",
             ["name", "id", "cas", "pubchem_cid", "pubchem_sid", "dtxsid", "inchikey", "definition",
              "parents", "trees", "parent_trees", "mesh_synonyms", "ctd_synonyms"])
    q["id"] = q.id.str.removeprefix("MESH:")

    m = pd.read_parquet(ids.OT / "drug_molecule", columns=["id", "inchiKey", "name", "synonyms", "tradeNames"])
    drugs = ids.drugs()
    by_inchi = {k: drugs[c][0] for c, k in zip(m.id, m.inchiKey) if isinstance(k, str) and c in drugs}
    by_name = {}
    for c, name, syn, trade in zip(m.id, m.name, m.synonyms, m.tradeNames):
        if c not in drugs:
            continue
        for x in [name, *(s["label"] for s in (syn if syn is not None else [])),
                  *(s["label"] for s in (trade if trade is not None else []))]:
            if isinstance(x, str):
                by_name.setdefault(x.lower(), drugs[c][0])
    # PrimeKG (DrugBank) names take precedence over Open Targets synonyms.
    for key, (node_type, name) in ids.all_primekg_keys().items():
        if node_type == "drug":
            by_name[str(name).lower()] = key

    mapping, via = {}, Counter()
    for qid, name, inchi, s1, s2 in zip(q.id, q["name"], q.inchikey, q.mesh_synonyms, q.ctd_synonyms):
        if isinstance(inchi, str) and inchi in by_inchi:
            mapping[qid] = by_inchi[inchi]; via["InChIKey"] += 1
        elif name.lower() in by_name:
            mapping[qid] = by_name[name.lower()]; via["name"] += 1
        else:
            for s in "|".join(x for x in (s1, s2) if isinstance(x, str)).split("|"):
                if s and s.lower() in by_name:
                    mapping[qid] = by_name[s.lower()]; via["synonym"] += 1
                    break
    print(f"CTD chemicals: {len(q):,}; as a graph drug: {len(mapping):,} {dict(via)}")
    return mapping


def direction(actions):
    """"decreases^expression|affects^binding" -> increases / decreases / other.

    Only the expression and activity of the gene count; if both
    directions are present, it is "other"."""
    signs = set()
    for a in actions.split("|"):
        sense, _, what = a.partition("^")
        if what in ("expression", "activity") and sense in ("increases", "decreases"):
            signs.add(sense)
    if signs == {"increases"}:
        return "increases"
    if signs == {"decreases"}:
        return "decreases"
    return "other"


def genes(mapping, nodes):
    g = read("CTD_chem_gene_ixns.tsv.gz",
             ["chemical", "chemical_id", "cas", "symbol", "gene_id", "forms", "organism", "org_id",
              "interaction", "actions", "pubmed"])
    g["x"] = g.chemical_id.map(mapping)
    g = g.dropna(subset=["x", "gene_id"])
    # Genes from other species have their own id: they are mapped to the human gene by
    # symbol (mouse Tnf -> human TNF). If there is none, out.
    h = pd.read_csv(ids.RAW / "hgnc" / "hgnc_complete_set.txt", sep="\t", low_memory=False,
                    usecols=["symbol", "entrez_id"]).dropna()
    human_keys = {f"NCBI:{int(e)}" for e in h.entrez_id}
    by_symbol = {s.upper(): f"NCBI:{int(e)}" for s, e in zip(h.symbol, h.entrez_id)}
    key = "NCBI:" + g.gene_id
    g["y"] = np.where(g.org_id == "9606", key, g.symbol.str.upper().map(by_symbol))
    g["y"] = g.y.where(g.y.isin(human_keys))
    print(f"Chemical-gene interactions with a drug: {len(g):,}; with a human gene: {g.y.notna().sum():,}")
    g = g.dropna(subset=["y"])
    g["display_relation"] = g.actions.fillna("").map(direction)
    g["human"] = g.org_id == "9606"
    g["n_articles"] = g.pubmed.fillna("").str.count(r"\|") + (g.pubmed.notna()).astype(int)
    agg = g.groupby(["x", "y", "display_relation"]).agg(
        n_articles=("n_articles", "sum"), n_interactions=("interaction", "size"), human=("human", "any")).reset_index()
    # If there are articles in both directions, the majority wins (tie -> other).
    pivot = agg.pivot_table(index=["x", "y"], columns="display_relation", values="n_articles", fill_value=0)
    for c in ("increases", "decreases"):
        if c not in pivot:
            pivot[c] = 0
    final_dir = np.where(pivot.increases > pivot.decreases, "increases",
                         np.where(pivot.decreases > pivot.increases, "decreases", "other"))
    final = pd.DataFrame({"x": pivot.index.get_level_values(0), "y": pivot.index.get_level_values(1),
                          "display_relation": final_dir})
    totals = agg.groupby(["x", "y"]).agg(n_articles=("n_articles", "sum"), human=("human", "any")).reset_index()
    final = final.merge(totals, on=["x", "y"])
    for key, symbol in g.drop_duplicates("y")[["y", "symbol"]].itertuples(index=False):
        nodes[key] = ("gene/protein", symbol, "CTD")
    return final.assign(relation="ctd_gene", score=final.n_articles.astype(float))


def diseases(mapping):
    d = read("CTD_chemicals_diseases.tsv.gz",
             ["chemical", "chemical_id", "cas", "disease", "disease_id", "evidence", "inference_gene",
              "inference_score", "omim", "pubmed"])
    d = d.dropna(subset=["evidence"])  # only the direct evidence; the inferred one comes from the genes
    d["x"] = d.chemical_id.map(mapping)
    to_mondo = ids.xref_to_mondo()
    d["y"] = [next(iter(m)) if len(m := to_mondo.get(ids._xref(i), ())) == 1 else None for i in d.disease_id]
    print(f"Direct chemical-disease: {len(d):,}; with a drug: {d.x.notna().sum():,}; "
          f"with a drug and MONDO: {(d.x.notna() & d.y.notna()).sum():,}")
    d = d.dropna(subset=["x", "y"])
    d["display_relation"] = d.evidence.map(
        lambda e: "therapeutic" if "therapeutic" in e else "marker/mechanism")
    d["n_articles"] = d.pubmed.fillna("").str.count(r"\|") + d.pubmed.notna().astype(int)
    d = d.groupby(["x", "y", "display_relation"], as_index=False).n_articles.sum()
    return d.assign(relation="ctd_disease")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    mapping = chemical_map()
    new = {}
    a = pd.concat([genes(mapping, new), diseases(mapping)], ignore_index=True).assign(source="CTD")
    primekg = ids.all_primekg_keys()
    nodes = pd.DataFrame([(k, *v) for k, v in new.items() if k not in primekg],
                         columns=["key", "node_type", "node_name", "node_source"])
    EXT.mkdir(parents=True, exist_ok=True)
    a.to_parquet(EXT / "ctd_edges.parquet", index=False)
    nodes.to_parquet(EXT / "ctd_nodes.parquet", index=False)
    print(f"\nEdges: {len(a):,}; new genes: {len(nodes):,}")
    print(a.groupby(["relation", "display_relation"]).size().to_string())


if __name__ == "__main__":
    main()
