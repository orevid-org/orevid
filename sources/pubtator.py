"""Relations from the whole of PubMed extracted by PubTator3 (NCBI, public domain) -> data/ext/.

    py -3 -m sources.pubtator

PubTator3 has automatically read the ~36 million PubMed abstracts and publishes
the relations between drugs, genes and diseases, each with its PMID
(https://ftp.ncbi.nlm.nih.gov/pub/lu/PubTator3/). It is automatic extraction (BioREx):
it has errors, so only the relations with at least MIN_ARTICLES
distinct articles are kept.

New relations:
    pubtator_drug_disease  drug -> disease (treats / causes)
        Does NOT enter the model or the validation: "X treats Y" is the answer the graph
        tries to guess. It is used to know what the literature says about a candidate.
    pubtator_drug_gene     drug -> gene (increases / decreases)
    pubtator_disease_gene  disease -> gene (associated / increases / decreases)

Attributes: n_articles, pmid_min, pmid_max (the PMID grows over time: the minimum indicates
when the first evidence was published).
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from sources import ids
from sources.ctd import chemical_map

RAW = ids.RAW / "pubtator" / "relation2pubtator3.gz"
EXT = Path(__file__).parent.parent / "data" / "ext"
MIN_ARTICLES = 2

# (PubTator type, entity 1, entity 2) -> (relation, subtype, needs flipping?)
RULES = {
    ("treat", "Chemical", "Disease"): ("pubtator_drug_disease", "treats", False),
    ("cause", "Chemical", "Disease"): ("pubtator_drug_disease", "causes", False),
    ("positive_correlate", "Chemical", "Gene"): ("pubtator_drug_gene", "increases", False),
    ("negative_correlate", "Chemical", "Gene"): ("pubtator_drug_gene", "decreases", False),
    ("stimulate", "Chemical", "Gene"): ("pubtator_drug_gene", "increases", False),
    ("inhibit", "Chemical", "Gene"): ("pubtator_drug_gene", "decreases", False),
    ("associate", "Disease", "Gene"): ("pubtator_disease_gene", "associated", False),
    ("stimulate", "Disease", "Gene"): ("pubtator_disease_gene", "increases", False),
    ("inhibit", "Disease", "Gene"): ("pubtator_disease_gene", "decreases", False),
    ("positive_correlate", "Disease", "Gene"): ("pubtator_disease_gene", "increases", False),
    ("negative_correlate", "Disease", "Gene"): ("pubtator_disease_gene", "decreases", False),
}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    chemicals = chemical_map()  # MeSH without prefix -> drug
    mondo = ids.xref_to_mondo()
    h = pd.read_csv(ids.RAW / "hgnc" / "hgnc_complete_set.txt", sep="\t", low_memory=False,
                    usecols=["entrez_id"]).dropna()
    human = {str(int(e)) for e in h.entrez_id}

    def entity(series):
        type_id = series.str.split("|", n=1, expand=True)
        return type_id[0], type_id[1]

    acc = []
    lines_read = 0
    for block in pd.read_csv(RAW, sep="\t", header=None, names=["pmid", "rel_type", "a", "b"], dtype=str,
                             chunksize=3_000_000, on_bad_lines="skip"):
        lines_read += len(block)
        ta, ia = entity(block.a)
        tb, ib = entity(block.b)
        block = block.assign(ta=ta, ia=ia, tb=tb, ib=ib)
        # Relations have no fixed order: each pair is put in the order of the rules.
        flip = block.ta.isin(["Gene", "Disease"]) & block.tb.isin(["Chemical", "Disease"]) & (block.ta != block.tb)
        flip &= ~((block.ta == "Disease") & (block.tb == "Gene"))
        block.loc[flip, ["ta", "ia", "tb", "ib"]] = block.loc[flip, ["tb", "ib", "ta", "ia"]].to_numpy()
        key = list(zip(block.rel_type, block.ta, block.tb))
        rule = [RULES.get(k) for k in key]
        block = block.assign(relation=[r[0] if r else None for r in rule],
                             display_relation=[r[1] if r else None for r in rule]).dropna(subset=["relation"])

        x = np.where(block.ta == "Chemical", block.ia.str.removeprefix("MESH:").map(chemicals),
                     [next(iter(m)) if len(m := mondo.get(f"MESH:{i.removeprefix('MESH:')}", ())) == 1 else None
                      for i in block.ia])
        y = np.where(block.tb == "Gene", ["NCBI:" + g if g in human else None for g in block.ib],
                     [next(iter(m)) if len(m := mondo.get(f"MESH:{i.removeprefix('MESH:')}", ())) == 1 else None
                      for i in block.ib])
        block = block.assign(x=x, y=y).dropna(subset=["x", "y"])
        block["pmid"] = pd.to_numeric(block.pmid, errors="coerce")
        acc.append(block.groupby(["x", "y", "relation", "display_relation"]).agg(
            n_articles=("pmid", "nunique"), pmid_min=("pmid", "min"), pmid_max=("pmid", "max")).reset_index())
        print(f"  {lines_read / 1e6:.0f} M lines; pairs so far: {sum(len(a) for a in acc):,}", flush=True)

    a = pd.concat(acc).groupby(["x", "y", "relation", "display_relation"]).agg(
        n_articles=("n_articles", "sum"), pmid_min=("pmid_min", "min"), pmid_max=("pmid_max", "max")).reset_index()
    a = a[a.n_articles >= MIN_ARTICLES]
    a = a.assign(score=np.log1p(a.n_articles), source="PubTator3")
    EXT.mkdir(parents=True, exist_ok=True)
    a.to_parquet(EXT / "pubtator_edges.parquet", index=False)
    print(f"\nEdges (≥ {MIN_ARTICLES} articles): {len(a):,}")
    print(a.groupby(["relation", "display_relation"]).size().to_string())


if __name__ == "__main__":
    main()
