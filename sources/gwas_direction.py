"""Direction from GWAS + QTL: activate or inhibit each gene in each disease -> data/ext/.

    py -3 -m sources.gwas_direction            # shared lead variant (cheap, ~110 MB)
    py -3 -m sources.gwas_direction --locus    # lead GWAS variant inside the QTL credible set (~1.5 GB)

Idea: if the allele that increases the risk of a disease (GWAS, beta > 0) also increases
the expression of a gene (eQTL/pQTL, beta > 0), more gene = more disease -> it must be inhibited.
If the risk allele lowers expression, it must be activated. Open Targets harmonises the GWAS and
QTL betas to the same allele (the alternative one of variantId), so the sign of
beta_GWAS × beta_QTL is enough when both come from the same variant.

Only the needed columns of the 200 credible_set files are read via HTTP
Range (sources/remote.py), without downloading them whole.

New relation: dir_gwas (disease -> gene, display activate / inhibit; n = number of GWAS-QTL
pairs that vote for it).
"""
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc

from sources import ids, remote
from sources.opentargets import Nodes, disease, gene

OT_URL = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output/credible_set/"
# Local caches of the remote reads; they keep their original file names in data/raw.
LIGHT = ids.OT / "credible_set_light.parquet"
LOCUS = ids.OT / "credible_set_locus_qtl.parquet"
EXT = Path(__file__).parent.parent / "data" / "ext"
QTL = ["eqtl", "pqtl", "sceqtl"]
COLUMNS = ["studyLocusId", "studyId", "studyType", "variantId", "beta", "pValueMantissa", "pValueExponent",
           "isTransQtl"]


def in_parallel(func, urls, threads=8):
    res, read_bytes = [], 0
    with ThreadPoolExecutor(threads) as ex:
        for i, (t, n) in enumerate(ex.map(func, urls)):
            res.append(t)
            read_bytes += n
            if (i + 1) % 20 == 0:
                print(f"  {i + 1}/{len(urls)} files, {read_bytes / 1e6:.0f} MB read", flush=True)
    return pa.concat_tables(res)


def light():
    if not LIGHT.exists():
        urls = remote.files(OT_URL)
        t = in_parallel(lambda u: remote.read(u, COLUMNS), urls)
        import pyarrow.parquet as pq
        pq.write_table(t, LIGHT)
    return pd.read_parquet(LIGHT)


def locus_qtl():
    """Variants (with beta) of the 95 % credible sets of the cis QTLs."""
    if not LOCUS.exists():
        def one(u):
            def keep(t):
                t = t.filter(pc.and_(pc.is_in(t["studyType"], pa.array(QTL)),
                                     pc.invert(pc.fill_null(t["isTransQtl"], False))))
                return t.select(["studyLocusId", "studyId", "locus"])
            return remote.read(u, ["studyLocusId", "studyId", "studyType", "isTransQtl", "locus"], keep)
        t = in_parallel(one, remote.files(OT_URL))
        rows = []
        for slid, sid, loc in zip(t["studyLocusId"].to_pylist(), t["studyId"].to_pylist(), t["locus"].to_pylist()):
            for v in loc or []:
                if v["is95CredibleSet"] and v["beta"] is not None:
                    rows.append((slid, sid, v["variantId"], v["beta"], v["posteriorProbability"]))
        pd.DataFrame(rows, columns=["studyLocusId", "studyId", "variantId", "beta", "pp"]).to_parquet(LOCUS, index=False)
    return pd.read_parquet(LOCUS)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    with_locus = "--locus" in sys.argv
    st = pd.read_parquet(ids.OT / "study.parquet", columns=["studyId", "studyType", "diseaseIds", "geneId"])
    cs = light()
    print(f"Credible sets: {len(cs):,} " + str(cs.studyType.value_counts().to_dict()))

    gw = cs[cs.studyType == "gwas"].merge(st[["studyId", "diseaseIds"]], on="studyId")
    gw = gw[gw.diseaseIds.map(lambda l: l is not None and len(l) > 0)].explode("diseaseIds")
    gw = gw.rename(columns={"diseaseIds": "diseaseId", "beta": "beta_g"})[["diseaseId", "variantId", "beta_g"]]
    gw = gw.dropna(subset=["beta_g"])

    if with_locus:
        q = locus_qtl()
        q = q[q.pp >= 0.01]
    else:
        q = cs[cs.studyType.isin(QTL) & ~cs.isTransQtl.fillna(False).astype(bool)][["studyId", "variantId", "beta"]]
    q = q.merge(st[["studyId", "geneId"]], on="studyId").dropna(subset=["geneId", "beta"])
    q = q.rename(columns={"beta": "beta_q"})[["geneId", "variantId", "beta_q"]]

    p = gw.merge(q, on="variantId")
    p["sign"] = np.sign(p.beta_g * p.beta_q)  # +: more gene, more disease -> inhibit
    p = p[p.sign != 0]
    print(f"GWAS-QTL pairs sharing a variant: {len(p):,}")

    v = p.groupby(["diseaseId", "geneId"])["sign"].agg(["sum", "size"]).reset_index()
    v = v[v["sum"] != 0]
    nodes = Nodes()
    dis = {x: disease(x, nodes) for x in v.diseaseId.unique()}
    gene_keys = {x: gene(x, nodes) for x in v.geneId.unique()}
    v["x"] = v.diseaseId.map(dis)
    v["y"] = v.geneId.map(gene_keys)
    v = v.dropna(subset=["x", "y"])
    v = v[v.x.str.startswith(("MONDO:", "OT:"))]
    v = v.groupby(["x", "y"], as_index=False).agg(vote=("sum", "sum"), n=("size", "sum"))
    v = v[v.vote != 0]
    v["display_relation"] = np.where(v.vote > 0, "inhibit", "activate")
    edges = v.assign(relation="dir_gwas", score=v.vote.abs().astype(float),
                     source="Open Targets (GWAS+QTL)")[["x", "y", "relation", "display_relation", "score", "n", "source"]]
    EXT.mkdir(parents=True, exist_ok=True)
    edges.to_parquet(EXT / "gwas_direction_edges.parquet", index=False)
    nodes.table().to_parquet(EXT / "gwas_direction_nodes.parquet", index=False)
    print(f"Disease-gene edges with direction: {len(edges):,} "
          + str(edges.display_relation.value_counts().to_dict()))


if __name__ == "__main__":
    main()
