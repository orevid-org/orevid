"""QTL colocalisations from Open Targets 26.06 (M1, docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026: M1»).

Streams the 200 parts of `colocalisation`, keeps the rows whose right side is a pQTL or eQTL with h4 ≥ 0.8 and
betaRatioSignAverage = ±1, adds the QTL gene, cis/trans and the GWAS disease ids, and writes one table of unique
(gene, type, cis/trans, diseases, sign) rows with counts. Each part is deleted after reading.

    py -3 -m sources.coloc       # -> data/ext/coloc_qtl.parquet
"""
import re
import sys
import urllib.request
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.06/output/colocalisation/"
TMP = ROOT / "data" / "raw" / "opentargets" / "colocalisation_tmp"
DEST = ROOT / "data" / "ext" / "coloc_qtl.parquet"
OT = ROOT / "data" / "raw" / "opentargets"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    with urllib.request.urlopen(BASE, timeout=120) as r:
        parts = sorted(set(re.findall(r'href="(part-[^"]+\.parquet)"', r.read().decode())))
    TMP.mkdir(parents=True, exist_ok=True)
    cs = pd.read_parquet(OT / "credible_set_light.parquet", columns=["studyLocusId", "studyId", "isTransQtl"])
    st = pd.read_parquet(OT / "study.parquet", columns=["studyId", "geneId", "diseaseIds"])
    loc_study = dict(zip(cs.studyLocusId, cs.studyId))
    loc_trans = dict(zip(cs.studyLocusId, cs.isTransQtl))
    gene = dict(zip(st.studyId, st.geneId))
    dis = {s: "|".join(v) for s, v in zip(st.studyId, st.diseaseIds) if v is not None and len(v)}
    del cs, st
    cols = ["gene", "rightStudyType", "trans", "disease_ids", "betaRatioSignAverage"]

    def one(name):
        """Downloads one part (to .tmp, then renamed), keeps the QTL rows as unique (gene, type, cis/trans, GWAS
        diseases, sign) with the number of colocalisations and the best h4, deletes the part."""
        f = TMP / name
        if not f.exists():
            urllib.request.urlretrieve(BASE + name, f.with_suffix(".tmp"))
            f.with_suffix(".tmp").rename(f)
        d = pd.read_parquet(f, columns=["leftStudyLocusId", "rightStudyLocusId", "rightStudyType", "h4",
                                        "betaRatioSignAverage"])
        n = len(d)
        d = d[d.rightStudyType.isin(["pqtl", "eqtl"]) & (d.h4 >= 0.8) & d.betaRatioSignAverage.abs().eq(1)]
        d = d.assign(gene=d.rightStudyLocusId.map(loc_study).map(gene), trans=d.rightStudyLocusId.map(loc_trans),
                     disease_ids=d.leftStudyLocusId.map(loc_study).map(dis).fillna(""))
        d = d.groupby(cols, dropna=False).agg(n=("h4", "size"), h4=("h4", "max")).reset_index()
        f.unlink()
        return n, d

    keep, n_rows = [], 0
    with ThreadPoolExecutor(8) as ex:  # one connection to the EBI FTP gives ~1 MB/s
        for k, (n, d) in enumerate(ex.map(one, parts)):
            n_rows += n
            keep.append(d)
            print(f"  {k + 1}/{len(parts)}: kept {len(d):,}", flush=True)
    c = pd.concat(keep, ignore_index=True)
    c = c.groupby(cols, dropna=False).agg(n=("n", "sum"), h4=("h4", "max")).reset_index()
    c.to_parquet(DEST, index=False)
    print(f"{n_rows:,} colocalisations read; {c.n.sum():,} QTL colocalisations kept in {len(c):,} unique rows "
          f"({c.n[c.rightStudyType == 'pqtl'].sum():,} pQTL) -> {DEST}")
    TMP.rmdir()
    return c


if __name__ == "__main__":
    build()
