"""PharmGKB (now ClinPGx) clinical annotations, CC BY-SA 4.0 (data/raw/pharmgkb/clinicalAnnotations.zip), for the
internal graph (X1, docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, evening»). Drugs -> graph by name (moieties,
synonyms), genes by symbol; one row per (drug, gene) with the best level of evidence and the phenotype categories.

    py -3 -m sources.pharmgkb    # -> data/ext/pharmgkb.parquet
"""
import sys
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "pharmgkb" / "clinicalAnnotations.zip"
DEST = ROOT / "data" / "ext" / "pharmgkb.parquet"
ORDER = ["1A", "1B", "2A", "2B", "3", "4"]


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.practical import moiety
    from sources.prism import drug_index
    g = Graph()
    by_name = drug_index(g)
    gene_node = {str(s).upper(): int(k) for k, s, t in zip(g.nodes.node_index, g.nodes.node_name, g.nodes.node_type)
                 if t == "gene/protein"}
    c = pd.read_csv(zipfile.ZipFile(RAW).open("clinical_annotations.tsv"), sep="\t")
    rows = []
    for r in c.itertuples():
        genes = [gene_node.get(x.strip().upper()) for x in str(r.Gene).split(";")]
        for name in str(r._11).split(";"):
            for f in by_name.get(moiety(name.strip()), set()):
                for gn in genes:
                    if gn is not None:
                        rows.append({"drug": f, "gene": gn, "level": r._4, "category": r._8, "phenotype": r._12, "url": r.URL})
    out = pd.DataFrame(rows)
    out["rank"] = out.level.map({k: i for i, k in enumerate(ORDER)})
    out = (out.sort_values("rank").groupby(["drug", "gene"])
           .agg(level=("level", "first"), category=("category", lambda x: "; ".join(sorted(set(map(str, x))))),
                n_annotations=("level", "size"), url=("url", "first")).reset_index())
    out.to_parquet(DEST, index=False)
    print(f"PharmGKB: {len(c):,} clinical annotations; {len(out):,} drug–gene pairs ({out.drug.nunique():,} drugs); "
          f"levels {out.level.value_counts().to_dict()} -> {DEST}")
    return out


if __name__ == "__main__":
    build()
