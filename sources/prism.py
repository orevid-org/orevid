"""PRISM Repurposing 19Q4 secondary screen (Corsello et al., Nat Cancer 2020; figshare 10.6084/m9.figshare.9393293,
CC BY 4.0): drug sensitivity of cancer cell lines (O2, docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, evening»).

Per drug and cell line the median AUC over screens (lower = more sensitive), then per drug and lineage (primary
tissue) the median over cell lines. Drugs -> graph by name through active moieties and synonyms (sources/practical.py).

    py -3 -m sources.prism       # -> data/ext/prism_lineage.parquet (drug node, lineage, median AUC, cell lines)
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "prism"
DEST = ROOT / "data" / "ext" / "prism_lineage.parquet"


def drug_index(g):
    """moiety -> graph drug nodes (names and synonyms)."""
    from sources.practical import moiety, synonyms
    syn = synonyms()
    by_key = defaultdict(set)
    for f in g.drugs:
        f = int(f)
        for k in ({moiety(g.name(f))} | {moiety(s) for s in syn.get(f, ())}) - {""}:
            by_key[k].add(f)
    return by_key


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.practical import moiety
    g = Graph()
    d = pd.read_csv(RAW / "secondary-screen-dose-response-curve-parameters.csv",
                    usecols=["broad_id", "depmap_id", "auc", "name"], low_memory=False).dropna(subset=["auc", "name"])
    cl = pd.read_csv(RAW / "secondary-screen-cell-line-info.csv", usecols=["depmap_id", "primary_tissue"])
    d = d.merge(cl.drop_duplicates("depmap_id"), on="depmap_id")
    d = d[d.primary_tissue != "fibroblast"]
    per_cell = d.groupby(["name", "depmap_id", "primary_tissue"]).auc.median().reset_index()
    lin = per_cell.groupby(["name", "primary_tissue"]).agg(auc=("auc", "median"), cell_lines=("auc", "size")).reset_index()
    by_key = drug_index(g)
    rows = []
    names = lin.name.unique()
    mapped = {n: by_key.get(moiety(n), set()) for n in names}
    for n, t, a, k in lin.itertuples(index=False):
        for f in mapped[n]:
            rows.append((f, n, t, a, k))
    out = pd.DataFrame(rows, columns=["drug", "prism_name", "lineage", "auc", "cell_lines"])
    out.to_parquet(DEST, index=False)
    print(f"PRISM: {len(names):,} compounds, {sum(bool(v) for v in mapped.values()):,} mapped to graph drugs "
          f"({out.drug.nunique():,} drugs), {out.lineage.nunique()} lineages, {len(out):,} drug–lineage rows -> {DEST}")
    return out


if __name__ == "__main__":
    build()
