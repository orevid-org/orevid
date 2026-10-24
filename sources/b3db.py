"""B3DB blood–brain-barrier permeability (theochem/B3DB, CC0) for N7 (docs/IMPROVEMENT_PLAN.md, «Registered 24-10-2026:
N5–N7»). Compounds -> graph drugs by name (moieties, synonyms); a drug with both labels is left out.

    py -3 -m sources.b3db          # -> data/ext/b3db.parquet (drug, bbb = +1 / −1)
"""
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.practical import moiety
    from sources.prism import drug_index
    g = Graph()
    b = pd.read_csv(ROOT / "data" / "raw" / "b3db" / "B3DB_classification.tsv", sep="\t", low_memory=False)
    by = drug_index(g)
    lab = {}
    for nm, v in zip(b.compound_name, b["BBB+/BBB-"]):
        for f in by.get(moiety(str(nm)), set()) if nm == nm else ():
            lab.setdefault(f, set()).add(1 if v == "BBB+" else -1)
    out = pd.DataFrame([{"drug": f, "bbb": next(iter(s))} for f, s in lab.items() if len(s) == 1])
    out.to_parquet(ROOT / "data" / "ext" / "b3db.parquet", index=False)
    print(f"B3DB: {len(out):,} graph drugs ({(out.bbb > 0).sum():,} BBB+, {(out.bbb < 0).sum():,} BBB−)")


if __name__ == "__main__":
    build()
