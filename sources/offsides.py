"""OFFSIDES (nSIDES, Tatonetti lab; FAERS-derived drug–adverse-event signals; licence to check) for L6
(docs/IMPROVEMENT_PLAN.md, «Registered 24-10-2026: L1–L6»). Drugs -> graph by name (moieties, synonyms); MedDRA
conditions -> graph diseases by exact normalised name or synonym, as OnSIDES (sources/onsides.py).

    py -3 -m sources.offsides     # -> data/ext/offsides.parquet (drug, disease, term, PRR) — the same columns D1 reads
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "offsides" / "OFFSIDES.csv.gz"
DEST = ROOT / "data" / "ext" / "offsides.parquet"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.drug_status import disease_names
    from sources.practical import moiety
    from sources.prism import drug_index
    g = Graph()
    d = pd.read_csv(RAW, usecols=["drug_concept_name", "condition_concept_name", "PRR", "A"], low_memory=False)
    by_name = drug_index(g)
    name_dis = defaultdict(set)
    for dis, names in disease_names(g).items():
        for n in names:
            name_dis[n].add(dis)
    drugs = {n: by_name.get(moiety(n), set()) for n in d.drug_concept_name.unique()}
    conds = {c: name_dis.get(str(c).strip().lower(), set()) for c in d.condition_concept_name.unique()}
    rows = [(f, dd, c, prr, a) for n, c, prr, a in zip(d.drug_concept_name, d.condition_concept_name, d.PRR, d.A)
            for f in drugs[n] for dd in conds[c]]
    out = pd.DataFrame(rows, columns=["drug", "disease", "term", "prr", "reports"]).sort_values("prr", ascending=False)
    out = out.drop_duplicates(["drug", "disease"])
    out.to_parquet(DEST, index=False)
    print(f"OFFSIDES: {len(d):,} signals; {len(out):,} graph drug–disease pairs ({out.drug.nunique():,} drugs, "
          f"{out.disease.nunique():,} diseases) -> {DEST}")


if __name__ == "__main__":
    build()
