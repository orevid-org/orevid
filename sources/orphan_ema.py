"""EMA orphan designations (medicines-output report, www.ema.europa.eu, reuse with acknowledgement) for L2
(docs/IMPROVEMENT_PLAN.md, «Registered 24-10-2026: L1–L6»). Active substance -> graph drug by name (moieties, synonyms);
the intended use («Treatment of …», «Prevention of …», «Diagnosis of …») -> graph disease by exact normalised name or
synonym (sources/drug_status.disease_names). Refused designations are left out.

    py -3 -m sources.orphan_ema   # -> data/ext/orphan_ema.parquet
"""
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "orphan_ema" / "orphan_designations.xlsx"
DEST = ROOT / "data" / "ext" / "orphan_ema.parquet"
LEAD = re.compile(r"^(treatment|prevention|diagnosis|treatment and prevention|prevention and treatment)\s+of\s+(the\s+)?", re.I)


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.drug_status import disease_names
    from sources.practical import moiety
    from sources.prism import drug_index
    g = Graph()
    head = pd.read_excel(RAW, header=None)
    h = next(i for i in range(30) if head.iloc[i].notna().sum() > 5)
    d = pd.read_excel(RAW, header=h)
    d = d[~d.Status.astype(str).str.lower().str.contains("refus")]
    by_name = drug_index(g)
    name_dis = defaultdict(set)
    for dis, names in disease_names(g).items():
        for n in names:
            name_dis[n].add(dis)
    rows, n_drug, n_dis = [], 0, 0
    for sub, use, date, status, num, url in zip(d["Active substance"], d["Intended use"], d["Date of designation / refusal"],
                                                d.Status, d["EU designation number"], d["Orphan designation URL"]):
        fs = set().union(*(by_name.get(moiety(s.strip()), set()) for s in re.split(r"[,;/]| and ", str(sub)))) if sub == sub else set()
        text = LEAD.sub("", str(use)).strip().rstrip(".").lower()
        ds = name_dis.get(text, set())
        n_drug += bool(fs)
        n_dis += bool(ds)
        for f in fs:
            for dd in ds:
                rows.append({"drug": f, "disease": dd, "date": pd.to_datetime(date, dayfirst=True, errors="coerce"),
                             "status": status, "designation": num, "intended_use": use, "url": url})
    out = pd.DataFrame(rows).drop_duplicates(["drug", "disease", "designation"])
    out.to_parquet(DEST, index=False)
    print(f"EMA orphan designations: {len(d):,} (not refused); drug mapped {n_drug:,}, disease mapped {n_dis:,}; "
          f"{len(out):,} graph pairs ({out.drug.nunique():,} drugs, {out.disease.nunique():,} diseases), "
          f"{(out.date >= '2023-01-01').sum():,} dated 2023 or later -> {DEST}")


if __name__ == "__main__":
    build()
