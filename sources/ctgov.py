"""Phase 2–3 trials with posted results from the ClinicalTrials.gov API v2 (public domain), for T1
(docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026, evening: T1»). AACT snapshots need an account; the API does not.

Kept: one drug / biological intervention that is not a placebo; conditions and interventions by their MeSH terms
(MeSH -> MONDO by SSSOM exact match -> graph disease; intervention names -> graph drugs by moieties and synonyms).
Outcome: «met» if any primary-outcome analysis has p < 0.05, «not met» if primary analyses exist and all are ≥ 0.05.

    py -3 -m sources.ctgov download   # -> data/raw/ctgov/studies.jsonl
    py -3 -m sources.ctgov            # -> data/ext/ctgov_trials.parquet
"""
import json
import re
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "ctgov" / "studies.jsonl"
DEST = ROOT / "data" / "ext" / "ctgov_trials.parquet"
API = "https://clinicaltrials.gov/api/v2/studies"
FIELDS = ("NCTId,Phase,OverallStatus,InterventionType,InterventionName,ConditionMeshTerm,ConditionMeshId,"
          "InterventionMeshTerm,OutcomeMeasureType,OutcomeAnalysisPValue,CompletionDate")


def download():
    RAW.parent.mkdir(parents=True, exist_ok=True)
    params = {"aggFilters": "results:with,studyType:int", "filter.advanced": "AREA[Phase](PHASE2 OR PHASE3)",
              "fields": FIELDS, "pageSize": 1000, "countTotal": "true"}
    n, token = 0, None
    with RAW.open("w", encoding="utf-8") as fh:
        while True:
            q = dict(params, **({"pageToken": token} if token else {}))
            for k in range(5):
                try:
                    with urllib.request.urlopen(f"{API}?{urllib.parse.urlencode(q)}", timeout=180) as r:
                        j = json.load(r)
                    break
                except Exception:
                    if k == 4:
                        raise
                    time.sleep(10)
            for s in j["studies"]:
                fh.write(json.dumps(s) + "\n")
            n += len(j["studies"])
            print(f"  {n:,} / {j.get('totalCount', '?')}", flush=True)
            token = j.get("nextPageToken")
            if not token:
                break
    (RAW.parent / "VERSION").write_text(f"ClinicalTrials.gov API v2, downloaded {time.strftime('%Y-%m-%d')}\n")


def pvalue(s):
    m = re.search(r"([0-9]*\.?[0-9]+(?:[eE]-?\d+)?)", str(s))
    if not m:
        return None
    v = float(m.group(1))
    return v - 1e-12 if str(s).strip().startswith("<") else v


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.ids import xref_to_mondo
    from sources.practical import moiety
    from sources.prism import drug_index
    g = Graph()
    to_mondo = xref_to_mondo()
    mondo_node = {}
    for k, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type):
        if t == "disease" and str(i).replace("_", "").isdigit():
            for p in str(i).split("_"):
                mondo_node.setdefault(f"MONDO:{int(p)}", int(k))
    by_name = drug_index(g)
    rows, n = [], 0
    for line in RAW.open(encoding="utf-8"):
        s = json.loads(line)
        n += 1
        ps, rs, ds = s.get("protocolSection", {}), s.get("resultsSection", {}), s.get("derivedSection", {})
        iv = [x for x in ps.get("armsInterventionsModule", {}).get("interventions", [])
              if x.get("type") in ("DRUG", "BIOLOGICAL") and "placebo" not in x.get("name", "").lower()]
        if len(iv) != 1:
            continue
        fs = by_name.get(moiety(iv[0]["name"]), set())
        for m in ds.get("interventionBrowseModule", {}).get("meshes", []):
            fs = fs or by_name.get(moiety(m.get("term", "")), set())
        meshes = [m["id"] for m in ds.get("conditionBrowseModule", {}).get("meshes", [])]
        dis = {mondo_node[x] for m in meshes for x in to_mondo.get(f"MESH:{m}", ()) if x in mondo_node}
        pv = [pvalue(a.get("pValue")) for o in rs.get("outcomeMeasuresModule", {}).get("outcomeMeasures", [])
              if o.get("type") == "PRIMARY" for a in o.get("analyses", []) if a.get("pValue")]
        pv = [x for x in pv if x is not None]
        outcome = "met" if any(x < 0.05 for x in pv) else "not met" if pv else "unknown"
        phases = "/".join(ps.get("designModule", {}).get("phases", []))
        nct = ps.get("identificationModule", {}).get("nctId")
        for f in fs:
            for d in dis:
                rows.append({"drug": f, "disease": d, "nct": nct, "phase": phases, "outcome": outcome,
                             "status": ps.get("statusModule", {}).get("overallStatus"), "min_p": min(pv) if pv else None,
                             "completion": ps.get("statusModule", {}).get("completionDateStruct", {}).get("date")})
    out = pd.DataFrame(rows).drop_duplicates(["drug", "disease", "nct"])
    out.to_parquet(DEST, index=False)
    print(f"ClinicalTrials.gov: {n:,} phase 2–3 studies with results; {out.nct.nunique():,} single-drug trials mapped to "
          f"{len(out):,} drug–disease rows; outcomes {out.drop_duplicates('nct').outcome.value_counts().to_dict()} -> {DEST}")
    return out


if __name__ == "__main__":
    download() if sys.argv[1:] == ["download"] else build()
