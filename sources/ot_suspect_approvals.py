"""Open Targets «approved» indications that are probably wrong, for the Open Targets team
(opentargets/issues#4579; docs/OPEN_TARGETS_FEEDBACK.md).

Works on the Open Targets 26.06 tables only (no PrimeKG needed for the list itself). For every
(drug, disease) pair at the approval stage in `clinical_indication`, it looks at the reports
behind it and flags two patterns found by hand:

    passenger     every approval-level label or regulatory report for the pair names several
                  molecules and none names the drug alone for the disease. Example: lincomycin for multiple
                  sclerosis, from the fingolimod lauryl sulfate label (FDA 214962), which lists
                  3 molecules. Fixed-dose combinations whose component is never approved alone for
                  the disease (clavulanate, cilastatin, ritonavir as booster) also match: review.
    acronym       every approval-level report for the pair gives the disease as a text of
                  ≤ 3 characters (the PMDA's «ad» = ADHD grounded as Alzheimer's disease).

Pairs supported only by curated resources (TTD, ATC, …) are not flagged here; they are counted
in the summary. Agreement with PrimeKG (DrugBank indications) is added when the prepared graph
is available, as a hint: an indication DrugBank also has is less likely to be an error.

    py -3 -m sources.ot_suspect_approvals     # -> data/ext/ot_suspect_approvals.csv + summary
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

from sources import ids

OT = ids.OT
ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "data" / "ext" / "ot_suspect_approvals.csv"
APPROVAL = {"APPROVAL", "PREAPPROVAL"}
OFFICIAL = {"DRUG_LABEL", "REGULATORY_AGENCY"}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ci = pd.read_parquet(OT / "clinical_indication", columns=["drugId", "diseaseId", "maxClinicalStage", "clinicalReportIds"])
    ci = ci[ci.maxClinicalStage.isin(APPROVAL)]
    r = pd.read_parquet(OT / "clinical_report", columns=["id", "clinicalStage", "type", "source", "drugs", "diseases", "url"])
    r = r[r.clinicalStage.isin(APPROVAL)]
    names = dict(zip(*pd.read_parquet(OT / "drug_molecule", columns=["id", "name"]).values.T))
    dnames = dict(zip(*pd.read_parquet(OT / "disease", columns=["id", "name"]).values.T))

    # Per (drug, disease): the approval-level reports that mention both.
    reports = defaultdict(list)
    single_official = defaultdict(bool)
    for rid, typ, src, drugs, diseases, url in zip(r.id, r.type, r.source, r.drugs, r.diseases, r.url):
        if drugs is None or diseases is None:
            continue
        mols = sorted({d["drugId"] for d in drugs if d["drugId"]})
        for b in diseases:
            if not b.get("diseaseId"):
                continue
            text = str(b.get("diseaseFromSource") or "").strip()
            for a in mols:
                reports[(a, b["diseaseId"])].append(dict(id=rid, type=typ, source=src, n=len(mols), mols=mols, text=text, url=url))
                if typ in OFFICIAL and len(mols) == 1:
                    single_official[(a, b["diseaseId"])] = True

    # How many official approval reports carry exactly the same set of molecules: a fixed-dose
    # combination recurs across labels and agencies; a drug merely mentioned in one label does not.
    same_set = defaultdict(set)
    for rid, typ, drugs in zip(r.id, r.type, r.drugs):
        if typ in OFFICIAL and drugs is not None:
            same_set[tuple(sorted({d["drugId"] for d in drugs if d["drugId"]}))].add(rid)
    pkg = primekg_indications()
    rows, counts = [], defaultdict(int)
    for a, b in zip(ci.drugId, ci.diseaseId):
        rep = reports.get((a, b), [])
        off = [x for x in rep if x["type"] in OFFICIAL]
        if not off:
            counts["curated resources only"] += 1
            continue
        kind = None
        if all(x["n"] > 1 for x in off) and not single_official[(a, b)]:
            kind = "passenger"
        if kind is None and all(len(x["text"]) <= 3 for x in rep):
            kind = "acronym"
        if kind is None:
            counts["not flagged"] += 1
            continue
        counts[kind] += 1
        x = off[0]
        rows.append({
            "pattern": kind, "drugId": a, "drug": names.get(a, ""), "diseaseId": b, "disease": dnames.get(b, ""),
            "disease_text_in_report": "; ".join(sorted({y["text"] for y in off})),
            "molecules_in_report": "; ".join(f"{m} ({names.get(m, '')})" for m in x["mols"]),
            "approved_alone_for_disease": "; ".join(f"{m} ({names.get(m, '')})" for m in sorted(
                {m for y in off for m in y["mols"] if m != a and single_official[(m, b)]})),
            "report_ids": "; ".join(sorted({y["id"] for y in off})),
            "report_sources": "; ".join(sorted({f"{y['type']}/{y['source']}" for y in off})),
            "url": x["url"] or "",
            "in_primekg_drugbank": pkg(a, b),
            "reports_with_same_molecules": min(len(same_set[tuple(y["mols"])]) for y in off),
        })
    df = pd.DataFrame(rows)
    # Priority for manual review: acronyms; passengers from a molecule set seen in a single report
    # and not a DrugBank indication; then the rest (mostly genuine fixed-dose combinations).
    df["priority"] = ["high" if k == "acronym" or (n == 1 and pk != "yes") else "low"
                      for k, n, pk in zip(df.pattern, df.reports_with_same_molecules, df.in_primekg_drugbank)]
    df = df.sort_values(["priority", "pattern", "drug", "disease"])
    DEST.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(DEST, index=False)
    print(f"Approval-stage pairs in clinical_indication: {len(ci):,}")
    for k, v in counts.items():
        print(f"  {k}: {v:,}")
    for k, e in df.groupby("pattern"):
        print(f"\n{k}: {len(e)} pairs; also a DrugBank indication in PrimeKG: {(e.in_primekg_drugbank == 'yes').sum()}"
              f" ; not in PrimeKG: {(e.in_primekg_drugbank == 'no').sum()}")
        h = e[e.priority == "high"]
        print(f"  high priority: {len(h)}")
        print(h.head(40)[["drug", "disease", "disease_text_in_report", "molecules_in_report"]].to_string(index=False))
    print(f"-> {DEST}")


def primekg_indications():
    """Function (ChEMBL id, disease id) -> 'yes'/'no'/'' (unknown) using PrimeKG's indications."""
    try:
        from sources.ot_approvals import node_index
        import numpy as np
        index = node_index()
        drug_map, (dis_map, _, _, _) = ids.drugs(), ids.diseases()
        e = pd.read_csv(ROOT / "data" / "raw" / "edges.csv", usecols=["relation", "x_index", "y_index"])
    except Exception:
        return lambda a, b: ""
    e = e[e.relation == "indication"]
    pairs = set(zip(e.x_index, e.y_index)) | set(zip(e.y_index, e.x_index))

    def f(a, b):
        x = index.get((drug_map.get(a) or (None,))[0])
        y = index.get(dis_map.get(b))
        if x is None or y is None:
            return ""
        return "yes" if (x, y) in pairs else "no"
    return f


if __name__ == "__main__":
    main()
