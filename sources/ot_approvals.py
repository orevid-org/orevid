"""Which Open Targets approvals are confirmed by a single-drug report (DECISIONS.md P17, P31).

Open Targets 26.06 turns into "approved" the substances mentioned in an FDA label: lincomycin
ends up "approved" for multiple sclerosis because it appears in the fingolimod label (a
report with 3 molecules). Here each approved (drug, disease) pair is classified:

    confirmed    some approval label or regulatory agency report with a single drug and with the
                 disease name written out (not an acronym of ≤ 3 letters)
    unconfirmed  only reports with several drugs or acronyms, or no approval report at all

Second error found: the PMDA writes "ad" (ADHD) and Open Targets reads it as Alzheimer's disease
(atomoxetine and methylphenidate "approved" for Alzheimer's). Third: curated sources without a label. Agreement
with PrimeKG of the pairs confirmed only by each type: label or agency 25.7 % (5,805), only ATC
class 4.2 % (2,482), only TTD 14.4 % (1,599). It is not an error measure (PrimeKG is incomplete), but the
ATC classes give generic approvals and TTD has errors: they do not confirm.

Legitimate combinations (two drugs in one product) also fall into "unconfirmed": few are lost and nothing
wrong is shown. The website shows the unconfirmed ones separately and no Open Targets approval counts as
evidence (the validated definition uses the PrimeKG indications).

    py -3 -m sources.ot_approvals     # -> data/ext/ot_approvals.parquet (indices of the prepared graph)
"""
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

from sources import ids
from sources.ids import primekg_keys

OT = ids.OT
ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "data" / "ext" / "ot_approvals.parquet"
APPROVAL = {"APPROVAL", "PREAPPROVAL"}


def node_index():
    """Canonical key -> node index of the prepared graph (like prepare.merge_external)."""
    n = pd.read_pickle(ROOT / "data" / "prep" / "nodes.pkl")
    index = {}
    for i, t, nid, origin in zip(n.node_index, n.node_type, n.node_id, n.origin):
        keys = primekg_keys(t, nid) if origin == "PrimeKG" else [nid]
        for c in keys:
            index.setdefault(c, int(i))
    return index


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    ci = pd.read_parquet(OT / "clinical_indication", columns=["drugId", "diseaseId", "maxClinicalStage"])
    ci = ci[ci.maxClinicalStage.isin(APPROVAL)]
    r = pd.read_parquet(OT / "clinical_report", columns=["clinicalStage", "type", "drugs", "diseases"])
    # Only labels and regulatory agencies: the curated sources (TTD, ATC classes) agree much less
    # with DrugBank (14 % and 4 %, against 26 %) and have errors like leflunomide for MS (TTD).
    r = r[r.clinicalStage.isin(APPROVAL) & r.type.isin(["DRUG_LABEL", "REGULATORY_AGENCY"])]
    single = defaultdict(bool)
    for drugs, diseases in zip(r.drugs, r.diseases):
        if drugs is None or diseases is None:
            continue
        one = len({d["drugId"] for d in drugs}) == 1
        for a in drugs:
            for b in diseases:
                # "ad" (ADHD at the PMDA) is read as Alzheimer's: an acronym of ≤ 3 letters does not confirm
                name = str(b.get("diseaseFromSource") or "").strip()
                single[(a["drugId"], b["diseaseId"])] |= one and len(name) > 3
    ci["status"] = ["confirmed" if single.get((a, b)) else "unconfirmed" for a, b in zip(ci.drugId, ci.diseaseId)]
    drug_map, (dis_map, _, _, _) = ids.drugs(), ids.diseases()
    index = node_index()
    ci["x_index"] = [index.get((drug_map.get(a) or (None,))[0]) for a in ci.drugId]
    ci["y_index"] = [index.get(dis_map.get(b)) for b in ci.diseaseId]
    ci = ci.dropna(subset=["x_index", "y_index"]).astype({"x_index": int, "y_index": int})
    # Several Open Targets ids can land on the same node pair: one confirmed is enough.
    res = (ci.assign(conf=ci.status == "confirmed").groupby(["x_index", "y_index"], as_index=False).conf.any()
             .assign(status=lambda d: d.conf.map({True: "confirmed", False: "unconfirmed"})).drop(columns="conf"))
    res.to_parquet(DEST, index=False)
    print(res.status.value_counts().to_string())
    print(f"-> {DEST}")


if __name__ == "__main__":
    main()
