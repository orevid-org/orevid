"""S1 (docs/IMPROVEMENT_PLAN.md, «Registered 08-11-2026: S1, R6b, R8b»): a drug is *approved somewhere* if any of
four sources says so — DrugBank group «approved» (or «withdrawn in some market»), DrugCentral (an approval by any
agency), ChEMBL maximum phase «approval» (via Open Targets) or a WHO ATC code."""
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def approved_somewhere(G):
    """{drug node: [sources saying it is approved]} for every drug approved by at least one source."""
    from predictor import drug_status
    from sources import ids
    status = drug_status(G)
    out = {}

    def add(f, s):
        out.setdefault(int(f), []).append(s)

    for f, s in status.items():
        if s in ("approved", "withdrawn in some market"):
            add(f, "DrugBank")
    for f in set(pd.read_parquet(ROOT / "data" / "ext" / "drugcentral_approvals.parquet").drug.astype(int)):
        add(f, "DrugCentral")
    node_of = {str(G.nodes.node_id[int(x)]): int(x) for x in np.flatnonzero(G.types == "drug")}
    for key, _, _, phase in ids.drugs().values():
        k = key.split(":", 1)[1] if key.startswith("DB:") else key
        if str(phase).upper() == "APPROVAL" and k in node_of and "ChEMBL" not in out.get(node_of[k], []):
            add(node_of[k], "ChEMBL")
    atc = G.descriptions["atc_2"].fillna("").astype(str)
    for i, v in atc.items():
        if v.strip():
            add(i, "ATC")
    return out
