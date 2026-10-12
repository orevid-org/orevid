"""Converts Open Targets into edges with canonical keys (data/ext/).

    py -3 -m sources.opentargets

New relations:
    moa              drug -> gene   (display: inhibits / activates / other)          ChEMBL
    clinical         drug -> disease (display: approved / phase 3 / 2 / 1 …)        ChEMBL, ClinicalTrials.gov, DailyMed…
    clinical_negative drug -> disease (display: lack of efficacy / safety)          stopped trials
    faers            drug -> adverse effect (statistical signal from FAERS)
    ot_assoc         disease -> gene (display: evidence type)
    disease_disease  new disease -> parent (to hang it in the ontology)

The "known_drug" associations are discarded: they are drug indications and
would contaminate any validation.
"""
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.dataset as ds

from sources import ids

OT = ids.OT
EXT = Path(__file__).parent.parent / "data" / "ext"

PHASE = {
    "APPROVAL": "approved", "PREAPPROVAL": "approved",
    "PHASE_3": "phase 3", "PHASE_2_3": "phase 3",
    "PHASE_2": "phase 2", "PHASE_1_2": "phase 2",
    "PHASE_1": "phase 1", "EARLY_PHASE_1": "phase 1",
    "IND": "preclinical", "PRECLINICAL": "preclinical", "UNKNOWN": "unknown",
}
PHASE_ORDER = ["unknown", "preclinical", "phase 1", "phase 2", "phase 3", "approved"]

INHIBITS = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "NEGATIVE MODULATOR",
            "INVERSE AGONIST", "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "ALLOSTERIC ANTAGONIST", "DEGRADER",
            "DISRUPTING AGENT", "GENE EDITING NEGATIVE MODULATOR"}
ACTIVATES = {"AGONIST", "PARTIAL AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR", "POSITIVE MODULATOR",
             "OPENER", "STABILISER", "EXOGENOUS PROTEIN", "EXOGENOUS GENE", "RELEASING AGENT"}

# Score threshold per evidence type (≈ upper quartile of each distribution)
THRESHOLDS = {
    "genetic_association": 0.3, "somatic_mutation": 0.3, "genetic_literature": 0.3,
    "affected_pathway": 0.3, "animal_model": 0.5, "literature": 0.25, "rna_expression": 0.2,
}


class Nodes:
    """Collects the nodes used by the edges, with type and name."""

    def __init__(self):
        self.rows = {}
        self.primekg = ids.all_primekg_keys()

    def add(self, key, node_type, name, source):
        if key not in self.primekg and key not in self.rows:
            self.rows[key] = (node_type, name, source)

    def table(self):
        return pd.DataFrame([(k, *v) for k, v in self.rows.items()],
                            columns=["key", "node_type", "node_name", "node_source"])


def drug(cid, nodes):
    f = ids.drugs().get(cid)
    if f is None:
        return None
    key, name, _, _ = f
    nodes.add(key, "drug", name, "ChEMBL")
    return key


def disease(did, nodes):
    mapping, names, _, _ = ids.diseases()
    key = mapping.get(did)
    if key is None:
        return None
    node_type = "effect/phenotype" if key.startswith("HPO:") else "disease"
    nodes.add(key, node_type, names.get(did, did), "Open Targets")
    return key


def gene(ensg, nodes):
    g = ids.genes().get(ensg)
    if g is None:
        return None
    nodes.add(g[0], "gene/protein", g[1], "HGNC")
    return g[0]


def mechanisms(nodes):
    m = pd.read_parquet(OT / "drug_mechanism_of_action", columns=["actionType", "chemblIds", "targets"])
    m = m.explode("chemblIds").explode("targets").dropna(subset=["chemblIds", "targets"])
    m["x"] = [drug(c, nodes) for c in m.chemblIds]
    m["y"] = [gene(t, nodes) for t in m.targets]
    m["display_relation"] = np.where(m.actionType.isin(INHIBITS), "inhibits",
                                     np.where(m.actionType.isin(ACTIVATES), "activates", "other"))
    m = m.dropna(subset=["x", "y"]).rename(columns={"actionType": "action"})
    return m.assign(relation="moa")[["x", "y", "relation", "display_relation", "action"]].drop_duplicates(["x", "y", "display_relation"])


def clinical(nodes):
    ci = pd.read_parquet(OT / "clinical_indication", columns=["drugId", "diseaseId", "maxClinicalStage"])

    r = pd.read_parquet(OT / "clinical_report")
    if "origin" in r:  # 26.09: type = INDICATION/SAFETY, origin = where it comes from
        r = r[(r.type == "INDICATION") & (r.origin == "CLINICAL_TRIAL")]
    else:              # 26.06: type = where it comes from
        r = r[r.type == "CLINICAL_TRIAL"]
    r = r[["drugs", "diseases", "trialStartDate", "trialStopReasonCategories"]]
    year = pd.to_datetime(r.trialStartDate).dt.year
    r = r.assign(
        year=year.where(year.between(1950, 2026)),
        negative=[c is not None and "Negative" in c for c in r.trialStopReasonCategories],
        safety=[c is not None and "Safety_Sideeffects" in c for c in r.trialStopReasonCategories],
    )
    r = r.explode("drugs").explode("diseases").dropna(subset=["drugs", "diseases"])
    r["drugId"] = [d["drugId"] for d in r.drugs]
    r["diseaseId"] = [d["diseaseId"] for d in r.diseases]
    agg = r.groupby(["drugId", "diseaseId"]).agg(
        year=("year", "min"), n_trials=("year", "size"),
        n_negative=("negative", "sum"), n_safety=("safety", "sum")).reset_index()

    c = ci.merge(agg, on=["drugId", "diseaseId"], how="left")
    c["x"] = [drug(d, nodes) for d in c.drugId]
    c["y"] = [disease(d, nodes) for d in c.diseaseId]
    c = c.dropna(subset=["x", "y"])
    c["display_relation"] = c.maxClinicalStage.map(PHASE).fillna("unknown")
    c["rank"] = c.display_relation.map(PHASE_ORDER.index)
    # Several molecules or terms can land on the same pair: the highest phase is kept.
    c = (c.sort_values("rank", ascending=False)
          .groupby(["x", "y"], as_index=False)
          .agg(display_relation=("display_relation", "first"), year=("year", "min"),
               n_trials=("n_trials", "sum"), n_negative=("n_negative", "sum"),
               n_safety=("n_safety", "sum")))
    clinical_edges = c.assign(relation="clinical")

    neg = []
    for reason, col in (("lack of efficacy", "n_negative"), ("safety", "n_safety")):
        s = c[(c[col] > 0) & (c.display_relation != "approved")]
        neg.append(s.assign(relation="clinical_negative", display_relation=reason))
    return pd.concat([clinical_edges, *neg], ignore_index=True)


def faers(nodes):
    f = pd.read_parquet(OT / "openfda_significant_adverse_drug_reactions")
    phenotypes = {name.lower(): k for k, (t, name) in nodes.primekg.items() if t == "effect/phenotype"}
    f["x"] = [drug(c, nodes) for c in f.chembl_id]
    ys = []
    for ev, code in zip(f.event, f.meddraCode):
        key = phenotypes.get(ev.lower()) or f"MEDDRA:{code or ev.lower()}"
        nodes.add(key, "effect/phenotype", ev.capitalize(), "MedDRA (FAERS)")
        ys.append(key)
    f["y"] = ys
    f = f.dropna(subset=["x"]).rename(columns={"llr": "score"})
    return f.assign(relation="faers", display_relation="adverse effect")[
        ["x", "y", "relation", "display_relation", "score"]].drop_duplicates(["x", "y"])


def associations(nodes):
    d = ds.dataset(OT / "association_by_datatype_direct")
    rows = []
    for evidence_type, threshold in THRESHOLDS.items():
        t = d.to_table(
            columns=["diseaseId", "targetId", "associationScore", "timeseries"],
            filter=(ds.field("aggregationValue") == evidence_type) & (ds.field("associationScore") >= threshold),
        ).to_pandas()
        # First year with evidence: will be used for validation over time.
        t["year"] = [min((p["year"] for p in s if p["yearlyEvidenceCount"] and p["year"] is not None), default=np.nan)
                     if s is not None else np.nan for s in t.timeseries]
        t["display_relation"] = evidence_type
        rows.append(t.drop(columns="timeseries"))
        print(f"    {evidence_type}: {len(t):,}", flush=True)
    a = pd.concat(rows, ignore_index=True)
    disease_keys = {d: disease(d, nodes) for d in a.diseaseId.unique()}
    gene_keys = {g: gene(g, nodes) for g in a.targetId.unique()}
    a["x"] = a.diseaseId.map(disease_keys)
    a["y"] = a.targetId.map(gene_keys)
    a = a.dropna(subset=["x", "y"]).rename(columns={"associationScore": "score"})
    a = a.sort_values("score", ascending=False).drop_duplicates(["x", "y", "display_relation"])
    return a.assign(relation="ot_assoc")[["x", "y", "relation", "display_relation", "score", "year"]]


def ontology(nodes):
    """Hangs the new diseases (without a PrimeKG node) from their parents."""
    mapping, names, parents, _ = ids.diseases()
    new = {k for k, (t, _, _) in nodes.rows.items() if t == "disease"}
    inverse = {}
    for did, key in mapping.items():
        inverse.setdefault(key, did)
    rows = []
    for key in new:
        did = inverse.get(key) if not key.startswith("OT:") else key[3:]
        for p in parents.get(did, []):
            cp = mapping.get(p)
            if cp and cp != key and (cp in nodes.primekg or cp in nodes.rows):
                rows.append((key, cp))
    return pd.DataFrame(rows, columns=["x", "y"]).assign(relation="disease_disease", display_relation="parent-child")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    EXT.mkdir(parents=True, exist_ok=True)
    nodes = Nodes()
    _, _, _, via = ids.diseases()
    print("Open Targets diseases ->", dict(via))
    parts = {}
    print("mechanisms…", flush=True); parts["moa"] = mechanisms(nodes)
    print("clinical…", flush=True); parts["clinical"] = clinical(nodes)
    print("FAERS…", flush=True); parts["faers"] = faers(nodes)
    print("associations…", flush=True); parts["assoc"] = associations(nodes)
    parts["ontology"] = ontology(nodes)
    # VERSION is the version file written by sources/download.py into data/raw.
    edges = pd.concat(parts.values(), ignore_index=True).assign(source=f"Open Targets {(OT / "VERSION").read_text().strip()}")
    edges.to_parquet(EXT / "opentargets_edges.parquet", index=False)
    nodes.table().to_parquet(EXT / "opentargets_nodes.parquet", index=False)

    print("\nEdges:")
    print(edges.groupby(["relation", "display_relation"]).size().to_string())
    print("\nNew nodes:", Counter((t, f) for t, _, f in nodes.rows.values()))


if __name__ == "__main__":
    main()
