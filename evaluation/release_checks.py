"""Release checks of G2 (genetic direction), X1 (CPIC) and R4 (Orphanet prevalence) (docs/IMPROVEMENT_PLAN.md,
«Registered 23-10-2026»): 50 random lines of `data/eval/annotations.jsonl` (`py -3 -m web.annotate --dump`) each, checked
against the source and kept with the answer. R4 is checked by opening orpha.net pages (browser); this script writes its
sample and reads the verdicts back.

    py -3 -m evaluation.release_checks g2|x1|r4-sample|report
"""
import json
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "data" / "eval"
N, SEED = 50, 62
OT_API = "https://api.platform.opentargets.org/api/v4/graphql"
LOWER = {"INHIBITOR", "ANTAGONIST", "BLOCKER", "NEGATIVE ALLOSTERIC MODULATOR", "INVERSE AGONIST", "DEGRADER",
         "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "DISRUPTING AGENT"}
RAISE = {"AGONIST", "ACTIVATOR", "POSITIVE ALLOSTERIC MODULATOR", "OPENER", "STABILISER", "PARTIAL AGONIST"}


def sample(kind):
    a = pd.read_json(OUT / "annotations.jsonl", lines=True)
    a = a[a.kind == kind]
    key = {"genetics": ["disease", "drug", "g"], "pgx": ["drug", "g"], "prevalence": ["disease"]}[kind]
    a = a.drop_duplicates(key)
    return a.sample(min(N, len(a)), random_state=SEED)


def ot(query):
    req = urllib.request.Request(OT_API, json.dumps({"query": query}).encode(), {"Content-Type": "application/json"})
    for k in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as r:
                j = json.load(r)
            if "errors" in j:
                raise RuntimeError(j["errors"])
            return j["data"]
        except Exception:
            if k == 3:
                raise
            time.sleep(5)


def coloc_rows(ensg, size=10):
    """Every QTL credible set of the gene with its GWAS colocalisations (small pages: the API refuses large ones)."""
    rows, page = [], 0
    while True:
        q = (f'query {{ target(ensemblId:"{ensg}") {{ credibleSets(page:{{index:{page},size:{size}}}) {{ count rows {{ '
             'studyType isTransQtl study { id } colocalisation(studyTypes:[gwas], page:{index:0,size:200}) { rows { h4 '
             'betaRatioSignAverage otherStudyLocus { study { id diseases { id } } } } } } } } }')
        cs = ot(q)["target"]["credibleSets"]
        rows += cs["rows"]
        page += 1
        if page * size >= cs["count"]:
            return rows


def g2():
    from sources import ids
    s = sample("genetics")
    chembl = {}
    for cid, (key, *_) in ids.drugs().items():
        chembl.setdefault(key, cid)
    from graph import Graph
    g = Graph()
    n = g.nodes
    node_key = {int(k): (f"DB:{i}" if src == "DrugBank" else f"CHEMBL:{i}") for k, i, t, src in
                zip(n.node_index, n.node_id, n.node_type, n.node_source) if t == "drug"}
    out = []
    for r in s.itertuples():
        types = {x.lower() for x in r.q.split("+")}
        ids_ = set(r.ot)
        hits = []
        for cs in coloc_rows(r.e):
            if cs["studyType"] not in types or cs["isTransQtl"]:
                continue
            for c in cs["colocalisation"]["rows"]:
                ds = {d["id"] for d in c["otherStudyLocus"]["study"]["diseases"]}
                if ds & ids_ and c["h4"] is not None and c["h4"] >= 0.8:
                    hits.append({"qtl_study": cs["study"]["id"], "gwas_study": c["otherStudyLocus"]["study"]["id"],
                                 "h4": round(c["h4"], 3), "sign": c["betaRatioSignAverage"]})
        signs = {h["sign"] for h in hits}
        genetics_ok = r.s in signs  # registered: a colocalisation (h4 ≥ 0.8) with the stated direction exists
        opposite = -r.s in signs
        cid = chembl.get(node_key.get(int(r.drug)), "")
        moa = ot(f'query {{ drug(chemblId:"{cid}") {{ mechanismsOfAction {{ rows {{ actionType targets {{ id }} }} }} }} }}')["drug"] if cid else None
        acts = {m["actionType"] for m in (moa or {}).get("mechanismsOfAction", {}).get("rows", []) if any(x["id"] == r.e for x in m["targets"])}
        drug_ok = bool(acts & (LOWER if r.a < 0 else RAISE))
        out.append({"disease": r.disease_name, "drug": r.drug_name, "gene": r.g, "ensg": r.e, "stated_sign": int(r.s),
                    "stated_action": "lowers" if r.a < 0 else "raises", "qtl": r.q, "colocalisations": hits,
                    "genetics_ok": genetics_ok, "opposite_also": opposite, "chembl": cid, "actions": sorted(acts), "drug_ok": drug_ok,
                    "wrong": not (genetics_ok and drug_ok)})
        print(f"  {r.drug_name} | {r.disease_name} | {r.g}: genetics {'ok' if genetics_ok else 'NO'} ({len(hits)} coloc, signs "
              f"{sorted(signs)}), drug {'ok' if drug_ok else 'NO'} {sorted(acts)}", flush=True)
    (OUT / "g2_check.jsonl").write_text("\n".join(json.dumps(x, ensure_ascii=False) for x in out) + "\n", "utf-8")
    return out


def x1():
    s = sample("pgx")
    pairs = json.loads(urllib.request.urlopen("https://api.cpicpgx.org/v1/pair", timeout=120).read())
    drugs = {d["drugid"]: d for d in json.loads(urllib.request.urlopen("https://api.cpicpgx.org/v1/drug", timeout=120).read())}
    from graph import Graph
    g = Graph()
    out = []
    for r in s.itertuples():
        f = int(r.drug)
        db = str(g.nodes.node_id[f])
        live = [p for p in pairs if p["genesymbol"] == r.g and drugs.get(p["drugid"], {}).get("drugbankid") == db and not p["removed"]]
        lvl = live[0]["cpiclevel"] if live else None
        name = drugs.get(live[0]["drugid"], {}).get("name") if live else None
        ok = bool(live) and lvl == r.l
        out.append({"drug": r.drug_name, "drugbank": db, "gene": r.g, "stated_level": r.l, "cpic_level": lvl,
                    "cpic_drug_name": name, "wrong": not ok})
    (OUT / "x1_check.jsonl").write_text("\n".join(json.dumps(x, ensure_ascii=False) for x in out) + "\n", "utf-8")
    return out


def r4_sample():
    s = sample("prevalence")
    rows = [{"disease": r.disease_name, "orpha": str(r.o), "orpha_name": r.n, "stated": r.c,
             "url": f"https://www.orpha.net/en/disease/detail/{r.o}"} for r in s.itertuples()]
    (OUT / "r4_sample.json").write_text(json.dumps(rows, ensure_ascii=False, indent=1), "utf-8")
    return rows


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    {"g2": g2, "x1": x1, "r4-sample": r4_sample}[sys.argv[1]]()
