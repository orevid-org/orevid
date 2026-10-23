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
             'betaRatioSignAverage otherStudyLocus { study { id traitFromSource nCases diseases { id } } } } } } } } }')
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


def g2b():
    """G2b: every genetic line of the pages, kept only if release 26.09 has a colocalisation with the stated direction and
    none with the opposite one (same QTL type, cis, h4 ≥ 0.8, a GWAS mapped to one of the line's diseases)."""
    a = pd.read_json(OUT / "annotations.jsonl", lines=True)
    a = a[a.kind == "genetics"].drop_duplicates(["disease", "drug", "g"])
    cache, out = {}, []
    for r in a.itertuples():
        if r.e not in cache:
            cache[r.e] = coloc_rows(r.e)
        types, ids_ = {x.lower() for x in r.q.split("+")}, set(r.ot)
        same, opposite = [], []
        for cs in cache[r.e]:
            if cs["studyType"] not in types or cs["isTransQtl"]:
                continue
            for c in cs["colocalisation"]["rows"]:
                st = c["otherStudyLocus"]["study"]
                if {d["id"] for d in st["diseases"]} & ids_ and c["h4"] is not None and c["h4"] >= 0.8:
                    dis = {d["id"] for d in st["diseases"]}
                    h = {"qtl_study": cs["study"]["id"], "gwas_study": st["id"], "trait": st["traitFromSource"],
                         "h4": round(c["h4"], 3), "cases": st.get("nCases") or 0, "only_page_disease": dis <= ids_}
                    if c["betaRatioSignAverage"] == r.s:
                        same.append(h)
                    elif c["betaRatioSignAverage"] == -r.s:
                        opposite.append(h)
        out.append({"disease_id": int(r.disease), "disease": r.disease_name, "drug_id": int(r.drug), "drug": r.drug_name,
                    "gene": r.g, "ensg": r.e, "sign": int(r.s), "action": int(r.a), "qtl": r.q, "support": same,
                    "opposite": opposite, "stable": bool(same) and not opposite})
    (OUT / "g2b_lines.jsonl").write_text("\n".join(json.dumps(x, ensure_ascii=False) for x in out) + "\n", "utf-8")
    ok = [x for x in out if x["stable"]]
    print(f"{len(out)} lines, {len(ok)} stable, {sum(bool(x['opposite']) for x in out)} with an opposite colocalisation, "
          f"{sum(not x['support'] for x in out)} without support in 26.09")
    import random
    random.Random(71).shuffle(ok)
    if not (OUT / "g2b_sample.json").exists():
        (OUT / "g2b_sample.json").write_text(json.dumps(ok[:N], ensure_ascii=False, indent=1), "utf-8")
    # G2c (P72): a stable line with a case–control supporting GWAS whose diseases are all the page's disease
    for x in out:
        x["g2c"] = x["stable"] and any(h["cases"] > 0 and h["only_page_disease"] for h in x["support"])
    (OUT / "g2b_lines.jsonl").write_text("\n".join(json.dumps(x, ensure_ascii=False) for x in out) + "\n", "utf-8")
    okc = [x for x in out if x["g2c"]]
    rows = []
    for x in okc:  # the lines web/annotate.py may show, with their best supporting study
        best = max((h for h in x["support"] if h["cases"] > 0 and h["only_page_disease"]), key=lambda h: h["h4"])
        rows.append({"disease": x["disease_id"], "drug": x["drug_id"], "gene": x["gene"], "gwas_study": best["gwas_study"],
                     "trait": best["trait"], "h4": best["h4"]})
    pd.DataFrame(rows).to_parquet(ROOT / "data" / "ext" / "genetic_lines_g2c.parquet", index=False)
    random.Random(72).shuffle(okc)
    (OUT / "g2c_sample.json").write_text(json.dumps(okc[:N], ensure_ascii=False, indent=1), "utf-8")
    print(f"G2c: {len(okc)} lines kept")
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
    {"g2": g2, "g2b": g2b, "x1": x1, "r4-sample": r4_sample}[sys.argv[1]]()
