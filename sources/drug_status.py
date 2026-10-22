"""N4: drug status audit against current regulatory sources (docs/IMPROVEMENT_PLAN.md, N4; display only).

1. Per drug: current US status (FDA Orange Book and Purple Book: a marketed product, or only discontinued ones) and EU
   status (EMA medicines table: authorised, or withdrawn / refused / suspended / expired with the year).
2. Per «approved» pair the website shows (PrimeKG indications and Open Targets approvals): is the disease named in the
   indication text of a current FDA label (openFDA drug labels, «indications and usage») or of the EMA medicine
   («therapeutic indication»), and is that sentence negated («not effective against…», «has not been established»)?
   Categories: confirmed · negated · only in a combination product · not found in a current label · no current label.
Matching of drugs by active moiety on the graph name and the ChEMBL synonyms (as sources/practical.py); of diseases by
the graph name, the names of its twin nodes and the Open Targets exact synonyms of the same MONDO term.

    py -3 -m sources.drug_status build     # -> data/ext/drug_status.parquet, data/ext/pair_label.parquet
"""
import json
import re
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
FDA = ROOT / "data" / "raw" / "openfda"
PRACT = ROOT / "data" / "raw" / "practical"
NEG = re.compile(r"\b(not|no|never)\b[^.;]{0,60}\b(effective|efficacy|indicated|established|shown|demonstrated|recommended|"
                 r"approved|active|be used|for use|intended)\b|\bineffective\b|\bhas not been\b|\bhave not been\b|"
                 r"\bis not\b|\bare not\b|\bshould not\b|\bdoes not\b|\bdo not\b", re.I)
SENT = re.compile(r"(?<=[.;:])\s+|\n+|•")


def moiety(n):
    from sources.practical import moiety as m
    return m(n)


# --- labels ------------------------------------------------------------------------------------

def fda_labels():
    """moiety -> list of (text, is_combination, title, effective_time) of current FDA labels."""
    out = defaultdict(list)
    for z in sorted(FDA.glob("drug-label-*.json.zip")):
        with zipfile.ZipFile(z) as zf:
            data = json.loads(zf.read(zf.namelist()[0]))
        for r in data["results"]:
            text = " ".join(r.get("indications_and_usage") or [])
            if not text:
                continue
            o = r.get("openfda") or {}
            subs = o.get("substance_name") or o.get("generic_name") or []
            names = {moiety(s) for s in subs} | {moiety(s) for s in (o.get("generic_name") or [])}
            names.discard("")
            gen = " ".join(o.get("generic_name") or [])
            combo = len(o.get("substance_name") or []) > 1 or bool(re.search(r"\band\b|,|/", gen))
            title = (o.get("brand_name") or [""])[0]
            for n in names:
                out[n].append((text, combo, title, r.get("effective_time", "")))
    return out


def ema_table():
    t = pd.read_excel(PRACT / "ema_medicines.xlsx", header=None, dtype=str)
    hdr = next(i for i in range(30) if t.iloc[i].astype(str).str.contains("Name of medicine", case=False).any())
    t.columns = t.iloc[hdr].astype(str).str.strip()
    t = t.iloc[hdr + 1:]
    return t[t["Category"].str.lower() == "human"]


def ema_labels(t):
    out = defaultdict(list)
    inn = "International non-proprietary name (INN) / common name"
    for name, ind, status, med in zip(t[inn], t["Therapeutic indication"], t["Medicine status"], t["Name of medicine"]):
        if not isinstance(name, str) or not isinstance(ind, str) or not str(status).lower().startswith("authorised"):
            continue
        parts = [p for p in re.split(r"\s*(?:/|,|\band\b)\s*", name) if p.strip()]
        combo = len(parts) > 1
        for p in parts:
            out[moiety(p)].append((ind, combo, med, "EMA"))
    return out


def ema_status(t):
    """moiety -> (EU status, year) from the single-substance human medicines of that INN."""
    inn = "International non-proprietary name (INN) / common name"
    wcol = "Withdrawal / expiry / revocation / lapse of marketing authorisation date"
    res = {}
    for name, status, wdate, rdate in zip(t[inn], t["Medicine status"], t[wcol], t["Refusal of marketing authorisation date"]):
        if not isinstance(name, str) or re.search(r"/|,|\band\b", name):
            continue
        m, st = moiety(name), str(status).strip().lower()
        year = next((str(x)[:4] for x in (wdate, rdate) if isinstance(x, str) and x[:4].isdigit()), "")
        rank = 0 if st.startswith("authorised") else 1
        cur = res.get(m)
        if cur is None or rank < cur[2]:
            res[m] = (st, year, rank)
    return {k: (v[0], v[1]) for k, v in res.items()}


def us_status():
    from sources.practical import RAW
    import io
    z = zipfile.ZipFile(RAW / "orange_book.zip")
    prod = pd.read_csv(io.BytesIO(z.read(next(f for f in z.namelist() if f.lower().startswith("products")))), sep="~",
                       dtype=str, encoding="latin-1")
    prod = prod[~prod.Ingredient.str.contains(";", na=False)]
    st = {}
    for ing, typ in zip(prod.Ingredient, prod.Type):
        m = moiety(ing)
        st[m] = st.get(m, False) or typ in ("RX", "OTC")
    lines = (RAW / "purple_book.csv").read_text("utf-8", errors="replace").splitlines()
    start = [i for i, l in enumerate(lines) if l.startswith("N/R/U,")][-1]
    pb = pd.read_csv(io.StringIO("\n".join(lines[start:])), dtype=str)
    pb = pb[pb["Proper Name"].notna()]
    for name, lic, mk in zip(pb["Proper Name"], pb["Licensure"].fillna(""), pb["Marketing Status"].fillna("")):
        m = moiety(re.sub(r"-[a-z]{4}$", "", name.lower().strip()))
        live = lic == "Licensed" and "disc" not in mk.lower()
        st[m] = st.get(m, False) or live
    return {k: ("marketed" if v else "discontinued") for k, v in st.items()}


# --- matching ------------------------------------------------------------------------------------

def disease_names(g):
    """graph disease node -> set of lower-case names (graph name, twin nodes, Open Targets exact synonyms)."""
    from evaluation.twins import twins
    from sources.ids import diseases
    names = defaultdict(set)
    tw = twins(g, "disease")
    for n in np.flatnonzero(g.types == "disease"):
        for t in tw.get(int(n), {int(n)}):
            names[int(n)].add(str(g.names[t]).lower())
    key_node = {}
    for n, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type):
        if t == "disease" and str(i).replace("_", "").isdigit():
            for p in str(i).split("_"):
                key_node.setdefault(f"MONDO:{int(p)}", int(n))
    mapping = diseases()[0] if isinstance(diseases(), tuple) else diseases()
    ot = pd.read_parquet(ROOT / "data" / "raw" / "opentargets" / "disease", columns=["id", "name", "synonyms"])
    for did, name, syn in zip(ot.id, ot.name, ot.synonyms):
        node = key_node.get(mapping.get(did, "")) if isinstance(mapping, dict) else None
        if node is None:
            continue
        names[node].add(str(name).lower())
        if syn is not None and syn.get("hasExactSynonym") is not None:
            names[node].update(str(s).lower() for s in syn["hasExactSynonym"])
    return {k: {s for s in v if len(s) >= 4} for k, v in names.items()}


def find(text, names):
    """(found, negated): the disease named in a sentence of `text`; negated if every such sentence is negated."""
    low = text.lower()
    hits = [s for s in names if re.search(r"(?<![a-z])" + re.escape(s) + r"(?![a-z])", low)]
    if not hits:
        return False, False
    neg_all = True
    for sent in SENT.split(text):
        sl = sent.lower()
        if any(re.search(r"(?<![a-z])" + re.escape(s) + r"(?![a-z])", sl) for s in hits):
            if not NEG.search(sent):
                neg_all = False
                break
    return True, neg_all


def classify(drug_keys, dnames, fda, ema):
    labels = [l for k in drug_keys for l in fda.get(k, [])] + [l for k in drug_keys for l in ema.get(k, [])]
    if not labels:
        return "no current label", ""
    single = [l for l in labels if not l[1]]
    combo = [l for l in labels if l[1]]
    neg_seen = None
    for group in (single, combo):
        for text, is_combo, title, src in group:
            found, neg = find(text, dnames)
            if found and not neg:
                return ("confirmed" if not is_combo else "only in a combination product"), f"{title} ({'EMA' if src == 'EMA' else 'FDA'})"
            if found and neg:
                neg_seen = neg_seen or f"{title} ({'EMA' if src == 'EMA' else 'FDA'})"
    if neg_seen:
        return "negated", neg_seen
    return "not found in a current label", ""


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import INDICATION, Graph
    from sources.practical import synonyms
    g = Graph()
    syn = synonyms()
    keys = {int(f): ({moiety(g.name(int(f)))} | {moiety(s) for s in syn.get(int(f), ())}) - {""} for f in g.drugs}
    fda = fda_labels()
    t = ema_table()
    ema = ema_labels(t)
    eu, us = ema_status(t), us_status()
    rows = []
    for f, ks in keys.items():
        u = next((us[k] for k in ks if k in us), "")
        e = next((eu[k] for k in ks if k in eu), ("", ""))
        rows.append({"node": f, "name": g.name(f), "us_status": u, "eu_status": e[0], "eu_year": e[1],
                     "fda_labels": int(sum(len(fda.get(k, [])) for k in ks)), "ema_labels": int(sum(len(ema.get(k, [])) for k in ks))})
    st = pd.DataFrame(rows)
    st.to_parquet(ROOT / "data" / "ext" / "drug_status.parquet", index=False)
    # pairs the website calls approved
    dn = disease_names(g)
    I = g.M(INDICATION).tocsr()
    pairs = {(int(f), int(d), "PrimeKG") for f in g.drugs for d in I.indices[I.indptr[f]:I.indptr[f + 1]]
             if g.types[d] == "disease"}
    ot = pd.read_parquet(ROOT / "data" / "ext" / "ot_approvals.parquet")
    for x, y, s in zip(ot.x_index.astype(int), ot.y_index.astype(int), ot.status):
        pairs.add((x, y, f"Open Targets ({s})"))
    out = []
    for f, d, src in sorted(pairs):
        cat, where = classify(keys.get(f, set()), dn.get(d, {str(g.name(d)).lower()}), fda, ema)
        out.append({"drug": f, "disease": d, "drug_name": g.name(f), "disease_name": g.name(d), "source": src,
                    "label": cat, "where": where})
    pl = pd.DataFrame(out)
    pl.to_parquet(ROOT / "data" / "ext" / "pair_label.parquet", index=False)
    from predictor import drug_status
    ds = drug_status(g)
    appr = st[st.node.map(lambda f: ds.get(f) == "approved")]
    print(f"{len(st):,} drugs; approved in DrugBank {len(appr):,}: US status {(appr.us_status != '').mean():.1%}, "
          f"EU status {(appr.eu_status != '').mean():.1%}, a current FDA label {(appr.fda_labels > 0).mean():.1%}, "
          f"an EMA indication {(appr.ema_labels > 0).mean():.1%}")
    print(pl.groupby(["source", "label"]).size().unstack(fill_value=0))


if __name__ == "__main__":
    {"build": build}[sys.argv[1]]()
