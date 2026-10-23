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
# Negation tied to efficacy or use only («patients who do not respond to…» is an ordinary indication)
NEG = re.compile(r"\b(not|no)\b[^.;]{0,40}\b(effective|efficacy|indicated|established|shown|demonstrated|recommended|"
                 r"approved|active against|be used|intended|evaluated)\b|\bineffective\b", re.I)
SENT = re.compile(r"(?<=[.;])\s+|\n+|•")  # not at «:» (a negation often introduces a list)


# esters and salts frequent in label substance names that evaluation.twins does not strip (triamcinolone acetonide)
EXTRA_SUFFIXES = ("hexacetonide", "acetonide", "diacetate", "aceponate", "nitrate", "pamoate", "oxalate", "lactate",
                  "gluconate", "mesilate", "besilate", "tosylate", "hemifumarate", "hemihydrate", "sesquihydrate",
                  "pentahydrate", "dihydrochloride", "hcl")


def moiety(n):
    from sources.practical import moiety as m
    k = m(n)
    changed = True
    while changed:
        changed = False
        for x in EXTRA_SUFFIXES:
            if k.endswith(" " + x) and len(k) - len(x) - 1 >= 4:
                k = k[: -len(x) - 1].strip(" ,")
                changed = True
    return k


# --- labels ------------------------------------------------------------------------------------

def fda_labels():
    """moiety -> list of (text, is_combination, title, effective_time) of current FDA labels. The 14 openFDA files are
    read as a stream (ijson) and the needed fields cached in data/ext/fda_labels.parquet."""
    cache = ROOT / "data" / "ext" / "fda_labels_v2.parquet"  # v2: with the SPL set id
    if not cache.exists():
        import ijson
        rows = []
        for z in sorted(FDA.glob("drug-label-*.json.zip")):
            with zipfile.ZipFile(z) as zf, zf.open(zf.namelist()[0]) as fh:
                for r in ijson.items(fh, "results.item"):
                    text = " ".join(r.get("indications_and_usage") or [])
                    if not text:
                        continue
                    o = r.get("openfda") or {}
                    subs = list(o.get("substance_name") or [])
                    gens = list(o.get("generic_name") or [])
                    rows.append({"substances": "|".join(subs), "generic": "|".join(gens),
                                 "title": (o.get("brand_name") or [""])[0], "time": r.get("effective_time", ""),
                                 "set_id": r.get("set_id", ""), "text": text[:6000]})
            print(f"  {z.name}: {len(rows):,} labels with indications so far", flush=True)
        pd.DataFrame(rows).to_parquet(cache, index=False)
    df = pd.read_parquet(cache)
    out = defaultdict(list)
    for subs, gens, title, sid, text in zip(df.substances, df.generic, df.title, df.set_id, df.text):
        sl = [x for x in subs.split("|") if x]
        gl = [x for x in gens.split("|") if x]
        names = {moiety(x) for x in sl + gl} - {""}
        combo = len(sl) > 1 or any(re.search(r"\band\b|,|/", x, re.I) for x in gl)
        for n in names:
            out[n].append((text, combo, title, f"https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={sid}"))
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
    for name, ind, status, med, url in zip(t[inn], t["Therapeutic indication"], t["Medicine status"], t["Name of medicine"],
                                           t["Medicine URL"]):
        if not isinstance(name, str) or not isinstance(ind, str) or not str(status).lower().startswith("authorised"):
            continue
        parts = [p for p in re.split(r"\s*(?:/|,|\band\b)\s*", name) if p.strip()]
        combo = len(parts) > 1
        for p in parts:
            out[moiety(p)].append((ind, combo, med, f"EMA {url}"))
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
        m = moiety(re.sub(r",?\s*(ultra)?microsize", "", ing, flags=re.I))
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


CUE = re.compile(r"\b(indicated|indication|treatment|treat|management|prevention|prevent|prophylaxis|relief|relieve|"
                 r"reduce the risk|therapy|adjunct|control of)\b", re.I)
QUAL = re.compile(r"\btype\s*(1|2|i|ii)\b", re.I)


def _pat(s):
    return r"(?<![a-z])" + re.escape(s) + r"(?:s|es)?(?![a-z])"


def _qualifier_ok(sentence, names):
    """A disease named with a type (type 1 diabetes) is not confirmed by a sentence about the other type only."""
    want = {m.group(1).lower().replace("ii", "2").replace("i", "1") for n in names for m in QUAL.finditer(n)}
    if not want:
        return True
    got = {m.group(1).lower().replace("ii", "2").replace("i", "1") for m in QUAL.finditer(sentence)}
    return not got or bool(want & got)


def find(text, names, broader=False):
    """(found, negated, confirming sentence): the disease named in `text`; negated if every sentence naming it is
    negated; confirmed only by a sentence where it follows an indication cue («indicated for», «treatment of»…) within
    120 characters. A broader term preceded by «forms of» («relapsing forms of multiple sclerosis») does not confirm."""
    low = text.lower()
    hits = [s for s in names if re.search(_pat(s), low)]
    if not hits:
        return False, False, ""
    neg_all, seen = True, False
    for sent in SENT.split(text):
        sl = sent.lower()
        spans = [(m.start(), s) for s in hits for m in re.finditer(_pat(s), sl)]
        if broader:
            spans = [(i, s) for i, s in spans if "forms of" not in sl[max(0, i - 12):i]]
        if not spans:
            continue
        seen = True
        if NEG.search(sent):
            continue
        neg_all = False
        cues = [m.end() for m in CUE.finditer(sent)]
        close = [s for i, s in spans if any(0 <= i - c <= 120 for c in cues)]
        if close and _qualifier_ok(sent, close):
            return True, False, sent.strip()[:300]
    if not seen:
        return False, False, ""
    return True, neg_all, ""


def classify(drug_keys, dnames, fda, ema, broader=None):
    """Category for a pair; if the disease is not named, the names of its direct broader terms are tried
    (a label for «relapsing forms of multiple sclerosis» covers relapsing-remitting MS), marked in `where`."""
    cat, where = _classify(drug_keys, dnames, fda, ema)
    if cat == "not found in a current label" and broader:
        c2, w2 = _classify(drug_keys, broader, fda, ema, broader=True)
        if c2 != "not found in a current label":
            return c2, w2 + " (broader term)"
    return cat, where


def _classify(drug_keys, dnames, fda, ema, broader=False):
    labels = [l for k in drug_keys for l in fda.get(k, [])] + [l for k in drug_keys for l in ema.get(k, [])]
    if not labels:
        return "no current label", ""
    single = [l for l in labels if not l[1]]
    combo = [l for l in labels if l[1]]
    neg_seen = None
    for group in (single, combo):
        for text, is_combo, title, src in group:
            found, neg, sent = find(text, dnames, broader)
            if sent:
                return ("confirmed" if not is_combo else "only in a combination product"), f"{title} | {src} | {sent}"
            if found and neg:
                neg_seen = neg_seen or f"{title} | {src}"
    if neg_seen:
        return "negated", neg_seen
    return "not found in a current label", ""


def dailymed_count(name, cache):
    """Number of SPLs DailyMed lists for a drug name (any kind: a conservative test for «no US product»)."""
    import time
    import urllib.parse
    import urllib.request
    key = name.lower().strip()
    if key not in cache:
        url = ("https://dailymed.nlm.nih.gov/dailymed/services/v2/spls.json?pagesize=1&drug_name="
               + urllib.parse.quote(key))
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (OREVID research; contact@orevid.org)"})
        try:
            cache[key] = int(json.loads(urllib.request.urlopen(req, timeout=60).read())["metadata"]["total_elements"])
        except Exception:
            cache[key] = -1  # unknown: then no claim is made
        time.sleep(0.2)
    return cache[key]


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
    dm_file = ROOT / "data" / "raw" / "dailymed" / "spl_counts.json"
    dm_file.parent.mkdir(parents=True, exist_ok=True)
    dm = json.loads(dm_file.read_text("utf-8")) if dm_file.exists() else {}
    rows = []
    for f, ks in keys.items():
        u = next((us[k] for k in ks if k in us), "")
        if any(fda.get(k) for k in ks):  # a current FDA label (SPL listing): a product on the US market
            u = "marketed"
        if u == "discontinued":  # P55: openFDA lacks some current labels; DailyMed must have none under the name
            names = {str(g.name(f))} | set(ks)
            if any(dailymed_count(n, dm) != 0 for n in names):
                u = ""
        e = next((eu[k] for k in ks if k in eu), ("", ""))
        rows.append({"node": f, "name": g.name(f), "us_status": u, "eu_status": e[0], "eu_year": e[1],
                     "fda_labels": int(sum(len(fda.get(k, [])) for k in ks)), "ema_labels": int(sum(len(ema.get(k, [])) for k in ks))})
    dm_file.write_text(json.dumps(dm), "utf-8")
    st = pd.DataFrame(rows)
    st.to_parquet(ROOT / "data" / "ext" / "drug_status.parquet", index=False)
    # pairs the website calls approved
    dn = disease_names(g)
    from sources.hierarchy import load as load_hierarchy
    h = load_hierarchy()
    parents = defaultdict(list)
    for a_, b_ in zip(h.parent.astype(int), h.child.astype(int)):
        parents[b_].append(a_)
    I = g.M(INDICATION).tocsr()
    pairs = {(int(f), int(d), "PrimeKG") for f in g.drugs for d in I.indices[I.indptr[f]:I.indptr[f + 1]]
             if g.types[d] == "disease"}
    ot = pd.read_parquet(ROOT / "data" / "ext" / "ot_approvals.parquet")
    for x, y, s in zip(ot.x_index.astype(int), ot.y_index.astype(int), ot.status):
        pairs.add((x, y, f"Open Targets ({s})"))
    out = []
    for f, d, src in sorted(pairs):
        broader = set().union(*(dn.get(p_, set()) for p_ in parents.get(d, ()))) if parents.get(d) else None
        cat, where = classify(keys.get(f, set()), dn.get(d, {str(g.name(d)).lower()}), fda, ema, broader)
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
