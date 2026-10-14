"""Practical repurposing information per drug (docs/IMPROVEMENT_PLAN.md, I3; display only).

Sources (downloaded to data/raw/practical/):
- FDA Orange Book (products, patents, exclusivity; US public domain): a US generic exists when a single-ingredient
  ANDA product of the active moiety is marketed (Rx or OTC, not discontinued); `exclusivity_until` is the latest
  patent or exclusivity date of a marketed NDA product of that moiety; routes from the dosage form.
- FDA Purple Book (licensed biologics; US public domain): a US biosimilar exists when a 351(k) product of the
  proper name (without the four-letter suffix) is licensed and not discontinued.
- EMA medicines table (© European Medicines Agency, reuse with attribution): an EU generic or biosimilar exists
  when an authorised human medicine of that INN is flagged as generic or biosimilar.
Matching: by active moiety (`evaluation.twins.active_moiety`: name without salts or esters; not by drug family,
so an enantiomer is not credited with its racemate's generics; a mineral salt keeps its full name), on the graph's
drug names and the Open Targets (ChEMBL) name and synonyms of the same molecule. Combination products are ignored. Match rates and a manual check are written to docs/PRACTICAL.md.

    py -3 -m sources.practical download
    py -3 -m sources.practical build      # -> data/ext/practical.parquet (node, generic_us, …)
"""
import io
import re
import sys
import zipfile
from datetime import date
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "practical"
DEST = ROOT / "data" / "ext" / "practical.parquet"
URLS = {
    "orange_book.zip": "https://www.fda.gov/media/76860/download",
    "purple_book.csv": "https://www.accessdata.fda.gov/drugsatfda_docs/PurpleBook/2026/purplebook-search-September-data-download.csv",
    "ema_medicines.xlsx": "https://www.ema.europa.eu/en/documents/report/medicines-output-medicines-report_en.xlsx",
}
HEADERS = {"User-Agent": "Mozilla/5.0 (OREVID research; contact@orevid.org)"}


def download():
    import urllib.request
    RAW.mkdir(parents=True, exist_ok=True)
    for name, url in URLS.items():
        data = urllib.request.urlopen(urllib.request.Request(url, headers=HEADERS), timeout=120).read()
        (RAW / name).write_bytes(data)
        print(f"{name}: {len(data):,} bytes")
    (RAW / "VERSION").write_text(f"downloaded {date.today().isoformat()}\n" + "\n".join(URLS.values()) + "\n", "utf-8")


# A mineral salt whose active moiety would be the bare ion keeps its full name (calcium citrate is not calcium chloride)
IONS = {"calcium", "magnesium", "potassium", "sodium", "zinc", "ferric", "ferrous", "iron", "copper", "cupric", "lithium",
        "aluminum", "aluminium", "ammonium", "chromium", "chromic", "chromous", "manganese", "selenium", "silver", "barium",
        "strontium", "bismuth", "gallium", "cobalt", "gold", "stannous", "sodium fluoride"}


def moiety(name):
    from evaluation.twins import active_moiety
    n = re.sub(r"\s+", " ", str(name).strip().lower())
    m = active_moiety(n)
    return n if m in IONS else m


def synonyms():
    """Graph drug node -> Open Targets name and ChEMBL synonyms, for INN/USAN differences (cromoglicic acid /
    cromolyn sodium, salbutamol / albuterol). Trade names are not used."""
    from sources import ids
    from graph import Graph
    g = Graph()
    key_to_node = {(f"DB:{i}" if s == "DrugBank" else f"CHEMBL:{i}"): int(n)
                   for n, i, s, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_source, g.nodes.node_type)
                   if t == "drug"}
    chembl = {c: v[0] for c, v in ids.drugs().items()}
    m = pd.read_parquet(ids.OT / "drug_molecule", columns=["id", "name", "synonyms"])
    out = {}
    for c, name, syn in zip(m.id, m.name, m.synonyms):
        node = key_to_node.get(chembl.get(c))
        if node is None:
            continue
        # ChEMBL synonyms only: those taken from trial registrations (AACT) include combination partners and bare ions
        labels = [name, *(s["label"] for s in (syn if syn is not None else []) if s["source"] == "ChEMBL")]
        out.setdefault(node, set()).update(x for x in labels if isinstance(x, str) and "/" not in x)
    return out


def orange_book():
    z = zipfile.ZipFile(RAW / "orange_book.zip")
    read = lambda n: pd.read_csv(io.BytesIO(z.read(next(f for f in z.namelist() if f.lower().startswith(n)))), sep="~",
                                 dtype=str, encoding="latin-1")
    prod, pat, exc = read("products"), read("patent"), read("exclusivity")
    prod = prod[~prod.Ingredient.str.contains(";", na=False)]  # single-ingredient products only
    prod = prod.assign(moiety=prod.Ingredient.map(moiety),
                       route=prod["DF;Route"].str.split(";").str[-1].str.replace(r"\s+", " ", regex=True).str.strip().str.lower())
    live = prod[prod.Type.isin(["RX", "OTC"])]
    generic = set(live.moiety[live.Appl_Type == "A"])
    nda = live[live.Appl_Type == "N"][["Appl_No", "Product_No", "moiety"]]
    dates = []
    for frame, col in ((pat, "Patent_Expire_Date_Text"), (exc, "Exclusivity_Date")):
        d = frame.merge(nda, on=["Appl_No", "Product_No"])
        dates.append(pd.DataFrame({"moiety": d.moiety, "until": pd.to_datetime(d[col], errors="coerce")}))
    until = pd.concat(dates).dropna().groupby("moiety").until.max()
    routes = live.groupby("moiety").route.agg(lambda s: sorted(set(s.dropna())))
    names = set(prod.moiety)
    return generic, until, routes, names


def purple_book():
    """The monthly file has a short list of changes and then every licensed product, with the same header."""
    lines = (RAW / "purple_book.csv").read_text("utf-8", errors="replace").splitlines()
    starts = [i for i, l in enumerate(lines) if l.startswith("N/R/U,")]
    pb = pd.read_csv(io.StringIO("\n".join(lines[starts[-1]:])), dtype=str)
    pb = pb[pb["Proper Name"].notna()]
    pb["base"] = pb["Proper Name"].str.lower().str.strip().str.replace(r"-[a-z]{4}$", "", regex=True)
    live = (pb["Licensure"].fillna("") == "Licensed") & ~pb["Marketing Status"].fillna("").str.contains("Disc", case=False)
    biosim = set(pb.base[live & pb["License Type"].fillna("").str.contains("351(k)", regex=False)])
    return biosim, set(pb.base)


def ema():
    t = pd.read_excel(RAW / "ema_medicines.xlsx", header=None, dtype=str)
    hdr = next(i for i in range(30) if t.iloc[i].astype(str).str.contains("Medicine name|Name of medicine", case=False).any())
    t.columns = t.iloc[hdr].astype(str).str.strip()
    t = t.iloc[hdr + 1:]
    col = lambda pat: next(c for c in t.columns if re.search(pat, c, re.I))
    cat, inn, status = col(r"^category"), col(r"non-proprietary|INN|common name"), col(r"medicine status|authorisation status")
    gen, bio = col(r"^generic"), col(r"^biosimilar")
    t = t[(t[cat].str.lower() == "human") & t[status].str.lower().str.startswith("authorised")]
    t = t[~t[inn].str.contains(r"/|,| and ", na=True)]  # single active substance
    t = t.assign(moiety=t[inn].map(moiety))
    yes = lambda c: t[c].astype(str).str.lower().isin(["yes", "y", "true", "1"])
    return set(t.moiety[yes(gen)]), set(t.moiety[yes(bio)]), set(t.moiety)


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    g = Graph()
    drugs = [int(f) for f in g.drugs]
    syn = synonyms()
    ob_generic, ob_until, ob_routes, ob_names = orange_book()
    pb_biosim, pb_names = purple_book()
    eu_generic, eu_biosim, eu_names = ema()
    rows = []
    for f in drugs:
        keys = {moiety(g.name(f))} | {moiety(s) for s in syn.get(f, ()) if s.strip()}
        keys.discard("")
        rows.append({
            "node": f, "name": g.name(f),
            "in_orange_book": bool(keys & ob_names), "in_purple_book": bool(keys & pb_names),
            "in_ema": bool(keys & eu_names),
            "generic_us": bool(keys & ob_generic), "biosimilar_us": bool(keys & pb_biosim),
            "generic_eu": bool(keys & eu_generic), "biosimilar_eu": bool(keys & eu_biosim),
            "exclusivity_until": max((ob_until[k] for k in keys if k in ob_until.index), default=pd.NaT),
            "routes": sorted({r for k in keys if k in ob_routes.index for r in ob_routes[k]}),
        })
    df = pd.DataFrame(rows)
    # the latest patent or exclusivity of a branded product matters only while there is no generic or biosimilar
    # (otherwise it belongs to a newer formulation, e.g. imatinib)
    until = pd.to_datetime(df.exclusivity_until).dt.date.astype(str)
    df["exclusivity_until"] = until.where(~until.isin(["NaT", "NaN", "nan", "None"]), "")
    df.loc[df.generic_us | df.biosimilar_us, "exclusivity_until"] = ""
    df["exclusivity_until"] = df.exclusivity_until.fillna("").replace({"NaT": "", "nan": "", "NaN": ""})
    df["routes"] = df.routes.map(lambda r: ", ".join(r))
    DEST.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(DEST, index=False)
    from predictor import drug_status
    st = drug_status(g)
    appr = df[df.node.map(lambda f: st.get(f) == "approved")]
    print(f"{len(df):,} drugs, {len(appr):,} approved in DrugBank")
    for c in ("in_orange_book", "in_purple_book", "in_ema", "generic_us", "biosimilar_us", "generic_eu", "biosimilar_eu"):
        print(f"  {c:16s} all {df[c].mean():6.1%}   approved {appr[c].mean():6.1%}")
    print(f"-> {DEST}")


if __name__ == "__main__":
    {"download": download, "build": build}[sys.argv[1]]()
