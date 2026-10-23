"""N4 second version: the three parts released one by one (docs/IMPROVEMENT_PLAN.md, N4, «Added 23-10-2026»), each with
its own check against the source pages (gate ≤ 5 % wrong, none of the 22-10 cases wrong):

  (a) «negated» pairs among the approvals the website shows — all of them;
  (b) drugs shown on the website that are withdrawn / revoked / refused / suspended / expired / lapsed in the EU, or
      have only discontinued US products and no current FDA label — a random sample of 50;
  (c) «unconfirmed» Open Targets approvals that a current label confirms — a random sample of 50.

Every row carries the URL of the source page to open (DailyMed SPL, EMA product page, or a DailyMed search) and the
sentence the classification relies on.

    py -3 -m evaluation.drug_status_parts      # -> data/eval/n4_part_{a,b,c}.jsonl
"""
import json
import sys
import urllib.parse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
EU_BAD = ("withdrawn", "refused", "revoked", "suspended", "expired", "lapsed")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.drug_status import ema_table, moiety
    g = Graph()
    pl = pd.read_parquet(ROOT / "data" / "ext" / "pair_label.parquet")
    st = pd.read_parquet(ROOT / "data" / "ext" / "drug_status.parquet")
    on_page = pl.disease.map(lambda d: g.origin[d] == "PrimeKG")
    shown = pl[on_page & (pl.source != "Open Targets (unconfirmed)")]
    split = lambda w: (w.split(" | ") + ["", "", ""])[:3]
    out = ROOT / "data" / "eval"
    # (a)
    a = shown[shown.label == "negated"].drop_duplicates(["drug_name", "disease_name"])
    with (out / "n4_part_a.jsonl").open("w", encoding="utf-8") as fh:
        for r in a.itertuples():
            title, url, _ = split(r.where)
            fh.write(json.dumps({"drug": r.drug_name, "disease": r.disease_name, "claim": "the current label says the drug "
                                 "is not effective / not indicated for this disease", "label": title,
                                 "url": url.replace("EMA ", ""), "source": r.source}, ensure_ascii=False) + "\n")
    # (b)
    t = ema_table()
    inn = "International non-proprietary name (INN) / common name"
    eu_url = {}
    for name, status, url in zip(t[inn], t["Medicine status"], t["Medicine URL"]):
        if isinstance(name, str) and str(status).strip().lower() in EU_BAD:
            eu_url.setdefault(moiety(name), url)
    in_site = set(pl[on_page].drug)
    bad = st[st.node.isin(in_site) & ((st.us_status == "discontinued") | st.eu_status.isin(EU_BAD))]
    b = bad.sample(min(50, len(bad)), random_state=3)
    with (out / "n4_part_b.jsonl").open("w", encoding="utf-8") as fh:
        for r in b.itertuples():
            claims = []
            if r.us_status == "discontinued":
                claims.append("US: no marketed product (only discontinued ones) and no current FDA label")
            if r.eu_status in EU_BAD:
                claims.append(f"EU: centrally authorised medicine {r.eu_status}" + (f" ({r.eu_year})" if r.eu_year else ""))
            fh.write(json.dumps({"drug": r.name, "claim": "; ".join(claims),
                                 "url_eu": eu_url.get(moiety(r.name), ""),
                                 "url_us": "https://dailymed.nlm.nih.gov/dailymed/search.cfm?query=" + urllib.parse.quote(r.name)},
                                ensure_ascii=False) + "\n")
    # (c)
    c = pl[on_page & (pl.source == "Open Targets (unconfirmed)") & (pl.label == "confirmed")].drop_duplicates(["drug_name", "disease_name"])
    cs = c.sample(min(50, len(c)), random_state=4)
    with (out / "n4_part_c.jsonl").open("w", encoding="utf-8") as fh:
        for r in cs.itertuples():
            title, url, sent = split(r.where)
            fh.write(json.dumps({"drug": r.drug_name, "disease": r.disease_name, "claim": "a current label lists this "
                                 "disease (or its direct broader term, if marked) as an indication of the drug",
                                 "label": title, "url": url.replace("EMA ", ""), "sentence": sent,
                                 "broader_term": "(broader term)" in r.where}, ensure_ascii=False) + "\n")
    print(f"(a) {len(a)} negated pairs; (b) {len(bad)} drugs, sample {len(b)}; (c) {len(c)} pairs, sample {len(cs)}")
    cases = pl[pl.drug_name.isin(["Griseofulvin", "Ublituximab", "Butoconazole", "Terconazole", "Zinc oxide", "Benzoic acid", "Olaratumab"])]
    print(cases[cases.disease_name.str.contains("candidiasis|multiple sclerosis|systemic sclerosis|sarcoma", case=False)]
          [["drug_name", "disease_name", "source", "label"]].drop_duplicates().to_string(index=False))
    print(st[st.name.isin(["Olaratumab", "Triamcinolone", "Griseofulvin"])][["name", "us_status", "eu_status", "eu_year"]].to_string(index=False))


if __name__ == "__main__":
    main()
