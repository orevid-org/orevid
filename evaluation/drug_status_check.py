"""N4 checks before release (docs/IMPROVEMENT_PLAN.md, N4): two random samples of 50 pairs for a manual check against the
label texts — 50 «approved» pairs of the website (any category) and 50 pairs whose category changes (an Open Targets
«unconfirmed» approval confirmed by a current label, or a PrimeKG / confirmed Open Targets approval flagged). For each
pair the evidence is written out: the label sentence that names the disease, or the start of the indication texts.

    py -3 -m evaluation.drug_status_check sample    # -> data/eval/drug_status_sample.jsonl (to check by hand)
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent


def evidence(row, fda, ema, keys, dn, parents):
    from sources.drug_status import SENT
    names = dn.get(row.disease, {row.disease_name.lower()})
    broader = set().union(*(dn.get(p, set()) for p in parents.get(row.disease, ())))
    labels = [l for k in keys for l in fda.get(k, [])] + [l for k in keys for l in ema.get(k, [])]
    out = []
    for text, combo, title, src in labels:
        for sent in SENT.split(text):
            low = sent.lower()
            if any(re.search(r"(?<![a-z])" + re.escape(s) + r"(?![a-z])", low) for s in names | broader):
                out.append(f"[{title}{' — combination' if combo else ''}] {sent.strip()[:400]}")
                break
        if len(out) >= 3:
            break
    if not out:
        out = [f"[{t}{' — combination' if c else ''}] {x[:500]}" for x, c, t, s in labels[:2]]
    return out


def sample():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.drug_status import disease_names, ema_labels, ema_table, fda_labels, moiety
    from sources.hierarchy import load as load_hierarchy
    from sources.practical import synonyms
    g = Graph()
    syn = synonyms()
    pl = pd.read_parquet(ROOT / "data" / "ext" / "pair_label.parquet")
    pl = pl[pl.disease.map(lambda d: g.origin[d] == "PrimeKG")]  # pairs that can appear on a page
    fda, ema = fda_labels(), ema_labels(ema_table())
    dn = disease_names(g)
    h = load_hierarchy()
    parents = {}
    for a, b in zip(h.parent.astype(int), h.child.astype(int)):
        parents.setdefault(b, []).append(a)
    shown = pl[pl.source != "Open Targets (unconfirmed)"]
    changed = pd.concat([pl[(pl.source == "Open Targets (unconfirmed)") & (pl.label == "confirmed")],
                         shown[shown.label != "confirmed"]])
    s1 = shown.sample(50, random_state=1).assign(sample="approved pairs")
    s2 = changed.sample(50, random_state=2).assign(sample="pairs that change category")
    out = ROOT / "data" / "eval" / "drug_status_sample.jsonl"
    with out.open("w", encoding="utf-8") as fh:
        for r in pd.concat([s1, s2]).itertuples():
            keys = ({moiety(g.name(r.drug))} | {moiety(x) for x in syn.get(r.drug, ())}) - {""}
            fh.write(json.dumps({"sample": r.sample, "drug": r.drug_name, "disease": r.disease_name, "source": r.source,
                                 "label": r.label, "where": r.where,
                                 "evidence": evidence(r, fda, ema, keys, dn, parents)}, ensure_ascii=False) + "\n")
    print(f"-> {out}")
    print(pd.concat([s1, s2]).groupby(["sample", "label"]).size())


if __name__ == "__main__":
    {"sample": sample}[sys.argv[1]]()
