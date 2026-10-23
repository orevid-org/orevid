"""OnSIDES adverse events from current drug labels (D1, docs/IMPROVEMENT_PLAN.md «Registered 23-10-2026: D1»).

OnSIDES v3.1.1 data (tatonetti-lab; v3.2.1 is a code-only release with the same data; data CC BY 4.0, cite
doi:10.1016/j.medj.2025.100642). Mapping fixed in the registration:
- ingredients: RxNorm ingredient names of single-ingredient products -> graph drugs by active moiety (as N4);
- adverse events: MedDRA preferred / lowest-level terms -> graph diseases. MONDO has only 3 MedDRA cross-references
  and OnSIDES has no UMLS ids, so, as registered, by exact normalised name of the disease, its twin nodes or its Open
  Targets exact synonyms.

    py -3 -m sources.onsides        # -> data/ext/onsides.parquet (drug node, disease node, term, sections)
"""
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "onsides" / "onsides-v3.1.1.zip"
DEST = ROOT / "data" / "ext" / "onsides.parquet"


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.drug_status import disease_names, moiety
    from sources.practical import synonyms
    z = zipfile.ZipFile(RAW)
    rd = lambda f, **k: pd.read_csv(z.open(f"csv/{f}"), **k)
    ing = rd("vocab_rxnorm_ingredient.csv", dtype=str)
    i2p = rd("vocab_rxnorm_ingredient_to_product.csv", dtype=str)
    l2p = rd("product_to_rxnorm.csv", dtype=str)
    med = rd("vocab_meddra_adverse_effect.csv", dtype=str)
    ae = rd("product_adverse_effect.csv", usecols=["product_label_id", "label_section", "effect_meddra_id"], dtype=str)
    # single-ingredient products only (a combination's adverse events are not given to each ingredient)
    n_ing = i2p.groupby("product_id").ingredient_id.nunique()
    single = i2p[i2p.product_id.isin(n_ing[n_ing == 1].index)]
    prod_ing = dict(zip(single.product_id, single.ingredient_id))
    label_ing = defaultdict(set)
    for lab, prod in zip(l2p.label_id, l2p.rxnorm_product_id):
        if prod in prod_ing:
            label_ing[lab].add(prod_ing[prod])
    ing_name = dict(zip(ing.rxnorm_id, ing.rxnorm_name))
    term = dict(zip(med.meddra_id, med.meddra_name.str.lower().str.strip()))
    # graph drugs by moiety
    g = Graph()
    syn = synonyms()
    by_key = defaultdict(set)
    for f in g.drugs:
        f = int(f)
        for k in ({moiety(g.name(f))} | {moiety(s) for s in syn.get(f, ())}) - {""}:
            by_key[k].add(f)
    ing_drugs = {i: by_key.get(moiety(n), set()) for i, n in ing_name.items()}
    # graph diseases by exact name
    dn = disease_names(g)
    name_dis = defaultdict(set)
    for d, names in dn.items():
        for n in names:
            name_dis[n].add(d)
    rows = defaultdict(set)
    n_pairs = 0
    for lab, sec, mid in zip(ae.product_label_id, ae.label_section, ae.effect_meddra_id):
        n_pairs += 1
        t = term.get(mid)
        if not t or lab not in label_ing:
            continue
        for d in name_dis.get(t, ()):
            for i in label_ing[lab]:
                for f in ing_drugs.get(i, ()):
                    rows[(f, d, t)].add(str(sec))
    out = pd.DataFrame([{"drug": f, "disease": d, "term": t, "sections": "|".join(sorted(s))} for (f, d, t), s in rows.items()])
    DEST.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(DEST, index=False)
    mapped_ing = sum(bool(v) for v in ing_drugs.values())
    terms_hit = {t for t in term.values() if t in name_dis}
    print(f"OnSIDES: {n_pairs:,} label–adverse-event rows, {len(label_ing):,} single-ingredient labels, "
          f"{mapped_ing:,} of {len(ing_name):,} ingredients mapped to graph drugs, {len(terms_hit):,} of {len(set(term.values())):,} "
          f"MedDRA terms matched to graph diseases; {len(out):,} drug–disease pairs ({out.drug.nunique():,} drugs, "
          f"{out.disease.nunique():,} diseases)")
    return out


if __name__ == "__main__":
    build()
