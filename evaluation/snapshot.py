"""Dated snapshot of the predictions, for a prospective validation (see the analysis plan it writes).

For every disease on the website it stores the model's top 100 drugs (in model order, with the
specificity score) and what the website shows for each one: direct extension, candidate (shared
genes), exploratory (no indication), or not shown. Hidden-but-ranked drugs are the control group.
It also copies the exact website pages (web/public/data) so that what was shown can be checked later.

    python -m evaluation.snapshot            # -> data/snapshot/<date>/  (zip it and upload to Zenodo)
"""
import datetime
import json
import shutil
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from graph import GENES_ALL, Graph
from predictor import Predictor

ROOT = Path(__file__).resolve().parent.parent
TOP = 100


def shown_sections():
    """{(disease id, drug name): section} from the published website data."""
    out = {}
    for f in sorted((ROOT / "web" / "public" / "data" / "d").glob("*.json")):
        for did, page in json.loads(f.read_text("utf-8")).items():
            for section, key in (("extension", "extensions"), ("candidate", "candidates"), ("exploratory", "no_use")):
                for c in page.get(key, []):
                    for name in [c["name"]] + c.get("also", []):
                        out[(int(did), name)] = section
    return out


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    today = datetime.date.today().isoformat()
    dest = ROOT / "data" / "snapshot" / today
    dest.mkdir(parents=True, exist_ok=True)
    shown = shown_sections()
    P = Predictor()
    g = P.g
    pkg = np.flatnonzero((g.types == "disease") & (g.origin == "PrimeKG"))
    with_data = pkg[(np.diff(g.M(GENES_ALL).indptr)[pkg] > 0)
                    | (np.diff(g.M("disease_phenotype_positive").indptr)[pkg] > 0)]
    ext_id = dict(zip(g.nodes.node_index, g.nodes.node_id.astype(str)))
    rows, t0 = [], time.time()
    for k, d in enumerate(with_data):
        drugs, score, _ = P.predict(int(d))
        order = np.argsort(-score)[:TOP]
        for rank, j in enumerate(order, 1):
            f = int(drugs[j])
            name = g.name(f)
            rows.append((int(d), g.name(int(d)), ext_id.get(int(d), ""), rank, f, name, ext_id.get(f, ""),
                         round(float(score[j]), 4), shown.get((int(d), name), "not shown")))
        if (k + 1) % 500 == 0:
            print(f"  {k + 1}/{len(with_data)} diseases ({time.time() - t0:.0f} s)", flush=True)
    df = pd.DataFrame(rows, columns=["disease_node", "disease", "disease_id", "rank", "drug_node", "drug", "drug_id",
                                     "score", "website_section"])
    df.to_csv(dest / "ranking_top100.csv.gz", index=False, compression="gzip")
    # shown pairs that fall outside the top 100 (extensions and exploratory lists come from deeper ranks)
    in_top = set(zip(df.disease_node, df.drug))
    extra = [(d, n, s) for (d, n), s in shown.items() if (d, n) not in in_top]
    pd.DataFrame(extra, columns=["disease_node", "drug", "website_section"]).to_csv(
        dest / "shown_outside_top100.csv", index=False)
    shutil.make_archive(str(dest / "website_data"), "zip", ROOT / "web" / "public" / "data")
    counts = df.website_section.value_counts().to_dict()
    (dest / "summary.json").write_text(json.dumps({"date": today, "diseases": int(len(with_data)), "top": TOP,
                                                   "sections_in_top100": counts, "shown_outside_top100": len(extra)},
                                                  indent=1), "utf-8")
    print(f"snapshot {today}: {len(with_data)} diseases, {len(df):,} rows, {counts} -> {dest}")


if __name__ == "__main__":
    main()
