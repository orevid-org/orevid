"""Compounds that bind a protein, according to ChEMBL (EBI public API; data CC BY-SA 3.0).

For a protein (UniProt): how many distinct compounds have potent measured activity on it
(pChEMBL ≥ 6, i.e. ≤ 1 µM; and ≥ 7, ≤ 100 nM), the most potent one and those that are already drugs.
They are chemical starting points even if they are not drugs (possible new derivatives). It is queried on demand (only the targets that are shown) and saved to
data/cache/chembl/ so queries are not repeated.

    py -3 -m sources.chembl_api P11473
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

API = "https://www.ebi.ac.uk/chembl/api/data"
CACHE = Path(__file__).parent.parent / "data" / "cache" / "chembl"
PAGES = 5  # at most 5,000 activities per target (if there are more, the counts are "at least")


def _get(path, **params):
    url = f"{API}/{path}.json?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "orevid/0.1", "Accept": "application/json"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                return json.loads(r.read().decode())
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2 * (attempt + 1))


def compounds(uniprot):
    """Dictionary with the summary (or None if ChEMBL does not have the protein)."""
    if not uniprot:
        return None
    f = CACHE / f"{uniprot}.json"
    if f.exists():
        cached = json.loads(f.read_text("utf-8"))
        # Entries cached before the English key names lack "complete": they are fetched again.
        if cached is None or "complete" in cached:
            return cached
    t = _get("target", target_components__accession=uniprot, target_type="SINGLE PROTEIN",
             organism="Homo sapiens", limit=1)
    res = None
    if t["targets"]:
        tid = t["targets"][0]["target_chembl_id"]
        best = {}
        offset, complete = 0, True
        for _ in range(PAGES):
            a = _get("activity", target_chembl_id=tid, pchembl_value__gte=6, limit=1000, offset=offset,
                     only="molecule_chembl_id,pchembl_value,molecule_pref_name")
            for x in a["activities"]:
                m, v = x["molecule_chembl_id"], float(x["pchembl_value"])
                if m not in best or v > best[m][0]:
                    best[m] = (v, x.get("molecule_pref_name"))
            if not a["page_meta"]["next"]:
                break
            offset += 1000
        else:
            complete = False
        top = sorted(best.items(), key=lambda kv: -kv[1][0])
        res = {"target_chembl_id": tid, "complete": complete,
               "n_1uM": len(best), "n_100nM": sum(v >= 7 for v, _ in best.values()),
               "best": [{"id": m, "pchembl": v, "name": n} for m, (v, n) in top[:5]],
               "named": [{"id": m, "pchembl": v, "name": n} for m, (v, n) in top if n][:5]}
    CACHE.mkdir(parents=True, exist_ok=True)
    f.write_text(json.dumps(res, ensure_ascii=False), "utf-8")
    return res


def summary(r):
    """One line for the output of the "targets" command."""
    if not r:
        return ""
    more = "" if r["complete"] else "≥"
    line = f"{more}{r['n_1uM']} compounds ≤ 1 µM ({more}{r['n_100nM']} ≤ 100 nM; ChEMBL {r['target_chembl_id']})"
    if r["best"]:
        m = r["best"][0]
        line += f"; most potent {m['name'] or m['id']} (pChEMBL {m['pchembl']:.1f})"
    return line


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    r = compounds(sys.argv[1])
    print(summary(r))
    print(json.dumps(r, ensure_ascii=False, indent=1)[:1500])
