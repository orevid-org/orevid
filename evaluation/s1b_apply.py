"""S1b (P92): rewrite the status shown on the built pages for drugs approved in some country by an agency
(DrugCentral) or by ChEMBL and the WHO ATC together, without rebuilding the pages. A full rebuild gives the same result
(`web/build.py` computes `status_shown` the same way).

    py -3 -m evaluation.s1b_apply
"""
import collections
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RELABEL = ("experimental", "investigational", "withdrawn", "unknown")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from predictor import drug_status
    from evaluation.approval import approved_somewhere
    G = Graph()
    status = drug_status(G)
    ok = {f for f, s in approved_somewhere(G).items() if "DrugCentral" in s or {"ChEMBL", "ATC"} <= set(s)}
    changed = collections.Counter()
    for fj in sorted((ROOT / "web" / "public" / "data" / "d").glob("*.json")):
        pages = json.loads(fj.read_text("utf-8"))
        dirty = False
        for v in pages.values():
            for sec in ("extensions", "candidates", "no_use"):
                for c in v.get(sec) or []:
                    f = c.get("node")
                    if f is not None and c.get("status") in RELABEL and status.get(int(f)) in RELABEL and int(f) in ok:
                        changed[(c["status"], c["name"])] += 1
                        c["status"] = "approved in some country"
                        dirty = True
        if dirty:
            fj.write_text(json.dumps(pages, ensure_ascii=False, separators=(",", ":"), allow_nan=False), "utf-8")
    by = collections.Counter()
    for (s, _), n in changed.items():
        by[s] += n
    print(f"entries relabelled: {sum(changed.values()):,} ({dict(by)}); drugs: {len({n for _, n in changed}):,}")
    print("most frequent:", [(n, k) for (s, n), k in changed.most_common(12)])


if __name__ == "__main__":
    main()
