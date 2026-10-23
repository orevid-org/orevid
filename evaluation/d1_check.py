"""D1 release check (docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026: D1», decision rule): 50 random pairs flagged by
the signal that passed (S1, the P15 rule) as the website would show them — website diseases, approved drugs not
indicated and not already contraindicated (PrimeKG) for the disease — each with the current US label sections
(contraindications, warnings and precautions, adverse reactions) from openFDA and the DailyMed URL to open.

A pair is right if a section of the current label names the disease, a direct broader term, or a condition the
disease consists of as a reason to avoid or monitor the drug, or as an adverse effect; wrong otherwise.

    py -3 -m evaluation.d1_check          # -> data/eval/d1_check.jsonl
"""
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
N, POOL, SEED = 50, 400, 23
SECTIONS = ("contraindications", "boxed_warning", "warnings_and_cautions", "warnings", "precautions", "adverse_reactions")


def label(name):
    q = urllib.parse.quote(f'openfda.generic_name:"{name}"')
    try:
        with urllib.request.urlopen(f"https://api.fda.gov/drug/label.json?search={q}&limit=1", timeout=30) as r:
            res = json.load(r)["results"][0]
    except Exception:
        return None
    set_id = res.get("set_id", "")
    return {"title": (res.get("openfda", {}).get("brand_name") or [""])[0],
            "url": f"https://dailymed.nlm.nih.gov/dailymed/drugInfo.cfm?setid={set_id}" if set_id else "",
            **{s: " ".join(res.get(s, []))[:20000] for s in SECTIONS}}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.direction import Direction, checkpoint_rule, risk_flags
    from evaluation.twins import twins
    from predictor import drug_status
    D = Direction()
    g = D.g
    status = drug_status(g)
    approved = np.array(sorted(f for f in g.drugs if status.get(int(f)) == "approved"), dtype=int)
    blocked, autoimmune = checkpoint_rule(g)
    tw = twins(g, "disease")
    C, A = D.contra.tocsr(), D.approved.tocsr()
    rng = np.random.default_rng(SEED)
    site = [int(d) for d in np.flatnonzero(g.types == "disease") if g.origin[d] == "PrimeKG"]
    pool = []
    for d in rng.choice(site, POOL, replace=False):
        d = int(d)
        same = tw.get(d, {d})
        drop = {int(x) for t in same for M in (C, A) for x in M.indices[M.indptr[t]:M.indptr[t + 1]]}
        cand = np.array([f for f in approved if int(f) not in drop], dtype=int)
        flag, nc, ni, pc = risk_flags(D, d, cand, blocked, autoimmune)
        pool += [(d, int(f), int(c), int(i), bool(p)) for f, c, i, p in zip(cand[flag], nc[flag], ni[flag], pc[flag])]
    print(f"{len(pool):,} flags in {POOL} random website diseases ({len(pool) / POOL:.0f} per disease)")
    pick = [pool[i] for i in rng.choice(len(pool), N, replace=False)]
    with (ROOT / "data" / "eval" / "d1_check.jsonl").open("w", encoding="utf-8") as fh:
        for d, f, c, i, p in pick:
            lab = label(g.name(f))
            time.sleep(0.3)
            fh.write(json.dumps({"disease": g.name(d), "drug": g.name(f), "contraindicated_related": c, "indicated_related": i,
                                 "checkpoint": p, "label": lab}, ensure_ascii=False) + "\n")
            print(f"{g.name(d)} | {g.name(f)} | {c}>{i}{' checkpoint' if p else ''} | {'label' if lab else 'NO LABEL'}")


if __name__ == "__main__":
    main()
