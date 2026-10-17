"""OREVID + RotatE ensemble (docs/IMPROVEMENT_PLAN.md, I5: «an OREVID + embedding ensemble follows the P33 rule»).

Score per disease = z(OREVID) + w · z(RotatE), with z the standardisation over the task's drugs for that disease
(benchmark prediction files; nothing is retrained). w ∈ {0.1, 0.25, 0.5, 1.0}, fixed before running, chosen on
zero-shot seeds 1-2 by mean Recall@100. P33 rule (docs/XGBOOST.md): adopted only if the chosen ensemble beats
OREVID on mean Recall@100 in all three confirmation groups (seeds 3-5, nine areas, external exam), by at least
0.01 on average over the three, and loses no more than 0.01 on any hard or non-obvious subset.

    py -3 -m evaluation.kge_ensemble      # -> docs/KGE.md section, data/eval/kge/ensemble.json
"""
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BENCH = ROOT / "benchmark"
WEIGHTS = (0.1, 0.25, 0.5, 1.0)
AREAS = ["cell_proliferation", "mental_health", "cardiovascular", "anemia", "adrenal_gland", "autoimmune",
         "metabolic_disorder", "diabetes", "neurodigenerative"]
GROUPS = {"selection (seeds 1-2)": [f"zero_shot_{s}" for s in (1, 2)],
          "seeds 3-5": [f"zero_shot_{s}" for s in (3, 4, 5)],
          "areas": [f"area_{a}" for a in AREAS], "external": ["external_ot2606"]}


def z(S):
    return (S - S.mean(axis=0)) / (S.std(axis=0) + 1e-9)


def evaluate(task_name, w):
    sys.path.insert(0, str(BENCH))
    from bench import evaluate_task, load_predictions
    task = json.loads((BENCH / "tasks" / f"{task_name}.json").read_text())
    O, _ = load_predictions(BENCH / "predictions" / "orevid_lambdarank" / f"{task_name}.parquet", task)
    R, _ = load_predictions(BENCH / "predictions" / "kge_rotate" / f"{task_name}.parquet", task)
    S = z(O) if w == 0 else z(O) + w * z(R)
    r = evaluate_task(task, S)
    subs = {"all": r["all"]["Recall@100"]}
    for k in ("hard", "non_obvious", "recent"):
        if k in r:
            subs[k] = r[k]["Recall@100"]
    return subs


def group_mean(names, w):
    rs = [evaluate(n, w) for n in names]
    return {k: float(np.mean([r[k] for r in rs if k in r])) for k in rs[0]}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    res = {g: {str(w): group_mean(n, w) for w in (0.0, *WEIGHTS)} for g, n in GROUPS.items()}
    sel = res["selection (seeds 1-2)"]
    w_best = max(WEIGHTS, key=lambda w: sel[str(w)]["all"])
    conf = ("seeds 3-5", "areas", "external")
    gains = {g: res[g][str(w_best)]["all"] - res[g]["0.0"]["all"] for g in conf}
    losses = {f"{g} {k}": res[g][str(w_best)][k] - res[g]["0.0"][k] for g in conf for k in ("hard", "non_obvious")
              if k in res[g]["0.0"]}
    adopted = bool(all(v > 0 for v in gains.values()) and np.mean(list(gains.values())) >= 0.01
                   and all(v >= -0.01 for v in losses.values()))
    out = {"results": res, "w": w_best, "gains": gains, "subset_changes": losses, "adopted": adopted}
    (ROOT / "data" / "eval" / "kge" / "ensemble.json").write_text(json.dumps(out, indent=1), "utf-8")
    L = ["| Group | OREVID | " + " | ".join(f"+ {w} · RotatE" for w in WEIGHTS) + " |", "|---|" + "---|" * (len(WEIGHTS) + 1)]
    for g in GROUPS:
        for k in res[g]["0.0"]:
            L.append(f"| {g}, {k} | {res[g]['0.0'][k]:.3f} | " + " | ".join(f"{res[g][str(w)][k]:.3f}" for w in WEIGHTS) + " |")
    L += ["", f"Chosen on seeds 1-2: w = {w_best}. Gains in the confirmation groups: "
              + ", ".join(f"{g} {v:+.3f}" for g, v in gains.items()) + f" (mean {np.mean(list(gains.values())):+.3f}); "
              + "hard / non-obvious changes: " + ", ".join(f"{k} {v:+.3f}" for k, v in losses.items()) + ".",
          "", f"**Decision (P33 rule):** {'adopted' if adopted else 'not adopted; production unchanged'}."]
    text = "\n".join(L)
    print(text)
    return text


if __name__ == "__main__":
    main()
