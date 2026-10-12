"""Progress table of the comparison with TxGNN: for each split that TxGNN has already finished,
runs our best model on it (extra sources + masked training) and writes docs/TXGNN_PROGRESS.md
with TxGNN, popularity and our model, and the state of the queue.

    py -3 -m evaluation.txgnn_progress            # one pass
    py -3 -m evaluation.txgnn_progress --watch    # every 5 minutes until the queue finishes
"""
import argparse
import re
import sys
import time
from datetime import datetime, timedelta

import pandas as pd

from evaluation.txgnn import ROOT, CACHE, AREAS, NO_NEAR_DUPLICATES, load_kg, node_types, evaluate, published

QUEUE = ["complex_disease 1", "mental_health 1", "diabetes 1", "neurodigenerative 1", "cardiovascular 1",
         "anemia 1", "autoimmune 1", "adrenal_gland 1", "metabolic_disorder 1", "cell_proliferation 1",
         "complex_disease 2", "complex_disease 3", "complex_disease 4", "complex_disease 5"]
NAMES = {"complex_disease": "Zero-shot", "cell_proliferation": "Cancer", "mental_health": "Mental health",
         "cardiovascular": "Cardiovascular", "anemia": "Anaemia", "adrenal_gland": "Adrenal",
         "autoimmune": "Autoimmune", "metabolic_disorder": "Metabolic", "diabetes": "Diabetes",
         "neurodigenerative": "Neurodegenerative"}
DEST = ROOT / "docs" / "TXGNN_PROGRESS.md"
PROGRESS = CACHE / "progress"
OURS = "gbm + masked training"
TX = "TxGNN (run here)"
HARD = NO_NEAR_DUPLICATES


def queue_state():
    """{«split seed»: (state, start, end)} from data/eval/txgnn_queue.txt (the output of
    evaluation/train_txgnn.sh, through evaluation/txgnn_queue.ps1)."""
    f = CACHE / "txgnn_queue.txt"
    res = {}
    if not f.exists():
        return res
    for line in f.read_text("utf-8", errors="ignore").splitlines():
        m = re.match(r"﻿?(\w+) (\d+): (starts|exit \d+|already done) ?(\S*)", line.strip())
        if not m:
            continue
        key, what, hour = f"{m.group(1)} {m.group(2)}", m.group(3), m.group(4)
        if what == "starts":
            res[key] = ("running", hour, None)
        elif what.startswith("exit"):
            res[key] = ("finished" if what == "exit 0" else "error", res.get(key, (None, None))[1], hour)
    return res


def results(kg, types, key):
    """Metrics of a split (with TxGNN), saved in data/eval/progress/."""
    split, seed = key.split()
    f = PROGRESS / f"{split}_{seed}.csv"
    if f.exists():
        return pd.read_csv(f, index_col=0)
    if not (CACHE / "txgnn_predictions" / f"TxGNN_{seed}_{split}.parquet").exists():
        return None
    r, c = evaluate(kg, types, split, int(seed), 1500, "extra", 5)
    d = pd.DataFrame(r).T
    d["diseases"], d["hard"] = c["_diseases"], c["_hard"]
    PROGRESS.mkdir(parents=True, exist_ok=True)
    d.to_csv(f)
    return d


def cell(d, row, m):
    return f"{d.loc[row, m]:.2f}" if row in d.index else "–"


def duel(d, m, pre=""):
    """«TxGNN → ours», with the winner in bold."""
    a, b = pre + TX, pre + OURS
    if a not in d.index or b not in d.index:
        return f"{cell(d, a, m)} → {cell(d, b, m)}"
    x, y = d.loc[a, m], d.loc[b, m]
    return f"{'**' if x > y else ''}{cell(d, a, m)}{'**' if x > y else ''} → "            f"{'**' if y > x else ''}{cell(d, b, m)}{'**' if y > x else ''}"


def write(kg, types):
    queue = queue_state()
    done = {c: results(kg, types, c) for c in QUEUE}
    done = {c: d for c, d in done.items() if d is not None}
    L = [f"# Progress: TxGNN against our model\n",
         f"Updated: {datetime.now():%d/%m %H:%M}. Regenerated automatically (`evaluation/txgnn_progress.py`).\n",
         "Same split and same metrics code for both (see `docs/TXGNN_EVALUATION.md`). "
         "«Ours» = gradient boosting with extra sources and masked training. "
         "«Hard» = test diseases without any near-identical sibling with known drugs.\n",
         "## Recall@100: out of every 100 good drugs, how many come out in the top 100\n",
         "| Split | Diseases (hard) | Popularity | TxGNN → ours | Hard: TxGNN → ours | TxGNN published |",
         "|---|---|---|---|---|---|"]
    for c, d in done.items():
        split, seed = c.split()
        pub = published(split)
        p = f"{pub.loc[('TxGNN', 'Recall@100'), 'mean']:.2f}" if pub is not None else "–"
        L.append(f"| {NAMES[split]} (seed {seed}) | {int(d['diseases'].iloc[0])} ({int(d['hard'].iloc[0])}) | "
                 f"{cell(d, 'popularity', 'Recall@100')} | {duel(d, 'Recall@100')} | {duel(d, 'Recall@100', HARD)} | {p} |")
    L += ["\n## Other metrics (TxGNN → ours; in bold, the best)\n",
          "AUPRC and AUROC: positives against the same number of random negatives among drugs with some indication. "
          "Precision@K: hits among the top K (K = good drugs). Avg Rank: 1 = at the top.\n",
          "| Split | AUPRC | AUROC | Precision@K | Avg Rank | AUPRC hard |", "|---|---|---|---|---|---|"]
    for c, d in done.items():
        split, seed = c.split()
        L.append(f"| {NAMES[split]} (seed {seed}) | {duel(d, 'AUPRC')} | {duel(d, 'AUROC')} | "
                 f"{duel(d, 'Precision@K')} | {duel(d, 'Avg Rank')} | {duel(d, 'AUPRC', HARD)} |")
    # queue state
    durations = []
    for c, (e, start, end) in queue.items():
        if e == "finished" and start and end:
            t = (datetime.strptime(end, "%H:%M:%S") - datetime.strptime(start, "%H:%M:%S")).seconds
            durations.append(t)
    missing = [c for c in QUEUE if queue.get(c, ("pending",))[0] != "finished"]
    L += ["\n## TxGNN queue\n", "| Split | State |", "|---|---|"]
    for c in QUEUE:
        e, start, end = queue.get(c, ("pending", None, None))
        e = {"running": f"⏳ running (since {start})", "finished": f"✅ finished ({end})", "error": "❌ error (will retry)",
             "pending": "pending"}.get(e, e)
        L.append(f"| {NAMES[c.split()[0]]} (seed {c.split()[1]}) | {e} |")
    if durations and missing:
        mean = sum(durations) / len(durations)
        L.append(f"\nEach split takes ~{mean / 60:.0f} min. {len(missing)} left: "
                 f"it would finish around {datetime.now() + timedelta(seconds=mean * len(missing)):%H:%M on %d/%m}.")
    DEST.write_text("\n".join(L) + "\n", "utf-8")
    print(f"{datetime.now():%H:%M} table with {len(done)} splits -> {DEST}", flush=True)
    return len(missing) == 0


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    p = argparse.ArgumentParser()
    p.add_argument("--watch", action="store_true")
    a = p.parse_args()
    kg = load_kg()
    types = node_types(kg)
    while True:
        finished = write(kg, types)
        if not a.watch or finished:
            break
        time.sleep(300)


if __name__ == "__main__":
    main()
