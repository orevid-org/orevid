"""Knowledge-graph embedding baselines on OREVID-bench (docs/IMPROVEMENT_PLAN.md, I5; protocol fixed there).

TransE, ComplEx and RotatE (PyKEEN) trained on each task's PrimeKG training edges (the directed edge list TxGNN
uses, `benchmark/bench.py` `load_kg` / `train_edges`), scoring (drug, indication, disease) for every drug of the
task and every test disease. Grid fixed in the plan: embedding dimension 128/256, learning rate 1e-3/1e-2,
100 epochs, negative sampler «basic», other settings PyKEEN's defaults (Adam, one negative per positive, each
model's default loss); chosen on zero_shot_1 only (Recall@100, all diseases), then run unchanged on seeds 2-5,
the nine areas and the external exam (trained on zero_shot_1, as every model). Batch size 4,096 (not in the plan;
PyKEEN's pipeline default of 256 would need ~16,000 steps per epoch on 4 M edges).

Steps (prepare and select on Windows; train in WSL with ~/pykeen-env: torch + CUDA, pykeen):

    py -3 -m evaluation.kge_bench prepare                 # -> data/eval/kge/kg.npz (needs data/raw/txgnn/kg.csv)
    bash evaluation/kge_bench.sh grid                     # WSL: 3 models x 4 settings on zero_shot_1
    py -3 -m evaluation.kge_bench select                  # -> data/eval/kge/selected.json, docs/KGE.md (grid)
    bash evaluation/kge_bench.sh final                    # WSL: chosen settings on every task
    py -3 benchmark/bench.py evaluate --pred benchmark/predictions/kge_<model>
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
TASKS = ROOT / "benchmark" / "tasks"
OUT = ROOT / "data" / "eval" / "kge"
MODELS = ["TransE", "ComplEx", "RotatE"]
GRID = [(d, lr) for d in (128, 256) for lr in (1e-3, 1e-2)]
EPOCHS, BATCH, SEED = 100, 4096, 0


def prepare():
    import pandas as pd
    # bench.load_kg's selection, reading only the columns it uses (the whole file needs ~4 GB of memory)
    kg = pd.read_csv(ROOT / "data" / "raw" / "txgnn" / "kg.csv", usecols=["relation", "x_index", "x_type", "y_index", "y_type"],
                     dtype={"relation": "category", "x_type": "category", "y_type": "category",
                            "x_index": np.int32, "y_index": np.int32})
    keep = []
    for rel in kg.relation.cat.categories:
        d = kg[kg.relation == rel]
        if "_" in rel and rel.split("_")[0] == rel.split("_")[1]:
            pair = np.minimum(d.x_index, d.y_index).astype(np.int64) * 10**6 + np.maximum(d.x_index, d.y_index)
            keep.append(d.index[~pair.duplicated().to_numpy()])
        else:
            keep.append(d.index[(d.x_type == d.x_type.iloc[0]).to_numpy()])
    kg = kg.loc[np.sort(np.concatenate(keep))].reset_index(drop=True)
    kg["relation"] = kg.relation.astype(str)
    s = np.load(TASKS / "zero_shot_1.split.npz")
    assert max(s["valid"].max(), s["test"].max()) < len(kg)
    rels = sorted(kg.relation.unique())
    r = pd.Categorical(kg.relation, categories=rels).codes.astype(np.int16)
    ind = kg[kg.relation == "indication"]
    drug_first = bool((ind.x_type == "drug").all())
    assert drug_first or bool((ind.y_type == "drug").all())
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(OUT / "kg.npz", h=kg.x_index.to_numpy(np.int32), r=r, t=kg.y_index.to_numpy(np.int32),
                        relations=np.array(rels), n_nodes=int(max(kg.x_index.max(), kg.y_index.max()) + 1),
                        drug_first=drug_first)
    print(f"{len(kg):,} directed edges, {len(rels)} relations; indication stored as "
          f"{'drug -> disease' if drug_first else 'disease -> drug'} -> {OUT / 'kg.npz'}")


def tasks_trained_on(split):
    """Tasks whose test diseases are scored with the model trained on `split`."""
    out = []
    for f in sorted(TASKS.glob("*.json")):
        t = json.loads(f.read_text())
        if t["name"] == split or t.get("train_split") == split:
            out.append(t)
    return out


def train(split, model_name, dim, lr, dest, epochs=EPOCHS):
    """Runs in WSL. Trains one model on the training edges of `split` and writes <dest>/<task>.parquet for
    every task that uses that training split."""
    import pandas as pd
    import torch
    from pykeen import models
    from pykeen.training import SLCWATrainingLoop
    from pykeen.triples import CoreTriplesFactory
    z = np.load(OUT / "kg.npz")
    s = np.load(TASKS / f"{split}.split.npz")
    mask = np.ones(len(z["h"]), bool)
    mask[s["valid"]] = False
    mask[s["test"]] = False
    triples = torch.as_tensor(np.stack([z["h"][mask], z["r"][mask].astype(np.int32), z["t"][mask]], 1), dtype=torch.long)
    rels = list(z["relations"])
    tf = CoreTriplesFactory(mapped_triples=triples, num_entities=int(z["n_nodes"]), num_relations=len(rels))
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = getattr(models, model_name)(triples_factory=tf, embedding_dim=dim, random_seed=SEED).to(dev)
    loop = SLCWATrainingLoop(model=model, triples_factory=tf, optimizer="Adam", optimizer_kwargs={"lr": lr},
                             negative_sampler="basic")
    t0 = time.time()
    losses = loop.train(triples_factory=tf, num_epochs=epochs, batch_size=BATCH, use_tqdm=False)
    print(f"{split} {model_name} d={dim} lr={lr}: {len(triples):,} triples, {time.time() - t0:.0f} s, "
          f"loss {losses[0]:.4f} -> {losses[-1]:.4f}", flush=True)
    model.eval()
    ind = rels.index("indication")
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    for task in tasks_trained_on(split):
        drugs = torch.as_tensor(task["drugs"], dtype=torch.long)
        rows = []
        with torch.no_grad():
            for r in task["diseases"]:
                d = torch.full_like(drugs, int(r["disease"]))
                h, t = (drugs, d) if bool(z["drug_first"]) else (d, drugs)
                hrt = torch.stack([h, torch.full_like(h, ind), t], 1).to(dev)
                sc = model.score_hrt(hrt).squeeze(-1).cpu().numpy()
                rows.append(pd.DataFrame({"disease": int(r["disease"]), "drug": task["drugs"], "score": sc}))
        pd.concat(rows).to_parquet(dest / f"{task['name']}.parquet", index=False)
        print(f"  -> {dest / (task['name'] + '.parquet')}", flush=True)


def grid_name(m, dim, lr):
    return f"{m}_d{dim}_lr{lr:g}"


def select():
    """Grid on zero_shot_1: best setting per model by Recall@100 (all test diseases)."""
    sys.path.insert(0, str(ROOT / "benchmark"))
    from bench import evaluate_task, load_predictions
    task = json.loads((TASKS / "zero_shot_1.json").read_text())
    rows, best = [], {}
    for m in MODELS:
        for dim, lr in GRID:
            f = OUT / "grid" / grid_name(m, dim, lr) / "zero_shot_1.parquet"
            if not f.exists():
                continue
            S, _ = load_predictions(f, task)
            r = evaluate_task(task, S)["all"]["Recall@100"]
            rows.append((m, dim, lr, r))
            if r > best.get(m, (None, None, -1))[2]:
                best[m] = (dim, lr, r)
    (OUT / "selected.json").write_text(json.dumps({m: {"dim": d, "lr": lr, "recall100_zero_shot_1": r}
                                                   for m, (d, lr, r) in best.items()}, indent=1))
    for m, dim, lr, r in rows:
        print(f"{m:8s} dim {dim:3d} lr {lr:g}: Recall@100 {r:.3f}{'  <- chosen' if best[m][:2] == (dim, lr) else ''}")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("cmd", choices=["prepare", "train", "select"])
    p.add_argument("--split")
    p.add_argument("--model")
    p.add_argument("--dim", type=int)
    p.add_argument("--lr", type=float)
    p.add_argument("--dest")
    p.add_argument("--epochs", type=int, default=EPOCHS)
    a = p.parse_args()
    if a.cmd == "prepare":
        prepare()
    elif a.cmd == "select":
        select()
    else:
        train(a.split, a.model, a.dim, a.lr, a.dest, a.epochs)
