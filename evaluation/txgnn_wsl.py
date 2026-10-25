"""Trains TxGNN with its official code on the same splits as evaluation/txgnn.py and saves its
per-disease predictions, so they can be evaluated with our metrics.

It runs inside WSL with the TxGNN environment (Python 3.10, PyTorch 1.13 + CUDA, DGL 0.9.1,
pandas 1.5; see docs/TXGNN_EVALUATION.md):

    python evaluation/txgnn_wsl.py --split complex_disease --seed 1 --data ~/txgnn_data \
        --splits <repo>/data/eval --out <repo>/data/eval/txgnn_predictions

- complex_disease: the split is made by TxGNN's code (it is deterministic); our module
  reproduces it identically and we check it when reading the predictions.
- areas: TxGNN makes it without a seed; here we give it ours (data/eval/split_<area>_1.npy)
  as the «split» column of its kg_directed.csv, and it makes the validation as usual.
Hyperparameters: those of TxGNN's reproduce/train.py.

Output: <out>/TxGNN_<seed>_<split>[_external|_test].parquet with columns
relation, disease_id, drug_id, score.
"""
import argparse
import os
import shutil
import sys
import time

import numpy as np
import pandas as pd

AREAS = ["cell_proliferation", "mental_health", "cardiovascular", "anemia", "adrenal_gland",
         "autoimmune", "metabolic_disorder", "diabetes", "neurodigenerative"]


def prepare_area(data, split, splits, repo):
    """kg_directed.csv of the area with our split (2 = test)."""
    dest = os.path.join(data, f"{split}_kg", "kg_directed.csv")
    if os.path.exists(dest):
        return
    base = os.path.join(data, "kg_directed.csv")
    if not os.path.exists(base):
        raise SystemExit("Missing kg_directed.csv: run complex_disease first (TxGNN creates it).")
    df = pd.read_csv(base)
    label = np.load(os.path.join(splits, f"split_{split}_1.npy"))
    assert len(df) == len(label), (len(df), len(label))
    df["split"] = np.where(label == 2, "test", "train")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    df.to_csv(dest, index=False)
    files = os.path.join(data, "disease_files")
    if not os.path.exists(files):
        shutil.copytree(os.path.join(repo, "data", "disease_files"), files)


def prepare_custom(data, folder, label_path, repo):
    """T3: kg_directed.csv with our own split column (label 2 = test, i.e. never trained on), in TxGNN's
    disease-area folder layout (<area>_kg_frac<test_size>), so TxGNN's area-split code reads it unchanged."""
    dest = os.path.join(data, folder, "kg_directed.csv")
    if os.path.exists(dest):
        return
    df = pd.read_csv(os.path.join(data, "kg_directed.csv"))
    label = np.load(label_path)
    assert len(df) == len(label), (len(df), len(label))
    df["split"] = np.where(label == 2, "test", "train")
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    df.to_csv(dest, index=False)
    files = os.path.join(data, "disease_files")
    if not os.path.exists(files):
        shutil.copytree(os.path.join(repo, "data", "disease_files"), files)


def _id(x):
    """TxGNN's convert2str."""
    x = str(x)
    try:
        return x if "_" in x else str(float(x))
    except ValueError:
        return x


def minimal_kg(data):
    """After the first split, kg.csv is only read by TxGNN's evaluation, and only for the
    drug and disease names; reading it whole (~4 GB in memory) brings WSL down. The original
    is saved as kg_full.csv and a kg.csv with one row per drug and disease is left."""
    kg = os.path.join(data, "kg.csv")
    if os.path.getsize(kg) < 100e6 or not os.path.exists(os.path.join(data, "kg_directed.csv")):
        return
    os.replace(kg, os.path.join(data, "kg_full.csv"))
    n = pd.read_csv(os.path.join(data, "node.csv"), sep="	", dtype=str)
    n = n[n.node_type.isin(["drug", "disease"])]
    pd.DataFrame({"relation": "name", "x_id": n.node_id, "x_type": n.node_type, "x_name": n.node_name,
                  "y_id": n.node_id, "y_type": n.node_type, "y_name": n.node_name}).to_csv(kg, index=False)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--split", default="complex_disease")
    p.add_argument("--seed", type=int, default=1)
    p.add_argument("--data", required=True, help="TxGNN data folder (with kg.csv)")
    p.add_argument("--splits", required=True, help="the repository's data/eval")
    p.add_argument("--out", required=True)
    p.add_argument("--repo", default=os.path.expanduser("~/TxGNN"))
    p.add_argument("--epochs", type=int, default=500)
    p.add_argument("--test", action="store_true", help="no pretraining and 1 epoch (to test the pipeline)")
    p.add_argument("--label", help="T3: .npy split column (2 = test) for a custom split, read through the area mechanism")
    p.add_argument("--test-size", type=float, default=0.05, help="T3: only names TxGNN's cache folders (distinct per cutoff)")
    p.add_argument("--tag", help="T3: model / output name (TxGNN_<tag>)")
    p.add_argument("--diseases", help="CSV (node_index,node_id): score these diseases with the already trained model "
                                      "(external exam, evaluation/external_exam.py) instead of the test ones")
    a = p.parse_args()
    sys.path.insert(0, a.repo)
    from txgnn import TxData, TxGNN, TxEval

    t0 = time.time()
    # TxData downloads edges.csv if it is missing, but only the areas with TxGNN's own split
    # use it (which we do not use here).
    open(os.path.join(a.data, "edges.csv"), "a").close()
    if a.label:  # T3: custom split through the area mechanism (any area name; the folder is set by test_size)
        a.split = "anemia"
        prepare_custom(a.data, f"{a.split}_kg_frac{a.test_size}", a.label, a.repo)
    elif a.split in AREAS:
        prepare_area(a.data, a.split, a.splits, a.repo)
    if a.label:
        # T3: TxGNN's own end-of-training test pass (finetune -> «Testing...») builds a graph from its test edges, which in a
        # custom split may lack an edge type and makes DGL fail after training; it only prints numbers we do not use (the
        # model is chosen on validation and T3 scores the diseases separately), so it is skipped when it fails.
        txgnn_module = sys.modules["txgnn.TxGNN"]  # the module (the package exports a class of the same name)
        _evaluate_fb = txgnn_module.evaluate_fb

        def evaluate_fb(*args, mode="valid", **kw):
            try:
                return _evaluate_fb(*args, mode=mode, **kw)
            except Exception as e:  # noqa: BLE001 — only the unused test pass
                if mode != "test":
                    raise
                print(f"test pass skipped ({e})", flush=True)
                return ({}, {}, float("nan"), float("nan"), float("nan"), float("nan")), float("nan"), None, None
        txgnn_module.evaluate_fb = evaluate_fb
    data = TxData(data_folder_path=a.data)
    data.prepare_split(split=a.split, seed=a.seed, no_kg=False, test_size=a.test_size)
    name = (f"TxGNN_{a.tag}" if a.tag else f"TxGNN_{a.seed}_{a.split}") + ("_test" if a.test else "")
    model = TxGNN(data=data, weight_bias_track=False, proj_name="TxGNN", exp_name=name, device="cuda:0")
    saved = os.path.join(a.data, "models", name)
    if os.path.exists(os.path.join(saved, "model.pt")):
        # Our own model saved by a previous run that crashed while evaluating.
        model.load_pretrained(saved)
        print("model loaded from", saved, flush=True)
    else:
        model.model_initialize(n_hid=100, n_inp=100, n_out=100, proto=True, proto_num=3,
                               attention=False, sim_measure="all_nodes_profile", agg_measure="rarity")
        if not a.test:
            model.pretrain(n_epoch=1, learning_rate=1e-3, batch_size=1024, train_print_per_n=20)
        model.finetune(n_epoch=1 if a.test else a.epochs, learning_rate=5e-4, train_print_per_n=5, valid_per_n=20)
        os.makedirs(os.path.dirname(saved), exist_ok=True)
        model.save_model(saved)
        print(f"trained in {(time.time() - t0) / 60:.0f} min; saved to {saved}", flush=True)
    import gc
    import torch
    gc.collect()
    torch.cuda.empty_cache()
    minimal_kg(a.data)

    if a.diseases:
        # TxGNN's internal index of each disease (convert2str of the ids, as TxGNN does).
        df = data.df
        id2idx = dict(zip(df[df.x_type == "disease"].x_id.map(_id), df[df.x_type == "disease"].x_idx))
        id2idx.update(zip(df[df.y_type == "disease"].y_id.map(_id), df[df.y_type == "disease"].y_idx))
        listing = pd.read_csv(a.diseases, dtype=str)
        idxs = [int(id2idx[_id(i)]) for i in listing.node_id if _id(i) in id2idx]
        print(f"external exam: {len(idxs)} of {len(listing)} diseases", flush=True)
        d = TxEval(model=model).eval_disease_centric(disease_idxs=idxs, relation="indication", show_plot=False,
                                                     verbose=False, save_result=False, return_raw=False,
                                                     simulate_random=False)
        res = {"rev_indication": d}
        name += "_external"
    else:
        res = TxEval(model=model).eval_disease_centric(disease_idxs="test_set", show_plot=False, verbose=False,
                                                       save_result=False, return_raw=False, simulate_random=False)
    rows = []
    for rel in res:
        d = res[rel]
        for disease, pred in zip(d.index, d["Prediction"]):
            ids = list(pred.keys())
            rows.append(pd.DataFrame({"relation": rel[4:], "disease_id": str(disease),
                                      "drug_id": [str(i) for i in ids],
                                      "score": np.array([pred[i] for i in ids], np.float32)}))
    os.makedirs(a.out, exist_ok=True)
    pd.concat(rows).to_parquet(os.path.join(a.out, f"{name}.parquet"), index=False)
    print(f"done in {(time.time() - t0) / 60:.0f} min: {name}", flush=True)


if __name__ == "__main__":
    main()
