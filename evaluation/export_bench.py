"""Exports the evaluation scores (data/eval/scores/*_extra_no_ctd_r5_rank.npz) to the OREVID-bench
format (benchmark/predictions/<model>/<task>.parquet: disease, drug, score).

    py -3 -m evaluation.export_bench
    py -3 -m evaluation.export_bench external    # external_ot2606 task (data/eval/external_scores.npz)
    py -3 -m evaluation.export_bench train-only  # OREVID trained without validation labels (…_rank_trainonly.npz)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from evaluation.txgnn import AREAS, CACHE, PRODUCTION

ROOT = Path(__file__).resolve().parent.parent
MODELS = {"LambdaRank_+_masked": "orevid_lambdarank",
          "graph_only_model": "orevid_no_popularity",
          "TxGNN_run_here": "txgnn"}


def main(suffix="", models=MODELS):
    sys.stdout.reconfigure(encoding="utf-8")
    for f in sorted((CACHE / "scores").glob(f"*_{PRODUCTION}_r5_rank{suffix}.npz")):
        split = next(s for s in ["complex_disease"] + AREAS if f.name.startswith(s))
        seed = f.name[len(split) + 1:].split("_")[0]
        task = f"zero_shot_{seed}" if split == "complex_disease" else f"area_{split}"
        z = np.load(f, allow_pickle=True)
        for key, name in models.items():
            if key not in z:
                continue
            S = z[key]
            dis, drugs = z["diseases"], z["drugs"]
            out = ROOT / "benchmark" / "predictions" / name
            out.mkdir(parents=True, exist_ok=True)
            pd.DataFrame({"disease": np.repeat(dis, len(drugs)).astype(np.int32),
                          "drug": np.tile(drugs, len(dis)).astype(np.int32),
                          "score": S.T.ravel().astype(np.float32)}).to_parquet(out / f"{task}.parquet", index=False)
        print(task)


if __name__ == "__main__" and len(sys.argv) == 1:
    main()
if __name__ == "__main__" and sys.argv[1:] == ["train-only"]:
    main("_trainonly", {"LambdaRank_+_masked": "orevid_lambdarank_train_only"})


def external():
    """Task «external_ot2606»: indications approved in Open Targets 26.06 that are not in PrimeKG.
    The models must be trained on the zero_shot_1 split (like TxGNN)."""
    import json
    z = np.load(CACHE / "external_scores.npz", allow_pickle=True)
    p = pd.read_parquet(CACHE / "external_positives.parquet")
    diseases, drugs = z["diseases"], z["drugs"]
    rows = []
    for d, known in zip(diseases, z["known"]):
        s = p[p.disease == d]
        rows.append({"disease": int(d), "positives": sorted(int(f) for f in s.drug),
                     "recent": sorted(int(f) for f in s[s["year"] >= 2018].drug),
                     "non_obvious": sorted(int(f) for f in s[~s.obvious].drug),
                     "exclude": sorted(int(f) for f in known)})
    task = {"name": "external_ot2606", "type": "external", "train_split": "zero_shot_1",
            "description": "Drug-disease pairs approved according to Open Targets 26.06 (ChEMBL) that are not in "
                           "PrimeKG. Candidates: all PrimeKG drugs except those PrimeKG already relates to the disease.",
            "drugs": [int(f) for f in drugs], "train_popularity": [int(x) for x in z["degree"]],
            "diseases": rows, "n_test_diseases": len(rows)}
    (ROOT / "benchmark" / "tasks" / "external_ot2606.json").write_text(json.dumps(task))
    for key, name in {"Ours_ranking": "orevid_lambdarank", "Ours_without_popularity": "orevid_no_popularity",
                      "TxGNN": "txgnn", "Popularity": "popularity"}.items():
        S = z[key]
        out = ROOT / "benchmark" / "predictions" / name
        out.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"disease": np.repeat(diseases, len(drugs)).astype(np.int32),
                      "drug": np.tile(drugs, len(diseases)).astype(np.int32),
                      "score": S.T.ravel().astype(np.float32)}).dropna().to_parquet(out / "external_ot2606.parquet", index=False)
    print(f"external_ot2606: {len(rows)} diseases, {len(p)} indications")


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "external":
    external()
