"""Chemical similarity between drugs (structure), from the ChEMBL SMILES (Open Targets 26.06).

Morgan fingerprints (radius 2, 2048 bits) and Tanimoto similarity. For each molecule the
20 most similar neighbours with Tanimoto ≥ 0.4 are kept (data/ext/structure_edges.parquet):
  relation "structure", display_relation "≥0.6" or "0.4-0.6", score = Tanimoto.
The keys are the graph ones (DB:… if it has DrugBank, otherwise CHEMBL:…); prepare.py links them.

Useful for drugs with no known use (they resemble drugs that work in similar
diseases) and to look for molecules similar to a drug (possible derivatives).

    py -3 -m sources.structure
"""
import time
from pathlib import Path

import numpy as np
import pandas as pd

from sources.ids import OT, drugs

DEST = Path(__file__).parent.parent / "data" / "ext" / "structure_edges.parquet"
NEIGHBOURS, MINIMUM = 20, 0.4


def main():
    import sys
    sys.stdout.reconfigure(encoding="utf-8")
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import rdFingerprintGenerator
    RDLogger.DisableLog("rdApp.*")
    t0 = time.time()
    m = pd.read_parquet(OT / "drug_molecule", columns=["id", "canonicalSmiles"]).dropna()
    keys = drugs()
    m["key"] = [keys[i][0] if i in keys else f"CHEMBL:{i}" for i in m.id]
    m = m.drop_duplicates("key")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps, ok = [], []
    for smi, key in zip(m.canonicalSmiles, m.key):
        mol = Chem.MolFromSmiles(smi)
        if mol is None or mol.GetNumHeavyAtoms() < 5:  # ions and molecules that are too small
            continue
        fps.append(gen.GetFingerprint(mol))
        ok.append(key)
    print(f"{len(ok)} molecules with a fingerprint ({time.time() - t0:.0f} s)", flush=True)
    rows = []
    for i, fp in enumerate(fps):
        s = np.array(DataStructs.BulkTanimotoSimilarity(fp, fps))
        s[i] = 0
        top = np.argpartition(-s, NEIGHBOURS)[:NEIGHBOURS]
        top = top[s[top] >= MINIMUM]
        rows.extend((ok[i], ok[j], s[j]) for j in top)
        if i % 2000 == 0:
            print(f"\r  {i}/{len(fps)}", end="", flush=True)
    d = pd.DataFrame(rows, columns=["x", "y", "score"])
    # a single row per pair
    a, b = np.minimum(d.x, d.y), np.maximum(d.x, d.y)
    d = d.assign(x=a, y=b).drop_duplicates(["x", "y"])
    d["relation"] = "structure"
    d["display_relation"] = np.where(d.score >= 0.6, "≥0.6", "0.4-0.6")
    d["source"] = "ChEMBL (structure)"
    DEST.parent.mkdir(parents=True, exist_ok=True)
    d[["x", "y", "relation", "display_relation", "score", "source"]].to_parquet(DEST, index=False)
    print(f"\n{len(d):,} similar pairs ({(d.score >= 0.6).sum():,} with Tanimoto ≥ 0.6) -> {DEST} "
          f"({time.time() - t0:.0f} s)")


if __name__ == "__main__":
    main()
