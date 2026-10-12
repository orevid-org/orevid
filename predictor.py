"""Production model for the «disease» command: the winner of the comparison with TxGNN
(docs/TXGNN_EVALUATION.md): LightGBM with a ranking objective (LambdaRank) over metapaths
(the «extra» set, without any Open Targets drug-disease relation), popularity and drug
degrees, with masked training. Scores are calibrated to a logit (Platt) so that the
adjustments of the «disease» command (failed trials, direction, literature) make sense on them.
Also, a model without popularity (logistic regression) for the list of drugs with no
previous use, where popularity would sink them (docs/NOVELTY.md).
Here it is trained with all the diseases of the prepared graph (data/prep).

    py -3 predictor.py              # trains and saves data/prep/predictor.pkl (~10 min)

The candidates are the PrimeKG drugs (the ones used in the evaluation).
"""
import pickle
import sys
import time

import numpy as np
import scipy.sparse as sp

from evaluation.txgnn import (PRODUCTION, FEATURE_SETS, TxGraph, fit_model, fit_ranker, masked_rows, drug_degrees, score_paths)
from graph import Graph, INDICATION, PREP

MODEL_FILE = PREP / "predictor.pkl"
FEATURE_SET, ROUNDS = PRODUCTION, 5


class ProductionGraph(TxGraph):
    """The whole prepared graph (data/prep) with the metapaths of the evaluation."""

    def __init__(self):
        Graph.__init__(self)
        self._sim = {}
        # Candidates: PrimeKG drugs (the evaluation ones); the Open Targets ones stay as graph nodes.
        self.drugs = np.flatnonzero((self.types == "drug") & (self.origin == "PrimeKG"))


def train(seed=0):
    t0 = time.time()
    g = ProductionGraph()
    paths = FEATURE_SETS[FEATURE_SET]
    drug_pos = np.full(g.N, -1)
    drug_pos[g.drugs] = np.arange(len(g.drugs))
    I = g.M(INDICATION).tocoo()
    is_drug = np.isin(I.row, g.drugs)
    x, y = I.row[is_drug], I.col[is_drug]  # drug, disease
    known = np.ones(len(x), bool)
    degree = np.diff(g.M(INDICATION).indptr)[g.drugs].astype(np.float64)
    gf = drug_degrees(g)
    diseases = np.unique(y)
    rng = np.random.default_rng(seed)
    rows_X, rows_y = [], []
    for i in range(0, len(diseases), 256):
        block = diseases[i:i + 256]
        F = score_paths(g, block, paths)
        for c, d in enumerate(block):
            lab = np.zeros(len(g.drugs), bool)
            lab[drug_pos[x[y == d]]] = True
            pop = degree - lab  # its own indication does not count
            sel = lab | (rng.random(len(lab)) < 0.1)
            rows_X.append(np.column_stack([F[sel, c, :], pop[sel], gf[sel]]))
            rows_y.append(lab[sel])
        print(f"\r  {min(i + 256, len(diseases))}/{len(diseases)} diseases", end="", flush=True)
    Xm, Ym, Gm = masked_rows(g, np.array([], int), known, x, y, drug_pos, paths, ROUNDS, seed)
    X = np.concatenate(rows_X + [Xm])
    Y = np.concatenate(rows_y + [Ym])
    rank = fit_ranker(X, Y, np.r_[[len(v) for v in rows_y], Gm], seed)
    # Platt: ranking score -> logit (over the training rows)
    from sklearn.linear_model import LogisticRegression
    platt = LogisticRegression().fit(rank(X)[:, None], Y)
    nc = len(paths)
    no_popularity, _ = fit_model(np.concatenate(rows_X)[:, :nc], np.concatenate(rows_y))
    MODEL_FILE.write_bytes(pickle.dumps({"rank": rank, "platt": (float(platt.coef_[0, 0]), float(platt.intercept_[0])),
                                         "no_popularity": no_popularity, "paths": paths, "feature_set": FEATURE_SET,
                                         "rounds": ROUNDS}))
    print(f"\nModel trained with {len(diseases)} diseases ({len(X):,} rows) in {time.time() - t0:.0f} s -> {MODEL_FILE}")


def developed_drugs(g):
    """Drugs approved or in clinical investigation according to DrugBank (not the compounds that are only
    «experimental», which are often crystallography ligands)."""
    if g.descriptions is None or "group" not in g.descriptions:
        return set(g.drugs)
    group = g.descriptions["group"].astype(str)
    return set(group.index[group.str.contains(APPROVED_HUMAN + "|investigational", regex=True).to_numpy()])


def diagnostic_drugs(g):
    """Contrast media and diagnostic agents (DrugBank category): they are not therapeutic candidates."""
    if g.descriptions is None or "category" not in g.descriptions:
        return set()
    cat = g.descriptions["category"].astype(str)
    return set(cat.index[cat.str.contains("Contrast Media|Diagnostic|Skin Test Antigen", case=False, regex=True).to_numpy()])


def drug_status(g):
    """Development status according to DrugBank: «approved», «investigational», «withdrawn»,
    «experimental» or «unknown» (to label the candidates on the web)."""
    res = {}
    if g.descriptions is None or "group" not in g.descriptions:
        return res
    import re
    for i, s_ in g.descriptions["group"].astype(str).items():
        approved = re.search(APPROVED_HUMAN, s_) is not None  # «vet_approved» is not approved in humans
        if "withdrawn" in s_:
            # «approved and withdrawn»: withdrawn in some market (e.g. thalidomide); without «approved»: withdrawn
            res[i] = "withdrawn in some market" if approved else "withdrawn"
        elif approved:
            res[i] = "approved"
        elif "investigational" in s_:
            res[i] = "investigational"
        elif "vet_approved" in s_:
            res[i] = "veterinary only"
        elif "experimental" in s_:
            res[i] = "experimental"
    return res


APPROVED_HUMAN = r"(?<!vet_)approved"


# Relations that count as «known use» for the list of drugs with no previous use (review: before,
# only PrimeKG, and drugs approved according to Open Targets showed up there).
KNOWN_USE = "indication|clinical:approved|off-label use"


class Predictor:
    def __init__(self):
        if not MODEL_FILE.exists():
            raise SystemExit("Model missing: run  py -3 predictor.py")
        d = pickle.loads(MODEL_FILE.read_bytes())  # our own file, generated by train()
        self.rank, self.platt, self.no_popularity, self.paths = d["rank"], d["platt"], d["no_popularity"], d["paths"]
        self.g = ProductionGraph()
        self.degree = np.diff(self.g.M(INDICATION).indptr)[self.g.drugs].astype(np.float64)
        self.gf = drug_degrees(self.g)

    def raw(self, disease):
        F = score_paths(self.g, np.array([disease]), self.paths)[:, 0, :]
        self.last_no_popularity = self.no_popularity(F)  # for the list of drugs with no previous use
        return self.rank(np.column_stack([F, self.degree, self.gf])), F

    def predict(self, disease):
        """(drugs, score, features per metapath) for a disease.

        Score = specificity: model score normalised per disease (z) minus half of the drug's
        mean over reference diseases. It removes the generic candidates (which show up
        everywhere) at a small cost in hits; weight chosen with a fixed criterion and confirmed out of
        sample (docs/DECISIONS.md, P13). It is not a probability."""
        s, F = self.raw(disease)
        z = (s - s.mean()) / (s.std() + 1e-9)
        if REFERENCE.exists():
            z = z - SPECIFICITY_WEIGHT * np.load(REFERENCE)
        return self.g.drugs, z, dict(zip(self.paths, F.T))


REFERENCE = PREP / "predictor_reference.npy"
SPECIFICITY_WEIGHT = 0.5


def reference(n=2000, seed=0):
    """Mean normalised score (z) of each drug over n random diseases with data."""
    from graph import GENES_ALL
    P = Predictor()
    g = P.g
    pkg = np.flatnonzero((g.types == "disease") & (g.origin == "PrimeKG"))
    with_data = pkg[(np.diff(g.M(GENES_ALL).indptr)[pkg] > 0) | (np.diff(g.M("disease_phenotype_positive").indptr)[pkg] > 0)]
    sample = np.random.default_rng(seed).choice(with_data, min(n, len(with_data)), replace=False)
    total = np.zeros(len(g.drugs))
    t0 = time.time()
    for k, d in enumerate(sample):
        s, _ = P.raw(d)
        total += (s - s.mean()) / (s.std() + 1e-9)
        if (k + 1) % 250 == 0:
            print(f"  {k + 1}/{len(sample)} ({time.time() - t0:.0f} s)", flush=True)
    np.save(REFERENCE, total / len(sample))
    print(f"specificity reference ({len(sample)} diseases) -> {REFERENCE}")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    if len(sys.argv) > 1 and sys.argv[1] == "reference":
        reference()
    else:
        train()
        print("Remember: py -3 predictor.py reference  (it must be redone if the model is retrained)")
