"""Reproduces the TxGNN evaluation (Huang et al., Nature Medicine 2024) and runs our model on it.

Same graph (the PrimeKG kg.csv that TxGNN uses, with grouped diseases), same disease splits and
same metrics as the notebook in their repository (reproduce/gather result-txgnn.ipynb), so the
results can be compared with the published figures (resources/txgnn/published_results.csv, from
reproduce/result_more_metrics.csv; MIT licence).

Splits:
- complex_disease («zero-shot»): 5 % of the diseases with drugs (shuffled with the seed) lose
  all their drug-disease edges; of the other relations, 5 % (test) and 11.875 % (validation)
  are removed. Exactly reproducible: it is TxGNN's code.
- disease areas (anemia, neurodigenerative…): all the diseases of a branch of the DOID
  ontology (lists in resources/txgnn/) lose their drug-disease edges, and on top of that random
  edges of the 2-hop neighbourhood are removed until 5 % of the graph is reached. TxGNN uses no
  seed there, so here it is an approximation (with a fixed seed).

Metrics (for the «indication» relation), as in the TxGNN notebook:
- AUPRC / AUROC: for each test disease, the positives and the same number of negatives drawn
  at random (seed 42, with replacement) among the drugs with some indication in training;
  all diseases are pooled.
- Recall@100, Precision@K (K = number of positives), Avg Rank (mean normalised rank of the
  positives; 1 = at the top): over all the drugs of the graph.

    py -3 -m evaluation.txgnn                       # zero-shot, seeds 1-5
    py -3 -m evaluation.txgnn --split anemia --seeds 1
    py -3 -m evaluation.txgnn --split areas         # the 9 areas
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import rankdata
from sklearn.metrics import average_precision_score, roc_auc_score

from graph import Graph, METAPATHS, INDICATION

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "txgnn"
CACHE = ROOT / "data" / "eval"
LISTS = ROOT / "resources" / "txgnn"
PUBLISHED = LISTS / "published_results.csv"

DD = ["contraindication", "indication", "off-label use"]
AREAS = ["cell_proliferation", "mental_health", "cardiovascular", "anemia", "adrenal_gland",
         "autoimmune", "metabolic_disorder", "diabetes", "neurodigenerative"]

# Metapaths that only use PrimeKG (see graph.METAPATHS).
METAPATHS.update({
    "target (all)": ("direct", "drug_protein", "disease_protein",
                     "drug → target, enzyme or transporter = disease gene"),
    "similar phenotype": ("similar", INDICATION, "disease_phenotype_positive",
                          "indicated for a disease that shares phenotypes"),
    "contra similar": ("similar", "contraindication", "disease_protein",
                       "contraindicated in a disease that shares genes"),
    "off-label similar": ("similar", "off-label use", "disease_protein",
                          "off-label use in a disease that shares genes"),
    "off-label ontology": ("ontology", "off-label use", None, "off-label use in a neighbouring disease"),
})
# New metapaths for the diseases with few shared genes (mental health, diabetes).
# «chain»: product of relations; «via drug»: similar drugs (by the given relation)
# to those scored by another metapath.
METAPATHS.update({
    "effect = phenotype": ("direct", "drug_effect", "disease_phenotype_positive",
                           "the drug's adverse effects are phenotypes of the disease"),
    "target → phenotype": ("ppi", "drug_protein", "disease_phenotype_positive",
                           "target associated with a phenotype of the disease", "phenotype_protein"),
    "target → pathway": ("chain", ["drug_protein", "pathway_protein", "pathway_protein", "disease_protein"],
                         None, "target in the same pathway (Reactome) as a disease gene"),
    "target → process": ("chain", ["drug_protein", "bioprocess_protein", "bioprocess_protein", "disease_protein"],
                         None, "target in the same biological process (GO) as a disease gene"),
    "ontology 2-step": ("chain", [INDICATION, "disease_disease", "disease_disease"],
                        None, "indicated for a disease 2 steps away in the ontology"),
    "drug (targets) · similar": ("via drug", "drug_protein", "similar",
                                 "shares targets with drugs indicated in diseases with genes in common"),
    "drug (effects) · similar": ("via drug", "drug_effect", "similar",
                                 "shares adverse effects with drugs indicated in diseases with genes in common"),
    "drug (effects) · phenotype": ("via drug", "drug_effect", "similar phenotype",
                                   "shares adverse effects with drugs indicated in diseases with phenotypes in common"),
    "drug (interactions) · similar": ("via drug", "drug_drug", "similar",
                                      "interacts with drugs indicated in diseases with genes in common"),
    "drug (targets) · ontology": ("via drug", "drug_protein", "ontology",
                                  "shares targets with drugs indicated in neighbouring diseases"),
})
NEW_PATHS = ["effect = phenotype", "target → phenotype", "target → pathway", "target → process", "ontology 2-step",
             "drug (targets) · similar", "drug (effects) · similar", "drug (effects) · phenotype",
             "drug (interactions) · similar", "drug (targets) · ontology"]

# Extra OREVID sources (data/prep). Biological knowledge only: no drug-disease relation
# (no trials, no CTD therapeutic, no PubTator «treats», no FAERS), because an equivalent
# unmapped disease could carry the answer into the graph.
EXTRA = ["moa", "ctd_gene:increases", "ctd_gene:decreases", "string", "string:≥900",
         "ot_assoc:genetic_association", "ot_assoc:genetic_literature", "ot_assoc:somatic_mutation",
         "ot_assoc:animal_model", "ot_assoc:rna_expression", "structure"]
GENES_EXTRA = ("disease_protein|ot_assoc:genetic_association|ot_assoc:genetic_literature|"
               "ot_assoc:somatic_mutation|ot_assoc:animal_model|ot_assoc:rna_expression")
METAPATHS.update({
    "target+ all": ("direct", "drug_protein|moa", GENES_EXTRA, "DrugBank + ChEMBL targets = associated gene (PrimeKG + Open Targets)"),
    "interaction STRING all": ("ppi", "drug_protein|moa", GENES_EXTRA, "target ↔ (STRING ≥ 700) associated gene", "string"),
    "effect CTD all": ("direct", "ctd_gene:increases|ctd_gene:decreases", GENES_EXTRA,
                       "the drug changes a gene associated with the disease (CTD)"),
    "similar all (OT)": ("similar", INDICATION, GENES_EXTRA, "indicated for a disease with genes in common (PrimeKG + OT)"),
    "drug (CTD) · similar all": ("via drug", "ctd_gene:increases|ctd_gene:decreases", "similar all (OT)",
                                 "affects the same genes (CTD) as drugs indicated in similar diseases"),
})
PATHS_EXTRA = ["target+ all", "interaction STRING all", "effect CTD all", "similar all (OT)",
               "drug (CTD) · similar all"]

# Drug degrees in other relations (drug features for the gradient boosting).
DEGREES = ["drug_protein", "drug_effect", "drug_drug", "contraindication", "off-label use"]
PATHS = ["target", "target (all)", "interaction", "similar", "similar phenotype", "ontology",
         "contra similar", "off-label similar", "off-label ontology"]
# Chemical structure (sources/structure.py): drugs that are chemically alike (Tanimoto ≥ 0.4).
METAPATHS.update({
    "drug (structure) · similar": ("via drug", "structure", "similar",
                                   "chemically similar to drugs indicated in diseases with genes in common"),
    "drug (structure) · phenotype": ("via drug", "structure", "similar phenotype",
                                     "chemically similar to drugs indicated in diseases with phenotypes in common"),
    "drug (structure) · ontology": ("via drug", "structure", "ontology",
                                    "chemically similar to drugs indicated in neighbouring diseases"),
})
PATHS_STRUCTURE = ["drug (structure) · similar", "drug (structure) · phenotype", "drug (structure) · ontology"]
# Production feature set (predictor, web, external exam, OREVID-bench). Without CTD: it adds nothing
# (ablation on the 14 splits and on the external exam) and it is the only source with a
# non-commercial licence (DECISIONS.md, P20).
PRODUCTION = "extra_no_ctd"
FEATURE_SETS = {"base": PATHS, "new": PATHS + NEW_PATHS, "extra": PATHS + NEW_PATHS + PATHS_EXTRA,
                "structure": PATHS + NEW_PATHS + PATHS_EXTRA + PATHS_STRUCTURE,
                # ablation: without CTD (drug-gene relations taken from the literature; docs/DECISIONS.md A14)
                "extra_no_ctd": PATHS + NEW_PATHS + [c for c in PATHS_EXTRA if "CTD" not in c]}

# Model score names (display names in the results tables) -> keys in the score npz files
# (data/eval/scores/), read by evaluation/novelty.py, evaluation/evidence.py and evaluation/export_bench.py.
SCORE_KEYS = {"popularity": "popularity",
              "model: graph only": "graph_only_model",
              "gbm: graph + popularity + degrees": "gbm_graph_+_popularity_+_degrees",
              "gbm + masked training": "gbm_+_masked_training",
              "ranking (LambdaRank) + masked": "LambdaRank_+_masked",
              "TxGNN (run here)": "TxGNN_run_here"}
# Prefix of the metrics computed only on the test diseases without a near-duplicate.
NO_NEAR_DUPLICATES = "[no near duplicates] "


# --- TxGNN graph ----------------------------------------------------------------------

def load_kg():
    """TxGNN's kg.csv, with one row per edge like its `kg_directed.csv`: in relations between
    nodes of the same type, a single direction per pair; in the rest, the direction of the
    relation's first record. Row order and row index are kept (the splits need them)."""
    f = CACHE / "kg_directed.parquet"
    if f.exists():
        return pd.read_parquet(f)
    if not (RAW / "kg.csv").exists():
        raise SystemExit("Missing data/raw/txgnn/kg.csv (Harvard Dataverse, file 7144484). See the README.")
    t = time.time()
    kg = pd.read_csv(RAW / "kg.csv", usecols=["relation", "display_relation", "x_index", "x_type",
                                               "x_name", "y_index", "y_type", "y_name"], low_memory=False)
    print(f"kg.csv: {len(kg):,} rows ({time.time() - t:.0f} s)")
    keep = []
    for rel in kg.relation.unique():
        d = kg[kg.relation == rel]
        if "_" in rel and rel.split("_")[0] == rel.split("_")[1]:
            pair = np.minimum(d.x_index, d.y_index).astype(np.int64) * 10**6 + np.maximum(d.x_index, d.y_index)
            keep.append(d.index[~pair.duplicated().to_numpy()])
        else:
            keep.append(d.index[(d.x_type == d.x_type.iloc[0]).to_numpy()])
    kg = kg.loc[np.sort(np.concatenate(keep))]
    CACHE.mkdir(parents=True, exist_ok=True)
    kg.to_parquet(f)
    return kg


def _sample(n, frac, seed):
    """Positions chosen by pandas.DataFrame.sample(frac=frac, random_state=seed)."""
    return np.random.RandomState(seed).choice(n, size=int(round(frac * n)), replace=False)


def split_edges(kg, split, seed):
    """Labels each row: 0 training, 1 validation, 2 test (like TxGNN's create_split)."""
    f = CACHE / f"split_{split}_{seed}.npy"
    if f.exists():
        return np.load(f)
    label = np.zeros(len(kg), np.int8)
    is_dd = kg.relation.isin(DD).to_numpy()
    rel = kg.relation.to_numpy()
    if split == "complex_disease":
        # complex_disease_fold: the diseases are shuffled in order of appearance.
        diseases = pd.unique(kg.y_index.to_numpy()[is_dd])
        np.random.seed(seed)
        np.random.shuffle(diseases)
        n = len(diseases)
        cuts = [int(0.83125 * n), int((0.83125 + 0.11875) * n)]
        val, test = set(diseases[cuts[0]:cuts[1]]), set(diseases[cuts[1]:])
        y = kg.y_index.to_numpy()
        label[is_dd & np.isin(y, list(val))] = 1
        label[is_dd & np.isin(y, list(test))] = 2
        for r in pd.unique(rel[~is_dd]):
            pos = np.flatnonzero(rel == r)
            p = _sample(len(pos), 0.05, seed)
            label[pos[p]] = 2
            rest = np.setdiff1d(np.arange(len(pos)), p)  # in the original order
            v = _sample(len(rest), 0.11875 / 0.95, 1)
            label[pos[rest[v]]] = 1
    elif split in AREAS:
        label = _split_area(kg, split, is_dd, rel)
    else:
        raise SystemExit(f"Unknown split: {split}")
    np.save(f, label)
    return label


def _id(x):
    """TxGNN's convert2str: numeric ids as float text («8840» -> «8840.0»)."""
    x = str(x)
    try:
        return x if "_" in x else str(float(x))
    except ValueError:
        return x


def area_diseases(split):
    """Indices of the diseases of an area. As in TxGNN's process_disease_area_split, they are
    matched by node_id (MONDO), because some lists carry indices from another version of the
    graph; a grouped node («8840_12166») counts for each of its ids."""
    nodes = pd.read_csv(RAW / "node.csv", sep="\t", dtype=str)
    nodes = nodes[nodes.node_type == "disease"]
    id2idx = {}
    for i, idx in zip(nodes.node_id, nodes.node_index.astype(int)):
        id2idx[_id(i)] = idx
        if "_" in i:
            for part in i.split("_"):
                id2idx.setdefault(_id(part), idx)
    d = pd.read_csv(LISTS / f"{split}.csv", dtype=str)
    return np.unique([id2idx[_id(i)] for i in d.node_id if _id(i) in id2idx])


def _split_area(kg, split, is_dd, rel, frac=0.05, seed=0):
    """Approximation of TxGNN's DataSplitter.get_test_kg_for_disease + random_fold."""
    nodes = area_diseases(split)
    x, y = kg.x_index.to_numpy(), kg.y_index.to_numpy()
    test = is_dd & (np.isin(x, nodes) | np.isin(y, nodes))
    # 2-hop neighbourhood; random edges between neighbourhood nodes up to 5 % of the graph
    # (TxGNN counts both directions: 5 % of 2 × edges, and picks from either direction).
    N = int(max(x.max(), y.max())) + 1
    A = sp.csr_matrix((np.ones(len(kg)), (x, y)), shape=(N, N))
    A = (A + A.T).tocsr()
    inside = np.zeros(N, bool)
    inside[nodes] = True
    for _ in range(2):
        inside |= np.asarray(A[inside].sum(axis=0)).ravel() > 0
    candidates = np.flatnonzero(inside[x] & inside[y] & ~test)
    n_dd = int(test.sum())
    n_random = int(round(2 * len(kg) * frac)) - 2 * n_dd
    # Each row here is one of TxGNN's two directions; only half of them fall in the direction
    # that is kept, so half as many are removed.
    pick = np.random.default_rng(seed).choice(candidates, min(len(candidates), max(n_random, 0) // 2), replace=False)
    test[pick] = True
    label = np.where(test, 2, 0).astype(np.int8)
    # random_fold(train_val, [0.875, 0.125, 0.0]): validation of 12.5 % of each relation (seed 1)
    for r in pd.unique(rel):
        pos = np.flatnonzero((rel == r) & (label == 0))
        label[pos[_sample(len(pos), 0.125, 1)]] = 1
    return label


# --- our model on the training graph -----------------------------------------------------

class TxGraph(Graph):
    """OREVID's Graph, but with the matrices built from the training edges.
    With extra=True the EXTRA relations of the prepared graph (data/prep) are added; it uses
    the same node indices for the PrimeKG nodes. The candidate drugs are always the PrimeKG
    ones."""

    def __init__(self, kg, train, types, extra=False):
        self.types = types
        self.drugs = np.flatnonzero(types == "drug")
        self.N = len(types)
        if extra:
            self.N = len(pd.read_pickle(ROOT / "data" / "prep" / "nodes.pkl"))
        self._w, self._total, self._attr = {}, None, None
        d = kg[train]
        self._m = {}
        for rel, s in d.groupby("relation"):
            self._m[rel] = self._symmetric(s)
        for disp, s in d[d.relation == "drug_protein"].groupby("display_relation"):
            self._m[f"drug_protein:{disp}"] = self._symmetric(s)
        if extra:
            import json
            meta = json.loads((ROOT / "data" / "prep" / "relations.json").read_text("utf-8"))["relations"]
            for rel in EXTRA:
                if rel in meta:
                    self._m[rel] = sp.load_npz(ROOT / "data" / "prep" / meta[rel]["file"]).tocsr()
        self.relations = dict.fromkeys(self._m)
        self._sim = {}

    def _symmetric(self, s):
        x, y = s.x_index.to_numpy(), s.y_index.to_numpy()
        m = sp.coo_matrix((np.ones(2 * len(s)), (np.r_[x, y], np.r_[y, x])), shape=(self.N, self.N)).tocsr()
        m.data[:] = 1
        return m

    def mask(self, diseases):
        """Copy of the graph without any drug-disease edge of these diseases."""
        import copy
        h = copy.copy(self)
        keep = np.ones(self.N)
        keep[diseases] = 0
        D = sp.diags(keep)
        h._m = {k: v for k, v in self._m.items() if "|" not in k}
        for r in DD:
            h._m[r] = (D @ self.M(r) @ D).tocsr()
            h._m[r].eliminate_zeros()
        h._w = {k: v for k, v in self._w.items() if not any(r in k for r in DD)}
        return h

    def drug_similarity(self, rel):
        """Drug × drug: shared targets/effects (normalised by degree) or direct relation."""
        if rel not in self._sim:
            f = self.drugs
            if rel in ("drug_drug", "structure"):
                S = self.W(rel)[f][:, f]
            else:
                W = self.W(rel)
                S = (W[f] @ W[:, f]).tocsr()
            S = S - sp.diags(S.diagonal())
            self._sim[rel] = S.tocsr()
        return self._sim[rel]

    def scores(self, diseases, paths=None, strict=False):
        diseases = np.asarray(diseases)
        own = [c for c in paths if METAPATHS[c][0] not in ("chain", "via drug")]
        base = [METAPATHS[c][2] for c in paths if METAPATHS[c][0] == "via drug"]
        res = super().scores(diseases, list(dict.fromkeys(own + base)), strict)
        for c in paths:
            kind, rel, rel2, *_ = METAPATHS[c]
            if kind == "chain":
                m = self.W(rel[-1])[:, diseases]
                for r in reversed(rel[1:-1]):
                    m = self.W(r) @ m
                # Without returning to the same disease (e.g. d -> neighbour -> d in «ontology 2-step»):
                # it would use its own indications (bug found in the independent review).
                m = m.tocoo()
                own_disease = m.row == diseases[m.col]
                m = sp.csr_matrix((m.data[~own_disease], (m.row[~own_disease], m.col[~own_disease])), shape=m.shape)
                m = self.W(rel[0]) @ m
                res[c] = np.asarray(m[self.drugs].todense())
            elif kind == "via drug":
                res[c] = np.asarray(self.drug_similarity(rel) @ res[rel2])
        return {c: res[c] for c in paths}


def node_types(kg):
    N = int(max(kg.x_index.max(), kg.y_index.max())) + 1
    t = np.empty(N, dtype=object)
    t[kg.x_index.to_numpy()] = kg.x_type.to_numpy()
    t[kg.y_index.to_numpy()] = kg.y_type.to_numpy()
    return t.astype(str)


def score_paths(g, diseases, paths, block=128):
    """Matrix (drugs × diseases × metapaths)."""
    res = np.zeros((len(g.drugs), len(diseases), len(paths)), np.float32)
    for i in range(0, len(diseases), block):
        p = g.scores(diseases[i:i + block], paths)
        for k, c in enumerate(paths):
            res[:, i:i + block, k] = p[c]
    return res


def fit_model(X, y):
    from sklearn.linear_model import LogisticRegression
    scale = np.array([X[X[:, j] > 0, j].mean() if (X[:, j] > 0).any() else 1.0 for j in range(X.shape[1])])
    clf = LogisticRegression(class_weight="balanced", max_iter=3000).fit(np.arcsinh(X / scale), y)
    return Logit(clf, scale), clf


class Logit:
    """Logit of the logistic regression on scaled features (can be saved with pickle)."""

    def __init__(self, clf, scale):
        self.clf, self.scale = clf, scale

    def __call__(self, Z):
        return self.clf.decision_function(np.arcsinh(Z / self.scale))


class Probability:
    """Probability of the positive class (can be saved with pickle, unlike a lambda)."""

    def __init__(self, clf):
        self.clf = clf

    def __call__(self, Z):
        return self.clf.predict_proba(Z)[:, 1]


def fit_gbm(X, y, seed=0):
    """Gradient boosting: it can learn interactions (e.g. that popularity only helps if there
    is a biological signal). Trees do not depend on the scale of the features."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                         min_samples_leaf=50, l2_regularization=1.0, random_state=seed)
    clf.fit(X, y)
    return Probability(clf)


def drug_degrees(g):
    """How well known each drug is in the training graph, apart from the indications."""
    return np.column_stack([np.diff(g.M(r).indptr)[g.drugs] if r in g.relations else np.zeros(len(g.drugs))
                            for r in DEGREES]).astype(np.float64)


def near_duplicates(g, diseases, positives, drug_pos):
    """For each disease, the maximum fraction of its drugs that a disease ≤ 2 steps away in
    the ontology already has as indications (in the training graph)."""
    O = g.M("disease_disease")
    O2 = (O + O @ O).tocsr()
    I = g.M(INDICATION)
    res = []
    for d, pos in zip(diseases, positives):
        neighbors = O2.indices[O2.indptr[d]:O2.indptr[d + 1]]
        neighbors = neighbors[neighbors != d]
        pos = set(g.drugs[pos])
        res.append(max((len(pos.intersection(I.indices[I.indptr[v]:I.indptr[v + 1]])) / len(pos)
                        for v in neighbors), default=0.0))
    return np.array(res)


# --- metrics from the TxGNN notebook ---------------------------------------------------

def metrics(scores, positives, pool_neg, tie_seed=0):
    """scores: (drugs × diseases); positives: list of arrays of drug positions;
    pool_neg: mask of drugs that can be negatives (some indication in training).
    Ties are broken at random (TxGNN gives continuous scores, without ties)."""
    rng = np.random.default_rng(tie_seed)
    pp, pn, rec, prk, rank = [], [], [], [], []
    n = scores.shape[0]
    for c, pos in enumerate(positives):
        s = scores[:, c] + rng.random(n) * 1e-9 * (np.abs(scores[:, c]).max() + 1e-9)
        np.random.seed(42)
        options = np.intersect1d(np.setdiff1d(np.arange(n), pos), np.flatnonzero(pool_neg))
        neg = np.random.choice(options, len(pos))
        pp += s[pos].tolist()
        pn += s[neg].tolist()
        order = np.argsort(s)[::-1]
        rec.append(len(np.intersect1d(order[:100], pos)) / len(pos))
        prk.append(len(np.intersect1d(order[:len(pos)], pos)) / len(pos))
        rank.append(np.mean(np.argsort(np.argsort(s))[pos]) / n)
    y = [1] * len(pp) + [0] * len(pn)
    return {"AUPRC": average_precision_score(y, pp + pn), "AUROC": roc_auc_score(y, pp + pn),
            "Recall@100": np.mean(rec), "Precision@K": np.mean(prk), "Avg Rank": np.mean(rank)}


# --- run -------------------------------------------------------------------------------

def masked_rows(g, test_diseases, known, x, y, drug_pos, paths, rounds, seed, seeds_per_round=40, maximum=400):
    """Training rows with masked neighbourhoods (see evaluate). No test disease ever enters."""
    rng = np.random.default_rng(1000 + seed)
    with_drugs = pd.unique(y[known])
    rows_X, rows_y = [], []
    for _ in range(rounds):
        round_seeds = rng.choice(with_drugs, seeds_per_round, replace=False)
        near = g._nearby("disease_protein", round_seeds).tocsc()
        ball, remove = set(), []
        for k in range(len(round_seeds)):
            v = near.indices[near.indptr[k]:near.indptr[k + 1]]
            v = v[g.types[v] == "disease"]
            if len(v) > maximum:  # nodes that are too general (e.g. «cancer»)
                continue
            remove.extend(v)
            ball.update(v)
        h = g.mask(np.array(sorted(ball | set(round_seeds))))
        labelled = np.array([d for d in with_drugs if d in ball or d in round_seeds])
        labelled = labelled[~np.isin(labelled, test_diseases)]
        degree = np.diff(h.M(INDICATION).indptr)[h.drugs].astype(np.float64)
        gf = drug_degrees(h)
        for i in range(0, len(labelled), 256):
            block = labelled[i:i + 256]
            F = score_paths(h, block, paths)
            for c, d in enumerate(block):
                lab = np.zeros(len(g.drugs), bool)
                lab[drug_pos[x[known & (y == d)]]] = True
                sel = lab | (rng.random(len(lab)) < 0.1)
                rows_X.append(np.column_stack([F[sel, c, :], degree[sel], gf[sel]]))
                rows_y.append(lab[sel])
    return np.concatenate(rows_X), np.concatenate(rows_y), np.array([len(v) for v in rows_y])


def fit_ranker(X, y, groups, seed=0):
    """LightGBM with a ranking objective (LambdaRank): it learns to order the drugs within each
    disease (group), instead of classifying each pair separately."""
    import lightgbm as lgb
    clf = lgb.LGBMRanker(objective="lambdarank", n_estimators=300, learning_rate=0.05, num_leaves=31,
                         min_child_samples=50, reg_lambda=1.0, random_state=seed, verbose=-1)
    clf.fit(X, y.astype(int), group=groups)
    return clf.predict


def evaluate(kg, types, split, seed, n_train=1500, feature_set="base", rounds=0, rank=False, train_only=False):
    t0 = time.time()
    paths = FEATURE_SETS[feature_set]
    label = split_edges(kg, split, seed)
    g = TxGraph(kg, label == 0, types, extra=(feature_set in ("extra", "structure", "extra_no_ctd")))
    drug_pos = np.full(g.N, -1)
    drug_pos[g.drugs] = np.arange(len(g.drugs))

    ind = kg.relation.to_numpy() == "indication"
    x, y = kg.x_index.to_numpy(), kg.y_index.to_numpy()
    test_ind = ind & (label == 2)
    if split in AREAS:
        test_ind &= np.isin(y, area_diseases(split))
    test_diseases = pd.unique(y[test_ind])
    test_pos = [drug_pos[x[test_ind & (y == d)]] for d in test_diseases]

    # Popularity: the drug's indications in the training graph.
    degree = np.diff(g.M(INDICATION).indptr)[g.drugs].astype(np.float64)
    pool_neg = degree > 0

    # Diseases to train the model: with known indications (training or validation). With train_only,
    # only the training edges are used as labels, as TxGNN does (it uses validation only for early
    # stopping): a fairer comparison (docs/FAIR_COMPARISON.md).
    known = ind & ((label == 0) if train_only else (label < 2)) & ~np.isin(y, test_diseases)
    train_diseases = pd.unique(y[known])
    train_diseases = np.random.default_rng(seed).choice(train_diseases, min(n_train, len(train_diseases)), replace=False)
    train_pos = {d: drug_pos[x[known & (y == d)]] for d in train_diseases}
    is_train = ind & (label == 0)

    F_test = score_paths(g, test_diseases, paths)
    rows_X, rows_y = [], []
    rng = np.random.default_rng(seed)
    gf = drug_degrees(g)
    for i in range(0, len(train_diseases), 256):
        block = train_diseases[i:i + 256]
        F_train = score_paths(g, block, paths)
        for c, d in enumerate(block):
            lab = np.zeros(len(g.drugs), bool)
            lab[train_pos[d]] = True
            pop = degree.copy()
            pop[drug_pos[x[is_train & (y == d)]]] -= 1  # its own indication does not count
            sel = lab | (rng.random(len(lab)) < 0.1)  # all the positives, 10 % of the negatives
            rows_X.append(np.column_stack([F_train[sel, c, :], pop[sel], gf[sel]]))
            rows_y.append(lab[sel])
    del F_train
    X, Y = np.concatenate(rows_X), np.concatenate(rows_y)
    nc = len(paths)
    graph_only, _ = fit_model(X[:, :nc], Y)
    full, clf = fit_model(X[:, :nc + 1], Y)
    gbm = fit_gbm(X, Y, seed)
    # Train the way it is tested: in each round the drug-disease edges of whole neighbourhoods
    # are hidden (diseases ≤ 2 steps away in the ontology or with near-identical genes) and these
    # diseases are added as examples. This way the model sees cases like the test ones: without
    # siblings with known drugs, and with drugs that have no other indication.
    if rounds:
        Xm, Ym, Gm = masked_rows(g, test_diseases, known, x, y, drug_pos, paths, rounds, seed)
        gbm_m = fit_gbm(np.concatenate([X, Xm]), np.concatenate([Y, Ym]), seed)
        if rank:
            groups = np.r_[[len(v) for v in rows_y], Gm]
            rank_m = fit_ranker(np.concatenate([X, Xm]), np.concatenate([Y, Ym]), groups, seed)

    def per_disease(f, cols):
        return np.column_stack([f(np.column_stack([F_test[:, c, :], degree, gf])[:, cols]) for c in range(len(test_diseases))])

    scores = {"popularity": np.repeat(degree[:, None], len(test_diseases), axis=1)}
    for k, c in enumerate(paths):
        scores[c] = F_test[:, :, k]
    scores["model: graph only"] = per_disease(graph_only, slice(0, nc))
    scores["model: graph + popularity"] = per_disease(full, slice(0, nc + 1))
    scores["gbm: graph + popularity + degrees"] = per_disease(gbm, slice(None))
    if rounds:
        scores["gbm + masked training"] = per_disease(gbm_m, slice(None))
        if rank:
            scores["ranking (LambdaRank) + masked"] = per_disease(rank_m, slice(None))
    tx = txgnn_predictions(split, seed, test_diseases, g.drugs)
    if tx is not None:
        scores["TxGNN (run here)"] = tx

    results = {k: metrics(v, test_pos, pool_neg) for k, v in scores.items()}
    # Scores for later analyses (evaluation/novelty.py).
    (CACHE / "scores").mkdir(parents=True, exist_ok=True)
    np.savez_compressed(CACHE / "scores" / f"{split}_{seed}_{feature_set}_r{rounds}{'_rank' if rank else ''}{'_trainonly' if train_only else ''}.npz",
                        diseases=test_diseases, drugs=g.drugs, degree=degree,
                        positives=np.array([np.asarray(p) for p in test_pos], dtype=object),
                        **{key: scores[k].astype(np.float32) for k, key in SCORE_KEYS.items() if k in scores})
    # The same metrics only with the test diseases that have no near-identical disease with
    # known drugs (ontology neighbour ≤ 2 steps away with < 50 % of their drugs).
    nd = near_duplicates(g, test_diseases, test_pos, drug_pos)
    hard = np.flatnonzero(nd < 0.5)
    for k in ["popularity", "model: graph + popularity", "gbm: graph + popularity + degrees", "gbm + masked training",
              "ranking (LambdaRank) + masked", "TxGNN (run here)"]:
        if k not in scores:
            continue
        results[f"{NO_NEAR_DUPLICATES}{k}"] = metrics(scores[k][:, hard], [test_pos[i] for i in hard], pool_neg)
    print(f"  {split} seed {seed}: {len(test_diseases)} test diseases ({len(hard)} without near duplicates), "
          f"{sum(map(len, test_pos))} indications; model trained on {len(train_diseases)} diseases ({time.time() - t0:.0f} s)")
    return results, {**dict(zip(paths + ["popularity"], clf.coef_[0])),
                     "_diseases": len(test_diseases), "_hard": len(hard)}


def txgnn_predictions(split, seed, diseases, drugs):
    """Predictions of TxGNN trained by evaluation/txgnn_wsl.py (drugs × diseases), or None.
    Checks that TxGNN has exactly the same test diseases."""
    f = CACHE / "txgnn_predictions" / f"TxGNN_{seed}_{split}.parquet"
    if not f.exists():
        return None
    d = pd.read_parquet(f)
    d = d[d.relation == "indication"]
    nodes = pd.read_csv(RAW / "node.csv", sep="	", dtype=str)
    idx = {(t, _id(i)): int(n) for n, i, t in zip(nodes.node_index, nodes.node_id, nodes.node_type)}
    d = d.assign(m=[idx.get(("disease", _id(i)), -1) for i in d.disease_id],
                 f=[idx.get(("drug", _id(i)), -1) for i in d.drug_id])
    theirs = set(d.m.unique())
    if theirs != set(diseases):
        print(f"  ⚠ TxGNN has {len(theirs)} test diseases and we have {len(set(diseases))} "
              f"({len(theirs & set(diseases))} in common): it is not the same split")
    pos_m = {m: j for j, m in enumerate(diseases)}
    pos_f = np.full(int(max(drugs.max(), d.f.max())) + 1, -1)
    pos_f[drugs] = np.arange(len(drugs))
    res = np.full((len(drugs), len(diseases)), -np.inf, np.float32)
    d = d[d.m.isin(list(pos_m)) & (d.f >= 0)]
    res[pos_f[d.f.to_numpy()], [pos_m[m] for m in d.m]] = d.score.to_numpy()
    res[np.isneginf(res)] = np.nanmin(res[np.isfinite(res)]) - 1
    return res


def published(split):
    if not PUBLISHED.exists():
        return None
    d = pd.read_csv(PUBLISHED)
    d = d[(d.Task == "indication") & (d.Split == split)]
    return d.groupby(["Method", "Metric Name"]).Metric.agg(["mean", "std"])


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--split", default="complex_disease", help="complex_disease, an area, or «areas»")
    p.add_argument("--seeds", default="1,2,3,4,5")
    p.add_argument("--n-train", type=int, default=1500)
    p.add_argument("--rounds", type=int, default=0, help="training rounds with masked neighbourhoods")
    p.add_argument("--rank", action="store_true", help="also try LightGBM with a ranking objective")
    p.add_argument("--train-only", action="store_true",
                   help="use only the training edges as labels (not the validation ones), like TxGNN")
    p.add_argument("--paths", choices=list(FEATURE_SETS), default="base",
                   help="base: PrimeKG; new: + new metapaths; extra: + OREVID biological sources")
    a = p.parse_args()
    kg = load_kg()
    types = node_types(kg)
    splits = AREAS if a.split == "areas" else [a.split]
    metr = ["AUPRC", "AUROC", "Recall@100", "Precision@K", "Avg Rank"]
    for split in splits:
        runs, coefs = [], []
        for s in [int(s) for s in a.seeds.split(",")]:
            r, c = evaluate(kg, types, split, s, a.n_train, a.paths, a.rounds, a.rank, a.train_only)
            runs.append(pd.DataFrame(r).T)
            coefs.append(c)
        mean = pd.concat(runs).groupby(level=0, sort=False).mean()
        print(f"\n== {split} (indication; mean of {len(runs)} seeds) ==")
        print(f"  {'':52s}" + "".join(f"{m:>12s}" for m in metr))
        pub = published(split)
        if pub is not None:
            for met in ["TxGNN", "RGCN", "HAN", "BioBERT"]:
                if met in pub.index.get_level_values(0):
                    print(f"  {'(published) ' + met:52s}" + "".join(f"{pub.loc[(met, m), 'mean']:12.3f}" for m in metr))
        for name, row in mean.iterrows():
            print(f"  {name:52s}" + "".join(f"{row[m]:12.3f}" for m in metr))
        print("  Weights of the full model (scaled features): " + ", ".join(
            f"{k} {np.mean([c[k] for c in coefs]):+.2f}" for k in coefs[0] if not k.startswith("_")))
        CACHE.mkdir(parents=True, exist_ok=True)
        pd.concat({s: t for s, t in zip(a.seeds.split(","), runs)}).to_csv(CACHE / f"results_{split}_{a.paths}{'_r' + str(a.rounds) if a.rounds else ''}.csv")


if __name__ == "__main__":
    main()
