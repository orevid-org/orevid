"""Direction of the effect: does a drug improve or worsen a disease?

Classifier «indication (improves) versus contraindication (worsens)» for (drug, disease)
pairs, with direction signals. Evaluation: PrimeKG indications and contraindications,
cross-validation over two halves of the diseases.

Signals (none uses the relation of the pair itself):
- genetics: targets where the drug mechanism (ChEMBL: activates/inhibits) agrees or not with
  what human genetics asks for the disease (dir_genetic, dir_gwas);
- analogues: drugs approved in neighbour diseases (ontology ≤ 2 steps) or similar ones (shared
  genes), without the same disease, that act on the same targets in the same direction or
  in the opposite one;
- adverse effects: phenotypes of the disease that the drug causes (SIDER, FAERS);
- neighbours: neighbour/similar diseases where the drug is indicated or contraindicated.

Drug counts «minus this pair» are NOT used (see docs/TARGETS.md: the model
recognises the drug and learns «lower than usual»).

    py -3 -m evaluation.direction          # evaluation
    py -3 -m evaluation.direction train    # production model (data/prep/direction_model.pkl)
"""
import pickle
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from graph import Graph, INDICATION, CONTRA

ROOT = Path(__file__).resolve().parent.parent
MODEL = ROOT / "data" / "prep" / "direction_model.pkl"
GROUPS = {
    "genetics": ["gene for", "gene against"],
    "analogues": ["analogues for", "analogues against"],
    "adverse effects": ["SIDER phenotypes", "SIDER fraction", "FAERS phenotypes", "FAERS fraction"],
    "neighbours": ["indicated in neighbours", "contraindicated in neighbours", "indicated in similar", "contraindicated in similar"],
}


def _directed(a, rel, N):
    """Disease × gene matrix: +1 activate, -1 inhibit."""
    s = a[a.relation == rel]
    v = np.where(s.display_relation == "activate", 1.0, -1.0)
    m = sp.csr_matrix((v, (s.x_index, s.y_index)), shape=(N, N))
    m.sum_duplicates()
    m.data = np.sign(m.data)
    return m


class Direction:
    def __init__(self):
        t0 = time.time()
        self.g = g = Graph()
        N = g.N
        a = g.attributes()
        moa = a[(a.relation == "moa") & a.display_relation.isin(["activates", "inhibits"])]
        v = np.where(moa.display_relation == "activates", 1.0, -1.0)
        self.moa = sp.csr_matrix((v, (moa.x_index, moa.y_index)), shape=(N, N))  # drug × gene ±1
        self.moa.sum_duplicates()
        self.moa.data = np.sign(self.moa.data)
        self.dir = (_directed(a, "dir_genetic", N) + _directed(a, "dir_gwas", N)).tocsr()
        self.dir.data = np.sign(self.dir.data)
        self.ind = g.M(INDICATION).astype(np.float64).tocsr()
        ap = a[(a.relation == "clinical") & (a.display_relation == "approved")]
        A = sp.csr_matrix((np.ones(len(ap)), (ap.y_index, ap.x_index)), shape=(N, N))
        self.approved = ((self.ind + A + A.T) > 0).astype(np.float64).tocsr()  # symmetric drug ↔ disease
        self.contra = g.M(CONTRA).astype(np.float64).tocsr()
        O = g.M("disease_disease").astype(np.float64)
        self.O2 = ((O + O @ O) > 0).astype(np.float64).tocsr()
        G = g.M("disease_protein").astype(np.float64).tocsr()
        self.G = G
        self.phen = g.M("disease_phenotype_positive").astype(np.float64).tocsr()
        self.sider = g.M("drug_effect").astype(np.float64).tocsr()
        self.faers = g.M("faers").astype(np.float64).tocsr() if "faers" in g.relations else self.sider * 0
        print(f"direction data loaded ({time.time() - t0:.0f} s)", flush=True)

    def neighbour_diseases(self, d):
        """(ontology neighbours ≤ 2 steps, 30 most similar by genes), without itself."""
        o = self.O2[d].indices
        o = o[o != d]
        sim = (self.G @ self.G[d].T).toarray().ravel()
        sim[d] = 0
        s = np.argsort(-sim)[:30]
        return o, s[sim[s] > 0]

    def features(self, d, drugs):
        """Matrix (drugs × features) for a disease."""
        drugs = np.asarray(drugs)
        o, s = self.neighbour_diseases(d)
        X = {}
        # genetics: drug targets where genetics says activate/inhibit
        needed = self.dir[d]
        M = self.moa[drugs]
        prod = (M.multiply(needed)).tocsr()  # +1 agrees, -1 opposite
        X["gene for"] = np.asarray((prod > 0).sum(axis=1)).ravel()
        X["gene against"] = np.asarray((prod < 0).sum(axis=1)).ravel()
        # analogues: action of the drugs approved in neighbours/similar on each gene
        neighbours = np.unique(np.r_[o, s]).astype(int)
        if len(neighbours):
            analogues = np.unique(self.approved[neighbours].indices)
            profile = np.asarray(self.moa[analogues].sum(axis=0)).ravel()  # >0 activated, <0 inhibited
            profile = np.sign(profile)
            P = M.multiply(profile[None, :]).tocsr()
            X["analogues for"] = np.asarray((P > 0).sum(axis=1)).ravel()
            X["analogues against"] = np.asarray((P < 0).sum(axis=1)).ravel()
        else:
            X["analogues for"] = X["analogues against"] = np.zeros(len(drugs))
        # adverse effects that are phenotypes of the disease
        fd = self.phen[d]
        nf = max(fd.nnz, 1)
        for name, m in [("SIDER", self.sider), ("FAERS", self.faers)]:
            c = np.asarray((m[drugs].multiply(fd)).sum(axis=1)).ravel()
            X[f"{name} phenotypes"] = c
            X[f"{name} fraction"] = c / nf
        # neighbours: indicated / contraindicated in neighbour or similar diseases
        for name, idx in [("neighbours", o), ("similar", s)]:
            if len(idx):
                X[f"indicated in {name}"] = np.asarray(self.approved[drugs][:, idx].sum(axis=1)).ravel()
                X[f"contraindicated in {name}"] = np.asarray(self.contra[drugs][:, idx].sum(axis=1)).ravel()
            else:
                X[f"indicated in {name}"] = X[f"contraindicated in {name}"] = np.zeros(len(drugs))
        names = [n for group in GROUPS.values() for n in group]
        return np.column_stack([X[n] for n in names]).astype(np.float64), names


def pairs(D):
    """(disease, drug, label): 1 indication, 0 contraindication (without those that are both)."""
    ind = sp.triu(D.ind).tocoo() if False else D.ind.tocoo()
    con = D.contra.tocoo()
    types = D.g.types
    pi = {(r, c) for r, c in zip(ind.row, ind.col) if types[r] == "disease" and types[c] == "drug"}
    pc = {(r, c) for r, c in zip(con.row, con.col) if types[r] == "disease" and types[c] == "drug"}
    both = pi & pc
    return pd.DataFrame([(d, f, 1) for d, f in pi - both] + [(d, f, 0) for d, f in pc - both],
                        columns=["disease", "drug", "y"])


def table(D, p):
    rows = []
    for d, s in p.groupby("disease"):
        X, names = D.features(d, s.drug.to_numpy())
        rows.append(pd.DataFrame(X, columns=names).assign(disease=d, drug=s.drug.to_numpy(), y=s.y.to_numpy()))
    return pd.concat(rows, ignore_index=True), names


def gbm(seed=0):
    from sklearn.ensemble import HistGradientBoostingClassifier
    return HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=30,
                                          l2_regularization=1.0, random_state=seed)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from sklearn.metrics import roc_auc_score, average_precision_score
    D = Direction()
    p = pairs(D)
    t, names = table(D, p)
    ind_degree = np.diff(D.ind.indptr)
    t["few uses"] = ind_degree[t.drug] <= 1
    print(f"{len(t):,} pairs ({t.y.mean():.1%} indications) in {t.disease.nunique()} diseases", flush=True)
    rng = np.random.default_rng(0)
    m = t.disease.unique()
    rng.shuffle(m)
    half = set(m[: len(m) // 2])
    fold = t.disease.isin(half).to_numpy()
    variants = {"genetics": GROUPS["genetics"], "analogues": GROUPS["analogues"],
                "adverse effects": GROUPS["adverse effects"],
                "mechanism (genetics + analogues)": GROUPS["genetics"] + GROUPS["analogues"],
                "mechanism + adverse effects": GROUPS["genetics"] + GROUPS["analogues"] + GROUPS["adverse effects"],
                "all (+ neighbours)": names}
    res = []
    for name, cols in variants.items():
        pred = np.zeros(len(t))
        for train_, test in [(fold, ~fold), (~fold, fold)]:
            clf = gbm().fit(t.loc[train_, cols], t.y[train_])
            pred[test] = clf.predict_proba(t.loc[test, cols])[:, 1]
        for sub, mask in [("all", np.ones(len(t), bool)), ("drugs with ≤ 1 indication", t["few uses"].to_numpy())]:
            y = t.y.to_numpy()[mask]
            res.append({"features": name, "subset": sub, "AUROC": roc_auc_score(y, pred[mask]),
                        "AUPRC (contraindication)": average_precision_score(1 - y, 1 - pred[mask]),
                        "n": int(mask.sum())})
    r = pd.DataFrame(res)
    print(r.pivot(index="features", columns="subset", values="AUROC").reindex(list(variants)).round(3).to_string())
    print()
    print(r.round(3).to_string(index=False))
    r.to_csv(ROOT / "data" / "eval" / "direction_results.csv", index=False)
    coverage = {g_: float((t[c].sum(axis=1) > 0).mean()) for g_, c in GROUPS.items()}
    print("\nFraction of pairs with some signal per group:", {k: round(v, 2) for k, v in coverage.items()})


def train():
    sys.stdout.reconfigure(encoding="utf-8")
    D = Direction()
    t, names = table(D, pairs(D))
    cols = GROUPS["genetics"] + GROUPS["analogues"] + GROUPS["adverse effects"] + GROUPS["neighbours"]
    clf = gbm().fit(t[cols], t.y)
    MODEL.write_bytes(pickle.dumps({"model": clf, "cols": cols}))
    print(f"direction model trained with {len(t):,} pairs -> {MODEL}")


def improvement_probability(D, d, drugs):
    """P(indication | indication or contraindication) for each drug, and the features."""
    m = pickle.loads(MODEL.read_bytes())  # our own file (train)
    X, names = D.features(d, drugs)
    Xc = pd.DataFrame(X, columns=names)[m["cols"]]
    return m["model"].predict_proba(Xc)[:, 1], pd.DataFrame(X, columns=names)


if __name__ == "__main__":
    {"train": train}.get(sys.argv[1] if len(sys.argv) > 1 else "", main)()


# --- explicit safety rule (not learned) ----------------------------------------------
# Drugs that block immune checkpoints activate the immune system and can trigger or
# worsen autoimmune diseases. The data (registered contraindications) do not cover this
# for new drugs, so it is flagged with a rule. Agonists of these targets
# (in development for autoimmune diseases) are left out.
CHECKPOINTS = {"PDCD1", "CD274", "PDCD1LG2", "CTLA4", "LAG3", "HAVCR2", "TIGIT"}
IMMUNE = "MONDO_0005046"  # immune system disorder (Open Targets therapeutic area)
AUTOIMMUNE = "MONDO_0007179"   # autoimmune disease


def checkpoint_rule(g):
    """(drugs that block checkpoints, immune diseases) as sets of nodes."""
    from sources.ids import OT
    a = g.attributes()
    genes = {i for i in np.flatnonzero(g.types == "gene/protein") if g.name(i) in CHECKPOINTS}
    targets = g.M("drug_protein|moa").tocsc()
    drugs = set()
    for gene in genes:
        drugs.update(targets.indices[targets.indptr[gene]:targets.indptr[gene + 1]])
    activates = a[(a.relation == "moa") & (a.display_relation == "activates") & a.y_index.isin(list(genes))]
    drugs -= set(activates.x_index)
    d = pd.read_parquet(OT / "disease", columns=["id", "therapeuticAreas", "ancestors"])
    # Only autoimmune diseases. The «immune system disorder» area includes lymphomas, leukaemias and
    # thymomas, where these drugs are indeed used (validation: 196 trials; error corrected).
    imm = d[[(AUTOIMMUNE in (an if an is not None else [])) or i == AUTOIMMUNE
             for i, an in zip(d.id, d.ancestors)]]
    mondo = {int(i.split("_")[1]) for i in imm.id if i.startswith("MONDO_")}
    n = g.nodes
    diseases = set()
    for i, t, nid in zip(n.node_index, n.node_type, n.node_id.astype(str)):
        if t == "disease" and not nid.startswith(("OT:", "MONDO:")):
            if any(p.isdigit() and int(p) in mondo for p in nid.split("_")):
                diseases.add(i)
    return drugs, diseases


def risk_flags(D, d, drugs, blockers, autoimmune):
    """Production risk flag (transparent rule, docs/DIRECTION.md): contraindicated in more
    neighbour or similar diseases than indicated (97 % of real contraindications; baseline 77 %),
    or a checkpoint blocker in an autoimmune disease.
    Returns (flag, contraindicated, indicated, checkpoint), arrays over `drugs`."""
    X, names = D.features(d, drugs)
    v = dict(zip(names, X.T))
    n_contra = v["contraindicated in neighbours"] + v["contraindicated in similar"]
    n_ind = v["indicated in neighbours"] + v["indicated in similar"]
    pc = np.isin(drugs, list(blockers)) & (d in autoimmune)
    return (n_contra > n_ind) | pc, n_contra, n_ind, pc
