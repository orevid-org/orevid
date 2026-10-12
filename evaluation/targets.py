"""Target prioritisation: for a disease, which genes/proteins to go after?

Evaluation: for each disease, the positives are the targets (ChEMBL mechanism of action) of the
drugs already approved for it (PrimeKG + Open Targets). Cross-validation over two halves of the
diseases: train on one and evaluate on the other.

Features per (disease, gene):
- evidence: PrimeKG gene-disease; Open Targets by type (genetic, somatic, animal model,
  expression, genetic literature); known genetic and GWAS+QTL direction;
- network: fraction of the gene's STRING neighbours associated with the disease; whether the gene is a target of drugs
  approved in diseases that share genes (without the same disease);
- gene: tractability (sources/targets.py), safety, LOEUF, membrane/secreted, probes, STRING
  degree, pleiotropy (associated diseases) and «target popularity» (drugs acting on it).

Variants: «biology» (without any feature that comes from drugs) and «full».
Open Targets literature (Europe PMC) is left out: it often mentions the target because there is already
a drug (leak).

    py -3 -m evaluation.targets
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from graph import Graph, INDICATION

ROOT = Path(__file__).resolve().parent.parent
OT_TYPES = ["genetic_association", "somatic_mutation", "animal_model", "rna_expression", "genetic_literature"]
GENE = ["tract_SM", "tract_AB", "tract_PR", "tract_OC", "loeuf", "n_safety", "membrane_or_secreted", "n_probes"]
# Essentiality (DepMap) and Open Targets prioritisation factors (sources/targets.py).
GENE_NEW = ["depmap_effect", "dependent_fraction", "essential", "isInMembrane", "isSecreted", "hasPocket", "hasLigand",
            "hasSmallMoleculeBinder", "geneticConstraint", "paralogMaxIdentityPercentage",
            "mouseOrthologMaxIdentityPercentage", "isCancerDriverGene", "mouseKOScore", "hasHighQualityChemicalProbes",
            "tissueSpecificity", "tissueDistribution", "celltypeSpecificity", "celltypeDistribution"]
FROM_DRUGS = {"target in similar diseases", "target in neighbour diseases", "hasLigand", "hasSmallMoleculeBinder",
              "n_probes", "n_safety", "hasHighQualityChemicalProbes", "tract_SM", "tract_AB", "tract_PR", "tract_OC", "target popularity"}


def score_matrix(a, rel, disp, N, node_types):
    """Disease × gene with the Open Targets score (oriented: rows = disease)."""
    s = a[(a.relation == rel) & (a.display_relation == disp)]
    x, y, v = s.x_index.to_numpy(), s.y_index.to_numpy(), s.score.fillna(1).to_numpy()
    flip = node_types[x] != "disease"
    x, y = np.where(flip, y, x), np.where(flip, x, y)
    m = sp.csr_matrix((v, (x, y)), shape=(N, N))
    m.sum_duplicates()
    return m


NEW_FEATURES = True  # essentiality and prioritisation features (False: as before, to compare)


class Targets:
    def __init__(self):
        t0 = time.time()
        self.g = g = Graph()
        self.N, node_types = g.N, g.types
        a = g.attributes()
        self.ot = {k: score_matrix(a, "ot_assoc", k, self.N, node_types) for k in OT_TYPES}
        self.pkg = g.M("disease_protein").astype(np.float64).tocsr()
        self.dir = {}
        for r in ["dir_genetic", "dir_gwas"]:
            m = score_matrix(a, r, "activate", self.N, node_types) + score_matrix(a, r, "inhibit", self.N, node_types)
            m.data[:] = 1
            self.dir[r] = m
        # Total evidence (without literature): defines the associated genes and the candidates.
        E = self.pkg.copy()
        for m in self.ot.values():
            E = E + (m > 0).astype(np.float64)
        E = E + self.dir["dir_gwas"]
        self.E = (E > 0).astype(np.float64).tocsr()
        # Approved drugs per disease (PrimeKG + Open Targets) and targets (ChEMBL).
        ind = g.M(INDICATION).astype(np.float64)
        ap = a[(a.relation == "clinical") & (a.display_relation == "approved")]
        A = sp.csr_matrix((np.ones(len(ap)), (ap.y_index, ap.x_index)), shape=(self.N, self.N))
        self.approved = ((ind + A) > 0).astype(np.float64).tocsr()  # disease × drug (and drug × disease)
        moa = a[a.relation == "moa"]
        self.moa = sp.csr_matrix((np.ones(len(moa)), (moa.x_index, moa.y_index)), shape=(self.N, self.N))
        self.moa.data[:] = 1
        self.action = moa.groupby(["x_index", "y_index"]).display_relation.first().to_dict()
        self.moa_c = self.moa.tocsc()
        self.T = ((self.approved @ self.moa) > 0).astype(np.float64).tocsr()  # disease × target gene
        self.Tc = self.T.tocsc()
        O = g.M("disease_disease").astype(np.float64)
        self.O2 = ((O + O @ O) > 0).astype(np.float64).tocsr()
        self.string = g.M("string").astype(np.float64).tocsr() if "string" in g.relations else g.M("protein_protein")
        self.string_degree = np.asarray(self.string.sum(axis=1)).ravel()
        self.pleio = np.asarray(self.E.sum(axis=0)).ravel()
        self.target_pop = np.asarray(self.moa.sum(axis=0)).ravel()
        d = pd.read_parquet(ROOT / "data" / "ext" / "targets.parquet").set_index("node_index")
        self.gene_columns = GENE + [c for c in GENE_NEW if c in d.columns] if NEW_FEATURES else GENE
        self.gene = pd.DataFrame(index=np.arange(self.N), columns=self.gene_columns, dtype=float)
        self.gene.loc[d.index, self.gene_columns] = d[self.gene_columns].to_numpy(float)
        self.annot = d
        self.genes = np.flatnonzero(node_types == "gene/protein")
        # Against leaks. The drugs of the same disease are NOT subtracted from the target popularity:
        # the model recognises each gene by its fixed traits (degree, pleiotropy…) and learns «lower
        # than usual ⇒ positive» (R@20 0.32 -> 0.62, false). «How many drugs act on it» is a fact of the
        # gene that does not depend on the disease. The «clinical» level of tractability is removed, though.
        self.subtract_own, self.cap = False, True
        self.biology_only = False
        print(f"target data loaded ({time.time() - t0:.0f} s)", flush=True)

    def similar(self, d):
        """Weight of each disease by gene similarity (without itself) and ontology neighbours (≤ 2 steps)."""
        sim = (self.E @ self.E[d].T).toarray().ravel()
        sim[d] = 0
        sim = sim / np.maximum(np.asarray(self.E.sum(axis=1)).ravel(), 1)
        ont = self.O2[d].toarray().ravel()
        ont[d] = 0
        return sim, ont

    def candidates(self, d, sim, ont):
        """Genes with evidence, their STRING neighbours and the targets of drugs approved in similar
        or neighbour diseases."""
        e = self.E[d].indices
        neighbors = self.string[e].indices if len(e) else np.array([], int)
        top = np.argsort(-sim)[:50]
        alt = np.flatnonzero(ont > 0)
        if self.biology_only:  # without candidates coming from drug targets (independent review)
            return np.unique(np.r_[e, neighbors]).astype(int)
        targets = self.T[np.r_[top[sim[top] > 0], alt]].indices
        return np.unique(np.r_[e, neighbors, targets]).astype(int)

    def features(self, d):
        """(candidate genes, feature matrix, names) for a disease."""
        sim, ont = self.similar(d)
        c = self.candidates(d, sim, ont)
        if not len(c):
            return c, np.zeros((0, 0)), []
        cols, names = [], []
        cols.append(self.E[d, c].toarray().ravel()); names.append("has evidence")
        cols.append(self.pkg[d, c].toarray().ravel()); names.append("PrimeKG gene-disease")
        for k, m in self.ot.items():
            cols.append(m[d, c].toarray().ravel()); names.append(f"OT {k}")
        for k, m in self.dir.items():
            cols.append(m[d, c].toarray().ravel()); names.append(k)
        ev = self.E[d].toarray().ravel()
        neighbors = self.string[c] @ ev
        cols.append(neighbors / np.maximum(self.string_degree[c], 1)); names.append("associated STRING neighbours")
        cols.append(self.Tc[:, c].T @ sim); names.append("target in similar diseases")
        cols.append(self.Tc[:, c].T @ ont); names.append("target in neighbour diseases")
        # Drugs acting on the gene (see subtract_own in __init__).
        own = (self.approved[d] @ self.moa[:, c]).toarray().ravel() if self.approved[d].nnz and self.subtract_own else 0
        cols.append(self.target_pop[c] - own); names.append("target popularity")
        cols.append(self.string_degree[c]); names.append("STRING degree")
        cols.append(self.pleio[c]); names.append("pleiotropy")
        G = self.gene.loc[c, self.gene_columns].to_numpy(float).copy()
        # Tractability without the «clinical» level: it could come from the drugs of this same disease.
        if self.cap:
            G[:, :4] = np.minimum(G[:, :4], 2)
        for j, k in enumerate(self.gene_columns):
            cols.append(G[:, j]); names.append(k)
        return c, np.column_stack(cols), names


def recall(s, cand, pos, k, n_total, rng):
    """Recall@k and rank of the positives; positives that are not candidates go to the end."""
    v = s + rng.random(len(s)) * 1e-9
    order = cand[np.argsort(-v)]
    rank = {g: i for i, g in enumerate(order)}
    r = np.array([rank.get(p, n_total) for p in pos])
    return {"Recall@20": np.mean(r < 20), "Recall@100": np.mean(r < 100), "Precision@K": np.mean(r < len(pos)),
            "coverage": np.mean(r < n_total)}


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    strict = "--strict" in sys.argv
    from sklearn.ensemble import HistGradientBoostingClassifier
    D = Targets()
    D.biology_only = strict
    node_types = D.g.types
    diseases = [d for d in np.flatnonzero(node_types == "disease") if D.T[d].nnz and D.E[d].nnz]
    print(f"{len(diseases)} diseases with approved drugs with a known target and associated genes", flush=True)
    rng = np.random.default_rng(0)
    rng.shuffle(diseases)
    halves = [diseases[: len(diseases) // 2], diseases[len(diseases) // 2:]]
    rows = {}
    for d in diseases:
        c, X, names = D.features(d)
        rows[d] = (c, X, D.T[d].indices)
    i_bio = [j for j, n in enumerate(names) if n not in FROM_DRUGS]
    res = []
    n_total = len(D.genes)
    for k_train, k_test in [(0, 1), (1, 0)]:
        Xs, Ys = [], []
        for d in halves[k_train]:
            c, X, pos = rows[d]
            y = np.isin(c, pos)
            sel = y | (rng.random(len(y)) < 0.2)
            Xs.append(X[sel]); Ys.append(y[sel])
        X, Y = np.concatenate(Xs), np.concatenate(Ys)
        models = {}
        without = lambda *drop: [j for j, n in enumerate(names) if n not in drop]
        variants = [("Model: biology only", i_bio), ("Model: full", list(range(X.shape[1]))),
                    ("Full without target popularity", without("target popularity")),
                    ("Full without neighbour diseases", without("target in neighbour diseases")),
                    ("Full without neighbour or similar", without("target in neighbour diseases", "target in similar diseases"))]
        for name, cols in variants:
            clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31,
                                                 min_samples_leaf=50, l2_regularization=1.0, random_state=0)
            clf.fit(X[:, cols], Y)
            models[name] = (clf, cols)
        for d in halves[k_test]:
            c, X, pos = rows[d]
            if not len(c):
                continue
            j = {n: i for i, n in enumerate(names)}
            scores = {
                "Target popularity": X[:, j["target popularity"]],
                "Genetics only (OT)": X[:, j["OT genetic_association"]],
                "Total evidence (sum)": X[:, [j["PrimeKG gene-disease"]] + [j[f"OT {t}"] for t in OT_TYPES]].sum(axis=1),
            }
            for name, (clf, cols) in models.items():
                scores[name] = clf.predict_proba(X[:, cols])[:, 1]
            for name, s in scores.items():
                res.append({"method": name, "disease": d, "n_pos": len(pos), "candidates": len(c),
                            **recall(s, c, pos, 20, n_total, rng)})
    r = pd.DataFrame(res)
    table = r.groupby("method")[["Recall@20", "Recall@100", "Precision@K", "coverage"]].mean()
    order = ["Target popularity", "Genetics only (OT)", "Total evidence (sum)", "Model: biology only",
             "Full without neighbour or similar", "Full without neighbour diseases", "Full without target popularity",
             "Model: full"]
    table = table.reindex(order)
    print(f"\nDiseases evaluated: {r.disease.nunique()}; positive targets per disease (mean): "
          f"{r.groupby('disease').n_pos.first().mean():.1f}; candidate genes (mean): {r.groupby('disease').candidates.first().mean():.0f}")
    print("Coverage: fraction of positive targets that are among the candidates (genes with some evidence).")
    print(table.round(3).to_string())
    (ROOT / "data" / "eval").mkdir(exist_ok=True)
    r.to_csv(ROOT / "data" / "eval" / "targets_results.csv", index=False)
    return table


if __name__ == "__main__" and (len(sys.argv) == 1 or sys.argv[1] == "--strict"):
    main()


# --- production: command «explore.py targets» ------------------------------------------------

TARGETS_MODEL = ROOT / "data" / "prep" / "targets_model.pkl"


def train_production(seed=0):
    """«Full» model trained with all the diseases with approved targets."""
    import pickle
    from sklearn.ensemble import HistGradientBoostingClassifier
    D = Targets()
    diseases = [d for d in np.flatnonzero(D.g.types == "disease") if D.T[d].nnz and D.E[d].nnz]
    rng = np.random.default_rng(seed)
    Xs, Ys = [], []
    for d in diseases:
        c, X, names = D.features(d)
        y = np.isin(c, D.T[d].indices)
        sel = y | (rng.random(len(y)) < 0.2)
        Xs.append(X[sel]); Ys.append(y[sel])
    clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05, max_leaf_nodes=31, min_samples_leaf=50,
                                         l2_regularization=1.0, random_state=seed).fit(np.concatenate(Xs), np.concatenate(Ys))
    TARGETS_MODEL.write_bytes(pickle.dumps({"model": clf, "names": names}))
    print(f"target model trained with {len(diseases)} diseases -> {TARGETS_MODEL}")


LEVELS = {3: "clinical", 2: "evidence", 1: "predicted"}


def prioritise(D, d, top=20, with_evidence=False):
    """Table of candidate targets for a disease, with how to go after them."""
    import pickle
    m = pickle.loads(TARGETS_MODEL.read_bytes())  # our own file (train_production)
    c, X, names = D.features(d)
    p = m["model"].predict_proba(X)[:, 1]
    j = {n: i for i, n in enumerate(names)}
    a = D.g.attributes()
    dirs = a[(a.relation.isin(["dir_genetic", "dir_gwas"])) & ((a.x_index == d) | (a.y_index == d))]
    direction = {}
    for r in dirs.itertuples():
        gene = r.y_index if r.x_index == d else r.x_index
        direction.setdefault(gene, set()).add(f"{r.display_relation} ({'genetic' if r.relation == 'dir_genetic' else 'GWAS+QTL'})")
    already = set(D.T[d].indices)
    order = np.argsort(-p)
    new = [k for k in order if c[k] not in already and (not with_evidence or X[k, j["has evidence"]] > 0)][:top]
    exploited = [k for k in order if c[k] in already][:top]
    rows = []
    for k in exploited + new:
        gene = c[k]
        ev = [n.replace("OT ", "") for n in ["PrimeKG gene-disease"] + [f"OT {t}" for t in OT_TYPES] if X[k, j[n]] > 0]
        an = D.annot.loc[gene] if gene in D.annot.index else None
        modal = []
        if an is not None:
            for mod, name in [("SM", "small molecule"), ("AB", "antibody"), ("PR", "PROTAC"), ("OC", "oligonucleotide")]:
                if an[f"tract_{mod}"] > 0:
                    modal.append(f"{name} ({LEVELS[int(an[f'tract_{mod}'])]})")
        acting = D.moa_c[:, gene].indices
        drugs = [D.g.name(f) for f in acting[:6]]
        actions = sorted({D.action.get((f, gene), "") for f in acting})
        rows.append({
            "gene": D.g.name(gene), "score": round(float(p[k]), 3),
            "evidence": ", ".join(ev) if ev else ("network neighbour" if X[k, j["associated STRING neighbours"]] > 0 else "target in similar diseases"),
            "direction": "; ".join(sorted(direction.get(gene, []))) or "unknown",
            "modalities": ", ".join(modal) or "none known",
            "safety": (an["safety"] if an is not None and an["n_safety"] > 0 else ""),
            "LOEUF": (round(float(an["loeuf"]), 2) if an is not None and an["loeuf"] == an["loeuf"] else None),
            "drugs acting on it": (", ".join(drugs) + (f" (+{len(acting) - 6})" if len(acting) > 6 else "")
                                   + (f" [{'/'.join(a_ for a_ in actions if a_)}]" if any(actions) else "")),
            "uniprot": (an["uniprot"] if an is not None and "uniprot" in an else ""),
            "essential": (int(an["essential"]) if an is not None and an.get("essential") == an.get("essential") else None),
            "dependent_fraction": (float(an["dependent_fraction"]) if an is not None and an.get("dependent_fraction") == an.get("dependent_fraction") else None),
            "tissue_specificity": (float(an["tissueSpecificity"]) if an is not None and an.get("tissueSpecificity") == an.get("tissueSpecificity") else None),
            "already_target_here": gene in already,
            "has_evidence": bool(X[k, j["has evidence"]] > 0),
        })
    return pd.DataFrame(rows)


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "train":
    train_production()
