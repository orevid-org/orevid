"""Access to the prepared graph and scoring of candidate drugs for a disease.

The score follows the Hetionet idea (DWPC): count the paths of a given
type ("metapath") between the drug and the disease, but each intermediate node weighs
less the more connected it is. That way a hub gene, which touches half the graph, does not dominate.

A relation can be written as a union: "drug_protein:target|moa" means the
DrugBank targets (PrimeKG) plus the ChEMBL ones (Open Targets).
"""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

PREP = Path(__file__).parent / "data" / "prep"

# Relations
TARGET_PKG = "drug_protein:target"
TARGET = "drug_protein:target|moa"
DISEASE_GENE = "disease_protein"
GENES_GENETIC = "ot_assoc:genetic_association|ot_assoc:genetic_literature|ot_assoc:somatic_mutation"
GENES_ANIMAL = "ot_assoc:animal_model"
GENES_LITERATURE = "ot_assoc:literature"
# Without "affected_pathway": it assigns identical blocks of hundreds of genes to different diseases.
GENES_ALL = f"{DISEASE_GENE}|{GENES_GENETIC}|{GENES_ANIMAL}|ot_assoc:rna_expression"
PPI = "protein_protein"
INDICATION = "indication"
OFF_LABEL = "off-label use"
CONTRA = "contraindication"
ONTOLOGY = "disease_disease"
EFFECT = "drug_effect"
CLINICAL = "clinical"
CLINICAL_ADVANCED = "clinical:approved|clinical:phase 3"
NEGATIVE = "clinical_negative"
FAERS = "faers"

# Effects of the drug on a gene: (relations that activate it, relations that inhibit it)
EFFECTS = {
    "moa": ("moa:activates", "moa:inhibits"),
    "all": ("moa:activates|ctd_gene:increases", "moa:inhibits|ctd_gene:decreases"),
    "all+lit": ("moa:activates|ctd_gene:increases|pubtator_drug_gene:increases",
                "moa:inhibits|ctd_gene:decreases|pubtator_drug_gene:decreases"),
}
DIR_ALL = "dir_genetic|dir_animal|dir_expression"

# name: (kind, drug-gene or drug-disease relation, disease-gene relation, description)
METAPATHS = {
    "target": ("direct", TARGET_PKG, DISEASE_GENE, "drug → target (DrugBank) = disease gene (PrimeKG)"),
    "interaction": ("ppi", TARGET_PKG, DISEASE_GENE, "drug → target (DrugBank) ↔ disease gene (PrimeKG)"),
    "similar": ("similar", INDICATION, DISEASE_GENE, "indicated for a disease that shares genes"),
    "ontology": ("ontology", INDICATION, None, "indicated for a neighbouring disease in the ontology"),
    # new layers (Open Targets)
    "target+": ("direct", TARGET, DISEASE_GENE, "DrugBank + ChEMBL targets = disease gene (PrimeKG)"),
    "target genetics": ("direct", TARGET, GENES_GENETIC, "target = gene with genetic evidence (GWAS, Mendelian, somatic)"),
    "interaction genetics": ("ppi", TARGET, GENES_GENETIC, "target ↔ gene with genetic evidence"),
    "target animal": ("direct", TARGET, GENES_ANIMAL, "target = gene with an animal model of the disease"),
    "target literature": ("direct", TARGET, GENES_LITERATURE, "target = gene associated in the literature (leakage risk)"),
    "target all": ("direct", TARGET, GENES_ALL, "target = gene associated by any evidence (without literature)"),
    "interaction all": ("ppi", TARGET, GENES_ALL, "target ↔ gene associated by any evidence"),
    "similar genetics": ("similar", INDICATION, GENES_GENETIC, "indicated for a disease with genetic genes in common"),
    "similar all": ("similar", INDICATION, GENES_ALL, "indicated for a disease with genes in common (all layers)"),
    "similar clinical": ("similar", CLINICAL_ADVANCED, GENES_ALL, "in phase 3 or approved for a disease with genes in common"),
    # STRING interactions instead of the PrimeKG ones (5th element: gene-gene relation)
    "interaction STRING": ("ppi", TARGET, GENES_ALL, "target ↔ (STRING ≥ 700) associated gene", "string"),
    "interaction STRING 900": ("ppi", TARGET, GENES_ALL, "target ↔ (STRING ≥ 900) associated gene", "string:≥900"),
    "interaction gen. STRING": ("ppi", TARGET, GENES_GENETIC, "target ↔ (STRING ≥ 900) genetic gene", "string:≥900"),
    # CTD: the drug changes the expression or activity of the gene (curated literature)
    "effect CTD": ("direct", "ctd_gene:increases|ctd_gene:decreases", GENES_ALL,
                   "the drug increases or decreases a disease gene (CTD)"),
    "effect CTD genetic": ("direct", "ctd_gene:increases|ctd_gene:decreases", GENES_GENETIC,
                           "the drug increases or decreases a genetic disease gene (CTD)"),
    "failed similar": ("similar", NEGATIVE, GENES_ALL, "trials stopped in diseases that share genes"),
    # Direction: the drug does to the gene what the disease needs (✓) or the opposite (✗).
    # 2nd element: drug effects (see EFFECTS); 3rd: what the disease needs (sources/direction.py).
    # QTL colocalisation direction (M1, P58; in the model by the author's decision, P61): targets the drug moves the
    # way human genetics asks minus those it moves the other way; 0 without a directed target
    "qtl direction": ("qtl_direction", None, None, "moves its targets in the direction human genetics (QTL colocalisation) asks"),
    # Orphanet genes (L3b, P79; in the model by the author's decision): log(1 + the drug's targets that are Orphanet genes
    # of the disease), from data/ext/orphadata_genes.parquet; no drug–disease edge
    "orphanet genes": ("orphanet_genes", "drug_protein|moa", None, "acts on a gene that Orphanet lists for the disease"),
    "dir genetics ✓": ("agrees", "all", "dir_genetic", "activates/inhibits the gene as genetics asks"),
    "dir genetics ✗": ("contradicts", "all", "dir_genetic", "does to the gene the opposite of what genetics asks"),
    "dir all ✓": ("agrees", "all", DIR_ALL, "does to the gene what genetics, mouse or expression ask"),
    "dir all ✗": ("contradicts", "all", DIR_ALL, "does to the gene the opposite of what genetics, mouse or expression ask"),
    "dir ChEMBL ✓": ("agrees", "moa", DIR_ALL, "ChEMBL mechanism in the direction the disease asks"),
    "dir ChEMBL ✗": ("contradicts", "moa", DIR_ALL, "ChEMBL mechanism in the opposite direction"),
    "dir lit ✓": ("agrees", "all+lit", DIR_ALL, "like \"dir all ✓\", with the PubTator effects"),
    "dir lit ✗": ("contradicts", "all+lit", DIR_ALL, "like \"dir all ✗\", with the PubTator effects"),
    "similar CTD": ("similar", "ctd_disease:therapeutic", GENES_ALL,
                    "therapeutic (CTD) in a disease that shares genes"),
}
# Metapaths of the model used by the "disease" command (see the README). Without
# "target literature" (it can leak drug knowledge) nor "ontology" (it copies
# near-identical diseases) nor popularity (it always pushes the same drugs).
COMBINED = ["target+", "target all", "target animal", "interaction", "interaction STRING",
            "similar all", "similar clinical", "effect CTD"]
MODEL = PREP / "model.json"
# Penalty (in logit) for the fraction of failed trials in the same disease. In Open Targets,
# with the same number of trials, having some stopped divides the probability of approval
# by 2-3, and if more than half were stopped, none was approved (README).
WEIGHT_FAILED = 3.0
# Direction adjustment (in logit), with the ChEMBL mechanism and the human genetic direction.
# Of the drug-disease pairs tested in the clinic where the gene must be activated, activators
# are approved 44 % of the time and inhibitors 28 % (log OR 0.71). But in general activators
# are already approved more than inhibitors (log OR 0.31), so the net effect is ~0.4.
# Knockout mice (IMPC) always say "activate" and their signal is exactly this
# bias (log OR 0.26): they are not used. CTD comes out inverted, PubTator and expression, null.
# See the README.
WEIGHT_DIRECTION = {"dir_genetic": 0.4}
# Penalty (in logit) for the fraction of PubTator articles that say the drug CAUSES
# the disease (relative to those that say it treats it), if there are enough (≥ MIN_ARTICLES_CAUSE).
# Of the tested pairs with > 50 articles, if most say "causes" 4.5 % are approved;
# if none says so, 43 %. See the README.
WEIGHT_CAUSE = 2.0
MIN_ARTICLES_CAUSE = 6


_ORPHA = {}


def orphanet_genes(n):
    """gene × disease matrix of Orphanet's disease genes (data/ext/orphadata_genes.parquet), built once."""
    if n not in _ORPHA:
        o = pd.read_parquet(PREP.parent / "ext" / "orphadata_genes.parquet")
        o = o[(o.disease < n) & (o.gene < n)]
        _ORPHA[n] = sp.csc_matrix((np.ones(len(o), np.float32), (o.gene, o.disease)), shape=(n, n))
    return _ORPHA[n]


_QTL = []


def qtl_direction():
    """The genetic direction of evaluation/m2_direction_feature.py, built once (data/ext/coloc_qtl.parquet)."""
    if not _QTL:
        from evaluation.m2_direction_feature import GeneticDirection
        _QTL.append(GeneticDirection())
    return _QTL[0]


class Graph:
    def __init__(self, prep=PREP):
        if not (prep / "nodes.pkl").exists():
            raise SystemExit("The graph is missing. Run:  py -3 download.py  and  py -3 prepare.py")
        self.prep = prep
        self.nodes = pd.read_pickle(prep / "nodes.pkl")
        meta = json.loads((prep / "relations.json").read_text("utf-8"))
        self.relations = meta["relations"]
        self.labels = {int(k): v for k, v in meta["labels"].items()}
        self.N = len(self.nodes)
        self.types = self.nodes.node_type.to_numpy()
        self.names = self.nodes.node_name.astype(str).to_numpy()
        self.origin = self.nodes.origin.to_numpy() if "origin" in self.nodes else np.full(self.N, "PrimeKG")
        self._names_lower = pd.Series(self.names).str.lower()
        self.drugs = np.flatnonzero(self.types == "drug")
        self._m, self._w, self._total, self._attr = {}, {}, None, None
        f = prep / "drugs.pkl"
        self.descriptions = pd.read_pickle(f).set_index("node_index") if f.exists() else None

    # --- matrices --------------------------------------------------------

    def M(self, rel):
        if rel not in self._m:
            if "|" in rel:
                # Pieces that do not exist (a source not downloaded) are skipped.
                pieces = [r for r in rel.split("|") if r in self.relations or "|" in r]
                if not pieces:
                    raise SystemExit(f"None of the relations in {rel!r} exists")
                m = sum(self.M(r) for r in pieces).tocsr()
                m.data[:] = 1
            elif rel in self.relations:
                m = sp.load_npz(self.prep / self.relations[rel]["file"]).tocsr()
            else:
                raise SystemExit(f"Unknown relation: {rel!r}. Available: {', '.join(self.relations)}")
            self._m[rel] = m
        return self._m[rel]

    def register(self, name, m):
        """Adds a computed relation (e.g. filtered by year)."""
        self._m[name] = m.tocsr()
        self._w.pop(name, None)

    def W(self, rel, damping=0.5):
        """Degree-normalised matrix: D^-damping · M · D^-damping."""
        if rel not in self._w:
            m = self.M(rel).astype(np.float64)
            g = np.asarray(m.sum(axis=1)).ravel()
            inv = np.zeros_like(g)
            inv[g > 0] = g[g > 0] ** -damping
            d = sp.diags(inv)
            self._w[rel] = (d @ m @ d).tocsr()
        return self._w[rel]

    def total(self):
        if self._total is None:
            self._total = sp.load_npz(self.prep / "all.npz").tocsr()
        return self._total

    def neighbors(self, i, rel):
        m = self.M(rel)
        return m.indices[m.indptr[i]:m.indptr[i + 1]]

    def neighbors_all(self, i):
        t = self.total()
        return t.indices[t.indptr[i]:t.indptr[i + 1]]

    def attributes(self):
        """Attributes of the external edges (phase, year, action, score…)."""
        if self._attr is None:
            f = self.prep / "attributes.parquet"
            self._attr = pd.read_parquet(f) if f.exists() else pd.DataFrame(
                columns=["relation", "display_relation", "x_index", "y_index"])
        return self._attr

    def matrix_until(self, relation, subtypes, cutoff_year):
        """Matrix of an external relation with only the edges from before the cutoff year."""
        a = self.attributes()
        a = a[(a.relation == relation) & a.display_relation.isin(subtypes) & (a["year"] < cutoff_year)]
        x, y = a.x_index.to_numpy(), a.y_index.to_numpy()
        m = sp.coo_matrix((np.ones(2 * len(a)), (np.r_[x, y], np.r_[y, x])), shape=(self.N, self.N)).tocsr()
        m.data[:] = 1
        return m

    # --- names -----------------------------------------------------------

    def resolve(self, text, node_type=None):
        """Indices of the matching nodes. "#123" is the direct index."""
        if text.startswith("#"):
            return [int(text[1:])]
        t = text.strip().lower()
        mask = np.ones(self.N, bool) if node_type is None else self.types == node_type
        exact = np.flatnonzero(mask & (self._names_lower.to_numpy() == t))
        if len(exact):
            return list(exact)
        return list(np.flatnonzero(mask & self._names_lower.str.contains(t, regex=False).to_numpy()))

    def name(self, i):
        return self.names[i]

    # --- scoring ---------------------------------------------------------

    def _nearby(self, genes_rel, diseases):
        """Diseases too similar to each disease (the same disease under another name):
        two steps or fewer in the ontology, or with a gene Jaccard index ≥ 0.5."""
        k = len(diseases)
        O = self.M(ONTOLOGY)
        near = O[:, diseases] + O @ O[:, diseases]
        G = self.M(genes_rel).astype(np.float64)
        shared = (G @ G[:, diseases]).tocoo()
        degree = np.asarray(G.sum(axis=1)).ravel()
        jac = shared.data / (degree[shared.row] + degree[diseases][shared.col] - shared.data)
        sel = jac >= 0.5
        near = near + sp.csr_matrix((np.ones(sel.sum()), (shared.row[sel], shared.col[sel])), shape=(self.N, k))
        near.data[:] = 1
        return near

    def scores(self, diseases, paths=None, strict=False):
        """For each metapath, a (drugs x diseases) matrix of DWPC scores.

        strict: the "similar" metapath cannot use near-identical
        diseases (see _nearby) and "ontology" is not computed.
        """
        diseases = np.asarray(diseases)
        k = len(diseases)
        paths = paths or list(METAPATHS)
        f = self.drugs
        res = {}
        for name in paths:
            kind, rel, genes_rel, _, *extra = METAPATHS[name]
            if kind == "direct":
                m = self.W(rel) @ self.W(genes_rel)[:, diseases]
            elif kind == "ppi":
                ppi = extra[0] if extra else PPI
                m = self.W(rel) @ (self.W(ppi) @ self.W(genes_rel)[:, diseases])
            elif kind == "similar":
                # Diseases that share genes, removing the disease itself
                # (otherwise the already indicated drugs would sneak in through the back door).
                U = (self.W(genes_rel) @ self.W(genes_rel)[:, diseases]).tocsc()
                remove = sp.csc_matrix((np.ones(k), (diseases, np.arange(k))), shape=U.shape)
                if strict:
                    remove = remove + self._nearby(genes_rel, diseases)
                    remove.data[:] = 1
                U = U - U.multiply(remove)
                m = self.W(rel) @ U
            elif kind in ("agrees", "contradicts"):
                act, inh = EFFECTS[rel]
                need_act = self.W("|".join(f"{b}:activate" for b in genes_rel.split("|")))[:, diseases]
                need_inh = self.W("|".join(f"{b}:inhibit" for b in genes_rel.split("|")))[:, diseases]
                if kind == "agrees":
                    m = self.W(act) @ need_act + self.W(inh) @ need_inh
                else:
                    m = self.W(act) @ need_inh + self.W(inh) @ need_act
            elif kind == "ontology":
                if strict:
                    continue
                m = self.W(rel) @ self.W(ONTOLOGY)[:, diseases]
            elif kind == "orphanet_genes":
                res[name] = np.log1p((self.M(rel)[f] > 0).astype(np.float32) @ orphanet_genes(self.N)[:, diseases]).toarray()
                continue
            elif kind == "qtl_direction":
                # no drug–disease edge: the same in every graph; 0 when no directed target (the logistic model needs
                # finite values; M2 tested NaN with LightGBM only)
                res[name] = np.nan_to_num(qtl_direction().column(f, diseases), nan=0.0)
                continue
            res[name] = np.asarray(m[f].todense())
        return res

    @staticmethod
    def model():
        return json.loads(MODEL.read_text("utf-8")) if MODEL.exists() else None

    @staticmethod
    def apply_model(scores, model):
        """Score of the learned model (logit) for each drug and disease."""
        z = model["intercept"]
        for name, coef, scale in zip(model["paths"], model["coef"], model["scale"]):
            z = z + coef * np.arcsinh(scores[name] / scale)
        return z

    @staticmethod
    def combine(scores, method="rank"):
        """Combines the metapaths into a single score per drug.

        max:  mean of each metapath scaled to its maximum (per disease).
        rank: mean of the percentiles within each metapath (zeros count as 0).
        """
        parts = []
        for p in scores.values():
            if method == "max":
                mx = p.max(axis=0, keepdims=True)
                parts.append(np.divide(p, mx, out=np.zeros_like(p), where=mx > 0))
            else:
                order = p.argsort(axis=0).argsort(axis=0).astype(np.float64)
                nz = (p > 0).sum(axis=0, keepdims=True)
                pos = order - (len(p) - nz)  # 0 for the worst non-zero
                parts.append(np.where(p > 0, (pos + 1) / np.maximum(nz, 1), 0.0))
        return np.mean(parts, axis=0)

    def direction_adjustment(self, disease):
        """Logit adjustment per drug depending on whether its mechanism (ChEMBL) acts on disease
        genes in the direction that human genetics or knockout mice ask for:
        +weight/2 if most targets go in favour, -weight/2 if they go against (per source).
        Also returns the detail per drug: [(source, gene, action, what is needed), ...]."""
        a = self.attributes()
        res = np.zeros(len(self.drugs))
        detail = {}
        for source, weight in WEIGHT_DIRECTION.items():
            need = a[(a.relation == source) & (a.x_index == disease)]
            if not len(need):
                continue
            need = dict(zip(need.y_index, need.display_relation))
            moa = a[(a.relation == "moa") & a.y_index.isin(list(need)) & a.display_relation.isin(["activates", "inhibits"])]
            for f, s in moa.groupby("x_index"):
                votes = [1 if (act == "activates") == (need[gn] == "activate") else -1
                         for gn, act in zip(s.y_index, s.display_relation)]
                res[np.searchsorted(self.drugs, f)] += weight / 2 * np.sign(sum(votes))
                detail.setdefault(f, []).extend((source, gn, act, need[gn]) for gn, act in zip(s.y_index, s.display_relation))
        return res, detail

    def literature(self, disease):
        """PubTator articles per drug: (treats, causes, pmid_min of treats) in this disease."""
        a = self.attributes()
        s = a[(a.relation == "pubtator_drug_disease") & (a.y_index == disease)]
        res = {}
        for f, t in s.groupby("x_index"):
            n = dict(zip(t.display_relation, t.n_articles))
            pm = t[t.display_relation == "treats"].pmid_min
            res[f] = (int(n.get("treats", 0)), int(n.get("causes", 0)), int(pm.iloc[0]) if len(pm) else None)
        return res

    def cause_fraction(self, disease):
        """Per drug, fraction of articles that say it causes the disease (0 if there are few)."""
        fr = np.zeros(len(self.drugs))
        for f, (t, c, _) in self.literature(disease).items():
            if t + c >= MIN_ARTICLES_CAUSE:
                fr[np.searchsorted(self.drugs, f)] = c / (t + c)
        return fr

    def failed_fraction(self, disease):
        """For each drug, fraction of its trials in this disease that were
        stopped for lack of efficacy or for safety (0 if there are none)."""
        a = self.attributes()
        c = a[(a.relation == CLINICAL) & (a.y_index == disease) & (a.n_trials > 0)]
        fr = np.zeros(len(self.drugs))
        if len(c):
            pos = np.searchsorted(self.drugs, c.x_index.to_numpy())
            fr[pos] = ((c.n_negative.fillna(0) + c.n_safety.fillna(0)) / c.n_trials).clip(0, 1).to_numpy()
        return fr

    # --- explanations ----------------------------------------------------

    def clinical_history(self, drug, disease):
        a = self.attributes()
        s = a[(a.x_index == drug) & (a.y_index == disease) & a.relation.isin([CLINICAL, NEGATIVE])]
        return s

    def action(self, drug, gene):
        a = self.attributes()
        s = a[(a.x_index == drug) & (a.y_index == gene) & (a.relation == "moa")]
        return "/".join(sorted(set(s.display_relation))) if len(s) else None

    def explain(self, drug, disease, genes_rel=GENES_ALL, limit=3):
        """Concrete paths that justify a drug's score."""
        genes = set(self.neighbors(disease, genes_rel))
        targets = self.neighbors(drug, TARGET)
        lines = []

        def target(g):
            a = self.action(drug, g)
            return f"{self.name(g)} ({a})" if a else self.name(g)

        _, detail = self.direction_adjustment(disease)
        for source, source_name in (("dir_genetic", "human genetics"),):
            d = [(gn, act, need) for f_, gn, act, need in detail.get(drug, []) if f_ == source]
            if d:
                lines.append(f"direction according to {source_name}: " + ", ".join(
                    f"{'✓' if (act == 'activates') == (need == 'activate') else '✗'} {act} {self.name(gn)} (needs {need})"
                    for gn, act, need in d[:5]))

        direct = [g for g in targets if g in genes]
        if direct:
            lines.append("target associated with the disease: " + ", ".join(target(g) for g in direct[:6]))

        a = self.attributes()
        ctd = a[(a.relation == "ctd_gene") & (a.x_index == drug) & a.y_index.isin(genes)
                & a.display_relation.isin(["increases", "decreases"])].sort_values("n_articles", ascending=False)
        if len(ctd):
            lines.append(f"changes disease genes (CTD, {len(ctd)}; does not indicate the therapeutic direction): " + ", ".join(
                f"{'↑' if r.display_relation == 'increases' else '↓'}{self.name(r.y_index)}"
                for r in ctd.head(6).itertuples()))

        indirect = []
        for t in targets:
            if t in genes:
                continue
            for g in self.neighbors(t, PPI):
                if g in genes:
                    indirect.append(f"{target(t)}↔{self.name(g)}")
        if indirect:
            lines.append(f"target that interacts with disease genes ({len(indirect)}): "
                         + ", ".join(indirect[:limit]))

        similar = []
        for d in self.neighbors(drug, f"{INDICATION}|{CLINICAL_ADVANCED}"):
            if d == disease:
                continue
            shared = len(genes.intersection(self.neighbors(d, genes_rel)))
            if shared:
                similar.append((shared, d))
        if similar:
            similar.sort(reverse=True)
            lines.append("indicated or in phase 3 for: " + ", ".join(
                f"{self.name(d)} ({c} genes in common)" for c, d in similar[:limit]))

        neighbours = set(self.neighbors(disease, ONTOLOGY)).intersection(self.neighbors(drug, INDICATION))
        if neighbours:
            lines.append("indicated for the neighbouring disease: " + ", ".join(self.name(d) for d in list(neighbours)[:limit]))

        for _, r in self.clinical_history(drug, disease).iterrows():
            if r.relation == CLINICAL:
                n = f", {int(r.n_trials)} trials" if r.n_trials == r.n_trials and r.n_trials else ""
                a = f", since {int(r['year'])}" if r["year"] == r["year"] else ""
                lines.append(f"already tested in this disease: {r.display_relation}{n}{a}")
                failed = np.nan_to_num(r.n_negative) + np.nan_to_num(r.n_safety)
                if failed and r.n_trials:
                    lines.append(f"⚠ {int(failed)} of {int(r.n_trials)} trials stopped for lack of efficacy or safety "
                                 f"(penalised)")
            else:
                lines.append(f"⚠ trial stopped for {r.display_relation}")

        t, c, pm = self.literature(disease).get(drug, (0, 0, None))
        if t or c:
            lines.append(f"{'⚠ ' if c >= t and c + t >= MIN_ARTICLES_CAUSE else ''}PubTator: {t} articles say "
                         f"it treats it, {c} that it causes it" + (f" (first: PMID {pm})" if pm else ""))
        attr = self.attributes()
        lit = attr[(attr.relation == "pubmed") & (attr.x_index == drug) & (attr.y_index == disease)]
        if len(lit):
            parts = []
            for rel, s in lit.groupby("display_relation"):
                cites = ", ".join(f"PMID {p} ({e})" for p, e in s[["pmid", "evidence"]].drop_duplicates().values[:4])
                parts.append(f"{rel}: {cites}")
            lines.append("literature (PubMed): " + "; ".join(parts))
        return lines
