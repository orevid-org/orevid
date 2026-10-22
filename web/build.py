"""Builds the data for the web explorer (web/public/data/): for each disease, approved drugs,
repurposing candidates (with explanation and direction), drugs with no prior use, and targets.

Explanations are structured (not sentences) so the web can show them in Catalan or in English.

    py -3 -m web.build --test             # a handful of diseases (prototype)
    py -3 -m web.build                    # every disease with genes or phenotypes (~1 h)
"""
import argparse
import json
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from graph import (CLINICAL, CONTRA, EFFECT, GENES_ALL, GENES_ANIMAL, GENES_GENETIC, INDICATION, NEGATIVE, OFF_LABEL,
                   PPI, TARGET)

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "web" / "public" / "data"
FRAGMENTS = 500  # files d/<id % 500>.json
TEST_DISEASES = ["multiple sclerosis", "Alzheimer disease", "Parkinson disease", "cystic fibrosis",
                 "type 2 diabetes mellitus", "asthma", "rheumatoid arthritis", "amyotrophic lateral sclerosis",
                 "Huntington disease", "psoriasis", "Crohn disease", "systemic lupus erythematosus",
                 "idiopathic pulmonary fibrosis", "glioblastoma", "Duchenne muscular dystrophy", "migraine disorder",
                 "major depressive disorder", "schizophrenia", "hypertension", "heart failure", "osteoarthritis",
                 "endometriosis", "fibromyalgia", "sickle cell anemia"]


# Evidence that justifies showing a candidate (abstention rule, docs/DECISIONS.md P27, which replaces
# P14/P21). Hit rate by evidence profile within the top 25 (evaluation/evidence.py, benchmark and
# confirmation on 9 areas): only «indicated for a neighbouring disease» (neighbour) and «indicated for a disease
# with ≥ 3 shared genes» (shared_genes) count: 18-33 % hits. «Acts on a disease gene» without either hits
# 2-3 %, like having no path at all; «disease 2 steps away» alone, 1-11 %. They are shown as context.
MIN_SHARED_GENES = 3
# In addition, the overlap must be significant (hypergeometric, p < 1e-5): overlaps with «hub» diseases
# (thousands of genes) hit 0-2.5 %, significant ones 23-26 % (DECISIONS.md P28). Since P40 the population N is
# the number of genes (before, every linked node, genes and diseases, with p < 1e-6; docs/HYPERGEOMETRIC_N.md).
# The genes shared are context, not evidence of a shared mechanism (P38).
MAX_P_SHARED = 1e-5
# «Hub» genes (associated with ≥ 500 diseases: IL6, TNF, HLA genes…) do not count for the shared-genes line, neither in
# the number of shared genes nor in the hypergeometric test (R1b, DECISIONS.md P48: the candidates they alone brought
# hit 1-3 %, like hidden ones).
HUB_DISEASES = 500


# «Specific» gene: associated with fewer than 100 diseases. With ≥ 3 specific shared genes, candidates
# hit 32-33 %; if the shared genes are mostly heavily studied genes, 9-11 % (P29). Shown, not filtered.
MAX_SPECIFIC_DEGREE = 100
MIN_SPECIFIC = 3


def valid_similar(m):
    """m = [disease, shared genes, names (specific first), p, specific genes]"""
    return m[1] >= MIN_SHARED_GENES and (len(m) < 4 or m[3] < MAX_P_SHARED)


def evidence_lines(expl):
    """Lines of evidence of a candidate: neighbour (neighbouring disease), shared_genes (≥ 3 shared genes),
    specific_genes (≥ 3 of them specific), target (acts on a disease gene)."""
    L = set()
    for e in expl:
        if e["t"] == "neighbour":
            L.add("neighbour")
        elif e["t"] == "similar" and any(valid_similar(m) for m in e["m"]):
            L.add("shared_genes")
            if any(valid_similar(m) and len(m) > 4 and m[4] >= MIN_SPECIFIC for m in e["m"]):
                L.add("specific_genes")
        elif e["t"] == "target":
            L.add("target")
    return L


def has_evidence(expl):
    return bool(evidence_lines(expl) & {"neighbour", "shared_genes"})


# Fc receptors and C1q: DrugBank lists them as «targets» of many antibodies (they bind through the Fc region),
# but they are not the mechanism. They are removed when the drug has other targets (expert review: sarilumab
# showed up as «FCGR1A, FCGR2A…» instead of IL6R).
FC = re.compile(r"^(FCGR\w*|FCER\w*|MS4A2|FCAR|FCAMR|FCRL\w*|FCGRT|C1Q[ABC]|C1R|C1S)$")


def mechanism_targets(g, targets):
    real = [t for t in targets if not FC.match(g.name(t))]
    return real if real else list(targets)


def explanation_text(e):
    """A structured explanation as plain text (English), as the web shows it (app.js)."""
    t = e["t"]
    if t == "target":
        return "acts on disease genes: " + ", ".join(g + (f" ({a})" if a else "") for g, a in e["g"])
    if t == "ctd":
        return (f"changes disease genes (CTD, {e['n']}; does not say whether it is therapeutic): "
                + ", ".join(x[0] for x in e["g"]))
    if t == "interaction":
        return f"its targets interact with disease genes ({e['n']}): " + ", ".join(f"{a}↔{b}" for a, b in e["g"])
    if t == "similar":
        return "indicated for diseases sharing genes: " + "; ".join(
            f"{m[0]} ({m[1]} genes, p = {m[3]:.0e}: {', '.join(m[2])})" for m in e["m"])
    if t == "neighbour":
        return "indicated for the neighbouring disease in the ontology: " + ", ".join(e["m"])
    if t == "neighbour2":
        return "indicated for a disease 2 steps away in the ontology: " + ", ".join(f"{m} (via {v})" for m, v in e["m"])
    if t == "phenotype":
        return "indicated for diseases sharing phenotypes: " + "; ".join(f"{m} ({c} phenotypes)" for m, c in e["m"])
    if t == "adverse_effect":
        return "recorded adverse effects (SIDER) that are phenotypes of the disease: " + ", ".join(e["f"])
    if t == "similar_drug":
        return "shares targets with drugs indicated for diseases sharing genes: " + "; ".join(
            f"{f2} → {m2} (via {d})" for f2, m2, d in e["g"])
    if t == "clinical":
        return (f"already tested in this disease: {e['phase']}" + (f", {e['n']} trials" if e.get("n") else "")
                + (f", since {e['year']}" if e.get("year") else "")
                + (f" — ⚠ {e['stopped']} stopped" if e.get("stopped") else ""))
    if t == "stopped":
        return f"⚠ trial stopped for {e['reason']}"
    if t == "literature":
        return f"PubTator: {e['treats']} articles say it treats it, {e['causes']} that it causes it"
    return str(e)


class Explainer:
    """Fast, structured explanations for a (drug, disease) pair."""

    def __init__(self, g):
        from evaluation.twins import twins
        self.g = g
        self.drug_twins = twins(g, "drug")  # the same drug in several nodes (P22): records may sit on either
        a = g.attributes()
        moa = a[a.relation == "moa"]
        self.action = moa.groupby(["x_index", "y_index"]).display_relation.first().to_dict()
        cl = a[a.relation.isin([CLINICAL, NEGATIVE])]
        self.clinical = defaultdict(list)
        for r in cl.itertuples():
            self.clinical[(r.x_index, r.y_index)].append(r)
        pt = a[a.relation == "pubtator_drug_disease"]
        self.lit = defaultdict(dict)
        for f, d, dr, n in zip(pt.x_index, pt.y_index, pt.display_relation, pt.n_articles.fillna(0)):
            self.lit[(f, d)][dr] = int(n)
        self.genes = g.M(GENES_ALL)
        self.targets = g.M(TARGET)
        self.targets_c = self.targets.tocsc()
        self.ppi = g.M(f"{PPI}|string")
        # PrimeKG (DrugBank) indications: the definition validated in evaluation/evidence.py. Open Targets
        # approvals do not count (partly wrong, P17; 3rd review, P31), nor does phase 3.
        self.ind = g.M(INDICATION)
        self.ont = g.M("disease_disease")
        O = self.ont.astype(np.float64)
        self.O2 = ((O + O @ O) > 0).astype(np.float64).tocsr()
        self.phen = g.M("disease_phenotype_positive")
        self.genes_c = self.genes.tocsc()
        self.n_genes = np.diff(self.genes.indptr)
        self.gene_degree = np.diff(self.genes_c.indptr)  # diseases per gene
        # genes with at least one disease (the matrix is symmetric: the diseases must not be counted; P40)
        self.N_genes = int(((self.gene_degree > 0) & (g.types == "gene/protein")).sum())
        # R1b (P48): hub genes are left out of the shared-genes line (count, population and gene-set sizes)
        self.hub = (self.gene_degree >= HUB_DISEASES) & (g.types == "gene/protein")
        self.n_genes_nh = self.n_genes - np.asarray(self.genes @ self.hub.astype(np.float64)).ravel().astype(int)
        self.N_genes_nh = int(((self.gene_degree > 0) & (g.types == "gene/protein") & ~self.hub).sum())
        # Adverse effects (SIDER) specific enough (≤ 5 % of drugs with effects): not a validated warning
        # (they do not separate contraindications from indications, DECISIONS.md P27), but an expert wants to see them.
        self.effects = g.M(EFFECT)
        has = np.diff(self.effects.indptr) > 0
        freq = np.asarray((self.effects[has] > 0).sum(axis=0)).ravel() / max(int(has.sum()), 1)
        self.specific_effect = freq <= 0.05
        self._d = None

    def _neighbors(self, m, i):
        return m.indices[m.indptr[i]:m.indptr[i + 1]]

    def _prepare(self, d):
        """Disease data shared by all the explanations."""
        if self._d == d:
            return
        self._d = d
        self.genes_d = set(self._neighbors(self.genes, d))
        self.phen_d = set(self._neighbors(self.phen, d))
        # diseases sharing genes (excluding itself), for «similar drugs»
        shared = (self.genes @ self.genes[d].T).tocoo()
        self.similar_d = {r for r, v in zip(shared.row, shared.data) if r != d and v >= 2}

    def explain(self, f, d, twins=None):
        """List of explanations; the evidence types are model paths (abstention rule)."""
        g = self.g
        self._prepare(d)
        twins = set(twins or {d})  # nodes that are the same disease (they do not count as «nearby»)
        genes = self.genes_d
        targets = mechanism_targets(g, self._neighbors(self.targets, f))
        e = []
        direct = [t for t in targets if t in genes]
        if direct:
            e.append({"t": "target", "g": [[g.name(t), self.action.get((f, t), "")] for t in direct[:6]]})
        # (no CTD explanations: the source is not used in production; DECISIONS.md P20)
        ind = []
        for t in targets:
            if t in genes:
                continue
            for x in self._neighbors(self.ppi, t):
                if x in genes:
                    ind.append([g.name(t), g.name(x)])
        if ind:
            e.append({"t": "interaction", "n": len(ind), "g": ind[:3]})
        sim = []
        genes_nh = {x for x in genes if not self.hub[x]}  # hub genes do not count (R1b, P48)
        for m in self._neighbors(self.ind, f):
            if m in twins:
                continue
            shared = genes_nh.intersection(self._neighbors(self.genes, m))
            if shared:
                # probability of sharing this many non-hub genes by chance (hypergeometric)
                p = float(hypergeom.sf(len(shared) - 1, self.N_genes_nh, self.n_genes_nh[m], len(genes_nh)))
                order = sorted(shared, key=lambda x: (self.gene_degree[x], g.name(x)))  # specific first
                spec = int(sum(self.gene_degree[x] < MAX_SPECIFIC_DEGREE for x in shared))
                sim.append((p, -len(shared), int(m), [g.name(x) for x in order], spec))
        if sim:
            sim.sort()
            e.append({"t": "similar", "m": [[g.name(m), -c, gs[:6], float(f"{p:.1e}"), spec]
                                            for p, c, m, gs, spec in sim[:3]]})
        ind_f = set(self._neighbors(self.ind, f)) - twins
        # neighbours of the disease or of any of its twins (like the «extensions» list in DiseasePages)
        neigh_d = {int(x) for b in twins for x in self._neighbors(self.ont, b)} - twins
        neighbours = sorted(neigh_d.intersection(ind_f))
        if neighbours:
            e.append({"t": "neighbour", "m": [g.name(m) for m in neighbours[:3]]})
        neigh2 = {int(x) for b in twins for x in self._neighbors(self.O2, b)}
        neighbours2 = sorted((neigh2 - set(neighbours) - twins).intersection(ind_f))
        if neighbours2:
            # with the intermediate disease, so the expert can see where the path goes through
            e.append({"t": "neighbour2",
                      "m": [[g.name(m), g.name(next(iter(sorted(neigh_d.intersection(self._neighbors(self.ont, m)))), m))]
                            for m in neighbours2[:3]]})
        if self.phen_d:
            phen = []
            for m in ind_f:
                c = len(self.phen_d.intersection(self._neighbors(self.phen, m)))
                if c:
                    phen.append((c, m))
            if phen:
                phen.sort(reverse=True)
                e.append({"t": "phenotype", "m": [[g.name(m), c] for c, m in phen[:3]]})
        ef = [x for x in self._neighbors(self.effects, f) if self.specific_effect[x] and x in self.phen_d]
        if ef:
            e.append({"t": "adverse_effect", "f": sorted(g.name(x) for x in ef)[:5]})
        if not sim and self.similar_d:
            # shares targets with a drug indicated for a disease sharing genes
            ex = []
            for t in targets[:20]:
                for f2 in self.targets_c.indices[self.targets_c.indptr[t]:self.targets_c.indptr[t + 1]]:
                    if f2 == f:
                        continue
                    m2 = self.similar_d.intersection(self._neighbors(self.ind, f2))
                    if m2:
                        ex.append([g.name(f2), g.name(next(iter(m2))), g.name(t)])
                        break
                if len(ex) >= 2:
                    break
            if ex:
                e.append({"t": "similar_drug", "g": ex})
        ft = self.drug_twins.get(f, {f})
        for r in [r for x in ft for m in twins for r in self.clinical.get((x, m), [])]:
            if r.relation == CLINICAL:
                stopped = int(np.nan_to_num(r.n_negative) + np.nan_to_num(r.n_safety))
                year = getattr(r, "year", None)
                e.append({"t": "clinical", "phase": r.display_relation,
                          "n": int(r.n_trials) if r.n_trials == r.n_trials else None,
                          "year": int(year) if year is not None and year == year else None,
                          "stopped": stopped})
            else:
                e.append({"t": "stopped", "reason": r.display_relation})
        lit = next((self.lit[(x, m)] for x in ft for m in twins if (x, m) in self.lit), None)
        if lit:
            e.append({"t": "literature", "treats": lit.get("treats", 0), "causes": lit.get("causes", 0)})
        return e


def external_ids(g):
    n = g.nodes
    return dict(zip(n.node_index, n.node_id.astype(str)))


class DiseasePages:
    """Builds the page of a disease (candidates, no prior use, targets, what is already known) with
    all the production rules. Used by the web (build) and by `explore.py disease`, so that they cannot
    drift apart (second independent review)."""

    def __init__(self, top=25, top_no_use=8, top_targets=12):
        from predictor import KNOWN_USE, Predictor, developed_drugs, diagnostic_drugs, drug_status
        from evaluation.direction import Direction, checkpoint_rule
        from evaluation.targets import TARGETS_MODEL, Targets
        from evaluation.twins import key, names_with_use, twins
        from graph import Graph
        t0 = time.time()
        self.top, self.top_no_use, self.top_targets = top, top_no_use, top_targets
        self.P = Predictor()
        self.g = self.P.g
        self.gc = Graph()  # full graph (every drug and disease)
        self.D = Direction()
        self.DI = Targets() if TARGETS_MODEL.exists() else None
        self.X = Explainer(self.g)
        self.ext = external_ids(self.g)
        gc = self.gc
        self.degree = np.diff(gc.M(KNOWN_USE).indptr)  # known use in any source (node index)
        self.blocked, self.autoimmune = checkpoint_rule(gc)
        self.diagnostic = diagnostic_drugs(gc)
        self.developed = list(developed_drugs(gc) - self.diagnostic)
        self.status = drug_status(gc)
        # Twins (same drug or disease in more than one node; DECISIONS.md P22): compared by name.
        self.twins = twins(gc, "disease")
        self.drug_twins = twins(gc, "drug")  # e.g. belimumab: DrugBank node + ChEMBL node holding its trials
        used_names = names_with_use(gc, KNOWN_USE)
        self.no_use_name = np.array([key(gc.names[f]) not in used_names for f in self.g.drugs])
        self.atc = atc_classes(gc)
        # Open Targets approvals confirmed by a single-drug label or agency report (sources/ot_approvals.py)
        f_ot = ROOT / "data" / "ext" / "ot_approvals.parquet"
        ot = pd.read_parquet(f_ot) if f_ot.exists() else pd.DataFrame(columns=["x_index", "y_index", "status"])
        ot = ot[ot.status == "confirmed"]
        self.ot_confirmed = set(zip(ot.x_index.astype(int), ot.y_index.astype(int)))
        self.targets_f = gc.M(TARGET)
        self.ind_f = gc.M(INDICATION)
        # Ontology with direction, for navigation only (sources/hierarchy.py; REVIEW_RESPONSE.md T6)
        from sources.hierarchy import load as load_hierarchy
        h = load_hierarchy()
        self.parents, self.children = defaultdict(list), defaultdict(list)
        for a_, b_ in zip(h.parent.astype(int), h.child.astype(int)):
            self.parents[b_].append(a_)
            self.children[a_].append(b_)
        self.n_genes = np.diff(self.g.M(GENES_ALL).indptr)
        self.markers = marker_rules(gc, self.children, self.targets_f)
        self.practical, self.trials_pair, self.papers_pair = practical_tables()
        pkg = (self.g.types == "disease") & (self.g.origin == "PrimeKG")
        self.has_page = pkg & ((self.n_genes > 0) | (np.diff(self.g.M("disease_phenotype_positive").indptr) > 0))
        print(f"models loaded ({time.time() - t0:.0f} s)", flush=True)

    def disease(self, d):
        from evaluation.twins import drug_family, twin_neighbors  # salts, esters, prodrugs, enantiomers
        from evaluation.direction import risk_flags
        from evaluation.targets import prioritise
        g, gc, X, status = self.g, self.gc, self.X, self.status
        d = int(d)
        dt = self.twins.get(d, {d})
        ind_pkg = twin_neighbors(gc, dt, INDICATION)
        offlabel, contra = twin_neighbors(gc, dt, OFF_LABEL), twin_neighbors(gc, dt, CONTRA)
        # Open Targets approvals that are not in DrugBank: shown separately (partly wrong, P17), without
        # those the graph itself says are contraindicated (contradictory data).
        names_contra = {x for i in contra for x in drug_family(g.name(i))}
        ind_ot = {i for i in twin_neighbors(gc, dt, "clinical:approved") - ind_pkg
                  if not (drug_family(g.name(i)) & names_contra)}
        indicated = ind_pkg | ind_ot  # to exclude candidates, all of them (cautious)
        ot_conf = {i for i in ind_ot if any((i, m) in self.ot_confirmed for m in dt)}
        names_ind = {x for i in indicated | offlabel for x in drug_family(g.name(i))}
        drugs, total, _ = self.P.predict(d)
        # The order is only that of the evaluated model (with specificity); failed trials and literature are
        # warnings inside the explanations (DECISIONS.md, P3).
        # Risk warning: a transparent rule (docs/DIRECTION.md), not a model; no probability.
        risk_f, n_contra, n_ind, checkpoint = risk_flags(self.D, d, drugs, self.blocked, self.autoimmune)
        extensions = twin_neighbors(gc, dt, "disease_disease") - dt
        broader, narrower = self.parents.get(d, []), self.children.get(d, [])
        equivalent = {m for m in set(broader) | set(narrower) | extensions if self.equivalent(d, m)}

        def candidate(j, score):
            f = int(drugs[j])
            reasons = {}
            if n_contra[j] > n_ind[j]:
                reasons = {"contraindicated": int(n_contra[j]), "indicated": int(n_ind[j])}
            if checkpoint[j]:
                reasons["checkpoint"] = 1
            expl = X.explain(f, d, twins=dt)
            return {"node": f, "name": g.name(f), "id": self.ext.get(f, ""), "score": round(float(score), 3),
                    "status": status.get(f, "unknown"), "atc": self.atc.get(f, ""),
                    "mechanism": self.mechanism(f), "indications": self.indications(f),
                    "lines": sorted(evidence_lines(expl)),
                    "risk": bool(risk_f[j]), "risk_reasons": reasons, "contra": f in contra, "expl": expl,
                    **self.practical_info(f, dt)}

        # Subtypes defined by an absent biomarker (resources/biomarker_subtypes.csv, REVIEW_RESPONSE.md T7):
        # drugs acting on that marker cannot work here (tamoxifen in triple-negative breast cancer).
        marker_block = {}
        for marker, status_, drugs_ in self.markers.get(d, []):
            if status_ == "absent":
                for f in drugs_:
                    marker_block.setdefault(f, marker)
        removed_marker = {}

        def known_or_unfit(f):
            # already approved (also by name: twins), diagnostic, contraindicated here, veterinary only
            fam = drug_family(g.name(f))
            return (f in indicated or f in offlabel or bool(fam & names_ind) or f in self.diagnostic
                    or f in contra or bool(fam & names_contra) or status.get(f) == "veterinary only")

        def excluded(f):
            return known_or_unfit(f) or f in marker_block

        # Abstention rule (validated, DECISIONS.md P27-P28, P31): of the top 25 non-excluded candidates
        # (direct extensions included, as in evaluation/evidence.py), only those with evidence.
        cand, ext_cands, hidden, seen = [], [], 0, 0
        for j in np.argsort(-total):
            if seen >= self.top:
                break
            f = int(drugs[j])
            if excluded(f):
                if f in marker_block and not known_or_unfit(f):
                    removed_marker.setdefault(g.name(f), marker_block[f])
                continue
            seen += 1
            in_ext = [m for m in g.neighbors(f, INDICATION) if m in extensions]
            if in_ext:
                # direct extension: already indicated for a neighbouring disease (the best hit rate; P27)
                c = candidate(j, total[j])
                c["neighbour_diseases"] = [g.name(m) for m in in_ext[:3]]
                same = [m for m in in_ext if m in equivalent]
                if same:  # indicated for a near-synonym (breast cancer / carcinoma): display only (T6)
                    c["equivalent"] = g.name(same[0])
                ext_cands.append(c)
                continue
            c = candidate(j, total[j])
            if has_evidence(c["expl"]):
                cand.append(c)
            else:
                hidden += 1
        new = np.flatnonzero((self.degree[drugs] == 0) & self.no_use_name & np.isin(drugs, self.developed)
                             & ~np.array([excluded(int(f)) for f in drugs]))
        new = new[np.argsort(-self.P.last_no_popularity[new])]
        risk = risk_f[new]  # those with a risk warning are discarded
        # With no indications they cannot have neighbour or shared_genes: an exploratory list, not validated
        # (Precision@K 0-2 %, docs/NOVELTY.md). At least they must act on a disease gene.
        no_use = [c for c in (candidate(j, self.P.last_no_popularity[j]) for j in new[~risk][:self.top_no_use])
                  if "target" in c["lines"]]
        excluded_names = [g.name(int(drugs[j])) for j in new[risk][:8]]
        targets = None
        if self.DI is not None and self.DI.E[d].nnz:
            t = prioritise(self.DI, d, self.top_targets, with_evidence=False)
            targets = {"exploited": t[t["already_target_here"]].gene.tolist(),
                       "new": [{"gene": r["gene"], "uniprot": r["uniprot"], "ev": r["evidence"], "dir": r["direction"],
                                "mod": r["modalities"], "safety": clean_safety(r["safety"]), "loeuf": r["LOEUF"],
                                "ess": r["essential"], "dep": r["dependent_fraction"], "tissue": r["tissue_specificity"],
                                "drugs": r["drugs acting on it"], "direct": bool(r["has_evidence"]),
                                "s": r["score"]}
                               for _, r in t[~t["already_target_here"]].iterrows()]}
        trials = {}
        for phase in ["phase 3", "phase 2", "phase 1"]:
            s_ = twin_neighbors(gc, dt, f"{CLINICAL}:{phase}")
            if len(s_):
                names = sorted({g.name(i) for i in s_})
                trials[phase] = {"n": len(names), "drugs": names[:30]}
        return {
            "id": d, "name": g.name(d), "mondo": self.ext.get(d, ""), "n_genes": int(len(g.neighbors(d, GENES_ALL))),
            "gene_sources": {k: int(len(g.neighbors(d, rel))) for k, rel in GENE_SOURCES.items()},
            "approved": sorted({g.name(i) for i in ind_pkg}), "approved_ot": sorted({g.name(i) for i in ot_conf})[:40],
            "approved_unconfirmed": sorted({g.name(i) for i in ind_ot - ot_conf})[:30],
            "offlabel": sorted({g.name(i) for i in offlabel})[:30],
            "contra": sorted({g.name(i) for i in contra})[:30], "trials": trials,
            "stopped": sorted({g.name(i) for i in twin_neighbors(gc, dt, NEGATIVE)})[:20],
            "genes": sorted(g.name(x) for x in g.neighbors(d, GENES_ALL))[:60],
            "candidates": group_salts(cand), "hidden": hidden, "little_info": len(cand) < 5,
            "no_use": no_use, "excluded": excluded_names, "extensions": group_by_class(group_salts(ext_cands)),
            "targets": targets,
            "broader": [self.link(m) for m in broader],
            "narrower": [self.link(m) for m in sorted(narrower, key=lambda m: -self.n_genes[m])],
            "n_descendants": self.descendants(d),
            "markers": sorted({(m, st) for m, st, _ in self.markers.get(d, [])}),
            "excluded_marker": [[n, m] for n, m in removed_marker.items()]}

    def practical_info(self, f, dt):
        """Generic / biosimilar status, routes, and how much the pair has been studied (I3; display only).
        «Underexplored»: no trial for the pair and at most 2 papers co-mentioning drug and disease (thresholds fixed in
        docs/IMPROVEMENT_PLAN.md)."""
        p = self.practical.get(f, {})
        ft = self.drug_twins.get(f, {f})  # the same drug in several nodes (P22): trials may sit on the ChEMBL one
        n_trials = max((self.trials_pair.get((x, m), 0) for x in ft for m in dt), default=0)
        n_papers = max((self.papers_pair.get((x, m), 0) for x in ft for m in dt), default=0)
        return {"generic": {k: p[k] for k in ("generic_us", "biosimilar_us", "generic_eu", "biosimilar_eu",
                                              "exclusivity_until", "routes") if p.get(k)},
                "n_trials": int(n_trials), "n_papers": int(n_papers),
                "underexplored": bool(n_trials == 0 and n_papers <= 2)}

    def link(self, m):
        return [int(m), self.g.name(m), bool(self.has_page[m])]

    def descendants(self, d):
        seen, stack = set(), list(self.children.get(d, []))
        while stack:
            m = stack.pop()
            if m not in seen:
                seen.add(m)
                stack.extend(self.children.get(m, []))
        return len(seen)

    _GENERIC = {"cancer", "carcinoma", "neoplasm", "neoplasms", "tumor", "tumour", "malignant", "disease",
                "disorder", "syndrome", "(disease)"}

    def equivalent(self, d, m):
        """Near-synonym nodes (breast cancer / breast carcinoma / breast neoplasm): same name once the
        generic words are removed, or associated genes almost identical (Jaccard ≥ 0.8). Display only."""
        a = set(self.g.name(d).lower().split()) - self._GENERIC
        b = set(self.g.name(m).lower().split()) - self._GENERIC
        if a and a == b:
            return True
        ga, gb = set(self.g.neighbors(d, GENES_ALL)), set(self.g.neighbors(m, GENES_ALL))
        return len(ga) >= 10 and len(gb) >= 10 and len(ga & gb) / len(ga | gb) >= 0.8

    def mechanism(self, f):
        """Drug targets with the action (ChEMBL) first; if there are none, those from DrugBank."""
        X, g = self.X, self.g
        targets = mechanism_targets(g, self.targets_f.indices[self.targets_f.indptr[f]:self.targets_f.indptr[f + 1]])
        with_action = [(g.name(t), X.action[(f, t)]) for t in targets if (f, t) in X.action]
        without = [(g.name(t), "") for t in targets if (f, t) not in X.action]
        return ([list(x) for x in (sorted(with_action) + sorted(without))[:4]]
                + ([["+", len(targets) - 4]] if len(targets) > 4 else []))

    def indications(self, f):
        """Current approved indications of the drug (to place it: «what it is used for today»)."""
        ms = sorted({self.g.name(m) for m in self.ind_f.indices[self.ind_f.indptr[f]:self.ind_f.indptr[f + 1]]})
        return ms[:4] + ([f"+{len(ms) - 4}"] if len(ms) > 4 else [])


def practical_tables():
    """Per drug: generic / biosimilar / exclusivity / routes (sources/practical.py). Per (drug, disease): number of
    registered trials (Open Targets clinical records) and of papers co-mentioning both (PubTator3)."""
    f = ROOT / "data" / "ext" / "practical.parquet"
    practical = {}
    if f.exists():
        d = pd.read_parquet(f)
        practical = {int(r.node): {k: getattr(r, k) for k in ("generic_us", "biosimilar_us", "generic_eu",
                                                            "biosimilar_eu", "exclusivity_until", "routes")}
                     for r in d.itertuples()}
    a = pd.read_parquet(ROOT / "data" / "prep" / "attributes.parquet",
                        columns=["relation", "x_index", "y_index", "n_trials", "n_articles"])
    c = a[a.relation == "clinical"]
    trials = {(int(x), int(y)): int(n) for x, y, n in zip(c.x_index, c.y_index, c.n_trials.fillna(0))}
    p = a[a.relation == "pubtator_drug_disease"].groupby(["x_index", "y_index"]).n_articles.sum()
    papers = {(int(x), int(y)): int(n) for (x, y), n in p.items()}
    return practical, trials, papers


def marker_rules(g, children, targets):
    """{disease node: [(marker, «absent»|«present», drug nodes acting on the marker)]} from
    resources/biomarker_subtypes.csv. An absent marker is also absent in the more specific forms
    (descendants). Drug rules: «endocrine» = ATC therapeutic group «endocrine therapy» (L02) except
    anti-androgens, or targets ESR1, ESR2 or CYP19A1 unless it is a cytotoxic antineoplastic (L01);
    «targets:G1,G2» = acts on one of those genes (DrugBank or ChEMBL); «drugs:A;B» = a closed list of drugs
    (e.g. antibacterials active only against gram-positive bacteria, on gram-negative infections)."""
    f = ROOT / "resources" / "biomarker_subtypes.csv"
    if not f.exists():
        return {}
    t = pd.read_csv(f, dtype={"drug_rule": str}).fillna("")
    drugs = np.flatnonzero(g.types == "drug")
    gene_of = {}
    for x in drugs:
        gene_of[int(x)] = {g.name(y) for y in targets.indices[targets.indptr[x]:targets.indptr[x + 1]]}
    desc = g.descriptions if g.descriptions is not None else pd.DataFrame()
    atc = lambda col: desc[col].fillna("").astype(str).str.lower() if col in desc else pd.Series(dtype=str)
    atc2, atc4 = atc("atc_2"), atc("atc_4")
    # L02 except anti-androgens (enzalutamide acts on AR, studied in AR-positive triple-negative breast
    # cancer); found by the manual check of T7.
    endocrine_atc = {int(i) for i, v in atc2.items() if "endocrine therapy" in v and "anti-androgens" not in atc4.get(i, "")}
    # Cytotoxic antineoplastic agents (L01) that also list ESR1 act through another mechanism
    # (estramustine: microtubules); found by the same check.
    cytotoxic = {int(i) for i, v in atc2.items() if "antineoplastic agents" in v}

    def acting(rule):
        if rule == "endocrine":
            genes = {"ESR1", "ESR2", "CYP19A1"}
            return {x for x, gs in gene_of.items() if gs & genes and x not in cytotoxic} | endocrine_atc
        if rule.startswith("targets:"):
            genes = set(rule.split(":", 1)[1].split(","))
            return {x for x, gs in gene_of.items() if gs & genes}
        if rule.startswith("drugs:"):  # a closed list of drug names (salts and esters included through drug_family)
            from evaluation.twins import drug_family
            listed = {k for n in rule.split(":", 1)[1].split(";") if n.strip() for k in drug_family(n.strip())}
            return {int(x) for x in drugs if drug_family(g.name(int(x))) & listed}
        return set()

    res = defaultdict(list)
    for node, marker, status, rule in zip(t.node, t.marker, t.status, t.drug_rule):
        node, ds = int(node), acting(rule) if status == "absent" else set()
        below, stack = set(), [node]
        while stack and status == "absent":  # an absent marker stays absent in the subtypes
            for c in children.get(stack.pop(), []):
                if c not in below:
                    below.add(c)
                    stack.append(c)
        for m in {node} | below:
            res[m].append((marker, status, ds))
    return dict(res)


# Where the associated genes come from (for the header; an expert wants to know whether it is genetics or literature)
GENE_SOURCES = {"genetic": GENES_GENETIC, "primekg": "disease_protein", "animal": GENES_ANIMAL,
                "expression": "ot_assoc:rna_expression"}

# Open Targets target safety contains ecotoxicology results («Decrease, Population growth rate»):
# they are not risks for humans (3rd review).
_ECOTOX = re.compile(r"^(decrease|decreased|increase|increased|promotion|n/a|impaired|impairment)\b.*,", re.I)


def clean_safety(text):
    if not text or text != text:
        return ""
    return "; ".join(x.strip() for x in str(text).split(";") if x.strip() and not _ECOTOX.match(x.strip()))


def group_by_class(cands):
    """Direct extensions of the same ATC class (e.g. 11 topical corticosteroids for SLE): one row with
    «also». Only for extensions, where the class is the hypothesis (3rd review)."""
    res, seen = [], {}
    for c in cands:
        # the first chunk of the class («corticosteroids, potent (group iii) and …» -> «corticosteroids»)
        k = re.split(r",| and ", c.get("atc") or "")[0].strip()
        if k and k in seen:
            seen[k].setdefault("also", []).extend([c["name"]] + c.get("also", []))
            continue
        if k:
            seen[k] = c
        res.append(c)
    return res


def atc_classes(g):
    """Level-4 ATC class (chemical and functional group) from DrugBank, as text."""
    if g.descriptions is None or "atc_4" not in g.descriptions:
        return {}
    res = {}
    for i, s_ in g.descriptions["atc_4"].dropna().astype(str).items():
        m = re.search(r" is (.+?)\.?$", s_)
        if m:
            res[int(i)] = m.group(1)
    return res


def group_salts(cands):
    """Merges salts and esters of the same active moiety (e.g. estradiol benzoate / dienanthate): the
    first one is shown with «also» (expert review, P27)."""
    from evaluation.twins import active_moiety
    res, seen = [], {}
    for c in cands:
        k = active_moiety(c["name"])
        if k in seen:
            seen[k].setdefault("also", []).append(c["name"])
            continue
        seen[k] = c
        res.append(c)
    return res


def build(diseases, **kw):
    P = DiseasePages(**kw)
    t0 = time.time()
    fragments = defaultdict(dict)
    for k, d in enumerate(diseases):
        fragments[int(d) % FRAGMENTS][str(int(d))] = P.disease(d)
        if (k + 1) % 50 == 0 or k + 1 == len(diseases):
            print(f"  {k + 1}/{len(diseases)} diseases ({time.time() - t0:.0f} s)", flush=True)
    return fragments


def clean(x):
    """NaN -> None (valid JSON), recursively; numpy types -> Python."""
    if isinstance(x, dict):
        return {k: clean(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [clean(v) for v in x]
    if isinstance(x, (float, np.floating)):
        return None if x != x else float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, np.bool_):
        return bool(x)
    return x


def patch():
    """Adds the navigation fields (REVIEW_RESPONSE.md T6) to the pages already built, and rebuilds in full the
    pages with a biomarker rule (T7), instead of rebuilding the 10,042 pages (4-5 h). It uses the same
    functions as DiseasePages.disease(), so the result is the same as a full rebuild with this version."""
    from evaluation.twins import twin_neighbors
    P = DiseasePages()
    g, gc = P.g, P.gc
    t0 = time.time()
    n = rebuilt = 0
    for f in sorted((DEST / "d").glob("*.json")):
        pages = json.loads(f.read_text("utf-8"))
        for k, page in pages.items():
            d = int(k)
            n += 1
            if d in P.markers:
                pages[k] = clean(P.disease(d))
                rebuilt += 1
                continue
            dt = P.twins.get(d, {d})
            extensions = twin_neighbors(gc, dt, "disease_disease") - dt
            broader, narrower = P.parents.get(d, []), P.children.get(d, [])
            equivalent = {m for m in set(broader) | set(narrower) | extensions if P.equivalent(d, m)}
            for c in page["extensions"]:
                c.pop("equivalent", None)
                same = [m for m in g.neighbors(c["node"], INDICATION) if m in extensions and m in equivalent]
                if same:
                    c["equivalent"] = g.name(same[0])
            page.update(clean({"broader": [P.link(m) for m in broader],
                               "narrower": [P.link(m) for m in sorted(narrower, key=lambda m: -P.n_genes[m])],
                               "n_descendants": P.descendants(d), "markers": [], "excluded_marker": []}))
        f.write_text(json.dumps(pages, ensure_ascii=False, separators=(",", ":"), allow_nan=False), "utf-8")
    print(f"patched {n} pages ({rebuilt} rebuilt in full for the biomarker rule) in {time.time() - t0:.0f} s")


def save(fragments, g, all_diseases):
    (DEST / "d").mkdir(parents=True, exist_ok=True)
    for fr, data in fragments.items():
        f = DEST / "d" / f"{fr}.json"
        old = json.loads(f.read_text("utf-8")) if f.exists() else {}
        old.update(clean(data))
        f.write_text(json.dumps(old, ensure_ascii=False, separators=(",", ":"), allow_nan=False), "utf-8")
    with_data = {int(k) for fr in fragments.values() for k in fr}
    idx = DEST / "index.json"
    old = {int(r[0]) for r in json.loads(idx.read_text("utf-8"))["diseases"] if r[2]} if idx.exists() else set()
    with_data |= old
    index = [[int(d), g.name(d), int(d) in with_data] for d in all_diseases]
    idx.write_text(json.dumps({"fragments": FRAGMENTS, "diseases": index}, ensure_ascii=False, separators=(",", ":")),
                   "utf-8")
    print(f"index: {len(index)} diseases ({len(with_data)} with data) -> {DEST}")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    p = argparse.ArgumentParser()
    p.add_argument("--test", action="store_true", help="only a handful of well-known diseases (prototype)")
    p.add_argument("--limit", type=int, help="only the first N diseases with data")
    p.add_argument("--sample", type=int, help="random sample (seed 0) of N diseases with data, for the panel")
    p.add_argument("--patch", action="store_true", help="add the navigation fields to the built pages and rebuild "
                                                         "the biomarker pages (instead of a full rebuild)")
    a = p.parse_args()
    if a.patch:
        return patch()
    from graph import Graph
    g = Graph()
    pkg = np.flatnonzero((g.types == "disease") & (g.origin == "PrimeKG"))
    with_data = pkg[(np.diff(g.M(GENES_ALL).indptr)[pkg] > 0)
                    | (np.diff(g.M("disease_phenotype_positive").indptr)[pkg] > 0)]
    if a.test:
        diseases = [g.resolve(n, "disease")[0] for n in TEST_DISEASES if g.resolve(n, "disease")]
    elif a.sample:
        diseases = np.random.default_rng(0).choice(with_data, a.sample, replace=False)
    else:
        diseases = with_data[:a.limit] if a.limit else with_data
    fr = build(diseases)
    save(fr, g, pkg)


if __name__ == "__main__":
    main()
