"""Display-only framing of a disease page (G3 display, F1; docs/IMPROVEMENT_PLAN.md, «Registered 09-11-2026 …»).
Nothing here adds, removes or reorders a candidate.

- genes_by_source: the disease's genes grouped by their strongest source (human genetics > PrimeKG curated >
  expression > mouse models only), so a reader sees which genes rest only on animal models.
- features (F1a): the diseases a proposal comes from (neighbour diseases of an extension, the diseases of a
  shared-genes line) that are a feature of the page's disease: the neighbour's MONDO term cross-references (MONDO
  `xref: HP:`) or names (name or exact synonym, case-insensitive) one of the page disease's current HPO phenotypes.
- trial_conditions (F1b): each listed trial drug with its trials (NCT id, phase) and the conditions as registered on
  ClinicalTrials.gov (API v2, `sources/ctgov.py conditions`), found through the NCT ids of the Open Targets clinical
  records the trial lists come from.
- merged (F1c): the MONDO terms a page node merges, when more than one.
- feature_pages (F1c): the page disease's HPO phenotypes that are diseases with their own page (NF1 -> MPNST); the
  relation is the HPO annotation itself.
"""
import re
from collections import defaultdict
from functools import cached_property
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SOURCE_ORDER = (("genetic", "ot_assoc:genetic_association|ot_assoc:genetic_literature|ot_assoc:somatic_mutation"),
                ("primekg", "disease_protein"), ("expression", "ot_assoc:rna_expression"), ("animal", "ot_assoc:animal_model"))
MAX_GENES, MAX_TRIALS = 60, 3
PHASE_RANK = {"APPROVAL": 7, "PHASE_4": 6, "PHASE_3": 5, "PHASE_2_3": 4, "PHASE_2": 3, "PHASE_1_2": 2, "PHASE_1": 1,
              "EARLY_PHASE_1": 0.5}


def _norm(s):
    return re.sub(r"\s+", " ", str(s).lower().replace("_", " ")).strip()


def _obo_terms(path, prefix):
    """{id: (name, [exact synonyms], [xrefs])} of an OBO file, for ids starting with `prefix`."""
    out, cur = {}, None
    for line in open(path, encoding="utf-8"):
        if line.startswith("[Term]") or line.startswith("[Typedef]"):
            cur = None
        elif line.startswith("id: "):
            i = line[4:].strip()
            cur = i if i.startswith(prefix) else None
            if cur:
                out[cur] = ["", [], []]
        elif cur and line.startswith("name: "):
            out[cur][0] = line[6:].strip()
        elif cur and line.startswith("synonym: ") and " EXACT" in line:
            out[cur][1].append(line.split('"')[1])
        elif cur and line.startswith("xref: "):
            out[cur][2].append(line[6:].split()[0])
    return out


class Framing:
    def __init__(self, g):
        self.g = g

    # --- genes by source (G3 display) ---
    @cached_property
    def _gene_mats(self):
        return [(k, self.g.M(rel)) for k, rel in SOURCE_ORDER]

    def genes_by_source(self, d):
        seen, out = set(), {}
        for k, m in self._gene_mats:
            xs = [int(x) for x in m.indices[m.indptr[d]:m.indptr[d + 1]] if int(x) not in seen]
            seen.update(xs)
            if xs:
                names = sorted(self.g.name(x) for x in xs)
                out[k] = {"n": len(names), "genes": names[:MAX_GENES]}
        return out

    # --- MONDO / HPO ---
    @cached_property
    def _mondo(self):
        return _obo_terms(ROOT / "data" / "raw" / "mondo" / "mondo.obo", "MONDO:")

    @cached_property
    def _hpo(self):
        return _obo_terms(ROOT / "data" / "raw" / "hpo" / "hp.obo", "HP:")

    def mondo_ids(self, d):
        return [f"MONDO:{int(p):07d}" for p in str(self.g.nodes.node_id[int(d)]).split("_") if p.isdigit()]

    def merged(self, d):
        ids = self.mondo_ids(d) if self.g.types[int(d)] == "disease" else []
        return [[i, self._mondo.get(i, ["", [], []])[0]] for i in ids] if len(ids) > 1 else []

    @cached_property
    def _phenotypes(self):
        """{disease node: {HP id}} from the current HPO annotations (sources/hpo_current.py)."""
        h = pd.read_parquet(ROOT / "data" / "ext" / "hpo_current.parquet", columns=["disease", "phenotype"])
        hp_of = {int(i): f"HP:{int(n):07d}" for i, n, t in zip(self.g.nodes.node_index, self.g.nodes.node_id, self.g.types)
                 if t == "effect/phenotype" and str(n).isdigit()}
        out = defaultdict(set)
        for d, p in zip(h.disease.astype(int), h.phenotype.astype(int)):
            if p in hp_of:
                out[d].add(hp_of[p])
        return out

    @cached_property
    def _hp_by_name(self):
        out = defaultdict(set)
        for i, (name, syn, _) in self._hpo.items():
            for s in [name, *syn]:
                out[_norm(s)].add(i)
        return out

    @cached_property
    def _disease_by_name(self):
        out = defaultdict(list)
        for x in np.flatnonzero(self.g.types == "disease"):
            out[self.g.name(int(x))].append(int(x))
        return out

    def hp_terms(self, n):
        """HP terms a disease node stands for: MONDO's HP cross-references, then name / exact-synonym matches."""
        xref, names = set(), {_norm(self.g.name(n))}
        for i in self.mondo_ids(n):
            name, syn, xr = self._mondo.get(i, ["", [], []])
            xref.update(x for x in xr if x.startswith("HP:"))
            names.update(_norm(s) for s in [name, *syn] if s)
        byname = set().union(*(self._hp_by_name.get(s, set()) for s in names))
        return xref, byname

    def features(self, dt, via_names):
        """{neighbour name: [HP id, HP name, route]} for the via diseases that are features of the page disease."""
        phen = set().union(*(self._phenotypes.get(int(m), set()) for m in dt))
        if not phen:
            return {}
        out = {}
        for name in sorted(set(via_names)):
            for n in self._disease_by_name.get(name, []):
                if n in dt:
                    continue
                xref, byname = self.hp_terms(n)
                hit = sorted(xref & phen)
                route = "MONDO xref"
                if not hit:
                    hit, route = sorted(byname & phen), "name"
                if hit:
                    out[name] = [hit[0], self._hpo.get(hit[0], [hit[0]])[0], route]
                    break
        return out

    @cached_property
    def _diseases_by_hp(self):
        out = defaultdict(set)
        for x in np.flatnonzero(self.g.types == "disease"):
            xref, byname = self.hp_terms(int(x))
            for h in xref | byname:
                out[h].add(int(x))
        return out

    def feature_pages(self, dt, has_page):
        """F1c: the page disease's HPO phenotypes that are diseases with their own page (MPNST for NF1), as
        [[node, name, HP id]]: the relation is stated by the HPO annotation."""
        phen = set().union(*(self._phenotypes.get(int(m), set()) for m in dt))
        out = {}
        for h in phen:
            for x in self._diseases_by_hp.get(h, ()):
                if x not in dt and has_page(x):
                    out.setdefault(self.g.name(x), [x, self.g.name(x), h])
        return [out[k] for k in sorted(out, key=str.lower)]

    # --- trials (F1b) ---
    @cached_property
    def _trials(self):
        """{(drug node, disease node): [(NCT id, phase, registered conditions)]} from Open Targets clinical reports."""
        from sources import ids
        from sources.ids import primekg_keys
        n = self.g.nodes
        index = {}
        for i, t, nid, o in zip(n.node_index, n.node_type, n.node_id, n.origin):
            for k in (primekg_keys(t, nid) if o == "PrimeKG" else [nid]):
                index.setdefault(k, int(i))
        dmap, _, _, _ = ids.diseases()
        fmap = ids.drugs()
        f = ROOT / "data" / "ext" / "ctgov_conditions.parquet"  # conditions as registered (F1b)
        registered = dict(zip(*pd.read_parquet(f).T.values)) if f.exists() else {}
        r = pd.read_parquet(ids.OT / "clinical_report", columns=["id", "type", "drugs", "diseases", "clinicalStage", "url"])
        r = r[(r.type == "CLINICAL_TRIAL") & r.diseases.notna() & r.drugs.notna()]
        out = defaultdict(list)
        for nct, drugs, dis, stage in zip(r.id, r.drugs, r.diseases, r.clinicalStage):
            conds = registered.get(str(nct).upper(), "")
            fs = {index.get(fmap[x["drugId"]][0]) for x in drugs if x.get("drugId") in fmap}
            ms = {index.get(dmap.get(x.get("diseaseId"))) for x in dis}
            for f in fs - {None}:
                for m in ms - {None}:
                    out[(f, m)].append((str(nct).upper(), str(stage or ""), conds))
        return out

    @cached_property
    def _trial_index(self):
        by_d = defaultdict(lambda: defaultdict(list))
        for (f, m), recs in self._trials.items():
            by_d[m][self.g.name(f)].extend(recs)
        return by_d

    def trial_conditions(self, dt, drug_names):
        out = {}
        want = set(drug_names)
        for m in dt:
            for name, recs in self._trial_index.get(int(m), {}).items():
                if name in want:
                    out.setdefault(name, []).extend(recs)
        res = {}
        for name, recs in out.items():
            uniq = sorted(dict.fromkeys(recs), key=lambda r: (-PHASE_RANK.get(r[1], 0), r[0]))  # later phases first
            res[name] = [list(r) for r in uniq[:MAX_TRIALS]] + ([["+", len(uniq) - MAX_TRIALS, ""]] if len(uniq) > MAX_TRIALS else [])
        return res


def via_names(page):
    """The diseases each shown proposal comes from: neighbour diseases of extensions, neighbours and shared-genes
    diseases of candidates."""
    from web.build import valid_similar
    out = set()
    for c in page.get("extensions") or []:
        if isinstance(c, dict):
            out.update(c.get("neighbour_diseases") or [])
    for c in page.get("candidates") or []:
        for e in c.get("expl") or []:
            if e.get("t") == "neighbour":
                out.update(e["m"])
            elif e.get("t") == "similar":
                out.update(m[0] for m in e["m"] if valid_similar(m))
    return out


def frame(F, d, dt, page, has_page):
    """The display fields added to a page (G3 display, F1)."""
    trial_drugs = [x for a in (page.get("trials") or {}).values() for x in (a if isinstance(a, list) else a["drugs"])]
    return {"genes_by_source": F.genes_by_source(d), "features": F.features(dt, via_names(page)),
            "trial_conditions": F.trial_conditions(dt, trial_drugs), "merged": F.merged(d),
            "feature_pages": F.feature_pages(dt, has_page)}
