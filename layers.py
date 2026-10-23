"""The two graphs: the internal general graph (every source we have, used or not) and the publishable one.

OREVID's ranking model and website read only `data/prep` (prepare.py) and a fixed list of relations
(`model_relations()`, frozen in resources/model_relations.txt and checked by tests/test_model_relations.py). The graphs
built here sit beside it and never change it:

- `data/graph/internal/`: the core graph (PrimeKG + the merged sources of data/prep, same node indices) plus every
  layer in LAYERS, including sources tested and not adopted; new node types (ATC codes, pathogen classes) are appended
  after the core nodes. Each layer is one file with its provenance in `manifest.json`.
- `data/graph/public/`: only the layers whose licence allows redistribution now, each in its own file under its own
  licence (aggregation, no merged file). Layers marked «check» or «no» need the creators' permission first.
  Nothing is uploaded anywhere by this script.

    py -3 layers.py build      # -> data/graph/internal/
    py -3 layers.py export     # -> data/graph/public/
    py -3 layers.py report     # -> docs/GRAPH_LAYERS.md
"""
import json
import re
import shutil
import sys
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
EXT = ROOT / "data" / "ext"
RAW = ROOT / "data" / "raw"
OUT = ROOT / "data" / "graph"

# --- the relations the production ranking reads ------------------------------------------------------

IMPLICIT = {"ppi": "protein_protein", "ontology": "disease_disease"}


def model_relations():
    """Relations read by the production ranking (feature set PRODUCTION: metapaths and drug degrees)."""
    import evaluation.txgnn as tx
    out = set()

    def add(rel):
        if isinstance(rel, (list, tuple)):
            for r in rel:
                add(r)
        elif isinstance(rel, str):
            if rel in tx.METAPATHS:
                resolve(rel)
            else:
                out.update(rel.split("|"))

    def resolve(name):
        kind, a, b, *rest = tx.METAPATHS[name]
        if kind == "qtl_direction":  # reads a layer of the internal graph, not data/prep (P61)
            out.add("layer:qtl_coloc")
            return
        add(a)
        add(b)
        if rest[1:]:
            add(rest[1])  # the network of a «ppi» path (STRING, phenotype_protein)
        elif kind in IMPLICIT:
            out.add(IMPLICIT[kind])

    for p in tx.FEATURE_SETS[tx.PRODUCTION]:
        resolve(p)
    out.update(tx.DEGREES)
    return out


# --- core sources (data/prep) -----------------------------------------------------------------------

CORE = {  # source: (version, licence, redistribution)
    "PrimeKG": ("Harvard Dataverse doi:10.7910/DVN/IXA7BM (TxGNN release)",
                "CC BY-NC-SA 4.0 as OREVID's data (aggregates DrugBank CC BY-NC 4.0, SIDER and DisGeNET CC BY-NC-SA 4.0)",
                "yes"),
    "Open Targets 26.06": ("26.06", "CC0 1.0", "yes"),
    "Open Targets (direction)": ("26.06, derived (sources/direction.py)", "CC0 1.0 (derived)", "yes"),
    "Open Targets (GWAS+QTL)": ("26.06, derived (sources/gwas_direction.py)", "CC0 1.0 (derived)", "yes"),
    "STRING 12.5": ("12.5", "CC BY 4.0", "yes"),
    "ChEMBL (structure)": ("Open Targets 26.06 molecules, Tanimoto ≥ 0.4", "CC0 1.0 (via Open Targets)", "yes"),
    "PubTator3": ("NCBI FTP relation2pubtator3", "Public domain (NCBI)", "yes"),
    "PubMed": ("pilot extraction (resources/pubmed_pilot_ms.jsonl)", "quoted sentences from abstracts", "check"),
    "CTD": ("CTD downloads", "CTD terms (non-commercial; redistribution needs permission)", "no"),
}


def core_layers():
    """One layer per (source, relation, subtype) of data/prep, with the attributes of the external edges."""
    e = pd.read_csv(RAW / "edges.csv", dtype={"relation": str, "display_relation": str})
    e = e.rename(columns={"x_index": "x", "y_index": "y"})
    e["source"] = "PrimeKG"
    a = pd.read_parquet(ROOT / "data" / "prep" / "attributes.parquet").rename(columns={"x_index": "x", "y_index": "y"})
    for df in (e, a):
        for (src, rel), d in df.groupby(["source", "relation"], observed=True):
            d = d.dropna(axis=1, how="all").drop(columns=["source"])
            lo, hi = np.minimum(d.x, d.y), np.maximum(d.x, d.y)
            d = d.assign(x=lo, y=hi).drop_duplicates(["x", "y", "display_relation"])
            yield src, rel, d


# --- added layers -----------------------------------------------------------------------------------

@dataclass
class Layer:
    name: str
    relation: str
    source: str
    version: str
    licence: str
    redistribution: str  # yes | check | no
    status: str
    load: Callable  # (graph, new_nodes) -> DataFrame x, y, relation, display_relation, attributes


def _graph():
    from graph import Graph
    return Graph()


def _named(df, x, y, rel, disp, **attrs):
    out = pd.DataFrame({"x": df[x].astype(int), "y": df[y].astype(int), "relation": rel,
                        "display_relation": disp if isinstance(disp, str) else df[disp].astype(str)})
    for k, v in attrs.items():
        out[k] = df[v].to_numpy()
    return out


def load_gencc(g, new):
    return _named(pd.read_parquet(EXT / "gencc.parquet"), "gene", "disease", "gene_disease_validity", "class",
                  submissions="submissions")


def load_onsides(g, new):
    return _named(pd.read_parquet(EXT / "onsides.parquet"), "drug", "disease", "label_adverse_event", "adverse event",
                  term="term", sections="sections")


def load_measured(source):
    def f(g, new):
        m = pd.read_parquet(EXT / "measured_targets.parquet")
        return _named(m[m.source == source], "drug", "gene", "measured_activity", "pX ≥ 6", px="pchembl_value")
    return f


def atc_classes():
    f = RAW / "chembl" / "atc_class.json"
    if not f.exists():
        base, url, rows = "https://www.ebi.ac.uk", "https://www.ebi.ac.uk/chembl/api/data/atc_class.json?limit=1000", []
        while url:
            with urllib.request.urlopen(url, timeout=120) as r:
                j = json.load(r)
            rows += j["atc"]
            url = base + j["page_meta"]["next"] if j["page_meta"]["next"] else None
        f.write_text(json.dumps(rows), "utf-8")
    return json.loads(f.read_text("utf-8"))


def load_atc(g, new):
    from sources.infections import atc_codes
    from sources import ids
    names = {}
    for r in atc_classes():
        for k in range(1, 5):
            names[r[f"level{k}"]] = r[f"level{k}_description"].capitalize()
        names[r["level5"]] = r["who_name"]
    n = g.nodes
    key_node = {(f"DB:{i}" if s == "DrugBank" else f"CHEMBL:{i}"): int(k)
                for k, i, t, s in zip(n.node_index, n.node_id, n.node_type, n.node_source) if t == "drug"}
    mol = ids.drugs()
    rows, parents = [], set()
    for cid, codes in atc_codes()["codes"].items():
        f = key_node.get(mol.get(cid, (None,))[0])
        if f is None:
            continue
        for code in codes:
            rows.append((f, new("ATC class", code, names.get(code, code)), "drug_atc", "ATC code"))
            for a, b in ((code[:5], code), (code[:4], code[:5]), (code[:3], code[:4]), (code[:1], code[:3])):
                parents.add((new("ATC class", a, names.get(a, a)), new("ATC class", b, names.get(b, b))))
    rows += [(a, b, "atc_parent", "parent of") for a, b in parents]
    return pd.DataFrame(rows, columns=["x", "y", "relation", "display_relation"]).drop_duplicates()


def load_infection_class(g, new):
    from sources.infections import CLASS_LABELS, disease_class_fn
    fn = disease_class_fn()
    n = g.nodes
    rows = []
    for i, t, nid, name in zip(n.node_index, n.node_type, n.node_id, n.node_name):
        if t == "disease":
            for c in fn(nid, name):
                rows.append((int(i), new("pathogen class", c, CLASS_LABELS[c]), "pathogen_class", c))
    return pd.DataFrame(rows, columns=["x", "y", "relation", "display_relation"])


def load_drugmechdb(g, new):
    from evaluation.drugmechdb import uniprot_to_node
    from sources.ids import xref_to_mondo
    paths = json.loads((RAW / "drugmechdb" / "indication_paths.json").read_text("utf-8"))
    n = g.nodes
    drug_node = {f"DB:{i}": int(k) for k, i, t in zip(n.node_index, n.node_id, n.node_type) if t == "drug"}
    mondo_node = {}
    for k, i, t in zip(n.node_index, n.node_id, n.node_type):
        if t == "disease" and str(i).replace("_", "").isdigit():
            for p in str(i).split("_"):
                mondo_node.setdefault(f"MONDO:{int(p)}", int(k))
    to_mondo, up = xref_to_mondo(), uniprot_to_node(g)
    rows = []
    for p in paths:
        f = drug_node.get(p["graph"].get("drugbank") or "")
        ds = {mondo_node[k] for k in to_mondo.get(p["graph"].get("disease_mesh") or "", ()) if k in mondo_node}
        if f is None or len(ds) != 1:
            continue
        d = ds.pop()
        pid = p["graph"].get("_id", "")
        rows.append((f, d, "curated_mechanism", "drug treats disease (DrugMechDB path)", pid))
        for node in p["nodes"]:
            if node["label"] == "Protein" and node["id"].split(":", 1)[1] in up:
                gn = up[node["id"].split(":", 1)[1]]
                rows += [(f, gn, "mechanism_path", "drug → protein in a curated path", pid),
                         (gn, d, "mechanism_path", "protein → disease in a curated path", pid)]
    return pd.DataFrame(rows, columns=["x", "y", "relation", "display_relation", "path_id"]).drop_duplicates()


def _gene_symbol_node(g):
    return {str(s).upper(): int(k) for k, s, t in zip(g.nodes.node_index, g.nodes.node_name, g.nodes.node_type)
            if t == "gene/protein"}


def load_lincs(g, new):
    from evaluation.lincs_reversal import drug_signatures
    sig, _ = drug_signatures(g)
    sym = _gene_symbol_node(g)
    rows = [(f, sym[s], "expression_signature", d) for f, (u, dn) in sig.items()
            for d, genes in (("up", u), ("down", dn)) for s in genes if s in sym]
    return pd.DataFrame(rows, columns=["x", "y", "relation", "display_relation"]).drop_duplicates()


def load_creeds(g, new):
    from evaluation.lincs_reversal import disease_signatures
    sig, _, _, _ = disease_signatures(g)
    sym = _gene_symbol_node(g)
    rows = {}
    for d, lst in sig.items():
        for u, dn in lst:
            for direction, genes in (("up", u), ("down", dn)):
                for s in genes:
                    if s in sym:
                        rows[(d, sym[s], direction)] = rows.get((d, sym[s], direction), 0) + 1
    return pd.DataFrame([(d, gn, "expression_signature", di, k) for (d, gn, di), k in rows.items()],
                        columns=["x", "y", "relation", "display_relation", "n_signatures"])


def load_coloc(g, new):
    from evaluation.m1_direction import gene_disease_signs
    out = []
    for qtl in ("pqtl", "eqtl"):
        signs, _, _, _ = gene_disease_signs(g, [qtl])
        out += [(gn, d, "qtl_colocalisation", f"{qtl}: higher level, {'higher' if s > 0 else 'lower'} risk", s)
                for (gn, d), s in signs.items()]
    return pd.DataFrame(out, columns=["x", "y", "relation", "display_relation", "sign"])


def load_labels(ema):
    def f(g, new):
        p = pd.read_parquet(EXT / "pair_label.parquet")
        is_ema = p["where"].str.contains("| EMA ", regex=False)  # where = "label title | URL (EMA …) | sentence"
        p = p[is_ema == ema]
        return _named(p, "drug", "disease", "current_label", "label", approval_source="source", where="where")
    return f


def load_verified_negations(g, new):
    v = pd.read_csv(ROOT / "resources" / "n4_verified_negations.csv")
    low = lambda t: {str(x).lower(): int(k) for k, x, ty in zip(g.nodes.node_index, g.nodes.node_name, g.nodes.node_type) if ty == t}
    dn, sn = low("drug"), low("disease")
    v["x"], v["y"] = v.drug.str.lower().map(dn), v.disease.str.lower().map(sn)
    v = v.dropna(subset=["x", "y"])
    return _named(v, "x", "y", "label_negation_verified", "not effective / not indicated (checked by hand)",
                  quote="quote", url="url", checked="checked")


def load_ot_approvals(g, new):
    return _named(pd.read_parquet(EXT / "ot_approvals.parquet"), "x_index", "y_index", "approval", "status")


def load_prism(g, new):
    p = pd.read_parquet(EXT / "prism_lineage.parquet")
    p["y"] = [new("cancer lineage", l, l.replace("_", " ")) for l in p.lineage]
    return _named(p, "drug", "y", "cell_line_sensitivity", "median AUC (PRISM)", auc="auc", cell_lines="cell_lines")


def load_civic(g, new):
    c = pd.read_parquet(EXT / "civic.parquet").dropna(subset=["drug", "disease"])
    return _named(c, "drug", "disease", "clinical_evidence", "significance", evidence_type="evidence_type",
                  direction="direction", level="level", gene="gene", molecular_profile="molecular_profile", url="url")


def load_hub(g, new):
    return _named(pd.read_parquet(EXT / "repurposing_hub.parquet"), "drug", "gene", "annotated_target", "moa", phase="phase")


def load_cpic(g, new):
    return _named(pd.read_parquet(EXT / "cpic.parquet"), "drug", "gene", "pharmacogenomic_pair", "cpic_level",
                  guideline="guideline", url="url")


def load_pharmgkb(g, new):
    return _named(pd.read_parquet(EXT / "pharmgkb.parquet"), "drug", "gene", "pharmacogenomic_annotation", "level",
                  category="category", n_annotations="n_annotations", url="url")


def load_trials(g, new):
    return _named(pd.read_parquet(EXT / "ctgov_trials.parquet"), "drug", "disease", "trial_result", "outcome",
                  nct="nct", phase="phase", status="status", min_p="min_p")


LAYERS = [
    Layer("gencc", "gene_disease_validity", "GenCC", "data/raw/gencc/VERSION", "CC0 1.0", "yes",
          "tested for the shared-genes line, not adopted (P54)", load_gencc),
    Layer("onsides", "label_adverse_event", "OnSIDES", "v3.1.1 data (exact disease names)", "CC BY 4.0 (cite doi:10.1016/j.medj.2025.100642)",
          "yes", "tested for «may worsen», not released (P56); mapping not checked by hand", load_onsides),
    Layer("chembl_activity", "measured_activity", "ChEMBL", "ChEMBL 37 web services", "CC BY-SA 3.0", "yes",
          "tested as a feature, gate failed (P59)", load_measured("ChEMBL")),
    Layer("gtopdb_activity", "measured_activity", "IUPHAR/BPS Guide to Pharmacology", "2026.3.1", "CC BY-SA 4.0", "yes",
          "tested as a feature, gate failed (P59)", load_measured("GtoPdb")),
    Layer("atc", "drug_atc|atc_parent", "ChEMBL ATC classification", "ChEMBL 37", "ChEMBL CC BY-SA 3.0; ATC © WHO Collaborating Centre",
          "check", "used by R3, rejected (P57)", load_atc),
    Layer("infection_class", "pathogen_class", "MONDO", "releases/2026-10-06", "CC BY 4.0", "yes",
          "used by R3, rejected (P57)", load_infection_class),
    Layer("drugmechdb", "curated_mechanism|mechanism_path", "DrugMechDB", "indication_paths.json", "CC0 1.0", "yes",
          "evaluation only: shared genes overlap curated mechanisms (P50)", load_drugmechdb),
    Layer("lincs", "expression_signature", "LINCS L1000 consensus signatures (Enrichr library)", "Enrichr download",
          "LINCS / Enrichr terms", "check", "tested, negative (P45)", load_lincs),
    Layer("creeds", "expression_signature", "CREEDS disease signatures from GEO (Enrichr library)", "Enrichr download",
          "CREEDS / Enrichr terms", "check", "tested, negative (P45)", load_creeds),
    Layer("qtl_coloc", "qtl_colocalisation", "Open Targets colocalisation", "26.06 (h4 ≥ 0.8, cis, consistent sign)",
          "CC0 1.0", "yes", "used by the ranking model, column «qtl direction» (author's decision, P61)", load_coloc),
    Layer("label_fda", "current_label", "openFDA / DailyMed", "openFDA labels 22-10-2026", "Public domain (US government)",
          "yes", "N4: status part shown (P55); pair labels not released", load_labels(False)),
    Layer("label_ema", "current_label", "EMA product information", "EMA medicines table 22-10-2026",
          "EMA reuse with acknowledgement", "check", "N4: not released (P55)", load_labels(True)),
    Layer("verified_negations", "label_negation_verified", "OREVID check of DailyMed / EMA", "resources/n4_verified_negations.csv",
          "quotes from DailyMed / EMA", "check", "curated list, shown only if the author agrees (P55)", load_verified_negations),
    Layer("prism", "cell_line_sensitivity", "PRISM Repurposing 19Q4 secondary screen", "figshare 10.6084/m9.figshare.9393293",
          "CC BY 4.0", "yes", "O2 (P63)", load_prism),
    Layer("civic", "clinical_evidence", "CIViC", "nightly 23-10-2026", "CC0 1.0", "yes", "O1 oncology exam (P64)", load_civic),
    Layer("repurposing_hub", "annotated_target", "Broad Drug Repurposing Hub", "20200324",
          "CC BY 4.0 on the website; «non-commercial use only» in the file", "check", "P1b (P65)", load_hub),
    Layer("cpic", "pharmacogenomic_pair", "CPIC", "API 23-10-2026", "CC0 1.0", "yes", "X1 display (P67)", load_cpic),
    Layer("pharmgkb", "pharmacogenomic_annotation", "PharmGKB / ClinPGx clinical annotations", "2025-07-05", "CC BY-SA 4.0",
          "yes", "internal graph (X1 uses CPIC for display)", load_pharmgkb),
    Layer("trials", "trial_result", "ClinicalTrials.gov API v2", "23-10-2026 (phase 2–3, posted results)", "Public domain",
          "yes", "T1 validation, claim not supported (P68)", load_trials),
    Layer("ot_approvals", "approval", "Open Targets 26.06 (classified, sources/ot_approvals.py)", "26.06", "CC0 1.0", "yes",
          "evaluation and display of approvals", load_ot_approvals),
]

NODE_ATTRIBUTES = {  # file: (source, licence, redistribution, status)
    "drug_status": ("openFDA, DailyMed, EMA", "US public domain; EMA reuse with acknowledgement", "check", "N4 part b shown (P55)"),
    "orphadata_prevalence": ("Orphadata product 9", "CC BY 4.0", "yes", "R4 display (P66)"),
    "practical": ("Orange Book, Purple Book, EMA", "US public domain; EMA reuse with acknowledgement", "check", "shown (P44)"),
}


# --- build ------------------------------------------------------------------------------------------

class NewNodes:
    def __init__(self, start):
        self.start, self.rows, self.index = start, [], {}

    def __call__(self, node_type, key, name):
        k = (node_type, key)
        if k not in self.index:
            self.index[k] = self.start + len(self.rows)
            self.rows.append({"node_index": self.index[k], "node_id": key, "node_type": node_type, "node_name": name,
                              "node_source": "OREVID layers", "origin": "layers"})
        return self.index[k]


def _slug(s):
    return re.sub(r"[^a-z0-9]+", "_", s.lower()).strip("_")


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    g = _graph()
    dest = OUT / "internal"
    if dest.exists():
        shutil.rmtree(dest)
    (dest / "layers").mkdir(parents=True)
    used = model_relations()
    manifest = {"core_nodes": int(len(g.nodes)), "layers": []}
    for src, rel, d in core_layers():
        name = f"core_{_slug(src)}_{_slug(rel)}"
        subs = sorted(d.display_relation.astype(str).unique())
        part = [s for s in subs if f"{rel}:{s}" in used]
        status = ("used by the ranking model" if rel in used else
                  f"subtypes used by the ranking model: {', '.join(part)}" if part else "in data/prep, not in the ranking model")
        ver, lic, red = CORE[src]
        d.to_parquet(dest / "layers" / f"{name}.parquet", index=False)
        manifest["layers"].append({"name": name, "relation": rel, "source": src, "version": ver, "licence": lic,
                                   "redistribution": red, "core": True, "edges": int(len(d)),
                                   "status": status})
        print(f"  {name}: {len(d):,}", flush=True)
    new = NewNodes(len(g.nodes))
    for L in LAYERS:
        try:
            d = L.load(g, new)
        except FileNotFoundError as e:
            print(f"  {L.name}: missing input ({e})")
            continue
        d.to_parquet(dest / "layers" / f"{L.name}.parquet", index=False)
        manifest["layers"].append({"name": L.name, "relation": L.relation, "source": L.source, "version": L.version,
                                   "licence": L.licence, "redistribution": L.redistribution, "core": False,
                                   "edges": int(len(d)), "status": L.status,
                                   "node_types": sorted({str(t) for t in _types(g, new, d)})})
        print(f"  {L.name}: {len(d):,}", flush=True)
    nodes = pd.concat([g.nodes[["node_index", "node_id", "node_type", "node_name", "node_source", "origin"]],
                       pd.DataFrame(new.rows)], ignore_index=True)
    nodes["node_id"] = nodes.node_id.astype(str)
    nodes.to_parquet(dest / "nodes.parquet", index=False)
    (dest / "node_attributes").mkdir()
    for f, (src, lic, red, status) in NODE_ATTRIBUTES.items():
        shutil.copy(EXT / f"{f}.parquet", dest / "node_attributes" / f"{f}.parquet")
        manifest.setdefault("node_attributes", []).append({"name": f, "source": src, "licence": lic,
                                                           "redistribution": red, "status": status})
    manifest["nodes"] = int(len(nodes))
    manifest["new_nodes"] = pd.DataFrame(new.rows).node_type.value_counts().to_dict() if new.rows else {}
    (dest / "manifest.json").write_text(json.dumps(manifest, indent=1, ensure_ascii=False), "utf-8")
    print(f"internal graph: {len(nodes):,} nodes, {sum(x['edges'] for x in manifest['layers']):,} edges in "
          f"{len(manifest['layers'])} layers -> {dest}")
    report()


def _types(g, new, d):
    t = dict(zip(range(len(g.nodes)), g.types))
    t.update({r["node_index"]: r["node_type"] for r in new.rows})
    return {f"{t.get(int(a))} – {t.get(int(b))}" for a, b in d[["x", "y"]].drop_duplicates().head(20000).itertuples(index=False)}


def export():
    sys.stdout.reconfigure(encoding="utf-8")
    src = OUT / "internal"
    man = json.loads((src / "manifest.json").read_text("utf-8"))
    dest = OUT / "public"
    if dest.exists():
        shutil.rmtree(dest)
    (dest / "layers").mkdir(parents=True)
    keep = [L for L in man["layers"] if L["redistribution"] == "yes"]
    used_nodes = set()
    for L in keep:
        d = pd.read_parquet(src / "layers" / f"{L['name']}.parquet")
        used_nodes |= set(d.x) | set(d.y)
        shutil.copy(src / "layers" / f"{L['name']}.parquet", dest / "layers" / f"{L['name']}.parquet")
    nodes = pd.read_parquet(src / "nodes.parquet")
    nodes[nodes.node_index.isin(used_nodes)].to_parquet(dest / "nodes.parquet", index=False)
    out = {**{k: v for k, v in man.items() if k not in ("layers", "node_attributes")}, "layers": keep,
           "left_out": [{"name": L["name"], "redistribution": L["redistribution"], "licence": L["licence"]}
                        for L in man["layers"] if L["redistribution"] != "yes"]}
    (dest / "manifest.json").write_text(json.dumps(out, indent=1, ensure_ascii=False), "utf-8")
    L = ["# OREVID graph, publishable layers (draft, not published)\n",
         "Each file in `layers/` is distributed under the licence of its source (below); the collection is an aggregation, "
         "not a merged database. `nodes.parquet` lists the nodes the layers use (PrimeKG node indices, then OREVID's).\n",
         "| Layer | Source | Licence | Edges |", "|---|---|---|---|"]
    L += [f"| {x['name']} | {x['source']} | {x['licence']} | {x['edges']:,} |" for x in keep]
    (dest / "README.md").write_text("\n".join(L) + "\n", "utf-8")
    print(f"public graph: {len(keep)} layers, {sum(x['edges'] for x in keep):,} edges; left out {len(out['left_out'])}")


def report():
    man = json.loads((OUT / "internal" / "manifest.json").read_text("utf-8"))
    used = sorted(model_relations())
    L = ["# Graph layers: the internal graph and the publishable one\n",
         "Generated by `py -3 layers.py report` from `data/graph/internal/manifest.json`. The ranking model and the "
         "website read only `data/prep` and the relations listed at the end (frozen in `resources/model_relations.txt`, "
         "checked by `tests/test_model_relations.py`). Every other layer is in the internal graph for exploration and "
         "future uses; a layer enters the model or the website only through a registered test (`docs/IMPROVEMENT_PLAN.md`).\n",
         "Redistribution: **yes** = in the publishable graph (`py -3 layers.py export`, each layer under its own licence); "
         "**check** = licence or terms to confirm with the creators before publishing; **no** = internal only.\n",
         f"Nodes: {man['nodes']:,} ({man['core_nodes']:,} core" + "".join(f", {v:,} {k}" for k, v in man["new_nodes"].items())
         + ").\n", "## Added layers (beyond data/prep)\n",
         "| Layer | Relation | Source (version) | Licence | Redistribution | Edges | Status |", "|---|---|---|---|---|---|---|"]
    for x in man["layers"]:
        if not x["core"]:
            L.append(f"| {x['name']} | {x['relation']} | {x['source']} ({x['version']}) | {x['licence']} | "
                     f"{x['redistribution']} | {x['edges']:,} | {x['status']} |")
    L += ["", "## Core layers (data/prep)\n", "| Layer | Source | Licence | Redistribution | Edges | Status |",
          "|---|---|---|---|---|---|"]
    for x in man["layers"]:
        if x["core"]:
            L.append(f"| {x['name']} | {x['source']} | {x['licence']} | {x['redistribution']} | {x['edges']:,} | {x['status']} |")
    L += ["", "## Node attributes\n", "| File | Source | Licence | Redistribution | Status |", "|---|---|---|---|---|"]
    L += [f"| {x['name']} | {x['source']} | {x['licence']} | {x['redistribution']} | {x['status']} |"
          for x in man.get("node_attributes", [])]
    L += ["", "## Not layers\n",
          "- Network proximity (N2, P51) and RotatE (P47) are computed from the graph, not sources.",
          "- Bioteque (N3): not downloaded (IRB websites offline). Drugs-to-avoid lists (porphyria, Brugada, long QT, "
          "G6PD, myasthenia gravis) and the Broad Repurposing Hub: licences to check, not downloaded.",
          "", "## Relations read by the ranking model\n", ", ".join(f"`{r}`" for r in used)]
    (ROOT / "docs" / "GRAPH_LAYERS.md").write_text("\n".join(L) + "\n", "utf-8")
    print("docs/GRAPH_LAYERS.md written")


if __name__ == "__main__":
    {"build": build, "export": export, "report": report}[sys.argv[1] if len(sys.argv) > 1 else "build"]()
