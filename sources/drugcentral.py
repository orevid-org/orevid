"""DrugCentral (CC BY-SA 4.0; PostgreSQL dump of 02-09-2026 in data/raw/drugcentral/) for L1 (docs/IMPROVEMENT_PLAN.md,
«Registered 24-10-2026: L1–L6»). The dump is read as text: only the COPY blocks of the tables needed.

- drug–disease relations (`omop_relationship`: indication, contraindication, off-label use, …), diseases by UMLS CUI /
  SNOMED CT -> MONDO (SSSOM exact) -> graph node;
- first approval per drug and agency (`approval`, with the orphan flag);
- target development level per gene (`target_component.tdl`, from TCRD / Pharos).
Drugs -> graph by DrugBank id (`identifier`), else by name (moieties, synonyms).

    py -3 -m sources.drugcentral    # -> data/ext/drugcentral_relations.parquet, _approvals.parquet, _tdl.parquet
"""
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw" / "drugcentral" / "drugcentral.dump.09022026.sql.gz"
EXT = ROOT / "data" / "ext"
TABLES = {"omop_relationship", "identifier", "approval", "structures", "target_component", "synonyms"}


def read_tables():
    out, cur, cols = defaultdict(list), None, None
    with gzip.open(RAW, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if cur is None:
                if line.startswith("COPY public."):
                    name = line[12:line.index(" (")]
                    if name in TABLES:
                        cur, cols = name, [c.strip() for c in line[line.index("(") + 1:line.index(")")].split(",")]
                        out[cur + "#cols"] = cols
                continue
            if line.startswith("\\."):
                cur = None
                continue
            vals = line.rstrip("\n").split("\t")
            if cur == "structures":  # only id and name (molfiles are large)
                vals = [vals[cols.index("id")], vals[cols.index("name")]]
            out[cur].append([None if v == "\\N" else v for v in vals])
    return {t: pd.DataFrame(out[t], columns=(["id", "name"] if t == "structures" else out[t + "#cols"])) for t in TABLES}


def build():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    from sources.ids import xref_to_mondo
    from sources.practical import moiety
    from sources.prism import drug_index
    g = Graph()
    T = read_tables()
    n = g.nodes
    db_node = {str(i): int(k) for k, i, t, s in zip(n.node_index, n.node_id, n.node_type, n.node_source)
               if t == "drug" and s == "DrugBank"}
    by_name = drug_index(g)
    ident = T["identifier"]
    struct_db = defaultdict(set)
    for sid, idv, typ in zip(ident.struct_id, ident.identifier, ident.id_type):
        if typ == "DRUGBANK_ID" and idv in db_node:
            struct_db[sid].add(db_node[idv])
    names = dict(zip(T["structures"].id, T["structures"].name))
    syn = T["synonyms"].dropna(subset=["id"]).groupby("id").name.apply(list).to_dict()

    def drugs_of(sid):
        fs = struct_db.get(sid)
        if fs:
            return fs, "DrugBank id"
        for nm in [names.get(sid)] + syn.get(sid, []):
            if nm and by_name.get(moiety(nm)):
                return by_name[moiety(nm)], "name"
        return set(), "none"
    to_mondo = xref_to_mondo()
    mondo_node = {}
    for k, i, t in zip(n.node_index, n.node_id, n.node_type):
        if t == "disease" and str(i).replace("_", "").isdigit():
            for p in str(i).split("_"):
                mondo_node.setdefault(f"MONDO:{int(p)}", int(k))
    o = T["omop_relationship"]
    rows, how = [], defaultdict(int)
    cache = {}
    for sid, rel, cname, cui, sct in zip(o.struct_id, o.relationship_name, o.concept_name, o.umls_cui, o.snomed_conceptid):
        if sid not in cache:
            cache[sid] = drugs_of(sid)
        fs, h = cache[sid]
        keys = ([f"UMLS:{cui}"] if cui else []) + ([f"SCTID:{sct}", f"SNOMEDCT:{sct}"] if sct else [])
        ds = {mondo_node[m] for k in keys for m in to_mondo.get(k, ()) if m in mondo_node}
        how[(h, bool(ds))] += 1
        for f in fs:
            for d in ds:
                rows.append({"drug": f, "disease": d, "relation": rel, "concept": cname, "struct_id": sid})
    rel = pd.DataFrame(rows).drop_duplicates(["drug", "disease", "relation"])
    rel.to_parquet(EXT / "drugcentral_relations.parquet", index=False)
    a = T["approval"].dropna(subset=["approval"])
    arows = []
    for sid, date, typ, orphan in zip(a.struct_id, a.approval, a.type, a.orphan):
        if sid not in cache:
            cache[sid] = drugs_of(sid)
        for f in cache[sid][0]:
            arows.append({"drug": f, "date": date, "agency": typ, "orphan": orphan == "t"})
    ap = pd.DataFrame(arows)
    ap.to_parquet(EXT / "drugcentral_approvals.parquet", index=False)
    tc = T["target_component"]
    tdl = tc.dropna(subset=["gene", "tdl"])[["gene", "tdl", "accession"]].drop_duplicates("gene")
    tdl.to_parquet(EXT / "drugcentral_tdl.parquet", index=False)
    print(f"DrugCentral: {len(o):,} drug–concept relations; mapping (drug how, disease mapped): "
          + ", ".join(f"{k}: {v:,}" for k, v in sorted(how.items(), key=str)) + f"; graph pairs {len(rel):,} "
          f"({rel.relation.value_counts().to_dict()}); approvals {len(ap):,} for {ap.drug.nunique():,} drugs; TDL for "
          f"{len(tdl):,} genes")


if __name__ == "__main__":
    build()
