"""Pathogen classes of diseases and anti-infective classes of drugs (R3, docs/IMPROVEMENT_PLAN.md, «Registered
23-10-2026: R3»).

- Diseases: the MONDO ids of each graph disease node, their is_a ancestors in the current MONDO release
  (data/raw/mondo/mondo.obo, VERSION file) and four class terms found by name or exact synonym.
- Drugs: ATC codes from ChEMBL (molecule.atc_classifications, cached in data/raw/chembl/atc.json), mapped to graph
  drugs through Open Targets' ChEMBL -> DrugBank cross-references (sources/ids.py); classes by code prefix only.

    py -3 -m sources.infections      # -> data/ext/infections.json (report of classes, drugs and diseases per class)
"""
import json
import re
import sys
import urllib.request
from collections import defaultdict
from functools import cache
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OBO = ROOT / "data" / "raw" / "mondo" / "mondo.obo"
ATC = ROOT / "data" / "raw" / "chembl" / "atc.json"
DEST = ROOT / "data" / "ext" / "infections.json"
CLASS_LABELS = {"bacterial": "bacterial infectious disease", "viral": "viral infectious disease",
                "fungal": "fungal infectious disease", "parasitic": "parasitic infection"}
PREFIXES = {"bacterial": ("J01", "J04"), "fungal": ("J02", "D01"), "viral": ("J05", "D06BB", "S01AD"),
            "parasitic": ("P01", "P02", "P03")}
# disease class -> drug classes that are out of place there (cross-class uses that are standard are allowed)
WRONG = {"viral": {"bacterial", "fungal", "parasitic"}, "bacterial": {"viral", "fungal"},
         "fungal": {"viral", "bacterial"}, "parasitic": {"viral"}}


@cache
def mondo():
    """(is_a parents by MONDO number, class root numbers by class, how each root was found)."""
    parents, names, exact = defaultdict(set), {}, defaultdict(set)
    cur = None
    for line in OBO.open(encoding="utf-8"):
        line = line.rstrip("\n")
        if line.startswith("[Term]"):
            cur = None
        elif line.startswith("id: MONDO:"):
            cur = int(line[10:])
        elif cur is None:
            continue
        elif line.startswith("name: "):
            names[line[6:].lower()] = cur
        elif line.startswith("synonym: ") and " EXACT" in line:
            exact[line.split('"')[1].lower()].add(cur)
        elif line.startswith("is_a: MONDO:"):
            parents[cur].add(int(line[12:19]))
    roots, how = {}, {}
    for cls, lab in CLASS_LABELS.items():
        if lab in names:
            roots[cls], how[cls] = names[lab], "name"
        elif len(exact.get(lab, ())) == 1:
            roots[cls], how[cls] = next(iter(exact[lab])), "exact synonym"
    return parents, roots, how


def mondo_ids(node_id):
    s = str(node_id)
    if re.fullmatch(r"[\d_]+", s):
        return {int(x) for x in s.split("_") if x}
    return {int(m) for m in re.findall(r"MONDO[:_](\d+)", s)}


def disease_class_fn():
    parents, roots, _ = mondo()
    root_cls = {v: k for k, v in roots.items()}

    @cache
    def ancestors_cls(m):
        seen, stack, out = {m}, [m], set()
        while stack:
            x = stack.pop()
            if x in root_cls:
                out.add(root_cls[x])
            for p in parents.get(x, ()):
                if p not in seen:
                    seen.add(p)
                    stack.append(p)
        return frozenset(out)

    def classes(node_id, name=""):
        if "pneumocyst" in str(name).lower():
            return set()
        return set().union(*(ancestors_cls(m) for m in mondo_ids(node_id))) if mondo_ids(node_id) else set()
    return classes


def atc_codes():
    """ChEMBL id -> ATC codes (downloaded once from the ChEMBL web services)."""
    if not ATC.exists():
        ATC.parent.mkdir(parents=True, exist_ok=True)
        base = "https://www.ebi.ac.uk"
        with urllib.request.urlopen(f"{base}/chembl/api/data/status.json", timeout=60) as r:
            version = json.load(r)["chembl_db_version"]
        url, out = (f"{base}/chembl/api/data/molecule.json?atc_classifications__isnull=false"
                    "&only=molecule_chembl_id,atc_classifications&limit=1000"), {}
        while url:
            with urllib.request.urlopen(url, timeout=120) as r:
                j = json.load(r)
            out.update({m["molecule_chembl_id"]: m["atc_classifications"] for m in j["molecules"]})
            url = base + j["page_meta"]["next"] if j["page_meta"]["next"] else None
        ATC.write_text(json.dumps({"version": version, "codes": out}), "utf-8")
    return json.loads(ATC.read_text("utf-8"))


def drug_classes(g):
    """graph drug node -> set of anti-infective classes (from ATC codes)."""
    from sources import ids
    key_node = {str(i): int(n) for n, i, t in zip(g.nodes.node_index, g.nodes.node_id, g.nodes.node_type) if t == "drug"}
    mol = ids.drugs()
    out = defaultdict(set)
    for cid, codes in atc_codes()["codes"].items():
        key = mol.get(cid, (None,))[0]
        f = key_node.get(key[3:]) if key and key.startswith("DB:") else None
        if f is None:
            continue
        for code in codes:
            for cls, pre in PREFIXES.items():
                if code.startswith(pre):
                    out[f].add(cls)
    return dict(out)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from graph import Graph
    g = Graph()
    _, roots, how = mondo()
    dcls = disease_class_fn()
    n = g.nodes
    per_cls, multi = defaultdict(int), 0
    for i, t, nid, name in zip(n.node_index, n.node_type, n.node_id, n.node_name):
        if t != "disease" or g.origin[int(i)] != "PrimeKG":
            continue
        c = dcls(nid, name)
        if len(c) == 1:
            per_cls[next(iter(c))] += 1
        elif len(c) > 1:
            multi += 1
    dc = drug_classes(g)
    drugs_per = defaultdict(int)
    for c in dc.values():
        for k in c:
            drugs_per[k] += 1
    rep = {"mondo": OBO.with_name("VERSION").read_text().strip(), "chembl": atc_codes()["version"],
           "roots": {k: f"MONDO:{v:07d} ({how[k]})" for k, v in roots.items()},
           "missing": [k for k in CLASS_LABELS if k not in roots],
           "website_diseases_one_class": dict(per_cls), "website_diseases_several_classes": multi,
           "drugs_per_class": dict(drugs_per), "drugs_with_a_class": len(dc),
           "drugs_several_classes": sum(len(c) > 1 for c in dc.values())}
    DEST.write_text(json.dumps(rep, indent=1), "utf-8")
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
