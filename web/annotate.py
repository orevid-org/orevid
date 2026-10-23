"""Adds information that does not change the ranking to the pages already built (no model run): G2 genetic direction,
X1 pharmacogenomic guidelines, R4 prevalence (docs/IMPROVEMENT_PLAN.md, «Registered 23-10-2026»). Each is applied only
once its release check has passed (RELEASED); `--dump` writes every line that would be shown, for the checks.

    py -3 -m web.annotate --dump     # -> data/eval/annotations.jsonl (all candidate lines, nothing changed)
    py -3 -m web.annotate            # writes the RELEASED fields into web/public/data/d/*.json
"""
import json
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
DEST = ROOT / "web" / "public" / "data" / "d"
RELEASED = {"pgx", "prevalence", "genetics"}  # release checks passed: X1 0/50, R4b 1/50, G2c 1/50 (P66, P67, P72)
LISTS = ("candidates", "no_use", "extensions")


class Annotations:
    def __init__(self):
        from evaluation.m1_direction import drug_actions, gene_disease_signs
        from evaluation.twins import twins
        from graph import Graph
        from sources import ids
        g = self.g = Graph()
        tw = twins(g, "disease")
        self.sign, self.qtl, self.ot_ids = {}, defaultdict(set), defaultdict(set)
        combined, ot_ids, _, _ = gene_disease_signs(g, ["pqtl", "eqtl"])
        for (gene, d), s in combined.items():
            for d2 in tw.get(d, {d}):
                self.sign[(gene, d2)] = s
                self.ot_ids[(gene, d2)] |= ot_ids[(gene, d)]
        for q in ("pqtl", "eqtl"):
            signs, _, _, _ = gene_disease_signs(g, [q])
            for (gene, d), s in signs.items():
                for d2 in tw.get(d, {d}):
                    self.qtl[(gene, d2)].add(q)
        self.acts = drug_actions(g)
        ncbi_ens = {v[0]: e for e, v in ids.genes().items()}
        n = g.nodes
        self.ensg = {int(k): ncbi_ens.get(f"NCBI:{i}", "") for k, i, t in zip(n.node_index, n.node_id, n.node_type)
                     if t == "gene/protein"}
        c = pd.read_parquet(ROOT / "data" / "ext" / "cpic.parquet")
        c = c[c.cpic_level.isin(["A", "B", "A/B"])]
        self.pgx = defaultdict(list)
        for r in c.itertuples():
            self.pgx[int(r.drug)].append({"g": r.gene_symbol, "l": r.cpic_level, "u": r.url})
        # G2c (P72): only lines whose colocalisation is stable in Open Targets 26.09 and rests on a case–control GWAS of
        # the page's disease (`py -3 -m evaluation.release_checks g2b` writes them; new candidates need that step again)
        v = ROOT / "data" / "ext" / "genetic_lines_g2c.parquet"
        self.verified = {(int(r.disease), int(r.drug), r.gene): (r.gwas_study, r.trait)
                         for r in pd.read_parquet(v).itertuples()} if v.exists() else {}
        p = pd.read_parquet(ROOT / "data" / "ext" / "orphadata_prevalence.parquet")
        self.prev = {int(r.disease): {"c": r.prevalence, "o": r.orpha, "n": r.orpha_name} for r in p.itertuples()}

    def genetics(self, f, d):
        out = []
        for gene, a in sorted(self.acts.get(int(f), {}).items()):
            s = self.sign.get((gene, int(d)))
            if s is not None:
                out.append({"g": self.g.name(gene), "e": self.ensg.get(gene, ""), "s": s, "a": a, "agree": s * -a > 0,
                            "q": "+".join(sorted(x.replace("qtl", "QTL") for x in self.qtl[(gene, int(d))])),
                            "ot": sorted(self.ot_ids[(gene, int(d))])})
        return out


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    A = Annotations()
    dump = "--dump" in sys.argv
    rows, n_pages = [], 0
    for f in sorted(DEST.glob("*.json")):
        pages = json.loads(f.read_text("utf-8"))
        for k, page in pages.items():
            d = int(k)
            n_pages += 1
            if d in A.prev:
                if dump:
                    rows.append({"kind": "prevalence", "disease": d, "disease_name": page["name"], **A.prev[d]})
                if "prevalence" in RELEASED:
                    page["prevalence"] = A.prev[d]
            for lst in LISTS:
                for c in page.get(lst) or []:
                    gl = A.genetics(c["node"], d)
                    if dump:
                        rows += [{"kind": "genetics", "disease": d, "disease_name": page["name"], "drug": c["node"],
                                  "drug_name": c["name"], "list": lst, **x} for x in gl]
                        rows += [{"kind": "pgx", "drug": c["node"], "drug_name": c["name"], **x} for x in A.pgx.get(c["node"], [])]
                    shown = []
                    for x in gl:
                        src = A.verified.get((d, int(c["node"]), x["g"]))
                        if src:
                            shown.append({k: x[k] for k in ("g", "e", "s", "a", "agree", "q")} | {"st": src[0], "tr": src[1]})
                    if "genetics" in RELEASED and shown:
                        c["genetics"] = shown
                    if "pgx" in RELEASED and A.pgx.get(c["node"]):
                        c["pgx"] = A.pgx[c["node"]]
        if not dump and RELEASED:
            f.write_text(json.dumps(pages, ensure_ascii=False, separators=(",", ":"), allow_nan=False), "utf-8")
    if dump:
        out = ROOT / "data" / "eval" / "annotations.jsonl"
        with out.open("w", encoding="utf-8") as fh:
            for r in rows:
                fh.write(json.dumps(r, ensure_ascii=False) + "\n")
        df = pd.DataFrame(rows)
        print(f"{n_pages:,} pages; lines that would be shown: " + ", ".join(
            f"{k} {len(v):,} ({v.disease.nunique() if 'disease' in v and v.disease.notna().any() else v.drug.nunique():,} "
            f"{'pages' if k != 'pgx' else 'drugs'})" for k, v in df.groupby("kind")) + f" -> {out}")
    else:
        print(f"{n_pages:,} pages annotated with {sorted(RELEASED) or 'nothing'}")


if __name__ == "__main__":
    main()
