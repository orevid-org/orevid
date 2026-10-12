"""Relations extracted from PubMed with Claude, each with its citation.

    py -3 -m sources.pubmed candidates "multiple sclerosis" [--top 15] [--per 5]
        For each candidate from the graph, searches "drug + disease" articles and extracts relations from them.
    py -3 -m sources.pubmed query "multiple sclerosis[MeSH] AND drug therapy[sh]" [--n 50]
        Any PubMed search.
    py -3 -m sources.pubmed import file.jsonl
        Imports extractions made elsewhere (same format as Claude's).
    py -3 -m sources.pubmed link
        Links the extracted relations to the graph nodes -> data/ext/pubmed_edges.parquet
        (prepare.py must be run again afterwards).

Needs an Anthropic API key in ANTHROPIC_API_KEY for "candidates" and "query".

Each relation carries the PMID and the literal sentence from the article. Relations whose
sentence does not appear in the abstract are discarded: the citation must be verifiable.

Caution: these relations do NOT enter the model or the validation. "This drug treats
this disease" is exactly the answer the graph tries to guess; using it as a
signal would inflate the results without helping to discover anything. They serve to know
what the literature says about each candidate (for, against, nothing).
"""
import argparse
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd

from sources import ids

DIR = ids.RAW / "pubmed"
ARTICLES = DIR / "articles"
EXTRACTED = DIR / "extracted"
EXT = Path(__file__).parent.parent / "data" / "ext"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MODEL = "claude-opus-5-5"

RELATIONS = ["treats", "no_improvement", "worsens", "increases", "decreases", "inhibits", "activates",
             "associated", "adverse_effect"]
TYPES = ["drug", "gene", "disease", "phenotype", "other"]
EVIDENCE = ["randomised_trial", "clinical_trial", "observational", "case_report", "animal",
            "in_vitro", "review", "hypothesis"]

SCHEMA = {
    "type": "object",
    "properties": {
        "relations": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "pmid": {"type": "string"},
                    "subject": {"type": "string"},
                    "subject_type": {"type": "string", "enum": TYPES},
                    "relation": {"type": "string", "enum": RELATIONS},
                    "object": {"type": "string"},
                    "object_type": {"type": "string", "enum": TYPES},
                    "evidence": {"type": "string", "enum": EVIDENCE},
                    "sentence": {"type": "string"},
                },
                "required": ["pmid", "subject", "subject_type", "relation", "object", "object_type",
                             "evidence", "sentence"],
                "additionalProperties": False,
            },
        }
    },
    "required": ["relations"],
    "additionalProperties": False,
}

INSTRUCTIONS = """You extract biomedical relations from PubMed abstracts for a drug-repurposing knowledge graph.

For each abstract, list the relations that the abstract states as findings (or, for reviews and hypotheses, as claims), between drugs, genes/proteins, diseases and phenotypes:
- treats: the drug improved / treated the disease or phenotype.
- no_improvement: the drug was tested and did NOT improve it (negative or null result). Do not drop these: negative results are as valuable as positive ones.
- worsens: the drug worsened it or triggered it.
- increases / decreases: the drug increased / decreased the level, expression or activity of a gene or protein.
- inhibits / activates: the drug is described as an inhibitor / activator (agonist) of a gene or protein.
- associated: a gene is associated with a disease (genetics, expression, biomarker).
- adverse_effect: the drug caused an adverse effect.

Rules:
- Use canonical English names: generic drug names (not brands), HGNC gene symbols, standard disease names.
- `sentence` must be one sentence copied verbatim from the abstract that supports the relation. Do not paraphrase.
- `evidence`: the kind of study the finding comes from. Use `hypothesis` when the abstract only proposes or speculates.
- Only include relations about the entities' effects, not background statements about other work, unless the abstract is a review.
- If an abstract has no relevant relation, return none for it."""


# --- PubMed ---------------------------------------------------------------

def _get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "orevid/0.1"})
    with urllib.request.urlopen(req) as r:
        data = r.read()
    time.sleep(0.4)  # NCBI: at most 3 requests per second without a key
    return data


def search(query, n):
    q = urllib.parse.urlencode({"db": "pubmed", "term": query, "retmax": n, "sort": "relevance", "retmode": "json"})
    return json.loads(_get(f"{EUTILS}/esearch.fcgi?{q}"))["esearchresult"]["idlist"]


def _cached(pmid):
    """True if the article is in the cache with the current (English) keys."""
    f = ARTICLES / f"{pmid}.json"
    return f.exists() and "abstract" in json.loads(f.read_text("utf-8"))


def articles(pmids):
    """Title, abstract and year of each PMID (cached in data/raw/pubmed/articles)."""
    ARTICLES.mkdir(parents=True, exist_ok=True)
    missing = [p for p in pmids if not _cached(p)]
    for i in range(0, len(missing), 100):
        q = urllib.parse.urlencode({"db": "pubmed", "id": ",".join(missing[i:i + 100]), "retmode": "xml"})
        root = ET.fromstring(_get(f"{EUTILS}/efetch.fcgi?{q}"))
        for art in root.iter("PubmedArticle"):
            pmid = art.findtext(".//MedlineCitation/PMID")
            parts = []
            for t in art.iterfind(".//Abstract/AbstractText"):
                text = "".join(t.itertext()).strip()
                parts.append(f"{t.get('Label')}: {text}" if t.get("Label") else text)
            year = art.findtext(".//Article/Journal/JournalIssue/PubDate/Year") or \
                (art.findtext(".//Article/Journal/JournalIssue/PubDate/MedlineDate") or "")[:4]
            data = {"pmid": pmid, "title": "".join(art.find(".//ArticleTitle").itertext()).strip(),
                    "abstract": "\n".join(parts), "year": int(year) if year.isdigit() else None}
            (ARTICLES / f"{pmid}.json").write_text(json.dumps(data, ensure_ascii=False), "utf-8")
    res = []
    for p in pmids:
        f = ARTICLES / f"{p}.json"
        if f.exists():
            res.append(json.loads(f.read_text("utf-8")))
    return res


# --- extraction -----------------------------------------------------------

def _normalise(text):
    return re.sub(r"\s+", " ", text).strip().lower()


def verify(relations, by_pmid):
    """Keeps the relations whose sentence is literally in the title or the abstract."""
    good = []
    for r in relations:
        a = by_pmid.get(r["pmid"])
        if a and _normalise(r["sentence"]).rstrip(".") in _normalise(a["title"] + " " + a["abstract"]):
            good.append(r)
    return good


def extract_batch(client, batch):
    import anthropic

    text = "\n\n".join(f'<article pmid="{a["pmid"]}">\nTitle: {a["title"]}\n{a["abstract"]}\n</article>' for a in batch)
    try:
        resp = client.beta.messages.create(
            model=MODEL,
            max_tokens=16000,
            betas=["server-side-fallback-2026-07-01"],
            fallbacks="default",
            output_config={"effort": "medium", "format": {"type": "json_schema", "schema": SCHEMA}},
            system=INSTRUCTIONS,
            messages=[{"role": "user", "content": text}],
        )
    except anthropic.RateLimitError:
        time.sleep(30)
        return extract_batch(client, batch)
    except anthropic.APIStatusError as e:
        print(f"  API error ({e.status_code}): {e.message}")
        return None
    if resp.stop_reason == "refusal":
        print("  the model declined this batch")
        return None
    if resp.stop_reason == "max_tokens":
        print("  truncated response (max_tokens): the batch is split")
        if len(batch) > 1:
            m = len(batch) // 2
            a, b = extract_batch(client, batch[:m]), extract_batch(client, batch[m:])
            return (a or []) + (b or [])
        return None
    txt = "".join(b.text for b in resp.content if b.type == "text")
    return json.loads(txt)["relations"]


def extract(arts, batch_size=8):
    """Extracts relations from the articles that have not been processed yet."""
    import anthropic

    EXTRACTED.mkdir(parents=True, exist_ok=True)
    pending = [a for a in arts if a["abstract"] and not (EXTRACTED / f"{a['pmid']}.json").exists()]
    if not pending:
        return
    client = anthropic.Anthropic()
    for i in range(0, len(pending), batch_size):
        batch = pending[i:i + batch_size]
        rels = extract_batch(client, batch)
        if rels is None:
            continue
        by_pmid = {a["pmid"]: a for a in batch}
        good = verify(rels, by_pmid)
        for a in batch:
            save(a["pmid"], [r for r in good if r["pmid"] == a["pmid"]], f"{MODEL} (API)")
        print(f"  {min(i + batch_size, len(pending))}/{len(pending)} articles: "
              f"{len(good)} relations verified out of {len(rels)}", flush=True)


def save(pmid, relations, extractor):
    (EXTRACTED / f"{pmid}.json").write_text(json.dumps(
        {"pmid": pmid, "extractor": extractor, "relations": relations}, ensure_ascii=False, indent=1), "utf-8")


# --- linking to the graph nodes -------------------------------------------

def name_index():
    """Lower-case name -> canonical key, per node type."""
    nodes = pd.read_pickle(ids.RAW.parent / "prep" / "nodes.pkl")
    idx = {t: {} for t in TYPES}
    graph_types = {"drug": "drug", "gene/protein": "gene", "disease": "disease", "effect/phenotype": "phenotype"}
    for t, nid, name, origin in zip(nodes.node_type, nodes.node_id, nodes.node_name, nodes.origin):
        tt = graph_types.get(t)
        if tt is None:
            continue
        key = ids.primekg_keys(t, nid)[0] if origin == "PrimeKG" else nid
        idx[tt].setdefault(str(name).lower(), key)
    # synonyms of drugs (Open Targets) and of diseases
    m = pd.read_parquet(ids.OT / "drug_molecule", columns=["id", "synonyms", "tradeNames"])
    drugs = ids.drugs()
    for c, syn, trade in zip(m.id, m.synonyms, m.tradeNames):
        if c in drugs:
            for s in [*(syn if syn is not None else []), *(trade if trade is not None else [])]:
                idx["drug"].setdefault(s["label"].lower(), drugs[c][0])
    mapping, _, _, _ = ids.diseases()
    d = pd.read_parquet(ids.OT / "disease", columns=["id", "exactSynonyms"])
    for did, syn in zip(d.id, d.exactSynonyms):
        if did in mapping and syn is not None:
            for s in syn:
                idx["disease"].setdefault(s.lower(), mapping[did])
    idx["phenotype"].update({k: v for k, v in idx["disease"].items() if k not in idx["phenotype"]})
    return idx


def link():
    idx = name_index()
    rows, unlinked = [], 0
    for f in sorted(EXTRACTED.glob("*.json")):
        e = json.loads(f.read_text("utf-8"))
        art = ARTICLES / f"{e['pmid']}.json"
        year = json.loads(art.read_text("utf-8")).get("year") if art.exists() else None
        for r in e["relations"]:
            x = idx.get(r["subject_type"], {}).get(r["subject"].lower())
            y = idx.get(r["object_type"], {}).get(r["object"].lower())
            if not (x and y):
                unlinked += 1
                continue
            rows.append({"x": x, "y": y, "relation": "pubmed", "display_relation": r["relation"],
                         "pmid": e["pmid"], "sentence": r["sentence"], "evidence": r["evidence"],
                         "year": year, "extractor": e["extractor"], "source": "PubMed"})
    a = pd.DataFrame(rows)
    EXT.mkdir(parents=True, exist_ok=True)
    a.to_parquet(EXT / "pubmed_edges.parquet", index=False)
    print(f"Relations linked to the graph: {len(a):,}; without a node: {unlinked:,}")
    if len(a):
        print(a.display_relation.value_counts().to_string())
    print("Now run: py -3 prepare.py")


# --- commands -------------------------------------------------------------

def candidates(a):
    sys.path.insert(0, str(Path(__file__).parent.parent))
    from graph import Graph

    g = Graph()
    found = g.resolve(a.disease, "disease")
    if len(found) != 1:
        raise SystemExit(f"\"{a.disease}\": {len(found)} matches; give an exact name.")
    d = found[0]
    model = g.model()
    scores = g.scores([d], model["paths"])
    total = g.apply_model(scores, model)[:, 0]
    already = set(g.neighbors(d, "indication")) | set(g.neighbors(d, "clinical:approved"))
    order = [g.drugs[j] for j in total.argsort()[::-1] if g.drugs[j] not in already][:a.top]
    pmids = []
    for f in order:
        q = f'"{g.name(f)}"[tiab] AND "{g.name(d)}"[tiab]'
        found = search(q, a.per)
        print(f"  {g.name(f)}: {len(found)} articles")
        pmids += found
    extract(articles(list(dict.fromkeys(pmids))))
    link()


def query(a):
    pmids = search(a.query, a.n)
    print(f"{len(pmids)} articles")
    extract(articles(pmids))
    link()


def import_file(a):
    """JSONL: one relation per line with the schema fields (plus an optional "extractor")."""
    rows = [json.loads(l) for l in Path(a.file).read_text("utf-8").splitlines() if l.strip()]
    by_pmid = {}
    for r in rows:
        by_pmid.setdefault(r["pmid"], []).append(r)
    arts = {x["pmid"]: x for x in articles(list(by_pmid))}
    EXTRACTED.mkdir(parents=True, exist_ok=True)
    total = 0
    for pmid, rels in by_pmid.items():
        good = verify(rels, arts)
        total += len(good)
        save(pmid, [{k: r[k] for k in SCHEMA["properties"]["relations"]["items"]["required"]} for r in good],
             rels[0].get("extractor", "imported"))
    print(f"Imported {total} of {len(rows)} relations (the others did not have the sentence in the abstract)")
    link()


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    p = argparse.ArgumentParser(description="PubMed relations with citations")
    s = p.add_subparsers(dest="command", required=True)
    c = s.add_parser("candidates"); c.add_argument("disease"); c.add_argument("--top", type=int, default=15)
    c.add_argument("--per", type=int, default=5)
    c = s.add_parser("query"); c.add_argument("query"); c.add_argument("--n", type=int, default=50)
    c = s.add_parser("import"); c.add_argument("file")
    s.add_parser("link")
    a = p.parse_args()
    {"candidates": candidates, "query": query, "import": import_file, "link": lambda _: link()}[a.command](a)


if __name__ == "__main__":
    main()
