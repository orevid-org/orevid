"""Integrity tests: the things that, if they fail, would give absurd results without any warning.

    python -m pytest tests -q

Most tests need the built graph (data/prep, from build_all.sh) and the website data; without them they
are skipped. The benchmark test only needs the committed benchmark files.
"""
import glob
import json
import math
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
WEB = ROOT / "web" / "public" / "data"


@pytest.fixture(scope="module")
def g():
    if not (ROOT / "data" / "prep" / "nodes.pkl").exists():
        pytest.skip("data/prep missing: build the graph first (build_all.sh)")
    from graph import Graph
    return Graph()


# --- filters and rules ---------------------------------------------------------------------

def test_developed_drugs_not_empty(g):
    from predictor import developed_drugs
    s = developed_drugs(g)
    names = {g.name(x).lower() for x in s}
    assert len(s) > 3000, "the approved/investigational filter has become almost empty"
    for known in ["metformin", "ibuprofen", "nivolumab"]:
        assert known in names, f"{known} should be «approved or investigational»"
    assert not any("mannopyranosyl" in n for n in names), "compounds that are only experimental get in"


def test_drug_status(g):
    from predictor import drug_status
    e = {g.name(k).lower(): v for k, v in drug_status(g).items()}
    assert e["metformin"] == "approved"
    assert e["rofecoxib"].startswith("withdrawn"), "a withdrawn drug cannot show up as approved"
    assert e["thalidomide"] == "withdrawn in some market"
    assert e["imidacloprid"] == "veterinary only", "«vet_approved» is not approved in humans"
    assert e["octylphenoxy polyethoxyethanol"] == "experimental"


def test_evidence():
    from web.build import has_evidence
    assert not has_evidence([{"t": "phenotype", "m": [["x", 5]]}]), "phenotypes alone are not evidence"
    assert not has_evidence([{"t": "similar", "m": [["x", 2, ["A", "B"]]]}]), "similar with < 3 genes is not evidence"
    assert not has_evidence([{"t": "interaction", "n": 3, "g": []}, {"t": "similar_drug", "g": []}])
    assert not has_evidence([{"t": "target", "g": [["G", ""]]}]), "a target alone hits 2-3 % (P27)"
    assert not has_evidence([{"t": "neighbour2", "m": [["x", "y"]]}]), "at 2 steps alone it is not confirmed (P27)"
    assert has_evidence([{"t": "similar", "m": [["x", 3, ["A", "B", "C"]]]}])
    assert has_evidence([{"t": "neighbour", "m": ["x"]}])


def test_diagnostics(g):
    from predictor import diagnostic_drugs
    names = {g.name(x).lower() for x in diagnostic_drugs(g)}
    assert "coccidioides immitis spherule" in names and "tuberculin purified protein derivative" in names
    assert not {"cisplatin", "carboplatin", "nitroprusside"} & names, "therapeutic drugs flagged as diagnostic"


def test_drug_family():
    from evaluation.twins import drug_family
    assert drug_family("Teriflunomide") & drug_family("Leflunomide")
    assert drug_family("Dexibuprofen") & drug_family("Ibuprofen")
    assert drug_family("Fosphenytoin") & drug_family("Phenytoin")
    assert not drug_family("Esmolol") & drug_family("Metoprolol")


def test_active_moiety():
    from evaluation.twins import active_moiety
    assert active_moiety("Prednisolone phosphate") == active_moiety("Prednisolone")
    assert active_moiety("Betamethasone sodium phosphate") == "betamethasone"
    assert active_moiety("Dimethyl fumarate") == "dimethyl fumarate"
    assert active_moiety("Sodium phosphate, dibasic") == "sodium phosphate, dibasic"


def test_twins(g):
    from evaluation.twins import twins, names_with_use
    from predictor import KNOWN_USE
    tw = twins(g, "disease")
    d = g.resolve("idiopathic pulmonary fibrosis", "disease")[0]
    assert len(tw[d]) >= 2, "IPF has an Open Targets twin"
    assert "exenatide" in names_with_use(g, KNOWN_USE), "exenatide has a known use (in some twin)"


def test_checkpoint_rule(g):
    from evaluation.direction import checkpoint_rule
    blockers, diseases = checkpoint_rule(g)
    nb = {g.name(x).lower() for x in blockers}
    assert {"nivolumab", "pembrolizumab", "ipilimumab"} <= nb
    assert not ({"abatacept", "belatacept", "peresolimab"} & nb), "agonists or fusion proteins flagged as blockers"
    nm = [g.name(x).lower() for x in diseases]
    assert any(n == "rheumatoid arthritis" for n in nm)
    for word in ["lymphoma", "leukemia", "myeloma", "thymoma"]:
        assert not any(word in n for n in nm), f"the immune rule includes a cancer ({word})"


# --- web data ----------------------------------------------------------------------

def _web_diseases():
    for f in glob.glob(str(WEB / "d" / "*.json")):
        yield from json.loads(Path(f).read_text("utf-8")).values()


def _no_nan(x):
    if isinstance(x, float):
        return not math.isnan(x)
    if isinstance(x, dict):
        return all(_no_nan(v) for v in x.values())
    if isinstance(x, list):
        return all(_no_nan(v) for v in x)
    return True


@pytest.mark.skipif(not (WEB / "index.json").exists(), reason="no web data")
def test_web_data(g):
    from predictor import KNOWN_USE, developed_drugs, diagnostic_drugs
    from evaluation.twins import key, names_with_use
    developed = developed_drugs(g) - diagnostic_drugs(g)
    degree = np.diff(g.M(KNOWN_USE).indptr)
    names_use = names_with_use(g, KNOWN_USE)
    n = 0
    for d in _web_diseases():
        n += 1
        assert _no_nan(d), f"{d['name']}: there are NaNs"
        from evaluation.twins import drug_family
        approved = {x for a in d["approved"] + d.get("approved_ot", []) for x in drug_family(a)}
        from web.build import has_evidence
        for c in d["candidates"] + d["extensions"]:
            assert has_evidence(c["expl"]), f"{d['name']}: {c['name']} is shown without enough evidence"
        for c in d["no_use"]:
            assert any(e["t"] == "target" for e in c["expl"]), f"{d['name']}: {c['name']} (exploratory) without a target"
        for c in d["candidates"] + d["extensions"]:
            assert not drug_family(c["name"]) & approved, f"{d['name']}: {c['name']} (or a salt, prodrug or enantiomer) is already approved"
            assert c["name"].lower() not in {x.lower() for x in d["contra"]}, f"{d['name']}: {c['name']} contraindicated here"
            assert c["status"] != "veterinary only", f"{d['name']}: {c['name']} is veterinary only"
            assert "dir_p" not in c, "no direction probability must be published"
            if c["risk"]:
                assert c["risk_reasons"], f"{d['name']}: {c['name']} with a risk flag without a reason"
        for c in d["no_use"]:
            assert degree[c["node"]] == 0, f"{d['name']}: {c['name']} has some known use and shows up as «no use»"
            assert key(c["name"]) not in names_use, f"{d['name']}: {c['name']} has a twin with a known use"
            assert c["node"] in developed, f"{d['name']}: {c['name']} is neither approved nor investigational, or is diagnostic"
            assert not c["risk"], f"{d['name']}: {c['name']} with a risk flag should not be shown"
    assert n > 0
    index = json.loads((WEB / "index.json").read_text("utf-8"))
    assert len(index["diseases"]) > 10000


# --- benchmark -----------------------------------------------------------------------

def test_benchmark_reproduces():
    t = json.loads((ROOT / "benchmark" / "tasks" / "zero_shot_1.json").read_text())
    assert t["n_test_diseases"] == 71
    assert sum(len(r["positives"]) for r in t["diseases"]) == 444
    assert sum(r["near_duplicate_overlap"] < 0.5 for r in t["diseases"]) == 29
    res = json.loads((ROOT / "benchmark" / "predictions" / "txgnn" / "results.json").read_text())
    rz = np.mean([res[f"zero_shot_{s}"]["all"]["Recall@100"] for s in range(1, 6)])
    assert abs(rz - 0.412) < 0.002, "TxGNN no longer reproduces the reference value"
