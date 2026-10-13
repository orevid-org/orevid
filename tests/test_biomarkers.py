"""Biomarker-defined subtypes (resources/biomarker_subtypes.csv; web/build.py `marker_rules`; REVIEW_RESPONSE.md T7).

    py -3 -m pytest tests -q
"""
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
WEB = ROOT / "web" / "public" / "data" / "d"
TRIPLE_NEGATIVE, ER_NEGATIVE, HER2_POSITIVE = 37655, 37773, 84130
ENDOCRINE = {"Tamoxifen", "Formestane", "Testolactone"}


@pytest.fixture(scope="module")
def pages():
    files = sorted(WEB.glob("*.json"))
    if not files:
        pytest.skip("website data missing")
    res = {}
    for f in files:
        res.update({int(k): v for k, v in json.loads(f.read_text("utf-8")).items()})
    if "markers" not in res.get(TRIPLE_NEGATIVE, {}):
        pytest.skip("website data built before the biomarker rule (T7)")
    return res


def shown(page):
    names = set()
    for k in ("candidates", "extensions"):
        for c in page[k]:
            names |= {c["name"], *c.get("also", [])}
    return names


@pytest.mark.parametrize("d", [TRIPLE_NEGATIVE, ER_NEGATIVE])
def test_no_endocrine_therapy_without_oestrogen_receptor(pages, d):
    assert not (shown(pages[d]) & ENDOCRINE)


def test_trastuzumab_still_approved_for_her2_positive(pages):
    assert "Trastuzumab" in pages[HER2_POSITIVE]["approved"]


def test_table_has_no_duplicates():
    import pandas as pd
    t = pd.read_csv(ROOT / "resources" / "biomarker_subtypes.csv")
    assert not t.duplicated(["node", "marker"]).any()
    assert set(t.status) <= {"absent", "present"}
