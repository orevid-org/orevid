"""Disease hierarchy used for the website's navigation (sources/hierarchy.py; REVIEW_RESPONSE.md T6).

    py -3 -m pytest tests -q
"""
from collections import defaultdict
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
HIERARCHY = ROOT / "data" / "prep" / "disease_hierarchy.parquet"
BREAST_CANCER, TRIPLE_NEGATIVE, HER2_POSITIVE = 36691, 37655, 84130


@pytest.fixture(scope="module")
def children():
    if not HIERARCHY.exists():
        pytest.skip("data/prep/disease_hierarchy.parquet missing (py -3 -m sources.hierarchy)")
    import pandas as pd
    h = pd.read_parquet(HIERARCHY)
    ch = defaultdict(set)
    for a, b in zip(h.parent.astype(int), h.child.astype(int)):
        ch[a].add(b)
    return ch


def descendants(ch, d):
    seen, stack = set(), list(ch[d])
    while stack:
        m = stack.pop()
        if m not in seen:
            seen.add(m)
            stack.extend(ch[m])
    return seen


def test_breast_cancer_subtypes_are_descendants(children):
    below = descendants(children, BREAST_CANCER)
    assert TRIPLE_NEGATIVE in below and HER2_POSITIVE in below


def test_no_cycles(children):
    state = {}

    def visit(u):  # iterative DFS: 1 = on the stack, 2 = done
        stack = [(u, iter(children[u]))]
        state[u] = 1
        while stack:
            node, it = stack[-1]
            nxt = next(it, None)
            if nxt is None:
                state[node] = 2
                stack.pop()
            elif state.get(nxt) == 1:
                return True
            elif nxt not in state:
                state[nxt] = 1
                stack.append((nxt, iter(children[nxt])))
        return False

    assert not any(visit(u) for u in list(children) if u not in state)


def test_pages_list_parents_and_children():
    """Every disease page lists as «broader» each parent that lists it as «narrower»."""
    import json
    files = sorted((ROOT / "web" / "public" / "data" / "d").glob("*.json"))
    if not files:
        pytest.skip("website data missing")
    pages = {}
    for f in files:
        pages.update(json.loads(f.read_text("utf-8")))
    with_h = {int(k): v for k, v in pages.items() if "narrower" in v}
    if not with_h:
        pytest.skip("website data built before the hierarchy (T6)")
    for d, p in with_h.items():
        for child, _, _ in p["narrower"]:
            c = with_h.get(int(child))
            if c is not None:
                assert d in {int(x[0]) for x in c["broader"]}
