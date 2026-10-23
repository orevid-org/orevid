"""The ranking model reads a fixed set of relations; new sources go to the internal graph (layers.py), not the model.

    py -3 -m pytest tests -q
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def frozen():
    lines = (ROOT / "resources" / "model_relations.txt").read_text("utf-8").splitlines()
    return {x.strip() for x in lines if x.strip() and not x.startswith("#")}


def test_model_relations_frozen():
    """A relation enters the ranking model only through a registered test: update the file in the same commit."""
    import layers
    assert layers.model_relations() == frozen()


def test_prepare_merges_a_fixed_list():
    """prepare.py merges a fixed list of data/ext files, so a new source cannot shift node indices by accident."""
    import prepare
    assert prepare.CORE_EXT == ["ctd", "direction", "gwas_direction", "opentargets", "pubmed", "pubtator", "string",
                                "structure"]


def test_layer_relations_not_in_core():
    """Added layers use relation names the model cannot read."""
    import layers
    added = {r for L in layers.LAYERS for r in L.relation.split("|")}
    assert not added & frozen()
    assert not added & {re.sub(r":.*", "", r) for r in frozen()}
