"""The shared-genes threshold must be the same in the page builder, the evaluation and the frontend (P40).

    py -3 -m pytest tests -q
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def test_same_threshold_everywhere():
    build = re.search(r"^MAX_P_SHARED = ([0-9e.-]+)", (ROOT / "web" / "build.py").read_text("utf-8"), re.M).group(1)
    evidence = re.search(r"^MAX_P = ([0-9e.-]+)", (ROOT / "evaluation" / "evidence.py").read_text("utf-8"), re.M).group(1)
    app = re.search(r"const validSimilar = .*?m\[3\] < ([0-9e.-]+)\)", (ROOT / "web" / "public" / "app.js").read_text("utf-8")).group(1)
    assert float(build) == float(evidence) == float(app)


def test_same_hub_threshold():
    """R1b (P48): the website and the evaluation leave out the same hub genes."""
    build = re.search(r"^HUB_DISEASES = (\d+)", (ROOT / "web" / "build.py").read_text("utf-8"), re.M).group(1)
    ci = re.search(r"^HUB_DISEASES = (\d+)", (ROOT / "evaluation" / "evidence_ci.py").read_text("utf-8"), re.M).group(1)
    assert int(build) == int(ci) == 500
