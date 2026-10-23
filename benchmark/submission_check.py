"""Evaluates a submission (SUBMITTING.md): reads submissions/<model>/submission.json, downloads the zip, checks its MD5,
extracts only the expected task files, then runs bench.py validate and evaluate. Used by
.github/workflows/bench-submission.yml and by maintainers.

    python submission_check.py submissions/<model>/submission.json [--out check_out]
"""
import argparse
import hashlib
import json
import re
import sys
import urllib.request
import zipfile
from pathlib import Path

from bench import HERE, evaluate, summary, validate

REQUIRED = ["model", "team", "contact", "date", "code_url", "predictions_doi", "predictions_url", "predictions_md5",
            "data_beyond_training_edges", "settings_chosen_on", "external_exam_eligible"]
MAX_BYTES = 2 << 30


def main():
    p = argparse.ArgumentParser()
    p.add_argument("submission")
    p.add_argument("--out", default="check_out")
    a = p.parse_args()
    sub = Path(a.submission)
    if not re.fullmatch(r"[a-z0-9_]+", sub.parent.name):
        sys.exit(f"model id «{sub.parent.name}»: use lower-case letters, digits and _")
    meta = json.loads(sub.read_text("utf-8"))
    if missing := [k for k in REQUIRED if k not in meta]:
        sys.exit(f"submission.json: missing fields {missing}")
    if meta.get("accepted"):
        sys.exit("submission.json: leave \"accepted\": false (set by a maintainer)")
    url = str(meta["predictions_url"])
    if not url.startswith("https://"):
        sys.exit("predictions_url must be https")
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    zpath = out / "predictions.zip"
    h, n = hashlib.md5(), 0
    with urllib.request.urlopen(url, timeout=600) as r, zpath.open("wb") as fh:
        while b := r.read(1 << 20):
            n += len(b)
            if n > MAX_BYTES:
                sys.exit("zip larger than 2 GB")
            h.update(b)
            fh.write(b)
    if h.hexdigest() != str(meta["predictions_md5"]).lower():
        sys.exit(f"MD5 mismatch: got {h.hexdigest()}")
    tasks = {json.loads(f.read_text())["name"] for f in (HERE / "tasks").glob("*.json")}
    pred = out / "predictions"
    pred.mkdir(exist_ok=True)
    with zipfile.ZipFile(zpath) as z:
        for info in z.infolist():
            name = Path(info.filename).name  # no folders, no paths outside pred/
            stem, _, ext = name.rpartition(".")
            if stem in tasks and ext in ("parquet", "csv"):
                (pred / name).write_bytes(z.read(info))
    if validate(HERE / "tasks", pred):
        sys.exit(1)
    res = evaluate(HERE / "tasks", pred)
    if not meta["external_exam_eligible"]:
        res = {k: v for k, v in res.items() if not k.startswith("external")}
        (pred / "results.json").write_text(json.dumps(res, indent=1))
    table = summary(res).round(3).to_string(index=False)
    (out / "summary.txt").write_text(f"{meta['model']} ({meta['team']})\n\n{table}\n", "utf-8")
    print(table)


if __name__ == "__main__":
    main()
