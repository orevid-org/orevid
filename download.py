"""Downloads PrimeKG (Harvard Dataverse, CC0) to data/raw/.

    py -3 download.py          # nodes, edges and drug descriptions (~405 MB)
    py -3 download.py --all    # + disease descriptions and TxGNN's kg.csv / node.csv (evaluation, benchmark)
"""
import hashlib
import sys
import urllib.request
from pathlib import Path

RAW = Path(__file__).parent / "data" / "raw"
URL = "https://dataverse.harvard.edu/api/access/datafile/{}"

# name: (Dataverse id, published md5, MB)
FILES = {
    "nodes.tab": (6180617, "7f9ab4109c54049e819ecd14e15a6038", 8.9),
    "edges.csv": (6180616, "5d4d211a22e88544b78fde2735e797bc", 386.6),
    "drug_features.tab": (6180619, "1eae9bb880675d4a26e9900c400dd4c8", 10.0),
}
OPTIONAL = {
    "disease_features.tab": (6180618, "8b38e95c5d31044deb207d6892130eb4", 113.5),
    # The PrimeKG release that TxGNN uses (grouped diseases): evaluation/txgnn.py and benchmark/bench.py
    "txgnn/kg.csv": (7144484, "aac8191d4fbc5bf09cdf8c3c78b4e75f", 982),
    "txgnn/node.csv": (7144482, "4924de04fb3deefa1e0a8dada424538e", 8.9),
}


def fetch(name, fid, md5, mb):
    dest = RAW / name
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists():
        print(f"  {name}: already there")
        return
    part = dest.with_suffix(dest.suffix + ".part")
    h = hashlib.md5()
    done, next_report = 0, 0.0
    req = urllib.request.Request(URL.format(fid), headers={"User-Agent": "orevid/0.1"})
    with urllib.request.urlopen(req) as r, open(part, "wb") as f:
        while block := r.read(1 << 20):
            f.write(block)
            h.update(block)
            done += len(block)
            if done / 1e6 >= next_report:
                print(f"\r  {name}: {done / 1e6:7.1f} / {mb} MB", end="", flush=True)
                next_report += max(mb / 20, 1)
    print(f"\r  {name}: {done / 1e6:7.1f} MB downloaded")
    if h.hexdigest() != md5:
        # Dataverse serves the .tab files in archival format; the published md5 is the original's.
        print(f"  warning: md5 {h.hexdigest()} does not match the published one ({md5})")
    part.rename(dest)


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    files = dict(FILES)
    if "--all" in sys.argv:
        files |= OPTIONAL
    for name, (fid, md5, mb) in files.items():
        fetch(name, fid, md5, mb)


if __name__ == "__main__":
    main()
