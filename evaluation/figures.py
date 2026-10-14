"""Figures of the preprint, drawn from the result files in the repository (no hand-entered results).

    Figure 1  method overview (diagram; the counts are those of the graph described in the README)
    Figure 2  Recall@100 of the popularity baseline, TxGNN and OREVID on four tests
              (benchmark/predictions/<model>/results.json)
    Figure 3  Recall@100 per zero-shot seed and per disease area (same files)
    Figure 4  hit rate among the 25 highest-ranked candidates, by website section, with disease-level
              intervals (table «T3» of docs/EVIDENCE_CI.md, written by evaluation/evidence_ci.py)

    py -3 -m evaluation.figures          # -> docs/preprint/figures/figure{1..4}.png and .pdf

Only needs matplotlib (no data/ directory), so it runs on any computer with a clone of the repository.
"""
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

ROOT = Path(__file__).resolve().parent.parent
PRED = ROOT / "benchmark" / "predictions"
OUT = ROOT / "docs" / "preprint" / "figures"

BLUE, TXGNN, POP, GRID, TEXT, MUTED = "#2878d6", "#c3c2b6", "#e4e3df", "#e6e6e6", "#1f1f1f", "#5f5f5f"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.edgecolor": "#bdbdbd",
                     "axes.labelcolor": MUTED, "xtick.color": MUTED, "ytick.color": TEXT,
                     "savefig.dpi": 300, "savefig.bbox": "tight"})

AREAS = {"adrenal_gland": "Adrenal gland", "anemia": "Anaemia", "autoimmune": "Autoimmune",
         "cardiovascular": "Cardiovascular", "cell_proliferation": "Cell proliferation", "diabetes": "Diabetes",
         "mental_health": "Mental health", "metabolic_disorder": "Metabolic disorder",
         "neurodigenerative": "Neurodegenerative"}
MODELS = [("popularity", "Popularity baseline", POP), ("txgnn", "TxGNN", TXGNN), ("orevid_lambdarank", "OREVID", BLUE)]


def results(model):
    return json.loads((PRED / model / "results.json").read_text(encoding="utf-8"))


def recall(res, task, subset="all"):
    return res[task][subset]["Recall@100"]


def save(fig, n):
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"figure{n}.{ext}")
    plt.close(fig)
    print(f"  figure{n}.png / .pdf")


def figure1():
    """Method overview: from the sources to what a disease page shows."""
    fig, ax = plt.subplots(figsize=(10, 5.6))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 56)
    ax.axis("off")
    w, h = 21, 9

    def box(x, y, title, sub, hl=False, height=h):
        ax.add_patch(FancyBboxPatch((x, y), w, height, boxstyle="round,pad=0,rounding_size=1.2", lw=1.3,
                                    ec=BLUE if hl else "#c9c9c2", fc="#e6f0fb" if hl else "white"))
        ax.text(x + w / 2, y + height * 0.66, title, ha="center", va="center", fontsize=10.5, fontweight="bold",
                color=TEXT)
        ax.text(x + w / 2, y + height * 0.3, sub, ha="center", va="center", fontsize=8.5, color=MUTED,
                linespacing=1.3)

    def arrow(p, q):
        ax.add_patch(FancyArrowPatch(p, q, arrowstyle="-|>", mutation_scale=9, lw=1.1, color="#b5b5ad",
                                     shrinkA=0, shrinkB=0))

    top, mid = 41, 20
    xs = [2, 26.5, 51, 75.5]
    box(xs[0], top, "Sources", "PrimeKG, Open Targets\n26.06, STRING 12.5")
    box(xs[1], top, "Knowledge graph", "157,822 nodes\n9,457,076 edges")
    box(xs[2], top, "Features", "22 metapaths and\ndrug degrees")
    box(xs[3], top, "Ranking model", "LambdaRank, trained\non masked diseases")
    for a, b in zip(xs, xs[1:]):
        arrow((a + w, top + h / 2), (b, top + h / 2))
    # ranking model -> specificity (elbow)
    ax.plot([xs[3] + w / 2, xs[3] + w / 2, xs[0] + w / 2], [top, top - 3, top - 3], color="#b5b5ad", lw=1.1)
    arrow((xs[0] + w / 2, top - 3), (xs[0] + w / 2, mid + h))
    box(xs[0], mid, "Specificity", "minus 0.5 × the\ndrug's mean score")
    box(xs[1], mid, "Top 25 per disease", "not yet approved,\nafter exclusions")
    box(xs[2], mid, "Evidence line?", "neighbour indication\nor shared genes", hl=True)
    for a, b in zip(xs[:2], xs[1:3]):
        arrow((a + w, mid + h / 2), (b, mid + h / 2))
    outs = [(mid + 9.5, "Direct extension", "shown"), (mid, "Shared genes", "shown, with p-value"),
            (mid - 9.5, "Hidden", "no line: not shown")]
    bx = xs[3]
    trunk = xs[2] + w + 1.8
    ax.plot([xs[2] + w, trunk], [mid + h / 2] * 2, color="#b5b5ad", lw=1.1)
    ax.plot([trunk, trunk], [outs[-1][0] + 3.5, outs[0][0] + 3.5], color="#b5b5ad", lw=1.1)
    for y, t, s in outs:
        box(bx, y, t, s, height=7)
        arrow((trunk, y + 3.5), (bx, y + 3.5))
    ax.text(2, 4, "Also on each page: a risk warning (contraindication rule), an exploratory list of drugs with no prior "
                  "indication\n(not validated) and literature relations (PubTator3) as context. None of them changes "
                  "the ranking or the evidence lines.", fontsize=8.5, color=MUTED, va="center", linespacing=1.4)
    save(fig, 1)


def figure2():
    """Recall@100 of three models on four tests (means over seeds or areas, as in Table 1)."""
    tests = [("Zero-shot", "zero_shot", "all"), ("Zero-shot, hard", "zero_shot", "hard"),
             ("Disease areas held out", "area_", "all"), ("External exam", "external", "all")]
    vals = {}
    for key, _, _ in MODELS:
        res = results(key)
        row = []
        for _, prefix, sub in tests:
            tasks = [t for t in res if t.startswith(prefix) and sub in res[t]]
            row.append(sum(recall(res, t, sub) for t in tasks) / len(tasks))
        vals[key] = row
    fig, ax = plt.subplots(figsize=(9, 3.9))
    width = 0.26
    for i, (key, label, colour) in enumerate(MODELS):
        xs = [j + (i - 1) * (width + 0.02) for j in range(len(tests))]
        bars = ax.bar(xs, vals[key], width, color=colour, label=label, zorder=2)
        for b, v in zip(bars, vals[key]):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.015, f"{v:.2f}", ha="center", va="bottom", fontsize=9,
                    color=TEXT if key == "orevid_lambdarank" else MUTED)
    ax.set_xticks(range(len(tests)), [t[0] for t in tests])
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.tick_params(axis="x", length=0, labelsize=10, labelcolor=TEXT)
    ax.set_title("Recall@100 of three models on four tests", loc="left", fontsize=13, color=TEXT, pad=34)
    ax.text(0, 1.1, "Recall@100: share of a disease's correct drugs ranked in the top 100 of 7,957 (higher is better)",
            transform=ax.transAxes, fontsize=9, color=MUTED)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 0.93), ncol=3, frameon=False, fontsize=9.5, handlelength=1,
              handleheight=1)
    save(fig, 2)


def figure3():
    """Recall@100 on each zero-shot seed and each disease area (dumbbell plot)."""
    res = {key: results(key) for key, _, _ in MODELS}
    rows = [(f"Seed {s}", f"zero_shot_{s}") for s in range(1, 6)]
    rows += [(f"{name} ({res['orevid_lambdarank'][f'area_{a}']['n_test_diseases']})", f"area_{a}")
             for a, name in AREAS.items()]
    ys, labels, headers, y = [], [], [], 0
    for i, (label, _) in enumerate(rows):
        if i in (0, 5):  # a header line before each group, and a gap before the second one
            y += 0.6 if i else 0
            headers.append(-y)
            y += 1
        ys.append(-y)
        labels.append(label)
        y += 1
    fig, ax = plt.subplots(figsize=(8.5, 7.4))
    ahead = 0
    for yy, (_, task) in zip(ys, rows):
        p, t, o = (recall(res[k], task) for k in ("popularity", "txgnn", "orevid_lambdarank"))
        ahead += o > t
        ax.plot([t, o], [yy, yy], color=GRID, lw=1.6, zorder=1)
        ax.scatter([p], [yy], s=55, facecolors="white", edgecolors=TXGNN, linewidths=1.5, zorder=2)
        ax.scatter([t], [yy], s=55, color=TXGNN, zorder=3)
        ax.scatter([o], [yy], s=70, color=BLUE, zorder=4)
        ax.text(o + 0.015, yy, f"{o:.2f}", va="center", fontsize=9, color=MUTED)
    ax.set_yticks(ys, labels)
    ax.set_xlim(0, 1)
    ax.set_xlabel("Recall@100 (share of correct drugs in the top 100)", loc="right")
    ax.grid(axis="x", color=GRID, zorder=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    for yy, header in zip(headers, ("Zero-shot, by seed", "Disease areas (test diseases)")):
        ax.text(-0.42, yy, header, transform=ax.get_yaxis_transform(), va="center", fontweight="bold", color=TEXT)
    ax.set_ylim(min(ys) - 0.7, headers[0] + 0.7)
    ax.set_title("Recall@100 per zero-shot seed and per disease area", loc="left", fontsize=13, color=TEXT,
                 pad=58, x=-0.42)
    ax.text(-0.42, 1.1, f"OREVID above TxGNN in {ahead} of {len(rows)} splits", transform=ax.transAxes,
            fontsize=9.5, color=MUTED)
    handles = [plt.Line2D([], [], marker="o", ls="", color=BLUE, ms=8, label="OREVID"),
               plt.Line2D([], [], marker="o", ls="", color=TXGNN, ms=8, label="TxGNN"),
               plt.Line2D([], [], marker="o", ls="", mfc="white", mec=TXGNN, mew=1.5, ms=8,
                          label="Popularity baseline")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.42, 1.01), ncol=3, frameon=False, fontsize=9.5)
    save(fig, 3)


def evidence_table():
    """{test: {section: (hit %, low, high, n)}} from the T3 table of docs/EVIDENCE_CI.md."""
    text = (ROOT / "docs" / "EVIDENCE_CI.md").read_text(encoding="utf-8")
    block = text.split("## T3.", 1)[1].split("\n\n", 2)[1]
    cell = re.compile(r"([\d.]+) \(([\d.]+)–([\d.]+); ([\d,]+)\)")
    rows = {}
    for line in block.splitlines()[2:]:
        parts = [p.strip() for p in line.strip("|").split("|")]
        vals = [cell.match(p).groups() for p in parts[1:4]]
        rows[parts[0]] = {s: (float(a), float(b), float(c), int(d.replace(",", "")))
                          for s, (a, b, c, d) in zip(("Direct extensions", "Shared genes", "Hidden"), vals)}
    return rows


def figure4():
    """Hit rate among the top 25 candidates by website section, with disease-level 95% intervals."""
    table = evidence_table()
    tests = ["Benchmark", "Benchmark, hard", "Areas", "Areas, hard", "External exam"]
    styles = {"Direct extensions": dict(color=TXGNN, mfc=TXGNN), "Shared genes": dict(color=BLUE, mfc=BLUE),
              "Hidden": dict(color=TXGNN, mfc="white")}
    offsets = {"Direct extensions": 0.25, "Shared genes": 0.0, "Hidden": -0.25}
    fig, ax = plt.subplots(figsize=(9, 5.4))
    for i, test in enumerate(tests):
        y0 = -i
        for section, off in offsets.items():
            v, lo, hi, _ = table[test][section]
            st = styles[section]
            ax.plot([lo, hi], [y0 + off] * 2, color=st["color"], lw=1.6, zorder=1)
            ax.plot([v], [y0 + off], marker="o", ms=8, color=st["color"], mfc=st["mfc"], mew=1.5, zorder=2)
            ax.text(hi + 1, y0 + off, f"{v:.1f}%", va="center", fontsize=9, color=MUTED)
    ax.set_yticks([-i for i in range(len(tests))], tests)
    ax.set_xlim(0, 72)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, _: f"{x:.0f}%"))
    ax.grid(axis="x", color=GRID, zorder=0)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=0)
    ax.set_title("Hit rate among the top 25 candidates, by website section", loc="left", fontsize=13, color=TEXT,
                 pad=48, x=-0.2)
    ax.text(-0.2, 1.1, "Percentage correct, 95% CI resampling whole diseases. Chance in the external exam is "
                       "about 0.07%.", transform=ax.transAxes, fontsize=9.5, color=MUTED)
    handles = [plt.Line2D([], [], marker="o", ls="", color=TXGNN, ms=8, label="Direct extensions"),
               plt.Line2D([], [], marker="o", ls="", color=BLUE, ms=8, label="Shared genes"),
               plt.Line2D([], [], marker="o", ls="", mfc="white", mec=TXGNN, mew=1.5, ms=8, label="Hidden")]
    ax.legend(handles=handles, loc="lower left", bbox_to_anchor=(-0.2, 1.0), ncol=3, frameon=False, fontsize=9.5)
    save(fig, 4)


def main():
    print(f"-> {OUT}")
    figure1()
    figure2()
    figure3()
    figure4()


if __name__ == "__main__":
    main()
