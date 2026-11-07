"""E1: blind expert rating of candidates (docs/IMPROVEMENT_PLAN.md, «Registered 07-11-2026: E1»).

    py -3 -m evaluation.e1_expert_sheet list [text]                       # eligible diseases (≥ 4 shown pairs)
    py -3 -m evaluation.e1_expert_sheet sheet E01 "disease A" "disease B"  # -> data/eval/e1/E01_sheet.xlsx + E01_key.csv
    py -3 -m evaluation.e1_expert_sheet analyse                          # data/eval/e1/returned/*.xlsx -> docs/E1_EXPERT_RATING.md

Pairs come from the frozen snapshot of 12-10-2026, with the groups of docs/PROSPECTIVE_PLAN.md. The sheet shows the
disease and the drug only, in random order; the key (pair -> group) stays in data/eval/e1/ (not public).
"""
import json
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SNAP = ROOT / "data" / "snapshot" / "2026-10-12"
OUT = ROOT / "data" / "eval" / "e1"
N_SHOWN, N_HIDDEN, N_LOWER, MIN_SHOWN = 4, 3, 3, 4
GROUPS = {"shown": "shown", "hidden": "ranked 1-25, not shown", "lower": "ranks 26-100"}
KNOWN = ["used in practice", "tried in a trial", "not that I know", "can't tell"]


def pages():
    """Archived disease pages of the snapshot: disease node -> page."""
    z = zipfile.ZipFile(SNAP / "website_data.zip")
    res = {}
    for f in z.namelist():
        if f.startswith("d/") and f.endswith(".json"):
            res.update({int(k): v for k, v in json.loads(z.read(f)).items()})
    return res


def _names(xs):
    return {(x["name"] if isinstance(x, dict) else str(x)).lower() for x in xs or []}


def excluded(page):
    """Drugs already approved, contraindicated or in a registered trial for the disease on the snapshot date."""
    out = _names(page.get("approved")) | _names(page.get("approved_ot")) | _names(page.get("approved_unconfirmed"))
    out |= _names(page.get("contra")) | _names(page.get("stopped"))
    for v in (page.get("trials") or {}).values():
        out |= _names(v.get("drugs"))
    return out


def pairs():
    """Snapshot top 100 with the E1 group of each pair; approved drugs only in the control groups."""
    from graph import Graph
    from predictor import drug_status
    status = drug_status(Graph())
    r = pd.read_csv(SNAP / "ranking_top100.csv.gz")
    r["group"] = np.where(r.website_section.isin(["extension", "candidate"]), "shown",
                          np.where(r["rank"] <= 25, "hidden", "lower"))
    r = r[(r.website_section != "exploratory")]
    approved = r.drug_node.map(lambda f: status.get(int(f)) == "approved")
    return r[(r.group == "shown") | approved]


def eligible(r, pg):
    rows = []
    for d, sub in r.groupby("disease_node"):
        ex = excluded(pg.get(int(d), {}))
        sub = sub[~sub.drug.str.lower().isin(ex)]
        n = sub.group.value_counts()
        if n.get("shown", 0) >= MIN_SHOWN and n.get("hidden", 0) >= N_HIDDEN and n.get("lower", 0) >= N_LOWER:
            rows.append((int(d), sub.disease.iloc[0], int(n["shown"]), int(n["hidden"]), int(n["lower"])))
    return pd.DataFrame(rows, columns=["disease_node", "disease", "shown", "hidden", "lower"])


def sheet(expert, diseases):
    import openpyxl
    from openpyxl.styles import Alignment, Font
    from openpyxl.worksheet.datavalidation import DataValidation
    pg = pages()
    r = pairs()
    el = eligible(r, pg)
    seed = int.from_bytes(expert.encode(), "little") % 2 ** 32
    rng = np.random.default_rng(seed)
    key = []
    for name in diseases:
        hit = el[el.disease.str.lower() == name.lower()]
        if hit.empty:
            sys.exit(f"not eligible or not found: {name!r} (see `list`)")
        d = int(hit.disease_node.iloc[0])
        sub = r[(r.disease_node == d) & ~r.drug.str.lower().isin(excluded(pg.get(d, {})))]
        for g, n in (("shown", N_SHOWN), ("hidden", N_HIDDEN), ("lower", N_LOWER)):
            pick = sub[sub.group == g].sample(n, random_state=int(rng.integers(2 ** 31)))
            key += [dict(disease_node=d, disease=x.disease, drug_node=int(x.drug_node), drug=x.drug, drug_id=x.drug_id,
                         rank=int(x.rank), section=x.website_section, group=g) for x in pick.itertuples()]
    key = pd.DataFrame(key)
    key = key.iloc[rng.permutation(len(key))].reset_index(drop=True)
    key.insert(0, "item", range(1, len(key) + 1))
    OUT.mkdir(parents=True, exist_ok=True)
    key.to_csv(OUT / f"{expert}_key.csv", index=False)

    wb = openpyxl.Workbook()
    ins = wb.active
    ins.title = "Instructions"
    text = [
        "OREVID — blind rating of drug–disease pairs (about one hour)",
        "",
        "OREVID is a non-profit, open research project (https://orevid.org). It ranks existing approved drugs as",
        "possible new uses for diseases. It is a hypothesis generator, not a clinical tool. Your help is voluntary and unpaid.",
        "",
        "Each row of the sheet «Rating» is a disease and a drug, in random order. Please rate each pair from your own",
        "knowledge, without opening the OREVID pages of these diseases until you have finished.",
        "",
        "Plausibility: could this drug benefit patients with this disease?",
        "   1 = no rationale · 2 = weak · 3 = some rationale · 4 = good rationale · 5 = strong, worth a study",
        "   Leave it empty if you cannot judge.",
        "Already known: «used in practice», «tried in a trial», «not that I know» or «can't tell».",
        "Safety concern: «yes» if the drug raises a concern specific to this disease.",
        "Comment: anything useful (mechanism, a known study, why not).",
        "",
        "Please send the file back to contact@orevid.org. Results are published only in aggregate,",
        "and your name appears only if you agree.",
    ]
    for i, t in enumerate(text, 1):
        ins.cell(i, 1, t).font = Font(bold=(i == 1))
    ins.column_dimensions["A"].width = 110
    ws = wb.create_sheet("Rating")
    head = ["Item", "Disease", "Drug", "Plausibility (1-5)", "Already known", "Safety concern", "Comment"]
    ws.append(head)
    for c in ws[1]:
        c.font = Font(bold=True)
    for x in key.itertuples():
        ws.append([x.item, x.disease, x.drug, None, None, None, None])
    n = len(key) + 1
    for col, opts in (("D", '"1,2,3,4,5"'), ("E", '"' + ",".join(KNOWN) + '"'), ("F", '"yes,no"')):
        dv = DataValidation(type="list", formula1=opts, allow_blank=True)
        ws.add_data_validation(dv)
        dv.add(f"{col}2:{col}{n}")
    for col, w in zip("ABCDEFG", (6, 40, 32, 18, 18, 15, 60)):
        ws.column_dimensions[col].width = w
    for row in ws.iter_rows(min_row=2):
        row[6].alignment = Alignment(wrap_text=True)
    f = OUT / f"{expert}_sheet.xlsx"
    wb.save(f)
    print(f"{f}: {len(key)} pairs, {len(diseases)} diseases; key {OUT / (expert + '_key.csv')}")


def boot(diff_by_disease, rng, n=10000):
    v = np.array(list(diff_by_disease.values()))
    m = rng.choice(v, (n, len(v))).mean(axis=1)
    return v.mean(), np.percentile(m, 2.5), np.percentile(m, 97.5)


def analyse():
    rows = []
    for f in sorted((OUT / "returned").glob("*_sheet.xlsx")):
        expert = f.name.split("_")[0]
        s = pd.read_excel(f, sheet_name="Rating")
        s.columns = ["item", "disease", "drug", "plaus", "known", "safety", "comment"]
        k = pd.read_csv(OUT / f"{expert}_key.csv")
        m = s.merge(k[["item", "disease_node", "group"]], on="item")
        m["expert"] = expert
        rows.append(m)
    if not rows:
        sys.exit("no returned sheets in data/eval/e1/returned/")
    a = pd.concat(rows)
    a["plaus"] = pd.to_numeric(a.plaus, errors="coerce")
    rng = np.random.default_rng(0)

    def contrast(df, g2):
        means = df.groupby(["expert", "disease_node", "group"]).plaus.mean().unstack()
        means = means.dropna(subset=["shown", g2])
        return boot(dict(zip(means.index, means["shown"] - means[g2])), rng), len(means)

    f3 = lambda x: f"{x[0]:+.2f} ({x[1]:+.2f} to {x[2]:+.2f})"
    (prim, nd), (sec, _) = contrast(a, "hidden"), contrast(a, "lower")
    unknown = a[~a.known.isin(KNOWN[:2])]
    (prim_u, nd_u) = contrast(unknown, "hidden")
    L = ["# E1. Blind expert rating\n",
         "Generated by `py -3 -m evaluation.e1_expert_sheet analyse`. Protocol: `docs/IMPROVEMENT_PLAN.md`, «Registered "
         "07-11-2026: E1». Mean plausibility (1–5) per group; differences paired by expert and disease, 95 % intervals "
         "resampling expert–disease units.\n",
         f"Experts: {a.expert.nunique()}; diseases: {a.disease_node.nunique()}; pairs rated: {a.plaus.notna().sum()} of "
         f"{len(a)}.\n",
         "| Group | Pairs | Mean plausibility | Already known | Safety concern |", "|---|---|---|---|---|"]
    for g, label in GROUPS.items():
        s = a[a.group == g]
        L.append(f"| {label} | {len(s)} | {s.plaus.mean():.2f} | {100 * s.known.isin(KNOWN[:2]).mean():.0f} % | "
                 f"{100 * (s.safety == 'yes').mean():.0f} % |")
    L += ["", f"Shown − ranked 1–25, not shown (primary): {f3(prim)} over {nd} expert–disease units.",
          f"Shown − ranks 26–100: {f3(sec)}.",
          f"Primary without pairs already known: {f3(prim_u)} over {nd_u} units.", ""]
    claim = prim[1] > 0 and a.disease_node.nunique() >= 10
    L.append(f"**Decision:** {'claim supported' if claim else 'no claim'} (rule: interval above 0 and ≥ 10 diseases; "
             f"{a.disease_node.nunique()} diseases rated).\n")
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "E1_EXPERT_RATING.md").write_text(text, "utf-8")
    print(text)


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    cmd = sys.argv[1] if len(sys.argv) > 1 else "list"
    if cmd == "sheet":
        sheet(sys.argv[2], sys.argv[3:])
    elif cmd == "analyse":
        analyse()
    else:
        el = eligible(pairs(), pages())
        if len(sys.argv) > 2:
            el = el[el.disease.str.contains(" ".join(sys.argv[2:]), case=False)]
        print(f"{len(el)} eligible diseases")
        print(el.sort_values("shown", ascending=False).head(60).to_string(index=False))
