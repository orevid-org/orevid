"""Evidence lines: disease-level confidence intervals, sensitivity of the shared-genes thresholds, a
permutation null, the split by gene-set size and the website's exact rules on the external exam
(docs/REVIEW_RESPONSE.md, tasks T3, T4, T5 and T8; decision rules written there before running).

Same candidates and sections as evaluation/evidence.py: the 25 highest-ranked candidates per test
disease; «neighbour» = indicated (training graph) for an ontology neighbour, «table» = indicated for a
disease that shares genes under the rule, «hidden» = the rest. Computed with matrices, so that every rule
variant and every permutation is cheap. The base rule (≥ 3 genes, p < 1e-6) must reproduce the published
counts of evaluation/evidence.py (checked at the start of the report).

    py -3 -m evaluation.evidence_ci          # -> data/eval/evidence_ci.json and docs/EVIDENCE_CI.md
"""
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import hypergeom

from evaluation.evidence import _transform

ROOT = Path(__file__).resolve().parent.parent
TOP, B, PERMUTATIONS = 25, 10_000, 100
SECTIONS = ("neighbour", "table", "hidden")
# Rule variants (minimum shared genes, p threshold, extra condition). «base» is production.
VARIANTS = {"base": (3, 1e-6, None)}
for k in (2, 3, 4, 5):
    for p in (1e-4, 1e-6, 1e-8):
        VARIANTS.setdefault(f"k{k}_p{p:.0e}", (k, p, None))
VARIANTS.update({"base+fold2": (3, 1e-6, "fold2"), "base+specific3": (3, 1e-6, "specific3"),
                 "base, N=genes only": (3, 1e-6, "ngenes")})
SIZE_BINS = ((0, 100, "< 100 genes"), (100, 500, "100–499 genes"), (500, 10 ** 9, "≥ 500 genes"))
_TYPES = None


def node_types_all():
    """Type of every node of the prepared graph (TxGraph.types only covers PrimeKG's nodes)."""
    global _TYPES
    if _TYPES is None:
        n = pd.read_pickle(ROOT / "data" / "prep" / "nodes.pkl").sort_values("node_index")
        _TYPES = n.node_type.to_numpy()
    return _TYPES


class Genes:
    """Shared-gene statistics between test diseases D and the diseases M where drugs are indicated."""

    def __init__(self, g, D, G=None):
        from evaluation.txgnn import GENES_EXTRA, INDICATION
        self.g, self.D = g, np.asarray(D)
        G = (g.M(GENES_EXTRA) if G is None else G).tocsr()
        I = g.M(INDICATION).tocsr()
        self.ng = np.diff(G.indptr)
        deg = np.diff(G.tocsc().indptr)
        self.N = int((deg > 0).sum())  # as evidence.py and the website (counts every linked node)
        self.N_genes = int(((deg > 0) & (node_types_all() == "gene/protein")).sum())
        is_drug = np.zeros(g.N, bool)
        is_drug[g.drugs] = True
        self.M = np.flatnonzero((self.ng > 0) & (np.diff(I.indptr) > 0) & ~is_drug)
        self.IndM = I[g.drugs][:, self.M].astype(np.float32).tocsr()  # drugs × M
        self.K = (G[self.D] @ G[self.M].T).toarray().astype(np.int32)
        specific = sp.diags(((deg < 100) & (deg > 0)).astype(np.float32))
        self.Ks = (G[self.D] @ specific @ G[self.M].T).toarray().astype(np.int32)
        posM = {m: j for j, m in enumerate(self.M)}
        for i, d in enumerate(self.D):  # the disease itself does not count
            if d in posM:
                self.K[i, posM[d]] = 0
                self.Ks[i, posM[d]] = 0
        self.ngD, self.ngM = self.ng[self.D][:, None], self.ng[self.M][None, :]

    def pvalues(self, N):
        p = np.ones(self.K.shape)
        m = self.K >= 2
        p[m] = hypergeom.sf(self.K[m] - 1, N, np.broadcast_to(self.ngM, self.K.shape)[m],
                            np.broadcast_to(self.ngD, self.K.shape)[m])
        return p

    def table(self, variant):
        """drugs × D: indicated for a disease that shares genes under the rule variant."""
        k, pv, extra = VARIANTS[variant]
        if not hasattr(self, "_p"):
            self._p = {}
        N = self.N_genes if extra == "ngenes" else self.N
        if N not in self._p:
            self._p[N] = self.pvalues(N)
        sig = (self.K >= k) & (self._p[N] < pv)
        if extra == "fold2":
            sig &= self.K * self.N >= 2 * self.ngD * self.ngM
        elif extra == "specific3":
            sig &= self.Ks >= 3
        return np.asarray(self.IndM @ sig.T.astype(np.float32)) > 0


def counts(top, V, T, hits):
    """Per disease: hits and candidates per section among `top` (list of drug-index arrays)."""
    out = {s: np.zeros((len(top), 2), int) for s in SECTIONS}
    for c, idx in enumerate(top):
        v, t, h = V[idx, c], T[idx, c], hits[c][idx]
        for s, m in (("neighbour", v), ("table", ~v & t), ("hidden", ~v & ~t)):
            out[s][c] = (int(h[m].sum()), int(m.sum()))
    return out


def bench_task(kg, types, split, seed):
    from evaluation.txgnn import CACHE, PRODUCTION, TxGraph, near_duplicates, score_paths, split_edges
    z = np.load(CACHE / "scores" / f"{split}_{seed}_{PRODUCTION}_r5_rank.npz", allow_pickle=True)
    S = _transform(z["LambdaRank_+_masked"].astype(float))
    D, pos = z["diseases"], [np.asarray(p, int) for p in z["positives"]]
    g = TxGraph(kg, split_edges(kg, split, seed) == 0, types, extra=True)
    drug_pos = np.full(g.N, -1)
    drug_pos[g.drugs] = np.arange(len(g.drugs))
    hard = near_duplicates(g, D, pos, drug_pos) < 0.5
    V = score_paths(g, D, ["ontology"])[:, :, 0] > 0
    hits = [np.isin(np.arange(len(g.drugs)), p) for p in pos]
    top = [np.argsort(-S[:, c])[:TOP] for c in range(len(D))]
    return dict(g=g, D=D, V=V, hits=hits, top=top, hard=hard)


def external_task(kg, types):
    from evaluation.external_exam import SEED, SPLIT
    from evaluation.txgnn import CACHE, TxGraph, score_paths, split_edges
    z = np.load(CACHE / "external_scores.npz", allow_pickle=True)
    p = pd.read_parquet(CACHE / "external_positives.parquet")
    D, drugs = z["diseases"], z["drugs"]
    g = TxGraph(kg, split_edges(kg, SPLIT, SEED) == 0, types, extra=True)
    assert np.array_equal(g.drugs, drugs)
    pos_f = {int(f): j for j, f in enumerate(drugs)}
    cand = np.ones((len(drugs), len(D)), bool)
    for c, known in enumerate(z["known"]):
        cand[[pos_f[int(f)] for f in known if int(f) in pos_f], c] = False
    hits, hits_new = [], []
    for d in D:
        e = p[p.disease == d]
        hits.append(np.isin(drugs, e.drug.to_numpy()))
        hits_new.append(np.isin(drugs, e[~e.obvious].drug.to_numpy()))
    raw = z["Ours_ranking"].astype(float)
    S = _transform(raw)
    top = [np.flatnonzero(cand[:, c])[np.argsort(-S[cand[:, c], c])][:TOP] for c in range(len(D))]
    V = score_paths(g, D, ["ontology"])[:, :, 0] > 0
    return dict(g=g, D=D, V=V, hits=hits, hits_new=hits_new, top=top, cand=cand, raw=raw)


# --- T5: the website's exact rules on the external exam ---------------------------------------

def website_top(kg, types, t):
    """Top 25 with the website's specificity (reference diseases) and exclusions (web/build.py)."""
    from evaluation.external_exam import SEED, SPLIT
    from evaluation.twins import drug_family, twin_neighbors, twins
    from evaluation.txgnn import CACHE, FEATURE_SETS, PRODUCTION, GENES_EXTRA, drug_degrees, score_paths
    from evaluation.xgboost_compare import CONFIGS, OUT, _fit
    from graph import CONTRA, INDICATION, OFF_LABEL, Graph
    from predictor import diagnostic_drugs, drug_status
    g = t["g"]
    # 1. specificity over 2,000 reference diseases of the seed-1 training graph, production model
    data = np.load(OUT / f"{SPLIT}_{SEED}.npz", allow_pickle=True)  # training rows (xgboost_compare)
    model = _fit("lgb_rank_n300_l31", data["X"], data["Y"], data["groups"], SEED)
    assert np.allclose(model(np.column_stack([data["F_ext"][:, 0, :].astype(np.float64), data["degree"], data["gf"]])),
                       t["raw"][:, 0], atol=1e-6), "the retrained model does not reproduce the published scores"
    paths = FEATURE_SETS[PRODUCTION]
    ng = np.diff(g.M(GENES_EXTRA).indptr)
    nph = np.diff(g.M("disease_phenotype_positive").indptr)
    from evaluation.external_exam import N_PKG
    dis = np.flatnonzero((node_types_all() == "disease") & ((ng > 0) | (nph > 0)) & (np.arange(g.N) < N_PKG))  # PrimeKG, as predictor.reference
    sample = np.random.default_rng(0).choice(dis, min(2000, len(dis)), replace=False)
    degree, gf = data["degree"], data["gf"]
    ref = np.zeros(len(g.drugs))
    for i in range(0, len(sample), 128):
        F = score_paths(g, sample[i:i + 128], paths)
        for c in range(F.shape[1]):
            s = model(np.column_stack([F[:, c, :], degree, gf]))
            ref += (s - s.mean()) / (s.std() + 1e-9)
    ref /= len(sample)
    z = (t["raw"] - t["raw"].mean(axis=0)) / (t["raw"].std(axis=0) + 1e-9)
    S = z - 0.5 * ref[:, None]
    # 2. exclusions (without Open Targets approvals: here they are the answer)
    gc = Graph()
    tw, diag, status = twins(gc, "disease"), diagnostic_drugs(gc), drug_status(gc)
    drugs = g.drugs
    fam = [drug_family(gc.name(int(f))) for f in drugs]
    cand = t["cand"].copy()
    for c, d in enumerate(t["D"]):
        dt = tw.get(int(d), {int(d)})
        known = twin_neighbors(gc, dt, INDICATION) | twin_neighbors(gc, dt, OFF_LABEL)
        contra = twin_neighbors(gc, dt, CONTRA)
        names_k = {x for i in known for x in drug_family(gc.name(i))}
        names_c = {x for i in contra for x in drug_family(gc.name(i))}
        for j, f in enumerate(drugs):
            f = int(f)
            if (f in known or f in contra or f in diag or status.get(f) == "veterinary only"
                    or fam[j] & names_k or fam[j] & names_c):
                cand[j, c] = False
    top = [np.flatnonzero(cand[:, c])[np.argsort(-S[cand[:, c], c])][:TOP] for c in range(len(t["D"]))]
    return top, cand


# --- permutation null (T4) ----------------------------------------------------------------------

def permuted_genes(g, G, rng, oversample=1.6):
    """Every disease keeps its number of genes; genes are redrawn with probability proportional to their
    number of diseases (draw with replacement, keep the first distinct ones)."""
    G = G.tocsr()
    is_gene = node_types_all() == "gene/protein"
    rows = np.flatnonzero((np.diff(G.indptr) > 0) & ~is_gene)
    deg = np.diff(G.tocsc().indptr).astype(float)
    genes = np.flatnonzero(is_gene & (deg > 0))
    prob = deg[genes] / deg[genes].sum()
    size = np.diff(G.indptr)[rows]
    want = (size * oversample + 10).astype(int)
    draw = genes[rng.choice(len(genes), int(want.sum()), p=prob)]
    owner = np.repeat(rows, want)
    df = pd.DataFrame({"r": owner, "c": draw}).drop_duplicates()
    df["k"] = df.groupby("r").cumcount()
    df = df[df.k.to_numpy() < np.repeat(size, df.groupby("r").size().reindex(rows).fillna(0).astype(int).to_numpy())]
    A = sp.csr_matrix((np.ones(len(df), np.float32), (df.r, df.c)), shape=G.shape)
    return ((A + A.T) > 0).astype(np.float32).tocsr()  # symmetric, as the graph's matrices


# --- bootstrap and report -------------------------------------------------------------------------

def bootstrap(parts, rng):
    """parts: list of (h, n) per-disease arrays, one per task; resamples diseases within tasks.
    Returns (rate, low, high) of sum(h) / sum(n)."""
    H = np.zeros(B)
    Nn = np.zeros(B)
    h_tot = n_tot = 0
    for h, n in parts:
        if not len(h):
            continue
        idx = rng.integers(0, len(h), (B, len(h)))
        H += h[idx].sum(axis=1)
        Nn += n[idx].sum(axis=1)
        h_tot += h.sum()
        n_tot += n.sum()
    if not n_tot:
        return None
    r = H / np.maximum(Nn, 1)
    return float(h_tot / n_tot), float(np.percentile(r, 2.5)), float(np.percentile(r, 97.5)), int(n_tot)


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    from evaluation.txgnn import AREAS, GENES_EXTRA, load_kg, node_types
    t0 = time.time()
    kg = load_kg()
    types = node_types(kg)
    rng = np.random.default_rng(0)
    groups = {"benchmark": [("complex_disease", s) for s in (1, 2, 3, 4, 5)], "areas": [(a, 1) for a in AREAS]}
    store = {}  # (group, subset, variant, section, size_bin) -> list of (h, n)

    def add(group, subset, variant, cnt, sizes):
        for s in SECTIONS:
            h, n = cnt[s][:, 0], cnt[s][:, 1]
            store.setdefault((group, subset, variant, s, "all sizes"), []).append((h, n))
            for lo, hi, name in SIZE_BINS:
                m = (sizes >= lo) & (sizes < hi)
                store.setdefault((group, subset, variant, s, name), []).append((h[m], n[m]))

    for group, tasks in groups.items():
        for split, seed in tasks:
            t = bench_task(kg, types, split, seed)
            gn = Genes(t["g"], t["D"])
            sizes = gn.ng[t["D"]]
            for v in VARIANTS:
                cnt = counts(t["top"], t["V"], gn.table(v), t["hits"])
                add(group, "all", v, cnt, sizes)
                hard = t["hard"]
                add(group, "hard", v, {s: c[hard] for s, c in cnt.items()}, sizes[hard])
            print(f"  {group} {split} {seed} ({time.time() - t0:.0f} s)", flush=True)

    t = external_task(kg, types)
    gn = Genes(t["g"], t["D"])
    sizes = gn.ng[t["D"]]
    tables = {v: gn.table(v) for v in VARIANTS}
    for v in VARIANTS:
        add("external", "all", v, counts(t["top"], t["V"], tables[v], t["hits"]), sizes)
        add("external", "non-obvious approvals", v, counts(t["top"], t["V"], tables[v], t["hits_new"]), sizes)
    # model only (top 25, no rule) and rules only (every candidate with a line, no ranking): check C
    extra = {}
    hits = t["hits"]
    extra["model only: top 25"] = [(np.array([hits[c][i].sum() for c, i in enumerate(t["top"])]),
                                    np.array([len(i) for i in t["top"]]))]
    shown = (t["V"] | tables["base"]) & t["cand"]
    extra["rules only: every candidate with a line"] = [(np.array([hits[c][shown[:, c]].sum() for c in range(len(t["D"]))]),
                                                         shown.sum(axis=0))]
    extra["no rule, no model: every candidate"] = [(np.array([hits[c][t["cand"][:, c]].sum() for c in range(len(t["D"]))]),
                                                    t["cand"].sum(axis=0))]
    print(f"  external ({time.time() - t0:.0f} s)", flush=True)

    # T5: website's exact rules
    web_top, web_cand = website_top(kg, types, t)
    add("external, website rules", "all", "base", counts(web_top, t["V"], tables["base"], t["hits"]), sizes)
    add("external, website rules", "non-obvious approvals", "base",
        counts(web_top, t["V"], tables["base"], t["hits_new"]), sizes)
    extra["website rules: model only, top 25"] = [(np.array([hits[c][i].sum() for c, i in enumerate(web_top)]),
                                                   np.array([len(i) for i in web_top]))]
    print(f"  website rules ({time.time() - t0:.0f} s)", flush=True)

    # T4: permutation null of the shared-genes section on the external exam
    G = t["g"].M(GENES_EXTRA)
    null = []
    prng = np.random.default_rng(1)
    for k in range(PERMUTATIONS):
        gp = Genes(t["g"], t["D"], G=permuted_genes(t["g"], G, prng))
        cnt = counts(t["top"], t["V"], gp.table("base"), t["hits"])
        h, n = cnt["table"][:, 0].sum(), cnt["table"][:, 1].sum()
        null.append((int(h), int(n)))
        if (k + 1) % 10 == 0:
            print(f"  permutation {k + 1}/{PERMUTATIONS} ({time.time() - t0:.0f} s)", flush=True)

    res = {"|".join(k): bootstrap(v, rng) for k, v in store.items()}
    res_extra = {k: bootstrap(v, rng) for k, v in extra.items()}
    out = {"sections": res, "checks": res_extra, "null_table_external": null}
    (ROOT / "data" / "eval" / "evidence_ci.json").write_text(json.dumps(out, indent=1), "utf-8")
    report(out)


def fmt(x):
    if not x:
        return "—"
    r, lo, hi, n = x
    return f"{100 * r:.1f} ({100 * lo:.1f}–{100 * hi:.1f}; {n:,})"


def report(out=None):
    out = out or json.loads((ROOT / "data" / "eval" / "evidence_ci.json").read_text("utf-8"))
    R = out["sections"]
    get = lambda *k: R.get("|".join(k))
    L = ["# Evidence lines: disease-level intervals, threshold sensitivity, permutation null and website rules\n",
         "Generated by `evaluation/evidence_ci.py` (docs/REVIEW_RESPONSE.md, T3, T4, T5, T8; decision rules fixed "
         "there before running). Hit rate in % among the 25 highest-ranked candidates per disease, with a 95% "
         "interval from a bootstrap that resamples whole diseases (10,000 resamples, within seeds or areas), and "
         "the number of candidates. Benchmark: zero-shot seeds 1–5; areas: nine areas, seed 1; external: Open "
         "Targets approvals absent from PrimeKG (training split zero-shot seed 1).\n",
         "## T3. Sections with disease-level intervals (production rule)\n",
         "| Test | Direct extensions | Shared genes | Hidden |", "|---|---|---|---|"]
    rows = [("Benchmark", "benchmark", "all"), ("Benchmark, hard", "benchmark", "hard"), ("Areas", "areas", "all"),
            ("Areas, hard", "areas", "hard"), ("External exam", "external", "all"),
            ("External, non-obvious approvals", "external", "non-obvious approvals"),
            ("External, website rules (T5)", "external, website rules", "all"),
            ("External, website rules, non-obvious", "external, website rules", "non-obvious approvals")]
    for name, gname, sub in rows:
        L.append(f"| {name} | " + " | ".join(fmt(get(gname, sub, "base", s, "all sizes")) for s in SECTIONS) + " |")
    L += ["", "Checks on the external exam (check C: model and rules measured apart):\n",
          "| Candidates | Hit rate |", "|---|---|"]
    L += [f"| {k} | {fmt(v)} |" for k, v in out["checks"].items()]
    L += ["", "## T4. Sensitivity of the shared-genes rule (hit rate of the «shared genes» section)\n",
          "| Variant | Benchmark | Areas | External |", "|---|---|---|---|"]
    for v in VARIANTS:
        L.append(f"| {v} | " + " | ".join(fmt(get(gname, "all", v, "table", "all sizes"))
                                          for gname in ("benchmark", "areas", "external")) + " |")
    null = out["null_table_external"]
    rates = np.array([h / n if n else 0.0 for h, n in null])
    real = get("external", "all", "base", "table", "all sizes")
    L += ["", f"**Permutation null** ({len(null)} permutations; every disease keeps its number of genes, genes redrawn "
          f"with probability proportional to their number of diseases; ranking and indications unchanged): "
          f"shared-genes hit rate on the external exam {100 * rates.mean():.2f} % on average "
          f"(99th percentile {100 * np.percentile(rates, 99):.2f} %; candidates per permutation "
          f"{int(np.mean([n for _, n in null])):,}), against {100 * real[0]:.2f} % observed "
          f"({'above' if real[0] > np.percentile(rates, 99) else 'NOT above'} the 99th percentile).\n",
          "## T8. Shared genes by the disease's gene-set size\n",
          "| Test | Size | Shared genes (production) | + fold enrichment ≥ 2 | + ≥ 3 specific genes | Hidden |",
          "|---|---|---|---|---|---|"]
    for gname in ("benchmark", "areas", "external"):
        for _, _, b in SIZE_BINS:
            L.append(f"| {gname} | {b} | " + " | ".join(fmt(get(gname, "all", v, s, b)) for v, s in
                                                       (("base", "table"), ("base+fold2", "table"),
                                                        ("base+specific3", "table"), ("base", "hidden"))) + " |")
    text = "\n".join(L) + "\n"
    (ROOT / "docs" / "EVIDENCE_CI.md").write_text(text, "utf-8")
    print(text)


def null_extra():
    """Exploratory, not pre-specified (added after the T4 null failed): the permutation null stored per
    disease, so that it can be compared within gene-set-size bins (the null favours large gene sets), for
    the production rule and the «≥ 3 specific genes» variant on the external exam (100 permutations), and
    for the production rule on the benchmark (zero-shot seeds 1-5, 20 permutations each).
    -> data/eval/evidence_null_extra.npz"""
    from evaluation.txgnn import GENES_EXTRA, load_kg, node_types
    t0 = time.time()
    kg = load_kg()
    types = node_types(kg)
    out = {}
    t = external_task(kg, types)
    G = t["g"].M(GENES_EXTRA)
    gn = Genes(t["g"], t["D"])
    out["ext_sizes"] = gn.ng[t["D"]]
    for v in ("base", "base+specific3"):
        out[f"ext_real_{v}"] = counts(t["top"], t["V"], gn.table(v), t["hits"])["table"]
    prng = np.random.default_rng(2)
    nb, ns = [], []
    for k in range(PERMUTATIONS):
        gp = Genes(t["g"], t["D"], G=permuted_genes(t["g"], G, prng))
        nb.append(counts(t["top"], t["V"], gp.table("base"), t["hits"])["table"])
        ns.append(counts(t["top"], t["V"], gp.table("base+specific3"), t["hits"])["table"])
        if (k + 1) % 20 == 0:
            print(f"  external permutation {k + 1} ({time.time() - t0:.0f} s)", flush=True)
    out["ext_null_base"], out["ext_null_base+specific3"] = np.array(nb), np.array(ns)
    rb, nbm, sz = [], [], []
    for seed in (1, 2, 3, 4, 5):
        t = bench_task(kg, types, "complex_disease", seed)
        G = t["g"].M(GENES_EXTRA)
        gn = Genes(t["g"], t["D"])
        sz.append(gn.ng[t["D"]])
        rb.append(counts(t["top"], t["V"], gn.table("base"), t["hits"])["table"])
        prng = np.random.default_rng(10 + seed)
        nbm.append(np.array([counts(t["top"], t["V"], Genes(t["g"], t["D"], G=permuted_genes(t["g"], G, prng)).table("base"),
                                    t["hits"])["table"] for _ in range(20)]))
        print(f"  benchmark seed {seed} ({time.time() - t0:.0f} s)", flush=True)
    out["bench_sizes"], out["bench_real_base"] = np.concatenate(sz), np.concatenate(rb)
    out["bench_null_base"] = np.concatenate(nbm, axis=1)
    np.savez(ROOT / "data" / "eval" / "evidence_null_extra.npz", **out)
    null_report()


def null_report():
    z = np.load(ROOT / "data" / "eval" / "evidence_null_extra.npz")
    L = ["| Test | Rule | Gene-set size | Observed | Null mean | Null 99th percentile | Observed above? |",
         "|---|---|---|---|---|---|---|"]
    bins = (("all sizes", 0, 10 ** 9),) + tuple((b, lo, hi) for lo, hi, b in SIZE_BINS)
    for test, rule, real, null, sizes in (
            ("external", "production", z["ext_real_base"], z["ext_null_base"], z["ext_sizes"]),
            ("external", "≥ 3 specific genes", z["ext_real_base+specific3"], z["ext_null_base+specific3"], z["ext_sizes"]),
            ("benchmark", "production", z["bench_real_base"], z["bench_null_base"], z["bench_sizes"])):
        for name, lo, hi in bins:
            m = (sizes >= lo) & (sizes < hi)
            h, n = real[m].sum(axis=0)
            if not n:
                continue
            rates = np.array([x[m, 0].sum() / max(x[m, 1].sum(), 1) for x in null])
            cands = int(np.mean([x[m, 1].sum() for x in null]))
            obs = h / n
            p99 = np.percentile(rates, 99)
            L.append(f"| {test} | {rule} | {name} | {100 * obs:.1f} % ({n:,}) | {100 * rates.mean():.1f} % ({cands:,}) | "
                     f"{100 * p99:.1f} % | {'yes' if obs > p99 else 'no'} |")
    text = "\n".join(L)
    print(text)
    return text


if __name__ == "__main__":
    if sys.argv[1:2] == ["null"]:
        sys.stdout.reconfigure(encoding="utf-8")
        null_report() if sys.argv[2:] == ["report"] else null_extra()
    elif sys.argv[1:] == ["report"]:
        sys.stdout.reconfigure(encoding="utf-8")
        report()
    else:
        main()
