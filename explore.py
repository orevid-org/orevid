"""Graph explorer.

    py -3 explore.py summary
    py -3 explore.py search "fibrosis"            [--type disease]
    py -3 explore.py drug "metformin"
    py -3 explore.py disease "cystic fibrosis"  [--top 20] [--paths "target all,similar all"]
    py -3 explore.py path "metformin" "alzheimer disease"
    py -3 explore.py targets "multiple sclerosis"     [--top 15]
    py -3 explore.py validate                       [--n 800] [--strict] [--paths all]
    py -3 explore.py temporal                     [--cutoff 2015] [--n 800]

An ambiguous name can be given by index: "#12345".
"""
import argparse
import json
import sys

import numpy as np
from scipy.stats import rankdata

from predictor import MODEL_FILE as PREDICTOR, Predictor
from graph import (CLINICAL, COMBINED, CONTRA, TARGET, EFFECT, FAERS, GENES_ALL, INDICATION, METAPATHS, MODEL, NEGATIVE, WEIGHT_CAUSE,
                   WEIGHT_FAILED, OFF_LABEL, Graph)


def one_node(g, text, node_type=None):
    found = g.resolve(text, node_type)
    if len(found) == 1:
        return found[0]
    if not found:
        raise SystemExit(f"Could not find «{text}»" + (f" of type {node_type}" if node_type else "") + ".")
    print(f"«{text}» is ambiguous ({len(found)} matches). Pick one with #index:")
    for i in found[:25]:
        print(f"  #{i:<7} {g.types[i]:20s} {g.name(i)}  [{g.origin[i]}]")
    raise SystemExit(1)


def name_list(g, idx, maximum=12):
    names = sorted(g.name(i) for i in idx)
    rest = f" … (+{len(names) - maximum})" if len(names) > maximum else ""
    return ", ".join(names[:maximum]) + rest if names else "—"


def read_paths(text, default):
    if not text:
        return list(default)
    if text == "all":
        return list(METAPATHS)
    paths = [c.strip() for c in text.split(",")]
    for c in paths:
        if c not in METAPATHS:
            raise SystemExit(f"Unknown metapath: {c!r}. Available: {', '.join(METAPATHS)}")
    return paths


def auroc(scores, positives):
    r = rankdata(scores)
    npos, nneg = positives.sum(), (~positives).sum()
    return (r[positives].sum() - npos * (npos + 1) / 2) / (npos * nneg)


# --- commands ---------------------------------------------------------------

def summary(g, a):
    print(f"Nodes: {g.N:,}\n")
    print(g.nodes.groupby(["node_type", "origin"]).size().to_string())
    print("\nRelations:")
    for name, r in sorted(g.relations.items(), key=lambda kv: -kv[1]["edges"]):
        print(f"  {name:45s} {r['edges']:>10,}")


def search(g, a):
    found = g.resolve(a.text, a.type)
    print(f"{len(found)} matches")
    for i in found[:a.top]:
        print(f"  #{i:<7} {g.types[i]:20s} {g.name(i)}  [{g.origin[i]}]")


def drug(g, a):
    f = one_node(g, a.name, "drug")
    print(f"== {g.name(f)}  (#{f}, {g.nodes.node_id[f]}, {g.origin[f]}) ==")
    if g.descriptions is not None and f in g.descriptions.index:
        d = g.descriptions.loc[f]
        for field in ("mechanism_of_action", "description"):
            if field in d and isinstance(d[field], str):
                print(f"\n{field.replace('_', ' ').capitalize()}:\n  {d[field][:600]}")
    atr = g.attributes()
    moa = atr[(atr.x_index == f) & (atr.relation == "moa")]
    print(f"\nTargets:           {name_list(g, g.neighbors(f, TARGET))}")
    for action, s in moa.groupby("display_relation"):
        print(f"  {action:16s} {name_list(g, s.y_index.unique())}")
    print(f"Indicated for:     {name_list(g, g.neighbors(f, INDICATION))}")
    print(f"Off-label use:     {name_list(g, g.neighbors(f, OFF_LABEL))}")
    print(f"Contraindicated in: {name_list(g, g.neighbors(f, CONTRA))}")
    clin = atr[(atr.x_index == f) & (atr.relation == CLINICAL)]
    if len(clin):
        print(f"\nClinical development ({len(clin)} diseases):")
        for phase in ["approved", "phase 3", "phase 2", "phase 1", "preclinical", "unknown"]:
            s = clin[clin.display_relation == phase]
            if len(s):
                print(f"  {phase:11s} {name_list(g, s.y_index, 10)}")
    neg = atr[(atr.x_index == f) & (atr.relation == NEGATIVE)]
    for reason, s in neg.groupby("display_relation"):
        print(f"  ⚠ trials stopped for {reason}: {name_list(g, s.y_index, 10)}")
    effects = g.neighbors(f, EFFECT)
    print(f"\nAdverse effects (SIDER, {len(effects)}): {name_list(g, effects, 15)}")
    fa = atr[(atr.x_index == f) & (atr.relation == FAERS)].sort_values("score", ascending=False)
    if len(fa):
        print(f"Pharmacovigilance signals (FAERS, {len(fa)}): "
              + ", ".join(g.name(i) for i in fa.y_index[:15]))
    print(f"Interactions with other drugs: {len(g.neighbors(f, 'drug_drug')):,}")


def disease(g, a):
    d = one_node(g, a.name, "disease")
    if a.paths or not PREDICTOR.exists():
        return disease_metapaths(g, a, d)
    # Production: the same page as the web (web/build.py, DiseasePages), so that they do not drift apart.
    from web.build import DiseasePages, explanation_text
    F = DiseasePages(top=a.top)
    m = F.disease(d)
    print(f"== {m['name']}  (#{d}, MONDO {m['mondo']}) · {m['n_genes']} associated genes ==")
    print(f"Already approved or indicated ({len(m['approved'])}): {', '.join(m['approved'][:15]) or '—'}")
    if m["offlabel"]:
        print(f"Off-label: {', '.join(m['offlabel'][:10])}")
    if m["contra"]:
        print(f"Contraindicated: {', '.join(m['contra'][:10])}")
    for phase, t in m["trials"].items():
        print(f"In {phase} ({t['n']}): {', '.join(t['drugs'][:8])}")
    if m["stopped"]:
        print(f"⚠ Stopped trials (efficacy or safety): {', '.join(m['stopped'])}")

    def risk(c):
        r = c["risk_reasons"]
        t = []
        if "contraindicated" in r:
            t.append(f"contraindicated in {r['contraindicated']} neighbour or similar diseases and indicated in {r['indicated']}")
        if "checkpoint" in r:
            t.append("blocks immune checkpoints in an autoimmune disease")
        return "⚠ possible risk: " + "; ".join(t) if t else ""

    print("\nNew candidates (LambdaRank ranking with specificity; not a probability):")
    if m["little_info"]:
        print("⚠ There is little information about this disease in the graph: the candidates are few and unreliable.")
    for i, c in enumerate(m["candidates"], 1):
        status = "" if c["status"] in ("approved", "unknown") else f"  [{c['status']}]"
        print(f"\n{i:>3}. {c['name']}  [{c['score']:.2f}]{status}")
        if risk(c):
            print(f"     {risk(c)}")
        for e in c["expl"]:
            print(f"     ↳ {explanation_text(e)}")
    if m["hidden"]:
        print(f"\n({m['hidden']} candidates of the first {a.top} hidden for lack of sufficient evidence.)")
    print("\nDrugs without a registered indication (model without popularity; low reliability):")
    for c in m["no_use"]:
        print(f"  · {c['name']}: " + " | ".join(explanation_text(e) for e in c["expl"][:2]))
    if m["excluded"]:
        print(f"    Excluded for possible risk: {', '.join(m['excluded'])}")
    if m["extensions"]:
        print("\nAlready indicated for a form or neighbour disease: " + "; ".join(f"{f} ({', '.join(ms)})" for f, ms in m["extensions"]))
    print("\nNo flag does not mean safe. These are hypotheses for research, not medical recommendations.")


def disease_metapaths(g, a, d):
    """Experimentation mode (--paths): score by metapaths, without the production rules."""
    paths = read_paths(a.paths, COMBINED)
    scores = g.scores([d], paths)
    total = g.combine(scores, "max")[:, 0]
    indicated = set(g.neighbors(d, INDICATION))
    print(f"== {g.name(d)}  (#{d}) · metapath mode: {', '.join(paths)} ==")
    n = 0
    for j in np.argsort(-total):
        f = g.drugs[j]
        if total[j] <= 0 or n >= a.top:
            break
        if f in indicated:
            continue
        n += 1
        print(f"{n:>3}. {g.name(f)}  [{total[j]:.3f}]  " + " | ".join(f"{k} {scores[k][j, 0]:.3f}" for k in scores if scores[k][j, 0] > 0))


def targets(g, a):
    """Which genes/proteins to go after for a disease, in which direction and with which modality."""
    from evaluation.targets import Targets, TARGETS_MODEL, prioritise
    if not TARGETS_MODEL.exists():
        raise SystemExit("Target model missing: py -3 -m evaluation.targets train")
    d = one_node(g, a.name, "disease")
    D = Targets()
    t = prioritise(D, d, a.top, a.with_evidence)
    print(f"== Targets for {g.name(d)}  (#{d}) ==")
    print("Ordered by the target model (docs/TARGETS.md). «already target here»: an approved drug for "
          "this disease already acts on it. Direction: what human genetics asks for (activate or inhibit the gene).")
    already = t[t["already_target_here"]]
    if len(already):
        print(f"\nAlready exploited (target of some approved drug for this disease), {len(already)} among the first: "
              + ", ".join(f"{g_} [{p_:.2f}]" for g_, p_ in zip(already.gene, already["score"])))
    t = t[~t["already_target_here"]]
    if a.with_evidence:
        t = t[t["has_evidence"]]
    print(f"\nNew targets (top {len(t)}{', only with direct evidence in the disease' if a.with_evidence else ''}):")
    for i, (_, r) in enumerate(t.iterrows(), 1):
        print(f"\n{i:>3}. {r['gene']}  [{r['score']:.2f}]")
        print(f"     evidence: {r['evidence']}")
        print(f"     direction: {r['direction']}")
        print(f"     can be targeted with: {r['modalities']}")
        if r["safety"]:
            print(f"     ⚠ known risks: {r['safety']}")
        if r["LOEUF"] is not None and r["LOEUF"] == r["LOEUF"] and r["LOEUF"] < 0.35:
            print(f"     ⚠ gene intolerant to loss of function (LOEUF {r['LOEUF']}): inhibiting it may be risky")
        if r["essential"] == 1 or (r["dependent_fraction"] is not None and r["dependent_fraction"] > 0.5):
            print(f"     ⚠ essential gene (DepMap): {r['dependent_fraction']:.0%} of cell lines depend on it; "
                  "inhibiting it may be toxic")
        if r["tissue_specificity"] is not None and r["tissue_specificity"] > 0.5:
            print("     expression specific to few tissues (less risk of whole-body effects)")
        if not a.no_chembl and r["uniprot"]:
            from sources.chembl_api import compounds, summary as chembl_summary
            try:
                q = chembl_summary(compounds(r["uniprot"]))
            except Exception as e:  # no network: not essential
                q = f"(could not query ChEMBL: {e.__class__.__name__})"
            if q:
                print(f"     chemistry: {q}")
        if r["drugs acting on it"]:
            print(f"     drugs acting on it: {r['drugs acting on it']}")
        dir_, drugs_ = r["direction"], r["drugs acting on it"]
        if dir_.startswith("activate") and "[inhibits]" in drugs_ or dir_.startswith("inhibit") and "[activates]" in drugs_:
            print("     ⚠ the drugs acting on it go in the opposite direction to what genetics asks for")
    print("\nThese are hypotheses to prioritise research, not medical recommendations.")


def path(g, a):
    x, y = one_node(g, a.source), one_node(g, a.target)
    T = g.total()
    degree = np.diff(T.indptr)
    et = lambda u, v: g.labels[int(T[u, v])].split(":")[-1]
    Nx, Ny = set(g.neighbors_all(x)) - {y}, set(g.neighbors_all(y)) - {x}
    print(f"{g.name(x)} ({g.types[x]}) → {g.name(y)} ({g.types[y]})\n")
    if T[x, y]:
        print(f"Direct: {g.name(x)} -[{et(x, y)}]- {g.name(y)}\n")

    shared = sorted(Nx & Ny, key=lambda m: degree[m])
    print(f"2-step paths: {len(shared)}")
    for m in shared[:a.top]:
        print(f"  {g.name(x)} -[{et(x, m)}]- {g.name(m)} ({g.types[m]}) -[{et(m, y)}]- {g.name(y)}")

    ax, ay = np.fromiter(Nx, int), np.fromiter(Ny, int)
    if len(ax) and len(ay):
        sub = T[ax][:, ay].tocoo()
        weight = 1.0 / (degree[ax[sub.row]] * degree[ay[sub.col]])
        best = np.argsort(-weight)[:a.top]
        print(f"\n3-step paths: {sub.nnz:,} (most specific first)")
        for k in best:
            m1, m2 = ax[sub.row[k]], ay[sub.col[k]]
            print(f"  {g.name(x)} -[{et(x, m1)}]- {g.name(m1)} ({g.types[m1]}) -[{et(m1, m2)}]- "
                  f"{g.name(m2)} ({g.types[m2]}) -[{et(m2, y)}]- {g.name(y)}")


def evaluate(g, sample, paths, positives_of, exclude_of, pool, strict=False, base=None, save=False):
    """AUROC per metapath over a sample of diseases.

    positives_of(d) / exclude_of(d): indices of positive drugs / drugs to ignore.
    pool: mask over g.drugs of the candidates that count.
    base(d): reference score (without the graph) to compare with.
    """
    drug_pos = {f: j for j, f in enumerate(g.drugs)}
    names = None
    auc, top20, best_rank, rows = {}, {}, [], []
    for i in range(0, len(sample), 64):
        block = sample[i:i + 64]
        scores = g.scores(block, paths, strict)
        scores["combined (rank)"] = g.combine({k: scores[k] for k in paths if k in scores}, "rank")
        scores["combined (max)"] = g.combine({k: scores[k] for k in paths if k in scores}, "max")
        if names is None:
            names = list(scores) + (["reference"] if base else [])
            auc = {k: [] for k in names}
            top20 = {k: [] for k in names}
        for c, d in enumerate(block):
            positives = np.zeros(len(g.drugs), bool)
            positives[[drug_pos[f] for f in positives_of(d)]] = True
            valid = pool.copy()
            valid[[drug_pos[f] for f in exclude_of(d)]] = False
            valid |= positives & pool
            y = positives[valid]
            if not y.any() or y.all():
                continue
            for k in names:
                p = (base(d) if k == "reference" else scores[k][:, c])[valid]
                auc[k].append(auroc(p, y))
                # pessimistic rank: ties count as last
                best = rankdata(-p, method="max")[y].min()
                top20[k].append(best <= 20)
                if k == "combined (rank)":
                    best_rank.append(best)
            if base:
                feats = [scores[k][valid, c] for k in paths if k in scores] + [base(d)[valid]]
                rows.append((np.column_stack(feats).astype(np.float32), y))
        print(f"\r  {min(i + 64, len(sample))}/{len(sample)}", end="", flush=True)

    print(f"\n\nOver {len(best_rank)} diseases. AUROC: 0.5 = random, 1 = perfect. "
          "Top 20: % of diseases with at least one positive among the first 20.")
    print(f"  {'':24s} {'AUROC':>6s} {'top 20':>7s}")
    for k in names:
        desc = METAPATHS[k][3] if k in METAPATHS else ""
        print(f"  {k:24s} {np.mean(auc[k]):6.3f} {np.mean(top20[k]):7.1%}   {desc}")
    print(f"\nCombined (rank): median rank of the best positive: {np.median(best_rank):.0f} of {pool.sum():,}")
    if rows:
        model_cv(rows, [k for k in paths if k in scores] + ["popularity"], save=save)


def model_cv(rows, names, seed=0, save=False):
    """Logistic regression over the metapaths + popularity, validated on two
    halves of the diseases (trained on one and evaluated on the other)."""
    from sklearn.linear_model import LogisticRegression

    rng = np.random.default_rng(seed)
    idx = rng.permutation(len(rows))
    halves = [idx[: len(idx) // 2], idx[len(idx) // 2:]]

    def fit(index, cols):
        X = np.concatenate([rows[i][0][:, cols] for i in index])
        y = np.concatenate([rows[i][1] for i in index])
        sel = y | (rng.random(len(y)) < 0.1)  # all positives, 10 % of the negatives
        X, y = X[sel], y[sel]
        scale = np.array([X[X[:, j] > 0, j].mean() if (X[:, j] > 0).any() else 1.0 for j in range(X.shape[1])])
        clf = LogisticRegression(class_weight="balanced", max_iter=2000).fit(np.arcsinh(X / scale), y)
        return clf, scale

    variants = {
        "model: graph + popularity": list(range(len(names))),
        "model: graph only": list(range(len(names) - 1)),
    }
    print("\nLearned model (logistic regression, trained on one half of the diseases and evaluated on the other):")
    for name, cols in variants.items():
        aucs, tops = [], []
        for train_, test in [(halves[0], halves[1]), (halves[1], halves[0])]:
            clf, scale = fit(train_, cols)
            for i in test:
                X, y = rows[i]
                p = clf.decision_function(np.arcsinh(X[:, cols] / scale))
                aucs.append(auroc(p, y))
                tops.append(rankdata(-p, method="max")[y].min() <= 20)
        print(f"  {name:24s} {np.mean(aucs):6.3f} {np.mean(tops):7.1%}")
    if save:
        clf, scale = fit(idx, variants["model: graph only"])
        model = {"paths": names[:-1], "coef": clf.coef_[0].tolist(),
                 "intercept": float(clf.intercept_[0]), "scale": scale.tolist()}
        MODEL.write_text(json.dumps(model, ensure_ascii=False, indent=1), "utf-8")
        print(f"«graph only» model saved to {MODEL}")
    clf, _ = fit(idx, variants["model: graph + popularity"])
    print("\nWeight of each signal in the full model (scaled features; larger = more useful):")
    for n_, c in sorted(zip(names, clf.coef_[0]), key=lambda t: -abs(t[1])):
        print(f"  {n_:24s} {c:+.2f}")


def validate(g, a):
    """No need to hide anything: the metapaths do not use the indications of the
    disease itself, so we can measure whether they recover them."""
    paths = read_paths(a.paths, COMBINED)
    rng = np.random.default_rng(a.seed)
    ind_degree = np.diff(g.M(INDICATION).indptr)
    gene_degree = np.diff(g.M(GENES_ALL).indptr)
    candidates = np.flatnonzero((g.types == "disease") & (ind_degree > 0) & (gene_degree > 0))
    sample = rng.choice(candidates, min(a.n, len(candidates)), replace=False)
    pool = np.ones(len(g.drugs), bool) if a.drugs == "all" else g.origin[g.drugs] == "PrimeKG"
    print(f"Diseases with indicated drugs and genes: {len(candidates):,}; validating {len(sample)}"
          f"{' (strict mode)' if a.strict else ''}. Candidates: {pool.sum():,} drugs.")

    # Reference without the graph: how «popular» the drug is (indications in other diseases).
    popular = ind_degree[g.drugs].astype(float)
    ind = g.M(INDICATION)

    def base(d):
        p = popular.copy()
        p[np.isin(g.drugs, g.neighbors(d, INDICATION))] -= 1
        return p

    evaluate(g, sample, paths,
             positives_of=lambda d: g.neighbors(d, INDICATION),
             exclude_of=lambda d: g.neighbors(d, OFF_LABEL),
             pool=pool, strict=a.strict, base=base, save=a.save)
    if not a.strict:
        print("\nCareful: PrimeKG has many nearly identical disease subtypes, and this inflates "
              "«similar» and «ontology». Try --strict.")


def temporal(g, a):
    """With what was known before the cutoff year, could the drugs that entered
    trials (phase ≥ 2) for each disease afterwards be predicted?"""
    T = a.cutoff
    g.register("genes<T", g.matrix_until("ot_assoc", ["genetic_association", "genetic_literature", "somatic_mutation"], T))
    g.register("clinical<T", g.matrix_until(CLINICAL, ["approved", "phase 3", "phase 2", "phase 1", "preclinical", "unknown"], T))
    METAPATHS.update({
        "target genetics<T": ("direct", TARGET, "genes<T", f"target = genetic gene known before {T}"),
        "interaction genetics<T": ("ppi", TARGET, "genes<T", f"target ↔ genetic gene known before {T}"),
        "similar clinical<T": ("similar", "clinical<T", "genes<T", f"tested before {T} in a disease with shared genes"),
        "ontology clinical<T": ("ontology", "clinical<T", None, f"tested before {T} in a neighbour disease"),
        "interaction STRING<T": ("ppi", TARGET, "genes<T", f"target ↔ (STRING) genetic gene known before {T}", "string"),
    })
    paths = ["target genetics<T", "interaction genetics<T", "interaction STRING<T", "similar clinical<T",
             "ontology clinical<T"]

    atr = g.attributes()
    clin = atr[atr.relation == CLINICAL]
    after = clin[(clin["year"] >= T) & clin.display_relation.isin(["phase 2", "phase 3", "approved"])]
    positives = after.groupby("y_index").x_index.apply(np.asarray).to_dict()
    all_clin = clin.groupby("y_index").x_index.apply(np.asarray).to_dict()

    # Drugs that existed before the cutoff: they have targets or some previous trial.
    target_degree = np.diff(g.M(TARGET).indptr)[g.drugs]
    clin_degree = np.diff(g.M("clinical<T").indptr)[g.drugs]
    pool = (target_degree > 0) | (clin_degree > 0)
    in_pool = set(g.drugs[pool])

    gene_degree = np.diff(g.M("genes<T").indptr)
    candidates = np.array([d for d, p in positives.items()
                           if gene_degree[d] > 0 and any(f in in_pool for f in p)])
    rng = np.random.default_rng(a.seed)
    sample = rng.choice(candidates, min(a.n, len(candidates)), replace=False)
    print(f"Cutoff {T}: {len(candidates):,} diseases with prior genetic genes and new drugs in phase ≥ 2 "
          f"afterwards; validating {len(sample)}. Candidates: {pool.sum():,} drugs.")

    popular = clin_degree.astype(float)  # how many diseases the drug had been tested in before the cutoff
    evaluate(g, sample, paths,
             positives_of=lambda d: positives.get(d, []),
             exclude_of=lambda d: all_clin.get(d, []),
             pool=pool, base=lambda d: popular)
    print("\nPositive = drug that started trials for the disease from the cutoff on and reached phase ≥ 2 there.\n"
          "Reference = number of diseases in which the drug had already been tested before (popularity).")


def main():
    sys.stdout.reconfigure(encoding="utf-8")
    p = argparse.ArgumentParser(description="Drug graph explorer")
    s = p.add_subparsers(dest="command", required=True)
    s.add_parser("summary")
    c = s.add_parser("search"); c.add_argument("text"); c.add_argument("--type"); c.add_argument("--top", type=int, default=30)
    c = s.add_parser("drug"); c.add_argument("name")
    c = s.add_parser("disease"); c.add_argument("name"); c.add_argument("--top", type=int, default=20); c.add_argument("--paths")
    c = s.add_parser("targets"); c.add_argument("name"); c.add_argument("--top", type=int, default=15)
    c.add_argument("--with-evidence", action="store_true", help="only genes with direct evidence in the disease")
    c.add_argument("--no-chembl", action="store_true", help="do not query ChEMBL compounds (no network)")
    c = s.add_parser("path"); c.add_argument("source"); c.add_argument("target"); c.add_argument("--top", type=int, default=12)
    c = s.add_parser("validate"); c.add_argument("--n", type=int, default=800); c.add_argument("--seed", type=int, default=0)
    c.add_argument("--paths"); c.add_argument("--strict", action="store_true")
    c.add_argument("--drugs", choices=["primekg", "all"], default="primekg")
    c.add_argument("--save", action="store_true", help="save the «graph only» model for the disease command")
    c = s.add_parser("temporal"); c.add_argument("--cutoff", type=int, default=2015); c.add_argument("--n", type=int, default=800)
    c.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    {"summary": summary, "search": search, "drug": drug, "disease": disease,
     "path": path, "targets": targets, "validate": validate, "temporal": temporal}[a.command](Graph(), a)


if __name__ == "__main__":
    main()
