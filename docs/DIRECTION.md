# Direction of the effect: does it improve or worsen the disease?

Code: `evaluation/direction.py`. Model: `py -3 -m evaluation.direction train`
(`data/prep/direction_model.pkl`). The `explore.py disease` command uses it to:

- label each candidate with «no risk signal», «uncertain» or «⚠ risk of worsening it» (with the reasons);
- **remove** from the list of drugs with no prior use those that would probably worsen the disease
  (e.g. PD-1 inhibitors for multiple sclerosis, which would activate the immune system), and
  say which ones were removed.

A high score **does not mean the drug works**: only that there are no risk signals.

## Evaluation

«Indication versus contraindication» classifier on the 39,889 PrimeKG pairs (23 %
indications), cross-validation on two halves of the diseases. AUROC (0.5 = chance):

| Signals | All | Drugs with ≤ 1 indication |
|---|---|---|
| Genetics (ChEMBL mechanism vs genetic direction) | 0.50 | 0.50 |
| Adverse effects (SIDER, FAERS) that are phenotypes of the disease | 0.50 | 0.50 |
| Analogues (drugs approved in neighbouring diseases on the same targets, same or opposite direction) | 0.62 | 0.53 |
| **All**, with «indicated / contraindicated in neighbouring or similar diseases» | **0.97** | **0.91** |

Reading:

- **Genetic direction exists for only 1 % of the pairs**: it is not yet useful in general.
- **Adverse effects do not discriminate**, because of confounding by indication: a disease's drugs
  often have its symptoms as an «adverse effect» (FAERS records the symptoms of the people taking them).
- **The strong signal is contraindications in neighbouring diseases**: contraindications tend to
  affect whole families of diseases. This is valid knowledge for warning about risks, but it means
  that for very new drugs (with no recorded contraindications) the model knows little.
- No variable uses the relation of the same pair or «minus this pair» counts (see
  the leak described in `docs/TARGETS.md`).

## General validation of the labels (all diseases, out of sample)

Fraction of real contraindications by label (baseline: 0.77):

| Label | Pairs | Contraindications |
|---|---|---|
| «risk» (P < 0.3 **and** a concrete reason) | 26,686 | **0.98** |
| P < 0.3 without any reason | 3,613 | 0.82 (almost the baseline: that is why it is not flagged as a risk) |
| «uncertain» | 1,597 | 0.54 |
| «no risk signal» (P ≥ 0.7) | 7,993 | 0.08 |

## Explicit rule: checkpoint blockers in autoimmune diseases

The data do not cover the known risk that immune checkpoint inhibitors (PD-1,
PD-L1, CTLA-4, LAG-3, TIM-3, TIGIT) trigger or worsen autoimmune diseases. It is flagged with
a rule (`checkpoint_rule`), validated on all diseases:

- **First version (wrong)**: Open Targets' «immune system disorder». It includes lymphomas, leukaemias,
  myeloma and thymoma, where these drugs are in fact used (196 trials). Fixed.
- **Current version**: diseases under «autoimmune disease» (123). No indication or
  contraindication in PrimeKG; the only trials are 6 phase 1 trials of nivolumab (probably patients
  with cancer and a pre-existing autoimmune disease) and one of unknown phase.
- Agonists of these targets (peresolimab, rosnilimab…) and fusion proteins such as
  abatacept are excluded from it.

Lesson: a rule built from a single case must be validated on all diseases before adopting it.

## Correction after the independent review (11-10): the model is no longer used in production

The direction model is, in practice, «copying the labels of neighbouring diseases»: a one-line
rule (indicated minus contraindicated in neighbouring and similar diseases) gives AUROC 0.953 (the model, 0.969), and
**in the pairs without any label for the drug in neighbouring diseases (9.7 %) the model gives 0.447, worse than
chance**. That is exactly the case of new candidates. Verified independently.

Decision: in production (website and `explore.py`) **no probability** of improvement is shown. Only a
**possible-risk warning** with a concrete reason (`risk_flags`):

| Rule | Pairs | Real contraindications |
|---|---|---|
| contraindicated in neighbouring/similar diseases > 0 | 74 % | 90.7 % |
| **contraindicated > indicated (in neighbouring/similar diseases)** (the one used) | 66 % | **97.4 %** |
| contraindicated > 0 and indicated = 0 | 53 % | 99.3 % |
| (baseline) | | 76.7 % |

Plus the checkpoint rule for autoimmune diseases. The absence of a warning does not mean safe.
