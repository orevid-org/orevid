# Prospective validation: analysis plan (written before any outcome is known)

Snapshot date: **12 October 2026**. Code version: OREVID v1.1.0 (git tag `v1.1.0`, Zenodo concept DOI
10.5281/zenodo.23323226). The snapshot itself is archived separately on Zenodo, doi:10.5281/zenodo.23326593
(https://doi.org/10.5281/zenodo.23326593), and is never modified. Any deviation from this plan will be reported as such.

## What was frozen

Produced by `python -m evaluation.snapshot` (files in the Zenodo snapshot record):

| File | Content |
|---|---|
| `ranking_top100.csv.gz` | For each of the 10,042 diseases on the website: the model's top 100 drugs in model order, the specificity score, and the website section of each drug (`extension`, `candidate`, `exploratory`, `not shown`) |
| `shown_outside_top100.csv` | Drugs shown on the website that rank below 100 (extensions and exploratory lists come from deeper ranks) |
| `website_data.zip` | The exact website data (every disease page as served on the snapshot date) |
| `summary.json` | Counts, date |

## Question

Do the candidates OREVID shows become real uses (new trials, new approvals) more often than drugs the
model ranked highly but did not show, and than drugs ranked lower?

## Groups (disease–drug pairs)

1. **Extensions**: website section `extension`.
2. **Candidates**: website section `candidate`; subgroups with and without ≥ 3 specific shared genes.
3. **Ranked, not shown (control 1)**: ranks 1–25 with section `not shown` (hidden for lack of evidence).
4. **Lower ranks (control 2)**: ranks 26–100.

Pairs already approved, contraindicated or in a registered trial for that disease on the snapshot date
(as listed on the archived disease page) are excluded from every group.

## Outcomes

- **Primary (12 and 24 months):** a new interventional clinical trial of phase 2 or later, first posted
  on ClinicalTrials.gov after 12 October 2026, with the drug (or the same active moiety) as an
  intervention and the disease (MONDO term or a direct child) as a condition. Trials are matched through
  the Open Targets `clinical_report` table of the first release after each time point, and spot-checked
  by hand on ClinicalTrials.gov.
- **Secondary (24 months):** a new approval of the drug for the disease in DrugBank or in Open Targets,
  counting Open Targets approvals only when confirmed as in `sources/ot_approvals.py` (single-drug label
  or regulatory report).

## Analysis

- Proportion of pairs with the outcome per group, with Wilson 95% confidence intervals.
- Primary comparison: Candidates against Ranked, not shown (control 1), as a risk ratio with a 95%
  confidence interval from a bootstrap over diseases (diseases resampled, 10,000 replicates), since pairs
  of the same disease are not independent.
- Secondary comparisons: Extensions against Candidates; specific against generic shared genes; every
  group against Lower ranks (control 2).
- All counts are reported, including null or negative results. The analysis code will be added to the
  repository before the outcome data are downloaded.

## Known limitations

Trials and approvals are slow and depend on factors unrelated to the biology (funding, patents,
commercial interest), so absolute rates will be low and the comparison between groups is what matters.
Matching drugs and diseases across sources can miss or mis-assign pairs; the same matching is used for
every group.

## Integrity

SHA-256 of the frozen files, so that any later copy can be checked against them:

- `ranking_top100.csv.gz`: `ce6bc819958f7b2a11da6480604a5a7c8a911ec80515c039ccfa7a20940664c8`
- `website_data.zip`: `df3a4a5d5cf6216fc660e21387f7c6049ab806a1d7e2ad4ee9c26e162a07f936`

Snapshot counts: 10,042 diseases; in the top 100, 16,209 extensions, 53,575 candidates, 5,587
exploratory and 928,829 not shown; 12,979 shown pairs rank below 100.
