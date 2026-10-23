# G2. Genetic direction shown on the pages — release check

Protocol: `docs/IMPROVEMENT_PLAN.md`, «Registered 23-10-2026, before any change to the website: G2». 50 random lines of
the rebuilt pages (365 lines on 141 pages with the `extra_qtl` model), checked with `evaluation/release_checks.py g2`
against the Open Targets GraphQL API (`data/eval/g2_check.jsonl`, queries and answers kept). **Deviation:** the API
serves release 26.09, while OREVID's data are 26.06 (26.09 is not used because of its multiple sclerosis merge).

| Check | Lines failing |
|---|---|
| (1) a colocalisation (h4 ≥ 0.8) of a cis QTL of the gene with a GWAS of the disease, with the stated direction | 9 of 50 (7 with no such colocalisation at all in 26.09; 2 with direction 0) |
| (2) the drug's mechanism on the gene is of the stated kind | 0 of 50 |

**Not released** (gate ≤ 2). The failing lines are eQTL colocalisations from a single study that are present in 26.06
and absent or undirected in 26.09 (e.g. MTTP – fatty liver disease, SLx-4090: no colocalisation at all with that
disease in 26.09; flecainide – Brugada syndrome, SCN5A, a correct «disagrees» biologically, also missing). Single-study
eQTL colocalisations are not stable between Open Targets releases, so a line on a page cannot cite them as a source.
The direction stays in the ranking model (author's decision, P61) and in the internal graph; showing it would need a
rule restricted to stable signals (e.g. pQTL, or support in two releases), registered and checked anew.

## G2b: only lines stable between 26.06 and 26.09 (24-10-2026)

Of the 365 lines of the pages, 262 are stable (26.09 has a colocalisation with the stated direction and none with the opposite one); 64 have an opposite colocalisation in 26.09 and 39 no support. Check of 50 stable lines read by the assistant (`data/eval/g2b_check.json`): the drug's mechanism right in 50/50, but in **7 of 50 the supporting GWAS is not the page's disease**: neuroticism mapped to «neurotic disorder» (twice), a medication questionnaire on the diabetes page, «asthma or reflux» on the reflux page, mood traits on the anxiety-disorder page, cardiac arrhythmias on the COPD page, and Barrett's oesophagus / adenocarcinoma on the adenosquamous-carcinoma page. **Not released** (gate ≤ 2): the weak point is now Open Targets' mapping of some GWAS traits (quantitative traits, questionnaires, combined phenotypes) to diseases.

## G2c: stable lines resting on a case–control GWAS of the page's disease (24-10-2026) — released

203 of the 365 lines pass (stable as G2b, and at least one supporting 26.09 GWAS has cases and only the page's disease ids). Check of 50 new lines read by the assistant (`data/eval/g2c_sample.json`): mechanism right 50/50; supporting trait the page's disease or a form of it in 49/50 (wrong: «mood disorders» on the anxiety-disorder page). **1 of 50 wrong → released**: 198 lines on the pages (some candidates share a line), each with the QTL type and a link to the supporting GWAS study. The verified set is `data/ext/genetic_lines_g2c.parquet`; after a rebuild with new candidates, `py -3 -m evaluation.release_checks g2b` must be run again before `web.annotate` (new lines are not shown until checked).
