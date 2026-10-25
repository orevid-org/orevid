# Preprint v3: submission to Preprints.org (new version of preprints-238326 / 202610.1918)

Files: `OREVID_preprint_v3.docx` (manuscript), `OREVID_preprint_v3.pdf` (A4), graphical abstract
`figures/figure1.png` (as for v1 and v2). Same subject, licence (CC BY 4.0), author and related links as v1.

## Changes in this version (for the «new version» field)

Version 3 reports the current model and adds three tests.

1. Model: two gene-level columns were added to the ranking model, the direction of human genetics (protein and
   expression QTLs colocalised with disease association studies, Open Targets 26.06) and drug targets among the
   disease genes listed by Orphanet. Neither met its pre-registered adoption rule; both were adopted by the author,
   and this is stated in Methods and Limitations. Every result was recomputed with this model (zero-shot
   Recall@100 0.812 against 0.412 for TxGNN; disease areas 0.627 against 0.237; external exam 0.355 against 0.121).
2. Temporal test corrected: drugs are now dated by their first approval (DrugCentral). The advantage over popularity
   at the 2015 cutoff is 0.102 (95% CI 0.078–0.126) instead of 0.148, because the first design left later
   indications in the training graph; the advantage on non-obvious uses is no longer shown.
3. New: TxGNN trained on the same dated 2015 graph ranks the later phase 3 or approved uses below popularity
   (0.117 against 0.181; OREVID 0.283).
4. New Results section «Independent exams»: CIViC oncology evidence, indications recorded only in DrugCentral or
   MED-RT, and EMA orphan designations (OREVID ahead of popularity in all three), and posted trial outcomes (the
   claim that the candidates shown meet trial end points more often was not supported).
5. Website section: prevalence (Orphanet), CPIC pharmacogenomic guidelines and verified genetic-direction lines.
6. Figures 2–4 and all tables regenerated; references 24–26 added (DrugCentral 2023, Orphadata, CIViC).

## Short abstract (if the form limits the length; the manuscript keeps the full one)

Knowledge-graph models can rank existing drugs for a disease, but their benchmarks mostly reward rediscovering known
treatments. OREVID ranks drugs for 10,042 diseases with an interpretable metapath model (LightGBM LambdaRank) on a
graph built from PrimeKG, Open Targets and STRING, and shows a candidate only when an evidence line supports it. On
TxGNN's own zero-shot splits, Recall@100 was 0.812 ± 0.038 for OREVID, 0.412 ± 0.042 for TxGNN (retrained with its
official code) and 0.358 for a popularity baseline; on nine held-out disease areas, 0.627 against 0.237; on 3,718
approvals absent from the graph, 0.355 against 0.121. Trained on the graph as known in 2015, with drugs dated by their
first approval, OREVID ranked later phase 3 or approved uses above popularity (0.283 against 0.181), while TxGNN
trained on the same graph fell below it (0.117). On labels the graph does not contain (CIViC, DrugCentral and MED-RT,
EMA orphan designations) OREVID was ahead of popularity in every case. The limits are reported with the same weight:
most of the signal comes from what already works in related diseases, the advantage on non-obvious uses is not
shown in the temporal test, the shared genes select no better than random genes of the same size, and the candidates
shown did not meet trial end points more often than hidden ones. Rankings were frozen on 12 October 2026 for a
prospective test. Code, data, the website and the benchmark (OREVID-bench) are open.
