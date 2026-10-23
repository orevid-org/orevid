# R4. Orphanet prevalence on rare-disease pages

Protocol: `docs/IMPROVEMENT_PLAN.md`, «Registered 23-10-2026, evening: R4» and «… night, after the R4 check: R4b».
Source: Orphadata product 9 (CC BY 4.0), point-prevalence class (worldwide and validated preferred), ORPHA → MONDO
(SSSOM exact) → graph node. Checked by reading the 50 Orphanet pages of each sample in the browser (orpha.net,
23-10-2026); samples in `data/eval/r4_sample.json` and `data/eval/r4b_sample.json`.

| Version | Pages with a value | Values that match orpha.net | On the wrong disease | Decision |
|---|---|---|---|---|
| R4 (all mapped nodes, known classes) | 3,320 | 50 / 50 | 3 (grouped PrimeKG nodes given a subtype's value: dystrophic epidermolysis bullosa ← nails-only form; «intellectual developmental disorder, autosomal recessive» ← a cataract–aciduria syndrome; «SEMD with joint laxity» ← Isidor type) | not released (gate ≤ 2) |
| **R4b** (single-MONDO nodes only) | 2,896 | 50 / 50 | 1 (the page «spastic quadriplegia», MONDO:0016215 = spastic quadriplegic cerebral palsy, receives the value of Orphanet's inherited congenital spastic tetraplegia; counted wrong because the page name reads as the common form) | **released** |

A first R4 sample, drawn before «Unknown» classes were hidden, also found two values the Orphanet page does not show
(Scheie syndrome, Unverricht–Lundborg syndrome: «-» on the page, a class in the file) and a broader node with a
narrower entity (proximal spinal muscular atrophy ← its autosomal dominant adult-onset form); R4b removes the latter
kind. Shown on the page as «Prevalence (Orphanet): class», with the ORPHA code and a link.
