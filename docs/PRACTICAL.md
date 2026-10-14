# Practical information per candidate (I3)

Protocol: `docs/IMPROVEMENT_PLAN.md` (I3). **Display only**: none of these fields enters the ranking, the website
sections or any evaluation. Code: `sources/practical.py` (sources → `data/ext/practical.parquet`) and
`web/build.py` `practical_tables` / `DiseasePages.practical_info` (per candidate on each page).

## Sources (downloaded 14-10-2026, `data/raw/practical/VERSION`)

| Source | What is taken | Licence |
|---|---|---|
| FDA Orange Book (products, patents, exclusivity) | US generic = a marketed (Rx or OTC, not discontinued) ANDA product of the same active moiety; latest patent or exclusivity date of a marketed NDA product; routes | US public domain |
| FDA Purple Book (September 2026 monthly file) | US biosimilar = a licensed, not discontinued 351(k) product of the same proper name (four-letter suffix removed) | US public domain |
| EMA medicines table | EU generic / biosimilar = an authorised human medicine of that INN flagged as generic / biosimilar | © European Medicines Agency, reuse with acknowledgement (legal page) |
| Open Targets 26.06 `clinical` rows (already in the graph) | number of trials for the drug and the disease (or its grouped twin terms), any phase | CC0 |
| PubTator3 (already in the graph) | number of articles mentioning both the drug and the disease | public domain |

**Underexplored** = shown candidate or direct extension with no trial and at most 2 co-mentioning articles.
`exclusivity_until` is left blank when a US generic or biosimilar exists (otherwise it is the date of a newer
formulation, e.g. imatinib).

## Matching

By active moiety (`evaluation.twins.active_moiety`: name without salts or esters, lowercase), on the graph's drug
name and the Open Targets name and ChEMBL synonyms of the same molecule (INN/USAN differences such as cromoglicic
acid / cromolyn sodium). Not by drug family, so an enantiomer is not credited with its racemate's generics.
Single-ingredient products only. Two problems found in the manual check below were fixed before release:

- mineral salts reduced to the bare ion (calcium citrate was credited with calcium chloride and gluconate
  injections): a salt whose moiety would be an ion now keeps its full name;
- synonyms taken from trial registrations (AACT) include combination partners and bare ions (tipiracil was
  credited with trifluridine eye drops): only ChEMBL synonyms are used.

Match rates (15,498 graph drugs; 2,596 approved in DrugBank):

| Field | All drugs | Approved |
|---|---|---|
| found in the Orange Book | 14.0 % | 59.5 % |
| found in the Purple Book | 3.5 % | 7.2 % |
| found in the EMA table | 6.9 % | 23.2 % |
| US generic | 7.1 % | 36.1 % |
| US biosimilar | 0.3 % | 0.7 % |
| EU generic (centralised) | 0.8 % | 4.2 % |
| EU biosimilar | 0.3 % | 0.9 % |

The EU figures are low by construction: the EMA table covers the centralised procedure only, and most generics of
old drugs are authorised country by country. **«No EU generic» does not mean there is none.**

## Manual check (50 approved drugs, random seed 14)

Each row compared with the source files and general knowledge of the product (after the two fixes above).

| Drug | Orange / Purple / EMA | Shown | Check |
|---|---|---|---|
| Acetyl sulfisoxazole | yes / – / – | nothing (discontinued) | correct |
| Iodine | – / – / – | nothing | correct (no single-ingredient product) |
| Hydroxocobalamin | yes / – / yes | US generic; injection | correct |
| Etafedrine | – / – / – | nothing | correct (not marketed in US or centrally in EU) |
| Lactose | – / – / – | nothing | correct (excipient) |
| Chromous sulfate | – / – / – | nothing | correct |
| Zinc sulfate | yes / – / – | US generic; intravenous | correct |
| Cefditoren | – / – / – | nothing | correct for display (US product discontinued; the pivoxil ester is not stripped, so it is not even found) |
| Iloperidone | yes / – / – | protected until 2031; oral | correct |
| Parecoxib | – / – / yes | nothing | correct (EU brand only) |
| Cromoglicic acid | yes / – / – | US generic; inhalation, nasal, ophthalmic, oral | correct (found through the USAN synonym) |
| Dexchlorpheniramine maleate | yes / – / – | US generic; oral | correct |
| Berberine | – / – / – | nothing | correct |
| Axitinib | yes / – / yes | US and EU generic; oral | consistent with the sources (marketed ANDA; EU generic authorised) |
| Dutasteride | yes / – / – | US generic; oral | correct (EU generics are national) |
| Papain | – / – / – | nothing | correct |
| Cefotaxime | yes / – / – | nothing; injection | correct (all US products discontinued) |
| Glimepiride | yes / – / – | US generic; oral | correct |
| Tamoxifen | yes / – / – | US generic; oral | correct |
| Reteplase | – / yes / – | nothing | correct (no biosimilar) |
| Ravulizumab | – / yes / yes | nothing | correct |
| Sodium phosphate, monobasic | – / – / – | nothing | correct (combination products only) |
| Dabrafenib | yes / – / yes | protected until 2038; oral | caveat: latest date of any marketed NDA product (includes the newer dispersible tablets) |
| Belimumab | – / yes / yes | nothing | correct |
| Ferric derisomaltose | yes / – / – | protected until 2036; intravenous | correct |
| Rocuronium | yes / – / – | US generic; injection, intravenous | correct |
| Tipiracil | – / – / – | nothing | correct (combination only) |
| Chlortetracycline | yes / – / – | nothing | correct (discontinued) |
| Topiramate | yes / – / – | US generic; oral | correct |
| Gluconic acid | – / – / – | nothing | correct |
| Selenomethionine Se-75 | yes / – / – | nothing | correct (discontinued) |
| Olodaterol | yes / – / – | protected until 2027; inhalation | correct |
| Estradiol dienanthate | yes / – / – | US generic; injection, oral, transdermal, vaginal | **misleading**: credited with estradiol products (same active moiety in the FDA sense, but no dienanthate product exists in the US) |
| Pentaerythritol tetranitrate | – / – / – | nothing | correct |
| Dyphylline | yes / – / – | nothing | correct (no marketed ANDA) |
| Calcium citrate | – / – / – | nothing | correct (after the ion fix) |
| Acalabrutinib | yes / – / yes | protected until 2036; oral | correct |
| Ferric sulfate | – / – / – | nothing | correct (after the ion fix) |
| Somatostatin | – / – / – | nothing | correct |
| Abatacept | – / yes / yes | nothing | correct |
| Oxetacaine | – / – / – | nothing | correct |
| Mecasermin | – / yes / yes | nothing | correct |
| Phendimetrazine | yes / – / – | US generic; oral | correct |
| Mizolastine | – / – / – | nothing | correct |
| Glucosamine | – / – / – | nothing | correct |
| Aspartic acid | – / – / – | nothing | correct |
| Mannitol | yes / – / yes | US generic; inhalation, injection | correct |
| Minoxidil | yes / – / – | US generic; oral, topical | correct |
| Abametapir | yes / – / – | nothing | correct (no marketed product) |
| Netarsudil | yes / – / yes | protected until 2034; ophthalmic | correct |

Result: 48 of 50 correct as displayed, 1 misleading (an ester credited with the parent's products), 1 caveat
(exclusivity date of a newer formulation). Before the two fixes, 3 more rows were wrong (calcium citrate, ferric
sulfate, tipiracil) and cromoglicic acid was not found.

## Known limitations

- Esters and salts are credited with the parent moiety's products (estradiol dienanthate); the routes column shows
  the routes of those products.
- The exclusivity date is the latest of any marketed brand product of the moiety and can belong to a newer
  formulation; it is a hint, not a legal opinion on freedom to operate.
- EU coverage is the centralised procedure only.
- Trial and article counts come from the graph sources (Open Targets clinical records, PubTator3 co-mentions) and
  miss trials registered only elsewhere and articles not annotated by PubTator3.
