# Data licence

The **code** in this repository is under the MIT licence (`LICENSE`).

The **derived data and results** (the website data in `web/public/data/`, the benchmark task files and
predictions, and any ranking produced by this code from the sources below) are released under
**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**,
because they are adaptations of sources with non-commercial and share-alike terms:

| Source | Licence |
|---|---|
| PrimeKG (Chandak et al., Scientific Data 2023): distributed CC0 by Harvard, but it aggregates sources with their own licences, including: | |
| · DrugBank (indications, targets, status, classes) | CC BY-NC 4.0 |
| · SIDER (adverse effects), DisGeNET (gene–disease associations) | CC BY-NC-SA 4.0 |
| Open Targets Platform 26.06 (includes ChEMBL mechanisms and molecules) | CC0 1.0 |
| STRING 12.5 | CC BY 4.0 |
| PubTator3 (NCBI) | Public domain |
| MONDO / HGNC | CC BY 4.0 / CC0 |
| HPO (Human Phenotype Ontology) | HPO licence (free use with attribution; https://hpo.jax.org/license) |

Not covered by the published data:

- **CTD** (Comparative Toxicogenomics Database, non-commercial) was used only in analyses; nothing from
  CTD is in the production model, the website data or the benchmark (`docs/DECISIONS.md`, P20).
- **ChEMBL web services** (CC BY-SA 3.0) are queried on demand by the command-line `explore.py targets`
  to list compounds binding a target; those results are not part of the published data.
- TxGNN files reused here are MIT-licensed; see `THIRD_PARTY_NOTICES.md`.

Not for commercial use. Not a medical device and not medical advice: research use only.
Please cite the project and the original sources.
