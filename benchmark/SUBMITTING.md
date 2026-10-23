# Submitting a model to OREVID-bench

Anyone can add a model to the [leaderboard](LEADERBOARD.md). The evaluation is always `bench.py evaluate` from this
repository's `main` branch, on the task files in `tasks/`, beside the popularity baseline. Rules fixed on
23-10-2026, before any outside submission.

## 1. Rules

1. **Training data.** For each task, use only the training edges of PrimeKG (`train_edges(kg, "tasks", task)`,
   README «Training edges») plus outside data that does **not** encode drug–disease therapeutic relations: no
   clinical-trial registries, drug labels, indication databases (DrugBank, ChEMBL, Open Targets indications…) or
   literature-mined «drug treats disease» relations. Gene, protein, pathway, chemistry and expression data are fine.
2. **One method, retrained per task.** The same code and settings for every task; the model is trained again on each
   task's training edges (5 zero-shot seeds, 9 areas, external exam).
3. **Settings** (hyperparameters, early stopping, feature choice) are chosen on the validation edges or on
   `zero_shot_1` only, and declared.
4. **External exam.** Its labels (`tasks/external_ot2606.json`) are public, so it is honour-based: a model that used
   any source with drug–disease relations newer than PrimeKG (2023) gets no external row.
5. **Open.** The code (public repository) and the prediction files (a record with a DOI, e.g. Zenodo or a Hugging
   Face dataset) are public. A paper is not needed.
6. Results are what `bench.py` reports; the maintainers do not round, select seeds or drop tasks. A submission that
   does not pass `bench.py validate` is not evaluated.

## 2. Prediction files

One file per task, named `<task>.parquet` (or `.csv`), columns `disease`, `drug` (PrimeKG `node_index`) and `score`
(higher = more likely to treat), with every drug of `tasks/<task>.json` → `drugs` scored for every test disease. All
15 task files go in one zip, without folders inside. Check them first:

```bash
python bench.py validate --pred my_model/
python bench.py evaluate --pred my_model/
```

## 3. How to submit

Open a pull request that adds **only** `submissions/<model_id>/submission.json` (lower-case letters, digits and `_`),
filled in from [`submissions/TEMPLATE.json`](submissions/TEMPLATE.json): model name, team, contact, date, code URL,
the DOI and direct download URL of the zip, its MD5, the data used beyond the training edges (with versions), and how
the settings were chosen. Leave `"accepted": false`.

The pull request runs `.github/workflows/bench-submission.yml`, which downloads the zip, checks the MD5, runs
`bench.py validate` and `bench.py evaluate` with the `main`-branch code and tasks, and attaches `results.json` and
the summary to the run. A maintainer then reads the declaration, adds `predictions/<model_id>/results.json`, sets
`"accepted": true`, regenerates `LEADERBOARD.md` (`python leaderboard.py`) and the disease-level intervals, and
merges. Questions: open an issue or write to contact@orevid.org.

## 4. Hidden-label track (planned)

Every published label can be learned from, so the strongest test is prospective: predictions deposited with a DOI
before a date and scored on uses that appear later. OREVID's own prospective snapshot (12-10-2026,
doi:10.5281/zenodo.23326593, `docs/PROSPECTIVE_PLAN.md`) follows this design; outside models will be able to join a
later round with the same plan.
