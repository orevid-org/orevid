#!/bin/bash
# Rebuilds everything from the downloaded raw data (data/raw): graph, models, evaluations, benchmark and website.
# Downloads first (once): python download.py --all && python -m sources.download (see README).
# TxGNN predictions are an input here (data/eval/txgnn_predictions, from evaluation/train_txgnn.sh on a GPU).
#   bash build_all.sh > data/eval/log_build.txt 2>&1
set -e
export PY="${PY:-python}"  # e.g. PY="py -3" on Windows
cd "$(dirname "$0")"
mkdir -p data/eval
step() { echo "== $1 $(date +%T)"; }

step "sources"
$PY -m sources.opentargets
$PY -m sources.string_db
$PY -m sources.ctd
$PY -m sources.direction
$PY -m sources.gwas_direction
$PY -m sources.pubtator
$PY -m sources.structure
step "graph"
$PY prepare.py
# the PubMed pilot links its relations to graph nodes, so it needs a first graph and a second pass
$PY -m sources.pubmed import resources/pubmed_pilot_ms.jsonl
$PY prepare.py
$PY -m sources.ot_approvals
$PY -m sources.targets           # target annotations, indexed by graph node
$PY -m sources.practical build   # generic / biosimilar status (download once: sources.practical download)
step "models"
$PY explore.py validate --strict --save
$PY predictor.py
$PY predictor.py reference
$PY -m evaluation.targets train
$PY -m evaluation.direction train
step "evaluations"
bash evaluation/rebuild_all.sh
$PY -m evaluation.evidence > data/eval/log_evidence.txt
$PY -m evaluation.temporal run 2012 extra_qtl && $PY -m evaluation.temporal run 2015 extra_qtl
$PY -m evaluation.temporal run 2012 extra_no_ot && $PY -m evaluation.temporal run 2015 extra_no_ot
$PY -m evaluation.temporal report
$PY -m evaluation.calibration
step "website"
$PY -m web.build
$PY -m web.annotate      # released information that does not change the ranking (P66, P67)
$PY -m web.panel
step "tests"
$PY -m pytest tests -q
step "done"
