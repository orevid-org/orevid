#!/bin/bash
# Rebuilds all the evaluations and the production model after a change to the model.
#   bash evaluation/rebuild_all.sh > data/eval/log_rebuild.txt 2>&1
set -e
export PY="${PY:-python}"  # e.g. PY="py -3" on Windows
cd "$(dirname "$0")/.."
FS=$($PY -c "from evaluation.txgnn import PRODUCTION; print(PRODUCTION)")  # production feature set
echo "== zero-shot (5 seeds) $(date +%T)"
$PY -m evaluation.txgnn --paths $FS --rounds 5 --rank > data/eval/log_rebuild_complex.txt
echo "== areas $(date +%T)"
$PY -m evaluation.txgnn --split areas --seeds 1 --paths $FS --rounds 5 --rank > data/eval/log_rebuild_areas.txt
echo "== external exam $(date +%T)"
$PY -m evaluation.external_exam evaluate > data/eval/log_rebuild_external.txt
echo "== OREVID-bench $(date +%T)"
$PY -m evaluation.export_bench
$PY -m evaluation.export_bench external
(cd benchmark && for m in popularity txgnn orevid_no_popularity orevid_lambdarank; do
   $PY bench.py evaluate --pred predictions/$m > /dev/null; done && $PY leaderboard.py > /dev/null)
echo "== production model $(date +%T)"
$PY predictor.py
$PY predictor.py reference
echo "== done $(date +%T)"
