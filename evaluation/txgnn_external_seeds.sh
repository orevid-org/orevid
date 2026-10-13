#!/bin/bash
# TxGNN on the external exam with the models already trained on zero-shot seeds 2-5
# (docs/REVIEW_RESPONSE.md T2). Inference only: txgnn_wsl.py loads ~/txgnn_data/models/TxGNN_<seed>_complex_disease.
# Usage (from Windows): wsl -d Ubuntu -- bash <repo>/evaluation/txgnn_external_seeds.sh [seeds...]
export PATH=$HOME/.local/bin:$PATH
REPO="$(cd "$(dirname "$0")/.." && pwd)"
source ~/txgnn-env/bin/activate
source ~/txgnn-env/cuda_libs.sh
P="$REPO/data/eval"
SEEDS="${@:-2 3 4 5}"
for s in $SEEDS; do
  if [ -f "$P/txgnn_predictions/TxGNN_${s}_complex_disease_external.parquet" ]; then
    echo "seed $s: already done"; continue
  fi
  echo "seed $s: starts $(date +%T)"
  python -u "$REPO/evaluation/txgnn_wsl.py" --split complex_disease --seed "$s" --data ~/txgnn_data \
    --splits "$P" --out "$P/txgnn_predictions" --diseases "$P/external_diseases.csv" \
    > "$P/log_txgnn_external_$s.txt" 2>&1
  echo "seed $s: exit $? $(date +%T)"
done
