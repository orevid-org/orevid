#!/bin/bash
# Usage: train_txgnn.sh [--test] <split> <seed> [<split> <seed> ...]
# Skips the splits that already have predictions (it can be resumed).
export PATH=$HOME/.local/bin:$PATH
REPO="$(cd "$(dirname "$0")/.." && pwd)"  # this repository, as seen from WSL (/mnt/c/...)
source ~/txgnn-env/bin/activate
source ~/txgnn-env/cuda_libs.sh
mkdir -p ~/txgnn_data
[ -f ~/txgnn_data/kg.csv ] || cp "$REPO/data/raw/txgnn/kg.csv" ~/txgnn_data/kg.csv
[ -f ~/txgnn_data/node.csv ] || cp "$REPO/data/raw/txgnn/node.csv" ~/txgnn_data/node.csv
P="$REPO/data/eval"
EXTRA=""
if [ "$1" = "--test" ]; then EXTRA="--test"; shift; fi
while [ $# -gt 0 ]; do
  SUF=""; [ -n "$EXTRA" ] && SUF="_test"
  if [ -f $P/txgnn_predictions/TxGNN_$2_$1$SUF.parquet ]; then
    echo "$1 $2: already done"
  else
    echo "$1 $2: starts $(date +%T)"
    python -u "$REPO/evaluation/txgnn_wsl.py" --split $1 --seed $2 --data ~/txgnn_data \
      --splits $P --out $P/txgnn_predictions $EXTRA > $P/log_txgnn_$1_$2$SUF.txt 2>&1
    echo "$1 $2: exit $? $(date +%T)"
  fi
  shift 2
done
