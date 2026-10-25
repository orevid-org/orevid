#!/bin/bash
# T3 (docs/IMPROVEMENT_PLAN.md): TxGNN on the dated temporal cutoffs, inside WSL.
#   wsl -e bash /mnt/c/dev/farmagraf/evaluation/t3_txgnn.sh
export PATH=$HOME/.local/bin:$PATH
REPO="$(cd "$(dirname "$0")/.." && pwd)"
source ~/txgnn-env/bin/activate
source ~/txgnn-env/cuda_libs.sh
P="$REPO/data/eval"
for C in 2015 2012; do
  TS=0.05${C:2:2}
  if [ -f $P/txgnn_predictions/TxGNN_temporal${C}_external.parquet ]; then echo "$C: already done"; continue; fi
  echo "$C: starts $(date +%T)"
  python -u "$REPO/evaluation/txgnn_wsl.py" --label $P/split_temporal${C}_1.npy --test-size $TS --tag temporal$C \
    --data ~/txgnn_data --splits $P --out $P/txgnn_predictions --diseases $P/temporal_${C}_diseases.csv \
    > $P/log_txgnn_temporal$C.txt 2>&1
  echo "$C: exit $? $(date +%T)"
done
