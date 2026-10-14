#!/bin/bash
# Usage (WSL): bash evaluation/kge_bench.sh grid | final        (see evaluation/kge_bench.py)
# Resumable: a run whose zero_shot / area prediction file exists is skipped.
REPO="$(cd "$(dirname "$0")/.." && pwd)"
PY=~/pykeen-env/bin/python
cd "$REPO"
OUT=data/eval/kge
mkdir -p $OUT
run() {  # split model dim lr dest
  if [ -f "$5/$1.parquet" ]; then echo "$1 $2 d$3 lr$4: already done"; return; fi
  $PY -u -m evaluation.kge_bench train --split $1 --model $2 --dim $3 --lr $4 --dest $5 >> $OUT/log_$2.txt 2>&1
  echo "$1 $2 d$3 lr$4: exit $? $(date +%T)"
}
if [ "$1" = "grid" ]; then
  for m in TransE ComplEx RotatE; do
    for d in 128 256; do
      for lr in 0.001 0.01; do
        name=$($PY -c "print(f'${m}_d${d}_lr{float(\"$lr\"):g}')")
        run zero_shot_1 $m $d $lr $OUT/grid/$name
      done
    done
  done
elif [ "$1" = "final" ]; then
  [ -f $OUT/selected.json ] || { echo "run: py -3 -m evaluation.kge_bench select"; exit 1; }
  for m in TransE ComplEx RotatE; do
    read d lr < <($PY -c "import json; s=json.load(open('$OUT/selected.json'))['$m']; print(s['dim'], s['lr'])")
    dest=benchmark/predictions/kge_$(echo $m | tr A-Z a-z)
    mkdir -p $dest
    src=$OUT/grid/$($PY -c "print(f'${m}_d${d}_lr{float(\"$lr\"):g}')")
    cp -n $src/zero_shot_1.parquet $src/external_ot2606.parquet $dest/  # same training split as the grid
    for s in 2 3 4 5; do run zero_shot_$s $m $d $lr $dest; done
    for a in cell_proliferation mental_health cardiovascular anemia adrenal_gland autoimmune metabolic_disorder diabetes neurodigenerative; do
      run area_$a $m $d $lr $dest
    done
  done
fi
