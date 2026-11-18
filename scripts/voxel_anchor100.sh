#!/bin/bash
# voxel : test seed_anchor_xyz (graine figee) contre base, zones z65a90 et z50, 4 graines, 30 generations, 1 thread (deterministe)
ulimit -Sn 750000
PY=~/vesuvius/gg/bin/python; C=~/vesuvius/code; B=~/vesuvius/anchor100; mkdir -p $B; OUT=$B/run.out
PRED=https://vesuvius-challenge-open-data.s3.amazonaws.com/PHerc0490A/representations/predictions/surfaces/20250521151210-surface-20260413222639-surface-m7-L0-th0.2.zarr
NG=~/vesuvius/data/0490A/normal_grids
sed -n 1828,1846p ~/vesuvius/villa/volume-cartographer/core/src/GrowPatch.cpp > $B/code_1828.txt
: > $B/tasks.txt
for z in z65a90 z50 z80a90; do while read n x y z3; do for arm in base anchor; do echo "$z $n $arm $x $y $z3" >> $B/tasks.txt; done; done < ~/vesuvius/seedbench/${z}_seeds.txt; done
run1() { z=$1; n=$2; arm=$3; x=$4; y=$5; zz=$6; W=$B/$z/$arm/$n; mkdir -p $W; cd $W
  A=""; [ $arm = anchor ] && A=", \"seed_anchor_xyz\": [$x, $y, $zz]"
  echo "{\"mode\": \"seed\", \"cache_size\": 8000000000, \"generations\": 100, \"thread_limit\": 1, \"min_area_cm\": 0.0, \"normal_grid_path\": \"$NG\"$A}" > p.json
  [ -d g100 ] || { vc_grow_seg_from_seed -v "$PRED" -t $W/tmp -p p.json --seed $x $y $zz --rng-seed 490 > g100.log 2>&1; mv "$(ls -d $W/tmp/*/ | head -1)" $W/g100 2>/dev/null; }
  echo "$z $arm $n : $(grep -m1 -o 'seed_anchor[^,]*' g100.log | head -1) $(tail -n 1 g100.log | cut -c1-80)" >> $OUT; }
export -f run1; export B PRED NG OUT
xargs -P 4 -L 1 bash -c 'run1 "$@"' _ < $B/tasks.txt
for z in z65a90 z50 z80a90; do for arm in base anchor; do
  $PY $C/seedcoh.py score $B/$z/$arm/{A,B,C,D}/g100 --json $B/${z}_$arm.json > /dev/null 2>&1
  echo "SCORE $z $arm $($PY -c "import json;d=json.load(open('$B/${z}_$arm.json'));print(d['pairs_starting_on_same_sheet'], d['mean'], [round(p['rings'].get('0-100',-1),2) for p in d['pairs']])")" >> $OUT
done; done
echo ANCHOR100_DONE >> $OUT
