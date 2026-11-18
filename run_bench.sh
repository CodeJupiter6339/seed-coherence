#!/bin/bash
# Grow 4 seeds around the seed of an existing GrowPatch mesh, to 100 then 140 generations, and score them.
#   run_bench.sh PRED_ZARR_URL NORMAL_GRIDS_DIR SOURCE_MESH OUT_DIR [THREADS]
# Needs vc_grow_seg_from_seed (VC3D, ScrollPrize/villa) on PATH and python3 with numpy scipy tifffile.
set -e
PRED=$1; NG=$2; SRC=$3; OUT=$4; TH=${5:-6}; HERE=$(cd "$(dirname "$0")" && pwd)
mkdir -p "$OUT"; python3 "$HERE/seedcoh.py" seeds "$SRC" > "$OUT/seeds.txt"; cat "$OUT/seeds.txt"
grow() { W=$OUT/$1; mkdir -p $W
  printf '{"mode": "seed", "cache_size": 20000000000, "generations": %s, "thread_limit": %s, "min_area_cm": 0.0, "normal_grid_path": "%s"}\n' 100 $TH "$NG" > $W/p100.json
  sed 's/"generations": 100/"generations": 140/' $W/p100.json > $W/p140.json
  vc_grow_seg_from_seed -v "$PRED" -t $W/tmp -p $W/p100.json --seed $2 $3 $4 --rng-seed 490 > $W/g100.log 2>&1
  mv "$(ls -d $W/tmp/*/ | head -1)" $W/g100
  vc_grow_seg_from_seed -v "$PRED" -t $W/g140 -p $W/p140.json --resume $W/g100 --resume-opt skip --segment-name g140 --rng-seed 490 > $W/g140.log 2>&1; }
while read n x y z; do grow $n $x $y $z & done < "$OUT/seeds.txt"; wait
for g in g100 g140; do python3 "$HERE/seedcoh.py" score $OUT/{A,B,C,D}/$g --json $OUT/score_$g.json; done
