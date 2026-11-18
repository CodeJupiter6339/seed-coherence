#!/bin/bash
# voxel : graine figee (seed_anchor_xyz) contre base sur 0841, jugee contre le segment humain w00 (g100, memes graines que seedbench2)
ulimit -Sn 750000
PY=~/vesuvius/gg/bin/python; C=~/vesuvius/code; B=~/vesuvius/anchor_w00; mkdir -p $B; OUT=$B/run.out
W=~/vesuvius/data/0841/w00
until grep -q SEEDBENCH2_DONE ~/vesuvius/seedbench2/run.out && grep -q ANCHOR100_DONE ~/vesuvius/anchor100/run.out; do sleep 120; done
echo "debut $(date +%T)" >> $OUT
grep '^w00' ~/vesuvius/seedbench2/tasks.txt > $B/tasks.txt
run1() { n=$1; x=$2; y=$3; z=$4; P=$5; NG=$6; Wd=$B/$n; mkdir -p $Wd; cd $Wd
  echo "{\"mode\": \"seed\", \"cache_size\": 15000000000, \"generations\": 100, \"thread_limit\": 5, \"min_area_cm\": 0.0, \"normal_grid_path\": \"$NG\", \"seed_anchor_xyz\": [$x, $y, $z]}" > p.json
  [ -d g100 ] || { vc_grow_seg_from_seed -v "$P" -t $Wd/tmp -p p.json --seed $x $y $z --rng-seed 490 > g100.log 2>&1; mv "$(ls -d $Wd/tmp/*/ | head -1)" $Wd/g100 2>/dev/null; }
  echo "$n fini $(date +%T) $(grep -m1 -o 'seed_anchor event=[^ ]*' g100.log)" >> $OUT; }
export -f run1; export B OUT
xargs -P 4 -L 1 bash -c 'run1 "$@"' _ < $B/tasks.txt
for n in $(cut -d/ -f1 $B/tasks.txt | sort -u); do
  $PY $C/truth.py $W $B/${n}_anchor_truth.json $B/$n/?/g100 >> $OUT 2>/dev/null
  $PY $C/seedcoh.py score $B/$n/?/g100 --json $B/${n}_anchor.json > /dev/null 2>&1
done
echo ANCHOR_W00_DONE >> $OUT
