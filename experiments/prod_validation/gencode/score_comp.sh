#!/usr/bin/env bash
PV=/SSD/logan/dev/pyfin/experiments/prod_validation
GC=/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results/gencode_full_sweep
TRUTH=$(readlink -f $PV/gencode/_ref/full/annotation.gtf)
IMG=quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0
score(){ local q="$1" lab="$2" s="$3"; [ -s "$q" ] || return
  local sd=$PV/gencode/$s/full/score_$lab; mkdir -p $sd
  cp "$q" "$sd/q.gtf"   # local copy so tmap is writable
  docker run --rm -u $(id -u):$(id -g) -v /SSD:/SSD -v /autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD:ro \
    -w "$sd" $IMG gffcompare -r "$TRUTH" -o gc "$sd/q.gtf" >/dev/null 2>&1
  local line=$(grep -A1 "Transcript level" $sd/gc.stats 2>/dev/null|tail -1)
  local sn=$(echo "$line"|grep -oE "[0-9]+\.[0-9]+"|head -1); local pr=$(echo "$line"|grep -oE "[0-9]+\.[0-9]+"|sed -n 2p)
  local nout=$(grep -c $'\ttranscript\t' "$q"); local f1=$(python3 -c "s=$sn or 0;p=$pr or 0;print(f'{2*s*p/(s+p):.1f}' if s+p else 'NA')" 2>/dev/null)
  echo -e "$s\t$lab\t$nout\t$sn\t$pr\t$f1" >> $PV/gencode/comparison.tsv
  echo "  $s/$lab: Sn=$sn Pr=$pr F1=$f1"
}
for s in SGNex_H9_directRNA_replicate2_run2 SGNex_H9_directRNA_replicate3_run1 SGNex_H9_directRNA_replicate3_run2 SGNex_H9_directRNA_replicate4_run2 SGNex_HEYA8_directRNA_replicate1_run2; do
  score "$GC/$s/full/assembly/isoquant.gtf" isoquant "$s"
  score "$GC/$s/full/assembly/espresso.gtf" espresso "$s"
done
