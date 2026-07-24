#!/usr/bin/env bash
DS=/SSD/logan/dev/pyfin/experiments/prod_validation/gencode
TRUTH=$DS/_ref/full/annotation.gtf
IMG=quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0
echo "sample	nout	Tx_Sn	Tx_Pr" > $DS/pyfin_vs_gencode.tsv
for g in $DS/*/full/pyfin.gtf; do
  s=$(basename $(dirname $(dirname $g)))
  sd=$(dirname $g)/scoring; mkdir -p $sd
  docker run --rm -u $(id -u):$(id -g) -v /SSD:/SSD -v /autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD:ro \
    -w "$sd" $IMG gffcompare -r "$(readlink -f $TRUTH)" -o gc "$(readlink -f $g)" >/dev/null 2>&1
  sn=$(grep -A1 "Transcript level" $sd/gc.stats 2>/dev/null | grep -oE "[0-9]+\.[0-9]+" | head -1)
  pr=$(grep -A1 "Transcript level" $sd/gc.stats 2>/dev/null | grep -oE "[0-9]+\.[0-9]+" | sed -n 2p)
  nout=$(grep -c $'\ttranscript\t' "$g")
  echo "$s	$nout	$sn	$pr" >> $DS/pyfin_vs_gencode.tsv
done
cat $DS/pyfin_vs_gencode.tsv
