#!/bin/bash
T=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid/expression_contract/_truth/SGNex_H9_directRNA_replicate2_run2.gtf
G=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
IMG=quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0
for arm in m2base m2off m2scram; do
  out=$G/$arm.expr.tsv
  echo -e "cond\ttool\tn_query\ttx_sn\ttx_pr\ttx_f1" > "$out"
  for d in $G/$arm/SGNex_H9_directRNA_replicate2_run2/*/; do
    c=$(basename "$d"); gtf="$d/out/assembly.gtf"
    [ -s "$gtf" ] || continue
    w="$d/gc_expr"; mkdir -p "$w"
    docker run --rm -v /SSD:/SSD $IMG gffcompare -r "$T" -o "$w/cmp" "$gtf" >/dev/null 2>&1
    s="$w/cmp.stats"
    sn=$(grep -m1 "Transcript level" "$s" | awk '{print $3}')
    pr=$(grep -m1 "Transcript level" "$s" | awk '{print $5}')
    nq=$(grep -m1 "Query mRNAs" "$s" | tr -dc '0-9 ' | awk '{print $1}')
    f1=$(python3 -c "sn=$sn;pr=$pr;print(0 if sn+pr==0 else 2*sn*pr/(sn+pr))")
    echo -e "$c\tpyfin_$arm\t$nq\t$sn\t$pr\t$f1" >> "$out"
  done
  echo "$arm: $(( $(wc -l < "$out") - 1 )) cells"
done
