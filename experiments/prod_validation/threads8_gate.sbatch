#!/bin/bash
#SBATCH -J t8gate2
#SBATCH -p mne4
#SBATCH -c 12
#SBATCH --mem=110G
#SBATCH -t 4:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/t8gate2_%j.out
# The layer-3 gate Codex asked for and I skipped: a REAL chr1/19/21 cell on the
# current tree, byte-compared against the pre-patch production output (attK1,
# 2026-09-11, K=1 on the att_land snapshot; the only assembly-path commits
# since are 07e1776 and 12fdf4d). Two cells: full (annotated) and p00 (de novo,
# no --gtf, so the reader-reuse path is exercised with no reader at all).
set -uo pipefail
FAIL=0
cd /SSD/logan/dev/pyfin
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset gencode_sub --workers 2 \
  --profile real-drna --extra="--m2-support-gate-attest-reads 1 --threads 8" \
  --samples SGNex_H9_directRNA_replicate2_run2 --conds full,p00 \
  --tag _t8gate2 --out $R/t8gate2 --tsv $R/t8gate2.tsv 2>&1 | tail -3
RC=${PIPESTATUS[0]}
# A failed runner leaves the PREVIOUS run's outputs in place; comparing those
# would pass the gate without having run the current code. Fail here first.
if [ "$RC" -ne 0 ]; then echo "GATE FAILED: runner rc=$RC $(date +%T)"; exit 1; fi
for c in full p00; do
  S=$(awk -F'\t' -v c=$c 'NR==1{for(i=1;i<=NF;i++) if($i=="status") k=i; next} $2==c{print $k}' $R/t8gate2.tsv 2>/dev/null | head -1)
  [ "${S#ok}" != "$S" ] || { echo "GATE FAILED: cell $c status=${S:-none} $(date +%T)"; exit 1; }
done
echo "=== BYTE IDENTITY vs attK1 (pre-patch production) ==="
for c in full p00; do for f in assembly.gtf scores.tsv abundance_refit.json; do
  A=$R/attK1/SGNex_H9_directRNA_replicate2_run2/$c/out/$f; B=$R/t8gate2/SGNex_H9_directRNA_replicate2_run2/$c/out/$f
  if cmp -s "$A" "$B"; then echo "IDENTICAL $c/$f"; else echo "DIFFERS   $c/$f  ($(wc -l < $A 2>/dev/null) vs $(wc -l < $B 2>/dev/null) lines)"; FAIL=1; fi
done; done
# SLURM COMPLETED is not a pass: the gate passes only if every comparison did.
if [ "$FAIL" -ne 0 ]; then echo "GATE FAILED $(date +%T)"; exit 1; fi
echo "GATE PASSED $(date +%T)"
