#!/bin/bash
#SBATCH -J gwcampB
#SBATCH -p mne4
#SBATCH -c 10
#SBATCH --mem=144G
#SBATCH -t 48:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/gwcampB_%A_%a.out
#SBATCH --array=1-325%4
# GENOME-WIDE campaign, PART B: the 325 (sample, condition) cells where
# stringtie3 has an assembly and part A (gw_campaign.sbatch, 65 cells chosen
# for espresso/isoquant coverage) did not run PyFIN. Same frozen source
# snapshot, same profile and flags, same output tree, so the two parts pool
# into one 390-cell head-to-head. Part A's rationale for --threads, --mem and
# the marker discipline is in gw_campaign.sbatch and is not repeated here.
#
# What differs from part A, and why:
#   manifest_b.tsv / _done_b/    part A's markers are keyed by index into
#                                manifest.tsv, which is pinned; a second
#                                manifest needs its own marker directory.
#   threads by library size      the two libraries over 1.9 / 1.6 GB (HEYA8
#                                rep3 run1, rep1 run1) OOM-killed at 8 workers
#                                inside 144G and finished at 4; 4 and 8
#                                workers produce byte-identical output.
#   a part-A collision guard     the runner rmtree()s its cell at start, so a
#                                manifest_b row that duplicated a part-A cell
#                                would destroy a finished result. The rows are
#                                disjoint by construction; this refuses anyway.
set -uo pipefail
cd /SSD/logan/dev/pyfin
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_gwcampaign
DONE=$R/_done_b
mkdir -p $DONE

[ -r "$R/manifest_b.tsv" ] || { echo "FATAL manifest_b unreadable"; exit 2; }
MANIFEST_MD5=f855a81005fff7f06cae161ee6d80a8a
[ "$(md5sum $R/manifest_b.tsv | cut -c1-32)" = "$MANIFEST_MD5" ] || { echo "FATAL manifest_b.tsv changed (md5 != $MANIFEST_MD5)"; exit 2; }
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/gw_campaign
[ -d "$PYFIN_REPO/fin" ] || { echo "FATAL frozen snapshot missing at $PYFIN_REPO"; exit 2; }

CELL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $R/manifest_b.tsv)
SAMPLE=$(printf '%s' "$CELL" | cut -f1)
COND=$(printf '%s' "$CELL" | cut -f2)
if [ -z "$SAMPLE" ] || [ -z "$COND" ] || [ "$SAMPLE" = "$COND" ]; then
  echo "FATAL array index ${SLURM_ARRAY_TASK_ID} has no valid manifest_b row: '$CELL'"; exit 2
fi
MARK=$DONE/${SLURM_ARRAY_TASK_ID}.ok
echo "CELL B${SLURM_ARRAY_TASK_ID}/325: $SAMPLE / $COND  START $(date +%T)"
# One attempt per cell at a time. The runner rmtree()s the cell when it
# starts, so a second submission of an element that is still running would
# delete the first attempt's output under it. The lock is taken BEFORE the
# marker is read and held until this shell exits; a second attempt fails
# (exit 4) rather than report success for work it did not do.
mkdir -p $R/_locks
exec 9> "$R/_locks/${SAMPLE}__${COND}.lock"
flock -n 9 || { echo "FATAL another attempt at ($SAMPLE, $COND) holds the cell lock; not starting"; exit 4; }
if [ -f "$MARK" ]; then echo "SKIP already scored: $(cat $MARK 2>/dev/null)"; exit 0; fi

# part-A collision guard: a finished part-A cell has this exact (sample, cond).
# The guard must fail CLOSED: an unreadable part-A manifest would make grep
# find nothing and let the task proceed.
[ -r "$R/manifest.tsv" ] || { echo "FATAL part-A manifest unreadable; cannot prove no collision"; exit 2; }
if grep -qxF "$(printf '%s\t%s' "$SAMPLE" "$COND")" $R/manifest.tsv; then
  echo "FATAL ($SAMPLE, $COND) is a part-A cell; refusing to touch its output"; exit 2
fi

case "$SAMPLE" in
  SGNex_HEYA8_directRNA_replicate3_run1|SGNex_HEYA8_directRNA_replicate1_run1) THREADS=4 ;;
  *) THREADS=8 ;;
esac

FREE_GB=$(df -BG --output=avail /SSD | tail -1 | tr -dc '0-9')
if [ "${FREE_GB:-0}" -lt 20 ]; then
  echo "FATAL only ${FREE_GB}G free on /SSD (shared volume); refusing to start"; exit 3
fi

TSV=$R/gw_b_${SLURM_ARRAY_TASK_ID}.tsv
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset gencode --workers 1 \
  --profile real-drna --extra="--m2-support-gate-attest-reads 1 --threads $THREADS" \
  --samples "$SAMPLE" --conds "$COND" \
  --tag _gw --out $R/gw --tsv $TSV > $R/gw_b_${SLURM_ARRAY_TASK_ID}.log 2>&1
RC=$?
tail -3 $R/gw_b_${SLURM_ARRAY_TASK_ID}.log

STATUS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="status") c=i; next} {print $c}' "$TSV" 2>/dev/null | head -1)
echo "CELL B$SLURM_ARRAY_TASK_ID rc=$RC status=${STATUS:-none} threads=$THREADS $(date +%T)"
case "${STATUS:-}" in
  pyfin_rc-9) echo "CELL B$SLURM_ARRAY_TASK_ID HINT: SIGKILL -- likely OOM; check 'Detected oom-kill' in this log" ;;
esac
if [ "$RC" -ne 0 ] || [ -z "${STATUS:-}" ] || [ "${STATUS#ok}" = "$STATUS" ]; then
  echo "CELL B$SLURM_ARRAY_TASK_ID FAILED"; exit 1
fi
printf '%s\t%s\tstatus=%s\tjob=%s\t%s\n' "$SAMPLE" "$COND" "$STATUS" \
  "${SLURM_ARRAY_JOB_ID:-na}_${SLURM_ARRAY_TASK_ID}" "$(date -Is)" > $MARK
echo "CELL B$SLURM_ARRAY_TASK_ID DONE $(date +%T)"
