#!/bin/bash
#SBATCH -J gwcamp
#SBATCH -p mne4
#SBATCH -c 10
#SBATCH --mem=150G
#SBATCH -t 16:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/gwcamp_%A_%a.out
#SBATCH --array=1-65%5
# GENOME-WIDE benchmark campaign. One array task per frozen (sample, condition)
# cell -- the 65 distinct cells where espresso or isoquant has an assembly.
#
# Why genome-wide: every competitor on this grid assembled from the full BAM,
# and only PyFIN ran on a chr1/19/21 subset. Scoring their genome-wide GTFs
# against a three-chromosome truth could not remove that asymmetry, only hide
# it -- their thresholds were computed over the whole genome either way.
#
# -t 16:00:00 deliberately: a cell measured 14.6 h BEFORE the endpoint-EM hoist
# (12fdf4d) and should now be ~5 h. An earlier probe of mine died at 55% because
# I set the limit from an estimate I had not measured.
set -uo pipefail
cd /SSD/logan/dev/pyfin
R=/SSD/logan/dev/pyfin/experiments/prod_validation/_gwcampaign
DONE=$R/_done
mkdir -p $DONE

# An empty SAMPLE or COND is NOT a harmless no-op: the runner reads empty
# --samples/--conds as "every sample, every condition", so one bad array index
# would launch the whole grid and its cell cleanup would delete other tasks'
# output. Refuse the row instead.
[ -r "$R/manifest.tsv" ] || { echo "FATAL manifest unreadable"; exit 2; }
CELL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $R/manifest.tsv)
SAMPLE=$(printf '%s' "$CELL" | cut -f1)
COND=$(printf '%s' "$CELL" | cut -f2)
if [ -z "$SAMPLE" ] || [ -z "$COND" ] || [ "$SAMPLE" = "$COND" ]; then
  echo "FATAL array index ${SLURM_ARRAY_TASK_ID} has no valid manifest row: '$CELL'"; exit 2
fi
MARK=$DONE/${SLURM_ARRAY_TASK_ID}.ok
echo "CELL $SLURM_ARRAY_TASK_ID/65: $SAMPLE / $COND  START $(date +%T)"

# Skip only on a marker written AFTER scoring succeeded. A nonempty
# assembly.gtf proves the assembler ran, not that the cell was scored -- a task
# interrupted between the two would otherwise be skipped forever.
# -f, not -s: markers written by an earlier revision of this script are empty,
# and a marker's VALUE is provenance, not the completion signal itself.
if [ -f "$MARK" ]; then echo "SKIP already scored: $(cat $MARK 2>/dev/null)"; exit 0; fi

TSV=$R/gw_${SLURM_ARRAY_TASK_ID}.tsv
python3 experiments/prod_validation/sirv4_grid_runner.py --dataset gencode --workers 1 \
  --profile real-drna --extra="--m2-support-gate-attest-reads 1" \
  --samples "$SAMPLE" --conds "$COND" \
  --tag _gw --out $R/gw --tsv $TSV > $R/gw_${SLURM_ARRAY_TASK_ID}.log 2>&1
RC=$?
tail -3 $R/gw_${SLURM_ARRAY_TASK_ID}.log

# The runner exits 0 even when a cell's metrics row records pyfin_rc.../
# gffcompare_rc..., so the row's own status is the completion test.
STATUS=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) if($i=="status") c=i; next} {print $c}' "$TSV" 2>/dev/null | head -1)
echo "CELL $SLURM_ARRAY_TASK_ID rc=$RC status=${STATUS:-none} $(date +%T)"
if [ "$RC" -ne 0 ] || [ -z "${STATUS:-}" ] || [ "${STATUS#ok}" = "$STATUS" ]; then
  echo "CELL $SLURM_ARRAY_TASK_ID FAILED"; exit 1
fi
printf '%s\t%s\tstatus=%s\tjob=%s\t%s\n' "$SAMPLE" "$COND" "$STATUS" \
  "${SLURM_ARRAY_JOB_ID:-na}_${SLURM_ARRAY_TASK_ID}" "$(date -Is)" > $MARK
echo "CELL $SLURM_ARRAY_TASK_ID DONE $(date +%T)"
