#!/bin/bash
#SBATCH -J gwcamp
#SBATCH -p mne4
#SBATCH -c 10
#SBATCH --mem=96G
#SBATCH -t 48:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/gwcamp_%A_%a.out
#SBATCH --array=1-65%5
# GENOME-WIDE benchmark campaign. One array task per frozen (sample, condition)
# cell -- the 65 distinct cells where espresso or isoquant has an assembly.
#
# --threads 8 (the runner's own default is 2). The interval loop is ~85% of a
# genome-wide cell and scales with workers; on a real chr1/19/21 cell 8 workers
# reproduced the 2-worker production output byte for byte (twice, so completion
# order is proven not to matter) and ran 2.4x faster. Preconditions Codex named,
# each verified: aggregation is canonicalised because
# post_selection_refit_effective is true in every campaign manifest; --no-gpu
# forces gpu_workers=0 so every worker is CPU; and RAM is ~11 GB parent + ~9 GB
# per worker, hence --mem=96G. Cells finished at 2 workers and cells run at 8
# have identical outputs and may be pooled; RUNTIMES must be stratified.
#
# manifest.tsv is ordered ROUND-ROBIN across samples, not grouped by sample.
# The per-sample evidence prepass is a one-time cost cached to the work dir, but
# grouping put five cells of the SAME sample in the concurrency window, so each
# recomputed that sample's prepass while competing with its own siblings for the
# node. Round-robin gives the window five different samples: each prepass runs
# once, and every later cell of that sample hits the warm cache.
#
# Why genome-wide: every competitor on this grid assembled from the full BAM,
# and only PyFIN ran on a chr1/19/21 subset. Scoring their genome-wide GTFs
# against a three-chromosome truth could not remove that asymmetry, only hide
# it -- their thresholds were computed over the whole genome either way.
#
# -t 48:00:00. The first campaign batch showed per-interval rate scales with
# LIBRARY size (35 iv/min on a 489 MB BAM, 16 iv/min on 920 MB, same node,
# same code), and the two samples the ~5 h estimate was calibrated on are the
# two SMALLEST of the fifteen. Extrapolated, the largest (1970 MB) is ~29 h per
# cell and seven samples exceed 16 h. Two probes of mine have already died to
# a wall limit set from an estimate; this one is set from the largest library.
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
# Completion markers are keyed by ARRAY INDEX, so the manifest must never
# change once a marker exists -- a reordered row would let a marker skip the
# wrong cell. Its checksum is pinned here; the file is also chmod a-w.
MANIFEST_MD5=f58a072bc4bb97223194f4fc1de788c3
[ "$(md5sum $R/manifest.tsv | cut -c1-32)" = "$MANIFEST_MD5" ] || { echo "FATAL manifest.tsv changed (md5 != $MANIFEST_MD5)"; exit 2; }
# Every cell executes the FROZEN source snapshot, not the live tree: queued
# tasks import whatever is on disk when they start, so an edit to fin/ during
# the campaign would otherwise split it across two versions of the code.
# The snapshot's _source_sha256 equals the one tasks 1-4 recorded
# (eedd9b2d5e2ed9e8) and the directory is read-only.
export PYFIN_REPO=/SSD/logan/dev/pyfin/experiments/prod_validation/_snap/gw_campaign
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
  --profile real-drna --extra="--m2-support-gate-attest-reads 1 --threads 8" \
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
