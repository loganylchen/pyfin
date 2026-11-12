#!/bin/bash
#SBATCH -J h2h
#SBATCH -p mne4
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 4:00:00
#SBATCH -o /SSD/logan/dev/pyfin/logs/slurm/h2h_%j.out
# Aggregate the genome-wide campaign under the EXPRESSION TRUTH CONTRACT
# (truth = MASTER GENCODE transcripts with NanoCount TPM > 0 for that sample),
# which is the contract every earlier comparison in this project used. The
# per-cell gw_<N>.tsv numbers score against the FULL annotation instead and are
# not comparable to them.
#
# Cells are cached by (sample, cond, tool) under h2h/cells; a rerun scores only
# what is new. No --force here: forcing would rescore all 64 cells for nothing.
set -uo pipefail
cd /SSD/logan/dev/pyfin
python3 experiments/prod_validation/gw_headtohead.py
RC=$?
echo "H2H rc=$RC $(date +%T)"
# A failed pass must fail the job: scores.tsv is then still the previous
# table (the scorer replaces it only when it finishes), and a job that
# exited 0 would present that old table as this pass's result.
[ "$RC" -eq 0 ] || { echo "H2H FAILED; scores.tsv was not replaced"; exit "$RC"; }
# The scorer removes its own scratch copies (query.gtf, cmp.annotated.gtf)
# before it releases the scoring lock; nothing is deleted from here, where
# another pass could already be working in the same tree.
echo "H2H done, /SSD free $(df -BG --output=avail /SSD | tail -1 | tr -dc '0-9')G"
