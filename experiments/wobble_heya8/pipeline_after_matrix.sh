#!/usr/bin/env bash
# Wait for the matrix run to finish, then score every cell and build per-sample tables.
# Idempotent — safe to re-run.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

LOG=run_logs/_pipeline.log
exec > >(tee -a "$LOG") 2>&1
echo "[$(date)] pipeline_after_matrix start (waiting on matrix)"

# Poll matrix log every 60s until 'MATRIX DONE' appears
while ! grep -q "MATRIX DONE" run_logs/_matrix.log 2>/dev/null; do
  sleep 60
done
echo "[$(date)] matrix finished. tail:"
tail -3 run_logs/_matrix.log

echo "[$(date)] starting scoring"
bash score_matrix.sh
echo "[$(date)] scoring done"

echo "[$(date)] building tables"
python3 build_tables.py
echo "[$(date)] tables built. files:"
ls -la tables/ 2>/dev/null

echo "[$(date)] scanning cassette-pair counts per cell"
python3 cassette_scan.py
echo "[$(date)] cassette scan done"

touch run_logs/_pipeline_done.flag
echo "[$(date)] PIPELINE COMPLETE"
