#!/usr/bin/env bash
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/lev2_fulltest
PY=/SSD/logan/miniforge3/envs/nanofusion/bin/python
while ! grep -q GENCODE_NOW_DONE "$ROOT/gencode_now.log" 2>/dev/null; do sleep 120; done
$PY "$ROOT/score_fulltest.py" > "$ROOT/SUMMARY.txt" 2>&1
echo "FINAL_SCORED $(date)" >> "$ROOT/watcher.log"
