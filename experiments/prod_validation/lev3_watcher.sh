#!/usr/bin/env bash
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation
PY=/SSD/logan/miniforge3/envs/nanofusion/bin/python
while ! grep -q LEV3_DONE "$ROOT/lev3_ablation/driver.log" 2>/dev/null; do sleep 60; done
$PY "$ROOT/lev3_ablation/score_lev3.py" > "$ROOT/lev3_ablation/SUMMARY.txt" 2>&1
echo "LEV3_SCORED $(date)" >> "$ROOT/lev3_watcher.log"
