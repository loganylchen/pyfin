#!/usr/bin/env bash
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation
PY=/SSD/logan/miniforge3/envs/nanofusion/bin/python
while ! grep -q LEV2_DONE "$ROOT/lev2_ablation/driver.log" 2>/dev/null; do sleep 60; done
$PY "$ROOT/lev2_ablation/score_lev2.py" > "$ROOT/lev2_ablation/SUMMARY.txt" 2>&1
echo "LEV2_SCORED $(date)" >> "$ROOT/lev2_watcher.log"
