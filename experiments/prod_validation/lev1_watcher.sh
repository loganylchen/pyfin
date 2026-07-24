#!/usr/bin/env bash
set -uo pipefail
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation
PY=/SSD/logan/miniforge3/envs/nanofusion/bin/python
# 1) heya8: wait for ABLATION_DONE, then score
while ! grep -q ABLATION_DONE "$ROOT/lev1_ablation_heya8/driver.log" 2>/dev/null; do sleep 60; done
$PY "$ROOT/lev1_ablation_heya8/score_ablation.py" > "$ROOT/lev1_ablation_heya8/SUMMARY.txt" 2>&1
echo "HEYA8_SCORED $(date)" >> "$ROOT/lev1_watcher.log"
# 2) gencode: wait for GENCODE_ON_DONE, then score
while ! grep -q GENCODE_ON_DONE "$ROOT/lev1_ablation_gencode/driver.log" 2>/dev/null; do sleep 60; done
$PY "$ROOT/lev1_ablation_gencode/score_on_vs_off.py" > "$ROOT/lev1_ablation_gencode/SUMMARY.txt" 2>&1
echo "GENCODE_SCORED $(date)" >> "$ROOT/lev1_watcher.log"
echo "ALL_SCORED $(date)" >> "$ROOT/lev1_watcher.log"
