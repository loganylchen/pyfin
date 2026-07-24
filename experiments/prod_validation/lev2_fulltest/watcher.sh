#!/usr/bin/env bash
ROOT=/SSD/logan/dev/pyfin/experiments/prod_validation/lev2_fulltest
PY=/SSD/logan/miniforge3/envs/nanofusion/bin/python
last=0
while true; do
  if grep -q FULLTEST_DONE "$ROOT/driver.log" 2>/dev/null; then break; fi
  n=$(grep -c DATASET_DONE "$ROOT/driver.log" 2>/dev/null||echo 0)
  if [ "$n" -gt "$last" ]; then
    $PY "$ROOT/score_fulltest.py" > "$ROOT/SUMMARY.txt" 2>&1
    echo "SCORED after $n datasets $(date)" >> "$ROOT/watcher.log"; last=$n
  fi
  sleep 120
done
$PY "$ROOT/score_fulltest.py" > "$ROOT/SUMMARY.txt" 2>&1
echo "FINAL_SCORED $(date)" >> "$ROOT/watcher.log"
