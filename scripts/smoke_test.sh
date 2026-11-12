#!/bin/sh
# Check an installation end to end: run PyFIN on the bundled signal fixture
# (91 reads of one SIRV locus) and confirm that models came out and that the
# raw signal was actually scored.
#
#   scripts/smoke_test.sh                 # uses the fixture next to this script
#   scripts/smoke_test.sh /data /tmp/out  # fixture directory, output directory
#
# Takes a few seconds on one CPU core. Exit status 0 means the installation works.
set -eu
DATA=${1:-"$(dirname "$0")/../tests/data/sirv1_mini"}
OUT=${2:-"$(mktemp -d)/pyfin_smoke"}

pyfin --version
pyfin --bam "$DATA/reads.bam" --genome "$DATA/reference.fa" --signal "$DATA/reads.blow5" \
      --gtf "$DATA/annotation.gtf" --output-dir "$OUT" --no-gpu

python3 - "$OUT/run_summary.json" <<'PY'
import json
import sys

summary = json.load(open(sys.argv[1]))
models, signal = summary["models"], summary["signal"]
if not models["total"]:
    sys.exit("FAILED: no transcript model was assembled")
if not signal["scored_contrasts"]:
    sys.exit("FAILED: the run finished but the signal was never scored (is krill installed?)")
print(f"OK: {models['total']} models ({models['annotated']} annotated, {models['novel']} novel); "
      f"{signal['tied_reads']} reads in an alignment-score tie, "
      f"{signal['scored_contrasts']} signal contrasts scored with the {signal['pore']} model")
PY
