#!/usr/bin/env bash
# gffcmp.sh <out_prefix> <ref.gtf> <query.gtf>
# Runs gffcompare (biocontainer) producing <out_prefix>.* (tmap/refmap/stats).
# All paths are made absolute and the common parent is mounted into the container.
set -euo pipefail
PREFIX="$(readlink -f "$1")"; REF="$(readlink -f "$2")"; QRY="$(readlink -f "$3")"
IMG="quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"
MNT=/SSD/logan/dev/pyfin
docker run --rm -u "$(id -u):$(id -g)" -v "$MNT":"$MNT" -w "$(dirname "$PREFIX")" \
  "$IMG" gffcompare -r "$REF" -o "$PREFIX" "$QRY"
echo "--- $(basename "$PREFIX").stats (Transcript level) ---"
grep -A1 "Transcript level" "${PREFIX}.stats" || true
