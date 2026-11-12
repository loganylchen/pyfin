#!/usr/bin/env python3
"""Lossless disk compaction of FINISHED genome-wide campaign cells.

/SSD is a shared volume and the campaign's per-cell footprint (~510 MB) would
exhaust it before part B finishes. Everything here is lossless and touches only
cells whose completion marker exists (part A: manifest.tsv + _done/, part B:
manifest_b.tsv + _done_b/), so a running or queued cell is never read or
modified. Idempotent: re-running does nothing new to cells already compacted.

  1. _sample_evidence_<hash>.json is byte-identical across every condition of
     one sample (same hash name, same md5). Byte-identical copies are replaced
     by hard links to one canonical copy. Every path still reads the same bytes;
     a later rmtree of one cell removes only that link.
  2. gw/<s>/<c>/run.log -> run.log.gz. The runner writes run.log only while the
     cell runs; nothing reads a finished cell's log.
  3. h2h/cells/<s>/<c>/<tool>/cmp.tracking and cmp.loci -> .gz. The scorer
     caches on cmp.stats + cache_key.json and never reads these two; a forced
     or stale rescore deletes cmp* (including the .gz) and regenerates them.

It must not run while the scorer does. The two share a lock (gw_cells): the
scorer holds it for a whole pass, and this script refuses to start when it is
held rather than wait, and holds it while it works.
"""
import argparse
import collections
import gzip
import hashlib
import os
import shutil
import sys
from pathlib import Path

sys.path.insert(0, "/SSD/logan/dev/pyfin/experiments/prod_validation")
from gw_cells import finished_cells, scoring_lock  # noqa: E402

CAMP = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_gwcampaign")


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_inplace(path):
    gz = path.with_name(path.name + ".gz")
    tmp = gz.with_name(gz.name + ".partial")
    with open(path, "rb") as src, gzip.open(tmp, "wb", compresslevel=6) as dst:
        shutil.copyfileobj(src, dst, 1 << 22)
    os.replace(tmp, gz)
    os.unlink(path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    try:
        with scoring_lock(CAMP, wait=False):
            compact(a.dry_run)
    except BlockingIOError:
        raise SystemExit("the h2h tree is locked (a scoring pass or another compaction "
                         "is running); refusing to compact") from None


def compact(dry_run):
    cells = {(s, c) for _, s, c in finished_cells(CAMP)}
    saved = collections.Counter()

    # 1. hard-link byte-identical sample evidence, per sample
    by_key = collections.defaultdict(list)
    for s, c in cells:
        for f in (CAMP / "gw" / s / c / "out").glob("_sample_evidence_*.json"):
            by_key[(s, f.name, f.stat().st_size)].append(f)
    for (s, name, size), files in by_key.items():
        files.sort()
        canon = files[0]
        canon_md5 = md5(canon)
        for f in files[1:]:
            if os.path.samefile(f, canon):
                continue
            if md5(f) != canon_md5:
                continue
            if not dry_run:
                tmp = f.with_name(f.name + ".linktmp")
                if tmp.exists():
                    tmp.unlink()
                os.link(canon, tmp)
                os.replace(tmp, f)
            saved["evidence_hardlink"] += size

    # 2. gzip finished cells' run.log
    for s, c in cells:
        log = CAMP / "gw" / s / c / "run.log"
        if log.exists():
            saved["run_log_gzip_input"] += log.stat().st_size
            if not dry_run:
                gzip_inplace(log)

    # 3. gzip h2h cmp.tracking / cmp.loci of finished cells
    for s, c in cells:
        for tooldir in (CAMP / "h2h" / "cells" / s / c).glob("*"):
            for name in ("cmp.tracking", "cmp.loci"):
                p = tooldir / name
                if p.exists():
                    saved["h2h_gzip_input"] += p.stat().st_size
                    if not dry_run:
                        gzip_inplace(p)

    for k, v in saved.items():
        print(f"{k:22s} {v / 1073741824:7.2f} GB")
    print("cells considered:", len(cells), "(dry run)" if dry_run else "")


if __name__ == "__main__":
    main()
