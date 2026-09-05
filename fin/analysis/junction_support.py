"""Always-on record of the interval junction registry: ``junction_support.tsv``.

Discovery collects, for every interval, each observed junction's raw support
and its indel-clean support (no I/D within 10 bp of either splice boundary,
MAPQ >= 1 primary reads). That evidence is recorded on every run so it can be
inspected and used later; nothing here decides anything. Worker shards are
merged deterministically into one sorted table at finalization.
"""
from __future__ import annotations

import shutil
from pathlib import Path
from typing import Dict, Optional, Tuple

TSV_FIELDS = ("chrom", "strand", "interval", "donor", "acceptor", "raw_reads", "indel_clean_reads")
SHARD_DIR = "_jsupport_shards"
OUT_NAME = "junction_support.tsv"


def write_support_shard(work_dir, interval, support: Dict[Tuple[int, int], Tuple[int, int]]) -> None:
    """Write one interval's ``{(donor, acceptor): (raw, clean)}`` as a shard."""
    if not support:
        return
    d = Path(work_dir) / SHARD_DIR
    d.mkdir(parents=True, exist_ok=True)
    name = f"{interval.chrom}_{interval.start}_{interval.end}_{interval.strand or 'u'}.tsv"
    tmp = d / (name + ".tmp")
    with tmp.open("w") as fh:
        for (donor, acceptor), (raw, clean) in sorted(support.items()):
            fh.write("\t".join(map(str, (
                interval.chrom, interval.strand or "", interval.region_string,
                donor, acceptor, raw, clean))) + "\n")
    tmp.replace(d / name)


def merge_support_shards(work_dir) -> Optional[str]:
    """Merge shards into a sorted ``junction_support.tsv``; remove the shard dir."""
    wd = Path(work_dir)
    shard_dir = wd / SHARD_DIR
    rows = []
    if shard_dir.is_dir():
        for shard in sorted(shard_dir.glob("*.tsv")):
            with shard.open() as fh:
                rows.extend(line.rstrip("\n") for line in fh if line.strip())
    if not rows:
        # Write the header even with no junctions: a missing file cannot be told
        # apart from "the recording step never ran", and this sidecar is
        # advertised as always present.
        wd.mkdir(parents=True, exist_ok=True)
        out = wd / OUT_NAME
        tmp = out.with_suffix(".tsv.tmp")
        with tmp.open("w") as fh:
            fh.write("\t".join(TSV_FIELDS) + "\n")
        tmp.replace(out)
        if shard_dir.is_dir():
            shutil.rmtree(shard_dir)
        return str(out)
    rows.sort(key=lambda l: (l.split("\t")[0], l.split("\t")[1], int(l.split("\t")[3]),
                             int(l.split("\t")[4]), l.split("\t")[2]))
    out = wd / OUT_NAME
    tmp = out.with_suffix(".tsv.tmp")
    with tmp.open("w") as fh:
        fh.write("\t".join(TSV_FIELDS) + "\n")
        for line in rows:
            fh.write(line + "\n")
    tmp.replace(out)
    shutil.rmtree(shard_dir)
    return str(out)
