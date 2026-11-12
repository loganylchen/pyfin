#!/usr/bin/env python3
"""The genome-wide campaign's finished cells, and the lock that guards scoring.

The campaign is two pinned manifests that index one output tree:

  part A   manifest.tsv     markers in _done/     idx = line number
  part B   manifest_b.tsv   markers in _done_b/   idx = line number + 1000

A marker `<N>.ok` means line N of its manifest finished and was scored by the
runner. Line N is the physical line, as `sed -n Np` in the sbatch scripts
reads it. Everything that turns markers back into (sample, condition) goes
through `finished_cells`, so the scorer and the compactor cannot disagree
about which cell a marker names, and none of them can proceed on a manifest
that is missing, edited or reordered.
"""
import contextlib
import fcntl
import hashlib
from pathlib import Path

CAMP = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_gwcampaign")
# (manifest, marker directory, idx offset, md5 the sbatch scripts pin)
PARTS = (
    ("manifest.tsv", "_done", 0, "f58a072bc4bb97223194f4fc1de788c3"),
    ("manifest_b.tsv", "_done_b", 1000, "f855a81005fff7f06cae161ee6d80a8a"),
)


class ManifestError(RuntimeError):
    pass


def _rows(path, md5):
    try:
        raw = path.read_bytes()
    except OSError as exc:
        raise ManifestError(f"{path}: cannot be read ({exc})") from exc
    got = hashlib.md5(raw).hexdigest()
    if got != md5:
        raise ManifestError(f"{path}: md5 {got}, pinned {md5}; markers no longer mean what they meant")
    rows = []
    for n, line in enumerate(raw.decode().splitlines(), 1):
        fields = line.split("\t")
        if len(fields) != 2 or not all(fields):
            raise ManifestError(f"{path}:{n}: expected 'sample<TAB>condition', got {line!r}")
        rows.append(tuple(fields))
    return rows


def finished_cells(camp=CAMP, parts=PARTS):
    """[(idx, sample, cond)] for every completion marker, sorted by idx.

    Raises ManifestError rather than return a subset: a manifest that is
    absent or changed, a marker outside its manifest, a marker whose name is
    not a line number, or one (sample, cond) reachable from two markers.
    """
    camp = Path(camp)
    cells, origin = [], {}
    for manifest, done_dir, offset, md5 in parts:
        rows = _rows(camp / manifest, md5)
        if len(rows) >= 1000 and offset == 0:
            raise ManifestError(f"{manifest}: {len(rows)} rows collide with the part-B idx offset")
        for marker in sorted((camp / done_dir).glob("*.ok")):
            if not marker.stem.isdigit() or not 1 <= int(marker.stem) <= len(rows):
                raise ManifestError(f"{marker}: not a line number of {manifest} (1..{len(rows)})")
            sample, cond = rows[int(marker.stem) - 1]
            if (sample, cond) in origin:
                raise ManifestError(f"{marker}: ({sample}, {cond}) is also {origin[(sample, cond)]}")
            origin[(sample, cond)] = marker
            cells.append((int(marker.stem) + offset, sample, cond))
    return sorted(cells)


@contextlib.contextmanager
def scoring_lock(camp=CAMP, *, wait):
    """Exclusive hold on the h2h tree. The scorer waits for it; the compactor
    refuses instead (``wait=False`` raises BlockingIOError when it is held)."""
    path = Path(camp) / "h2h" / ".lock"
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fcntl.flock(fh, fcntl.LOCK_EX | (0 if wait else fcntl.LOCK_NB))
        try:
            yield
        finally:
            fcntl.flock(fh, fcntl.LOCK_UN)
