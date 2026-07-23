"""Evidence layer: per-interval read evidence computed ONCE and reused.

The m2_em pipeline used to re-open the BAM several times per interval to rebuild the same
observed-junction map (cluster-recheck GTF guard, Lever-2 per-junction support gate,
junction-dominance pre-gate, guided junction-support gate). This module centralizes that
computation so each interval fetches its evidence once.

Scope (M1, behavior-preserving): only ``observed_junctions`` is consolidated here — it is the
one map that was recomputed identically by multiple gates. Read genomic spans and polyA remain
sourced by their current callers (they use slightly different read sets; unifying them is a
separate, explicitly-validated step). See PIPELINE reorganization plan, contract 5.
"""
from __future__ import annotations

import logging
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

from fin.candidates.intron_chains import extract_intron_chain
from fin.io.interval_manager import GenomicInterval

logger = logging.getLogger(__name__)

# strand -> Counter{(donor, acceptor): n_reads}
ObservedJunctions = Dict[str, "Counter"]


def compute_observed_junctions(
    bam_path: str, interval: GenomicInterval
) -> Optional[ObservedJunctions]:
    """Strand-keyed ``{strand: Counter{(donor, acceptor): n_reads}}`` of intron junctions
    directly observed in the interval's primary-read CIGARs.

    Skips secondary/supplementary/unmapped reads. pysam region strings are 1-based;
    ``interval.start`` is 0-based, so use ``max(start+1, 1)`` (``region_string`` would yield
    ``chr:0-end`` for contig-start intervals, which pysam rejects). Returns None on empty or
    unreadable BAM so every support gate self-disables (fail-open) rather than treating "no
    evidence" as "zero support".

    CAVEAT (unchanged from the original ``_observed_junctions``): ``get_reads_in_region``
    swallows a fetch exception and returns whatever it accumulated first, so a mid-iteration
    failure yields an INCOMPLETE non-empty map (support counts can under-count). Inherent to the
    shared source, preserved here verbatim.
    """
    if not bam_path or not Path(bam_path).exists():
        return None
    from fin.io.io_bam import BamReader

    observed: dict = defaultdict(Counter)
    region = f"{interval.chrom}:{max(interval.start + 1, 1)}-{interval.end}"
    try:
        with BamReader(bam_path) as bam:
            for rd in bam.get_reads_in_region(region):
                if (rd.get("is_secondary") or rd.get("is_supplementary")
                        or not rd.get("is_mapped", True)):
                    continue
                ct = rd.get("cigartuples")
                if not ct:
                    continue
                strand = "-" if rd.get("is_reverse") else "+"
                ic = extract_intron_chain(ct, rd["reference_start"])
                for intr in ic.introns:
                    observed[strand][intr] += 1
    except Exception as exc:  # unreadable/invalid BAM -> recall-safe no-op
        logger.warning(
            "observed-junction build failed for %s (%s); support gates self-disable "
            "for this interval", interval.region_string, exc,
        )
        return None
    return observed if observed else None


@dataclass
class IntervalBundle:
    """Per-interval evidence computed once and shared by every consumer in an interval.

    M1 carries ``observed_junctions`` only; grows to hold read spans / polyA in later steps.
    ``observed_junctions`` is lazily computed on first access and memoized."""

    interval: GenomicInterval
    bam_path: str
    _observed_junctions: Optional[ObservedJunctions] = None
    _observed_built: bool = False

    @property
    def observed_junctions(self) -> Optional[ObservedJunctions]:
        if not self._observed_built:
            self._observed_junctions = compute_observed_junctions(
                self.bam_path, self.interval)
            self._observed_built = True
        return self._observed_junctions
