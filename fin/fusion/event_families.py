"""Which reported models describe one junction geometry.

Three BAG6>SLC44A4 rows were published from the whole K562 library as
4 / 2 / 2 supporting molecules with three different `alternative_group_id`s
and `support_compatible` 0 on all three. Two of them carry byte-identical
breakpoints; the third sits 6 bp away at the acceptor. Nothing in the output
linked them, so a reader could not tell one event with eight molecules and
three reported structures from three independent weak events.

The existing group (`CAG_`) is built from SHARED MOLECULES. Exact membership
does not in general stop a molecule from fitting several envelopes -- that is
exactly what `compatible_read_ids` records -- but in these three BAG6 rows no
molecule did, so each group degenerated to its own model. Shared molecules are
the right grouping for "which models compete for the same reads"; they are the
wrong one for "which models describe the same junction".

This module adds the second grouping and nothing else. It does not change
which models are emitted, does not touch the support gate, does not move a
molecule between models, and never mutates a model. Run it after
`build_fused_transcript_models` -- both merge passes and the endpoint
recomputation must already have decided which models exist.

What the family does NOT claim
------------------------------
* Not a shared candidate set. The name of the work item said "shared candidate
  set"; this is shared REPORTING. Discovery and sensitivity are unchanged.
* Geometry, not biology. Junctions within tolerance MAY be grouped -- the
  bounded-diameter partition deliberately separates some nearby pairs, and it
  groups without regard to whether they arose from one molecular event. The
  arms' gene assignments are deliberately not part of the key (three of the
  six duplicate breakpoint strings in the reference run are at unannotated
  loci, and a gene-based key would not group them at all, while chr16 holds
  two families 14-15 bp apart that the tolerance keeps separate). The
  defensible sentence is "three reported
  models in one tolerance-defined junction family, supported by eight distinct
  molecules" -- not "one event, three isoforms".
* The diameter bound is on the REPRESENTATIVE breakpoints each model reports.
  A model's retained `coordinate_alternatives` are not membership edges: the
  two BAG6 `...784` models hold an alternative at 31865787, so across every
  represented placement that family spans 9 bp, not 6. Joining on any matching
  alternative would give ambiguity a second route to bridge families.
* Family support is a union of molecules over the family's EMITTED models, via
  their unique or compatible assignments. Molecules supporting a structure that
  failed the gate are absent, and a compatible molecule may support models in
  more than one family, so support is neither exclusive nor additive across
  families. Union the ids; never sum the counts.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from typing import Dict, Iterable, Tuple

# Bumped when the partition rule or the key changes, so a family id computed
# under a different rule can never collide with one computed under this rule.
ALGORITHM_VERSION = "1"


@dataclass(frozen=True)
class EventFamily:
    family_id: str
    model_ids: Tuple[str, ...]
    topology: Tuple[Tuple[str, str, str, str], ...]
    breakpoint_vectors: Tuple[Tuple[int, ...], ...]
    support_read_ids: Tuple[str, ...]

    @property
    def model_count(self) -> int:
        return len(self.model_ids)

    @property
    def support(self) -> int:
        return len(self.support_read_ids)


def _topology(model) -> Tuple[Tuple[str, str, str, str], ...]:
    """Ordered contig/strand identity of every junction, both sides.

    `model.breakpoints` entries are (chromL, exitL, strandL, chromR, entryR,
    strandR). Only the identity half is topology; coordinates are compared
    separately and with a tolerance.
    """
    return tuple((left, strand_left, right, strand_right)
                 for left, _, strand_left, right, _, strand_right in model.breakpoints)


def _vector(model) -> Tuple[int, ...]:
    """Exit and entry coordinate of every junction, in order.

    Exit and entry are carried as separate elements so they are compared
    independently: a donor that moved 1 bp and an acceptor that moved 5 must
    not average into one distance.
    """
    return tuple(coordinate
                 for _, exit_position, _, _, entry_position, _ in model.breakpoints
                 for coordinate in (exit_position, entry_position))


def _distance(first: Tuple[int, ...], second: Tuple[int, ...]) -> int:
    return max((abs(left - right) for left, right in zip(first, second)), default=0)


def _family_id(topology, vectors, tolerance) -> str:
    payload = json.dumps(
        [ALGORITHM_VERSION, tolerance, [list(edge) for edge in topology],
         [list(vector) for vector in sorted(vectors)]],
        separators=(",", ":"), sort_keys=True,
    )
    return "EFAM_" + hashlib.sha256(payload.encode()).hexdigest()[:16]


def assign_event_families(models: Iterable, *, breakpoint_tolerance: int = 6
                          ) -> Tuple[EventFamily, ...]:
    """Partition reported models by junction geometry, bounded diameter.

    Bounded diameter (every pair in a family within tolerance) does not by
    itself pick a partition -- at tolerance 6, {0,6}|{12} and {0}|{6,12} both
    satisfy it -- so the algorithm is specified, not merely the property:

    1. bucket by ordered topology (junction count, contigs and strands);
    2. take the DISTINCT breakpoint vectors of the bucket, lexicographically;
    3. admit each vector only to families where it is within tolerance of
       every current member;
    4. among the eligible, take the smallest maximum distance, breaking ties
       on the family's smallest vector;
    5. attach every model holding that vector together.

    Partitioning distinct vectors rather than models is what makes step 4
    independent of how many models happen to sit at a geometry: publishing
    another isoform at an existing junction cannot repartition the bucket.
    """
    if breakpoint_tolerance < 0:
        raise ValueError("breakpoint tolerance must not be negative")
    models = list(models)
    by_topology: Dict[tuple, Dict[tuple, list]] = {}
    for model in models:
        if not model.breakpoints:
            # A reported model always spans at least two arms, so it always has
            # a junction. Were that ever to change, an empty vector must not
            # collapse every such model into one family.
            raise ValueError(f"model {model.model_id} has no junction to group on")
        by_topology.setdefault(_topology(model), {}).setdefault(_vector(model), []).append(model)
    families = []
    for topology in sorted(by_topology):
        placements = by_topology[topology]
        groups: list[list[tuple]] = []
        for vector in sorted(placements):
            eligible = [(max(_distance(vector, member) for member in members), min(members), members)
                        for members in groups
                        if all(_distance(vector, member) <= breakpoint_tolerance for member in members)]
            if eligible:
                min(eligible, key=lambda candidate: candidate[:2])[2].append(vector)
            else:
                groups.append([vector])
        for vectors in groups:
            members = [model for vector in sorted(vectors) for model in placements[vector]]
            members.sort(key=lambda model: model.model_id)
            support = {read_id for model in members
                       for read_id in model.unique_read_ids + model.compatible_read_ids}
            families.append(EventFamily(
                _family_id(topology, vectors, breakpoint_tolerance),
                tuple(model.model_id for model in members),
                topology, tuple(sorted(vectors)), tuple(sorted(support)),
            ))
    return tuple(sorted(families, key=lambda family: family.family_id))


def families_by_model(families: Iterable[EventFamily]) -> Dict[str, EventFamily]:
    """model_id -> its family. Every reported model belongs to exactly one."""
    out: Dict[str, EventFamily] = {}
    for family in families:
        for model_id in family.model_ids:
            if model_id in out:
                raise ValueError(f"model {model_id} is in more than one event family")
            out[model_id] = family
    return out
