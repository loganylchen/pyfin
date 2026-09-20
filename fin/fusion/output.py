"""Deterministic fusion artifacts and atomic, failure-visible publication."""

from __future__ import annotations

import csv
import hashlib
import json
from dataclasses import asdict
from pathlib import Path
import tempfile

from fin.fusion.chains import coordinate_key
from fin.fusion.validation import stitched_template


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def correction_id(observed):
    """Identify a correction record by the observed path it was derived from.

    The single expression, called from both places that publish an id here and
    from the junction-evidence exporter, so the two artifacts cannot drift into
    disagreeing about what an id is. The callers still pass their own observed
    path -- the corrections loop its correction's, the evaluated loop the
    item's -- because those are genuinely different arguments.
    """
    return hashlib.sha256(_json(asdict(observed)).encode()).hexdigest()[:20]


def _write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(_json(value) + "\n")
    temporary.replace(path)


def _hash_file(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


class FusionPublication:
    """Switch one directory pointer; a failed rerun cannot expose old calls."""

    def __init__(self, output_dir):
        self.root = Path(output_dir)
        self.lock = None
        self.work_dir = None
        self.state_dir = None
        self.completed = False

    def _point_to(self, directory):
        pointer = self.root / "fusion"
        if pointer.exists() and not pointer.is_symlink():
            raise ValueError("fusion output path exists and is not a managed symlink")
        temporary = self.root / (".fusion-pointer-" + directory.name)
        try:
            temporary.symlink_to(directory.relative_to(self.root), target_is_directory=True)
            temporary.replace(pointer)
        finally:
            temporary.unlink(missing_ok=True)

    def __enter__(self):
        import fcntl

        self.root.mkdir(parents=True, exist_ok=True)
        self.lock = (self.root / ".fusion.lock").open("a")
        try:
            fcntl.flock(self.lock.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            runs = self.root / ".fusion-runs"
            runs.mkdir(exist_ok=True)
            self.state_dir = Path(tempfile.mkdtemp(prefix="status-", dir=runs))
            self.work_dir = Path(tempfile.mkdtemp(prefix="run-", dir=runs))
            _write_json(self.state_dir / "fusion_manifest.json", {"schema_version": 1, "status": "running"})
            self._point_to(self.state_dir)
            return self
        except BaseException:
            self.lock.close()
            raise

    def complete(self, manifest):
        payload = dict(manifest, schema_version=1, status="complete")
        payload["files"] = {path.name: _hash_file(path)
                            for path in sorted(self.work_dir.iterdir()) if path.is_file()}
        _write_json(self.work_dir / "fusion_manifest.json", payload)
        self._point_to(self.work_dir)
        self.completed = True

    def __exit__(self, exc_type, exc_value, traceback):
        try:
            if not self.completed:
                _write_json(self.state_dir / "fusion_manifest.json", {
                    "schema_version": 1, "status": "failed",
                    "error_type": exc_type.__name__ if exc_type else "IncompleteFusionRun",
                    "error": str(exc_value) if exc_value else "fusion results were not published",
                })
        finally:
            if self.lock is not None:
                self.lock.close()


def write_fusion_outputs(analysis, genome, output_dir):
    """Emit complete chains, all N-1 edges, coordinate alternatives and audits."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    reads = {read.read_id: read for read in analysis.reads}
    fields = ("chain_id", "status", "event_class", "gene_path", "arm_count", "full_read_support",
              "full_read_ids", "partial_read_ids", "edge_support", "coordinate_status", "coordinate_alternatives",
              "signal_status",
              # appended: the support gate on its own, for EVERY event class.
              # `status` is a legacy composite -- it reports the gate for
              # chains that fail it and geometry or coordinate resolution for
              # chains that pass. The gate is still recoverable from it, but
              # only by knowing which other values it can take; this column is
              # the stable thing to filter on. Note bedpe's `score` is CAPPED
              # EDGE support and counts partial evidence, so it is not a
              # substitute for full-chain support -- join on chain_id instead.
              "chain_support_status")
    with (output_dir / "fusions.tsv").open("w") as table, (output_dir / "fusions.bedpe").open("w") as bedpe, \
            (output_dir / "fusion_transcripts.fa").open("w") as fasta, \
            (output_dir / "fusion_arms.jsonl").open("w") as arms_output:
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        bed_writer = csv.writer(bedpe, delimiter="\t", lineterminator="\n")
        bed_writer.writerow(("#chrom1", "start1", "end1", "chrom2", "start2", "end2", "name", "score",
                             "strand1", "strand2", "chain_id", "edge_index", "status", "coordinate_status"))
        for chain in analysis.chains:
            path = chain.representative
            coordinates = {}
            geometries = {}
            for read_id in chain.full_support_read_ids:
                for alternative in reads[read_id].selected:
                    coordinates.setdefault(coordinate_key(alternative), set()).add(read_id)
                    geometry = tuple((arm.chrom, arm.strand, arm.gene_ids, arm.exons) for arm in alternative.arms)
                    geometries.setdefault(geometry, alternative)
            coordinate_alternatives = [{"positions": positions, "read_ids": sorted(read_ids)}
                                       for positions, read_ids in sorted(coordinates.items())]
            signal_statuses = [audit["status"] for read_id in chain.full_support_read_ids
                               for audit in reads[read_id].signal if audit.get("retained_alternative", True)]
            signal_status = ("not_scored" if not signal_statuses else
                             "signal_supported" if set(signal_statuses) == {"supported"} else
                             "partially_signal_supported" if "supported" in signal_statuses else
                             "sequence_only_signal_abstained")
            writer.writerow((chain.chain_id, chain.status, path.event_class, _json(path.gene_path), len(path.arms),
                             len(chain.full_support_read_ids), _json(chain.full_support_read_ids),
                             _json(chain.partial_support_read_ids),
                             _json([len(read_ids) for read_ids in chain.edge_support_read_ids]),
                             chain.coordinate_status, _json(coordinate_alternatives), signal_status,
                             chain.support_status))
            for index, (left, right) in enumerate(zip(path.arms, path.arms[1:])):
                left_base = left.exit_position - 1 if left.strand == "+" else left.exit_position
                right_base = right.entry_position if right.strand == "+" else right.entry_position - 1
                bed_writer.writerow((left.chrom, left_base, left_base + 1, right.chrom, right_base, right_base + 1,
                                     f"{chain.chain_id}.edge{index + 1}",
                                     min(1000, len(chain.edge_support_read_ids[index])), left.strand, right.strand,
                                     chain.chain_id, index + 1, chain.status, chain.coordinate_status))
            emitted_templates = set()
            for geometry, alternative in sorted(geometries.items()):
                sequence, _ = stitched_template(alternative, genome)
                if not sequence:
                    raise ValueError("fusion output template is unavailable")
                template_id = chain.chain_id + ".template_" + hashlib.sha256(sequence.encode()).hexdigest()[:16]
                geometry_id = chain.chain_id + ".path_" + hashlib.sha256(_json(geometry).encode()).hexdigest()[:16]
                if template_id not in emitted_templates:
                    fasta.write(f">{template_id} status={chain.status} coordinates={chain.coordinate_status}\n")
                    for offset in range(0, len(sequence), 80):
                        fasta.write(sequence[offset:offset + 80] + "\n")
                    emitted_templates.add(template_id)
                arms_output.write(_json({"chain_id": chain.chain_id, "path_id": geometry_id,
                                         "template_id": template_id, "status": chain.status,
                                         "coordinate_status": chain.coordinate_status,
                                         "arms": [{"arm_index": index + 1, "chrom": arm.chrom,
                                                   "strand": arm.strand, "gene_ids": arm.gene_ids,
                                                   "exons": arm.exons}
                                                  for index, arm in enumerate(alternative.arms)]}) + "\n")
    with (output_dir / "fusion_reads.jsonl").open("w") as read_output, \
            (output_dir / "junction_corrections.jsonl").open("w") as correction_output, \
            (output_dir / "junction_signal.jsonl").open("w") as signal_output:
        for read in analysis.reads:
            for audit in read.signal:
                signal_output.write(_json({"read_id": read.read_id, **audit}) + "\n")
            evaluated = []
            seen_corrections = set()
            for correction in read.corrections:
                identifier = correction_id(correction.observed)
                if identifier not in seen_corrections:
                    correction_output.write(_json({"read_id": read.read_id, "correction_id": identifier,
                                                   **asdict(correction)}) + "\n")
                    seen_corrections.add(identifier)
            for item in read.evaluated:
                identifier = None
                if item.correction is not None:
                    identifier = correction_id(item.observed)
                    if identifier not in seen_corrections:
                        correction_output.write(_json({"read_id": read.read_id, "correction_id": identifier,
                                                       **asdict(item.correction)}) + "\n")
                        seen_corrections.add(identifier)
                evaluated.append({"observed": asdict(item.observed), "proposed": asdict(item.proposed),
                                  "correction_id": identifier, "validation": asdict(item.validation),
                                  "min_mapping_quality": item.min_mapping_quality})
            read_output.write(_json({"read_id": read.read_id, "status": read.status, "reason": read.reason,
                                     "path_count": read.path_count, "search_windows": read.search_windows,
                                     "correction_ids": sorted(seen_corrections),
                                     "signal": read.signal, "assembly_junctions": read.assembly_junctions,
                                     "sa_supersessions": list(read.sa_supersessions),
                                     "coverage_domain": list(read.coverage_domain),
                                     "selected": [asdict(path) for path in read.selected], "evaluated": evaluated}) + "\n")
