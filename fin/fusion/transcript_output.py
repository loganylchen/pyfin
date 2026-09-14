"""Fused-transcript products: model table, two FASTAs, exon table, evidence.

Added beside the chain-level artifacts, never replacing them, so the six
existing biological outputs stay byte-identical. Two FASTAs are deliberate:
the evidence-constrained one contains only structure a molecule supports; the
annotated-completion one extends outer arms to compatible annotated ends and
says so in every header. No downstream user should be able to mistake one
for the other, so they never share a record id.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

from dataclasses import asdict

from fin.fusion.transcript_models import (
    FusedTranscriptModel, completion_exons, model_sequence,
)


def _json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _token(value):
    """Header-safe: no whitespace, `=` or `;` inside a metadata value."""
    return str(value).replace(" ", "_").replace("=", ":").replace(";", ",")


def _gene_label(arm):
    names = arm.resolution.gene_names or arm.resolution.gene_ids or ()
    return "|".join(names) if names else f"{arm.chrom}:{arm.exons[0][0]}-{arm.exons[-1][1]}"


def _breakpoint_label(model):
    return ",".join(f"{cl}:{el}:{sl}>{cr}:{er}:{sr}" for cl, el, sl, cr, er, sr in model.breakpoints)


def _junction_names(model):
    """Conventional exon-rank names per breakpoint, e.g. BCR:e13>ABL1:e2."""
    names = []
    for left, right in zip(model.arms, model.arms[1:]):
        left_ranks = {rank for _, rank in left.resolution.breakpoint_exon_ranks}
        # the right arm's entry exon rank: rank of the exon containing entry
        right_rank = _entry_rank(right)
        left_label = "e" + "/".join(str(r) for r in sorted(left_ranks)) if left_ranks else "?"
        right_label = f"e{right_rank}" if right_rank is not None else "?"
        names.append(f"{_gene_label(left)}:{left_label}>{_gene_label(right)}:{right_label}")
    return names


def _entry_rank(arm):
    values = {rank for _, rank in arm.resolution.entry_exon_ranks}
    return "/".join(str(r) for r in sorted(values)) if values else None


def _header(record_id, model, *, kind, length, extension=None):
    tx = ">".join(
        f"{_gene_label(arm)}:[{','.join(arm.resolution.compatible_transcript_ids) or '-'}]"
        for arm in model.arms)
    fields = [
        f"genes={_token('>'.join(_gene_label(a) for a in model.arms))}",
        f"junction={_token(';'.join(_junction_names(model)))}",
        f"tx={_token(tx)}",
        f"bp={_token(_breakpoint_label(model))}",
        f"arms={len(model.arms)}",
        f"model={model.model_id}",
        f"parent_chain={model.primary_chain_id}",
        f"support_unique={len(model.unique_read_ids)}",
        f"support_compatible={len(model.compatible_read_ids)}",
        f"support_partial={len(model.partial_read_ids)}",
        f"nearby_not_counted={len(model.nearby_read_ids)}",
        f"envelope={model.envelope_status}",
        f"end5={model.five_prime_status}",
        f"end3={model.three_prime_status}",
        f"coordinates={model.coordinate_status}",
        f"coordinate_alternatives={len(model.coordinate_alternatives)}",
        f"class={model.event_class}",
        f"kind={kind}",
        f"length={length}",
        "sequence_source=reference_projected",
    ]
    if extension is not None:
        fields.append(f"extension_bases_5p={extension[0]}")
        fields.append(f"extension_bases_3p={extension[1]}")
        fields.append("extension_observed=false")
    return f">{record_id} " + " ".join(fields)


def write_transcript_outputs(models, collapsed, genome, output_dir):
    """Write fusion_models.tsv / .fa / .annotated.fa / _exons.tsv / _evidence.jsonl / _collapsed.jsonl.

    Any sequence that cannot be projected from the reference -- observed OR
    annotated completion -- raises, so the publication is marked failed rather
    than quietly missing a record.
    """
    output_dir = Path(output_dir)
    fields = ("model_id", "sequence_id", "primary_chain_id", "source_chain_ids", "alternative_group_id",
              "arm_count", "gene_path", "gene_names", "junction_names", "breakpoints",
              "compatible_transcript_ids", "arm_match_classes",
              "support_unique", "support_compatible", "support_partial", "chain_edge_support",
              "unique_read_ids", "compatible_read_ids", "partial_read_ids", "nearby_read_ids",
              "envelope_status", "five_prime_observed_reads", "three_prime_observed_reads",
              "coordinate_status", "coordinate_alternatives", "event_class", "five_prime_status", "three_prime_status",
              "observed_length", "annotated_completion_id", "annotated_completion_length",
              "extension_bases_5p", "extension_bases_3p", "completion_transcript_ids")
    exon_fields = ("model_id", "sequence_id", "record_kind", "arm_index", "exon_index_in_arm",
                   "exon_index_in_model", "chrom", "strand", "start0", "end0", "tx_start0", "tx_end0",
                   "gene_ids", "compatible_transcript_ids", "arm_match_class", "provenance")
    with (output_dir / "fusion_models.tsv").open("w") as table, \
            (output_dir / "fusion_models.fa").open("w") as observed_fa, \
            (output_dir / "fusion_models.annotated.fa").open("w") as annotated_fa, \
            (output_dir / "fusion_model_exons.tsv").open("w") as exon_table, \
            (output_dir / "fusion_model_evidence.jsonl").open("w") as evidence, \
            (output_dir / "fusion_model_collapsed.jsonl").open("w") as collapsed_out:
        for item in collapsed:
            collapsed_out.write(_json(asdict(item)) + "\n")
        writer = csv.writer(table, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        exon_writer = csv.writer(exon_table, delimiter="\t", lineterminator="\n")
        exon_writer.writerow(exon_fields)
        for model in models:
            sequence, rows = model_sequence(model.arms, genome)
            if not sequence:
                raise ValueError("fused transcript sequence is unavailable from the reference")
            sequence_id = model.model_id + ".obs"
            observed_fa.write(_header(sequence_id, model, kind="evidence_constrained", length=len(sequence)) + "\n")
            for offset in range(0, len(sequence), 80):
                observed_fa.write(sequence[offset:offset + 80] + "\n")
            model_exon_index = 0
            for arm_index, exon_index, chrom, strand, start, end, tx_start, tx_end in rows:
                model_exon_index += 1
                arm = model.arms[arm_index - 1]
                exon_writer.writerow((model.model_id, sequence_id, "evidence_constrained", arm_index, exon_index,
                                      model_exon_index, chrom, strand, start, end, tx_start, tx_end,
                                      _json(list(arm.resolution.gene_ids)),
                                      _json(list(arm.resolution.compatible_transcript_ids)),
                                      arm.resolution.match_class, "observed"))
            completion = completion_exons(model)
            completion_id = ""
            completion_length = ""
            ext5 = ext3 = ""
            completion_tx = ""
            if completion is not None and completion != tuple(arm.exons for arm in model.arms):
                full, full_rows = model_sequence(model.arms, genome, exons_override=completion)
                if not full:
                    raise ValueError("annotated completion sequence is unavailable from the reference")
                if full:
                    completion_id = model.model_id + ".ann"
                    completion_length = len(full)
                    first, last = model.arms[0], model.arms[-1]
                    ext5 = sum(e - s for s, e in completion[0]) - sum(e - s for s, e in first.exons)
                    ext3 = sum(e - s for s, e in completion[-1]) - sum(e - s for s, e in last.exons)
                    if ext5 < 0 or ext3 < 0:
                        # the resolver only offers a completion that contains
                        # every observed exon; a negative extension means it
                        # did not, and that is a defect, not a record
                        raise ValueError("annotated completion does not contain the observed exons")
                    completion_tx = _json({
                        "first_arm": list(first.resolution.completion_transcript_ids),
                        "last_arm": list(last.resolution.completion_transcript_ids)})
                    annotated_fa.write(_header(completion_id, model, kind="annotated_completion",
                                               length=len(full), extension=(ext5, ext3)) + "\n")
                    for offset in range(0, len(full), 80):
                        annotated_fa.write(full[offset:offset + 80] + "\n")
                    observed_exons = {(arm_index, exon) for arm_index, arm in enumerate(model.arms, start=1)
                                      for exon in arm.exons}
                    model_exon_index = 0
                    for arm_index, exon_index, chrom, strand, start, end, tx_start, tx_end in full_rows:
                        model_exon_index += 1
                        arm = model.arms[arm_index - 1]
                        provenance = ("observed" if (arm_index, (start, end)) in observed_exons else
                                      "annotated_extension" if not any(
                                          s < end and start < e for s, e in arm.exons) else
                                      "observed_partial_annotated_extension")
                        exon_writer.writerow((model.model_id, completion_id, "annotated_completion", arm_index,
                                              exon_index, model_exon_index, chrom, strand, start, end, tx_start,
                                              tx_end, _json(list(arm.resolution.gene_ids)),
                                              _json(list(arm.resolution.compatible_transcript_ids)),
                                              arm.resolution.match_class, provenance))
            writer.writerow((
                model.model_id, sequence_id, model.primary_chain_id, _json(list(model.source_chain_ids)),
                model.alternative_group_id, len(model.arms),
                _json([list(g) for g in model.gene_path]),
                _json([list(arm.resolution.gene_names) for arm in model.arms]),
                _json(_junction_names(model)), _breakpoint_label(model),
                _json([list(arm.resolution.compatible_transcript_ids) for arm in model.arms]),
                _json([arm.resolution.match_class for arm in model.arms]),
                len(model.unique_read_ids), len(model.compatible_read_ids), len(model.partial_read_ids),
                _json(list(model.chain_edge_support)),
                _json(list(model.unique_read_ids)), _json(list(model.compatible_read_ids)),
                _json(list(model.partial_read_ids)), _json(list(model.nearby_read_ids)),
                model.envelope_status,
                len(model.five_prime_observed_read_ids), len(model.three_prime_observed_read_ids),
                model.coordinate_status,
                _json([json.loads(alt) for alt in model.coordinate_alternatives]),
                model.event_class, model.five_prime_status, model.three_prime_status,
                len(sequence), completion_id, completion_length, ext5, ext3, completion_tx,
            ))
            evidence.write(_json({
                "model_id": model.model_id, "alternative_group_id": model.alternative_group_id,
                "primary_chain_id": model.primary_chain_id, "source_chain_ids": list(model.source_chain_ids),
                "unique_read_ids": list(model.unique_read_ids),
                "compatible_read_ids": list(model.compatible_read_ids),
                "partial_read_ids": list(model.partial_read_ids),
                "nearby_read_ids": list(model.nearby_read_ids),
                "envelope_status": model.envelope_status,
                "chain_edge_support": list(model.chain_edge_support),
                "five_prime_observed_read_ids": list(model.five_prime_observed_read_ids),
                "three_prime_observed_read_ids": list(model.three_prime_observed_read_ids),
                "coordinate_alternatives": [json.loads(alt) for alt in model.coordinate_alternatives],
                "nearby_observed_placements": [json.loads(q) for q in model.nearby_observed_placements],
                "arms": [{"arm_index": i + 1, "chrom": arm.chrom, "strand": arm.strand,
                          "exons": [list(e) for e in arm.exons],
                          "match_class": arm.resolution.match_class,
                          "gene_ids": list(arm.resolution.gene_ids),
                          "gene_names": list(arm.resolution.gene_names),
                          "compatible_transcript_ids": list(arm.resolution.compatible_transcript_ids),
                          "completion_exons": [list(e) for e in arm.resolution.completion_exons],
                          "breakpoint_exon_ranks": [list(r) for r in arm.resolution.breakpoint_exon_ranks]}
                         for i, arm in enumerate(model.arms)],
                "interpretation": ("evidence-constrained structure is reference-projected from molecule "
                                   "alignments; annotated completion is a hypothesis, not an observation"),
            }) + "\n")
