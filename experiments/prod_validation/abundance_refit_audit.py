#!/usr/bin/env python3
"""Mechanical identity and mass audit for paired refit-on/off runs."""
from __future__ import annotations

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


GTF_QUANT_ATTR = re.compile(
    r'(?:abundance|confidence|num_reads) "[^"]*"; ?'
)
TSV_QUANT_COLUMNS = {"abundance", "confidence", "num_reads", "tpm", "max_R"}


def _normalized_gtf_records(path: Path) -> list[str]:
    records = []
    for line in path.read_text().splitlines():
        if "\ttranscript\t" in line:
            line = GTF_QUANT_ATTR.sub("", line)
        records.append(line)
    return records


def _read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or ()), list(reader)


def audit(on_dir: Path, off_dir: Path, legacy_dir: Path | None = None) -> dict:
    on_gtf = on_dir / "assembly.gtf"
    off_gtf = off_dir / "assembly.gtf"
    on_tsv = on_dir / "scores.tsv"
    off_tsv = off_dir / "scores.tsv"
    diagnostics_path = on_dir / "abundance_refit.json"
    for path in (on_gtf, off_gtf, on_tsv, off_tsv, diagnostics_path):
        if not path.exists():
            raise FileNotFoundError(path)

    on_gtf_records = _normalized_gtf_records(on_gtf)
    off_gtf_records = _normalized_gtf_records(off_gtf)
    gtf_order_identity = on_gtf_records == off_gtf_records
    gtf_order_differences = sum(
        left != right
        for left, right in zip(on_gtf_records, off_gtf_records)
    ) + abs(len(on_gtf_records) - len(off_gtf_records))
    gtf_structural_identity = (
        Counter(on_gtf_records) == Counter(off_gtf_records)
    )
    if not gtf_structural_identity:
        raise AssertionError(
            "normalized refit-on/off GTF record multisets differ structurally"
        )

    on_header, on_rows = _read_tsv(on_tsv)
    off_header, off_rows = _read_tsv(off_tsv)
    if on_header != off_header:
        raise AssertionError("refit-on/off TSV headers differ")
    if len(on_rows) != len(off_rows):
        raise AssertionError("refit-on/off TSV row counts differ")

    invariant_columns = [
        column for column in on_header if column not in TSV_QUANT_COLUMNS
    ]
    changed = {column: 0 for column in sorted(TSV_QUANT_COLUMNS)}
    for index, (on_row, off_row) in enumerate(zip(on_rows, off_rows)):
        for column in invariant_columns:
            if on_row[column] != off_row[column]:
                raise AssertionError(
                    f"TSV invariant changed at row {index}, column {column}: "
                    f"{off_row[column]!r} -> {on_row[column]!r}"
                )
        for column in changed:
            if on_row[column] != off_row[column]:
                changed[column] += 1

    diagnostics = json.loads(diagnostics_path.read_text())
    if diagnostics.get("effective") is not True:
        raise AssertionError("refit diagnostics do not record effective=true")
    if float(diagnostics.get("mass_balance_error", float("inf"))) > 1e-6:
        raise AssertionError("refit total mass does not close")
    if float(diagnostics.get("max_per_read_conservation_error", float("inf"))) > 1e-7:
        raise AssertionError("refit per-read mass does not close")

    legacy_gtf_identity = None
    legacy_tsv_identity = None
    if legacy_dir is not None:
        legacy_gtf_identity = off_gtf.read_bytes() == (legacy_dir / "assembly.gtf").read_bytes()
        legacy_tsv_identity = off_tsv.read_bytes() == (legacy_dir / "scores.tsv").read_bytes()
        if not (legacy_gtf_identity and legacy_tsv_identity):
            raise AssertionError("refit-off output is not byte-identical to legacy")

    report = {
        "schema_version": 1,
        "refit_on": str(on_dir),
        "refit_off": str(off_dir),
        "legacy": str(legacy_dir) if legacy_dir is not None else None,
        "gtf_comparison": "normalized_record_multiset",
        "gtf_structural_identity": gtf_structural_identity,
        "gtf_record_multiset_identity": gtf_structural_identity,
        "gtf_line_order_identity": gtf_order_identity,
        "gtf_line_order_differences": gtf_order_differences,
        "tsv_rows": len(on_rows),
        "tsv_invariant_columns": invariant_columns,
        "tsv_quant_columns_changed": changed,
        "legacy_gtf_byte_identity": legacy_gtf_identity,
        "legacy_tsv_byte_identity": legacy_tsv_identity,
        "mass_balance_error": diagnostics["mass_balance_error"],
        "max_per_read_conservation_error": diagnostics[
            "max_per_read_conservation_error"
        ],
        "selection_orphaned_reads": diagnostics["selection_orphaned_reads"],
        "unassigned_mass": diagnostics["unassigned_mass"],
        "rescued_mass": diagnostics["rescued_mass"],
    }
    diagnostics["structural_identity"] = "passed"
    tmp = diagnostics_path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(diagnostics, indent=2, sort_keys=True) + "\n")
    tmp.replace(diagnostics_path)
    report_path = on_dir / "refit_identity.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--on-dir", required=True, type=Path)
    parser.add_argument("--off-dir", required=True, type=Path)
    parser.add_argument("--legacy-dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.on_dir, args.off_dir, args.legacy_dir), indent=2))


if __name__ == "__main__":
    main()
