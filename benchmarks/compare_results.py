#!/usr/bin/env python
"""Compare benchmark results across tools into a flat TSV.

Output schema (tab-separated, with header row):
  tool    metric    value

Metrics collected per tool:
  - All key/value pairs from result.json (e.g. status, dataset)
  - num_transcripts: count of 'transcript' feature lines in any *.gtf found
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def collect_tool_rows(tool_dir: Path) -> list[tuple[str, str, str]]:
    """Return (tool, metric, value) rows for a single tool output directory."""
    rows: list[tuple[str, str, str]] = []
    name = tool_dir.name

    result_file = tool_dir / "result.json"
    if result_file.exists():
        try:
            data = json.loads(result_file.read_text())
            for k, v in data.items():
                rows.append((name, k, str(v)))
        except json.JSONDecodeError:
            rows.append((name, "error", "invalid json"))

    for gtf in sorted(tool_dir.glob("*.gtf")):
        n = sum(
            1
            for line in gtf.read_text().splitlines()
            if line and not line.startswith("#") and "\ttranscript\t" in line
        )
        rows.append((name, "num_transcripts", str(n)))

    return rows


def main() -> None:
    p = argparse.ArgumentParser(
        description="Aggregate per-tool benchmark outputs into a comparison TSV."
    )
    p.add_argument("--input-dir", required=True, help="Directory containing per-tool subdirs")
    p.add_argument("--output", required=True, help="Path for output TSV file")
    args = p.parse_args()

    in_dir = Path(args.input_dir)
    if not in_dir.is_dir():
        p.error(f"--input-dir does not exist: {in_dir}")

    rows: list[tuple[str, str, str]] = [("tool", "metric", "value")]
    for tool_dir in sorted(in_dir.iterdir()):
        if not tool_dir.is_dir():
            continue
        rows.extend(collect_tool_rows(tool_dir))

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join("\t".join(r) for r in rows) + "\n")
    print(f"Wrote {args.output} with {len(rows) - 1} data rows")


if __name__ == "__main__":
    main()
