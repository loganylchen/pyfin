#!/usr/bin/env python3
"""Reproducible SIRV/real-dRNA profile ablation runner.

The primary metric is NanoCount-T3 honest transcript F1. Standard gffcompare
transcript sensitivity/precision/F1 are retained for compatibility with older
reports. Every pyfin run writes its own run_manifest.json.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
CODE_ROOT = Path(os.environ.get("PYFIN_CODE_ROOT", REPO)).resolve()
PROD = REPO / "experiments" / "prod_validation"
SIF = PROD / "_img" / "pyfin_gpu_e268c9b.sif"
GFFCOMPARE_IMAGE = "quay.io/biocontainers/gffcompare:0.12.10--h9948957_0"
ENST = re.compile(r"ENST\d+")
NANORNA_ROOT = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans")
REAL_SAMPLE_SHEET = NANORNA_ROOT / "benchmark/sgnex/batches/gencode_full_sweep/samples.tsv"


SIRV_VARIANTS = {
    "baseline": [],
    "overlap_locus": ["--isoform-fraction-locus", "overlap"],
    "legacy_profile": ["--no-floor-gtf-abundance", "--min-polya5p-reads", "1"],
    "m2_sum": ["--m2-metric", "summed_llr"],
    "m2_off": ["--m2-metric", "off"],
    "abundance_1": ["--min-abundance", "1"],
    "abundance_2": ["--min-abundance", "2"],
    "finalize_off": ["--min-fulllen-fraction", "0", "--min-polya5p-reads", "0"],
    "fulllen_off": ["--min-fulllen-fraction", "0"],
    "polya_off": ["--min-polya5p-reads", "0"],
    "junction_1": ["--novel-junction-min-reads", "1"],
    "containment_off": ["--no-containment-cluster"],
    "recheck_off": ["--no-m2-cluster-recheck"],
    "canonical_off": ["--no-canonical-gate"],
    "softmass_off": ["--max-soft-mass-ratio", "0"],
    "isoform_fraction_off": ["--min-isoform-fraction", "0"],
    "floor_gtf_on": ["--floor-gtf-abundance"],
    "best_base": ["--min-polya5p-reads", "0", "--floor-gtf-abundance"],
    "best_mean": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "mean"],
    "best_m2_off": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "off"],
    "best_finalize_off": ["--min-fulllen-fraction", "0", "--min-polya5p-reads", "0", "--floor-gtf-abundance"],
    "best_sum_m1_f4": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "1", "--m2-summed-llr-flank", "4"],
    "best_sum_m1_f6": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "1", "--m2-summed-llr-flank", "6"],
    "best_sum_m1_f8": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "1", "--m2-summed-llr-flank", "8"],
    "best_sum_m2_f4": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "2", "--m2-summed-llr-flank", "4"],
    "best_sum_m2_f6": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "2", "--m2-summed-llr-flank", "6"],
    "best_sum_m2_f8": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "2", "--m2-summed-llr-flank", "8"],
    "best_sum_m3_f4": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "3", "--m2-summed-llr-flank", "4"],
    "best_sum_m3_f6": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "3", "--m2-summed-llr-flank", "6"],
    "best_sum_m3_f8": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "3", "--m2-summed-llr-flank", "8"],
    "best_sum_m4_f4": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "4", "--m2-summed-llr-flank", "4"],
    "best_sum_m4_f6": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "4", "--m2-summed-llr-flank", "6"],
    "best_sum_m4_f8": ["--min-polya5p-reads", "0", "--floor-gtf-abundance", "--m2-metric", "summed_llr", "--m2-summed-llr-margin", "4", "--m2-summed-llr-flank", "8"],
}

REAL_VARIANTS = {
    "baseline": [],
    "overlap_locus": ["--isoform-fraction-locus", "overlap"],
    "mono3": ["--drop-mono-exon-novel", "--min-mono-exon-reads", "3"],
    "mono5": ["--drop-mono-exon-novel", "--min-mono-exon-reads", "5"],
    "mono5_snap6": ["--drop-mono-exon-novel", "--min-mono-exon-reads", "5", "--junction-snap", "--junction-snap-tolerance", "6", "--junction-snap-min-support", "2", "--junction-snap-min-ratio", "2"],
    "mono5_r2": ["--drop-mono-exon-novel", "--min-mono-exon-reads", "5", "--min-novel-reads", "2"],
    "mono5_r2_snap6": ["--drop-mono-exon-novel", "--min-mono-exon-reads", "5", "--min-novel-reads", "2", "--junction-snap", "--junction-snap-tolerance", "6", "--junction-snap-min-support", "2", "--junction-snap-min-ratio", "2"],
    "mono5_r2_ab15": ["--drop-mono-exon-novel", "--min-mono-exon-reads", "5", "--min-novel-reads", "2", "--min-abundance", "1.5"],
    "m2_sum": ["--m2-metric", "summed_llr"],
    "sum_m1_f4": ["--m2-metric", "summed_llr", "--m2-summed-llr-margin", "1", "--m2-summed-llr-flank", "4"],
    "sum_m1_f6": ["--m2-metric", "summed_llr", "--m2-summed-llr-margin", "1", "--m2-summed-llr-flank", "6"],
    "sum_m1_f8": ["--m2-metric", "summed_llr", "--m2-summed-llr-margin", "1", "--m2-summed-llr-flank", "8"],
    "sum_m2_f4": ["--m2-metric", "summed_llr", "--m2-summed-llr-margin", "2", "--m2-summed-llr-flank", "4"],
    "sum_m2_f6": ["--m2-metric", "summed_llr", "--m2-summed-llr-margin", "2", "--m2-summed-llr-flank", "6"],
    "sum_m2_f8": ["--m2-metric", "summed_llr", "--m2-summed-llr-margin", "2", "--m2-summed-llr-flank", "8"],
    "sum_m3_f4": ["--m2-metric", "summed_llr", "--m2-summed-llr-margin", "3", "--m2-summed-llr-flank", "4"],
    "sum_m3_f6": ["--m2-metric", "summed_llr", "--m2-summed-llr-margin", "3", "--m2-summed-llr-flank", "6"],
    "sum_m3_f8": ["--m2-metric", "summed_llr", "--m2-summed-llr-margin", "3", "--m2-summed-llr-flank", "8"],
    "m2_off": ["--m2-metric", "off"],
    "abundance_0": ["--min-abundance", "0"],
    "abundance_2": ["--min-abundance", "2"],
    "abundance_3": ["--min-abundance", "3"],
    "finalize_sirv": ["--min-fulllen-fraction", "0.1", "--min-polya5p-reads", "1"],
    "fulllen_on": ["--min-fulllen-fraction", "0.1"],
    "polya_on": ["--min-polya5p-reads", "1"],
    "junction_1": ["--novel-junction-min-reads", "1"],
    "containment_off": ["--no-containment-cluster"],
    "recheck_off": ["--no-m2-cluster-recheck"],
    "canonical_off": ["--no-canonical-gate"],
    "softmass_off": ["--max-soft-mass-ratio", "0"],
    "isoform_fraction_off": ["--min-isoform-fraction", "0"],
}


@dataclass(frozen=True)
class Dataset:
    domain: str
    sample: str
    guide: str
    bam: Path
    genome: Path
    fastq: Path
    signal: Path
    annotation: Path | None
    truth: Path
    nanocount: Path


def sirv_datasets(samples: list[str], guides: list[str]) -> list[Dataset]:
    root = PROD / "sirv4"
    out = []
    for sample in samples:
        stage = root / sample / "stage"
        for guide in guides:
            annotation = root / "_ref" / guide / "annotation.gtf"
            transcript_count = 0
            with annotation.open(errors="ignore") as handle:
                for line in handle:
                    if "\ttranscript\t" in line:
                        transcript_count += 1
                        if transcript_count >= 2:
                            break
            if transcript_count < 2:
                annotation = None
            out.append(Dataset(
                domain="sirv",
                sample=sample,
                guide=guide,
                bam=stage / "input.bam",
                genome=stage / "genome.fa",
                fastq=stage / "reads.fq.gz",
                signal=stage / "signal.blow5",
                annotation=annotation,
                truth=root / "_ref" / "full" / "annotation.gtf",
                nanocount=stage / "nanocount.tsv",
            ))
    return out


def _real_raw_inputs(sample: str) -> tuple[Path, Path]:
    """Return full raw inputs when mapped-subset stage links are unavailable."""
    with REAL_SAMPLE_SHEET.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("sample") == sample:
                return (
                    (NANORNA_ROOT / row["fastq"]).resolve(),
                    (NANORNA_ROOT / row["signal_path"]).resolve(),
                )
    raise ValueError(f"sample {sample!r} is absent from {REAL_SAMPLE_SHEET}")


def real_dataset(sample: str = "gencode_p00") -> Dataset:
    root = PROD / "gencode"
    stage = root / "_p00val" / "stage" if sample == "gencode_p00" else root / sample / "stage"
    nanocount_sample = (
        "SGNex_H9_directRNA_replicate2_run2"
        if sample == "gencode_p00" else sample
    )
    fastq = stage / "reads.fq.gz"
    signal = stage / "signal.blow5"
    if sample != "gencode_p00" and (not fastq.exists() or not signal.exists()):
        fastq, signal = _real_raw_inputs(sample)
    return Dataset(
        domain="real-drna",
        sample=sample,
        guide="p00",
        bam=stage / "input.bam",
        genome=stage / "genome.fa",
        fastq=fastq,
        signal=signal,
        annotation=None,
        truth=root / "_p00val" / "_gc" / "truth.gtf",
        nanocount=Path(
            "/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/"
            f"results/gencode_full_sweep/{nanocount_sample}/"
            "full/assembly/nanocount.tsv"
        ),
    )


def normalized_ref_id(value: str) -> str:
    match = ENST.search(value)
    if match:
        return match.group(0)
    return value.rsplit("|", 1)[-1].split(".", 1)[0]


def expressed_truth(path: Path, threshold: float = 3.0) -> set[str]:
    keep = set()
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                count = float(row["est_count"])
            except (KeyError, TypeError, ValueError):
                continue
            if count >= threshold:
                keep.add(normalized_ref_id(row["transcript_name"]))
    return keep


def parse_metrics(
    out_dir: Path,
    truth_by_threshold: dict[int, set[str]],
    runtime: float,
) -> dict[str, object]:
    stats = (out_dir / "gc.stats").read_text(errors="ignore")
    match = re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)", stats)
    if not match:
        raise RuntimeError(f"missing transcript metrics in {out_dir / 'gc.stats'}")
    sensitivity, precision = map(float, match.groups())
    standard_f1 = (
        2 * sensitivity * precision / (sensitivity + precision)
        if sensitivity + precision else 0.0
    )
    matched = set()
    with (out_dir / "gc.tracking").open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 4 and fields[3] == "=":
                matched.add(normalized_ref_id(fields[2]))
    honest = {}
    for threshold, truth in sorted(truth_by_threshold.items()):
        expressed_matched = matched & truth
        honest_precision = (
            precision * len(expressed_matched) / len(matched) if matched else 0.0
        )
        corrected_recall = (
            100.0 * len(expressed_matched) / len(truth) if truth else 0.0
        )
        honest_f1 = (
            2 * honest_precision * corrected_recall
            / (honest_precision + corrected_recall)
            if honest_precision + corrected_recall else 0.0
        )
        honest.update({
            f"expressed_truth_t{threshold}": len(truth),
            f"expressed_matched_t{threshold}": len(expressed_matched),
            f"honest_precision_t{threshold}": honest_precision,
            f"corrected_recall_t{threshold}": corrected_recall,
            f"honest_f1_t{threshold}": honest_f1,
        })
    output_gtf = out_dir / "assembly.gtf"
    transcripts = sum("\ttranscript\t" in line for line in output_gtf.open())
    return {
        "transcripts": transcripts,
        "tx_sensitivity": sensitivity,
        "tx_precision": precision,
        "tx_f1": standard_f1,
        "matched": len(matched),
        **honest,
        "runtime_seconds": runtime,
    }


def run_one(
    dataset: Dataset,
    profile: str,
    variant: str,
    flags: list[str],
    root: Path,
    *,
    score: bool = True,
) -> dict[str, object]:
    out_dir = root / variant / dataset.sample / dataset.guide
    out_dir.mkdir(parents=True, exist_ok=True)
    runtime_file = out_dir / "runtime_seconds.txt"
    if not (out_dir / "assembly.gtf").exists():
        bind_dirs = {Path("/SSD")}
        for input_path in (
            dataset.bam,
            dataset.genome,
            dataset.fastq,
            dataset.signal,
            dataset.annotation,
        ):
            if input_path is None:
                continue
            resolved = input_path.resolve(strict=False)
            if str(resolved).startswith("/autofs/"):
                # /autofs is a nested automount: binding its top-level mount
                # hides mounted children. Bind each realized input directory.
                bind_dirs.add(resolved.parent)
        binds = [item for path in sorted(bind_dirs) for item in ("-B", str(path))]
        cli_module = os.environ.get("PYFIN_CLI_MODULE", "fin.cli")
        command = [
            "singularity", "exec", "--nv", *binds, str(SIF),
            "env", f"PYTHONPATH={CODE_ROOT}", "/usr/bin/python3.10", "-m", cli_module,
            "--profile", profile,
            "--bam", str(dataset.bam),
            "--genome", str(dataset.genome),
            "--fastq", str(dataset.fastq),
            "--signal", str(dataset.signal),
            "--signal-format", "slow5",
            "--output-dir", str(out_dir),
            "--quant-mode", "m2_em",
        ]
        if dataset.domain == "sirv":
            command.append("--gpu")
        else:
            command.extend((
                "--no-gpu",
                "--threads", os.environ.get("PYFIN_REAL_THREADS", "8"),
            ))
        if dataset.annotation is not None:
            command.extend(("--gtf", str(dataset.annotation)))
        command.extend(flags)
        start = time.monotonic()
        with (out_dir / "run.log").open("w") as log:
            subprocess.run(command, cwd=REPO, stdout=log, stderr=subprocess.STDOUT, check=True)
        runtime_file.write_text(f"{time.monotonic() - start:.6f}\n")
    runtime = float(runtime_file.read_text()) if runtime_file.exists() else 0.0
    if not score:
        return {
            "domain": dataset.domain,
            "profile": profile,
            "variant": variant,
            "sample": dataset.sample,
            "guide": dataset.guide,
            "flags": " ".join(flags),
            "runtime_seconds": runtime,
        }
    if not (out_dir / "gc.stats").exists():
        subprocess.run([
            "docker", "run", "--rm", "-u", f"{os.getuid()}:{os.getgid()}",
            "-v", "/SSD:/SSD", GFFCOMPARE_IMAGE,
            "gffcompare", "-r", str(dataset.truth), "-o", str(out_dir / "gc"),
            str(out_dir / "assembly.gtf"),
        ], cwd=REPO, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True)
    metrics = parse_metrics(
        out_dir,
        {
            1: expressed_truth(dataset.nanocount, threshold=1.0),
            3: expressed_truth(dataset.nanocount, threshold=3.0),
        },
        runtime,
    )
    return {
        "domain": dataset.domain,
        "profile": profile,
        "variant": variant,
        "sample": dataset.sample,
        "guide": dataset.guide,
        "flags": " ".join(flags),
        **metrics,
    }


def write_results(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else []
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--domain", choices=("sirv", "real-drna"), required=True)
    parser.add_argument(
        "--profile",
        choices=("sirv", "real-drna", "real-drna-precision", "custom"),
    )
    parser.add_argument("--samples")
    parser.add_argument("--guides", default="p00,full")
    parser.add_argument("--variants", default="all")
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--run-only", action="store_true")
    args = parser.parse_args()

    if args.domain == "sirv":
        profile = args.profile or "sirv"
        variants = SIRV_VARIANTS
        samples = (args.samples or (
            "SGNex_H9_directRNA_replicate4_run2,"
            "SGNex_H9_directRNA_replicate2_run2"
        )).split(",")
        datasets = sirv_datasets(samples, args.guides.split(","))
        root = args.output_root or PROD / "sirv4" / "_goal_opt" / "profile_sweep"
    else:
        profile = args.profile or "real-drna"
        variants = REAL_VARIANTS
        samples = (args.samples or "gencode_p00").split(",")
        datasets = [real_dataset(sample) for sample in samples]
        root = args.output_root or PROD / "gencode" / "_goal_opt" / "profile_sweep"

    root = root.resolve()
    selected = list(variants) if args.variants == "all" else args.variants.split(",")
    rows = []
    results_name = "results.tsv" if len(selected) > 1 else f"results_{selected[0]}.tsv"
    results_path = root / results_name
    for variant in selected:
        if variant not in variants:
            raise SystemExit(f"unknown variant {variant!r}")
        for dataset in datasets:
            print(f"RUN {args.domain} {variant} {dataset.sample}/{dataset.guide}", flush=True)
            row = run_one(
                dataset,
                profile,
                variant,
                variants[variant],
                root,
                score=not args.run_only,
            )
            rows.append(row)
            if args.run_only:
                print(f"  complete runtime={row['runtime_seconds']:.1f}s", flush=True)
            else:
                write_results(results_path, rows)
                print(
                    f"  honestF1={row['honest_f1_t3']:.3f} "
                    f"corrRec={row['corrected_recall_t3']:.3f} "
                    f"honPr={row['honest_precision_t3']:.3f} "
                    f"stdF1={row['tx_f1']:.3f}",
                    flush=True,
                )


if __name__ == "__main__":
    main()
