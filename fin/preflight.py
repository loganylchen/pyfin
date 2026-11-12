"""Input checks that run before any expensive work.

Every check here corresponds to an input contract in ``docs/USAGE.md`` section 1.
A violated contract that would otherwise surface hours into a run -- or never
surface, and quietly produce a result that looks like biology -- is reported
up front, with the fix. Nothing in this module changes what the pipeline
computes; it only decides whether to start.

``PreflightError`` means the run cannot produce a meaningful result.
Warnings mean it can, but the user probably did not intend what they asked for.
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

logger = logging.getLogger(__name__)

POD5_CONVERSION = "blue-crab p2s reads.pod5 -o reads.blow5      # pip install blue-crab"
KRILL_INSTALL = (
    "pip install krill --no-deps "
    "--index-url https://loganylchen.github.io/krill-dist/simple/"
)
KRILL_UPGRADE = KRILL_INSTALL.replace("pip install krill", "pip install --upgrade krill")
_SIGNAL_SUFFIXES = {".blow5": "slow5", ".slow5": "slow5", ".pod5": "pod5"}
# How many leading primary alignments are inspected. Enough to tell "these are
# the same molecules" from "these are different files"; small enough to be free.
SAMPLE_READS = 200


class PreflightError(Exception):
    """An input contract is violated in a way that makes the run meaningless."""


@dataclass
class PreflightReport:
    signal_format: str
    warnings: List[str] = field(default_factory=list)
    sampled_reads: int = 0
    spliced_reads: int = 0
    signal_hits: Optional[int] = None
    pore: str = "rna002"
    pore_evidence: str = ""

    def warn(self, message: str) -> None:
        self.warnings.append(message)
        logger.warning("preflight: %s", message)


def resolve_signal_format(signal_path: str, requested: str) -> str:
    """Return ``slow5`` or ``pod5``; ``auto`` decides from the file name.

    An explicit format that contradicts an unambiguous file name is an error:
    opening a POD5 with the SLOW5 reader fails later with a message that does
    not name the cause.
    """
    path = Path(signal_path)
    by_name = "pod5" if path.is_dir() and any(path.glob("*.pod5")) else _SIGNAL_SUFFIXES.get(path.suffix.lower())
    if requested == "auto":
        if by_name is None:
            raise PreflightError(
                f"cannot tell the signal format of '{signal_path}' from its name "
                "(expected .blow5, .slow5 or .pod5). Pass --signal-format slow5 or "
                "--signal-format pod5."
            )
        return by_name
    if by_name is not None and by_name != requested:
        raise PreflightError(
            f"--signal-format {requested} was given, but '{signal_path}' is a "
            f"{by_name.upper()} file. Pass --signal-format {by_name}, or omit the "
            "option to let it be detected."
        )
    return requested


def refuse_pod5(signal_path: str, signal_format: str) -> None:
    """POD5 cannot be scored: say so before the backend is handed the file.

    The signal backend opens the file itself and reads SLOW5/BLOW5 only. Handed
    a POD5 it does not raise -- it takes the process down -- so this is checked
    whether or not the other input checks are enabled.
    """
    if signal_format == "pod5":
        raise PreflightError(
            f"'{signal_path}' is POD5. PyFIN's signal backend reads SLOW5/BLOW5 only.\n"
            f"  Convert once (lossless):   {POD5_CONVERSION}\n"
            "  then pass the .blow5 file to --signal."
        )


def check_backend(m2_metric: str) -> None:
    """The signal backend must import whenever signal scoring is enabled."""
    if m2_metric == "off":
        return
    try:
        import krill  # noqa: F401
    except ImportError as exc:
        raise PreflightError(
            f"signal scoring is enabled (m2_metric={m2_metric}) but the krill backend "
            f"cannot be imported ({exc}).\n"
            f"  Install it:   {KRILL_INSTALL}\n"
            "  krill is NOT the package of the same name on PyPI; use the index above.\n"
            "  To run without it, pass --m2-metric off: splice junctions are then decided by\n"
            "  sequence alone, and without krill there is no poly(A) evidence from the signal."
        ) from exc


def _polya_refusal(pore: str) -> Optional[str]:
    """Why the backend will not estimate poly(A) tails for ``pore``, or None.

    Builds the aligner the run itself builds (``fin.scoring.polya``), so the
    answer comes from the backend and not from comparing version strings.
    """
    try:
        import krill

        krill.Aligner(pore=pore, use_gpu=False, hmm_confidence=True, polya=True)
    except ImportError:
        return None  # not installed: check_backend has decided whether that is fatal
    except Exception as exc:  # noqa: BLE001 - krill raises broad errors
        # A backend that is installed and fails -- importing or building the
        # aligner -- gives the run no poly(A) calls either way. The run treats
        # both as "no calls", so neither may stop it here.
        return str(exc) or type(exc).__name__
    return None


def check_polya_backend(pore: str, report: PreflightReport) -> None:
    """Warn when the installed backend has no poly(A) model for this pore.

    krill before 0.1.6 has none for RNA004 and refuses to build the aligner.
    PyFIN treats that as "no poly(A) calls" (``fin.scoring.polya``): the run
    completes, every read abstains on the 3' endpoint channel, endpoint
    refinement falls back to end sharpness, and the polyA+5' filter no-ops.
    That is worth saying before the run.
    """
    refusal = _polya_refusal(pore)
    if refusal is None:
        return
    from fin.versions import library_versions

    report.warn(
        f"the installed krill ({library_versions().get('krill')}) cannot estimate "
        f"poly(A) tails for pore '{pore}' ({refusal}). The run will complete, but "
        "without poly(A) calls: every read abstains on the 3' endpoint evidence, and "
        "where they are enabled, endpoint refinement judges 3' ends by end sharpness "
        "alone and the polyA+5' filter (off in the default profile) is skipped. "
        "krill 0.1.6 and later carry an "
        "RNA004 poly(A) model.\n"
        f"  Upgrade:   {KRILL_UPGRADE}"
    )


def _primary_records(bam):
    for rec in bam.fetch(until_eof=True):
        if not (rec.is_unmapped or rec.is_secondary or rec.is_supplementary):
            yield rec


def _bam_contigs_and_sample(bam_path: str, report: PreflightReport) -> Tuple[Dict[str, int], List[str], List[str]]:
    import pysam

    try:
        bam = pysam.AlignmentFile(bam_path, "r")  # BAM or CRAM, detected from the file
    except (OSError, ValueError) as exc:
        raise PreflightError(f"cannot open the BAM '{bam_path}': {exc}") from exc
    with bam:
        if not bam.has_index():
            raise PreflightError(
                f"the BAM '{bam_path}' has no index. PyFIN reads it by region.\n"
                f"  Fix:   samtools index {bam_path}\n"
                "  (the BAM must be coordinate-sorted first: samtools sort)"
            )
        sort_order = bam.header.to_dict().get("HD", {}).get("SO")
        if sort_order not in (None, "coordinate"):
            report.warn(f"the BAM header declares sort order '{sort_order}', expected 'coordinate'")
        lengths = dict(zip(bam.references, bam.lengths))
        used = [s.contig for s in bam.get_index_statistics() if s.mapped > 0]

        ids: List[str] = []
        sampled_contigs: Dict[str, None] = {}
        spliced = 0
        without_seq = 0
        for rec in _primary_records(bam):
            ids.append(rec.query_name)
            sampled_contigs.setdefault(rec.reference_name)
            if rec.query_sequence is None:
                without_seq += 1
            if rec.cigartuples and any(op == 3 for op, _ in rec.cigartuples):
                spliced += 1
            if len(ids) >= SAMPLE_READS:
                break
    if not ids:
        raise PreflightError(f"the BAM '{bam_path}' contains no mapped primary alignments")
    if not used:
        # An index without per-contig counts (CRAM's .crai keeps none): the
        # contigs the sampled records lie on are what can be checked.
        used = list(sampled_contigs)
    report.sampled_reads = len(ids)
    report.spliced_reads = spliced
    if without_seq == len(ids):
        raise PreflightError(
            f"the first {len(ids)} primary alignments in '{bam_path}' carry no read "
            "sequence (SEQ is '*'). PyFIN takes read sequences from the BAM; re-align "
            "keeping the sequences on primary records."
        )
    if len(ids) >= 50 and spliced == 0:
        looks_transcriptomic = sum(
            1 for name in used[:200] if name.startswith(("ENST", "NM_", "NR_", "XM_", "XR_")) or "|" in name
        ) > len(used[:200]) / 2
        report.warn(
            f"none of the first {len(ids)} primary alignments is spliced. "
            + ("The reference names look like transcripts: PyFIN needs a GENOME alignment. "
               if looks_transcriptomic else "")
            + "Direct-RNA reads should be aligned with a splice-aware preset "
            "(minimap2 -ax splice -uf -k14); an unspliced BAM yields a mono-exon result "
            "that is an input artefact, not biology."
        )
    return lengths, used, ids


def _any_read_with_signal(bam_path: str, available: set, skip: int) -> bool:
    """Whether ANY primary record after the first ``skip`` is held by the signal file.

    Reads the BAM to its end if it has to: a bounded search could only say "not
    found so far", and that is not grounds for refusing a run.
    """
    import pysam

    logger.info("preflight: no sampled read has signal; searching the whole BAM for one that does")
    with pysam.AlignmentFile(bam_path, "r") as bam:
        for n, rec in enumerate(_primary_records(bam)):
            if n >= skip and rec.query_name in available:
                return True
    return False


def _genome_contigs(genome_path: str) -> Dict[str, int]:
    fai = Path(genome_path + ".fai")
    if not fai.exists():
        import pysam

        try:
            with pysam.FastaFile(genome_path):  # writes the .fai next to the FASTA
                pass
        except (OSError, ValueError) as exc:
            raise PreflightError(
                f"the genome '{genome_path}' has no .fai index and one could not be "
                f"created ({exc}).\n  Fix:   samtools faidx {genome_path}"
            ) from exc
    contigs: Dict[str, int] = {}
    with open(fai) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                contigs[parts[0]] = int(parts[1])
    return contigs


def _examples(names: Sequence[str], n: int = 3) -> str:
    return ", ".join(list(names)[:n])


def check_genome(genome_path: str, bam_lengths: Dict[str, int], used: Sequence[str], report: PreflightReport) -> None:
    """Same contig names AND lengths: the genome is read by contig name."""
    genome = _genome_contigs(genome_path)
    missing = [c for c in used if c not in genome]
    if len(missing) == len(used):
        raise PreflightError(
            "no contig with mapped reads in the BAM exists in the genome FASTA.\n"
            f"  BAM contigs:    {_examples(used)}\n"
            f"  genome contigs: {_examples(list(genome))}\n"
            "  The BAM was aligned to a different reference, or the two use different "
            "contig naming (e.g. 'chr1' vs '1'). Use the FASTA the reads were aligned to."
        )
    wrong_length = [c for c in used if c in genome and genome[c] != bam_lengths[c]]
    if wrong_length:
        c = wrong_length[0]
        raise PreflightError(
            f"{len(wrong_length)} contig(s) have a different length in the BAM header and "
            f"in the genome FASTA (e.g. {c}: BAM {bam_lengths[c]}, genome {genome[c]}). "
            "They are different assembly builds; use the FASTA the reads were aligned to."
        )
    if missing:
        report.warn(
            f"{len(missing)} contig(s) with mapped reads are absent from the genome FASTA "
            f"(e.g. {_examples(missing)}); reads on them cannot be assembled"
        )


def check_annotation(gtf_path: Optional[str], bam_lengths: Dict[str, int], report: PreflightReport) -> None:
    """A GTF whose contig names never match the BAM annotates nothing, silently.

    The file is read until a contig of the BAM appears -- normally the first
    feature line -- so only an annotation that matches nowhere is read to the end.
    """
    if not gtf_path:
        return
    import gzip

    opener = gzip.open if gtf_path.endswith(".gz") else open
    seen: Dict[str, None] = {}
    features = set()
    matched = False
    with opener(gtf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split("\t", 3)
            if len(parts) < 3:
                continue
            seen.setdefault(parts[0])
            features.add(parts[2])
            if parts[0] in bam_lengths:
                matched = True
                if "exon" in features:
                    break
    if not seen:
        report.warn(f"the annotation '{gtf_path}' has no feature lines; running without annotation")
        return
    if not matched:
        raise PreflightError(
            "the annotation and the BAM share no contig name, so no transcript would "
            "ever be matched and every model would be reported as novel.\n"
            f"  annotation contigs: {_examples(list(seen))}\n"
            f"  BAM contigs:        {_examples(list(bam_lengths))}\n"
            "  Use an annotation for the same assembly and naming as the alignment "
            "(e.g. 'chr1' vs '1')."
        )
    if "exon" not in features:
        report.warn(
            f"no 'exon' feature in '{gtf_path}'; PyFIN builds annotated models from "
            "exon records"
        )


class _Slow5Signal:
    """A SLOW5/BLOW5 file: which reads it holds, and what every read group says about its run.

    A merged file carries one header per read group, and the groups need not
    share a chemistry -- so every group is reported, not group zero.
    """

    def __init__(self, signal_path: str):
        import pyslow5

        self._handle = pyslow5.Open(str(signal_path), "r")
        ids, _count = self._handle.get_read_ids()
        self.ids = set(ids)

    def close(self) -> None:
        self._handle.close()

    def infos(self) -> List[Dict[str, object]]:
        """Header facts of every read group in the file."""
        out = []
        for group in range(int(self._handle.get_num_read_groups())):
            header = self._handle.get_all_headers(group)
            out.append({
                "read_group": group,
                "sequencing_kit": header.get("sequencing_kit"),
                "sample_rate": header.get("sample_frequency"),
                "experiment_type": header.get("experiment_type"),
            })
        return out


def _open_signal(signal_path: str, signal_format: str) -> _Slow5Signal:
    return _Slow5Signal(signal_path)


# Sampling rate of each direct-RNA chemistry, used only when the kit is not named.
_PORE_SAMPLE_RATE = {"rna002": 3012.0, "rna004": 4000.0}


def detect_pore(info: Dict[str, object]) -> Tuple[Optional[str], str]:
    """``(pore, evidence)`` from a signal header; pore is None when it does not say.

    The kit name decides when present. The sampling rate is the fallback: the
    two chemistries are run at different rates, and a re-basecalled or
    converted file often keeps the rate after losing the kit.
    """
    kit = str(info.get("sequencing_kit") or "").lower()
    for pore in _PORE_SAMPLE_RATE:
        if pore in kit:
            return pore, f"sequencing kit '{kit}'"
    if "rna001" in kit:
        return "rna002", f"sequencing kit '{kit}' (same pore model as RNA002)"
    try:
        rate = float(info.get("sample_rate"))  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return None, "the signal header names neither a kit nor a sampling rate"
    for pore, expected in _PORE_SAMPLE_RATE.items():
        if abs(rate - expected) <= 0.01 * expected:
            return pore, f"sampling rate {rate:g} Hz"
    return None, f"unrecognised kit '{kit}' and sampling rate {rate:g} Hz"


def resolve_pore(requested: str, infos: Sequence[Dict[str, object]], report: PreflightReport) -> str:
    """Decide the pore model; ``auto`` follows the signal file's read-group headers.

    ``infos`` holds one entry per read group of the signal file. Nothing here
    guesses: a wrong pore model does not fail, it quietly degrades signal
    scoring, so every case that is not unambiguous is either refused or said
    out loud.

    * every group is a non-RNA experiment        -> refused
    * groups name different chemistries          -> ``auto`` is refused (the
      file does not say which the BAM's reads are); an explicit model is kept
      with a warning
    * some groups name a chemistry, some nothing -> the named one, with a warning
    * no group names a chemistry                 -> ``auto`` falls back to
      rna002, what every earlier release assumed, with a warning
    """
    rna = [info for info in infos if "rna" in str(info.get("experiment_type") or "rna").lower()]
    if infos and not rna:
        raise PreflightError(
            f"the signal file was recorded as experiment type "
            f"'{infos[0].get('experiment_type')}'. PyFIN assembles direct-RNA runs."
        )
    found: Dict[str, str] = {}
    silent: List[str] = []
    for info in rna:
        pore, evidence = detect_pore(info)
        if pore is None:
            silent.append(evidence)
        else:
            found.setdefault(pore, evidence)
    if len(found) > 1:
        report.pore_evidence = "; ".join(f"{p}: {e}" for p, e in found.items())
        if requested == "auto":
            raise PreflightError(
                f"the signal file holds read groups of different chemistries ({report.pore_evidence}) "
                "and does not say which the reads of this BAM belong to. One run scores every "
                "read with one pore model: pass --pore rna002 or --pore rna004, and split the "
                "data by chemistry if the BAM mixes them."
            )
        report.warn(
            f"--pore {requested} was given, and the signal file holds read groups of different "
            f"chemistries ({report.pore_evidence}); reads of the other chemistry are scored with "
            "the wrong model"
        )
        return requested
    detected = next(iter(found), None)
    report.pore_evidence = found[detected] if detected else (silent[0] if silent else "no read group")
    if detected is None:
        if requested == "auto":
            report.warn(
                f"could not tell the pore model ({report.pore_evidence}); assuming rna002. "
                "Pass --pore rna004 if this is an RNA004 run: a wrong model does not "
                "fail, it quietly degrades signal scoring."
            )
            return "rna002"
        return requested
    if silent:
        report.warn(
            f"{len(silent)} read group(s) of the signal file name neither a kit nor a known "
            f"sampling rate; their reads are scored as {requested if requested != 'auto' else detected} "
            "like the rest. Pass --pore explicitly if that is wrong."
        )
    if requested == "auto":
        return detected
    if detected != requested:
        report.warn(
            f"--pore {requested} was given, but the signal file looks like {detected} "
            f"({report.pore_evidence}). Scoring with the wrong pore model degrades the result."
        )
    return requested


def check_signal(signal_path: str, signal_format: str, read_ids: Sequence[str],
                 report: PreflightReport, bam_path: Optional[str] = None) -> List[Dict[str, object]]:
    """The signal file must hold molecules of the BAM. Returns its read-group header facts.

    A signal file may legitimately hold only part of the BAM's reads, and that
    part need not be at the start of a coordinate-sorted BAM. So finding none
    of the sampled reads is not proof of unrelated files: the whole BAM is then
    searched for one read the signal file does hold, and the run is refused
    only when there is none.
    """
    try:
        signal = _open_signal(signal_path, signal_format)
    except PreflightError:
        raise
    except Exception as exc:  # noqa: BLE001 - reader libraries raise broadly
        raise PreflightError(
            f"cannot read '{signal_path}' as {signal_format.upper()}: {exc}"
        ) from exc
    try:
        hits = sum(1 for rid in read_ids if rid in signal.ids)
        report.signal_hits = hits
        if read_ids and hits == 0:
            if not (bam_path and _any_read_with_signal(bam_path, signal.ids, len(read_ids))):
                raise PreflightError(
                    f"no read of the BAM is present in the signal file '{signal_path}' "
                    f"({len(signal.ids)} reads). The alignment and the signal come from "
                    "different runs, or the read ids were rewritten (e.g. a prefix).\n"
                    f"  BAM read id example:    {read_ids[0]}\n"
                    f"  signal read id example: {next(iter(signal.ids), '(file is empty)')}"
                )
            report.warn(
                f"none of the first {len(read_ids)} BAM reads is present in the signal file, "
                "but later reads are: the signal file covers only part of the alignment, and "
                "reads without signal are assembled from sequence alone"
            )
        elif hits < len(read_ids):
            report.warn(
                f"only {hits} of the first {len(read_ids)} BAM reads are present in the signal "
                "file; reads without signal are assembled from sequence alone"
            )
        return signal.infos()
    finally:
        signal.close()


def run_preflight(
    *,
    bam_path: str,
    genome_path: str,
    signal_path: str,
    signal_format: str,
    m2_metric: str,
    pore: str = "auto",
    gtf_path: Optional[str] = None,
    fusion: bool = False,
    fastq_path: Optional[str] = None,
) -> PreflightReport:
    """Run every check; raise ``PreflightError`` on the first fatal one."""
    report = PreflightReport(signal_format=resolve_signal_format(signal_path, signal_format))
    refuse_pod5(signal_path, report.signal_format)
    check_backend(m2_metric)
    if fusion and not fastq_path:
        report.warn(
            "--fusion without --fastq: the sequence of hard-clipped supplementary "
            "alignments cannot be recovered, so fusion junctions on such reads cannot "
            "be validated. Pass the reads with --fastq."
        )
    for label, path in (("BAM", bam_path), ("genome", genome_path), ("signal", signal_path)):
        if not os.access(path, os.R_OK):
            raise PreflightError(f"the {label} '{path}' is not readable")
    lengths, used, ids = _bam_contigs_and_sample(bam_path, report)
    check_genome(genome_path, lengths, used, report)
    check_annotation(gtf_path, lengths, report)
    infos = check_signal(signal_path, report.signal_format, ids, report, bam_path=bam_path)
    report.pore = resolve_pore(pore, infos, report)
    check_polya_backend(report.pore, report)
    logger.info(
        "preflight ok: signal=%s, pore=%s (%s), %d reads sampled (%d spliced, %s with signal), "
        "%d warning(s)",
        report.signal_format, report.pore, report.pore_evidence or "as requested",
        report.sampled_reads, report.spliced_reads, report.signal_hits, len(report.warnings),
    )
    return report
