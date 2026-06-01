"""Fake-FASTQ signal metric: resolve M1 ties using f5c eventalign match quality.

For ambiguous reads (R_best < threshold), run f5c eventalign with fake FASTQ
(candidate cDNA at mapped region) to force ALL signal to align against each
candidate. Compute match score:

  matched = events with |z| < 2 AND model_kmer != NNNNN
  not_matched = NNNNN events + bad z events
  gap_events = gap_positions × events_per_base (estimated)
  score = matched / (total_events + gap_events)

The correct candidate scores ~0.87; wrong candidates score ~0.58-0.73.
Use scores as soft weights to redistribute R proportionally.
"""
from __future__ import annotations

import logging
import math
import os
import subprocess
import tempfile
import shutil
from typing import Dict, List, Optional, Tuple

import mappy
import numpy as np
import pysam

from fin.candidates.dataclasses import TranscriptCandidate

logger = logging.getLogger(__name__)


def _run_f5c_match_score(
    read_name: str,
    candidate_seq: str,
    signal_path: str,
    f5c_path: str = "f5c",
) -> float:
    """Run f5c eventalign with fake FASTQ/BAM and compute match score.

    Returns score in [0, 1]. 0 = total fail, 1 = perfect match.
    """
    ref_len = len(candidate_seq)
    if ref_len < 20:
        return 0.0

    total_ref_pos = ref_len - 4  # 5-mer positions
    work = tempfile.mkdtemp()
    try:
        # Write fake reference
        ref_fa = os.path.join(work, "ref.fa")
        with open(ref_fa, "w") as f:
            f.write(">candidate\n" + candidate_seq + "\n")
        pysam.faidx(ref_fa)

        # Write fake FASTQ (candidate sequence as "basecalled read")
        fq = os.path.join(work, "fake.fq")
        qual = "I" * ref_len
        with open(fq, "w") as f:
            f.write("@" + read_name + "\n" + candidate_seq + "\n+\n" + qual + "\n")

        # Write fake BAM (all M)
        bp = os.path.join(work, "fake.bam")
        hdr = pysam.AlignmentHeader.from_dict({
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": "candidate", "LN": ref_len}],
        })
        with pysam.AlignmentFile(bp, "wb", header=hdr) as out:
            seg = pysam.AlignedSegment(hdr)
            seg.query_name = read_name
            seg.query_sequence = candidate_seq
            seg.flag = 0
            seg.reference_id = 0
            seg.reference_start = 0
            seg.mapping_quality = 60
            seg.cigar = [(0, ref_len)]
            seg.query_qualities = pysam.qualitystring_to_array(qual)
            out.write(seg)
        pysam.sort("-o", bp, bp)
        pysam.index(bp)

        # f5c index
        res = subprocess.run(
            [f5c_path, "index", "--slow5", signal_path, fq],
            capture_output=True, text=True, timeout=30,
        )
        if res.returncode != 0:
            return 0.0

        # f5c eventalign
        res = subprocess.run(
            [f5c_path, "eventalign", "--rna",
             "--reads", fq, "--bam", bp, "--genome", ref_fa,
             "--min-mapq", "0", "--signal-index", "--scale-events",
             "--print-read-names", "--slow5", signal_path],
            capture_output=True, text=True, timeout=60,
        )
        if res.returncode != 0:
            return 0.0

        # Parse events and compute score
        total_ev = 0
        matched = 0
        positions = set()

        for line in res.stdout.split("\n"):
            if not line or line.startswith("contig"):
                continue
            parts = line.split("\t")
            if len(parts) < 13:
                continue
            try:
                pos = int(parts[1])
                model_kmer = parts[9]
                m_stdv = float(parts[11])
                ev_mean = float(parts[6])
                m_mean = float(parts[10])
            except (ValueError, IndexError):
                continue

            total_ev += 1
            positions.add(pos)

            if model_kmer == "NNNNN" or m_stdv <= 0:
                pass  # not matched
            else:
                z = (ev_mean - m_mean) / m_stdv
                if abs(z) < 2.0:
                    matched += 1

        if total_ev == 0:
            return 0.0

        covered = len(positions)
        gap_pos = max(0, total_ref_pos - covered)
        epb = total_ev / covered if covered > 0 else 0
        gap_ev = gap_pos * epb
        denom = total_ev + gap_ev
        return matched / denom if denom > 0 else 0.0

    except Exception as e:
        logger.debug("f5c match score failed for %s: %s", read_name, e)
        return 0.0
    finally:
        shutil.rmtree(work, ignore_errors=True)


def fake_fastq_tiebreak(
    R: np.ndarray,
    read_ids: List[str],
    read_seqs: Dict[str, str],
    candidates: List[TranscriptCandidate],
    signal_path: str,
    f5c_path: str = "f5c",
    ambig_threshold: float = 0.90,
) -> np.ndarray:
    """Resolve M1 ties using fake-FASTQ f5c match score.

    For each ambiguous read:
      1. Align to each tied candidate with mappy to get (r_st, r_en)
      2. Run f5c with candidate[r_st:r_en] as fake FASTQ + all-M BAM
      3. Compute match score (matched events / (total + gap estimate))
      4. Redistribute R proportionally by match scores

    Returns updated R matrix (copy).
    """
    R_new = R.copy()
    n_reads, n_cands = R.shape
    if n_cands < 2:
        return R_new

    # Build mappy aligners
    aligners = {}
    for j, c in enumerate(candidates):
        if c.sequence and len(c.sequence) > 20:
            a = mappy.Aligner(seq=c.sequence, preset="map-ont")
            if a:
                aligners[j] = a

    n_ambiguous = 0
    n_tiebroken = 0

    for i in range(n_reads):
        best_R = float(R[i].max())
        if best_R >= ambig_threshold:
            continue

        n_ambiguous += 1
        rid = read_ids[i]
        rseq = read_seqs.get(rid, "")
        if not rseq:
            continue

        # Find tied candidates
        tied_indices = np.where(R[i] > 0.1 * best_R)[0]
        if len(tied_indices) < 2:
            continue

        # Get mappy hits
        cand_hits = {}
        for j in tied_indices:
            if j not in aligners:
                continue
            best_hit = None
            for hit in aligners[j].map(rseq):
                if best_hit is None or hit.mlen > best_hit.mlen:
                    best_hit = hit
            if best_hit:
                cand_hits[j] = best_hit

        if len(cand_hits) < 2:
            continue

        # Run f5c match score for each tied candidate
        scores = {}
        for j, hit in cand_hits.items():
            fake_seq = candidates[j].sequence[hit.r_st:hit.r_en]
            score = _run_f5c_match_score(rid, fake_seq, signal_path, f5c_path)
            scores[j] = score

        if not scores or max(scores.values()) <= 0:
            continue

        # Check if scores actually differ
        score_vals = list(scores.values())
        if max(score_vals) - min(score_vals) < 0.02:
            # All candidates score similarly → don't tiebreak
            continue

        # Soft redistribution by match score, preserving original row sum
        total_s = sum(scores.values())
        if total_s <= 0:
            continue

        # Preserve the original total R for tied candidates
        tied_set = set(tied_indices)
        original_tied_sum = sum(float(R[i, j]) for j in tied_set)

        for j in tied_indices:
            if j in scores:
                R_new[i, j] = original_tied_sum * (scores[j] / total_s)
            else:
                R_new[i, j] = 0.0
        # Non-tied candidates keep their original R values (unchanged)

        n_tiebroken += 1

    logger.info(
        "Fake-FASTQ tiebreak: %d ambiguous, %d tiebroken",
        n_ambiguous, n_tiebroken,
    )
    return R_new
