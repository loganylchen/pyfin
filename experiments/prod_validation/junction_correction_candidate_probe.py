#!/usr/bin/env python3
"""Candidate-level evaluation of step-5a junction correction on SIRV -- BEFORE EM.

The question the owner cares about: after per-read junction correction, of the
candidate intron chains that are NEW (present after, absent before), how many are
TP (a real SIRV chain) and how many are FP? And what was lost?

For each cell (sample x condition) we run discovery twice per interval -- with
`junction_correction=None` and with the runner's `require` callback -- and diff
the multi-exon candidate chain sets. Truth = SIRV_Set4 catalog chains (structural
truth) and the per-sample NanoCount TPM>0 subset (expression truth). Both are
evaluation-only inputs. Read-level moves are cross-tabulated too (own chain
TP/FP -> corrected chain TP/FP).
"""
from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from junction_clean_probe import truth_junctions  # noqa: E402


def truth_chains(gtf_path: Path):
    """{chain tuple} of every multi-exon transcript, 0-based half-open introns."""
    exons = collections.defaultdict(list)
    for line in open(gtf_path):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] != "exon":
            continue
        tid = None
        for kv in f[8].split(";"):
            kv = kv.strip()
            if kv.startswith("transcript_id"):
                tid = kv.split('"')[1]
                break
        if tid:
            exons[tid].append((f[0], int(f[3]), int(f[4])))
    chains = set()
    for tid, ex in exons.items():
        ex.sort(key=lambda e: e[1])
        if len(ex) < 2:
            continue
        chains.add((ex[0][0], tuple((a[2], b[1] - 1) for a, b in zip(ex, ex[1:]))))
    return chains


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--cond", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--fq", required=True)
    ap.add_argument("--blow5", required=True)
    ap.add_argument("--gtf", default="")
    ap.add_argument("--genome", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--expressed", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tol", type=int, default=15)
    ap.add_argument("--flank", type=int, default=4)
    ap.add_argument("--margin", type=float, default=10.0)
    ap.add_argument("--runner-up-margin", type=float, default=10.0)
    ap.add_argument("--min-reads", type=int, default=2)
    ap.add_argument("--gpu", type=int, default=1)
    args = ap.parse_args()

    from fin.candidates.discovery import discover_candidates
    from fin.io.interval_manager import generate_isolated_intervals
    from fin.pipeline.config import PipelineConfig
    from fin.pipeline.runner import PipelineRunner

    outdir = Path(args.out)
    outdir.mkdir(parents=True, exist_ok=True)
    work = outdir / f"_work_{args.sample}_{args.cond}"
    work.mkdir(parents=True, exist_ok=True)
    cfg = PipelineConfig(
        bam_path=args.bam, genome_fasta_path=args.genome,
        fastq_path=args.fq, signal_path=args.blow5, work_dir=str(work),
        use_gpu=bool(args.gpu), krill_pore="rna002", gtf_path=(args.gtf or None),
        junction_correction_mode="require", junction_correction_tol_bp=args.tol,
        junction_correction_flank_bp=args.flank, junction_correction_margin=args.margin,
        junction_correction_runner_up_margin=args.runner_up_margin,
        junction_correction_new_chain_min_reads=args.min_reads,
    )
    runner = PipelineRunner(cfg)
    runner.setup()
    T_cat = truth_chains(Path(args.truth))
    T_exp = truth_chains(Path(args.expressed))

    intervals = generate_isolated_intervals(cfg.bam_path, gtf_path=cfg.gtf_path or None,
                                            max_gap=cfg.max_gap, max_reads=cfg.max_reads)["intervals"]
    before, after = {}, {}          # chain -> set(read_ids)   (multi-exon only)
    moves = collections.Counter()   # (own_TP?, new_TP?) read-level
    n_rechained = 0
    for iv in intervals:
        chrom_seq = runner._genome_fasta[iv.chrom] if runner._genome_fasta and iv.chrom in runner._genome_fasta else ""
        common = dict(interval=iv, bam_path=cfg.bam_path, gtf_reader=runner._gtf_reader,
                      genome_fasta=chrom_seq, threshold=cfg.three_prime_threshold,
                      min_novel_reads=cfg.min_novel_reads, chain_cluster=True,
                      chain_cluster_wobble_bp=cfg.chain_cluster_wobble_bp,
                      chain_cluster_cassette_max_exon_bp=cfg.chain_cluster_cassette_max_exon_bp,
                      chain_cluster_fold_monoexon=False, chain_cluster_fold_span_guard=True,
                      clustering="families")
        cs0 = discover_candidates(**common, junction_correction=None)
        cb = runner._make_junction_correction(iv, chrom_seq)
        captured = {}

        def wrapped(read_chains, indels_by_read, _cb=cb, _cap=captured):
            res = _cb(read_chains, indels_by_read)
            _cap["proposals"] = res.proposals
            return res

        cs1 = discover_candidates(**common, junction_correction=wrapped if cb else None,
                                  junction_correction_apply=runner._junction_correction_apply())
        for cs, store in ((cs0, before), (cs1, after)):
            for c in cs.candidates:
                ch = tuple(c.intron_chain.introns)
                if ch:
                    store.setdefault((iv.chrom, ch), set()).update(c.supporting_read_ids)
        for p in captured.get("proposals", []):
            n_rechained += 1
            moves[((iv.chrom, p.original_chain) in T_cat, (iv.chrom, p.corrected_chain) in T_cat)] += 1

    def cls(ch):
        return ("TP_cat" if ch in T_cat else "FP_cat", "TP_exp" if ch in T_exp else "FP_exp")

    new = {ch: r for ch, r in after.items() if ch not in before}
    lost = {ch: r for ch, r in before.items() if ch not in after}
    summ = dict(
        sample=args.sample, cond=args.cond, tol=args.tol, flank=args.flank, margin=args.margin,
        runner_up_margin=args.runner_up_margin, min_reads=args.min_reads,
        n_intervals=len(intervals), reads_proposing=n_rechained,
        chains_before=len(before), chains_after=len(after),
        before_TP_cat=sum(1 for c in before if c in T_cat), before_TP_exp=sum(1 for c in before if c in T_exp),
        after_TP_cat=sum(1 for c in after if c in T_cat), after_TP_exp=sum(1 for c in after if c in T_exp),
        new_chains=len(new), new_TP_cat=sum(1 for c in new if c in T_cat), new_TP_exp=sum(1 for c in new if c in T_exp),
        new_FP=sum(1 for c in new if c not in T_cat),
        new_FP_reads=sum(len(r) for c, r in new.items() if c not in T_cat),
        new_TP_reads=sum(len(r) for c, r in new.items() if c in T_cat),
        lost_chains=len(lost), lost_TP_cat=sum(1 for c in lost if c in T_cat), lost_TP_exp=sum(1 for c in lost if c in T_exp),
        lost_FP=sum(1 for c in lost if c not in T_cat),
        proposal_moves={f"own_{'TP' if a else 'FP'}->proposed_{'TP' if b else 'FP'}": n for (a, b), n in sorted(moves.items())},
        # truth chains gained/lost at the whole-cell level
        truth_exp_covered_before=sum(1 for c in T_exp if c in before),
        truth_exp_covered_after=sum(1 for c in T_exp if c in after),
        n_truth_exp=len(T_exp),
    )
    (outdir / f"{args.sample}.{args.cond}.summary.json").write_text(json.dumps(summ, indent=1))
    with open(outdir / f"{args.sample}.{args.cond}.new_chains.tsv", "w") as fh:
        fh.write("chrom\tchain\tcat\texp\tn_reads\n")
        for (chrom, ch), r in sorted(new.items()):
            fh.write(f"{chrom}\t{';'.join(f'{d}-{a}' for d, a in ch)}\t{cls((chrom, ch))[0]}\t{cls((chrom, ch))[1]}\t{len(r)}\n")
    with open(outdir / f"{args.sample}.{args.cond}.lost_chains.tsv", "w") as fh:
        fh.write("chrom\tchain\tcat\texp\tn_reads\n")
        for (chrom, ch), r in sorted(lost.items()):
            fh.write(f"{chrom}\t{';'.join(f'{d}-{a}' for d, a in ch)}\t{cls((chrom, ch))[0]}\t{cls((chrom, ch))[1]}\t{len(r)}\n")
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
