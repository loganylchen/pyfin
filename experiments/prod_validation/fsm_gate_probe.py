#!/usr/bin/env python3
"""Offline probe: FSM-style complete-splice-chain support as a survival gate.

Hypothesis under test (user-proposed): a multi-exon transcript should be
supported by at least T reads that cover its COMPLETE intron chain - i.e. the
read reaches the 5'-most splice junction - while the read's 5' terminus may
fall anywhere inside/before the first exon (5' degradation tolerated).

For every emitted transcript in the audited P0 cells this script counts, from
the raw sample BAM:

- n_fsm_exact: primary-alignment reads whose intron chain EQUALS the model
  chain exactly (same junction count, same coordinates);
- n_fsm_tol:   same but each junction may deviate by <= TOL bp (alignment
  jitter around splice sites);
- n_ism3p_tol: reads whose chain is a strict 3' SUFFIX of the model chain
  (missing >=1 of the 5'-most junctions) - degradation-shaped partial support,
  reported for contrast, never counted as FSM.

Truth labels (expression contract) are joined OFFLINE from the rescore tmap.
NanoCount/labels never feed any PyFIN run. Read-only: no fin/ code touched.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import pysam

BENCH = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/"
             "results/sirv4_full_sweep_ratio_grid_v2")
W = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid")
SAMPLES = ["SGNex_H9_directRNA_replicate3_run1",
           "SGNex_H9_directRNA_replicate4_run1"]
CONDS = ["full", "p99", "c_trunc3p50", "c_trunc5p50",
         "c_merge10", "c_ir50", "c_skip50"]
TOL = 6


def read_chains(bam_path: Path):
    """Per-chrom list of intron chains for primary alignments."""
    chains = defaultdict(list)
    with pysam.AlignmentFile(str(bam_path)) as bam:
        for aln in bam.fetch():
            if aln.is_secondary or aln.is_supplementary or aln.is_unmapped:
                continue
            pos = aln.reference_start
            chain = []
            for op, ln in aln.cigartuples:
                if op in (0, 7, 8, 2):      # M/=/X/D consume reference
                    pos += ln
                elif op == 3:               # N: intron
                    chain.append((pos, pos + ln))
                    pos += ln
            chains[aln.reference_name].append(tuple(chain))
    return chains


def model_chain(exons):
    ex = sorted(exons)
    return tuple((ex[i][1], ex[i + 1][0] - 1) for i in range(len(ex) - 1))


def close(a, b, tol):
    return abs(a[0] - b[0]) <= tol and abs(a[1] - b[1]) <= tol


def chain_match(read, model, tol):
    return len(read) == len(model) and all(
        close(r, m, tol) for r, m in zip(read, model))


def suffix_match(read, model, tol):
    """read is a strict 3'-suffix of model on + strand (5' junctions missing).
    Strand-aware: for '-' models the genomic PREFIX is the 3' side."""
    n, m = len(read), len(model)
    if n == 0 or n >= m:
        return None  # not a strict partial
    for offset in range(1, m - n + 1):
        if all(close(r, model[offset + i], tol) for i, r in enumerate(read)):
            return "genomic_suffix"
        if all(close(r, model[i], tol) for i, r in enumerate(read)):
            return "genomic_prefix"
    return None


def main() -> int:
    out_rows = []
    for sample in SAMPLES:
        bam = BENCH / sample / "align" / f"{sample}.sorted.bam"
        chains = read_chains(bam)
        for cond in CONDS:
            gtf = W / "p0_strict" / sample / cond / "out" / "assembly.gtf"
            tmap = (W / "expression_score" / "cells" / sample / cond /
                    "pyfin_p0" / "cmp.query.gtf.tmap")
            lab = {r["qry_id"]: ("TP" if r["class_code"] == "=" else "FP")
                   for r in csv.DictReader(tmap.open(), delimiter="\t")}
            models = {}
            for line in gtf.open():
                f = line.rstrip("\n").split("\t")
                if len(f) < 9 or f[2] != "exon":
                    continue
                tid = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
                models.setdefault(tid, {"chrom": f[0], "strand": f[6],
                                        "exons": [], "source": "?"})
                models[tid]["exons"].append((int(f[3]), int(f[4])))
            for line in gtf.open():
                f = line.split("\t")
                if len(f) > 3 and f[2] == "transcript":
                    tid = re.search(r'transcript_id "([^"]+)"', line).group(1)
                    src = re.search(r'transcript_source "([^"]+)"', line)
                    if tid in models and src:
                        models[tid]["source"] = src.group(1)
            for tid, m in models.items():
                mc = model_chain(m["exons"])
                if not mc:
                    continue  # mono-exon: criterion undefined
                pool = chains.get(m["chrom"], [])
                nfe = sum(1 for rc in pool if chain_match(rc, mc, 0))
                nft = sum(1 for rc in pool if chain_match(rc, mc, TOL))
                nism = 0
                for rc in pool:
                    sm = suffix_match(rc, mc, TOL)
                    if sm is None:
                        continue
                    deg_side = ("genomic_suffix" if m["strand"] == "+"
                                else "genomic_prefix")
                    if sm == deg_side:
                        nism += 1
                out_rows.append({
                    "sample": sample, "cond": cond, "tid": tid,
                    "source": m["source"], "n_junc": len(mc),
                    "label": lab.get(tid, "NA"),
                    "n_fsm_exact": nfe, "n_fsm_tol": nft,
                    "n_ism3p_tol": nism,
                })
    out = W / "fsm_gate_probe.tsv"
    with out.open("w") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(out_rows[0]), delimiter="\t")
        wr.writeheader()
        wr.writerows(out_rows)
    print(f"wrote {out} rows={len(out_rows)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
