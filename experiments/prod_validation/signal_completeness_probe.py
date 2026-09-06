#!/usr/bin/env python3
"""Can the RAW SIGNAL say whether a read's 5' end is the MOLECULE's 5' end?

This is the one question the basecalled alignment cannot answer, and it is the
question the 5' channel needs: a read that starts upstream of an annotated TSS
is evidence about that transcript's start only if the molecule really began
there. In direct RNA the strand translocates 3'->5', so the signal is
adapter -> polyA -> transcript body, and krill reports the polyA boundary in
SAMPLES (`polya_sample_start/end`) plus the translocation rate (`read_rate`,
samples per base). Samples after the polyA, divided by the rate, estimate the
molecule's transcript length independently of what aligned.

If a read's signal-derived length is much longer than what it aligned, the
molecule extended past the alignment and its alignment start is NOT a
transcript start. This measures that gap for reads that WOULD trigger a 5'
contradiction against a correct annotation (false triggers by construction)
against reads that would not.
"""
import argparse, bisect, collections, random, statistics as st

import pysam


def annotated_tss(gtf, chroms):
    out = collections.defaultdict(list)
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript" or f[0] not in chroms:
                continue
            s, e = int(f[3]) - 1, int(f[4])
            out[(f[0], f[6])].append(s if f[6] == "+" else e)
    return {k: sorted(set(v)) for k, v in out.items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bam", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--signal", required=True)
    ap.add_argument("--chroms", default="chr1,chr19,chr21")
    ap.add_argument("--n", type=int, default=600, help="reads per group")
    ap.add_argument("--seed", type=int, default=20260906)
    a = ap.parse_args()
    import krill
    import pyslow5
    from fin.scoring.krill_aligner import make_krill_aligner

    chroms = set(a.chroms.split(","))
    tss = annotated_tss(a.gtf, chroms)
    trig, cons = [], []
    with pysam.AlignmentFile(a.bam, "rb") as bam:
        for r in bam.fetch():
            if (r.is_unmapped or r.is_secondary or r.is_supplementary
                    or r.reference_name not in chroms or not r.query_sequence):
                continue
            strand = "-" if r.is_reverse else "+"
            five = r.reference_start if strand == "+" else r.reference_end
            pos = tss.get((r.reference_name, strand))
            if not pos:
                continue
            i = bisect.bisect_left(pos, five)
            near = [pos[j] for j in (i - 1, i) if 0 <= j < len(pos)]
            if not near:
                continue
            d = min((five - t for t in near), key=abs)
            u = -d if strand == "+" else d
            rec = (r.query_name, r.query_sequence, r.query_alignment_length,
                   len(r.query_sequence), u)
            (trig if u > 50 else cons).append(rec)
    rng = random.Random(a.seed)
    rng.shuffle(trig); rng.shuffle(cons)
    groups = {"trigger": trig[:a.n], "consistent": cons[:a.n]}
    print({k: len(v) for k, v in groups.items()}, flush=True)

    aligner, gpu = make_krill_aligner(krill, "rna002", False,
                                      hmm_confidence=True, polya=True)
    s5 = pyslow5.Open(a.signal, "r")
    out = {}
    for name, recs in groups.items():
        vals = []
        for i in range(0, len(recs), 128):
            chunk = recs[i:i + 128]
            try:
                res = krill.align_reads_variants(
                    a.signal, {c[0]: {"self": c[1]} for c in chunk},
                    aligner=aligner, use_gpu=False)
            except Exception as exc:
                print("krill batch failed:", exc, flush=True)
                continue
            by = {c[0]: c for c in chunk}
            for rid, rl in res.items():
                if not rl:
                    continue
                r = rl[0]
                if r.get("status", -1) != 0 or r.get("read_rate", 0) <= 0:
                    continue
                si = r.get("start_idx")
                if si is None or len(si) == 0:
                    continue
                # DRS translocates 3'->5', so the signal index runs from the 3'
                # end upward and start_idx DECREASES along the basecalled read:
                # si[0] is the 5'-most event. Signal recorded ABOVE it is
                # molecule that kept translocating after the basecall stopped --
                # i.e. the read's 5' end is NOT the molecule's 5' end.
                rec = s5.get_read(rid, pA=False)
                if rec is None:
                    continue
                total = len(rec["signal"])
                beyond = max(total - int(si[0]), 0) / float(r["read_rate"])
                _, _, aln_len, read_len, u = by[rid]
                vals.append((beyond, aln_len, read_len, u))
        out[name] = vals
        if vals:
            beyond = [v[0] for v in vals]
            print(f"{name}: scored n={len(vals)}/{len(recs)} "
                  f"median 5'-beyond-basecall bases={st.median(beyond):.0f} "
                  f"mean={st.mean(beyond):.0f} "
                  f">50b {100*sum(1 for x in beyond if x > 50)/len(beyond):.1f}% "
                  f">200b {100*sum(1 for x in beyond if x > 200)/len(beyond):.1f}% "
                  f"| median aligned={st.median([v[1] for v in vals]):.0f}",
                  flush=True)


if __name__ == "__main__":
    main()
