#!/usr/bin/env python3
"""Do the read-extent siblings' TSSs REPRODUCE in independent libraries?

The end-aware contract (`gffcompare --strict-match`) costs us on human because
`.ext` models carry a read-derived 5' end that differs from the annotated TSS.
This benchmark defines the annotation as truth, so it cannot say who is right.
The strongest INTERNAL evidence available is replication: if the same start
position is used by reads of the same transcript in an independent library,
the position is reproducible structure rather than a one-library artefact.

PRE-REGISTERED (written before any result was looked at):

  Population: the `.ext` models of the discovery library that gffcompare
  matched `=` to a truth transcript AND whose 5' end is more than
  `--fail-bp` (25) from the annotated TSS -- exactly the models that fail the
  end-aware contract.

  For each, `t` = annotated TSS, `s` = the model's claimed start.

  In every replicate library, over primary MAPQ>=1 reads whose intron chain is
  the transcript's chain or a strand-aware 3' suffix of it (the same
  compatibility rule the attestation exemption uses):
      n_s = reads whose 5' end is within `--tol` (10) bp of s
      n_t = reads whose 5' end is within `--tol` bp of t
  The claimed TSS REPRODUCES in that library iff n_s >= `--min-reads` (3) and
  n_s > n_t.

  NULL CONTROL: the mirror position t' = t + (t - s), the same distance from
  the annotated TSS on the opposite side, scored by the identical rule. Read
  5' ends disperse tens of bases around a TSS, so a "cluster" that is only
  dispersion will reproduce at the mirror too.

  DECISION, fixed in advance:
    reproduction rate >> mirror rate  -> the claimed starts are reproducible
        structure; the annotation is what disagrees, and suppressing them
        would be fitting the benchmark rather than the biology.
    reproduction rate ~ mirror rate   -> they are library-specific dispersion;
        the right rule is a reproducibility requirement, not a bp floor.

Replication demonstrates a reproducible transcription/alignment position. It
does NOT demonstrate cap-defined initiation, and a shared systematic artefact
(the same aligner, the same protocol) can replicate.
"""
import argparse, collections, csv, re, sys
from pathlib import Path

BENCH = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/results"
             "/gencode_full_sweep_ratio_grid_v2")


def transcript_geometry(gtf, keep):
    ex = collections.defaultdict(list)
    meta = {}
    with open(gtf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m or m.group(1) not in keep:
                continue
            ex[m.group(1)].append((int(f[3]) - 1, int(f[4])))
            meta[m.group(1)] = (f[0], f[6])
    out = {}
    for t, e in ex.items():
        e.sort()
        c, s = meta[t]
        chain = tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1))
        out[t] = dict(chrom=c, strand=s, chain=chain,
                      tss=e[0][0] if s == "+" else e[-1][1],
                      lo=e[0][0], hi=e[-1][1])
    return out


def read_starts(bam_path, targets, pad=2000):
    """{tx: [read 5' ends]} over reads whose chain is the tx chain or its 3' suffix."""
    import pysam
    out = collections.defaultdict(list)
    bam = pysam.AlignmentFile(bam_path, "rb")
    refs = set(bam.references)
    for t, g in targets.items():
        if g["chrom"] not in refs or not g["chain"]:
            continue
        for aln in bam.fetch(g["chrom"], max(0, g["lo"] - pad), g["hi"] + pad):
            if (aln.is_unmapped or aln.is_secondary or aln.is_supplementary
                    or aln.mapping_quality < 1):
                continue
            if ("-" if aln.is_reverse else "+") != g["strand"]:
                continue
            pos, chain = aln.reference_start, []
            for op, ln in aln.cigartuples or ():
                if op in (0, 7, 8, 2):
                    pos += ln
                elif op == 3:
                    chain.append((pos, pos + ln))
                    pos += ln
            chain = tuple(chain)
            if not chain:
                continue
            n, full = len(chain), g["chain"]
            if n > len(full):
                continue
            seg = full[len(full) - n:] if g["strand"] == "+" else full[:n]
            if seg != chain:
                continue
            out[t].append(aln.reference_start if g["strand"] == "+" else aln.reference_end)
    bam.close()
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ext-tmap", required=True, help="discovery library's cmp.query.gtf.tmap")
    ap.add_argument("--ext-gtf", required=True, help="its query.gtf")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--replicates", required=True,
                    help="comma list of SGNex sample names to test in")
    ap.add_argument("--fail-bp", type=int, default=25)
    ap.add_argument("--tol", type=int, default=10)
    ap.add_argument("--min-reads", type=int, default=3)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sys.path.insert(0, "/SSD/logan/dev/pyfin")

    ext = {r["qry_id"]: r["ref_id"]
           for r in csv.DictReader(open(a.ext_tmap), delimiter="\t")
           if r["class_code"] == "=" and ".ext" in r["qry_id"]}
    truth = transcript_geometry(a.truth, set(ext.values()))
    qry = transcript_geometry(a.ext_gtf, set(ext))
    targets = {}
    for q, ref in ext.items():
        if ref not in truth or q not in qry:
            continue
        t, s = truth[ref]["tss"], qry[q]["tss"]
        if abs(s - t) <= a.fail_bp:
            continue                      # this one already passes the contract
        targets[ref] = dict(truth[ref], claimed=s, mirror=t + (t - s), qry=q)
    print(f"population: {len(targets)} `.ext` models whose 5' end is >{a.fail_bp} bp "
          f"from the annotated TSS (of {len(ext)} matched `.ext`)", flush=True)

    reps = [s for s in a.replicates.split(",") if s]
    rows = []
    per_lib = {}
    for lib in reps:
        bam = BENCH / lib / "align" / f"{lib}.sorted.bam"
        if not bam.exists():
            print(f"  {lib}: no BAM, skipped", flush=True)
            continue
        starts = read_starts(str(bam), targets)
        rep = mir = both0 = 0
        for t, g in targets.items():
            v = starts.get(t, [])
            n_s = sum(1 for x in v if abs(x - g["claimed"]) <= a.tol)
            n_t = sum(1 for x in v if abs(x - g["tss"]) <= a.tol)
            n_m = sum(1 for x in v if abs(x - g["mirror"]) <= a.tol)
            r = n_s >= a.min_reads and n_s > n_t
            m = n_m >= a.min_reads and n_m > n_t
            rep += r
            mir += m
            both0 += (len(v) == 0)
            rows.append(dict(tx=t, library=lib, n_compatible=len(v), n_claimed=n_s,
                             n_annotated=n_t, n_mirror=n_m,
                             reproduces=int(r), mirror_reproduces=int(m)))
        per_lib[lib] = (rep, mir, len(targets), both0)
        print(f"  {lib}: claimed TSS reproduces {rep}/{len(targets)} "
              f"({100*rep/max(1,len(targets)):.0f}%)  |  MIRROR control "
              f"{mir}/{len(targets)} ({100*mir/max(1,len(targets)):.0f}%)  "
              f"| no compatible read at all: {both0}", flush=True)
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    by_tx = collections.defaultdict(list)
    for r in rows:
        by_tx[r["tx"]].append(r)
    libs = [l for l in per_lib]
    n2 = sum(1 for t, v in by_tx.items() if sum(x["reproduces"] for x in v) >= 2)
    m2 = sum(1 for t, v in by_tx.items() if sum(x["mirror_reproduces"] for x in v) >= 2)
    print(f"\nreproduces in >= 2 of {len(libs)} libraries: {n2}/{len(targets)} "
          f"({100*n2/max(1,len(targets)):.0f}%)   MIRROR control: {m2}/{len(targets)} "
          f"({100*m2/max(1,len(targets)):.0f}%)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
