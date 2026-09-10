#!/usr/bin/env python3
"""Discovery-level fate of truth transcripts we did not emit.

Replays production interval generation and candidate discovery (same
parameters as the run manifest of the cell) for the intervals that contain
the target transcripts, and reports for each target whether its intron chain
is a surviving candidate after discovery, a folded shadow (and of what), or
absent altogether. Everything downstream of discovery (EM, interval gates,
global selection) is NOT replayed here; a target that survives discovery and
is still not emitted died there.
"""
import argparse, collections, csv, json, re, sys
from pathlib import Path


def chains_from_gtf(path, keep):
    exons = collections.defaultdict(list)
    meta = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m or m.group(1) not in keep:
                continue
            exons[m.group(1)].append((int(f[3]) - 1, int(f[4])))
            meta[m.group(1)] = (f[0], f[6])
    out = {}
    for t, ex in exons.items():
        ex.sort()
        out[t] = (meta[t][0], meta[t][1],
                  tuple((ex[i][1], ex[i + 1][0]) for i in range(len(ex) - 1)),
                  ex[0][0], ex[-1][1])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", required=True, help="file with one transcript_id per line")
    ap.add_argument("--truth", required=True)
    ap.add_argument("--annotation", required=True, help="the --gtf the cell was run with")
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--manifest", required=True, help="run_manifest.json of the cell")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sys.path.insert(0, "/SSD/logan/dev/pyfin")
    import pysam
    from fin.candidates.discovery import discover_candidates
    from fin.io.interval_manager import generate_isolated_intervals
    from fin.io.io_gtf import GTFReader

    cfg = json.load(open(a.manifest))["config"]
    targets = [t for t in open(a.targets).read().split() if t]
    chains = chains_from_gtf(a.truth, set(targets))
    by_chrom = collections.defaultdict(list)
    for t, (c, s, ch, st, en) in chains.items():
        by_chrom[c].append(t)

    intervals = generate_isolated_intervals(
        a.bam, gtf_path=a.annotation, max_gap=cfg["max_gap"],
        max_reads=cfg.get("max_reads"))["intervals"]
    gtf_reader = GTFReader(a.annotation)
    gtf_reader.open()
    gtf_reader.parse()
    fasta = pysam.FastaFile(a.genome)
    seqs = {}
    # locate each target's interval(s): same chrom, overlapping span, same strand or unstranded
    iv_of = collections.defaultdict(list)
    for iv in intervals:
        for t in by_chrom.get(iv.chrom, ()):
            c, s, ch, st, en = chains[t]
            if st < iv.end and en > iv.start and (iv.strand is None or iv.strand == s):
                iv_of[t].append(iv)
    todo = {}
    for t, ivs in iv_of.items():
        for iv in ivs:
            todo.setdefault(iv.region_string, (iv, []))[1].append(t)

    mono_fold = (cfg["chain_cluster_fold_monoexon"]
                 and not (cfg["mono_resolve_post_em"] and cfg["quant_mode"] == "m2_em"))
    rows = []
    for key, (iv, ts) in sorted(todo.items()):
        if iv.chrom not in seqs:
            seqs[iv.chrom] = fasta.fetch(iv.chrom).upper()
        cs = discover_candidates(
            iv, a.bam, gtf_reader, seqs[iv.chrom],
            threshold=cfg["three_prime_threshold"], min_novel_reads=cfg["min_novel_reads"],
            chain_cluster=cfg["chain_cluster_discovery"], clustering=cfg["clustering"],
            chain_cluster_wobble_bp=cfg["chain_cluster_wobble_bp"],
            chain_cluster_cassette_max_exon_bp=cfg["chain_cluster_cassette_max_exon_bp"],
            chain_cluster_fold_monoexon=mono_fold,
            chain_cluster_fold_span_guard=cfg["chain_cluster_fold_span_guard"])
        by_id = {c.candidate_id: c for c in cs.candidates}
        by_chain = collections.defaultdict(list)
        for c in cs.candidates:
            by_chain[tuple(c.intron_chain.introns)].append(c)
        shadow_of = {}
        for container, shadows in (cs.shadows or {}).items():
            for chain_t, rids in shadows:
                shadow_of.setdefault(tuple(chain_t), []).append((container, len(rids)))
        for t in ts:
            c, s, ch, st, en = chains[t]
            same = by_chain.get(ch, [])
            fate, detail = "absent", ""
            if t in by_id:
                cand = by_id[t]
                fate = "candidate_self"
                detail = f"src={cand.source} reads={len(cand.supporting_read_ids)}"
            elif same:
                fate = "candidate_same_chain"
                detail = ";".join(f"{x.candidate_id}:{x.source}:{len(x.supporting_read_ids)}" for x in same)
            elif ch in shadow_of:
                fate = "folded_shadow"
                detail = ";".join(f"into {k} ({n} reads); container_src="
                                  f"{by_id[k].source if k in by_id else '?'}" for k, n in shadow_of[ch])
            else:
                # is a GTF candidate with this id in the reader's region at all?
                names = [x.transcript_id for x in gtf_reader.get_transcripts_in_region(iv.chrom, iv.start, iv.end)]
                detail = ("gtf_in_region" if t in names else "gtf_NOT_in_region") + \
                         f"; n_cands={len(cs.candidates)} n_reads={len(cs.read_ids)}"
            rows.append(dict(tx=t, interval=key, strand=s, n_introns=len(ch),
                             fate=fate, detail=detail))
    for t in targets:
        if t not in iv_of:
            rows.append(dict(tx=t, interval="-", strand=chains.get(t, ("", "?"))[1],
                             n_introns=len(chains[t][2]) if t in chains else -1,
                             fate="no_interval", detail=""))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=["tx", "interval", "strand", "n_introns", "fate", "detail"],
                           delimiter="\t")
        w.writeheader(); w.writerows(rows)
    fates = collections.Counter(r["fate"] for r in rows)
    print(f"targets={len(targets)} rows={len(rows)}")
    for k, v in fates.most_common():
        print(f"  {k:22s} {v}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
