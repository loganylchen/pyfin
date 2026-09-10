#!/usr/bin/env python3
"""Replay discovery for the intervals holding gate-dropped annotated candidates
and report, per candidate, the in-pipeline compatible_read_depth, an
unconstrained (no span / no indel) count, and whether the PRODUCTION veto
(jitter + novel-junction, MAPQ>=1 registry, annotation universe over every
annotated column) would fire. Cross with truth.

Kept in step with production: it passes ``collect_compatible_depths=True``
because that collection is gated on the opt-in flag, and it calls
``unjittered_gtf_columns``, not the historical ``uncontradicted_gtf_columns``.
"""
import argparse, collections, csv, json, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from human_recall_decomposition import chains_from_gtf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dropped", required=True, help="gate_drop_precision TSV")
    ap.add_argument("--annotation", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--k", type=int, default=3)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sys.path.insert(0, "/SSD/logan/dev/pyfin")
    import pysam
    from fin.candidates.discovery import discover_candidates, chain_compatible_read_depths
    from fin.io.interval_manager import generate_isolated_intervals
    from fin.io.io_gtf import GTFReader
    from fin.scoring.m2_junction_nll import unjittered_gtf_columns

    cfg = json.load(open(a.manifest))["config"]
    dropped = {r["candidate_id"]: r for r in csv.DictReader(open(a.dropped), delimiter="\t")}
    chroms = {"chr1", "chr19", "chr21"}
    anno = chains_from_gtf(a.annotation, chroms=chroms)
    targets = {cid for cid in dropped if cid in anno and anno[cid][2]}
    intervals = generate_isolated_intervals(a.bam, gtf_path=a.annotation, max_gap=cfg["max_gap"],
                                            max_reads=cfg.get("max_reads"))["intervals"]
    gtf_reader = GTFReader(a.annotation); gtf_reader.open(); gtf_reader.parse()
    fasta = pysam.FastaFile(a.genome); seqs = {}
    # locate intervals per target by the annotation span
    span = {}
    for tx in gtf_reader.transcripts.values():
        if tx.transcript_id in targets:
            span[tx.transcript_id] = (tx.chrom, tx.strand, tx.start, tx.end)
    todo = collections.defaultdict(list)
    for iv in intervals:
        for t, (c, s, st, en) in span.items():
            if c == iv.chrom and st < iv.end and en > iv.start and (iv.strand is None or iv.strand == s):
                todo[iv.region_string].append(t)
    iv_by_key = {iv.region_string: iv for iv in intervals}
    mono_fold = (cfg["chain_cluster_fold_monoexon"]
                 and not (cfg["mono_resolve_post_em"] and cfg["quant_mode"] == "m2_em"))
    rows = []
    for key, ts in sorted(todo.items()):
        iv = iv_by_key[key]
        if iv.chrom not in seqs:
            seqs[iv.chrom] = fasta.fetch(iv.chrom).upper()
        cs = discover_candidates(
            iv, a.bam, gtf_reader, seqs[iv.chrom],
            threshold=cfg["three_prime_threshold"], min_novel_reads=cfg["min_novel_reads"],
            chain_cluster=cfg["chain_cluster_discovery"], clustering=cfg["clustering"],
            chain_cluster_wobble_bp=cfg["chain_cluster_wobble_bp"],
            chain_cluster_cassette_max_exon_bp=cfg["chain_cluster_cassette_max_exon_bp"],
            chain_cluster_fold_monoexon=mono_fold,
            chain_cluster_fold_span_guard=cfg["chain_cluster_fold_span_guard"],
            collect_compatible_depths=True)   # the production gate is opt-in now
        by_id = {c.candidate_id: c for c in cs.candidates}
        # the production veto: the MAPQ>=1 registry, the jitter+novel rules,
        # and the annotation universe from every annotated column
        obs = cs.junction_support or {}
        cands = list(cs.candidates)
        universe = {j for c in cands if getattr(c, "source", "") == "gtf"
                    for j in tuple(c.intron_chain.introns)}
        unc = unjittered_gtf_columns(cands, obs, window_bp=10,
                                     novel_junction_veto=True,
                                     annotated_junctions=universe)
        unc_ids = {cands[j].candidate_id for j in unc}
        # unconstrained recount: no span checks, no indel policy, huge tolerances
        gm = collections.defaultdict(list)
        for c in cands:
            if c.source == "gtf" and c.intron_chain.introns:
                gm[tuple(c.intron_chain.introns)].append(c)
        rj = getattr(cs, "read_junctions", None) or {}
        rs = getattr(cs, "read_strands", None) or {}
        sp = getattr(cs, "read_spans", None) or {}
        reads = [(q, rs.get(q, "+"), rj[q], sp[q][0], sp[q][1]) for q in rj if q in sp]
        loose = chain_compatible_read_depths(gm, reads, {}, end_tol_bp=10**9, flank_bp=10**9) if reads else {}
        for t in ts:
            c = by_id.get(t)
            if c is None:
                rows.append(dict(candidate_id=t, interval=key, in_candidates=0, depth=-1, loose=-1,
                                 uncontradicted=-1, truth=dropped[t]["truth"], suffix_reads=dropped[t]["suffix_reads"]))
                continue
            rows.append(dict(candidate_id=t, interval=key, in_candidates=1,
                             depth=int(getattr(c, "compatible_read_depth", 0)),
                             loose=int(loose.get(t, 0)) if reads else -2,
                             uncontradicted=int(t in unc_ids), truth=dropped[t]["truth"],
                             suffix_reads=dropped[t]["suffix_reads"]))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)
    k = a.k
    n = len(rows)
    print(f"targets {len(targets)} rows {n}  (reads available for loose recount: {sum(1 for r in rows if r['loose']>=0)})")
    off = [r for r in rows if int(r["suffix_reads"]) >= k]
    print(f"offline suffix>={k}: {len(off)} (truth {sum(int(r['truth']) for r in off)})")
    dep = [r for r in off if r["depth"] >= k]
    print(f"  of which in-pipeline depth>={k}: {len(dep)} (truth {sum(int(r['truth']) for r in dep)})")
    lo = [r for r in off if r["loose"] >= k]
    print(f"  of which loose (no span/indel) depth>={k}: {len(lo)}")
    ex = [r for r in dep if r["uncontradicted"] == 1]
    print(f"  depth>={k} AND uncontradicted (= exempt): {len(ex)} (truth {sum(int(r['truth']) for r in ex)})")
    vet = [r for r in dep if r["uncontradicted"] == 0]
    print(f"  depth>={k} but VETOED by contradiction: {len(vet)} (truth {sum(int(r['truth']) for r in vet)})")
    allv = [r for r in rows if r["uncontradicted"] == 0 and r["in_candidates"] == 1]
    print(f"vetoed among all gate-dropped annotated multi-exon: {len(allv)}/{sum(1 for r in rows if r['in_candidates']==1)} (truth {sum(int(r['truth']) for r in allv)})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
