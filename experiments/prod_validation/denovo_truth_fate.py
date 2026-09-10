#!/usr/bin/env python3
"""De-novo (p00) fate of every expressed truth transcript, by intron chain.

On a de-novo cell every candidate id is a random novel uuid, so a truth
transcript cannot be looked up in the selection audit by name. This replays
production discovery for the intervals that contain the truth transcripts,
maps each truth intron chain to the candidate id(s) carrying it (exactly, or as
a folded shadow), and then reads that id's fate out of the audit.
"""
import argparse, collections, csv, json, re, sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from human_recall_decomposition import chains_from_gtf, read_chains

AUDIT_COLS = ("interval", "stage", "candidate_id", "source", "n_reads",
              "action", "reason", "fold_into")
ORDER = ["discovered", "canonical_gate", "junction_dominance_gate", "entered_em",
         "interval_cascade", "interval_survivor", "extent_unification",
         "global_cascade", "junction_snap", "post_refit", "post_snap", "final"]


def load_audit(out_dir):
    rows = []
    for shard in sorted((Path(out_dir) / "selection_audit").glob("*.tsv")):
        with open(shard) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == len(AUDIT_COLS):
                    rows.append(dict(zip(AUDIT_COLS, f)))
    g = Path(out_dir) / "selection_audit_global.tsv"
    if g.exists():
        with open(g) as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--bam", required=True)
    ap.add_argument("--genome", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--audit-dir", required=True)
    ap.add_argument("--chroms", default="chr1,chr19,chr21")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    sys.path.insert(0, "/SSD/logan/dev/pyfin")
    import pysam
    from fin.candidates.discovery import discover_candidates
    from fin.io.interval_manager import generate_isolated_intervals

    cfg = json.load(open(a.manifest))["config"]
    chroms = set(a.chroms.split(","))
    truth = chains_from_gtf(a.truth, chroms=chroms)
    reads = read_chains(a.bam, chroms)
    rows = load_audit(a.audit_dir)
    by_cid = collections.defaultdict(list)
    for r in rows:
        by_cid[r["candidate_id"]].append(r)
    final = {r["candidate_id"] for r in rows if r["stage"] == "final"}
    rank = {s: i for i, s in enumerate(ORDER)}

    def fate(cid, depth=0):
        recs = sorted(by_cid.get(cid, []), key=lambda r: rank.get(r["stage"], 99))
        if not recs:
            return ("no_audit_record", "", "")
        if cid in final:
            return ("emitted", "final", "")
        dec = [r for r in recs if r["action"] in ("drop", "fold")]
        if dec:
            d = dec[-1]
            if d["action"] == "fold" and d["fold_into"] and depth < 8:
                f = fate(d["fold_into"], depth + 1)
                if f[0].startswith("emitted"):
                    return ("emitted_via_fold", d["stage"], d["reason"])
            return (d["action"], d["stage"], d["reason"])
        return ("vanished_after", recs[-1]["stage"], "")

    intervals = generate_isolated_intervals(
        a.bam, gtf_path=None, max_gap=cfg["max_gap"],
        max_reads=cfg.get("max_reads"))["intervals"]
    fasta = pysam.FastaFile(a.genome)
    seqs = {}
    by_chrom = collections.defaultdict(list)
    for t, (c, s, ch) in truth.items():
        by_chrom[c].append(t)
    span = {}
    for t, (c, s, ch) in truth.items():
        if ch:
            span[t] = (c, s, ch[0][0], ch[-1][1])
    todo = collections.defaultdict(list)
    iv_by_key = {}
    for iv in intervals:
        iv_by_key[iv.region_string] = iv
        for t in by_chrom.get(iv.chrom, ()):
            if t not in span:
                continue
            c, s, lo, hi = span[t]
            if lo < iv.end and hi > iv.start and (iv.strand is None or iv.strand == s):
                todo[iv.region_string].append(t)
    mono_fold = (cfg["chain_cluster_fold_monoexon"]
                 and not (cfg["mono_resolve_post_em"] and cfg["quant_mode"] == "m2_em"))
    out = []
    seen = set()
    for key, ts in sorted(todo.items()):
        iv = iv_by_key[key]
        if iv.chrom not in seqs:
            seqs[iv.chrom] = fasta.fetch(iv.chrom).upper()
        cs = discover_candidates(
            iv, a.bam, None, seqs[iv.chrom],
            threshold=cfg["three_prime_threshold"], min_novel_reads=cfg["min_novel_reads"],
            chain_cluster=cfg["chain_cluster_discovery"], clustering=cfg["clustering"],
            chain_cluster_wobble_bp=cfg["chain_cluster_wobble_bp"],
            chain_cluster_cassette_max_exon_bp=cfg["chain_cluster_cassette_max_exon_bp"],
            chain_cluster_fold_monoexon=mono_fold,
            chain_cluster_fold_span_guard=cfg["chain_cluster_fold_span_guard"])
        by_chain = collections.defaultdict(list)
        for c in cs.candidates:
            by_chain[tuple(c.intron_chain.introns)].append(c.candidate_id)
        shadow_of = collections.defaultdict(list)
        for container, shadows in (cs.shadows or {}).items():
            for chain_t, rids in shadows:
                shadow_of[tuple(chain_t)].append(container)
        for t in ts:
            if t in seen:
                continue
            c, s, ch = truth[t]
            exact = reads.get((c, s, ch), 0)
            cids = by_chain.get(ch, [])
            if cids:
                fs = [fate(cid) for cid in cids]
                best = next((f for f in fs if f[0].startswith("emitted")), fs[0])
                out.append(dict(tx=t, interval=key, n_introns=len(ch), exact_reads=exact,
                                where="candidate", fate=best[0], stage=best[1], reason=best[2]))
            elif ch in shadow_of:
                fs = [fate(cid) for cid in shadow_of[ch]]
                best = next((f for f in fs if f[0].startswith("emitted")), fs[0])
                out.append(dict(tx=t, interval=key, n_introns=len(ch), exact_reads=exact,
                                where="folded_shadow", fate=best[0], stage=best[1], reason=best[2]))
            else:
                out.append(dict(tx=t, interval=key, n_introns=len(ch), exact_reads=exact,
                                where="never_generated", fate="-", stage="-", reason="-"))
            seen.add(t)
    for t, (c, s, ch) in truth.items():
        if t not in seen:
            out.append(dict(tx=t, interval="-", n_introns=len(ch),
                            exact_reads=reads.get((c, s, ch), 0),
                            where="no_interval", fate="-", stage="-", reason="-"))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(out[0]), delimiter="\t")
        w.writeheader(); w.writerows(out)
    print(f"expressed truth transcripts: {len(out)}")
    for k, n in collections.Counter(r["where"] for r in out).most_common():
        print(f"  {k:20s} {n:6d}")
    print("\nfate of the ones that WERE generated:")
    gen = [r for r in out if r["where"] in ("candidate", "folded_shadow")]
    for k, n in collections.Counter((r["fate"], r["reason"]) for r in gen).most_common(12):
        print(f"  {k[0]:18s} {k[1]:26s} {n:6d}")
    print("\nby read evidence (exact full-chain reads):")
    for tier, f in (("0", lambda r: r["exact_reads"] == 0), ("1", lambda r: r["exact_reads"] == 1),
                    (">=2", lambda r: r["exact_reads"] >= 2)):
        g = [r for r in out if f(r)]
        em = sum(1 for r in g if str(r["fate"]).startswith("emitted"))
        ng = sum(1 for r in g if r["where"] in ("never_generated", "no_interval"))
        print(f"  exact={tier:3s} n={len(g):6d}  emitted {em:5d}  never generated {ng:5d}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
