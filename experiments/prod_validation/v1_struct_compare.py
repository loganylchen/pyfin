#!/usr/bin/env python3
"""Structural identity of two PyFIN assemblies.

A change that is not supposed to move the result must leave ``assembly.gtf``
identical in everything but the SPELLING of the identifiers minted for novel
models and novel genes (``novel_<hex>``), which are random per run.

Spelling is all that is ignored. Each transcript is reduced to (chromosome,
strand, span, exon chain with every exon's own chromosome, strand and
attributes, all transcript attributes with a minted transcript id masked); a
novel gene is identified by the multiset of its member transcripts, so
regrouping the same transcripts into different novel genes is a difference; an
annotated gene keeps its id; every gene line keeps its span and attributes.
The two files are compared as multisets of gene and transcript records.

Exit status 0 when identical, 1 when different, 2 when a file is malformed.

    v1_struct_compare.py A/assembly.gtf B/assembly.gtf
"""
import collections
import re
import sys

_ATTR = re.compile(r'(\S+) "([^"]*)"')
_MINTED = re.compile(r"^novel_[0-9a-f]+(\.ext)?$")


class Malformed(Exception):
    pass


def _mask(value):
    return _MINTED.sub(lambda m: "novel" + (m.group(1) or ""), value)


def _other(attrs):
    """Every attribute except the two identifiers, which are handled structurally."""
    return tuple(sorted((k, v) for k, v in attrs.items() if k not in ("gene_id", "transcript_id")))


def records(path):
    genes = {}                                # gene_id -> (chrom, strand, start, end, attrs)
    tx = {}                                   # transcript_id -> (chrom, strand, start, end, gene_id, attrs)
    exons = collections.defaultdict(list)     # transcript_id -> [(chrom, strand, start, end, attrs, gene_id)]
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            attrs = dict(_ATTR.findall(f[8]))
            span = (f[0], f[6], int(f[3]), int(f[4]))
            if f[2] == "gene":
                if attrs["gene_id"] in genes:
                    raise Malformed(f"{path}: duplicate gene_id {attrs['gene_id']}")
                genes[attrs["gene_id"]] = span + (_other(attrs),)
            elif f[2] == "transcript":
                tid = attrs["transcript_id"]
                if tid in tx:
                    raise Malformed(f"{path}: duplicate transcript_id {tid}")
                tx[tid] = span + (attrs["gene_id"], attrs)
            elif f[2] == "exon":
                exons[attrs["transcript_id"]].append(span + (_other(attrs), attrs["gene_id"]))
            else:
                raise Malformed(f"{path}: unexpected feature type {f[2]!r}")
    if set(exons) - set(tx):
        raise Malformed(f"{path}: exon lines without a transcript line")

    structure = {}                            # transcript_id -> hashable structure
    members = collections.defaultdict(list)   # gene_id -> [structure]
    for tid, (chrom, strand, start, end, gid, attrs) in tx.items():
        for exon in exons[tid]:
            e_gid = exon[-1]
            if e_gid != gid:
                raise Malformed(f"{path}: exon of {tid} names gene {e_gid}, transcript names {gid}")
        kept = tuple(sorted(
            (k, _mask(v) if k == "transcript_id" else v)
            for k, v in attrs.items() if k != "gene_id"
        ))
        chain = tuple(sorted(e[:5] for e in exons[tid]))
        structure[tid] = (chrom, strand, start, end, chain, kept)
        members[gid].append(structure[tid])

    def gene_key(gid):
        if _MINTED.match(gid):
            return ("novel-gene", tuple(sorted(members[gid], key=repr)))
        return ("gene", gid)

    out = collections.Counter()
    for gid in set(genes) | set(members):
        out[("G", genes.get(gid), gene_key(gid))] += 1
    for tid, (_c, _s, _b, _e, gid, _attrs) in tx.items():
        out[("T", structure[tid], gene_key(gid))] += 1
    return out


def _describe(rec):
    if rec[0] == "G":
        return f"gene {rec[1]} {rec[2][0]} {rec[2][1] if rec[2][0] == 'gene' else f'{len(rec[2][1])} member(s)'}"
    chrom, strand, start, end, chain, kept = rec[1]
    return f"transcript {chrom}:{start}-{end}({strand}) exons={len(chain)} {dict(kept)}"


def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    try:
        a, b = records(a_path), records(b_path)
    except Malformed as exc:
        print(f"MALFORMED {exc}")
        return 2
    only_a, only_b = a - b, b - a
    print(f"A: {sum(a.values())} records   B: {sum(b.values())} records   "
          f"only-A: {sum(only_a.values())}   only-B: {sum(only_b.values())}")
    if not only_a and not only_b:
        print("IDENTICAL")
        return 0
    for name, side in (("only-A", only_a), ("only-B", only_b)):
        for rec in sorted(side, key=repr)[:5]:
            print(f"  {name}: {_describe(rec)}")
    print("DIFFERENT")
    return 1


if __name__ == "__main__":
    sys.exit(main())
