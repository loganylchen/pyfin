import sys, gzip, pysam
bam_path, out_path = sys.argv[1], sys.argv[2]
b = pysam.AlignmentFile(bam_path, "rb")
seen = set(); n = 0
with gzip.open(out_path, "wt") as fh:
    for r in b.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        if r.query_name in seen: continue
        seen.add(r.query_name)
        seq = r.get_forward_sequence() or r.query_sequence
        q = r.get_forward_qualities() or r.query_qualities
        if seq is None: continue
        qstr = "".join(chr(qq+33) for qq in (q if q is not None else [30]*len(seq)))
        fh.write(f"@{r.query_name}\n{seq}\n+\n{qstr}\n"); n += 1
print(f"wrote {n} reads -> {out_path}")
