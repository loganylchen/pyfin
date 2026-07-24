import sys,gzip,pysam
b=pysam.AlignmentFile(sys.argv[1],"rb");seen=set();n=0
with gzip.open(sys.argv[2],"wt") as fh:
 for r in b.fetch(until_eof=True):
  if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
  if r.query_name in seen: continue
  seen.add(r.query_name); s=r.get_forward_sequence() or r.query_sequence
  q=r.get_forward_qualities() or r.query_qualities
  if s is None: continue
  qs="".join(chr(c+33) for c in (q if q is not None else [30]*len(s)))
  fh.write(f"@{r.query_name}\n{s}\n+\n{qs}\n");n+=1
print(n)
