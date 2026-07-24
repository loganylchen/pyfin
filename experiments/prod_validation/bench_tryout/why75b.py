import sys, json, pysam
sys.path.insert(0,"/SSD/logan/dev/pyfin")
from fin.io.interval_manager import is_fusion_read
from fin.candidates.intron_chains import extract_intron_chain
from collections import Counter
BAM=sys.argv[1]
cases=json.load(open("missed340.json")); b=pysam.AlignmentFile(BAM,"rb")
def rd_of(r):
    return {"query_name":r.query_name,"cigartuples":r.cigartuples,"reference_start":r.reference_start,
            "reference_end":r.reference_end,"is_supplementary":r.is_supplementary,"is_secondary":r.is_secondary}
V=Counter(); fl_softclip=[]
for case in cases:
    ch=case["chrom"]; s=case["start"]; e=case["end"]
    truth_chain=tuple((a,bb) for a,bb in case["chain"])
    # this transcript's full-length reads (exact chain + span)
    flreads=[]
    for r in b.fetch(ch,max(0,s-200),e+200):
        if r.is_secondary or r.is_supplementary or r.is_unmapped: continue
        if tuple(x for x in [(p,p2) for p,p2 in []]) : pass
        introns=[]; pos=r.reference_start
        for op,ln in r.cigartuples:
            if op in (0,2,7,8): pos+=ln
            elif op==3: introns.append((pos,pos+ln)); pos+=ln
        if tuple(introns)==truth_chain and r.reference_start<=s+20 and r.reference_end>=e-20:
            flreads.append(r)
    if not flreads: continue
    # of these full-length reads, how many are fusion-flagged (softclip>=50)?
    fusion=[r for r in flreads if is_fusion_read(rd_of(r))]
    if len(flreads)-len(fusion) < 2:   # <2 non-fusion full-length reads survive
        V["all_or_most_fulllen_reads_fusion_flagged(softclip>=50)"]+=1
        # record max softclip
        for r in flreads[:1]:
            c=r.cigartuples; sc=max(c[0][1] if c[0][0]==4 else 0, c[-1][1] if c[-1][0]==4 else 0)
            fl_softclip.append(sc)
    else:
        V["enough_nonfusion_fulllen_reads(other cause)"]+=1
b.close()
print("=== of the missed transcripts, are their full-length reads fusion-flagged? ===")
for k,v in V.most_common(): print(f"  {v:4}  {k}")
if fl_softclip:
    fl_softclip.sort()
    print(f"  sample max-softclip of fusion-flagged full-len reads: min={fl_softclip[0]} median={fl_softclip[len(fl_softclip)//2]} max={fl_softclip[-1]}")
