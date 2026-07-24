#!/usr/bin/env bash
set -uo pipefail
OUT=/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/p00val
while squeue -h -j 211278 2>/dev/null | grep -q .; do sleep 30; done
echo "job 211278 finished $(date)"
for i in $(seq 1 30); do [ -f "$OUT/pyfin_m2a.gtf" ] && break; sleep 10; done
python3 - <<'PY'
import re
def recs(path):
    tx={}; ex={}
    for ln in open(path):
        if ln.startswith('#'): continue
        f=ln.rstrip('\n').split('\t')
        if len(f)<9: continue
        m=re.search('transcript_id "([^"]+)"',f[8])
        if not m: continue
        t=m.group(1)
        if f[2]=='transcript':
            nr=re.search('num_reads "([0-9]+)"',f[8]); ab=re.search('abundance "([0-9.]+)"',f[8]); sr=re.search('transcript_source "([^"]+)"',f[8])
            tx[t]=(f[0],int(f[3]),int(f[4]),f[6],nr.group(1) if nr else '',ab.group(1) if ab else '',sr.group(1) if sr else '')
        elif f[2]=='exon': ex.setdefault(t,[]).append((int(f[3]),int(f[4])))
    return [tx[t]+(tuple(sorted(ex.get(t,[]))),) for t in tx]
from collections import Counter
b=Counter(recs('/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/p00val/pyfin_ms.gtf'))
n=Counter(recs('/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout/p00val/pyfin_m2a.gtf'))
print(f"baseline tx={sum(b.values())}  m2a tx={sum(n.values())}")
print("M2 STEP A STRUCTURALLY BYTE-IDENTICAL:", b==n)
if b!=n:
    print("  base-only:",sum((b-n).values()),"  m2a-only:",sum((n-b).values()))
    for r,c in list((b-n).items())[:3]: print("   BASE:",r)
    for r,c in list((n-b).items())[:3]: print("   M2A :",r)
PY
echo "DONE $(date)"
