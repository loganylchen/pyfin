import re, collections
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"
S="SGNex_H9_directRNA_replicate2_run2"
NC=f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
TRK=f"gc_full/{S}__full__pyfin.tracking"
GTF=f"{B}/refs/host/gencode.v44.primary_assembly.annotation.gtf"
ENST=re.compile(r"(ENST\d+)")

# 1) expression (base ENST -> est_count)
est={}
with open(NC) as fh:
    next(fh,None)
    for ln in fh:
        p=ln.rstrip("\n").split("\t")
        if len(p)<3: continue
        m=ENST.search(p[0])
        if m:
            try: est[m.group(1)]=float(p[2])
            except ValueError: pass
expr={k for k,v in est.items() if v>=3}

# 2) pyfin matched (class '=') base ENSTs
matched=set()
for ln in open(TRK):
    c=ln.rstrip("\n").split("\t")
    if len(c)<4 or c[3]!="=": continue
    m=ENST.search(c[2])
    if m: matched.add(m.group(1))

# 3) exon count per expressed ENST from GTF
nex=collections.Counter()
tid_re=re.compile(r'transcript_id "([^".]+)')
with open(GTF) as fh:
    for ln in fh:
        if "\texon\t" not in ln: continue
        m=tid_re.search(ln)
        if m: nex[m.group(1)]+=1

def bucket_expr(v):
    if v<5: return "3-5"
    if v<10: return "5-10"
    if v<50: return "10-50"
    return "50+"
def bucket_exon(n):
    if n<=1: return "mono(1)"
    if n<=2: return "2"
    if n<=5: return "3-5"
    return "6+"

print(f"# pyfin recall diagnostic — {S} | expressed(est>=3)={len(expr)} | pyfin matched(expr)={len(expr&matched)} ({100*len(expr&matched)/len(expr):.1f}%)")
def strat(name, keyfn):
    tot=collections.Counter(); hit=collections.Counter()
    for e in expr:
        k=keyfn(e); tot[k]+=1
        if e in matched: hit[k]+=1
    print(f"\n## by {name}")
    print(f"  {'stratum':10} {'expressed':>10} {'recovered':>10} {'recall%':>8} {'missed':>8}")
    order=sorted(tot, key=lambda k: (len(k),k))
    for k in order:
        r=100*hit[k]/tot[k] if tot[k] else 0
        print(f"  {k:10} {tot[k]:>10} {hit[k]:>10} {r:>7.1f}% {tot[k]-hit[k]:>8}")
strat("expression", lambda e: bucket_expr(est.get(e,0)))
strat("exon count", lambda e: bucket_exon(nex.get(e,0)))
# cross: mono vs multi within low expr
print("\n## missed profile (the ~32% gap)")
missed=expr-matched
mono_missed=sum(1 for e in missed if nex.get(e,0)<=1)
low_missed=sum(1 for e in missed if est.get(e,0)<5)
print(f"  total missed        : {len(missed)}")
print(f"  missed & mono-exon  : {mono_missed} ({100*mono_missed/len(missed):.1f}%)")
print(f"  missed & est<5      : {low_missed} ({100*low_missed/len(missed):.1f}%)")
