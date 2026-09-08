#!/usr/bin/env python3
"""Are the count feature and the signal feature INDEPENDENT, and is there a
calibrated way to combine them?

The owner's principle: two ~90%-accurate features with independent errors
should combine to beat either alone. That holds for a majority of three, or for
two features with per-case calibrated confidence; with two hard decisions that
disagree you must follow the more accurate one, and then the pair is exactly as
good as its better member. This checks (1) independence, via the 2x2 of
count-right x signal-right on one-annotated pairs where both decide, and (2)
whether any per-case confidence separates the signal's fixes from its breaks.

Reads the rows produced by snap_signal_audit.py (same arguments).
"""
import sys, io, contextlib, collections
src=open("/SSD/logan/dev/pyfin/experiments/prod_validation/snap_signal_audit.py").read()
src=src.replace('if __name__ == "__main__":\n    raise SystemExit(main())','')
src=src.replace("    return 0\n\n\n","    return rows\n\n\n")
ns={}
sys.argv=["x","--out","/SSD/logan/dev/pyfin/experiments/prod_validation/_gencodegrid/gcsigaudit/SGNex_H9_directRNA_replicate2_run2/p00/out",
          "--gtf","/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/refs/host/gencode.v44.primary_assembly.annotation.gtf"]
with contextlib.redirect_stdout(io.StringIO()):
    exec(src, ns); rows=ns["main"]()
both=[r for r in rows if r["status"]=="one" and r["count_right"] is not None and r["sig_right"] is not None]
c=collections.Counter((r["count_right"], r["sig_right"]) for r in both); n=len(both)
print(f"one-annotated pairs where BOTH count and signal decide: {n}\n")
print(f"{'':22s} {'signal right':>13s} {'signal wrong':>13s}")
print(f"{'count right':22s} {c[(True,True)]:13d} {c[(True,False)]:13d}")
print(f"{'count wrong':22s} {c[(False,True)]:13d} {c[(False,False)]:13d}")
pc=(c[(True,True)]+c[(True,False)])/n; ps=(c[(True,True)]+c[(False,True)])/n
print(f"\ncount accuracy {100*pc:.2f}%   signal accuracy {100*ps:.2f}%")
cw=c[(False,True)]+c[(False,False)]
print(f"P(signal wrong | count wrong) = {c[(False,False)]/max(1,cw):.2f}   vs P(signal wrong) = {1-ps:.4f}")
print(f"expected joint-wrong if INDEPENDENT: {n*(1-pc)*(1-ps):.3f} pairs;  observed: {c[(False,False)]}")
print(f"\ncount's mistakes the signal would fix:  {c[(False,True)]}")
print(f"count's correct calls the signal breaks: {c[(True,False)]}")
print(f"-> signal as a second opinion, net: {c[(False,True)]-c[(True,False)]:+d} on {n} pairs")
# and the abstentions: where count's lean (the higher-count one) and the signal both speak
ab=[r for r in rows if r["status"]=="one" and r["decision"]=="abstain" and r["sig_pick"] is not None]
agree=sum(1 for r in ab if (r["sig_pick"]==(r["x"] if r["rx"]>=r["ry"] else r["y"])))
print(f"\nabstentions where the signal calls: {len(ab)}; the signal agrees with the count's LEAN on {agree}")

print("\n=== can per-case confidence separate the 5 fixes from the 12 breaks?")
def conf(r):
    n=r["vx"]+r["vy"]; return max(r["vx"],r["vy"])/n if n else 0, n
def cratio(r):
    hi,lo=max(r["rx"],r["ry"]),min(r["rx"],r["ry"]); return hi/max(lo,1)
fix=[r for r in both if not r["count_right"] and r["sig_right"]]
brk=[r for r in both if r["count_right"] and not r["sig_right"]]
for lab,sel in (("signal FIXES count (5)",fix),("signal BREAKS count (12)",brk)):
    print(f"  {lab}:")
    for r in sorted(sel,key=lambda r:-conf(r)[0]):
        f,n=conf(r); print(f"     signal supermajority {f:.2f} on {n:3d} ballots | count {r['rx']}:{r['ry']} (ratio {cratio(r):.1f})")
