set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  n=0; for r in full p10 p00; do [ -s "prodfull_wob6r5/$r/pyfin.gtf" ] && n=$((n+1)); done
  squeue -u "$USER" -h -n pyfin_wr5 2>/dev/null | grep -q . || break
  [ "$n" -ge 3 ] && break
  sleep 300
done
python3 eval_honest.py prodfull_wob6r5 gc_wob6r5 > eval_honest_WOB6R5.txt 2>/dev/null
echo WOB6R5_EVAL_DONE
python3 - <<'PY'
import re,os
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"; S="SGNex_H9_directRNA_replicate2_run2"
def sp(scdir,pdir,c):
    trk=f"{scdir}/{S}__{c}__pyfin.tracking"; g=f"{pdir}/{c}/pyfin.gtf"
    if not os.path.exists(trk): return None
    n=sum(1 for ln in open(trk) if len(ln.split('\t'))>=4 and ln.split('\t')[3]=="=")
    tx=sum(1 for l in open(g) if "\ttranscript\t" in l); return n,tx,100*n/tx
def hf(f,c):
    for ln in open(f):
        p=ln.split()
        if p and p[0]==c: return float(p[6]),float(p[5])
    return None,None
print("cond | struct_prec OFF->R5 | corrRec OFF->R5 | honF1 OFF->R5  (iso p00 honF1=42.1)")
for c in ["full","p10","p00"]:
    o=sp("gc_full","prodfull",c); n=sp("gc_wob6r5","prodfull_wob6r5",c)
    ohf,orc=hf("eval_honest_OFF.txt",c); nhf,nrc=hf("eval_honest_WOB6R5.txt",c)
    if o and n: print(f"{c:5} | {o[2]:.1f}->{n[2]:.1f} | {orc}->{nrc} | {ohf}->{nhf}")
PY
