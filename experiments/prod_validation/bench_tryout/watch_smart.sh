set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  n=0; for a in 3 5 8; do [ -s "prodfull_smart$a/p00/pyfin.gtf" ] && n=$((n+1)); done
  squeue -u "$USER" -h -n pyfin_smart 2>/dev/null | grep -q . || break
  [ "$n" -ge 3 ] && break
  sleep 300
done
echo "=== SMART sweep honest metrics (p00) ==="
printf "%-10s %8s %7s %7s %8s %7s\n" cond tx orphan honPr corrRec honF1
for a in 3 5 8; do
  python3 eval_honest.py prodfull_smart$a gc_smart$a 2>/dev/null | grep -E "^p00" | sed "s/^p00/smartA$a  /"
done
grep "^p00" eval_honest_OFF.txt | sed 's/^p00/OFF        /'
echo "(competitors p00 honF1: isoquant 42.1  bambu 39.4  stringtie3 38.7)"
echo "=== recovery of the 1387 iso-only missed transcripts ==="
python3 - <<'PY'
import os, re
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"; S="SGNex_H9_directRNA_replicate2_run2"
HERE="/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
NC=f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"; ENST=re.compile(r"(ENST\d+)")
def meq(t):
    s=set()
    if not os.path.exists(t): return s
    for ln in open(t):
        c=ln.rstrip("\n").split("\t")
        if len(c)>=4 and c[3]=="=":
            m=ENST.search(c[2]);
            if m: s.add(m.group(1))
    return s
est={}
for ln in open(NC):
    p=ln.rstrip("\n").split("\t")
    if len(p)>=3:
        m=ENST.search(p[0])
        if m:
            try: est[m.group(1)]=float(p[2])
            except: pass
expr3={k for k,v in est.items() if v>=3}
iso=meq(f"{HERE}/gc_frontier/{S}__p00__isoquant.tracking")
pyf=meq(f"{HERE}/gc_frontier/{S}__p00__pyfin.tracking")
iso_only=(iso-pyf)&expr3
for a in (3,5,8):
    sm=meq(f"{HERE}/gc_smart{a}/{S}__p00__pyfin.tracking")
    rec=iso_only & sm
    print(f"  smartA{a}: recovered {len(rec)}/{len(iso_only)} iso-only  corrRec {100*len(pyf&expr3)/len(expr3):.1f}->{100*len(sm&expr3)/len(expr3):.1f}")
PY
echo SMART_EVAL_DONE
