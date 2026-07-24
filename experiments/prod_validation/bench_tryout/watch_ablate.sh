set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  [ -s "prodfull_ablate/p00/pyfin.gtf" ] && break
  squeue -u "$USER" -h -n pyfin_ablate 2>/dev/null | grep -q . || break
  sleep 300
done
[ -s "prodfull_ablate/p00/pyfin.gtf" ] || { echo "ABLATE FAILED / no output"; exit 1; }
echo "=== honest metrics: OFF vs ABLATE (p00) ==="
python3 eval_honest.py prodfull_ablate gc_ablate 2>/dev/null | grep -E "^cond|^p00"
grep "^p00" eval_honest_OFF.txt | sed 's/$/   [OFF]/'
echo "=== recovery of the 1387 iso-only missed transcripts ==="
python3 - <<'PY'
import os, re
B="/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex"; S="SGNex_H9_directRNA_replicate2_run2"
HERE="/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
NC=f"{B}/results/gencode_full_sweep/{S}/full/assembly/nanocount.tsv"
ENST=re.compile(r"(ENST\d+)")
def meq(trk):
    s=set()
    if not os.path.exists(trk): return s
    for ln in open(trk):
        c=ln.rstrip("\n").split("\t")
        if len(c)>=4 and c[3]=="=":
            m=ENST.search(c[2])
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
abl=meq(f"{HERE}/gc_ablate/{S}__p00__pyfin.tracking")
iso_only=(iso-pyf)&expr3
recovered=iso_only & abl
print(f"iso-only missed by OFF pyfin (expressed): {len(iso_only)}")
print(f"  now RECOVERED by ablation run:          {len(recovered)} ({100*len(recovered)/len(iso_only):.0f}%)")
print(f"ablate expressed '=' total: {len(abl&expr3)}  (OFF was {len(pyf&expr3)})")
print(f"=> corrRec {100*len(pyf&expr3)/len(expr3):.1f}% -> {100*len(abl&expr3)/len(expr3):.1f}%")
PY
echo ABLATE_EVAL_DONE
