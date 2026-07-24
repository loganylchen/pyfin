set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  n=0; for r in full p10 p00; do [ -s "prodfull_graphb/$r/pyfin.gtf" ] && n=$((n+1)); done
  squeue -u "$USER" -h -n pyfin_graphb 2>/dev/null | grep -q . || break
  [ "$n" -ge 3 ] && break
  sleep 300
done
python3 eval_honest.py prodfull_graphb gc_graphb > eval_honest_GRAPHB.txt 2>/dev/null
echo GRAPHB_EVAL_DONE
python3 - <<'PY'
import os
S="SGNex_H9_directRNA_replicate2_run2"
def sp(scdir,pdir,c):
    trk=f"{scdir}/{S}__{c}__pyfin.tracking"; g=f"{pdir}/{c}/pyfin.gtf"
    if not os.path.exists(trk) or not os.path.exists(g): return None
    n=sum(1 for ln in open(trk) if len(ln.split('\t'))>=4 and ln.split('\t')[3]=="=")
    tx=sum(1 for l in open(g) if "\ttranscript\t" in l); return n,tx,100*n/tx if tx else 0
def hf(f,c):
    if not os.path.exists(f): return None,None
    for ln in open(f):
        p=ln.split()
        if p and p[0]==c and len(p)>=7:
            try: return float(p[6]),float(p[5])
            except: return None,None
    return None,None
print("cond | struct_prec OFF->GRAPH->GRAPHB | corrRec OFF->GRAPH->GRAPHB | honF1 OFF->GRAPH->GRAPHB  (iso p00 honF1=42.1)")
for c in ["full","p10","p00"]:
    o=sp("gc_full","prodfull",c); g=sp("gc_graph","prodfull_graph",c); b=sp("gc_graphb","prodfull_graphb",c)
    ohf,orc=hf("eval_honest_OFF.txt",c); ghf,grc=hf("eval_honest_GRAPH.txt",c); bhf,brc=hf("eval_honest_GRAPHB.txt",c)
    op=f"{o[2]:.1f}" if o else "NA"; gp=f"{g[2]:.1f}" if g else "NA"; bp=f"{b[2]:.1f}" if b else "NA"
    print(f"{c:5} | prec {op}->{gp}->{bp} | rec {orc}->{grc}->{brc} | honF1 {ohf}->{ghf}->{bhf}")
PY
