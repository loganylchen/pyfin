set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  [ -s "prodfull_cs0/p00/pyfin.gtf" ] && break
  squeue -u "$USER" -h -n pyfin_cs0 2>/dev/null | grep -q . || break
  sleep 240
done
[ -s "prodfull_cs0/p00/pyfin.gtf" ] || { echo "cs0 no output"; exit 1; }
echo "cs0 tx: $(grep -c $'\ttranscript\t' prodfull_cs0/p00/pyfin.gtf)"
python3 eval_honest.py prodfull_cs0 gc_cs0 >/dev/null 2>&1
echo "=== SHADOW FACTORY OFF (canonical-search-bp 0) vs ON (=4) — frontier ==="
python3 frontier_of.py prodfull_cs0 gc_cs0 "cs0_OFF" 2>&1 | tail -11
echo
python3 frontier_of.py prodfull_ladder_m2_em gc_ladder_m2_em "cs4_ON" 2>&1 | tail -11
echo "=== real honest-F1 for cs0 at floors 3/4/5 ==="
python3 - <<'PY'
import re, os
HERE="/SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout"
SC=f"{HERE}/prodfull_cs0/p00/work/scores.unfiltered.tsv"; GTF=f"{HERE}/prodfull_cs0/p00/pyfin.gtf"
feat={}
with open(SC) as fh:
    h=fh.readline().rstrip("\n").split("\t"); ci={c:i for i,c in enumerate(h)}
    for ln in fh:
        c=ln.rstrip("\n").split("\t"); feat[c[ci["candidate_id"]]]=float(c[ci["abundance"]])
tid=re.compile(r'transcript_id "([^"]+)"'); lines=open(GTF).readlines()
for T in (3,4,5):
    keep={c for c,a in feat.items() if a>=T}
    d=f"{HERE}/prodfull_cs0f{T}/p00"; os.makedirs(d,exist_ok=True)
    with open(f"{d}/pyfin.gtf","w") as o:
        for ln in lines:
            if ln.startswith("#"): o.write(ln); continue
            m=tid.search(ln)
            if m and m.group(1) not in keep: continue
            o.write(ln)
PY
for T in 3 4 5; do python3 eval_honest.py prodfull_cs0f$T gc_cs0f$T 2>/dev/null | grep "^p00" | sed "s/^p00/cs0_floor$T/"; done
echo "(cols tx F1std orphan honPr corrRec honF1) refs: smart5 36.8 | bambu 39.4 | isoquant 42.1"
echo CS0_DONE
