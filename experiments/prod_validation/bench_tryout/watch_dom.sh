set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  n=0; for r in full p10 p00; do [ -s "prodfull_dom/$r/pyfin.gtf" ] && n=$((n+1)); done
  q=$(squeue -u "$USER" -h -n pyfin_dom 2>/dev/null | grep -c .)
  echo "[$(date +%H:%M)] dom gtfs=$n/3 in_queue=$q"
  { [ "$n" -ge 3 ] || [ "$q" -eq 0 ]; } && break
  sleep 300
done
echo "=== dominance-filter eval $(date) ==="
python3 eval_honest.py prodfull_dom gc_dom > eval_honest_DOM.txt 2>/dev/null
echo DOM_EVAL_DONE
echo "===== honest F1: OFF -> dominance-filter (full/p10/p00) ====="
python3 - <<'PY'
def load(f):
    d={}
    for ln in open(f):
        p=ln.split()
        if len(p)>=7 and p[0] not in ("#","cond","eval_honest"):
            try: d[p[0]]=(float(p[3]),float(p[4]),float(p[5]),float(p[6]))
            except: pass
    return d
off=load("eval_honest_OFF.txt"); dom=load("eval_honest_DOM.txt")
for c in ["full","p10","p00"]:
    if c in off and c in dom:
        print(f"{c:6}  honF1 {off[c][3]:.1f}->{dom[c][3]:.1f}  honPr {off[c][1]:.1f}->{dom[c][1]:.1f}  corrRec {off[c][2]:.1f}->{dom[c][2]:.1f}  (isoquant p00 honF1=42.1)")
PY
