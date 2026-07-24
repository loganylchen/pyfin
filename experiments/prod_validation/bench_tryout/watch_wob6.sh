set -uo pipefail
cd /SSD/logan/dev/pyfin/experiments/prod_validation/bench_tryout
end=$(( $(date +%s) + 8*3600 )); sleep 120
while [ "$(date +%s)" -lt "$end" ]; do
  n=0; for r in full p10 p00; do [ -s "prodfull_wob6/$r/pyfin.gtf" ] && n=$((n+1)); done
  squeue -u "$USER" -h -n pyfin_wob6 2>/dev/null | grep -q . || break
  [ "$n" -ge 3 ] && break
  sleep 300
done
python3 eval_honest.py prodfull_wob6 gc_wob6 > eval_honest_WOB6.txt 2>/dev/null
echo WOB6_EVAL_DONE
python3 - <<'PY'
def load(f):
    d={}
    for ln in open(f):
        p=ln.split()
        if len(p)>=7 and p[0] not in ("#","cond","eval_honest"):
            try: d[p[0]]=(float(p[1]),float(p[3]),float(p[4]),float(p[5]),float(p[6]))
            except: pass
    return d
off=load("eval_honest_OFF.txt"); w=load("eval_honest_WOB6.txt")
print("cond    tx OFF->WOB  honPr OFF->WOB  corrRec OFF->WOB  honF1 OFF->WOB  (isoquant p00 honF1=42.1)")
for c in ["full","p10","p00"]:
    if c in off and c in w:
        o=off[c]; n=w[c]
        print(f"{c:6} {o[0]:.0f}->{n[0]:.0f}  {o[1]:.1f}->{n[1]:.1f}  {o[2]:.1f}->{n[2]:.1f}  {o[4]:.1f}->{n[4]:.1f}")
PY
