#!/usr/bin/env python3
"""Score the full Lever-2 matrix: on vs off per cell, mean delta per (ds, ratio).

sirv/heya8: honest F1@3/Sn@3/Pr@3 (nanocount expressed-truth, est>=3).
gencode: raw transcript-level Sn/Pr vs full GENCODE (diluted; DELTA is the signal).
Always scored against the dataset's FULL truth (ratio only changed the INPUT GTF).
"""
import os, re, glob, subprocess
REPO = "/SSD/logan/dev/pyfin"
OUT = f"{REPO}/experiments/prod_validation/lev2_fulltest"
IMG = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"
DSROOT = {
    "sirv":   f"{REPO}/experiments/prod_validation/sirv4",
    "heya8":  f"{REPO}/experiments/wobble_heya8/matrix",
    "gencode": f"{REPO}/experiments/prod_validation/gencode",
}
TRUTH = {k: os.path.realpath(f"{v}/_ref/full/annotation.gtf") for k, v in DSROOT.items()}
NANO = {  # nanocount path template; None => raw scoring (gencode)
    "sirv": f"{DSROOT['sirv']}/{{s}}/stage/nanocount.tsv",
    "heya8": f"{DSROOT['heya8']}/{{s}}/stage/nanocount.tsv",
    "gencode": None,
}

def parse_tx(p):
    s = set()
    if not os.path.exists(p): return s
    for ln in open(p):
        if ln.startswith("#"): continue
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 9 and c[2] == "transcript":
            m = re.search(r'transcript_id "([^"]+)"', c[8])
            if m: s.add(m.group(1))
    return s

def nano(p):
    e = {}
    if not p or not os.path.exists(p): return e
    fh = open(p); next(fh, None)
    for ln in fh:
        c = ln.rstrip("\n").split("\t")
        if len(c) >= 3:
            try: e[c[0]] = float(c[2])
            except ValueError: pass
    return e

def gff(query, sd, truth):
    os.makedirs(sd, exist_ok=True)
    tm = glob.glob(f"{sd}/gc.*.tmap")
    if not tm:
        if not (os.path.exists(query) and os.path.getsize(query) > 0): return None, None, None
        import shutil; shutil.copyfile(query, f"{sd}/in.gtf")
        subprocess.run(["docker","run","--rm","-u",f"{os.getuid()}:{os.getgid()}",
            "-v",f"{REPO}:{REPO}","-v","/autofs/mnemosyne3_SSD:/autofs/mnemosyne3_SSD:ro",
            "-w",sd,IMG,"gffcompare","-r",truth,"-o","gc","in.gtf"], capture_output=True)
        tm = glob.glob(f"{sd}/gc.*.tmap")
    if not tm: return None, None, None
    matched = set()
    fh = open(tm[0]); h = next(fh).rstrip("\n").split("\t"); ci = {x:i for i,x in enumerate(h)}
    for ln in fh:
        f = ln.rstrip("\n").split("\t")
        if f[ci["class_code"]] == "=": matched.add(f[ci["ref_id"]])
    sn = pr = None
    st = f"{sd}/gc.stats"
    if os.path.exists(st):
        for ln in open(st):
            if ln.strip().startswith("Transcript level"):
                n = re.findall(r"[0-9]+\.[0-9]+", ln)
                if len(n) >= 2: sn, pr = float(n[0]), float(n[1])
    return matched, sn, pr

def f1(a, b): return 2*a*b/(a+b) if (a+b) else 0.0
ttx = {k: parse_tx(v) for k, v in TRUTH.items()}
rows = []  # (ds, samp, ratio, arm, nout, Sn, Pr, F1)
for ds in ("sirv", "heya8", "gencode"):
    for d in sorted(glob.glob(f"{OUT}/{ds}/SGNex_*/*/")):
        parts = d.rstrip("/").split("/"); samp, ratio = parts[-2], parts[-1]
        e3 = None
        if NANO[ds]:
            est = nano(NANO[ds].format(s=samp))
            e3 = {t for t, v in est.items() if v >= 3 and t in ttx[ds]}
        for arm in ("off", "on"):
            q = f"{d}/{arm}.gtf"
            if not os.path.exists(q): continue
            nout = len(parse_tx(q))
            mref, sn, pr = gff(q, f"{d}/score_{arm}", TRUTH[ds])
            if e3 is not None and mref is not None:
                m3 = len(mref & e3)
                Sn = 100*m3/max(len(e3), 1); Pr = 100*m3/max(nout, 1); F = f1(Sn, Pr)
            else:
                Sn, Pr, F = (sn or 0.0), (pr or 0.0), 0.0
            rows.append((ds, samp, ratio, arm, nout, Sn, Pr, F))

by = {(r[0], r[1], r[2], r[3]): r for r in rows}
md = f"{OUT}/SUMMARY.md"
with open(md, "w") as fh:
    fh.write("# Lever-2 FULL test (--novel-junction-min-reads 2), on vs off, all datasets/ratios\n\n")
    fh.write("sirv/heya8: honest Sn3/Pr3/F1@3 (nanocount expressed-truth). gencode: raw Tx_Sn/Pr.\n\n")
    fh.write("## Mean delta (on - off) by (dataset, ratio)\n\n")
    fh.write("| ds | ratio | n | dNout | dSn | dPr | dF1@3 |\n|---|---|---|---|---|---|---|\n")
    seen = {}
    for ds, samp, ratio, arm, *_ in rows:
        seen.setdefault((ds, ratio), set()).add(samp)
    for (ds, ratio), samps in sorted(seen.items()):
        ds_rows = []
        for s in samps:
            o = by.get((ds, s, ratio, "off")); n = by.get((ds, s, ratio, "on"))
            if o and n: ds_rows.append((n[4]-o[4], n[5]-o[5], n[6]-o[6], n[7]-o[7]))
        if not ds_rows: continue
        k = len(ds_rows)
        fh.write(f"| {ds} | {ratio} | {k} | {sum(x[0] for x in ds_rows)/k:+.0f} | "
                 f"{sum(x[1] for x in ds_rows)/k:+.1f} | {sum(x[2] for x in ds_rows)/k:+.1f} | "
                 f"{sum(x[3] for x in ds_rows)/k:+.1f} |\n")
    # dataset-level rollup
    fh.write("\n## Mean delta by dataset (all ratios)\n\n| ds | cells | dSn | dPr | dF1@3 |\n|---|---|---|---|---|\n")
    for ds in ("sirv", "heya8", "gencode"):
        pairs = [(by.get((ds, s, r, "off")), by.get((ds, s, r, "on")))
                 for (d, s, r, a, *_) in rows if d == ds and a == "on"]
        ds_rows = [(n[5]-o[5], n[6]-o[6], n[7]-o[7]) for o, n in pairs if o and n]
        if not ds_rows: continue
        k = len(ds_rows)
        fh.write(f"| {ds} | {k} | {sum(x[0] for x in ds_rows)/k:+.1f} | "
                 f"{sum(x[1] for x in ds_rows)/k:+.1f} | {sum(x[2] for x in ds_rows)/k:+.1f} |\n")

    # --- ABSOLUTE pyfin off/on global means + competitor ranking ---------
    def comp_means(ds):
        """Parse existing tables/SUMMARY.md 'Global mean F1@3' per-tool table."""
        path = {"sirv": f"{DSROOT['sirv']}/tables/SUMMARY.md",
                "heya8": f"{DSROOT['heya8']}/tables/SUMMARY.md"}.get(ds)
        out = {}
        if not path or not os.path.exists(path): return out
        grab = False
        for ln in open(path):
            if "Global mean F1@3" in ln: grab = True; continue
            if grab and ln.startswith("|"):
                c = [x.strip() for x in ln.strip("|\n").split("|")]
                if len(c) >= 2 and c[0] not in ("tool", "---") and not c[0].startswith("---"):
                    try: out[c[0]] = float(c[1])
                    except ValueError: pass
        return out

    fh.write("\n## pyfin-ON vs competitors (global mean honest F1@3)\n\n")
    for ds in ("sirv", "heya8"):
        offs = [r[7] for r in rows if r[0] == ds and r[3] == "off"]
        ons = [r[7] for r in rows if r[0] == ds and r[3] == "on"]
        if not offs: continue
        pyfin_off = sum(offs)/len(offs); pyfin_on = sum(ons)/len(ons)
        comp = comp_means(ds)
        rank = dict(comp); rank["pyfin_OFF(prod)"] = pyfin_off; rank["pyfin_ON(lever2)"] = pyfin_on
        # drop the table's own pyfin_prod (we recompute it as pyfin_OFF)
        rank.pop("pyfin_prod", None)
        fh.write(f"### {ds}  (pyfin OFF measured {pyfin_off:.1f} vs table pyfin_prod {comp.get('pyfin_prod','?')})\n\n")
        fh.write("| rank | tool | F1@3 |\n|---|---|---|\n")
        for i, (t, v) in enumerate(sorted(rank.items(), key=lambda kv: -kv[1]), 1):
            star = " ⬅" if t.startswith("pyfin") else ""
            fh.write(f"| {i} | {t} | {v:.1f}{star} |\n")
        fh.write("\n")
print(open(md).read())
