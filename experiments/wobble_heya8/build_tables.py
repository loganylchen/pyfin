#!/usr/bin/env python3
"""Build per-sample comparison TSVs + a printable markdown summary.

For each (sample, ratio, tool):
  reads <scoring_dir>/<tool>.<query>.tmap  (exact-match '=' set)
  reads <query.gtf>                        (output transcript count)
  reads matrix/<sample>/stage/nanocount.tsv  (expressed-truth definition)

Honest metrics:
  expr1 = truth tx with nanocount est_count >= 1   (always full annotation truth)
  expr3 = truth tx with nanocount est_count >= 3
  Sn_k  = #(exact-match recovered ∩ expr_k) / |expr_k|
  Pr_k  = #(exact-match recovered ∩ expr_k) / |tool's reported transcripts|
  F1_k  = 2*Sn*Pr/(Sn+Pr)
"""
import glob, os, re, sys
from collections import defaultdict

ROOT = os.path.dirname(os.path.abspath(__file__))
MATRIX = os.path.join(ROOT, "matrix")
SAMPLES = [
    "SGNex_HEYA8_directRNA_replicate1_run1",
    "SGNex_HEYA8_directRNA_replicate1_run2",
    "SGNex_HEYA8_directRNA_replicate2_run1",
    "SGNex_HEYA8_directRNA_replicate2_run2",
    "SGNex_HEYA8_directRNA_replicate3_run1",
]
RATIOS = ["full", "p00", "c_jitter10bp", "c_skip10", "c_spurious5",
          "c_merge5", "c_flip5", "c_ir10"]
TOOLS = ["pyfin_bp10", "pyfin_bp20", "pyfin_bp20cas", "pyfin_v2", "pyfin_v3", "pyfin_prodfinal", "bambu", "espresso",
         "flair", "isoquant", "isotools", "lafite", "stringtie3", "talon"]
TRUTH = os.path.join(MATRIX, "_ref", "full", "annotation.gtf")

def parse_tx_ids(path):
    tx = set()
    if not path or not os.path.exists(path): return tx
    for line in open(path):
        if line.startswith("#"): continue
        p = line.rstrip("\n").split("\t")
        if len(p) < 9 or p[2] != "transcript": continue
        m = re.search(r'transcript_id "([^"]+)"', p[8])
        if m: tx.add(m.group(1))
    return tx

def parse_nanocount(path):
    est = {}
    if not os.path.exists(path): return est
    fh = open(path); next(fh, None)
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) >= 3:
            try: est[p[0]] = float(p[2])
            except ValueError: pass
    return est

def find_tool_gtf(sample, ratio, tool):
    if tool.startswith("pyfin_"):
        cfg = tool.split("_", 1)[1]
        p = os.path.join(MATRIX, sample, ratio, f"out_{cfg}", "assembly.gtf")
        return p if os.path.exists(p) else None
    p = os.path.join(MATRIX, "_other_tools", sample, ratio, f"{tool}.gtf")
    return p if os.path.exists(p) and os.path.getsize(p) > 0 else None

def find_tmap(sample, ratio, tool):
    # gffcompare writes <prefix>.<query>.tmap *next to the query GTF*, not in -o dir.
    gtf = find_tool_gtf(sample, ratio, tool)
    if not gtf: return None
    qdir = os.path.dirname(os.path.realpath(gtf))
    pat = os.path.join(qdir, f"{tool}.*.tmap")
    g = glob.glob(pat)
    return g[0] if g else None

def matched_refs(tmap_path):
    out = set()
    if not tmap_path or not os.path.exists(tmap_path): return out
    with open(tmap_path) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        ci = {h: i for i, h in enumerate(hdr)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[ci["class_code"]] == "=":
                out.add(f[ci["ref_id"]])
    return out

truth_tx = parse_tx_ids(TRUTH)
def f1(s, p): return 2*s*p/(s+p) if (s+p) else 0.0

all_rows = []  # for cross-sample md summary
out_root = os.path.join(ROOT, "tables")
os.makedirs(out_root, exist_ok=True)

for sample in SAMPLES:
    est = parse_nanocount(os.path.join(MATRIX, sample, "stage", "nanocount.tsv"))
    expr1 = {t for t, v in est.items() if v >= 1 and t in truth_tx}
    expr3 = {t for t, v in est.items() if v >= 3 and t in truth_tx}
    tsv = os.path.join(out_root, f"{sample}.tsv")
    with open(tsv, "w") as fh:
        fh.write("ratio\ttool\tnout\tSn1\tPr1\tF1@1\tSn3\tPr3\tF1@3\n")
        for ratio in RATIOS:
            for tool in TOOLS:
                gtf = find_tool_gtf(sample, ratio, tool)
                tmap = find_tmap(sample, ratio, tool)
                if not gtf or not tmap:
                    fh.write(f"{ratio}\t{tool}\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n")
                    continue
                nout = len(parse_tx_ids(gtf))
                m_ref = matched_refs(tmap)
                m1 = len(m_ref & expr1); m3 = len(m_ref & expr3)
                sn1 = 100*m1/max(len(expr1),1); pr1 = 100*m1/max(nout,1)
                sn3 = 100*m3/max(len(expr3),1); pr3 = 100*m3/max(nout,1)
                fh.write(f"{ratio}\t{tool}\t{nout}\t{sn1:.1f}\t{pr1:.1f}"
                         f"\t{f1(sn1,pr1):.1f}\t{sn3:.1f}\t{pr3:.1f}\t{f1(sn3,pr3):.1f}\n")
                all_rows.append((sample, ratio, tool, nout, sn1, pr1, f1(sn1,pr1),
                                 sn3, pr3, f1(sn3,pr3)))
    print(f"wrote {tsv}  (expr1={len(expr1)} expr3={len(expr3)})")

# Per-sample markdown summaries with pyfin highlighted
md = os.path.join(out_root, "SUMMARY.md")
with open(md, "w") as fh:
    fh.write("# Per-sample comparison (heya8 dual-spike, current dev pyfin + 8 tools)\n\n")
    fh.write("Honest metrics: expressed-truth denominator = nanocount est_count >= K (always vs full truth, even when input annotation is corrupted).\n\n")
    for sample in SAMPLES:
        fh.write(f"## {sample}\n\n")
        rows = [r for r in all_rows if r[0] == sample]
        for ratio in RATIOS:
            rrows = [r for r in rows if r[1] == ratio]
            if not rrows: continue
            fh.write(f"### {ratio}\n\n")
            fh.write("| tool | nout | Sn@1 | Pr@1 | F1@1 | Sn@3 | Pr@3 | F1@3 |\n")
            fh.write("|---|---|---|---|---|---|---|---|\n")
            # sort by F1@3 desc with NA last
            def key(r):
                try: return -float(r[9])
                except: return 1e9
            for r in sorted(rrows, key=key):
                _,_,tool,nout,sn1,pr1,f11,sn3,pr3,f13 = r
                marker = "**" if tool.startswith("pyfin") else ""
                fh.write(f"| {marker}{tool}{marker} | {nout} | {sn1:.1f} | {pr1:.1f} | {f11:.1f} | {sn3:.1f} | {pr3:.1f} | {f13:.1f} |\n")
            fh.write("\n")
print(f"wrote {md}")
