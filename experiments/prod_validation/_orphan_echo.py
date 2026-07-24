#!/usr/bin/env python3
"""Orphan / annotation-echo rate from prod_validation tmaps + nanocount.

orphan = a truth transcript the tool 'matched' (class '=') that is NOT expressed
(nanocount est_count < 1 read). A high orphan count means the tool re-emits
annotation transcripts that have zero read support => pure annotation echo.

For each dataset/tool at the FULL ratio (clean annotation), averaged over samples:
  matched_total   = # truth tx matched (class '=')
  matched_expr1   = matched & est_count>=1
  orphan          = matched & est_count<1   (echo)
  orphan_frac     = orphan / matched_total
"""
import os, glob, re
from collections import defaultdict

DS_ROOT = "/autofs/mnemosyne4_SSD/logan/dev/pyfin/experiments/prod_validation"
DATASETS = ["sirv4", "heya8", "sequin"]
TOOLS = ["pyfin_prod", "bambu", "espresso", "flair", "isoquant",
         "isotools", "lafite", "stringtie3", "talon"]
RATIO = "full"


def _resolve(path):
    """heya8 nanocount is a symlink to /SSD/logan/dev/pyfin/... (mnemosyne4's
    local mount name, absent on this host). Remap to the autofs path."""
    if os.path.exists(path):
        return path
    if os.path.islink(path):
        tgt = os.readlink(path)
        alt = tgt.replace("/SSD/logan/dev/pyfin",
                          "/autofs/mnemosyne4_SSD/logan/dev/pyfin")
        if os.path.exists(alt):
            return alt
    return path


def nanocount(path):
    e = {}
    path = _resolve(path)
    if not os.path.exists(path):
        return e
    fh = open(path)
    next(fh, None)
    for ln in fh:
        p = ln.rstrip("\n").split("\t")
        if len(p) >= 3:
            try:
                e[p[0]] = float(p[2])
            except ValueError:
                pass
    return e


def matched_refs(tmap):
    out = set()
    if not tmap or not os.path.exists(tmap):
        return None
    with open(tmap) as fh:
        h = next(fh).rstrip("\n").split("\t")
        ci = {x: i for i, x in enumerate(h)}
        for ln in fh:
            f = ln.rstrip("\n").split("\t")
            if f[ci["class_code"]] == "=":
                out.add(f[ci["ref_id"]])
    return out


print(f"{'dataset':8s} {'tool':11s} | {'matched':>7s} {'expr>=1':>7s} {'orphan':>6s} {'orphan%':>7s}")
print("-" * 60)
summary = {}
for ds in DATASETS:
    dsdir = os.path.join(DS_ROOT, ds)
    samples = sorted(os.path.basename(os.path.dirname(p))
                     for p in glob.glob(os.path.join(dsdir, "*", "stage")))
    for tool in TOOLS:
        m_tot, m_e1, orph, n = [], [], [], 0
        for s in samples:
            est = nanocount(os.path.join(dsdir, s, "stage", "nanocount.tsv"))
            sc = os.path.join(dsdir, s, RATIO, "scoring")
            g = glob.glob(os.path.join(sc, f"gc_{tool}.*.tmap"))
            if not g:
                continue
            mr = matched_refs(g[0])
            if mr is None:
                continue
            tot = len(mr)
            e1 = sum(1 for t in mr if est.get(t, 0) >= 1)
            o = tot - e1
            m_tot.append(tot); m_e1.append(e1); orph.append(o); n += 1
        if n:
            mt = sum(m_tot) / n; me = sum(m_e1) / n; mo = sum(orph) / n
            of = 100 * mo / mt if mt else 0.0
            summary[(ds, tool)] = (mt, me, mo, of)
            print(f"{ds:8s} {tool:11s} | {mt:7.1f} {me:7.1f} {mo:6.1f} {of:6.1f}%")
    print()

# write tsv
out = os.path.join(DS_ROOT, "_orphan_echo.tsv")
with open(out, "w") as fh:
    fh.write("dataset\ttool\tmatched\texpr_ge1\torphan\torphan_pct\n")
    for (ds, tool), (mt, me, mo, of) in summary.items():
        fh.write(f"{ds}\t{tool}\t{mt:.2f}\t{me:.2f}\t{mo:.2f}\t{of:.2f}\n")
print("wrote", out)
