#!/usr/bin/env python3
"""Rescore already-scored cells under gffcompare's END-AWARE contract.

The project's whole result ladder uses gffcompare's default `=`, which is an
INTRON-CHAIN match: "matching transcripts are allowed to differ on the length
of the first and last exons" (gffcompare docs). Every claim about precision is
therefore blind to whether a model's TSS/TES are right. `--strict-match` makes
`=` additionally require both terminal exons to agree within `-e` (default
100 bp); models that match the chain but not the ends become `~`.

This reruns gffcompare with `--strict-match` on the SAME query and truth GTFs
the expression contract already staged (`<cells>/<sample>/<cond>/<tool>/`), so
nothing is re-assembled and the comparison is exactly paired.

CONVENTION WARNING, measured on gffcompare 0.12.6 (2026-09-10): under
`--strict-match` the `.tmap` CLASS CODE and the `.stats` TRANSCRIPT LEVEL
disagree. The class code demands EXACT terminal coordinates -- every one of the
2167 `=` survivors on a real human cell had max |terminal offset| = 0, and a
19 bp shift was demoted to `~` even at `-e 25` -- while the stats honour `-e`
(the same run reported 2557 blind, 2468 matching at `-e 100`, 2333 at `-e 25`,
with byte-identical tmap codes). THIS SCRIPT REPORTS THE STATS, i.e. the
tolerance convention. Do not attribute a stats-level loss with tmap class
codes; recompute terminal offsets from the coordinates instead.
"""
import argparse, collections, csv, re, statistics, subprocess
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

IMAGE = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"


def parse_stats(path):
    if not path.exists():
        return {}
    txt = path.read_text()
    out = {}
    m = re.search(r"Transcript level:\s*([\d.]+)\s*\|\s*([\d.]+)", txt)
    if m:
        sn, pr = float(m.group(1)), float(m.group(2))
        out.update(tx_sn=sn, tx_pr=pr,
                   tx_f1=0.0 if sn + pr == 0 else 2 * sn * pr / (sn + pr))
    for pat, key in ((r"Query mRNAs\s*:\s*(\d+)", "n_query"),
                     (r"Matching transcripts:\s*(\d+)", "n_match")):
        m = re.search(pat, txt)
        if m:
            out[key] = int(m.group(1))
    return out


def score(job):
    cell, e_bp, force, truth_override, flat = job
    cell = Path(cell).resolve()      # docker needs absolute paths
    query = cell / "query.gtf"
    if flat:
        # flat layout: <root>/cells/<cond>/<tool>, one sample, truth given
        sample, cond = "-", cell.parent.name
    else:
        # contract layout: <work>/expression_contract/cells/<sample>/<cond>/<tool>
        sample, cond = cell.parents[1].name, cell.parent.name
    truth = (Path(truth_override).resolve() if truth_override
             else cell.parents[3] / "_truth" / f"{sample}.gtf")
    if not query.exists() or not truth.exists():
        return None
    out = cell / f"strict{e_bp}"
    if force or not (out.with_suffix(".stats")).exists():
        p = subprocess.run(["docker", "run", "--rm", "-v", "/SSD:/SSD", IMAGE, "gffcompare",
                            "--strict-match", "-e", str(e_bp), "-r", str(truth),
                            "-o", str(out), str(query)], capture_output=True)
        if p.returncode != 0:
            # a silent failure would look like a missing cell, and a missing
            # cell in a paired comparison is a wrong answer, not a gap
            raise RuntimeError(f"gffcompare failed for {cell}: "
                               f"{p.stderr.decode()[:400]}")
    st = parse_stats(out.with_suffix(".stats"))
    if not st:
        return None
    return dict(sample=sample, cond=cond, tool=cell.name, **st)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cells", required=True, help="<grid>/expression_contract/cells")
    ap.add_argument("--tools", required=True, help="comma list of cell dir names to score")
    ap.add_argument("--conds", default="")
    ap.add_argument("--e-bp", type=int, default=100)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--truth", default="", help="explicit truth GTF (flat layouts)")
    ap.add_argument("--flat", action="store_true",
                    help="cells are <root>/<cond>/<tool> with no sample level")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    tools = set(a.tools.split(","))
    conds = set(c for c in a.conds.split(",") if c)
    jobs = []
    pattern = "*/*" if a.flat else "*/*/*"
    for cell in sorted(Path(a.cells).glob(pattern)):
        if cell.name not in tools or (conds and cell.parent.name not in conds):
            continue
        jobs.append((str(cell), a.e_bp, a.force, a.truth, a.flat))
    print(f"scoring {len(jobs)} cells with --strict-match -e {a.e_bp}", flush=True)
    rows = []
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        for f in as_completed([ex.submit(score, j) for j in jobs]):
            r = f.result()
            if r:
                rows.append(r)
    if len(rows) != len(jobs):
        raise SystemExit(f"scored {len(rows)} of {len(jobs)} cells -- a paired "
                         f"comparison with missing cells is not a result")
    if not rows:
        raise SystemExit("nothing scored")
    rows.sort(key=lambda r: (r["sample"], r["cond"], r["tool"]))
    with open(a.out, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    per = collections.defaultdict(list)
    for r in rows:
        per[r["tool"]].append(r)
    print(f"{'tool':14s} {'cells':>5s} {'Sn':>6s} {'Pr':>6s} {'F1':>6s}")
    for t, g in sorted(per.items()):
        m = lambda k: statistics.mean(x[k] for x in g)
        print(f"{t:14s} {len(g):5d} {m('tx_sn'):6.2f} {m('tx_pr'):6.2f} {m('tx_f1'):6.2f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
