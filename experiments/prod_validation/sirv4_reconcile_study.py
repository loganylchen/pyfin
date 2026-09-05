#!/usr/bin/env python3
"""Offline study: can annotation + annotation-independent discovery be
reconciled into a better assembly than either alone?

Motivation. On the SIRV4 ratio grid PyFIN is top in the annotation-poor cells
(p00-p20) and far behind in the annotation-rich ones, because it barely uses
the supplied catalog: transcript-level Sn moves only 39.7 -> 47.5 going from
"no annotation" to "the complete truth as annotation". Worse, a corrupted
catalog drags it BELOW its own de novo baseline (c_trunc3p50 F1 44.71 < p00
54.50) while 19 truth intron chains that were present in the corrupted
catalog were emitted at `full` and dropped at `c_trunc3p50`.

This script does NOT touch production. It replays fixed reconciliation
policies over artifacts that already exist:

  * annotation  A = <bench>/_ref/<sample>/<cond>/annotation.gtf
  * de novo     D = PyFIN's own p00 output for that sample (catalog-free)
  * current     C = PyFIN's output for that cell

and scores every policy with the same gffcompare/truth contract as the
benchmark. No policy may branch on the condition name -- a policy that knows
it is in `c_trunc3p50` is not a method, it is a lookup table.
"""
from __future__ import annotations

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

BENCH = Path("/autofs/mnemosyne3_SSD/logan/NanoRNATrans/benchmark/sgnex/"
             "results/sirv4_full_sweep_ratio_grid_v2")
WORK = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_sirv4grid")
TRUTH = WORK / "_ref" / "SIRV_Set4.gtf"
GFFCOMPARE_IMAGE = "quay.io/biocontainers/gffcompare:0.12.6--h9f5acd7_0"


# ---------------------------------------------------------------- GTF model
def load_gtf(p: Path) -> dict:
    """transcript_id -> {chrom, strand, exons, gene}"""
    ex: dict = defaultdict(list)
    meta: dict = {}
    if not p.exists():
        return {}
    with open(p) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or 'transcript_id "' not in f[8]:
                continue
            tid = f[8].split('transcript_id "')[1].split('"')[0]
            gid = (f[8].split('gene_id "')[1].split('"')[0]
                   if 'gene_id "' in f[8] else tid)
            ex[tid].append((int(f[3]), int(f[4])))
            meta[tid] = (f[0], f[6], gid)
    return {t: {"chrom": meta[t][0], "strand": meta[t][1], "gene": meta[t][2],
                "exons": tuple(sorted(v))} for t, v in ex.items()}


def chain_key(m: dict):
    """Structure identity. Chromosome and strand are part of the key: two
    identical intron patterns on different contigs or strands are different
    transcripts."""
    e = m["exons"]
    if len(e) == 1:
        return (m["chrom"], m["strand"], "MONO", e[0][0], e[0][1])
    return (m["chrom"], m["strand"],
            tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))


def is_multi(m: dict) -> bool:
    return len(m["exons"]) > 1


def write_gtf(models: dict, path: Path) -> None:
    with open(path, "w") as fh:
        for tid, m in sorted(models.items(),
                             key=lambda kv: (kv[1]["chrom"], kv[1]["exons"][0][0])):
            s, e = m["exons"][0][0], m["exons"][-1][1]
            attr = f'gene_id "{m["gene"]}"; transcript_id "{tid}";'
            fh.write(f'{m["chrom"]}\trec\ttranscript\t{s}\t{e}\t.\t'
                     f'{m["strand"]}\t.\t{attr}\n')
            for a, b in m["exons"]:
                fh.write(f'{m["chrom"]}\trec\texon\t{a}\t{b}\t.\t'
                         f'{m["strand"]}\t.\t{attr}\n')


# ------------------------------------------------------------------ policies
def pol_anno_only(A, D, C):
    return dict(A)


def pol_denovo_only(A, D, C):
    return dict(D)


def pol_union_dedup(A, D, C):
    """Catalog, plus any de novo structure the catalog does not contain."""
    out = dict(A)
    seen = {chain_key(m) for m in A.values()}
    for t, m in D.items():
        k = chain_key(m)
        if k not in seen:
            out[f"dn_{t}"] = m
            seen.add(k)
    return out


def pol_chain_reconcile(A, D, C):
    """Union, but when catalog and reads agree on structure the READ-derived
    endpoints win. Rationale: dRNA endpoints are observed, catalog endpoints
    may be truncated/jittered."""
    dn_by_key = {}
    for t, m in D.items():
        dn_by_key.setdefault(chain_key(m), (t, m))
    out, seen = {}, set()
    for t, m in A.items():
        k = chain_key(m)
        if k in dn_by_key:
            out[t] = dict(dn_by_key[k][1], gene=m["gene"])
        else:
            out[t] = m
        seen.add(k)
    for t, m in D.items():
        k = chain_key(m)
        if k not in seen:
            out[f"dn_{t}"] = m
            seen.add(k)
    return out


def pol_evidence_gated(A, D, C):
    """Keep the whole catalog. Repair a catalog entry ONLY where the reads
    independently produced the same intron chain (positive evidence); never
    delete a catalog entry, and add de novo structures the catalog lacks."""
    dn_by_key = {}
    for t, m in D.items():
        dn_by_key.setdefault(chain_key(m), (t, m))
    out, seen = {}, set()
    for t, m in A.items():
        k = chain_key(m)
        seen.add(k)
        if k in dn_by_key and is_multi(m):
            rd = dn_by_key[k][1]
            out[t] = {**m, "exons": rd["exons"]}
        else:
            out[t] = m
    for t, m in D.items():
        k = chain_key(m)
        if k not in seen:
            out[f"dn_{t}"] = m
            seen.add(k)
    return out


def pol_current(A, D, C):
    return dict(C)


def pol_current_plus_catalog(A, D, C):
    """What PyFIN emits today, plus catalog entries it dropped."""
    out = dict(C)
    seen = {chain_key(m) for m in C.values()}
    for t, m in A.items():
        k = chain_key(m)
        if k not in seen:
            out[f"cat_{t}"] = m
            seen.add(k)
    return out


def pol_oracle(A, D, C):
    return load_gtf(TRUTH)


POLICIES = {
    "anno_only": pol_anno_only,
    "denovo_only": pol_denovo_only,
    "union_dedup": pol_union_dedup,
    "chain_reconcile": pol_chain_reconcile,
    "evidence_gated": pol_evidence_gated,
    "current": pol_current,
    "current_plus_catalog": pol_current_plus_catalog,
    "oracle": pol_oracle,
}


# ------------------------------------------------------------------- scoring
def parse_stats(p: Path) -> dict:
    if not p.exists():
        return {}
    txt = p.read_text()
    out = {}
    for label, key in (("Transcript level", "tx"), ("Intron chain level", "ic")):
        m = re.search(rf"{re.escape(label)}:\s*([\d.]+)\s*\|\s*([\d.]+)", txt)
        if m:
            sn, pr = float(m.group(1)), float(m.group(2))
            out[f"{key}_sn"], out[f"{key}_pr"] = sn, pr
            out[f"{key}_f1"] = 0.0 if sn + pr == 0 else 2 * sn * pr / (sn + pr)
    m = re.search(r"Query mRNAs\s*:\s*(\d+)", txt)
    if m:
        out["n_query"] = int(m.group(1))
    return out


def score(models: dict, tmpdir: Path, tag: str) -> dict:
    if not models:
        return {"n_query": 0, "tx_f1": 0.0, "tx_sn": 0.0, "tx_pr": 0.0}
    g = tmpdir / f"{tag}.gtf"
    write_gtf(models, g)
    out = tmpdir / tag
    rc = subprocess.run(
        ["docker", "run", "--rm", "-v", "/SSD:/SSD", GFFCOMPARE_IMAGE,
         "gffcompare", "-r", str(TRUTH), "-o", str(out), str(g)],
        capture_output=True).returncode
    if rc != 0:
        return {}
    return parse_stats(Path(str(out) + ".stats"))


def run_cell(job) -> list:
    sample, cond, base_root = job
    A = load_gtf(BENCH / "_ref" / sample / cond / "annotation.gtf")
    D = load_gtf(Path(base_root) / sample / "p00" / "out" / "assembly.gtf")
    C = load_gtf(Path(base_root) / sample / cond / "out" / "assembly.gtf")
    rows = []
    with tempfile.TemporaryDirectory(dir="/SSD/logan/dev/pyfin/experiments/"
                                         "prod_validation/_sirv4grid") as td:
        tmp = Path(td)
        for name, fn in POLICIES.items():
            try:
                models = fn(A, D, C)
            except Exception as exc:  # noqa: BLE001
                rows.append({"sample": sample, "cond": cond, "policy": name,
                             "status": f"err:{exc}"})
                continue
            st = score(models, tmp, name)
            r = {"sample": sample, "cond": cond, "policy": name,
                 "status": "ok" if st else "score_fail"}
            r.update(st)
            rows.append(r)
    return rows


COLS = ["sample", "cond", "policy", "status", "n_query",
        "tx_sn", "tx_pr", "tx_f1", "ic_sn", "ic_pr", "ic_f1"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", default=str(WORK / "base"),
                    help="root of the PyFIN baseline grid run")
    ap.add_argument("--tsv", default=str(WORK / "reconcile_study.tsv"))
    ap.add_argument("--workers", type=int, default=8)
    a = ap.parse_args()

    samples = sorted(p.name for p in BENCH.glob("SGNex_*") if p.is_dir())
    gc = BENCH / "benchmark" / "gffcompare"
    first = next(iter(sorted(gc.glob("SGNex_*"))), None)
    conds = sorted(p.name for p in first.glob("*") if p.is_dir())
    jobs = [(s, c, a.base) for s in samples for c in conds]
    print(f"cells={len(jobs)} policies={len(POLICIES)} "
          f"gffcompare_runs={len(jobs)*len(POLICIES)}", flush=True)

    rows, done = [], 0
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        futs = [ex.submit(run_cell, j) for j in jobs]
        for f in as_completed(futs):
            rows.extend(f.result())
            done += 1
            if done % 20 == 0:
                print(f"  [{done}/{len(jobs)}]", flush=True)
    with open(a.tsv, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in sorted(rows, key=lambda z: (z["sample"], z["cond"],
                                             z["policy"])):
            w.writerow(r)
    print(f"wrote {a.tsv} ({len(rows)} rows)", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
