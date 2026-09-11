#!/usr/bin/env python3
"""Who owns the recall loss NOW, after the attestation exemption became default?

The previous audit put 62% of expressed truth on m2_support_gate. That was the
input to promoting the exemption. This re-runs the attribution against the new
default so the NEXT lever is chosen from measurement rather than from the
memory of what used to be true.

Truth is the expression contract (annotated transcripts with NanoCount TPM>0).
A truth transcript is "lost at stage X" when a candidate carrying its exact
intron chain was dropped there and no surviving model carries that chain.
"""
import argparse, collections, csv, re, sys


def chains(gtf, keep=None):
    ex, meta = collections.defaultdict(list), {}
    with open(gtf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m or (keep is not None and m.group(1) not in keep):
                continue
            ex[m.group(1)].append((int(f[3]) - 1, int(f[4])))
            meta[m.group(1)] = (f[0], f[6])
    out = {}
    for t, e in ex.items():
        e.sort()
        c, s = meta[t]
        out[t] = (c, s, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--truth", required=True)
    ap.add_argument("--assembly", required=True)
    ap.add_argument("--audit-dir", required=True)
    ap.add_argument("--tmap", required=True)
    a = ap.parse_args()

    truth = chains(a.truth)
    truth_multi = {t: v for t, v in truth.items() if v[2]}
    matched = {r["ref_id"] for r in csv.DictReader(open(a.tmap), delimiter="\t")
               if r["class_code"] == "="}
    missed = {t for t in truth_multi if t not in matched}

    # every audited action, keyed by candidate id -- annotated candidates carry
    # the transcript_id, so a missed truth is traceable by name
    acts = collections.defaultdict(list)
    import glob, os
    for path in glob.glob(os.path.join(a.audit_dir, "*")):
        if os.path.isdir(path):
            continue
        # shards are HEADERLESS; the global file carries the header. Columns:
        # interval, stage, candidate_id, source, n_reads, action, reason, fold_into
        with open(path) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 7 or f[1] == "stage":
                    continue
                acts[f[2]].append((f[5], f[6]))

    owner = collections.Counter()
    never = 0
    for t in missed:
        seq = acts.get(t)
        if not seq:
            never += 1
            continue
        drops = [r for act, r in seq if act in ("drop", "fold")]
        owner[drops[-1] if drops else "survived selection, lost later"] += 1

    n = len(missed)
    print(f"expressed multi-exon truth: {len(truth_multi)}   matched `=`: "
          f"{len(truth_multi) - n}   MISSED: {n}\n")
    print(f"{'owner of the loss':<38}{'n':>7}{'share':>8}")
    print(f"{'never a candidate (generation)':<38}{never:>7}{100*never/n:>7.1f}%")
    for k, v in owner.most_common():
        print(f"{k or '(kept, but no = match)':<38}{v:>7}{100*v/n:>7.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
