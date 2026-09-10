#!/usr/bin/env python3
"""Rule 4 of the attestation promotion test: do the models added by the
exemption that could have been attested by the SAME reads pay for themselves?

The failure hypothesis is that attestation credits every compatible annotated
EXTENSION of one shared 3' suffix rather than recovering independently
supported isoforms. If that were happening, the added models sitting in a
sharing group should be markedly less likely to recover a distinct truth
transcript than the ones added alone, and the false positives should
concentrate there.

WHAT THIS DOES AND DOES NOT SHOW. It measures the OUTCOME of sharing groups in
the finished assembly. It does not observe the actual supporting read IDs, so
it cannot show that shared-read ambiguity was absent or biologically resolved
-- and it compares different sets of loci rather than the same models with and
without sharing, so a bucket difference is an association, not a measured cost.
Because it only sees models that survived the whole selection cascade,
a clean result is also consistent with SELECTION ENRICHMENT: the GTF floor and
later stages may already have removed the weak members of ambiguous families.
That is a good pipeline outcome, not proof that compatibility alone identifies
expression. Read the buckets as "did the vulnerable regime cost us", not as
"the vulnerability does not exist".

It also counts new recoveries without counting baseline truths the arm lost, so
its break-even comparison is an additive diagnostic rather than a per-bucket
counterfactual F1. The end-to-end grid numbers remain the promotion evidence.

Three defects of the first version, all found by Codex, are fixed here:

* it counted query-side ``=`` labels, so several added models matching the SAME
  reference all scored as wins. It now counts DISTINCT reference transcripts
  that the baseline did not already match, which is the quantity in the
  marginal-precision argument;
* its sharing predicate required one whole chain to be a suffix of the other,
  which misses the ordinary case: ``[A,C]`` and ``[B,C]`` are both attestable
  by reads carrying ``[C]`` while neither contains the other. It now links two
  models when they share ANY common strand-aware 3' suffix;
* its grouping joined a model to the first group it matched and never merged
  the others, so ``[A,C]``, ``[B,C]``, ``[C]`` could come out as a pair plus a
  singleton. It now uses union-find, so groups are true connected components.
"""
import argparse
import collections
import csv
import re


def chains(gtf):
    ex, meta = collections.defaultdict(list), {}
    with open(gtf) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon":
                continue
            m = re.search(r'transcript_id "([^"]+)"', f[8])
            if not m:
                continue
            ex[m.group(1)].append((int(f[3]) - 1, int(f[4])))
            meta[m.group(1)] = (f[0], f[6])
    out = {}
    for t, e in ex.items():
        e.sort()
        c, s = meta[t]
        out[t] = (c, s, tuple((e[i][1], e[i + 1][0]) for i in range(len(e) - 1)))
    return out


def suffixes(chain, strand):
    """Every strand-aware contiguous 3' suffix of an intron chain.

    dRNA reads are 3'-anchored, so a read attesting a model carries one of
    these. Two models that share one are STRUCTURALLY capable of being attested
    by a single read -- necessary, not sufficient: production additionally
    applies span and TES-overshoot checks, and a nonempty shared suffix means
    only that the two chains share a terminal intron. So these are structural
    ambiguity groups, not observed read-sharing groups.
    """
    n = len(chain)
    if strand == "+":
        return {chain[n - k:] for k in range(1, n + 1)}
    return {chain[:k] for k in range(1, n + 1)}


class Union:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[rb] = ra


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--base-gtf", required=True)
    ap.add_argument("--arm-gtf", required=True)
    ap.add_argument("--arm-tmap", required=True,
                    help="the arm's tmap scored against the EXPRESSED truth")
    ap.add_argument("--base-tmap", required=True,
                    help="the baseline's tmap against the same truth, so a "
                         "reference already matched without the lever is not "
                         "counted as recovered by it")
    a = ap.parse_args()

    base, arm = chains(a.base_gtf), chains(a.arm_gtf)
    seen = set(base.values())
    added = {t: v for t, v in arm.items() if v not in seen and v[2]}

    def matched(tmap):
        return {r["ref_id"]: r["qry_id"]
                for r in csv.DictReader(open(tmap), delimiter="\t")
                if r["class_code"] == "="}

    arm_hit, base_hit = matched(a.arm_tmap), matched(a.base_tmap)
    ref_of = {r["qry_id"]: r["ref_id"]
              for r in csv.DictReader(open(a.arm_tmap), delimiter="\t")
              if r["class_code"] == "="}
    newly = set(arm_hit) - set(base_hit)      # distinct truth the lever recovered

    # link added models that share ANY strand-aware 3' suffix
    u = Union()
    by_suffix = collections.defaultdict(list)
    for t, (c, s, ch) in added.items():
        u.find(t)
        for suf in suffixes(ch, s):
            by_suffix[(c, s, suf)].append(t)
    for members in by_suffix.values():
        for other in members[1:]:
            u.union(members[0], other)

    size = collections.Counter(u.find(t) for t in added)
    buckets = collections.defaultdict(lambda: [0, set()])
    for t in added:
        k = size[u.find(t)]
        b = "alone (1)" if k == 1 else ("pair (2)" if k == 2 else "group (>=3)")
        buckets[b][0] += 1
        r = ref_of.get(t)
        if r in newly:
            buckets[b][1].add(r)

    n_all = len(added)
    tp_all = len({r for v in buckets.values() for r in v[1]})
    print(f"added models: {n_all}   distinct NEWLY-matched truth: {tp_all}   "
          f"marginal precision {100*tp_all/max(1, n_all):.0f}%\n")
    print(f"{'sharing bucket':<16}{'models':>8}{'new truth':>11}{'marginal Pr':>13}")
    for b in ("alone (1)", "pair (2)", "group (>=3)"):
        n, refs = buckets[b]
        if n:
            print(f"{b:<16}{n:>8}{len(refs):>11}{100*len(refs)/n:>12.0f}%")
    print("\nThe failure hypothesis predicts the sharing buckets fall well "
          "below the alone bucket. Break-even for F1 to improve is F1_old/2, "
          "not the current precision.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
