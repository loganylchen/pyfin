#!/usr/bin/env python3
"""Can raw signal decide the junction pairs the count-based snap cannot?

`junction_snap` moves a novel junction to a neighbour within 6 bp when that
neighbour has >= 2 reads and >= 2x the support; otherwise it ABSTAINS and both
survive. Signal can only help the abstentions. This measures, on a cell where
the annotation is the truth proxy for "which coordinate is real" (`full`):

  * every pair of registry junctions within 6 bp (both coordinates, same strand)
  * whether snap's ratio rule decides it or abstains
  * the annotation status of the pair: BOTH annotated (a genuine tandem site --
    signal must not touch), ONE (a wobble -- signal could decide), NEITHER
  * on the ONE-annotated abstentions, whether an AGGREGATE signal vote built
    from the per-read junction_correction audit records agrees with the
    annotation more often than the count ratio does

The aggregate vote is deliberately not the failed per-read ballot: a read
counts once per pair, abstains unless it has events on both flanks, and the
pair needs a minimum electorate before it is called.

Thresholds below are fixed BEFORE the numbers are seen (see --help).
"""
import argparse, collections, csv, re


def load_annotated_junctions(gtf, chroms):
    ex = collections.defaultdict(list); meta = {}
    with open(gtf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "exon" or f[0] not in chroms:
                continue
            t = re.search(r'transcript_id "([^"]+)"', f[8]).group(1)
            ex[t].append((int(f[3]) - 1, int(f[4]))); meta[t] = (f[0], f[6])
    out = set()
    for t, e in ex.items():
        e.sort(); c, s = meta[t]
        for i in range(len(e) - 1):
            out.add((c, s, e[i][1], e[i + 1][0]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True, help="<cell>/out directory")
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--tol", type=int, default=6, help="snap tolerance (bp, each coord)")
    ap.add_argument("--min-support", type=int, default=2)
    ap.add_argument("--min-ratio", type=float, default=2.0)
    ap.add_argument("--min-electorate", type=int, default=5,
                    help="reads with a valid signal ballot before a pair is called")
    ap.add_argument("--supermajority", type=float, default=0.75)
    ap.add_argument("--vote-margin", type=float, default=0.0,
                    help="per-read |own_nll - alt_nll| below this = abstain")
    a = ap.parse_args()

    # --- the SAME support snap uses: primary BAM junction observations ------
    # (not the registry's per-interval counts; snap rescans the BAM itself)
    import json, math, sys
    sys.path.insert(0, "/SSD/logan/dev/pyfin")
    from fin.pipeline.junction_snap import collect_junction_support
    bam = json.load(open(f"{a.out}/run_manifest.json"))["config"]["bam_path"]
    support = collect_junction_support(bam) or {}
    reg = {}
    chroms = set()
    for (chrom, strand), cnt in support.items():
        for (d, acc), n in cnt.items():
            reg[(chrom, strand, d, acc)] = (n, n)
            chroms.add(chrom)
    anno = load_annotated_junctions(a.gtf, chroms)

    # --- neighbouring pairs within tol at both coordinates -----------------
    by_cs = collections.defaultdict(list)
    for k in reg:
        by_cs[(k[0], k[1])].append(k)
    pairs = []
    for cs, ks in by_cs.items():
        ks.sort(key=lambda k: k[2])
        for i, x in enumerate(ks):
            for y in ks[i + 1:]:
                if y[2] - x[2] > a.tol:
                    break
                if abs(y[3] - x[3]) <= a.tol and y != x:
                    pairs.append((x, y))

    # --- per-read signal records, keyed by unordered junction pair ---------
    ballots = collections.defaultdict(dict)   # pair -> read -> (+1 for x, -1 for y, 0)
    n_rec = n_valid = 0
    with open(f"{a.out}/junction_correction.tsv") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            n_rec += 1
            try:
                own = tuple(int(v) for v in r["own_junction"].split("-"))
                alt = tuple(int(v) for v in r["alt_junction"].split("-"))
                own_nll, alt_nll = float(r["own_nll"]), float(r["alt_nll"])
            except (ValueError, KeyError):
                continue
            # a record with no scored events on either hypothesis carries no
            # ballot; "abstain / own_margin_below_threshold" DOES carry one --
            # it only means the read did not switch under the correction margin
            if int(r.get("n_events_own") or 0) == 0 or int(r.get("n_events_alt") or 0) == 0:
                continue
            # chrom/strand come from the interval string
            m = re.match(r"([^:]+):", r["interval"])
            if not m:
                continue
            chrom, strand = m.group(1), r["strand"]
            kx = (chrom, strand) + own; ky = (chrom, strand) + alt
            key = tuple(sorted([kx, ky]))
            d = alt_nll - own_nll          # > 0: own fits better
            if abs(d) < a.vote_margin:
                vote = 0
            else:
                vote = (1 if d > 0 else -1) if kx == key[0] else (-1 if d > 0 else 1)
            # one read, one vote per pair: keep the strongest-margin record;
            # remember which coordinate the read's OWN chain carried (its origin)
            origin = kx
            prev = ballots[key].get(r["read_id"])
            if prev is None or abs(d) > abs(prev[1]):
                ballots[key][r["read_id"]] = (vote, d, origin)
            n_valid += 1

    # --- classify ----------------------------------------------------------
    C = collections.Counter()
    rows = []
    for x, y in pairs:
        rx, ry = reg[x][0], reg[y][0]
        ax, ay = x in anno, y in anno
        status = "both" if ax and ay else ("one" if ax or ay else "neither")
        hi, lo = (x, y) if rx >= ry else (y, x)
        rhi, rlo = reg[hi][0], reg[lo][0]
        # snap's rule verbatim: the target needs max(min_support, floor(cur*ratio)+1)
        # -- STRICTLY more than ratio x current, so 2 vs 4 abstains, 2 vs 5 moves
        required = max(a.min_support, math.floor(rlo * a.min_ratio) + 1)
        snap_moves = rhi >= required and rlo > 0
        decision = "moved" if snap_moves else "abstain"
        C[(decision, status)] += 1
        # count's pick = hi; annotation's pick (when exactly one) = the annotated one
        anno_pick = x if (ax and not ay) else (y if (ay and not ax) else None)
        count_right = (anno_pick == hi) if anno_pick else None
        # signal aggregate
        key = tuple(sorted([x, y]))
        b = ballots.get(key, {})
        vx = sum(1 for v, _d, _o in b.values() if v == 1)
        vy = sum(1 for v, _d, _o in b.values() if v == -1)
        # Codex's decisive test: among reads whose OWN chain carried the count
        # LOSER, how many does the signal move to the winner? A genuine
        # alternative keeps a coherent subset loyal; a wobble's reads defect.
        lo_key = lo
        loyal = sum(1 for v, _d, o in b.values()
                    if o == lo_key and v != 0 and (key[0] if v == 1 else key[1]) == lo_key)
        defect = sum(1 for v, _d, o in b.values()
                     if o == lo_key and v != 0 and (key[0] if v == 1 else key[1]) == hi)
        n = vx + vy
        sig_pick = None
        if n >= a.min_electorate:
            if vx / n >= a.supermajority: sig_pick = key[0]
            elif vy / n >= a.supermajority: sig_pick = key[1]
        sig_right = (sig_pick == anno_pick) if (anno_pick and sig_pick) else None
        rows.append(dict(x=x, y=y, rx=rx, ry=ry, status=status, decision=decision,
                         count_right=count_right, n_ballots=n, vx=vx, vy=vy,
                         sig_pick=sig_pick, sig_right=sig_right,
                         loyal=loyal, defect=defect))

    print(f"registry junctions {len(reg)}   pairs within {a.tol} bp: {len(pairs)}")
    print(f"per-read signal records {n_rec}, usable {n_valid}, pairs with ballots {len(ballots)}\n")
    print(f"{'':10s} {'both anno':>10s} {'one anno':>9s} {'neither':>8s}")
    for d in ("moved", "abstain"):
        print(f"{d:10s} {C[(d,'both')]:10d} {C[(d,'one')]:9d} {C[(d,'neither')]:8d}")

    def acc(rs, k):
        v = [r[k] for r in rs if r[k] is not None]
        return (100 * sum(v) / len(v), len(v)) if v else (float("nan"), 0)

    print("\n--- (b) moved-but-shouldn't: pairs snap MOVED where BOTH coordinates are annotated")
    both_moved = [r for r in rows if r["decision"] == "moved" and r["status"] == "both"]
    print(f"    {len(both_moved)}  (signal cannot help these; each is a real site collapsed)")

    print("\n--- (a) the population signal could address: ONE-annotated ABSTENTIONS")
    one_abs = [r for r in rows if r["decision"] == "abstain" and r["status"] == "one"]
    print(f"    n = {len(one_abs)}")
    pc, nc = acc(one_abs, "count_right")
    print(f"    count ratio picks the annotated one: {pc:.1f}% of {nc} (ratio < {a.min_ratio}, "
          f"so it did NOT act on these)")
    called = [r for r in one_abs if r["sig_pick"] is not None]
    ps, ns = acc(one_abs, "sig_right")
    print(f"    signal calls {len(called)} of {len(one_abs)} "
          f"(electorate >= {a.min_electorate}, supermajority {a.supermajority}); "
          f"agrees with annotation on {ps:.1f}% of {ns}")

    print("\n--- calibration check: ONE-annotated pairs snap MOVED (count already decided)")
    one_moved = [r for r in rows if r["decision"] == "moved" and r["status"] == "one"]
    pc, nc = acc(one_moved, "count_right"); ps, ns = acc(one_moved, "sig_right")
    print(f"    n = {len(one_moved)}   count right {pc:.1f}% of {nc}   "
          f"signal right {ps:.1f}% of {ns} (called {sum(1 for r in one_moved if r['sig_pick'])})")

    print("\n--- loser-origin reads (Codex test): do the count LOSER's own reads defect to the winner?")
    for lab, sel in (("one-anno abstentions", one_abs), ("both-anno abstentions",
                     [r for r in rows if r["decision"] == "abstain" and r["status"] == "both"])):
        withb = [r for r in sel if r["loyal"] + r["defect"] >= 3]
        if not withb:
            print(f"    {lab}: no pair with >= 3 loser-origin ballots"); continue
        fr = [r["defect"] / (r["loyal"] + r["defect"]) for r in withb]
        print(f"    {lab}: {len(withb)} pairs with >= 3 loser-origin ballots; "
              f"median defect fraction {sorted(fr)[len(fr)//2]:.2f}; "
              f"pairs with >= 3 LOYAL loser reads (a real minority): "
              f"{sum(1 for r in withb if r['loyal'] >= 3)}")

    print("\n--- the flip: could loser-loyalty VETO the count moves that collapse real tandem sites?")
    for lab, sel in (("both-anno MOVED (real site collapsed)",
                      [r for r in rows if r["decision"] == "moved" and r["status"] == "both"]),
                     ("one-anno MOVED (count was right 99.9%)",
                      [r for r in rows if r["decision"] == "moved" and r["status"] == "one"])):
        withb = [r for r in sel if r["loyal"] + r["defect"] >= 3]
        veto = sum(1 for r in withb if r["loyal"] >= 3)
        print(f"    {lab}: {len(sel)} pairs, {len(withb)} with >= 3 loser-origin ballots, "
              f"veto (>= 3 loyal) would fire on {veto}"
              + (f"  <- {100*veto/max(1,len(withb)):.0f}% of testable real sites protected"
                 if "both" in lab else
                 f"  <- {100*veto/max(1,len(withb)):.0f}% of testable CORRECT moves blocked"))

    print("\n--- the danger: BOTH-annotated ABSTENTIONS where the signal would force a choice")
    both_abs = [r for r in rows if r["decision"] == "abstain" and r["status"] == "both"]
    forced = sum(1 for r in both_abs if r["sig_pick"] is not None)
    print(f"    n = {len(both_abs)}; signal reaches a supermajority on {forced} "
          f"({100*forced/max(1,len(both_abs)):.1f}%) -- each would collapse a real tandem site")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
