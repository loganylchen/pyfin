#!/usr/bin/env python3
"""Merge a run's selection-audit shards and attribute a target list to its fate.

Reads <out>/selection_audit/*.tsv (interval stages) and
<out>/selection_audit_global.tsv, writes one merged table, and for each target
candidate id prints the LAST decisive record: the stage at which it was
dropped or folded, or "emitted" if it is in the final set.
"""
import argparse, collections, csv
from pathlib import Path

COLS = ("interval", "stage", "candidate_id", "source", "n_reads", "action", "reason", "fold_into")
ORDER = ["discovered", "canonical_gate", "junction_dominance_gate", "entered_em",
         "interval_cascade", "interval_survivor", "extent_unification",
         "global_cascade", "junction_snap", "post_refit", "post_snap", "final"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True, help="pyfin output dir (work_dir)")
    ap.add_argument("--targets", help="file of candidate ids to attribute")
    ap.add_argument("--merged", help="write the merged table here")
    ap.add_argument("--fates", help="write per-target fate table here")
    a = ap.parse_args()
    out = Path(a.out_dir)
    rows = []
    for shard in sorted((out / "selection_audit").glob("*.tsv")):
        with open(shard) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) == len(COLS):
                    rows.append(dict(zip(COLS, f)))
    g = out / "selection_audit_global.tsv"
    if g.exists():
        with open(g) as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    if a.merged:
        with open(a.merged, "w") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t")
            w.writeheader(); w.writerows(rows)
    by_cid = collections.defaultdict(list)
    for r in rows:
        by_cid[r["candidate_id"]].append(r)
    rank = {s: i for i, s in enumerate(ORDER)}

    final_ids = {r["candidate_id"] for r in rows if r["stage"] == "final"}

    def fate(cid, depth=0):
        recs = sorted(by_cid.get(cid, []), key=lambda r: rank.get(r["stage"], 99))
        if not recs:
            return ("never_discovered", "", "", "")
        if cid in final_ids:
            return ("emitted", "final", "", "")
        decisive = [r for r in recs if r["action"] in ("drop", "fold")]
        if decisive:
            d = decisive[-1]
            if d["action"] == "fold" and d["fold_into"] and depth < 8:
                # a fold is not a loss if the survivor it points to (transitively)
                # is emitted: the folded chain rides on that model
                f = fate(d["fold_into"], depth + 1)
                if f[0] in ("emitted", "emitted_via_fold"):
                    return ("emitted_via_fold", d["stage"], d["reason"], d["fold_into"])
            return (d["action"], d["stage"], d["reason"], d["fold_into"])
        last = recs[-1]
        return ("vanished_after", last["stage"], "", "")

    print(f"records={len(rows)} candidates={len(by_cid)}")
    stage_counts = collections.Counter((r["stage"], r["action"], r["reason"]) for r in rows)
    for (s, act, reason), n in sorted(stage_counts.items(), key=lambda x: (rank.get(x[0][0], 99), -x[1])):
        print(f"  {s:24s} {act:8s} {reason:32s} {n}")
    if a.targets:
        targets = [t for t in open(a.targets).read().split() if t]
        fates = [(t, *fate(t)) for t in targets]
        if a.fates:
            with open(a.fates, "w") as fh:
                fh.write("candidate_id\tfate\tstage\treason\tfold_into\n")
                for f in fates:
                    fh.write("\t".join(f) + "\n")
        print(f"\ntargets={len(targets)}")
        for k, n in collections.Counter((f[1], f[2], f[3]) for f in fates).most_common():
            print(f"  {k[0]:16s} {k[1]:20s} {k[2]:32s} {n}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
