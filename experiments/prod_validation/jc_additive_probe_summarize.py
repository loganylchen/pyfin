#!/usr/bin/env python3
"""Aggregate the additive junction-correction candidate probe grid.

Per config (margin M, runner-up R, min voters K): sums over cells of new
chains (TP by catalog / expressed truth, FP), lost chains, and the read-level
proposal moves. Hard constraints: lost_chains == 0 and lost_TP == 0 (additive
semantics). Metric: new_TP / (new_TP + new_FP).
"""
import collections
import json
import sys
from pathlib import Path

root = Path(sys.argv[1] if len(sys.argv) > 1 else
            "experiments/prod_validation/_sirv4grid")
rows = []
for d in sorted(root.glob("jc_additive_m*_r*_k*")):
    if not d.is_dir():
        continue
    for f in sorted(d.glob("*.summary.json")):
        s = json.loads(f.read_text())
        s["config"] = d.name
        rows.append(s)
by = collections.defaultdict(list)
for s in rows:
    by[s["config"]].append(s)
keys = ("chains_before", "chains_after", "new_chains", "new_TP_cat", "new_TP_exp",
        "new_FP", "new_TP_reads", "new_FP_reads", "lost_chains", "lost_TP_cat",
        "reads_proposing", "truth_exp_covered_before", "truth_exp_covered_after")
print("config\tcells\t" + "\t".join(keys) + "\tprecision_cat\tprecision_exp")
for cfg, ss in sorted(by.items()):
    tot = {k: sum(s[k] for s in ss) for k in keys}
    p_cat = tot["new_TP_cat"] / max(1, tot["new_chains"])
    p_exp = tot["new_TP_exp"] / max(1, tot["new_chains"])
    print(f"{cfg}\t{len(ss)}\t" + "\t".join(str(tot[k]) for k in keys) +
          f"\t{p_cat:.3f}\t{p_exp:.3f}")
print("\nper-condition (new_TP_cat / new_FP / lost) by config:")
for cfg, ss in sorted(by.items()):
    bc = collections.defaultdict(lambda: [0, 0, 0, 0])
    for s in ss:
        b = bc[s["cond"]]
        b[0] += s["new_TP_cat"]; b[1] += s["new_FP"]; b[2] += s["lost_chains"]; b[3] += 1
    print(cfg, {c: f"TP{v[0]}/FP{v[1]}/lost{v[2]} n={v[3]}" for c, v in sorted(bc.items())})
moves = collections.defaultdict(collections.Counter)
for s in rows:
    moves[s["config"]].update(s.get("proposal_moves", {}))
print("\nproposal moves (read-level, own chain TP/FP -> proposed TP/FP):")
for cfg, m in sorted(moves.items()):
    print(cfg, dict(sorted(m.items())))
