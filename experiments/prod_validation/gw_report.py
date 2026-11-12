#!/usr/bin/env python3
"""Detailed head-to-head report over every scored genome-wide cell.

Reads h2h/scores.tsv (expression truth contract) and h2h/passthrough.tsv (B1)
only; computes nothing that the scorer did not already measure. Sections:

  1. coverage          cells scored per tool
  2. overall           each tool against pyfin on the SAME cells (paired)
  3. vs stringtie3     per condition: pyfin and stringtie3 Sn/Pr/F1, paired
                       deltas, win counts, and a library-block bootstrap CI
                       on the mean F1 delta (libraries resampled whole)
  4. per library       pyfin vs stringtie3 aggregated over each library's cells
  5. condition groups  intact / reduced / de novo / corrupted
  6. library size      Spearman of library size against the Pr gap, per condition
  7. baselines         bambu vs B1 passthrough

Sn and Pr are reported everywhere F1 is. F1 means are means of per-cell F1,
never F1 of mean Sn and mean Pr.

Two means are reported, because they answer different questions and differ
whenever libraries contribute unequal numbers of cells (a competitor that
finished only part of the grid, or a campaign still running):

  每格平均   every cell counts once. The headline; its interval resamples whole
             libraries, so a library's cells stay together.
  文库等权   each library's cells are averaged first, then the libraries. A
             library with many finished cells cannot outweigh one with few.

On a complete balanced grid the two coincide.
"""
import collections
import csv
import math
import os
import random
import sys
from pathlib import Path

CAMP = Path("/SSD/logan/dev/pyfin/experiments/prod_validation/_gwcampaign")
TOOLS = ("pyfin", "stringtie3", "bambu", "isotools", "isoquant", "espresso", "lafite", "talon", "flair")
ORDER = ["full", "p99", "p90", "p80", "p50", "p20", "p10", "p05", "p01", "p00"]
GROUPS = {
    "注释完整 (full, p99, p90)": {"full", "p99", "p90"},
    "注释缺失 (p80 .. p01)": {"p80", "p50", "p20", "p10", "p05", "p01"},
    "从头组装 (p00)": {"p00"},
}


def lab(s):
    return (s.replace("SGNex_", "").replace("_directRNA_", " ")
             .replace("replicate", "r").replace("_run", "/r"))


def f(r, k):
    return float(r[k])


def load():
    cells = collections.defaultdict(dict)
    for r in csv.DictReader(open(CAMP / "h2h" / "scores.tsv"), delimiter="\t"):
        if r["status"].startswith("ok"):
            cells[(r["sample"], r["cond"])][r["tool"]] = r
    for r in csv.DictReader(open(CAMP / "h2h" / "passthrough.tsv"), delimiter="\t"):
        if r["status"].startswith("ok"):
            cells[(r["sample"], r["cond"])]["B1"] = r
    return cells


def lib_bootstrap(pairs, n=2000, seed=7):
    """pairs: list of (library, delta). Resample LIBRARIES with replacement."""
    by = collections.defaultdict(list)
    for lib_, d in pairs:
        by[lib_].append(d)
    libs = sorted(by)
    rng = random.Random(seed)
    means = []
    for _ in range(n):
        pick = [rng.choice(libs) for _ in libs]
        vals = [d for l in pick for d in by[l]]
        means.append(sum(vals) / len(vals))
    means.sort()
    return means[int(0.025 * n)], means[int(0.975 * n) - 1]


def library_sizes(samples):
    """BAM size in MB per library, nan where it cannot be read.

    The size is a covariate of sections 5 and 6 only, so an input that has
    moved costs those two columns and never the rest of the report.
    """
    size = dict.fromkeys(samples, float("nan"))
    try:
        sys.path.insert(0, str(Path(__file__).resolve().parent))
        os.environ["PYFIN_GRID_DATASET"] = "gencode"
        import sirv4_grid_runner as G
    except Exception:  # noqa: BLE001 - the runner resolves cluster paths at import
        return size
    for s in samples:
        try:
            size[s] = os.path.getsize(G.sample_inputs(s)[0]) / 1e6
        except Exception:  # noqa: BLE001
            pass
    return size


def lib_mean(pairs):
    """Mean over libraries of each library's mean. pairs: (library, value)."""
    by = collections.defaultdict(list)
    for lib_, d in pairs:
        by[lib_].append(d)
    return sum(sum(v) / len(v) for v in by.values()) / len(by) if by else float("nan")


def ranks(v):
    """Ranks with ties sharing their average rank."""
    order = sorted(range(len(v)), key=lambda i: v[i])
    r = [0.0] * len(v)
    k = 0
    while k < len(order):
        j = k
        while j + 1 < len(order) and v[order[j + 1]] == v[order[k]]:
            j += 1
        for i in order[k:j + 1]:
            r[i] = (k + j) / 2.0
        k = j + 1
    return r


def spearman(x, y):
    """Spearman rho; pairs with a non-finite value are left out, and the
    result is nan when fewer than three remain or either side is constant."""
    keep = [(a, b) for a, b in zip(x, y) if math.isfinite(a) and math.isfinite(b)]
    if len(keep) < 3:
        return float("nan")
    rx, ry = ranks([a for a, _ in keep]), ranks([b for _, b in keep])
    n = len(keep)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    den = (sum((a - mx) ** 2 for a in rx) * sum((b - my) ** 2 for b in ry)) ** 0.5
    return num / den if den else float("nan")


def main():
    cells = load()
    py = {k: t for k, t in cells.items() if "pyfin" in t}
    conds = [c for c in ORDER if any(c == k[1] for k in py)]
    conds += sorted({k[1] for k in py} - set(ORDER))

    print(f"## 1. 覆盖\n\npyfin 已评分 {len(py)} 格(样本 {len({k[0] for k in py})} 个,条件 {len(conds)} 个)\n")
    print("| 工具 | 与 pyfin 同格 |\n|---|---:|")
    for t in TOOLS[1:]:
        print(f"| {t} | {sum(1 for x in py.values() if t in x)} |")

    print("\n## 2. 总体(同格配对,每格平均)\n")
    print("最后一列是文库等权的 ΔF1:先在每个文库内平均,再对文库平均。\n")
    print("| 对手 | 格数 | 对手 Sn | 对手 Pr | 对手 F1 | pyfin Sn | pyfin Pr | pyfin F1 | ΔF1 | pyfin 领先 | ΔF1 文库等权 |")
    print("|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
    for t in TOOLS[1:]:
        pr = [(k[0], x[t], x["pyfin"]) for k, x in py.items() if t in x]
        n = len(pr)
        if not n:   # the tool shares no scored cell with pyfin: say so, compute nothing
            print(f"| {t} | 0 | — | — | — | — | — | — | — | — | — |")
            continue
        m = lambda i, k: sum(f(p[i], k) for p in pr) / n
        w = sum(1 for _, a, b in pr if f(b, "tx_f1") > f(a, "tx_f1"))
        eq = lib_mean([(s, f(b, "tx_f1") - f(a, "tx_f1")) for s, a, b in pr])
        print(f"| {t} | {n} | {m(1,'tx_sn'):.1f} | {m(1,'tx_pr'):.1f} | {m(1,'tx_f1'):.1f} | "
              f"{m(2,'tx_sn'):.1f} | {m(2,'tx_pr'):.1f} | {m(2,'tx_f1'):.1f} | "
              f"{m(2,'tx_f1')-m(1,'tx_f1'):+.2f} | {w}/{n} | {eq:+.2f} |")

    print("\n## 3. 按条件:pyfin vs stringtie3\n")
    print("ΔF1 的 95% 置信区间用文库分块 bootstrap(整库重抽样);格数 < 3 不给区间。\n")
    print("| 条件 | 格 | pyfin Sn/Pr/F1 | stringtie3 Sn/Pr/F1 | ΔSn | ΔPr | ΔF1 | 95% CI | F1 领先 | Sn 领先 | Pr 领先 |")
    print("|---|---:|---|---|---:|---:|---:|---|---:|---:|---:|")
    allpairs = []
    for c in conds:
        rows = [(k[0], t) for k, t in py.items() if k[1] == c and "stringtie3" in t]
        if not rows:
            continue
        n = len(rows)
        m = lambda tool, k: sum(f(t[tool], k) for _, t in rows) / n
        d = [(s, f(t["pyfin"], "tx_f1") - f(t["stringtie3"], "tx_f1")) for s, t in rows]
        allpairs += d
        ci = lib_bootstrap(d) if n >= 3 else None
        ws = lambda k: sum(1 for _, t in rows if f(t["pyfin"], k) > f(t["stringtie3"], k))
        cis = f"[{ci[0]:+.2f}, {ci[1]:+.2f}]" if ci else "—"
        print(f"| {c} | {n} | {m('pyfin','tx_sn'):.1f}/{m('pyfin','tx_pr'):.1f}/{m('pyfin','tx_f1'):.1f} | "
              f"{m('stringtie3','tx_sn'):.1f}/{m('stringtie3','tx_pr'):.1f}/{m('stringtie3','tx_f1'):.1f} | "
              f"{m('pyfin','tx_sn')-m('stringtie3','tx_sn'):+.1f} | {m('pyfin','tx_pr')-m('stringtie3','tx_pr'):+.1f} | "
              f"{m('pyfin','tx_f1')-m('stringtie3','tx_f1'):+.2f} | {cis} | {ws('tx_f1')}/{n} | {ws('tx_sn')}/{n} | {ws('tx_pr')}/{n} |")
    if allpairs:
        lo, hi = lib_bootstrap(allpairs)
        w = sum(1 for _, d in allpairs if d > 0)
        print(f"\n**合计 {len(allpairs)} 格:pyfin 领先 {w}/{len(allpairs)},平均 ΔF1 "
              f"{sum(d for _, d in allpairs)/len(allpairs):+.2f},文库分块 95% CI [{lo:+.2f}, {hi:+.2f}];"
              f"文库等权 ΔF1 {lib_mean(allpairs):+.2f}**")
    else:
        print("\n**没有 pyfin 与 stringtie3 同时评分的格。**")

    print("\n## 4. 按条件组\n")
    print("| 组 | 格 | pyfin Sn/Pr/F1 | stringtie3 Sn/Pr/F1 | ΔF1 | 95% CI | 领先 |")
    print("|---|---:|---|---|---:|---|---:|")
    groups = dict(GROUPS)
    groups["注释腐蚀 (c_*)"] = {c for c in conds if c.startswith("c_")}
    for g, cs in groups.items():
        rows = [(k[0], t) for k, t in py.items() if k[1] in cs and "stringtie3" in t]
        if not rows:
            continue
        n = len(rows)
        m = lambda tool, k: sum(f(t[tool], k) for _, t in rows) / n
        d = [(s, f(t["pyfin"], "tx_f1") - f(t["stringtie3"], "tx_f1")) for s, t in rows]
        lo, hi = lib_bootstrap(d)
        print(f"| {g} | {n} | {m('pyfin','tx_sn'):.1f}/{m('pyfin','tx_pr'):.1f}/{m('pyfin','tx_f1'):.1f} | "
              f"{m('stringtie3','tx_sn'):.1f}/{m('stringtie3','tx_pr'):.1f}/{m('stringtie3','tx_f1'):.1f} | "
              f"{m('pyfin','tx_f1')-m('stringtie3','tx_f1'):+.2f} | [{lo:+.2f}, {hi:+.2f}] | "
              f"{sum(1 for _, x in d if x > 0)}/{n} |")

    print("\n## 5. 按文库:pyfin vs stringtie3\n")
    size = library_sizes({k[0] for k in py})
    unknown = sorted(s for s, v in size.items() if not math.isfinite(v))
    if unknown:
        print(f"{len(unknown)} 个文库的 BAM 大小读不到(记为 nan,不进入第 6 节):"
              f"{', '.join(lab(s) for s in unknown)}\n")
    print("| 文库 | BAM MB | 格 | pyfin Sn/Pr/F1 | stringtie3 Sn/Pr/F1 | ΔSn | ΔPr | ΔF1 | 领先 |")
    print("|---|---:|---:|---|---|---:|---:|---:|---:|")
    for s in sorted(size, key=lambda s: (0, size[s], s) if math.isfinite(size[s]) else (1, 0.0, s)):
        rows = [t for k, t in py.items() if k[0] == s and "stringtie3" in t]
        n = len(rows)
        if not n:
            print(f"| {lab(s)} | {size[s]:.0f} | 0 | — | — | — | — | — | — |")
            continue
        m = lambda tool, k: sum(f(t[tool], k) for t in rows) / n
        w = sum(1 for t in rows if f(t["pyfin"], "tx_f1") > f(t["stringtie3"], "tx_f1"))
        print(f"| {lab(s)} | {size[s]:.0f} | {n} | {m('pyfin','tx_sn'):.1f}/{m('pyfin','tx_pr'):.1f}/{m('pyfin','tx_f1'):.1f} | "
              f"{m('stringtie3','tx_sn'):.1f}/{m('stringtie3','tx_pr'):.1f}/{m('stringtie3','tx_f1'):.1f} | "
              f"{m('pyfin','tx_sn')-m('stringtie3','tx_sn'):+.1f} | {m('pyfin','tx_pr')-m('stringtie3','tx_pr'):+.1f} | "
              f"{m('pyfin','tx_f1')-m('stringtie3','tx_f1'):+.2f} | {w}/{n} |")

    print("\n## 6. 文库大小与精确率差(Spearman,每条件 ≥ 8 格才算)\n")
    print("并列的值取平均秩;大小未知的文库不参与。\n")
    print("| 条件 | 格 | 库大小 vs ΔPr | vs pyfin Pr | vs stringtie3 Pr | 库大小 vs ΔF1 |")
    print("|---|---:|---:|---:|---:|---:|")
    for c in conds:
        rows = [(size[k[0]], t) for k, t in py.items()
                if k[1] == c and "stringtie3" in t and math.isfinite(size[k[0]])]
        if len(rows) < 8:
            continue
        x = [r[0] for r in rows]
        dp = [f(t["pyfin"], "tx_pr") - f(t["stringtie3"], "tx_pr") for _, t in rows]
        pp = [f(t["pyfin"], "tx_pr") for _, t in rows]
        sp = [f(t["stringtie3"], "tx_pr") for _, t in rows]
        df = [f(t["pyfin"], "tx_f1") - f(t["stringtie3"], "tx_f1") for _, t in rows]
        print(f"| {c} | {len(rows)} | {spearman(x, dp):+.2f} | {spearman(x, pp):+.2f} | "
              f"{spearman(x, sp):+.2f} | {spearman(x, df):+.2f} |")

    print("\n## 7. 基线:bambu vs B1 透传(把输入 GTF 原样输出)\n")
    both = {k: t for k, t in cells.items() if "bambu" in t and "B1" in t}
    print(f"B1 只在 part A 的格上跑过:这里是 {len(both)} 格,不是整个 campaign。"
          "Δ 是 bambu 减 B1,带符号;F1 接近不代表 Sn 和 Pr 都接近,所以三者都列。\n")
    print("| 条件组 | 格 | bambu Sn/Pr/F1 | B1 Sn/Pr/F1 | ΔSn | ΔPr | ΔF1 | 平均 \\|ΔF1\\| | \\|ΔF1\\| ≤ 0.5 |")
    print("|---|---:|---|---|---:|---:|---:|---:|---:|")
    bg = collections.defaultdict(list)
    for k, t in both.items():
        c = k[1]
        g = ("p00–p05" if c in ("p00", "p01", "p05") else
             "p10–p50" if c in ("p10", "p20", "p50") else "完整 / 近完整 / c_*")
        bg[g].append(t)
    for g in ("完整 / 近完整 / c_*", "p10–p50", "p00–p05"):
        v = bg[g]
        if not v:
            continue
        m = lambda tool, k: sum(f(t[tool], k) for t in v) / len(v)
        ad = [abs(f(t["bambu"], "tx_f1") - f(t["B1"], "tx_f1")) for t in v]
        print(f"| {g} | {len(v)} | {m('bambu','tx_sn'):.1f}/{m('bambu','tx_pr'):.1f}/{m('bambu','tx_f1'):.1f} | "
              f"{m('B1','tx_sn'):.1f}/{m('B1','tx_pr'):.1f}/{m('B1','tx_f1'):.1f} | "
              f"{m('bambu','tx_sn')-m('B1','tx_sn'):+.1f} | {m('bambu','tx_pr')-m('B1','tx_pr'):+.1f} | "
              f"{m('bambu','tx_f1')-m('B1','tx_f1'):+.2f} | {sum(ad)/len(ad):.2f} | "
              f"{sum(1 for x in ad if x <= 0.5)}/{len(ad)} |")


if __name__ == "__main__":
    main()
