# 人类 recall 差距的根因与"兼容性作证"豁免（2026-09-09/10）

状态：开发中，**默认关闭的实验性特性**；`real-drna` profile 不改。
当前实现 = **v6**（快照 `experiments/prod_validation/_snap/attest6`，加上提交前的
门控收集/几何校验两处修复 = `_snap/attest7`）。
§3 的 v1/v2/v3 和 §4 的 v4 面板都是**历史版本**的数字，保留是为了说明每一步为什么这样改，
不要当成当前实现的结果；当前实现的四个面板在 §4d。
前置记录：`HANDOFF.md` §4b5（人类 26 格头对头）、§4o（本轮）；`docs/LITERATURE_REVIEW_2026-09-09.md`。
每一步都经 Codex（gpt-6-astra）审查；两轮意见原文摘录在文末。

## 1. 问题：人类 GENCODE 上我们输给 stringtie3 的是什么

H9 rep2 run2 chr1/19/21，表达合同（真值 = GENCODE v44 中 NanoCount TPM>0 的转录本，gffcompare `=`）。
竞品网格已跑完（stringtie3/bambu/isotools 各 15 样本 × 26 条件），我们的生产默认（`GCSMC`）补齐了全部 26 个条件：

- 26 格均值 +1.41，**13/26 胜**，交叉点在 p80（+0.01）。
- 注释 ≥90% 完整、或只有 10% 轻度腐败：stringtie3 每格领先 0.8–2.0（full 44.68 vs 42.75）。
- 注释 ≤50% 或重度腐败：我们领先 2–10（c_jitter50 +10.4、c_flip50 +8.6）。
- **每一格我们 precision 都更高；差的全是 recall。**

### 1.1 `full` 上的差距解剖（`human_gap_dissect.py`）

| 集合 | n | 精确全链 read=0 | 只有 3' 后缀 read | NanoCount 中位数 |
|---|---:|---:|---:|---:|
| 双方都命中 | 1938 | 11% | 176 | 3.3 |
| stringtie3 独有 | 554 | 51% | 156 | 1.0（345 个 <1） |
| pyfin 独有 | 355 | 94% | 259 | 1.0 |

stringtie3 独有的 554 个里有 272 个在 BAM 里至少一条 read 精确携带全链（152 个 ≥2 条），
但我们在 471 个上一个模型都没有输出。我们独有的 355 个则几乎全靠 GTF 先验（5' 降解 read + 注释）。

### 1.2 它们死在哪一层（`--selection-audit`，新工具）

生产代码原本不记录 `SelectionOutcome`。加了 opt-in 的 `--selection-audit`：每个 interval 一个 shard
（stage：discovered / canonical_gate / junction_dominance_gate / entered_em / interval_cascade /
interval_survivor），主进程写 extent_unification / global_cascade / junction_snap / post_refit / post_snap / final。
GTF 输出逐字节不变（真实单元验证：模型多重集相同）。`selection_audit_collect.py` 合并 shard，
按目标 id 归因（fold 会沿 redirect 追到被输出的幸存者）。

554 个 stringtie3 独有转录本的归因：**475 个被 `m2_support_gate` 丢掉**，28 个 extent 折叠，
24 个 mono_resolve，17 个 cluster_recheck，8 个 min_abundance。
全部 7398 个表达真值：**4566 个（62%）死于 `m2_support_gate`**。该门在这个单元丢掉了 31,418 个进入 EM 的候选中的 20,294 个。
554 个在 discovery 层全部存活为 GTF 候选，没有一个被 fold。

### 1.3 被 gate 丢掉的注释候选，按 read 证据分层（`gate_drop_precision.py`）

| 证据层（后缀计数含精确） | n | 真值率 |
|---|---:|---:|
| 精确全链 read ≥3 | 122 | 75% |
| 精确 =2 | 71 | 65% |
| 精确 =1 | 347 | 45% |
| 无精确、3' 后缀 ≥3 | 1367 | 77.5% |
| 无精确、后缀 1–2 | 2160 | 59% |
| 没有任何链 read | 15,929 | 12% |

精确∪后缀 ≥3：1595 个候选，77.0% 真值——高于当时 67.9% 的总体精度。
Codex 的提醒：这个"真值率"是与 NanoCount 伪真值的一致率（NanoCount 本身按兼容性摊 read），不是独立的表达证据。

## 2. 豁免规则（`--m2-support-gate-attest-reads K`，默认 0 = 逐字节惰性）

注释的多外显子候选在 `m2_support_gate` 处豁免，当且仅当：

1. `compatible_read_depth ≥ K`：不同 query name 的 MAPQ≥1 primary read，其 intron 链**等于**候选链、或是候选链的
   链方向 3' **后缀**（'+' 取链尾，'-' 取链头，intron 按基因组升序）；read 跨度不越过后缀上游那条模型内含子（容差 10 bp，
   否则是保留内含子而非降解片段）；read 3' 端不超出注释 3' 端 50 bp（只限制超出，不限制提前终止）。
   在 `discovery.chain_compatible_read_depths` 里总是计算；`.ext` 兄弟继承其父注释的深度（链级支持）。
2. 否决：候选的任一**未观测** junction 若在 **两端都** ±W（默认 10 bp）内有观测 junction，视为该 junction 的 jitter，
   整个候选不豁免（`m2_junction_nll.unjittered_gtf_columns`）。观测集合用 `CandidateSet.junction_support`
   （MAPQ≥1，与深度同一 read 策略）。集合为空则 fail-closed。
3. 不要求 EM 质量；其它门（abundance 地板、isoform-fraction、canonical、containment）照常。

不做的事（都测过或被 Codex 否决）：
- **不**要求 read 自身所有 junction indel-clean（v1）：在 dRNA 上这是 read 长度过滤器，把 1079 砍到 210 而精度只从 90% 到 95%；
  HANDOFF 里 complete-chain rescue 当年就否决过"全 read AND"。
- **不**用"共用任一剪接位点 ±6 bp 即矛盾"的 OR 形否决（v1）：在人类上否决的 129 个里 128 个是表达的——那是可变剪接。
- **不**用"精确匹配某注释链的 read 只为该链作证"（v3）：Codex：精确链匹配不等于完整分子（长异构体的截断 read 可以恰好等于短注释链）；
  在腐败阶梯上它从不触发（真链总被删掉），只在人类花 0.3。
- **不**加 TSS 尖峰阻断：离线看 M≥5 会丢 181 个真值；Codex：5' 端集中既不证明完整分子也不证明长异构体零贡献，不能"peak ⇒ complete"。

## 3. 哨兵结果（一库一格，配对 vs 生产默认；单库 CI 不是生物学重复）

| | v1（B\* 原样） | v2（去 indel 策略 + jitter 形否决） | v3（v2 + exclusive） |
|---|---|---|---|
| 人类 full | +0.30 | **Sn 31.2→34.8, Pr 67.9→69.3, F1 42.75→46.33 (+3.58)** | +3.29 |
| 人类 c_jitter10 / 50 | – | +0.83 / +0.06 | +0.72 / +0.06 |
| 人类 p00 | 0 | 0 | 0 |
| SIRV4 lib（26） | −0.05 | −0.12 [−0.28, 0.00]；ir50/merge50/skip50 各 −0.9/−0.9/−1.3 | 同 v2 |
| HEYA8 lib（26） | +0.01 | −0.03；skip50 −0.68 | 同 v2 |

v2 在人类 `full` 上超过 stringtie3（44.68），两轴同升。K=0 对照与生产模型多重集完全相同（v1、v2 均验证）。

**spike-in 的 FP 是什么**：c_merge50 的 SIRV409 是"SIRV409 5' 部分 + SIRV410 链"的合并注释，而 SIRV410 自己的条目被删掉；
SIRV410 的 35 条 read 等于任何注释链都不等于 → 按规则是"片段" → 为嵌合体作证。c_skip50 的 SIRV403 同形；
SIRV303/306 被腐败成与 SIRV302 同链（重复链条目）；c_ir50 是保留内含子条目，缺失的 junction 是沉默。
**结论：腐败 junction 落在 read 覆盖之外时（merge/skip/ir），链级兼容性无法区分它和人类的降解片段。**

**Codex 的判定**：预先冻结的停止规则（腐败家族均值 ≥ −0.5）在 SIRV4 哨兵的 c_skip（−0.67）上被触发；
后续 v2 全面板是**表征**，不是未中断的晋级测试。修复后的 arm 需要新冻结的协议：人类均值 ΔF1 ≥ +1 且配对下界 > 0；
其它每个数据集下界 > −0.25；每个预定义腐败家族均值 ≥ −0.5；缺格阻断决策；损失不能跨数据集抵消；bootstrap 按文库分块。

## 4. 表征面板（快照 `_snap/attest4` = 修复后的 v2，jitter 否决）

| 数据集 | 格数/库数 | ΔSn | ΔPr | ΔF1 [95% CI] | 升/降/同 | 越界家族 |
|---|---|---:|---:|---|---|---|
| 人类 | 26 / 1 | **+2.21** | **+0.64** | **+2.19 [+1.70, +2.67]** | **25 / 0 / 1** | 无（全家族为正） |
| SIRV4 | 156 / 6 | +0.55 | −0.49 | −0.00 [−0.19, +0.27] | 51 / 37 / 68 | **c_skip −0.62** |
| Sequin | 104 / 4 | +0.00 | −0.19 | −0.08 [−0.11, −0.06] | 0 / 7 / 97 | **c_merge −0.77** |
| HEYA8 | 130 / 5 | +0.20 | −0.20 | +0.02 [−0.08, +0.18] | 17 / 13 / 100 | 无（c_merge −0.31） |

K=0 同快照对照在 `full` 上与生产逐格相同。

**人类 26 格对 stringtie3**：生产 13/26 胜 → v4 **25/26 胜**（均值 +3.60）。唯一输的是 p00（−0.30），
那是 de novo 单元，没有注释，本杠杆按构造不生效。逐格差（v4 − stringtie3）：full +1.63、p99 +1.88、p90 +2.44、
p80 +3.00、p50 +4.25、p20 +5.59、c_merge10 +1.54、c_trunc3p50 +1.17、c_jitter50 +10.46、c_flip50 +10.52。

## 4b. v5：加"未注释 junction"否决（快照 `_snap/attest5`）

发现：OR 形否决在人类扔掉的 652 个里 598 个（92%）按合同是真值——比放进来的还高；它是可变剪接探测器。
但腐败场景有人类没有的特征：**与候选未观测 junction 共用位点的那个观测 junction，在腐败后的注释里不属于任何注释链**
（merge/skip 把真链替换掉了；IR 把剪接掉的 junction 变成了外显子内部）。人类的可变受体/供体属于另一条注释异构体。

v5 否决 = jitter 规则 ∪ 未注释规则：
(a) 一个**不属于本 interval 任何注释链**的观测 junction 与候选未观测 junction 共用一个剪接位点（±10 bp）；
(b) 这样的 junction 整个落在候选的一个外显子内（长度 > 10 bp）。
注释集合取自 interval 内**所有** GTF 候选（含已被前面的门 mask 掉的），观测集合 = `junction_support`（MAPQ≥1）。
离线（人类 full，1075 个深度≥3 候选）：(a)∪(b) 否决 356（333 真值），保留 719/638；OR 否决保留 426/373。

v5 哨兵（配对 vs 生产）：

| | v4（jitter 否决） | **v5（jitter ∪ 未注释）** |
|---|---|---|
| 人类 full | +3.58 | **+2.76**（Sn 31.2→34.0, Pr 67.9→68.8；stringtie3 44.68） |
| 人类 c_merge50 | – | **+1.54**（37.58→39.13；stringtie3 37.68） |
| 人类 p50 | – | **+1.46**（Sn +1.3, Pr +0.8） |
| SIRV4 lib（26） | −0.12，3 格降 | **26/26 相同** |
| Sequin lib（26） | −0.12，3 格降 | **26/26 相同** |

### 4b.0 v5 面板（快照 `_snap/attest5`）与 v6 的 P1 修复

v5 spike-in 面板：SIRV4 156/6 **+0.09 [−0.11, +0.34]**（ΔSn +0.36 ΔPr −0.15，37 升/16 降/103 同，全家族 ≥ +0.07）；
Sequin 104/4 **−0.01 [−0.03, +0.00]**（0 升/1 降/103 同，仅 c_ir −0.15）；HEYA8 130/5 **+0.05 [+0.01, +0.12]**
（12 升/2 降/116 同，全家族 ≥ 0）。人类 13/26 格已出，`full` +2.76、c_merge50 +1.54、p50 +1.46、p00 0.00。

**Codex 第三轮（提交前审查）发现一个 P1 集成 bug**：rule 2 的"已注释 junction 全集"是从 `live`
（被前面门 mask 掉的列置为 `None`）构建的，于是**一个已被丢弃的注释，其 junction 不再算作注释**——
"这条 read junction 是不是新的"变成了依赖选择结果的问题。修复：`unjittered_gtf_columns` 增加显式参数
`annotated_junctions`，由调用方从**完整的 `cand_list`**（含 masked）构建；helper 的默认值仅供直接调用/测试。
同时按 Codex 要求：畸形几何（span 不包住自己的链、intron 乱序/重叠）拒绝豁免；审计目录清理失败改为 fatal；
修正 docstring 里写反的 polyA 侧。因此 v4/v5 的面板数字**属于修复前的实现**，v6（`_snap/attest6`）重跑全部四个面板。

### 4b.1 v5/v6 冻结协议（跑面板前写下）

- 快照 `_snap/attest5`（= 工作树 fin/，`git diff` 于 HEAD 30f14b3 之上）；arm `--m2-support-gate-attest-reads 3`（veto 默认 `novel`，W=10）；
  基线 = 生产默认（GCSMC / SMC 各网格已有）。
- 面板：人类 26（`gcatt5`）、SIRV4 156（`_sirv4grid/att5`）、Sequin 104（`_sequingrid/att5`）、HEYA8 125（`_heya8grid/att5`）。
  哨兵格重跑（同快照、同参数，按 runner 的幂等规则会覆盖），不单独抵扣。
- 晋级条件（Codex 第二轮原文，全部满足才晋级）：人类均值 ΔF1 ≥ +1 且配对下界 > 0（单库：格级 bootstrap，
  明示为"库内基准结论"）；SIRV4 / Sequin / HEYA8 每个数据集按库分块 bootstrap 下界 > −0.25；每个预定义腐败家族
  （c_jitter/c_merge/c_skip/c_flip/c_ir/c_spurious/c_trunc3p/c_trunc5p）在每个数据集上的均值 ≥ −0.5；缺格阻断；
  损失不跨数据集抵消；Sn/Pr 分开报。
- 敏感性（晋级后补做，不影响判定）：NanoCount est_count ≥1 / ≥3 的真值重评分；冗余模型计数。
- 承认：v5 的否决形状是在看过 v2 面板和 SIRV 失败案例之后选定的（post-hoc），不是事前注册；这在文档里明说。
  Codex 的定性：v5 面板是**选择后的开发表征**，通过阈值不能恢复确证地位，bootstrap 区间没有计入模型选择；
  要变成晋级，需要在**未被用于设计规则的数据**上事前冻结确认（两个有基因组比对但未 stage 的人类文库是可行选项；
  已进入面板汇总的库不能事后当作 held-out）。
- 本次提交的目标：**默认关闭（K=0，逐字节惰性）的实验性特性 + 审计工具**；`real-drna` 默认不改。
  Codex：默认开启的晋级需要修复后实现的回归结果、在更多人类文库上事前冻结的确认、注释不完整压力测试，
  以及更严真值/独立支持检查。`sirv` profile 不动。

## 4c. p00（de novo）为什么不在这条杠杆的射程内——以及它的天花板在哪

人类 26 格里唯一还输给 stringtie3 的是 p00（24.67 vs 24.97），那是没有注释的 de novo 单元，
作证豁免按构造不生效（它只豁免 ANNOTATED 候选）。对 p00 做同样的审计（`gcaud00`，
`denovo_truth_fate.py` 按内含子链把真值映到随机 uuid 候选再查审计）：

| 7398 条表达真值 | n |
|---|---:|
| **从未生成**（没有任何 read 的 CIGAR 链等于它） | **5239（71%）** |
| 生成为候选 | 1647 |
| 生成为折叠 shadow | 351 |
| 单外显子（度量口径外，见下） | 157 |

生成出来的 1998 条里 **1380 条被输出（69%）**；丢弃分布：`novel_junction_support` 533、
canonical_gate 24、cluster_recheck 24、min_abundance 23、其它 < 15。
门的规模也和 `full` 相反：`m2_support_gate` 在 p00 只丢 162（没有注释可门），
最大的是 `novel_junction_support` 3450。

**`novel_junction_support` 这一层没有杠杆**：3450 个丢弃里只有 533（15.4%）携带表达真值，
且 87%（2997）是单 read 候选；F1 的边际平衡点约 22%，我们当前 p00 精度 56.2%——
放进来两轴都掉。

**所以 p00 的 −0.30 不是选择层问题**：71% 的表达真值在这个文库里没有任何一条 read 携带完整链
（dRNA 5' 降解），那是数据上限；生成出来的已经输出了 69%；唯一的大门在这个数据集上是负价值的。

（口径说明：`read_chains` 把空链的单外显子 read 只按染色体/链归组，会把互不相关的无剪接 read 池在一起，
所以"157 条单外显子真值有 ≥2 条精确 read"是这个池化造成的假象，不是生成 bug。上表把它们单列。）

## 4d. 当前实现（v6）的四个面板

配对 vs 生产默认；>1 库时按库分块 bootstrap，单库时是格级 bootstrap 并明示"库内结论"。

| 数据集 | 格数/库数 | ΔSn | ΔPr | ΔF1 [95% CI] | 升/降/同 | 最差腐败家族 |
|---|---|---:|---:|---|---|---|
| 人类 | 26/1 | +1.65 | +0.42 | **+1.63 [+1.25, +2.02]** | **24 / 0 / 2** | c_jitter +0.34 |
| SIRV4 | 156/6 | +0.36 | −0.15 | +0.09 [−0.11, +0.34] | 37 / 16 / 103 | c_trunc3p +0.07 |
| Sequin | 104/4 | +0.00 | −0.03 | −0.01 [−0.03, +0.00] | 0 / 1 / 103 | **c_ir −0.15** |
| HEYA8 | 130/5 | +0.10 | −0.02 | +0.05 [+0.01, +0.12] | 12 / 2 / 116 | c_jitter +0.05 |

人类 26 格对 stringtie3（端点盲合同）：生产 13/26 胜（均值 +1.41）→ **v6 25/26 胜（均值 +3.04）**，
唯一还输的是 p00（−0.30），见 §4c。端点感知合同下的同一比较见
`docs/PRODUCTION_QUANT_ACCURACY_2026-09-10.md` §8——那里我们并不领先。

K=0 同快照对照在人类 `full` 上与生产 **逐字节相同**；把提交前的两处修复
（门控收集、几何校验前置）打进快照后，K=3 在 `full` 和 `c_merge50` 上与已跑面板的产物
**逐字节相同**（`_snap/attest8` vs `_snap/attest6`），所以面板数字可以归属到被提交的这份代码。

**提交前验证抓到的一个真 bug**（记下来，因为它说明为什么协议要求"验证"而不是"论证"）：
`collect_compatible_depths` 加进了外层和内层签名，但**转发那一行在一次脚本异常里丢了**，
于是门控版本的豁免恒不触发（输出退回生产的 3377 个模型）。
当时的单元测试只验了 runner 传参、没验端到端，所以没抓到；现在补了一个检查转发存在的测试。

**这些数字不是确证。** 否决规则的形状（rule 2）是在看过 v2 面板和 SIRV 失败案例之后选定的，
所以这是**选择后的开发表征**；bootstrap 区间没有计入模型选择；人类只有一个文库，
`+1.21` 是库内的格级下界，不是种群级的。三个数据集上 precision 都有小幅下降（虽在容差内）。
不要写"四个数据集验证通过"、"通过晋级合同"、"无回归"或"最优"。

## 5. Codex 对本杠杆的定性（原文）

> "This annotation-prior rescue exempts multi-exon models from the support gate when at least three distinct eligible primary
> reads have an exact or 3′-suffix intron chain and pass the implemented span checks, provided no unobserved model junction has
> an observed neighbor within 10 bp at both boundaries. On the inspected H9 full-annotation cell, v2 improved F1 by 3.58 points,
> increased sensitivity and precision, and exceeded StringTie3 on that cell under the NanoCount expression contract. It can also
> rescue incorrect annotated structures whose distinguishing upstream region is unobserved, as demonstrated by merge/skip/IR
> sentinel false positives."

不要说"片段 read"是既定身份、不要说"junction 未被 jitter"是保证、不要说"先验在真实 dRNA 上 90% 正确"（那是离线选定队列里与伪真值的一致率）。

## 6. 文件

- 代码：`fin/candidates/discovery.py`（`chain_compatible_read_depths`，`compatible_read_depth` 字段）、
  `fin/scoring/m2_junction_nll.py`（`unjittered_gtf_columns`）、`fin/pipeline/selection.py`（豁免块 + 审计 helper）、
  `fin/pipeline/runner.py`（审计）、`fin/pipeline/config.py`、`fin/cli.py`。
- 测试：`tests/unit/test_support_gate_attest.py`（34）、`tests/unit/test_selection_audit.py`（11）；全套 1647 通过。
- 脚本：`human_gap_dissect.py`、`human_lost_funnel.py`、`selection_audit_collect.py`、`gate_drop_precision.py`、
  `attest_diagnose.py`、`arm_compare.py`（按文库分块 bootstrap、缺格清单、显式家族映射）。
- 产物：`_gencodegrid/_headtohead/`（gap / funnel / fates / gate_drop_precision / attest_diagnose 表）、
  `_gencodegrid/gcaudit`（带审计的 full 单元）、`gcatt3*` / `gcatt0*`、`_sirv4grid/att3*`、`_heya8grid/att3*`。
