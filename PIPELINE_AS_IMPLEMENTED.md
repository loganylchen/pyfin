# PyFIN 现行算法：按实际执行顺序的逐步说明（2026-09-03）

本文只描述**代码里真正在跑的东西**。每一条都有 `file:line` 出处（以 2026-09-03 的 live 工作区为准，HEAD = ca82319 + 未提交改动）。
由 Codex（gpt-5.6-sol, xhigh）逐行核对、我再抽查关键论断后整理；原始核对稿在 `.omc/codex_walkthrough_raw_2026-09-03.txt`。

## 0. 标签约定与两条命令

每个组件都归入且只归入下面一类：

| 标签 | 含义 |
|---|---|
| **[A] 默认路径** | `--profile sirv --no-gpu` 什么都不加就会跑 |
| **[B] 已进入算法，需 opt-in** | 代码在、有 flag、默认关；会注明"当前最佳命令"是否打开 |
| **[C] 已落地但未接线** | 模块和测试都在，但 runner 里没有任何调用 |
| **[D] 只在实验脚本里** | `experiments/` 下，`fin/` 从不 import |
| **[E] 只是设计/讨论** | 没写成代码 |

两条参考命令：

- **plain SIRV**：`--profile sirv --no-gpu`
- **当前最佳 SIRV**（HANDOFF 里 91.87 F1 那条）：plain 之上加 `--min-gtf-abundance 0 --no-floor-gtf-abundance --require-gtf-positive-mass --complete-chain-rescue-reads 2 --junction-election-mode require --junction-election-burden defense`。**它没有打开** containment-evidence、junction-correction、endpoint-refine、fusion、catalog-mode。

SIRV profile 解析出的工作点（`config.py:16-33`）：`min_abundance=3.0`（含等号）、`min_gtf_abundance=1.0`、`floor_gtf_abundance=True`、`min_isoform_fraction=0.01`、`post_selection_refit=True`、`max_soft_mass_ratio=2.0`、`min_fulllen_fraction=0.1`、`min_polya5p_reads=0`、`canonical_gate=True`、`novel_junction_min_reads=2`、`containment_cluster=True`、`m2_cluster_recheck=True`、`m2_metric=auto`。`auto` 在 GTF 有 ≥2 条转录本时解析为 `mean`，否则（p00）为 `summed_llr`（`config.py:74-95`）。

## 1. 一页总览：真实执行顺序

```
BAM(+GTF) ──> intervals（按 chrom+strand 单链合并；GTF 转录本各成一个 interval 再与 read interval 合并）
   └─ 每个 interval：
      discovery
        ├─ 抓 read，CIGAR 的 N 直接变成 intron chain（无任何校正）
        ├─ [A] evidence scan（indel/junction 计数；无条件执行，写 junction_support.tsv；只记录不决策）
        ├─ [B, 未提交] junction_correction：逐 read 用 krill 改 chain（默认关）
        ├─ cluster_families：按"完全相同的 intron chain"分组 → wobble/cassette/containment 单链接 union-find
        ├─ [A, 2026-09-03/04] GTF 多外显子 chain 作为普通 0-read variant key：先只用 read chain 做 union-find，再把纯注释 chain 挂到**至多一个**已有分量（不可桥接）；可作 fold 容器，**不作 foldee**
        ├─ collapse：exact contiguous subchain 折进最长容器（span_guard 现在按 CIGAR 判断 read-through）
        ├─ mono read 留在 holding bucket（m2_em 下延后到 EM 之后处理）
        └─ _emit：key 若是注释 chain → 发射该 chain 上的全部 GTF 模型（注释几何 + key 的 read 作为支持）；否则 novel 候选，start/end = 成员 read span 的并集
      pre-EM gates：canonical gate（默认开，只删 novel）；junction dominance（默认关）
      assignment：每个候选一个 mappy aligner → AS 矩阵 → 严格 AS 并列集 → 只在并列集里算 M2 junction NLL
                  → d_tx（0 / −log p / MISSING）→ EM（sigma=1, beta=0, tol=1e-4 收敛）→ R
                  → hard assignment 由**显式规则**决定（`assignment.hard_assign`：最大后验 → atol=1e-12 并列集 →
                    优先 annotated → 再按 candidate_id；**不依赖候选列表顺序**）
      select_m2_interval（顺序执行，所有统计都读同一个 pre-selection R；不重跑 EM）：
        [B] containment-evidence → [B] junction election → cluster recheck → novel junction support
        → [B] guided support → m2 support gate（[B] rescue 豁免）→ containment cluster
        → [B] lever-1 collapse → mono_resolve_post_em（硬分给 EM 丰度最高宿主）
   跨 interval 聚合（read 去重、丰度按重复观测修正）
   [B] finalize_containment / merge_election_shards / merge_correction_shards（只合并审计表）
   select_global：A4 丰度 floor（[B] P0 严格 >0）→ isoform fraction → soft/hard 比 ([B] rescue 豁免)
                  → mono 规则 → fulllen fraction → polya5p → [B] junction snap
   [B] endpoint refine → refit（每条 read 在幸存者上重新归一，一次）→ [B] refit 后 P0 二次清扫
   [B] catalog completion → gene id → writers
```

## 2. 逐步说明

### 2.1 interval 构造 — `fin/io/interval_manager.py`

- **[A]** 全 BAM 扫描；跳过 unmapped 和 secondary，**不看 MAPQ**，不看 duplicate/QC flag（`:224-233`）。supplementary 先留着用于识别 fusion。
- **[A]** fusion read = supplementary 或首/尾 soft-clip ≥ 250 bp；不参与普通 interval 构造，存入 fusion side list（`:135-159, :243`）。
- **[A]** 普通 read 按 `(chrom, strand)` 分桶，strand 直接取 `is_forward`；每条 read 是 `[reference_start, reference_end)`（`:174-191`）。
- **[A]** 桶内按 start 排序，下一条 start ≤ 当前 end + `max_gap`（普通配置 `max_gap=0`）就合并，边界取并集（`:313-319`）。
- **[A]** 有 GTF 时，**每个 GTF 转录本单独生成一个同链 interval**，再与 read interval 合并（`:370-446`）——注释可以扩张或连接 read interval，但不能跨链。合并后 `read_count==0` 的纯 GTF interval 被删除（`:465`）。
- **[B]** fusion 子流程默认关（`fusion_enabled=False`，`config.py:743`），当前最佳命令也没开。

### 2.2 discovery — `fin/candidates/discovery.py::_chain_cluster_candidates`

**read 与 chain**
- **[A]** 对 interval 做 region fetch，跳过 unmapped/secondary，read strand 必须等于 interval strand（`:626-649`）。read 序列只从 primary 非 supplementary 比对取（`:659`）。
- **[A]** intron chain = CIGAR 里的 `N` 段 `(pos, pos+len)`，M/D/=/X 消耗参考，S/I/H/P 不消耗（`intron_chains.py:18`）。**没有任何 snap、canonical 校正或端点约束**。read span 记为 `[reference_start, reference_end)`（`:331`）。

**[A] evidence scan**（`:340-366`；2026-09-03 起**无条件**执行，此前只在 rescue / containment / election / correction 任一打开时才做）
- 只收 MAPQ ≥ 1 的 read；D 记参考区间、I 记参考点、N 记 junction。junction "indel-clean" = donor 和 acceptor 两侧各 10 bp 内没有任何 I/D。
- `junction_clean_depths`（`:199-227`）：对候选的 exact complete-chain read（成员 read 减去 folded 的 read）逐 junction 数 clean 支持；`weakest depth` = 最小值，记在候选的 `complete_chain_weakest_indel_clean_depth` 上。
- interval 级注册表 `{(donor, acceptor): (raw, clean)}` 记在 `CandidateSet.junction_support`，并在每次运行结束时写成 **`junction_support.tsv`**（`fin/analysis/junction_support.py`）。
- **只记录，不决策**：谁用它仍然各自 opt-in（rescue [B]、election [B]、correction [B]）。已验证 `full/p50/p00` 三个默认格子的 `assembly.gtf`/`scores.tsv`/`abundance_refit.json` 与 `_sirv4grid/base` 字节一致。

**[B, 未提交] junction correction**（`fin/analysis/junction_correction.py`，hook 在 `:368`，`cluster_families` 之前）
- runner 通过 `_make_junction_correction` 传入 callback（`runner.py:1682`），krill 上下文每个 worker 建一次。默认 `off`，当前最佳命令**未开**。
- registry：每个 junction 的 raw 数、clean 数、read id（`:107`）。提名：备选与本 junction 双坐标 ≤ 15 bp，clean ≥ max(3, 3×own_clean)（`:134-157`）；**只是资格筛查，不决定结果**。每 interval 最多处理 400 条被提名 read（`:58`）。
- ballot：read 自己的完整 chain vs 只换这一个 junction 的 chain，两条序列从基因组拼出（负链 revcomp）（`:187`）；先对 own 跑一次拿 `init_scale/init_shift`，再用**同一固定缩放**跑两个假设（`:273-277`）；krill 不返回这两个字段就报错（`:310`）。打分窗口 = 变化坐标 ±8 bp 在两个假设下映射到的采样区间并集；每 event 代价 `0.5·z² + log σ`，`σ = model_stdv·√var`，不可建模罚 6.0；每假设 ≥ 3 个 event（`:215-230, :304`）。
- 判定：`margin = own_nll − alt_nll > junction_correction_margin`（默认 0），多个备选取 margin 最大；一条 read 的所有胜出替换一起应用，组合 chain 必须有序、不重叠、无重复（`:407-421`）。`audit` 只写表、返回原 chain；`require` 才返回校正 chain。校正后的候选**不享受任何门的豁免**。
- 实测（表达契约，6 样本）：jitter/flip 24 格 dF1 −0.17（CI 跨零）；p00/full/p50 30 格 +0.02（CI 跨零）。**没有增益**。

**cluster_families**（`fin/candidates/chain_cluster.py`）
- **[A]** 默认 `chain_cluster_discovery=True`、`clustering="families"`（`config.py:149-150`）。
- **[A]** 第一步按**完全相等的 intron tuple** 分组；5' 端、3' 端、比对端点一概不参与（`:473`）。空 chain 的 read 进 mono bucket。
- **[A]** `_related(a,b)`：同 intron 数且每个 donor/acceptor 差 ≤ 6 bp → wobble；长链比短链多一个 intron → cassette（被跳过外显子 0 < 长度 < 70）；否则 containment（`:144-198`）。
- **[A]** **`_contains` 接受长链中任意 offset 的连续子窗口**，每个 junction ±6 bp，**不要求是链向 3' 后缀**（`:152`）。
- **[A]** 单链接 union-find：A~B、B~C 即把 A、B、C 放进同一 family，即使 A、C 不相关（`:486-498`）。family 的 read_pool 是各 exact-chain read 集的并集；clustering 阶段不重新分配 read（`:508`）。

**GTF 作为普通 variant（2026-09-03 改，默认、无开关）**
- **[A]** interval 内同链的 GTF 转录本先按 chain 分组成 `gtf_models`（同一 chain 可有多个模型：不同 TSS/TES）；多外显子 chain 在 `cluster_families` 里作为 **0-read variant key** 与 read chain 进入同一个 union-find（`chain_cluster.py::cluster_families(gtf_variants=…)`）。read chain 与注释 chain 相同时就是同一个 key（不再有"第一个 GTF 赢"的事后合并）。
- **[A]** 注释 chain **不能桥接**两个 read family（2026-09-04 改）：union-find 第一遍只跑 read-derived chain，第二遍把每条纯注释 chain 挂到至多一个已有分量（相关 variant 最多者，再按坐标），注释与注释之间也不连边。理由：零读数（甚至被 corrupt 的）假设不应改变**数据**产生的图连通性——那会连带改写 collapse 的容器集合、EM 竞争范围和 isoform-fraction 分母。它仍是普通 variant：可作 collapse 的容器（read 子链折进注释 chain，read 成为它的 discovery 支持）；**但注释 key 永远不是 foldee**：fold 是"这些 read 是容器的 5' 降解"这一 read 级判断，真正的 5' 截短注释异构体的 read 永远无法反驳它（它们必然从多出的 intron 下游开始），0-read 注释 key 更无从判断——所以注释 key 一律保留为成员，存在与否交给 EM + 最终 GTF 门（P0 要求正质量）。这是对 Codex "纯 GTF foldee 无条件折叠" 规则的**有意偏离**，由 156 格实测决定：按 Codex 规则 plain default 85.58 → 82.55 F1（Sn −3.8，137 down / 3 up，损失随注释比例单调增大，p00 +0.15），丢失的全部是 SIRV303/503/606/608/609/615 这类"更长注释异构体的精确子链"（记录：`_sirv4grid/gv_off_foldgtf`）。
- **[A]** `_emit` 对幸存的注释 key 发射该 chain 上的**每个** GTF 模型（保留注释 start/end/序列，`supporting_read_ids` = 该 key 的 read）；family_id 由 family 全部 variant chain（含注释）计算。mono GTF 仍按 genomic overlap 与 mono locus 匹配（`discovery.py::_best_overlap_gtf`）。
- 影响：guided 条件（full/p50/c_*）输出会变，需要重跑 grid 基线；p00 无 GTF，字节不变。

**collapse**（`chain_cluster.py:533-621`）
- **[A]** 每个 family 单独 `collapse(span_guard=True)`；只折叠**严格 exact contiguous subchain**，不折 wobble/cassette/非精确 containment。
- **[A]** 链按 intron 数降序、read 数降序、坐标排序；短链折进第一个包含它的已保留长链。
- **[A]** span_guard（2026-09-03 改为 CIGAR 级，`_cigar_exonically_covers`）：对容器多出来的每个 intron `[d,a)`，取窗口 `W=[d−10, a+10)`；短链 read 只有在 **(1)** M/=/X/D 足迹覆盖整个 W、**(2)** W 内没有 N、**(3)** W 内单个 D ≤ 3 bp 且 D 总长 ≤ 6 bp、**(4)** 上下游 10 bp 外显子 flank 各有 ≥ 8 bp M/=/X、**(5)** intron 体内 M/=/X ≥ (a−d)−6 时，才算"真的读穿了这个 intron"而保留在短候选；否则折入。旧规则只看 `reference_start/end` 是否夹住 intron，大 deletion / soft clip / N 都会被误判为读穿。图示：`experiments/prod_validation/_sirv4grid/fig/span_guard_explainer.png`。
- **[A]** 被折叠的 read 加进容器 `read_ids`；原短链作为 `GenCandidate.folded` 影子来源保存，不增加 read mass；后续仅供 5'-TSS recovery / ledger 用（`:614-621`）。

**mono read**
- **[A]** `mono_resolve_post_em=True` 且 `quant_mode=m2_em` 时，runner 让 generation 时的 `chain_cluster_fold_monoexon` 不生效；mono read 留在 holding bucket，EM 之后处理（`runner.py:584`）。
- **[B]** 非 m2_em 或关掉 post-EM 解析时，才在 generation 时把完全落在某候选外显子内的 mono read 折入确定性最佳宿主（`discovery.py:387-391`）。当前最佳命令不走这支。

**_emit**（`discovery.py:460-537`）
- **[A]** 候选在列表中的顺序**不再承载任何语义**（2026-09-04）：hard assignment 由 `assignment.hard_assign` 的显式规则决定。此前 `R.argmax` 用最小列号打破精确并列，而 discovery 恰好先追加 GTF，于是"并列时判给注释"这条 assignment 策略被藏在 discovery 的列表顺序里（实测：仅改变发射顺序就使 SIRV4 156 格 Pr −0.93）。
- **[A]** 多外显子成员的 novel 候选 **start/end = 成员 read 的 `min(reference_start)` 和 `max(reference_end)`（并集）**（`:503-511`）；不是中位数、不是代表 read、不是 3' 锚定共识。
- **[A]** 成员 key 是注释 chain → 该 chain 上的每个 GTF 模型都发射，read 作为支持并入，GTF 的 start/end/序列不改；0-read 的注释 key 也发射（interval 内同链的每个 GTF 转录本都一定成为候选）。
- **[A]** 否则用并集 start/end + chain 从基因组拼序列（负链 revcomp）；空序列不生成候选（`:143-156, :481`）。
- **[A]** mono bucket 先按 span 重叠单链接分 locus，每 locus 取并集 span，优先 overlap 匹配一个 mono GTF，否则生成 novel mono（`:58, :533`）。
- **[A]** `min_novel_reads` 只删支持 read 数不够的 novel；默认 1，所以不删单 read 候选；GTF 不受影响（`:546`；`config.py:141`）。

### 2.3 pre-EM gates — `runner.process_interval`

- **[B]** fusion augmentation：`fusion_enabled=False`，两条命令都没开（`runner.py:623`）。
- **[A]** canonical gate（`canonical_gate=True`）：只删**novel** 多外显子中含非 canonical junction 的模型，允许 GT-AG / GC-AG / AT-AC；GTF、fusion、mono 豁免；缺基因组序列时 fail-open（`selection.py:840-850`）。删候选但**不删 read**，这些 read 仍在 assignment 里竞争其他候选（`:843`）。
- **[B]** junction dominance：默认关，两条命令都没开（`config.py:671`）。

### 2.4 assignment — `fin/pipeline/assignment.py::Assigner.assign` + `tie_nll`

- **[A]** 候选顺序固定；read 轴 = 候选集 read id 排序后且有 primary 序列的 read；无序列的 read 在任何 mappy/EM 之前就被删（`:583`）。
- **[A]** **每个候选单独一个 `mappy.Aligner(seq=candidate.sequence, preset=map-ont)`**，不是 multi-index（`:598-604`）。每个 read×候选格取该 aligner 所有 hit 中 `score_hit` 最大者；`score_hit = 2·mlen − 4·mismatch − Σ min(4+2k, 24+k)`，任一单个 I/D > 50 bp（`M1_MAX_INDEL_BP`）该 hit 作废（`mappy_score.py:15-31`）。
- **[A]** read 至少要有一个格 AS > 0 且 ≥ `MAPPY_R1_MIN_AS`（默认 0），否则整条 read 不进 EM（`:625`）。EM 前删 read 的原因只有：无 primary 序列，或对所有候选都没有 AS > 0 的比对。
- **[A]** **tie set = 该行所有有限 AS 中 `AS ≥ best − 1e-9` 的格**，即严格最佳并列（`:334`）。
- **[A]** **M2 只对 `len(tie) ≥ 2` 的并列格构造信号 payload**；不在 tie 里的候选永远不会因 M2 得到这条 read 的似然（`:336-357`）。对每个 tied 格取最高 AS hit，只 eventalign 候选序列的 `[r_st:r_en]` 切片（`:361-371`）。任一 tied 假设切片含非 ACGTU 字符 → 整条 read 的 M2 弃权，保留完整 tie（`:46, :375`）。
- **[A]** `_m2_score_region`：`mean`（有注释）用 class junction window，flank 2 bp、每侧骨架 k=10（`:71-91`）；`summed_llr`（p00）要求 tied 候选 intron 数相同否则弃权，窗口是 differing-junction windows，SIRV profile 把 flank 覆盖为 4（`:78-86`；`config.py:31`）。event 代价 `0.5·z² + log(model_stdv)`，mean 除以有效 event 数，summed 不除（`m2_junction_nll.py:1204-1238`）。**主流程的 M2 不用 `√var` 二次缩放**（junction_correction 用了）。
- **[A]** d_tx（`m2_diff_cover_gate=True`，`config.py:474`）：唯一 AS 最佳 → 该格距离 0；tie 行若两假设 NLL 差 ≥ margin（`mean` 用 `m2_diff_cover_margin=0.5`，`summed` 用 profile 的 `m2_summed_llr_margin=1.0`；**不是 `m2_tiebreak_margin=1e-9`**，那个给 election 用）→ 只给最低 NLL 候选距离 0，其余 tied 格 MISSING；若 read 事件覆盖了全部 diff span，该 winner 得一张 `covered_vote`（`:682-712`）。margin 不足或无信号的 tie read **不删**：按 tie 候选的 covered-vote 比例给 `d_tx = −log p`，零票候选 MISSING；整组无票则全部设 0（=1/K）（`:719-726`）。
- **[A]** EM：`d_tx` 外的格是 MISSING；read×read 距离全 0；调用通用引擎时**硬编码 `sigma=1.0, beta=0.0`**（配置里的 `em_sigma` 主路径不读）（`:747-753`）。初始 `R_ij = exp(−(d_ij − row_min)/σ)` 行归一化（`assignments.py:196`）；**每轮检查 `max|R − R_new| < tol` 收敛，`max_iter=1000`、`tol=1e-4`，不是固定次数**（`assignments.py:240`；`config.py:370-373`）。hard assignment = 行 argmax；丰度 = 列和；hard read 数 = argmax 命中数（`assignments.py:269`；`quantification.py:112`）。
- **[A]** EM 后的 `krill_tiebreak` 二次覆盖默认关（`config.py:429`）。

### 2.5 interval 级选择 — `select_m2_interval`（`selection.py:57-438`）

**共同语义 [A]**：R、tie set、NLL、quant_results 全部来自同一次 assignment；各门按顺序把要删的列累积进 `drop_cols`，需要时用 `exclude=drop_cols` 让已删候选不参与后续关系，但**从不删 R 列重跑 EM**；最后按 id 从 quant_results 删掉。也就是"顺序 gate + 原始责任矩阵统计"。

按执行顺序：

1. **[B] containment-evidence layer**（`:86-130`）：默认关，当前最佳未开。枚举 short-in-long 对（own_tes / tss_only / equal_chain / mono_span，TES 距离 25 bp），parent 的独有 junction 最小 clean depth ≥ 2 或 ≥ 2 条延伸到 short 5' 上游的 read 才允许反驳；结构性 TES 判定可即时 drop/exempt，hazard 判定等全局 FDR；弃权不动。
2. **[B] junction election**（`:132-165`；`junction_election.py`）：默认关，**当前最佳 = require + signal + defense**。提名：accused clean ≤ 1，competitor 双坐标 ≤ 15 bp 且 clean ≥ max(3, 3×max(accused_clean,1))；**不检查 raw 支持**（`:74, :168`）。每个 junction 的代表 = 携带它的候选中丰度最大者。共享 champion 的 accused 合并成 K-way contest，每 interval ≤ 50。窗口 = donor 范围 ±8 ∪ acceptor 范围 ±8；代表在窗口内除自身 junction 外的 exon 边界指纹必须一致，否则 `window_impure`。**选民只取每个代表候选的 `assigned_read_ids`（P1 bug：漏掉同 junction 其他携带者的 read）**，span 须覆盖全部窗口，按 SHA1(read id) 排序截到 30（`:448-461`）。每张票 = `read_cand_mean_nll`，两假设各 ≥ 3 event，runner-up − winner ≥ 1e-9 才有效。quorum 5、supermajority 0.75：accused 自己是 supermajority winner → `upheld`；别人 winner 且 accused ≤ 0.25 → `defeated`；有效票 < 5 → `no_quorum`；其余 `split`（`:496-509`）。**prosecution 只删 `defeated`；defense 除 `upheld` 外全删——所以当前最佳命令下 `no_quorum`、`split`、`no_ballots`、`window_impure` 都会删掉被告候选**（`:515`）。
3. **[A] cluster recheck**（`:167-226`；`m2_junction_nll.py:420-601`）：所有多外显子候选按同 intron 数、每 junction ≤ 20 bp 或 cassette < 70 bp 做连通分量；分量内原始 `R.sum(axis=0)` 最大者为 anchor；`候选/anchor < 0.15` 才进入 shadow 处理；novel sibling 直接删，fusion 不删，GTF sibling 只在允许 novel 顶替且其区分 junction 没有 ≥ 1 条精确 read（tol 0）时才删。
4. **[A] novel junction support**（`:228-245`）：SIRV 门槛 2；只查 novel 多外显子；任一 junction 在同链 CIGAR 观测图中 ±2 bp 支持 < 2 → 整条删。GTF/fusion/mono 豁免。
5. **[B] guided support gate**（`:247-268`）：`guided_junction_min_reads=0`，两条命令都没开。
6. **[A] m2 support gate**（`:270-319`；`m2_junction_nll.py:366`）：GTF 和 novel 多外显子都查；候选必须至少是某条 read 的 M1 sole-best，或某条 AS-tie read 的最低 M2 NLL（`tie_ok=True`，并列都算）。**[B] rescue 豁免**（当前最佳 K=2）：GTF weakest clean depth ≥ 2；novel ≥ 2 junction 时 ≥ 3；novel 单 junction ≥ 10；fusion/mono 不救（`:460-485`）。
7. **[A] containment cluster**（`:321-351`；`m2_junction_nll.py:812-899`）：只删 novel 多外显子 short；**short 可匹配长链任意连续子窗口**（每 junction ±6），不要求 3' 后缀；同时满足 `short_EM ≤ 0.3·parent_EM`、`short_reads ≤ 0.3·parent_reads`、short raw read ≤ 10 才删。**动作是 drop，不把 mass 折给 parent**。
8. **[B] lever-1 containment collapse**（`:353-408`）：默认关。这是**唯一**实现了"严格链向 3' 后缀 + 3' 端差 ≤ 20"的地方；动作是把 short 的 read 和 soft 丰度加到 parent 再删 short。
9. **[A] mono_resolve_post_em**（`:410-432`；`m2_junction_nll.py:741-803`）：每条 mono hard-assigned read，找同链、span 落在候选 `start/end ±10` 内且不与任何 intron 相交的多外显子宿主；一个宿主全量折入；**多个宿主取原始 EM 丰度最大者**（tie 取 id）；折入的 read 给宿主丰度直接 +1；未覆盖的留在 mono，剩余 < 2 条则删，否则 mono 丰度重置为未覆盖 hard read 数。**不是按似然或比例软分配。**
10. **[B] endpoint refine**：默认关，而且**不在这里**——它在跨 interval 聚合和全局选择之后、snap/refit 之前（`runner.py:265`）。

### 2.6 跨 interval 聚合 — `quantification.aggregate_across_intervals`（`:544-680`）

- **[A]** 同 id 的 hard read id 做并集；丰度先跨观测求和，再乘 `unique_reads / Σ max(interval_reads,1)`（按重复观测修正，不是简单相加）。
- **[A]** confidence 等按 hard read 数加权；max_R、fulllen_frac、complete-chain weakest depth 取最大；containment verdict 按 `contradicted > supported > unidentifiable`（`:16`）。

### 2.7 全局选择 — `select_global`（`selection.py:565-835`），按顺序

1. **[A] A4 floor**：plain SIRV novel ≥ 3（含等号），GTF floor = max(1, 3) = 3；fusion 豁免（`:590-626`）。**[B] 当前最佳**：GTF 数值 floor 为 0，但 **P0 要求 abundance > 0**（`:605-631`；`gtf_zero_mass` 原因）；novel 仍 ≥ 3。
2. **[A] isoform fraction** 0.01：只删 novel 多外显子；分母 = 同 family 内所有来源的最大丰度（无 family 时退回同链重叠）；自己是最大者绝不删（`quantification.py:160-221`）。
3. **[A] soft-mass ratio** 2.0：只查 novel 多外显子；hard=0 必删；`soft/hard ≥ 2.0` 删（`quantification.py:240`）。**[B]** rescue（同 2/3/10 规则）和 containment `supported` 可豁免（`:710`）。
4. **[B] mono 规则**：`drop_mono_exon_novel=False`，两条命令都不执行。
5. **[A] fulllen fraction** 0.1：只查 novel 且 ≥ 2 intron；read 的 5' 和 3' 都与候选差 ≤ 25 bp 才算全长；≥ 4 条有 span 的 read 才计算，否则 −1 豁免；fraction < 0.1 删（`quantification.py:341-374`）。
6. **[B] polya5p**：SIRV 解析为 0，关。
7. **[B] junction snap**：默认关，两条命令都没开（`config.py:570`）。
8. **[A] refit**（`abundance_refit.py:143-327`）：**所有存在性决定都冻结在 refit 前的丰度上，refit 不重跑任何门**（`:338`）。ledger 里同一 read/候选跨 interval 的权重先相加；删掉非幸存者列后，每条普通 read 对幸存权重做 `p_j = w_j / Σ w`；没有幸存列的 read 成为显式 orphaned mass（不会硬塞给谁）。mono 强制 read 不按比例归一，而是选幸存目标中丰度最大者、概率 1（`:218`）。总 mass 守恒 `assignable = assigned + orphaned`，否则报错（`:284-327`）。
9. **[B] refit 后 P0 二次清扫**（当前最佳开）：refit 把已过 A4 的 GTF 归零则再删（`runner.py:353`；`selection.py:512`）。
10. **[B] catalog completion**：已接线（`finalize.py:73`），默认 off，当前最佳未开。
11. **[A]** gene id 解析 + writers（`finalize.py:94-103`），不再删候选。

### 2.8 runner 收尾的 sidecar

- **[B]** `finalize_containment`（mode 非 off 才跑；当前最佳不跑）：合并 shard、hazard 行做 grouped FDR、require 下执行最终 drop/exempt。**已知 bug**：最后用 `{short_id: verdict}` 字典推导回写，同一 short 多行时排序最后一行获胜，没有用聚合阶段的优先级（`containment_evidence.py:489`）。
- **[B]** `merge_election_shards`（当前最佳跑）：只合并 `junction_election.tsv`，drop 早在 interval 阶段完成。
- **[B, 未提交]** `merge_correction_shards`：只合并审计表。

## 3. 未接线 / 只在实验里 / 只是设计

**[C] 已落地但 runner 不调用**
- `fin/scoring/candidate_align.py`：`align_reads_to_candidates`（整个 family 一个 multi-sequence mappy index，输出 AS / q_cov / event identity / soft clip / 3' 距离）和 `is_compatible`（q_cov ≥ 0.80、identity ≥ max(0.80, best−0.05)、3' 距离 ≤ 50）。只有 `tests/unit/test_candidate_align.py` import。
- `fin/candidates/diff_regions.py::cluster_diff_regions`：候选间差异区域 + 共享 context。只有单测和实验 import；主流程用 `m2_junction_nll.py` 自己的窗口函数。
- `fin/candidates/explore.py::explore_family_paths`：family 内 intron DAG 的 beam 路径枚举（beam 256、最多新增 512）。只有测试调用。

**[D] 只在 `experiments/` 里**
- `junction_clean_probe.py`：TRUE / NEAR_FAKE / FAR_FAKE junction 的 raw/clean 统计；**"raw×3 ≤ competitor raw" 守卫只在这里，production 的 `nominate_contests` 没有**。
- `junction_signal_probe.py` + `junction_signal_plot.py`：同一 read 两个 junction 假设的 krill 对比；**按采样点的似然 `Σ n·[log σ + ((m−μ)²+s²)/(2σ²)]`、HMM 模式（K/B 行过滤）、run 依赖性审计（dominant-run fraction、leave-one-run-out）都只在这里**。
- `experiments/wobble_heya8/*`（collect_evidence、validate_signal、read_vote 等）：全部离线研究，方向是 experiments→fin。

**[E] 只是设计**
- `.omc/codex_consult_step5_resolver_2026-09-02.txt` 及其中文版：conflict component、CONFIRMED/UNRESOLVED/REFUTED 三态、K=3 witness、固定 25 次 EM 的统一 resolver。live `fin/` 没有这些状态或执行块。registry 保存 read id 的要求：live `CandidateSet.junction_support` 只存 `(raw, clean)` 计数（`discovery.py:568`）。

## 4. 与你的五步心智模型逐条对照

**你的第 1 步"isolated intervals，GTF 进入 interval"** —— 基本成立，但有两点不同：
1. GTF 不只是"进入 interval 提供候选"，它**先参与 interval 构造**（每个转录本一个 interval，可以扩张/连接 read interval），之后在 discovery 里**作为普通 0-read variant** 与 read chain 一起聚类/折叠/发射（2026-09-03 起；之前是事后"挂载"，不桥接不折叠）。
2. interval 构造不看 MAPQ；evidence scan 才看 MAPQ ≥ 1。

**你的第 2 步"junction 枚举（interval 级注册表）"** —— **不是一个独立阶段**。默认路径直接按 read 的 raw chain 聚类，没有 junction registry；interval 级的 raw/clean 注册表只在 rescue/containment/election/correction 任一打开时作为 evidence scan 出现，只有 correction 和 election 消费它（`discovery.py:340`；`runner.py:594`）。

**你的第 3 步"聚类：wobble 合并、少 5' exon 可合并、少 3' exon 不可合并"** —— wobble/cassette 成立；**"3' 必须完整"没有实现**：
3. `_contains` 接受长链**任意 offset 的连续子窗口**，没有链向 3' 后缀检查（`chain_cluster.py:152`）；exact-chain 分组也**完全忽略端点**（`:473`）。
4. 单链接 union-find 会把不直接相关的链拉进同一 family（`:486-498`）。
5. novel 候选的 start/end 是成员 read span 的**并集**（`discovery.py:503-511`），不是 3' 锚定，也不是代表 read。

**你的第 4 步"cluster 内 collapse"** —— 成立，但：
6. 折叠的是**任意** exact contiguous subchain（前缀、后缀、中间窗口都折），不只是 5' 降解梯子；只有默认关闭的 lever-1 collapse 实现了严格 3' 后缀（`chain_cluster.py:123`；`m2_junction_nll.py:680`）。
7. ~~span_guard 只看 genomic span~~ → 2026-09-03 起按 CIGAR 足迹判断读穿（见 2.2 collapse）。

**你的第 5 步"补全候选集 → mappy+M2 全体分配 → 过滤 → 重新分配被删候选的 read → 识别 5' 子集"** —— **这是差异最大的地方**：
8. **候选集没有"补全"阶段**：候选 = read 实际产生的 chain + GTF 假设；`explore_family_paths` 没接线。唯一已接入的补全机制是 opt-in 的 per-read junction correction（只替换被提名的 junction，不做路径补全），且实测无增益（`discovery.py:381`；`explore.py:80`）。
9. **"mappy + M2 对完整候选集评分"只有一半成立**：mappy 确实评估所有格；但 **M2 只在严格 AS 并列集里算**，AS 非最佳的候选即使信号更好也拿不到这条 read 的 M2 似然（`assignment.py:334-357, :638`）。M2 不是通用 per-read 似然：唯一 AS 赢家直接距离 0，非 tie 格 MISSING。
10. 主分配的判决 margin 是 0.5（mean）/ 1.0（summed），**不是 `m2_tiebreak_margin=1e-9`**（那个给 election 用）（`assignment.py:682`）。
11. **EM 用 tol=1e-4 收敛（上限 1000），不是固定次数**；sigma、beta 在主路径硬编码为 1 和 0（`assignment.py:747-753`；`assignments.py:240`）。
12. **"过滤 → 重新分配被删候选的 read → 再过滤"这个循环不存在**：interval 门只是累积 drop 列并持续读**同一个** pre-selection R；全局门也一样；只在**全部**门结束后做一次 refit 把每条 read 在幸存者上重新归一（`selection.py:71, :438`；`runner.py:336`；`abundance_refit.py:338`）。
13. **丰度大量参与存在性判定**：cluster recheck 用 `<0.15` 丰度比删 shadow，containment cluster 用 0.3 丰度比，A4 要求 novel/GTF 达 floor，isoform fraction 和 soft/hard 比也直接决定候选是否存在（`selection.py:184, :590`；`m2_junction_nll.py:901`；`quantification.py:235`）。
14. **mono read 不是按似然分配**：硬分给覆盖它的幸存宿主中原始 EM 丰度最高者，宿主丰度直接 +1（`m2_junction_nll.py:795-800`）。
15. **election（当前最佳的 defense 模式）把弃权当有罪**：`no_quorum`、`split`、`no_ballots`、`window_impure` 全部 drop（`junction_election.py:515`）。今天的探针证实了一个真实受害者：SIRV3 8966/8975 这对真实可变供体在 ~46 个格子里因 `no_quorum` 被删。
16. election 的选民只取代表候选的 read（P1 bug）；提名只看 clean 不看 raw（raw×3 守卫只在实验脚本里）（`junction_election.py:450, :170`）。
17. "识别 5' 子集"在 live 代码里对应的是：collapse 的 folded 影子 + **[B] `--folded-tss-recovery-mode audit|require`**（2026-09-03 新增，默认 off）：m2_em 第一遍 EM 后，对每个 parent 的每个纯 5' 截断 shadow（+ 链 exact suffix / − 链 exact prefix，≥3 read，≥80% read 的 3' 端在 parent 3' ±25 bp 内），用 shadow read 的 5' 端中位数提出 TSS，对第一遍 `R[read,parent]≥0.05` 的 read 5' offset 做 `evaluate_internal_tss`（bin 25、local background、Bonferroni 到 parent 全长 bin 数）；`supported` 才以 shadow 自身 chain + bin 内 read 5' 中位数 + parent 3' 端重新发射为 novel 候选，然后整个 assignment/EM **重跑一次**，选择只看第二遍。记录在 `tss_shadow_recovery.tsv`。
    - **为什么还需要一个 5' 位置似然项**：shadow read 对 parent 和恢复出的短模型的序列似然完全相同（同 AS、无 junction 差异），第二遍 EM 只按丰度分 → 短模型拿不到任何 hard read，soft/hard 比门随即删除它（实测 SHARP 例：supported，p=0.028，二次 EM 后短模型 soft 7.0 / hard 0）。按 Codex 的方案 2（`.omc/codex_consult_tss_5p_likelihood_2026-09-03.md`）实现了 `assignment.apply_start_likelihood`：**只在第二遍**、只对"AS 并列集恰好是 {parent, recovered}"的 read，按 read 5' 端在 parent 转录本坐标中的 bin 加一项冻结的位置 NLL——落在 supported bin 内：parent 格 += −log h_P(b0)（h_P 是第一遍测出的 parent 在该 bin 的终止 hazard，SHARP 例 0.019 → 3.96），短模型格 +0；起点在 TSS 上游：短模型格 = MISSING（短转录本不可能产生比自己 TSS 更 5' 的 read）；下游：弃权。该项在 junction sigma2 归一化**之后**叠加（不被 median 归一化压扁）。同一对比也镜像进 `nlls_by_read`，使 M2/M1 support gate 按完整证据判断"这条 read 最像哪个模型"（不是豁免：恢复出的模型仍必须是某条 read 的最佳）。SHARP 例结果：hard 14/14、短模型 soft ≈13.7、通过全部门。第一遍 assignment、off/audit 模式与旧路径逐字节相同（`test_first_pass_is_byte_identical_with_recovery_on`）。
18. 统一 resolver（conflict component、三态、K=3、固定 EM）**只是设计文本**。

**读 mass 守恒的边界**：无序列或无 AS > 0 格的 read 在 EM 前就消失；refit 只对有 ledger 行的 read 守恒；幸存列全被删的 read 成为 orphaned mass（`assignment.py:585, :625`；`abundance_refit.py:250`）。

## 5. 本轮（2026-09-02/03）新增的实证，及各自的状态

| 内容 | 状态 | 结论 |
|---|---|---|
| raw/clean junction 探针（6 样本） | [D] | 假 junction 的 read 89–97% 在 ±10 bp 内带 indel，真 junction 的也有 50%；SIRV3 8966（50 read/1 clean）被 clean-only 规则误起诉；加 `raw×3` 守卫后真 junction 误起诉 0/6 样本，假 junction 覆盖保留 97% |
| per-read 信号探针（6 样本，固定初始缩放） | [D] | 假 junction 的 read 78% 偏向真；真 junction 的 read 85% 守住；真–真对照 91% 守住自己。**但单条 read 的分数常被一个重复 k-mer run 主导**（例：+28.5 里 82% 来自一个 12-event run），需要 run 审计；HMM 模式掉到随机（原因待查：两假设 event 数不一致） |
| 按采样点似然 vs 按 event 求和 | [D] | 无差别（77.9% vs 78.0%） |
| per-read junction correction（5a） | **[B]，未提交** | 已接线、audit 字节惰性；require 在 jitter/flip 24 格 **−0.17**、p00/full/p50 30 格 +0.02，CI 都跨零，**无增益** |
| GTF 作为普通 variant key（默认、无开关，不可桥接） | **[A]，未提交** | 2026-09-03/04 实现。第一版按 Codex 规则折叠纯 GTF 子链：156 格 plain default **−3.03 F1**（Sn 83.32→79.49，Pr 88.16→86.04，137 down/3 up；p00 +0.15）→ 已改为"注释 key 不作 foldee"，156 格 `gv_off`/`gv_best`/`gv_tss_*` 重跑中（SLURM 211628） |
| CIGAR 级 collapse span_guard | **[A]，未提交** | 2026-09-03 实现（Codex 规则 1-5 + 8 类边界单测）；所有条件的默认输出都可能变，同上重跑 |
| `--folded-tss-recovery-mode` + 5' 位置似然 | **[B]，未提交** | 2026-09-03 实现；audit 字节惰性与 require 行为在 `gv_tss_audit`/`gv_tss_req` 156 格中验证；重点看 c_trunc5p*/c_merge* |
| 加法式 junction correction 候选级探针（margin 10/20 × K 2/3，margin 0 参照；6 样本 × 4 条件，结果与条件无关因为校正在 GTF 之前） | [D] **负结果** | lost_chains=0（加法语义成立）。新增候选几乎全是 FP：m10/K2 新增 128 条中 TP 12 / FP 116（每样本集 3 TP/29 FP，TP 全是同一条 SIRV6 三 intron 链，2–3 票；FP 有 5–7 票的）；m10/K3 4/52；m20/K2 **0/64**；m0 参照 44/1612。read 级：own_TP→proposed_FP 612 vs own_FP→proposed_TP 448——ballot 把真 read 改错的次数多于把错 read 改对。结论：per-read K-way krill ballot 作为候选补全机制不成立；汇总 `_sirv4grid/jc_additive_probe_summary.txt` |
| 显式 hard-assignment tie-break（去掉列表顺序依赖） | **[A]，未提交** | 2026-09-04；Sequin 104 格 ΔF1 −0.01（CI [−0.01,+0.00]），中性；目的是移除隐式策略而非提性能 |
| 第二套数据 Sequin（4×26=104 格 + 独立浓度真值） | 环境/协议 | 同一份代码只换 `--profile`：Sn 8.85→16.05、complete ρ 0.465→0.585，而 **matched ρ 不变** → SIRV 丰度门只是删转录本，不是存在性判据。见 `CANDIDATE_EXISTENCE_DESIGN.md` |
| krill 0.1.6 overlay（`_img/krill016_{cu122,cpu}`） | 环境 | `init_scale/init_shift` 可用；SIF 自带 0.1.5 不支持；grid runner 加 `PYFIN_PYPATH_PREFIX` |

## 6. 提交状态

- 已提交 `ca82319`：P0 + complete-chain rescue + catalog completion。
- 未提交：junction election + containment evidence（Codex review 3 条待修：election 选民集 P1、containment 回写优先级 P2、每 contest 重建 aligner P2）；junction correction（新，无增益）；三个探针脚本和绘图；`sirv4_grid_runner.py` 的 `PYFIN_PYPATH_PREFIX`；`sirv4_expression_rescore.py` 新增的 jc 模式。
