# 统一候选存在性判据（candidate existence test）— 设计与预注册

状态：**设计冻结中，未实现**。目标是替换掉现在这一串按丰度判存在的门
（`min_abundance` / `max_soft_mass_ratio` / `min_fulllen_fraction` /
`m2_support_gate` + `complete_chain_rescue`），换成一个基于分子证据的判据。

来源：Codex gpt-5.6-sol xhigh 咨询（`.omc/codex_consult_principled_sequin_2026-09-04.md`），
以及 2026-09-04 的双数据集实测。

---

## 0. 为什么现在做这件事（证据）

同一份代码、同一批数据，只换操作点（`--profile`），Sequin 104 格：

| 操作点 | Sn | Pr | matched ρ | complete ρ | complete NRMSE | 召回 /164 |
|---|---|---|---|---|---|---|
| `sirv`（六样本 H9 拟合） | 8.85 | 91.33 | 0.727 | 0.465 | 0.327 | 14.6 |
| `real-drna` | **16.05** | **92.35** | 0.718 | **0.585** | **0.278** | **26.5** |

读法：**matched ρ 几乎不变**（0.727 → 0.718）。也就是说 SIRV 的丰度门
**没有提高"已召回转录本的定量质量"**，它只是删掉了转录本；precision 还略微上升。
它不是一个存在性判据，它是一个针对 SIRV 深度分布的截断。

Sequin 的深度现实：chrIS 上只有 190–575 条 primary read，165 条转录本里每样本只有
28–42 条真正表达，`est_count ≥ 3` 的只有 24 条。所以 `min_abundance=3` 在这套数据上
把 recall 上限直接压到 ~57%。这就是"补丁只适用于一套数据"的量化形式。

结论：不要再为 Sequin 找一套"更合适的常数"。要让常数退出存在性决策。

---

## 1. 判据（一句话）

> 一个候选存在，当且仅当：把它加入模型，能对**其他候选无法解释**的 read 证据，
> 在**留出的 read** 上产生可重复的预测改善，且该改善在样本级 FDR 下显著。

形式化：对候选 `C`，在其 conflict component 内构造两个嵌套模型
`H1(C)`（含 C）与 `H0(C)`（移除 C，其余候选重新拟合丰度），统计量

```
T_C = 2 * Σ_r [ log P(x_r | H1_C^{(-fold(r))}) − log P(x_r | H0_C^{(-fold(r))}) ]
```

`fold(r)` 由 read id 的稳定哈希决定（5 折）：每条 read 只在**没有参与拟合**的那一折里
作为证据计入，因此"用同一批 read 生成候选、拟合丰度、再验证"的选择偏差被切断。

**保留规则是三态的，不是"不显著就删"**（Codex 审查修正，2026-09-04）：

| 状态 | 含义 | 动作 |
|---|---|---|
| `supported` | 检验有效、held-out ΔLL > 0、通过样本级 FDR | 保留，且已获存在性认证 |
| `contradicted` | 通过**独立的、预注册的**负证据/等价性判据 | 唯一允许删除的情形 |
| `unresolved` | 包含 `unidentifiable`、不显著、p 无效、bootstrap 分辨率不足 | **保留**，报表注明未获认证 |

现有检验的方向只能证明"加入 C 有预测收益"。`p > 0.05` 只是**没有证明存在**，不是证明不存在 ——
把它当作删除依据会直接违反红线"证据不足不等于不存在证据"。在有可靠的 contradiction 判据之前，
"大量 unresolved 被保留"是正确且透明的结果，不应伪装成判据成功。

`fold_salt`、bootstrap 次数 `B`、archetype、null 构造在 Sequin 揭盲前冻结。
**`B` 必须与 BY 分辨率绑定**：用 `(b+1)/(B+1)` 的经验 p，family 大小 `m` 下 BY 最强假设需要
`p <= 0.05/(m·H_m)`（m=40 时约 2.9e-4），因此 `B` 至少要能表示该量级（m=40 需 ≳3400 次）。

---

## 2. 三个必须分开的问题

现在的 cascade 把三件事混在一起判：

1. **存在性**：这个结构模型该不该在输出里？→ 本判据。
2. **归属**：这条 read 更像哪个模型？→ per-read likelihood（M1 AS + M2 junction NLL
   + 5' hazard + 3'/polyA），EM 的 E 步。
3. **丰度**：模型该有多少质量？→ EM 的 M 步。

"丰度低所以不存在"是用 (3) 回答 (1)，这正是不可迁移的根源：丰度分布是数据集属性，
而"证据是否支持这个结构"不是。

---

## 3. 证据通道（缺失 = 弃权，不是反对）

对 read `r` 和候选 `c`：

```
ℓ(r,c) = ℓ_align + ℓ_junction + ℓ_5' + ℓ_3'/polyA
```

- **ℓ_align / ℓ_junction**：现有 mappy AS + M2 tie NLL。存在性检验不问"是否至少赢过
  一次"，而是累加所有**真正覆盖 discrimination region** 的 read-level LLR。不覆盖
  差异区的 read 贡献 0，而不是分数支持。
- **indel-clean complete-chain depth**（`discovery.junction_clean_depths`）：不再是
  `K=2/3/10` 的 rescue 类别边界，而是给每个 junction observation 一个**样本内自适应**的
  error probability——indel-clean 且完整穿过的 read 似然比高，带局部 indel 的权重低。
  于是 "2 条 read" 和 "10 条 read" 只是累积 log evidence 的不同数值，不是类别。
- **ℓ_5'**：仅 5' 不同的 parent/short 用条件终止 hazard（`tss_evidence.evaluate_internal_tss`）。
  TSS 上游的 read 对 short 是**强负证据**；TSS bin 聚集是正证据；普通下游降解 read 近似弃权。
- **ℓ_3'/polyA**：polyA PASS 且在候选 TES 附近结束 = 正证据；polyA read 读穿 claimed TES
  = 强负证据；没有 polyA call = 缺失（不是负证据）。注意现在的 `min_polya5p_reads` 把
  polyA 和 **5' 邻近**绑在一起（`quantification.py`），存在性模型里必须解绑：polyA 是 3' 证据。

---

## 4. 四类候选各自必须拿出什么（预注册）

| 原型 | 必须提供 | 不再要求 |
|---|---|---|
| 全长 novel isoform | read 覆盖它与每个主要替代模型的 distinguishing junction/block；这些区域的 M2/sequence likelihood 支持它；junction observation 与样本内 indel-error null 不符 | 固定 ≥3 reads、fulllen fraction、hard-assignment winner |
| 5' 截断 shadow | parent 内部窄 5' peak 高于条件降解 hazard；存在约束 at-risk 群体的上游 read；上游 read 对 short 给负似然 | 额外 3' 证据（与 parent 共享 TES 时） |
| wobble sibling | read 真正跨过 wobble 差异区；M2 NLL 明确偏向它；不是局部高 indel 区产生的位移 junction；在留出 read 上可重复 | 从 anchor 借来的 EM soft mass、posterior tie、canonical motif |
| 被破坏的 annotation | 与 novel **同一套**分子证据：jitter/flip 要有覆盖改动 junction 的 read 级证据；skip/merge 要有支持其独有 junction 组合的 read；trunc5p 过 internal-TSS hazard；trunc3p 要在 claimed TES 有 endpoint/polyA 聚集（read-through 是负证据） | annotation 来源本身（不提供任何存在性证据） |

---

## 5. 统计与过拟合控制

### 5.1 选择偏差：只对丰度做 cross-fitting **不够**

候选的生成、family 聚类、collapse、shadow recovery 都用了**全部** read。即使 H0/H1 的丰度只在
其余四折拟合，留出的 read 仍可能已经参与了"促成该候选被生成 / 触发或阻止 collapse /
决定 discrimination region"。普通的 predictive cross-fitting 校正不了这一层。

因此 **v1 只检验外部预先指定的 GTF multi-exon 假设**：它们不是由样本 read 提议的，
按构造就是 selection-honest（`selection_honest=True`）。这一片同时也能把 likelihood、H0、
bootstrap、FDR 整条管线先验证一遍。

要扩展到 novel candidate，只有两条正确路径（都要另外实现）：
- 一次诚实划分：discovery/collapse/丰度拟合只用 train 折，test 折只评分（最省算力，损失 power）；
- 或每折都在其余四折上**重跑 discovery + collapse + assignment**（≈5 倍 discovery 成本）。

仅按全数据候选的 `supporting_read_ids` 事后过滤 **不正确**：其他留出 read 仍可能通过
union-find / collapse 间接改变了候选集合。

### 5.2 存在性检验用的是**它自己的模型**（Codex 审查修正）

第一版曾写成"生产默认 `beta=0` + `abundance_feedback=False`，所以删除 C 只需逐行重新归一化，
是精确解"。**这是错的**，必须纠正：生产的逐 read softmax **没有任何可拟合的参数**，
因此也没有东西可以 cross-fit；而 assignment 的后验本身并不是"候选是否存在"的生成似然。

存在性检验因此显式定义为**自己的模型**：在 conflict component 上的一个
**共享权重有限混合模型**（mixture weights 每折拟合）。H0 = 删掉 C 后**重新拟合**其余权重
—— 让替代模型重新吸收质量，正是嵌套比较诚实的关键；把其余权重固定会人为削弱 H0、
使检验偏向 H1、p 值反保守。

推论：单条留出 read 上出现**负** Δℓ 是正常的（训练折分给 C 的质量稀释了该 read 更匹配的
替代模型），它表示该折的样本外证据不支持 C，而不是实现有 bug。

**证据矩阵必须是原始的，不能是 EM 消费的那个。** 生产在 `m2_diff_cover_gate=True`（默认）
下把 `d_tx` 写成**决策矩阵**：达到 margin 的判别 read 变成 one-hot（只给最优候选 0，
其余 MISSING），而 fuzzy read 携带的是整个 component 的 `covered_vote` 比例，不再是该 read
自己的 M2 似然。用它评分会同时（a）跳过最强的证据 read，（b）把别的 read 的投票当成这条
read 的似然喂回来。因此 `QuantOutput` 另外保留 `structural_d_tx_raw`（gate 变换**之前**、
raw-AS tie 集上的 per-read junction NLL 对比），存在性模块只读这个矩阵。

### 5.2b 折的划分必须按 component 自适应

全局 `hash(read) mod 5` 不能保证一个小 component 的每折非空 —— 在 Sequin 的深度下
（每条染色体 190–575 read）大多数折没有测试 read，候选因纯技术原因被丢弃。
改为 **按 component 自适应**：`K = min(5, n_component_reads)`，用
`hash(component_id, read_id, salt)` 排序后 round-robin 发牌。于是每折非空、每条 read
恰好被评分一次且不在自己的训练集中；`n ≤ 5` 时自然退化为 LOO（正确的极限），
`n < 2` 没有诚实划分 → `unresolved`。**同一 component 的所有候选共用一份折划分**：
给每个候选单独优化划分会是第二重选择效应。

注意 LOO 的各折拟合高度相关：这不破坏 honesty，但意味着不能把各 read 的贡献当独立量套
正态/χ² null，bootstrap 必须整套重跑所有折拟合，而不是对已算好的 `fold_delta_ll` 重采样。

### 5.3 scope 与 read 集合

- refit **只在 conflict component 内**，前提是 component 按"共享有限 likelihood support"构建：
  若存在 read 对 component 内外候选都有有限 likelihood，两个集合不独立，必须合并。
- 使用**所有对该 component 至少一个候选有有限原始 likelihood 的 interval read**；
  不能只用 `R>0`、hard-assigned 或对 C 有非零 responsibility 的 read —— 那些集合已由 H1 结果决定
  （post-selection conditioning），而且会漏掉"删掉 C 后本应改判给替代模型"的 read。
- **不要**用 discovery 的 `family_id` 当统计 component：生产 m2_em 是 interval-wide assignment。

### 5.4 identifiability（必须保留，不得删除）

对候选 C 与折 k，held-out discriminating set：

```
D(C,k) = { read i : fold(i)=k
           且 C 可由 training reads(−k) 独立生成，或是外部预指定假设
           且 C 对 i 有有限 likelihood
           且 component 内至少一个 alternative 对 i 有有限 likelihood
           且 i 真正覆盖 C 与该 alternative 的 discrimination region }
```

满足任一条即 `unidentifiable`：所有折都无法在不用 held-out read 的情况下生成 C；所有 `D(C,k)` 为空；
在全部 eligible held-out read 上 H1 与 H0 的预测分布数值相同；C 的 likelihood 列与某个
alternative（或其混合）在可观测 read 上完全等价；active channel 全缺失；bootstrap null 退化。
"bootstrap 次数不足以达到 BY 所需分辨率"单独标记为 `untestable_resolution`，同样按 unresolved 保留。

### 5.5 null 与 FDR

- **null**：不用渐近 χ²（mixture weight 在边界 π=0、候选共享 read、选择过程复杂）。
  用局部 parametric bootstrap，在 `H0(C)` 下固定深度/长度/信号缺失模式/样本级 error 率。
- **FDR 是分层的**：component 内先做 omnibus 检验、再用 Bonferroni/Holm；跨 component 用 BY。
  必须记录候选数、identifiable 数、component 数与最小可达 p。
- **必须排除在 FDR family 之外**：`identifiable=False`、`selection_honest=False`、`p_valid=False`、
  bootstrap 退化、channel 未激活、p 分辨率不足。`fdr_family_size` 记录**实际进入校正的检验数**，
  不是 TSV 行数。
- FDR 是**样本级**：`_quant_m2_em` 是逐 interval 的，所以 interval 只写 shard，
  合并与 BY 放在 runner 的 finalize（与其他 audit artifact 同层）。
  **不按染色体做** —— chrIS 只有一条染色体，按染色体做会在迁移到真实多染色体数据时突然改变负担。

## 6. 落地方式（分阶段，audit 优先）

新模块 `fin/analysis/candidate_existence.py`，只暴露：
`build_conflict_components` / `evaluate_candidate_existence` / `apply_sample_fdr`。

**conflict component v1 的形式定义**：候选 `a,b` 之间有边，当存在 read `i` 使得
`a,b ∈ ties_by_read[i]` **且** `a,b ∈ nlls_by_read[i]`。
**coverage 不是建边条件**（Codex 审查修正）：component 界定的是统计**依赖**——哪些候选
在争夺质量、H0 refit 的范围到哪里——而一条 read 支持的两个候选无论能否被它区分都是耦合的。
coverage 决定的是这条 read 能否**贡献证据**，因此作用在评分那一步（逐 read）。
连通分量按 `candidate_id` 排序，`component_id` 由排序后的 id 稳定哈希生成。
**不能用 `R[i,j] > ε` 建边** —— 那会让图依赖已拟合的 H1 后验，且 ε 会改变检验 family。

**v1 切片（本周）**：archetype 只做 `gtf_multi_exon_m2` —— 外部 GTF 多外显子候选，
其 raw-AS tie 内至少有一个多外显子 alternative，且 M2 对至少两个候选有有限 NLL。
证据通道只做：M1 raw-AS tie 作 applicability mask、M2 有限 NLL contrast、
`cover_by_read is True` 作 eligibility；junction clean depth 只作 audit 协变量，不进统计量。
**不要**同时上 5' hazard 或 polyA，否则校准失败时无法归因。
（polyA 通道在该调用点根本拿不到：`CandidateSet`/`QuantOutput` 都没有 per-read polyA call，
需要另外新增 evidence collection。）

调用点：interval 的最终 assignment 完成后、ledger/selection 之前（`runner._quant_m2_em`）
只写 shard；sample 级合并 + BY 放在 runner finalize。

产物 `candidate_existence.tsv`（最低字段）：
`sample_id, interval, chrom, strand, candidate_id, source, family_id, archetype,
component_id, component_size, dominant_alternative, selection_honest,
identifiable, identifiability_reason, n_candidate_folds, n_tested_folds,
n_discriminating_reads, n_covered_reads, n_finite_contrasts,
delta_predictive_ll, fold_delta_ll, structural_llr, tss_llr, tes_polya_llr,
null_method, bootstrap_replicates, bootstrap_exceedances, p_value, p_valid,
q_value, fdr_family_size, fdr_method, status, decision,
h1_converged, h0_converged, software_version, fold_salt, bootstrap_seed`
（v1 的 `tss_llr` / `tes_polya_llr` **留空，不填 0**。）

**接线任何 drop 之前，audit 必须先证明**：
1. assembly/quant 输出逐字节不变，audit 只写 artifact；
2. read 顺序、候选顺序、worker shard 顺序都不改变 component / T / p / q；
3. 全数据候选没有被错误标成 `selection_honest=True`；
4. null/control 分层的 p 值不过度集中于 0（报 QQ/直方图与 bootstrap 分辨率）；
5. 同一分子证据下 GTF/novel 来源不改变统计量（source 只能是报表字段）；
6. `p_valid=False`、`unidentifiable`、分辨率不足的行完全不进 FDR family；
7. component 行先生成，样本级 FDR 只在 finalize 执行；
8. SIRV 开发集上 TP/FP 的 T 分布有方向性，但**不得用该结果回头调 `q`**；
9. Sequin 揭盲前冻结 fold salt、B、archetype、null。

判据稳定后才谈删除或降级为 compatibility mode 的门：m2 support gate、complete-chain rescue、
abundance floor、isoform fraction、soft/hard ratio、fulllen fraction。
`canonical_gate` 可保留为低成本序列合法性检查，但不能代替存在性证据。

## 7. 跨数据集协议（对 HANDOFF.md §5 的扩展）

- **SIRV4 = 开发集**（6 样本 × 26 条件 = 156 格，多异构体/窄动态范围，结构判别的压力测试）。
- **Sequin = 留出泛化集**（4 样本 × 26 条件 = 104 格，2 异构体/基因、6.9 个数量级浓度、
  低深度；**独立浓度真值**）。**不在 Sequin 上调任何参数**；它只用来揭盲。
- 两套数据**分别**做 paired bootstrap，**不合并**（设计、深度、复杂度都不同，合并的
  平均数没有意义）。
- 定量轴必须同时报 `matched_*`（条件于召回的定量质量）和 `complete_*`（未召回记 0 的
  联合指标）；只报前者可以靠"只输出最容易的几条"刷高。
- "改动可迁移"的量化定义：在两套数据上**方向一致**，且在留出集（Sequin）上
  Sn 与 Pr 均不显著变差（paired bootstrap 95% CI 下界 > −0.5），
  同时 complete Spearman 不下降。


---

## 8. 实测结论（2026-09-05，E1–E4 完成后的审计）

**null 的校准**（`existence_null_calibration.py`，三种 H0 构造，B=255–511）：
within-read label-permutation bootstrap 的 `p+`、`p−` 在 0.025/0.05/0.10 三个水平上均 ≤ 名义值
（例：dense k3 `p+≤.05` 0.017、`p−≤.05` 0.042），H1 power 在 n=6/effect 0.5 时 0.80、
n=30/effect 0.2 时 0.87，重复列（duplicate）退化为 unresolved。曾经试过的两种 null 被数据否决：
centred pairs bootstrap（`p−` 在 0.05 水平上 0.15，因为冗余候选的 cross-fitted T 期望为负），
duplicate-column null（零方差）。

**真实数据上的行为**（SIRV4 8 条件 × 6 样本，B=1023）：识别出的注释候选 465 条中 83% 有 T>0，
但 **`p+ ≤ 0.025` 为 0%**。原因是结构性的：几乎每个注释候选在 component 里都有一个 wobble 孪生
（junction 相差 ≤6 bp 的 read 派生变体），与近乎相同的伙伴交换标签不改变统计量，null 复现观测值。
对"这个候选 id 是否存在"而言这是**正确**的判决——它无法与孪生区分。

**推论**：候选级存在性认证不能成为生产决策（否则"unresolved 一律保留"会退化成保留全部）。
Part A 实际蕴含的分解是：**存在性在 component（locus）级别由 omnibus 认证；哪个结构变体承载质量
是 identifiability / 后验问题，由 assignment 的混合权重决定。** 因此本判据保持 audit-only；
面向生产的工作转到把各证据通道放进 assignment 后验（Part A8）——`--assignment-channels tes`
已实现（3' polyA 项进入 d_tx；单外显子位点保留 read 自身的 extent 作为独立假设），
5' 项同理。这也解释了三套数据上一致的短板（trunc3p/trunc5p/merge 的 precision）：
它们是端点问题，决定权在 assignment 而不在门。

## 9. 端点通道落地：注释模型的 read-extent 假设（2026-09-05 下午）

**为什么 trunc 条件的 precision 会掉**（HEYA8 `c_trunc3p50`，SIRV107，可复现：
`experiments/prod_validation/_heya8grid/_dbg_sirv107/`）：注释被截短 100 bp；reads 的主导
junction 链是注释（非经典位点）的 2–3 bp minimap2 位移，被 canonical gate 删掉；mappy 比对随即
偏向**带完整 span 的 wobble novel**（每条 read 都比截短注释多出 100 bp 可比对序列），注释饿死
（0 read），三个 wobble `j` 假阳性拿走 678 reads。也就是说，fold 把 read 支持的**extent**
抹掉了，只保留了注释的 extent——这是一个 extent 参数被错误固定的问题，不是 junction 问题。

**设计**（Part A：extent 是模型参数，注释只是先验值）：当 sample nuisance 模型判定 reads
与注释端点矛盾（polyA-PASS read 越过所声称的 TES 且 narrow_posterior<0.5，或 read 起点在
所声称 TSS 上游且 start.narrow_posterior<0.5；**排除**能被同链/wobble 链上其它注释模型的端点
解释的 read——注释 APA / 可选 TSS 不是矛盾），discovery 额外发出一个 **read-extent 兄弟**：
同一 intron 链，reads 的 span，`candidate_id = <gtf_id>.ext`，`source="gtf"`（继承 junction
信任：canonical gate、gtf 丰度下限、mono gate 豁免），`extent_of = gtf_id`。兄弟与注释是
**同一模型的两个 extent 估计**，不是两条转录本：selection 的第一步 `resolve_extent_pairs`
保留后验质量更大的那个（平局归注释），把另一个的 reads / 软质量 / tie 成员 / NLL 并入赢家
（后续所有门看到的是合并后的一列）。多外显子键的矛盾证据集 = family 中与该链 wobble 等价的
全部 read（它们就是该注释链声称要解释的 read，无论其变体是否活过后面的门）。

**为什么不是 novel 兄弟**：第一版把兄弟作为 novel 发出，Sequin 上 1–3 read 的截短注释（gffcompare
仍记 `=`）的 reads 被兄弟拿走后，兄弟死于 novel 专属下限（`novel > 1`、mono `min_reads=5`），
注释饿死，位点整体消失：trunc 条件 Sn −2.7…−5.8。`source="gtf"` + pair resolution 使幸存者
面对的门与注释本来会面对的完全一致。

**Codex 审查（gpt-5.6-sol xhigh）后修正**：同链 APA 的候选特异性矛盾判据；合并后的 R/ties/nlls
对后续门可见；catalog completion 把 `.ext` 视为覆盖其父结构；existence audit 的
selection_honest 与 containment n_bins 排除 `.ext`（其端点是从数据里选的）；重复匹配的兄弟
合并 span/reads 而非丢弃；outcome ledger 每候选一条决定性结果；work dir 的 sample evidence
校验版本与输入身份。默认路径（`assignment_channels=none`）逐字节不变（SIRV1 复现：54/54 行相同）。

**TSS require**（`folded_tss_recovery_mode=require`，三套网格）：SIRV +0.05（13↑/0↓），
Sequin ±0，HEYA8 +0.03（CI [+0.01,+0.07]，6↑/1↓）——三套非负，符合 maximin 默认准入。
`ends` 通道（含 read-extent 假设）的三网格结果见 endsarm4（冻结快照 `_snap/endsarm4` 运行）。
