# 生产流水线逐步演算：read 怎样变成最终转录本

## 0. 这份文档描述什么

本文解释 `dev` 分支 HEAD `19d5399` 的 `--profile real-drna`、不额外覆盖参数时，数据怎样经过生产流水线；图是路线图，本文是逐步演算，规则以源码为准。
`docs/diagrams/README.md` 开头标的是较早的 `aa21414`，因此不能把它的版本标记当成本次 HEAD；下面重新核对了实际调用与默认配置。
生产组合是 families discovery、`m2_em`、junction 与 ends 证据、folded-TSS recovery、区间选择、全样本选择、snap、merge、refit；不是把所有实验开关同时打开。依据：`fin/pipeline/config.py:38`、`fin/pipeline/config.py:174`、`fin/pipeline/assignment.py:1086`。

**先纠正一个会让例子走不通的前提：注释模型不能与 novel 模型在 snap 阶段合并。**
注释模型的 `.ext` sibling 也仍然是 `source="gtf"`；snap 和它的两种结构合并只处理 novel 多外显子。
所以本例保留一个注释模型 A，再放入两个 novel wobble 变体 B、C：A 展示注释端点修正，B、C 展示 snap 后合并，不能把两种身份变化混为一谈。依据：`fin/candidates/discovery.py:520`、`fin/pipeline/junction_snap.py:217`、`fin/pipeline/junction_snap.py:244`。

### 贯穿全文的 locus 与记账口径

所有基因组坐标除明确写作 GTF 输出者外，都是 **0-based、左闭右开**；例如 `[1043000,1043457)` 长度是 457 bp。
这与 CIGAR 提取和内部 exon 表示一致；GTF 写出时 start 加 1，end 不变。依据：`fin/candidates/intron_chains.py:18`、`fin/io/io_gtf.py:583`。
以下 read、AS、计数和 posterior 是教学构造，不是实测 benchmark；明确标“实测”的核参数才来自仓库记录。
所有模型位于 `chr1` 正链，第二个 intron 都是 `(1043900,1044500)`。

| 教学简称 | source | 第一个 intron | TSS | TES | 内部 exon 列表 |
| --- | --- | --- | ---: | ---: | --- |
| A | gtf | (1043457,1043600) | 1043000 | 1045000 | (1043000,1043457), (1043600,1043900), (1044500,1045000) |
| B | novel | (1043454,1043600) | 1043010 | 1045183 | (1043010,1043454), (1043600,1043900), (1044500,1045183) |
| C | novel | (1043451,1043600) | 1043030 | 1045183 | (1043030,1043451), (1043600,1043900), (1044500,1045183) |

A、B、C 是解释用的短名字，不声称 novel 输出 ID 就是单个字母；实际 novel ID 由 discovery 生成。依据：`fin/candidates/discovery.py:565`。
先跟踪六条 read；为了演示生产支持阈值，再给出同一区间其余八条 read 的分组，总共 14 条。
这不是让算法先只看六条再补数据：真实 discovery 一次看到全部 14 条，我们只是先放大六条观察。

| read | 原始链 | reference_start | reference_end | polyA length | QC |
| --- | --- | ---: | ---: | ---: | --- |
| r1 | A | 1043000 | 1045000 | 46 | PASS |
| r2 | A | 1043020 | 1045180 | 42 | PASS |
| r3 | B | 1043010 | 1045180 | 51 | PASS |
| r4 | B | 1043025 | 1045183 | 37 | PASS |
| r5 | C | 1043030 | 1045183 | 44 | PASS |
| r6 | 仅 (1043900,1044500) | 1043700 | 1045181 | 39 | PASS |

其余八条是 r7–r9 的 3 条 A 链、r10–r13 的 4 条 B 链、r14 的 1 条 C 链。
为固定几何，补充 read 不超出表中各组端点；其中若干 A 链 read 的 3′ 端也达到 1045180–1045183。
完整的原始链计数因此是 A=5、B=6、C=2、短链=1，总数 `5+6+2+1=14`。
**原始链支持数不等于最终 abundance**：前者来自 BAM 的 splice，后者来自重新比对后的 responsibility 列和。依据：`fin/pipeline/junction_snap.py:19`、`fin/analysis/quantification.py:94`。

## 1. 输入与样本级预处理

### 1.1 先决定哪些 read 在同一个工作区间

入口需要 BAM、genome、FASTQ、signal、输出目录；GTF 是指导注释，不是表达真值。CLI 的缺项检查包括 `--signal`，不能把 runner 内部的降级能力当成 CLI 默认允许无信号。依据：`fin/cli.py:420`。
`run()` 先调用 `generate_isolated_intervals()`，把 read 覆盖产生的区间与注释 span 一起考虑，然后按区间处理。依据：`fin/pipeline/runner.py:193`。
区间合并按 chromosome 和 strand 分桶；桶内按 start 排序，重叠或间隔不超过 `max_gap` 就合并，默认 `max_gap=0`。依据：`fin/io/interval_manager.py:403`、`fin/pipeline/config.py:165`。

例如先有正链 read 覆盖岛 `[1043000,1043500)`、`[1044500,1045183)`。
两个岛相距 `1044500−1043500=1000 bp`，仅凭 `max_gap=0` 不会相接。
加入 A 的注释 span `[1043000,1045000)` 后，它同时覆盖两岛的一部分，合并结果成为 `[1043000,1045183)`。
同坐标的负链区间仍是另一个工作单元；“isolated”不意味着正负链在基因组坐标上绝不重叠。依据：`fin/io/interval_manager.py:403`。

### 1.2 polyA 判定是 signal 证据，端点位置来自 BAM

生产的 `assignment_channels="ends"` 触发样本级预遍历；默认 `threads=1` 时后续区间串行处理。依据：`fin/pipeline/config.py:70`、`fin/pipeline/runner.py:202`。
预遍历按 read ID 去重，取 primary alignment 的链向、3′ 和 5′ 基因组端点；正链 3′=reference_end，5′=reference_start，负链互换。依据：`fin/pipeline/runner.py:1616`、`fin/pipeline/runner.py:1715`。
krill 用 whole-read basecalled sequence 的 `self` variant 做 `polya=True` eventalign，返回 tail length 和 QC；PyFIN 接受 caller 的判定，而不是自己数 FASTQ 末尾有几个 A。依据：`fin/scoring/polya.py:18`、`fin/scoring/polya.py:59`。
默认端点通道只接受 `PASS` 且可用长度 `>0`；缺失、失败和默认不接纳的 QC 是弃权，不是给某个转录本一张反对票。依据：`fin/analysis/molecule_evidence.py:51`。

对 r2，krill 返回 `(42,"PASS")`，因此允许它提供 3′ 位置证据。
它的 3′ 基因组端点是 `1045180`，**不是** `1045180+42`；tail length 不是 genomic extension 长度。
相对 A 的 TES，signed distance 为 `1045180−1045000=+180 bp`，正号表示 read 越过候选宣称的终点。
如果同一 read 返回 `(0,"PASS")` 或没有调用结果，这条 3′ 证据就不使用；普通 5′ 通道仍有独立的可用性检查。依据：`fin/analysis/molecule_evidence.py:51`、`fin/pipeline/assignment.py:617`。

### 1.3 缓存的是共同证据，不是某个候选的最后判决

缓存路径由输入身份确定；有匹配身份、可复用状态的 artifact 时跳过 krill，发布采用 create-only 语义，避免 worker 读到不同版本。依据：`fin/pipeline/runner.py:1616`、`fin/pipeline/runner.py:1818`。
例如第一次读入 14 条本 locus 的 primary 记录，得到 14 组端点；第二次相同输入可重用这些值，不再对这 14 条重复运行 whole-read polyA。
这不等于整次流水线只调用一次 krill：后面 junction tie 的 signal 评分是另一个工作。依据：`fin/scoring/polya.py:59`、`fin/pipeline/assignment.py:240`。
缓存不伪装失败：有 signal 但 polyA 未成功时，不应把结果当成“样本没有 polyA”的可复用答案；不可用记录让对应通道弃权。依据：`fin/pipeline/runner.py:1629`。

### 1.4 nuisance：从同一样本拟合“端点本来就会抖多少”

3′ 模型是两个 **离散 Laplace** 的混合，不是单个高斯，也不是固定的“距离超过 25 bp 就删除”。依据：`fin/analysis/sample_nuisance.py:97`。
令 `q=exp(−1/b)`，整数距离的概率质量是 `L(d;b)=(1−q)/(1+q) × q^|d|`。
混合为 `P(d)=(1−eta)L(d;b_narrow)+eta L(d;b_broad)`，`loglik=ln(P)`，`NLL=−ln(P)`；窄项表示 jitter，宽项给大偏差留下非零概率。依据：`fin/analysis/sample_nuisance.py:97`。
拟合用样本内端点相对输入注释 TES 的距离；排除当前区间后还有至少两个这种距离，就走注释路径，否则走 read 端点模式的 split-half 距离路径。依据：`fin/analysis/sample_nuisance.py:477`。
split-half 在一半 read 上建 mode，另一半量距离，再交换角色；避免一个 read 用自己帮助定义的 mode 证明自己端点很准。依据：`fin/analysis/sample_nuisance.py:159`。
拟合结果至少需要 `n_fit≥50` 且 `b_narrow≤25 bp` 才是可用核；“两个距离能进入拟合”不等于“两条 read 足够启动通道”。依据：`fin/analysis/sample_nuisance.py:74`、`fin/analysis/sample_nuisance.py:97`。

仓库人类样本的**实测**注释拟合为 `b_narrow=1.8686806129384694 bp`、`eta=0.20867754234594615`，约为 1.9 bp、0.21。依据：`HANDOFF.md:841`。
下面演算用四舍五入的 `b_narrow=1.9`、`eta=0.21`，另设教学用 `b_broad=100 bp`、`n_fit=10000`；**100 bp 和 10000 不是上述实测记录**。
`q_narrow=exp(−1/1.9)≈0.590778`；`q_broad=exp(−1/100)≈0.990050`。
把 r2 的 180 bp 代入得到 `P(180)=0.0001735624`，`loglik(180)=−8.658973`。
3 bp 则是 `P(3)=0.0429223347`，`loglik(3)=−3.148363`；因此同一 read 更容易由距离 3 bp 的 TES 解释，差 `8.658973−3.148363=5.510610 NLL`。
这仍是概率比较，不是宣称 r2 生物学上绝不可能来自 A。

5′ 另拟合 one-sided start mixture；只惩罚 read 出现在候选 TSS 上游的情况，候选 TSS 下游的 read 在这个普通通道弃权，不因降解而自动受罚。依据：`fin/analysis/sample_nuisance.py:235`、`fin/pipeline/assignment.py:684`。
例如教学 start 核也暂用 `(1.9,100,0.21)`，read 在 TSS 上游 20 bp 时附加 `NLL(20)−NLL(0)=7.052653−1.588288=5.464365`。
read 在 TSS 下游 700 bp 时，该通道增加 0，不增加“700 bp 的终止距离代价”。

**assignment 用 leave-one-interval-out，discovery 的 contradiction hook 用 pooled。**
假设样本有 10014 条可用端点观测，本区间有 14 条，assignment 排除这 14 条后拟合；discovery 可用全部观测。
代码实际排除的是 `candidate_set.read_ids`，不是根据候选最后有多少 hard reads 来排除；discovery 则没有 interval exclusion。依据：`fin/pipeline/runner.py:1974`、`fin/pipeline/runner.py:1884`。

## 2. Discovery：单区间内怎样生成候选

### 2.1 CIGAR 变成 intron chain

提取器沿 reference 坐标走，`N` 产生 intron，`M/D/=/X` 消耗 reference，`I/S/H/P` 不推进 reference。依据：`fin/candidates/intron_chains.py:18`。
r1 从 1043000 开始，教学 CIGAR 为 `457M143N300M600N500M`。
`457M` 到 1043457；`143N` 记录 `(1043457,1043600)`；`300M` 到 1043900；`600N` 记录 `(1043900,1044500)`；最后 `500M` 到 1045000。
输出链是 `((1043457,1043600),(1043900,1044500))`；外侧 start/end 另存，不塞进 chain key。
r2 对应 `437M143N300M600N680M`，虽然左右端都与 r1 不同，得到的两个 intron 完全一样。
r6 则是 `200M600N681M`，只得到第二个 intron；“少一个 intron”不等于已经证明存在一个短转录本。

### 2.2 精确分组：端点先不参与身份判断

family discovery 先按精确 intron tuple 汇集 read，端点不用于拆开相同链；原始 spans、CIGAR footprint、junction 支持另存。依据：`fin/candidates/discovery.py:337`。
六条放大 read 的精确组是 `{r1,r2}`、`{r3,r4}`、`{r5}`、`{r6}`，即 2、2、1、1 条。
完整 14 条数据对应 5、6、2、1 条；r1 的 TES=1045000 与 r2 的 TES=1045180 不会因此分成两个 chain key。
这一层只回答“观察到哪些精确剪接组合”，不回答“每个组合最后有几种 TSS/TES”。

### 2.3 family：把需要比较的近邻放在一起，不把坐标改成共识

wobble 要求 intron 数相等，而且逐个 donor、acceptor 都在 6 bp 内；family 是关系图的连通分组，不是投票选出一个统一链。依据：`fin/candidates/chain_cluster.py:223`、`fin/candidates/chain_cluster.py:473`、`fin/pipeline/config.py:176`。
A 与 B 第一 donor 差 3 bp，B 与 C 差 3 bp，A 与 C 差 6 bp；acceptor 差均为 0，第二 intron 也相同，因此进入同一 family F。
r6 的短链是长链的连续子链，containment 关系也可把它接入 F；family containment 可以容忍逐 junction wobble，但下一步 collapse 的 exact-subchain 要求更严格。依据：`fin/candidates/chain_cluster.py:231`、`fin/candidates/chain_cluster.py:672`。
例如两内含子 `(1043457,1043600),(1043660,1044500)` 与一个跨过小 exon 的 `(1043457,1044500)` 可形成 cassette 关系。
中间 exon 长 `1043660−1043600=60 bp`，满足 `0<60<70`；恰好 70 bp 不满足，代码是严格小于，不是“≤70”。依据：`fin/candidates/chain_cluster.py:243`。
这些只是候选的关系；此时 donor 1043451、1043454、1043457 仍然是三个不同坐标。

### 2.4 collapse：精确子链折叠成 shadow

family 内按链长优先、支持数及坐标排序；短链尝试进入最长的 eligible retained container，注释链永远不能当 foldee。依据：`fin/candidates/chain_cluster.py:729`。
完整数据中 B 有 6 条、A 有 5 条、C 有 2 条，长链长度相同，所以本例短链 r6 先找到 B 作为容器。
r6 的唯一 intron 与 B 的第二 intron 完全相等；它 5′ 从 1043700 才开始，没有走过 B 的第一个 intron，于是折叠到 B。
B 的生成支持变成 `6+1=7`，但其“原始完整 B 链 read 数”仍是 6；r6 以 dormant shadow 留在 B 的 shadow 记录中。依据：`fin/candidates/chain_cluster.py:747`。
A、B、C 都有两个 intron，彼此不是严格更短的 exact subchain，因此不会在这一步把 wobble 变体折成一个。
若 r6 的短链本身有 GTF 注释，它即使看起来像降解也不能在 discovery 被折掉；后续 assignment 和 gate 再判断。依据：`fin/candidates/chain_cluster.py:755`。
纯 annotation-only 链也不能无条件吸走有 read 的子链；容器资格还检查它是否属于未被 read 矛盾的注释集合。依据：`fin/candidates/chain_cluster.py:729`。

span guard 检查的是 **CIGAR 是否外显子式地穿过额外 intron**，不是仅看 read 的总 span 是否覆盖那里。依据：`fin/candidates/chain_cluster.py:163`。
对 B 第一 intron `(1043454,1043600)`，窗口为 `[1043444,1043610)`，长度 `146+20=166 bp`。
另一个教学 read 如果在这个窗口里是连续 166M，则左右 flank 各匹配 10 bp、body 匹配 146 bp、无 N、无 D，说明它读穿了 B 所谓的 intron，应保留为不同结构。
若是 `10M146N10M`，则有 N，并不是读穿；若中间有单个 10D，也超过单个 deletion≤3 的上限，不能作为 guard 的保留证据。
完整条件还包括总 deletion≤6、每侧至少 8/10 matched bp、body 至少 `146−6=140` matched bp。依据：`fin/candidates/chain_cluster.py:163`。

### 2.5 emit：生成序列与端点假设

exact chain 命中 GTF 时发出该链的注释模型，保留各自注释几何；没有 exact GTF 对应的 B、C 则用其 read span 的 min(start)、max(end) 生成 novel，再从 genome 拼接 cDNA。依据：`fin/candidates/discovery.py:520`、`fin/candidates/discovery.py:583`。
本例 A 先保留 `[1043000,1045000)`；B 是 `[1043010,1045183)`；C 是 `[1043030,1045183)`。
r6 已折到 B，但其 start=1043700 比 B 的 1043010 晚、end=1045181 比 1045183 早，不改变 B 的外侧边界。

对 A，pooled endpoint 核发现 PASS read r2 越过 TES 180 bp，而且窄分量 posterior 约 `8.4×10^−39<0.5`。
如果同链或 wobble-equivalent 的另一个注释模型能用其窄核解释这个端点，就不把该 read 当成 A 的矛盾；本例没有这样的另一注释。依据：`fin/candidates/discovery.py:618`、`fin/pipeline/runner.py:1934`。
因此只重估 A 的 TES：矛盾 read 最大 end 是 1045183，`max(1045000,1045183)=1045183`；TSS 没被矛盾，仍为 1043000。
发出 A.ext，exons 是 `(1043000,1043457),(1043600,1043900),(1044500,1045183)`，`extent_of=A`、`source=gtf`。依据：`fin/candidates/discovery.py:299`、`fin/candidates/discovery.py:546`。
这一步不是把 A 原地改长，也不是看到任意更长 read 就取全体 read 的最外端；它保留两种 extent 假设，只用矛盾 read、只动被矛盾的端、只向外。
相反，若 PASS read 停在 1044990，距离为 −10 bp，不满足 TES contradiction 的 `d>0`，不能让 A 向内缩 10 bp。依据：`fin/pipeline/runner.py:1951`。
发现上游 TSS contradiction 还需样本中邻近 start 支持达到默认 3 条；这项 cluster 要求属于 discovery，不属于普通逐 read 的 5′ assignment penalty。依据：`fin/pipeline/runner.py:1925`、`fin/pipeline/runner.py:1956`。

单外显子 read 进入独立 holding bucket，按空间重叠拆组，再匹配单外显子注释或生成 novel；生产中不在生成阶段直接折进 multi，而留给 post-EM mono resolution。依据：`fin/candidates/discovery.py:658`、`fin/pipeline/config.py:180`。
例如两条 mono read `[1044600,1044800)`、`[1044700,1044900)` 可组成 span `[1044600,1044900)`；另一个远处 `[1050000,1050200)` 不应仅因“空链相同”与它们合成一个模型。
本例主线的 14 条都有至少一个 intron，因此这里没有 mono 候选；上面的 mono 数字是独立分支演示。

## 3. 区间内 pre-assignment gate

### 3.1 canonical 门会删候选，但不会删掉它的 read 竞争资格

生产 canonical gate 检查 novel 的每个 intron；GTF、fusion 豁免，mono 无 junction 自然通过，无 genome sequence 时弃权。依据：`fin/pipeline/selection.py:1016`。
默认可接受 motif 为 `GT-AG`、`GC-AG`、`AT-AC`；正链 donor 查 `genome[d:d+2]`，acceptor 查 `genome[a−2:a]`，负链按 transcript orientation 取反向互补。依据：`fin/pipeline/config.py:485`、`fin/candidates/canonical.py:8`。
本例假设 B、C 的第一 donor 位置都读出 GT，acceptor 1043600 前两碱基为 AG，第二 intron 也为 GT-AG，因此两个 novel 都过门。
如果 C 的第二个 intron 是 `GT-AA`，即便第一个 intron 合格，也因“不是每个 junction 都 canonical”被删除。
A 与 A.ext 都是 GTF source，不经过 novel motif 否决；这不是声称 GTF 一定正确，而是该门的明确豁免规则。
候选删除不改 `candidate_set.read_ids`：例如 C 被删，r5、r14 仍可重新比对竞争 A/A.ext/B。依据：`fin/pipeline/selection.py:1016`。

### 3.2 junction registry 是记录，不是另一层投票淘汰

discovery 收集 raw/indel-clean junction depth，runner 写出 registry shard；没有 gate 回读这个 TSV 决定去留。依据：`fin/candidates/discovery.py:363`、`fin/pipeline/runner.py:723`。
若本例所有边界 flanks 都 indel-clean，第一 junction 三行 raw/clean 分别是 A `5/5`、B `6/6`、C `2/2`，第二 junction 则是 `14/14`。
第二 junction 的 14 条包括只带这一条 junction 的 r6；它对该 junction 有证据，但不因此变成完整 B 链证据。
registry“只记录”不意味着 junction 不参与别处：assignment 用 signal，区间支持门看观测 junction，final snap 另读全样本支持。
`junction_dominance_filter`、junction election、read-level junction correction 默认关闭，不插入本例流程。依据：`fin/pipeline/config.py:204`、`fin/pipeline/config.py:208`、`fin/cli.py:243`。

## 4. Assignment：先限定竞争集合，再用额外证据分配 read

### 4.1 重新对候选 cDNA 比对，生成 read×candidate AS 矩阵

`Assigner.assign()` 为各候选 spliced sequence 建 mappy aligner，将有序列的 read 对各候选比对，保留各自最佳合法 hit；没有 hit 的 cell 是缺失，不是零分命中。依据：`fin/pipeline/assignment.py:825`。
评分重建为 `2×matches−4×mismatches−Σgap_cost`，gap 长 k 的代价是 `min(4+2k,24+k)`；单个 I/D 超过默认 50 bp 的 hit 不接受。依据：`fin/scoring/mappy_score.py:15`。
例如 610 matched bases、5 mismatches、无 gap，AS=`1220−20=1200`；再有一个 3 bp gap 且其他统计不变，扣 `min(10,27)=10`，得到 1190。
教学 AS 可以是局部比对结果；不能仅凭 genomic span 长度直接推算匹配碱基数。

对 r2 假设 mappy 合法 hit 给出以下分数；这是一行具体教学输入，不声称由前面的简化 CIGAR 唯一推出。

| read | A AS | A.ext AS | B AS | C AS | best-AS tie set |
| --- | ---: | ---: | ---: | ---: | --- |
| r2 | 1200 | 1200 | 1194 | 1188 | A, A.ext |

best-AS tie 是每个 read 自己的集合，容差为 `1e−9`；B 差 6 分不在 tie，后续端点证据不能把它拉回。依据：`fin/pipeline/assignment.py:54`。
即使 B 的 TES 看起来很接近 r2，也不能跨过 AS 这道候选资格限制。

### 4.2 Junction：summed-LLR 产生 mask 或一次性 covered-vote prior

生产用 `summed_llr`、差异窗口 flank=8、区分 margin=1.0；intron 数不等的 tie 或无差异窗口的 tie，junction contrast 弃权。依据：`fin/pipeline/config.py:58`、`fin/pipeline/assignment.py:71`。
r2 的 A/A.ext intron chain 完全一样，没有 junction 差异窗口，所以它这一行不靠 junction 选端点。
同链 A/A.ext 在 junction 层属于同一种结构，不能互当 runner-up 把真正的结构 margin 冲成 0。依据：`fin/pipeline/assignment.py:975`。

另看 r3 的一个教学 tie `{B,C}`：krill 差异窗口 summed NLL 为 B=120.0、C=122.4。
差 `122.4−120.0=2.4≥1.0`，B 所在结构获准，C cell 留在 MISSING mask；**不是**把 2.4 作为连续代价后留两者做 softmax。依据：`fin/pipeline/assignment.py:1023`。
只有明确覆盖且能区分的 read 给 winning structure 的 covered_vote 加 1；缺少明确 coverage 不参与 prior 计票。
若某处 B 得 3 票、C 得 1 票，一条模糊 tie 的 structural prior 就是 B=3/4、C=1/4，对应代价 `−ln(0.75)=0.287682`、`−ln(0.25)=1.386294`。
这类代价是**票数 prior 的编码**，不是原始 junction NLL 加权项；没有任何 covered votes 时退回 flat tie。依据：`fin/pipeline/assignment.py:1039`。
如果 B 的结构有两个 endpoint sibling，0.75 再平分为各 0.375，保持结构总份额；有票的 tie 中零票结构保持 mask。依据：`fin/pipeline/assignment.py:1052`。

### 4.3 polyA 3′：只在仍 eligible 的 tie cell 加 NLL

端点通道要求 raw tie 至少两个候选，而且剩余 `d_tx<1e5` 的 eligible cell 至少两个；仅一个没被 junction mask 掉时，不再做端点竞争。依据：`fin/pipeline/assignment.py:657`。
r2 的 A/A.ext 都还 eligible，PASS、length=42，端点核可用，因此二者分别得到距离 +180、−3 的 NLL。
离散 Laplace 对绝对距离对称，`P(−3)=P(3)`；两项是 8.658973、3.148363。
代码减去本行新增端点代价中的最小值，得到 A `5.510610`、A.ext `0`；这只是共同平移，softmax 结果不变。依据：`fin/pipeline/assignment.py:670`。
这里没有把长度 42 当成额外权重，也没有要求 tail 比 40 nt 更长；默认 admission 只要求它的 PASS 长度为正。依据：`fin/analysis/molecule_evidence.py:51`。

### 4.4 upstream 5′：普通通道逐 read 计算，不要求 start 簇

r2 的 start=1043020，比 A/A.ext 的 TSS=1043000 晚 20 bp，所以普通 5′ 项为 0，不改变上一小节的两个代价。依据：`fin/pipeline/assignment.py:684`。
作为同 locus 的条件演算：某条 read 从 1043010 开始，若 B/C 都在 eligible tie，B 的 `u=0`、C 的 `u=1043030−1043010=20`。
教学 start 核给 C 额外 5.464365，给 B 0；即便只有这一条 read 从 1043010 开始，这个普通 assignment penalty 也照算。
若 read 从 1043700 开始，对这两个 TSS 都在下游，两者增加 0；不能据此说它一定是一个新 TSS。
这正是 discovery 的“至少三条 start 支持”与 assignment 的“每条 read 上游矛盾”两个规则的区别。依据：`fin/pipeline/runner.py:1956`、`fin/pipeline/assignment.py:684`。

### 4.5 softmax 与 abundance：这里的 EM 不是 abundance 迭代

生产调用 generic engine 时 `sigma=1`、`beta=0`，`abundance_feedback=False`；即固定 cost 的逐行 softmax，不迭代 theta 把模糊 read 越推越偏。依据：`fin/pipeline/assignment.py:1082`、`fin/cli.py:261`。
r2 的责任度为 `R(A)=exp(−5.510610)/(1+exp(−5.510610))=0.004027353`。
`R(A.ext)=0.995972647`，B、C 不在 tie，所以为 0；一行和为 1。
covered-vote prior 可能来自别的 read，但它只构造一次 cost，不是从当前 abundance 反复更新的反馈循环。依据：`fin/pipeline/assignment.py:1039`。

为了把账带到最后，设完整 14 条 read 的 assignment 汇总如下；其余行的 AS/信号假设使所示 winner 合法。

| read 集合 | read 数 | 每条给 A | 每条给 A.ext | 每条给 B | 每条给 C |
| --- | ---: | ---: | ---: | ---: | ---: |
| r1 | 1 | 1 | 0 | 0 | 0 |
| r2 | 1 | 0.004027353 | 0.995972647 | 0 | 0 |
| r6,r7,r8,r9,r10 | 5 | 0 | 1 | 0 | 0 |
| r3,r4,r11,r12,r13 | 5 | 0 | 0 | 1 | 0 |
| r5,r14 | 2 | 0 | 0 | 0 | 1 |

列和是 A=1.004027353、A.ext=5.995972647、B=5、C=2，总和 14；这是 abundance，允许是小数。依据：`fin/analysis/quantification.py:94`。
hard winner 数则是 1、6、5、2；例如 r2 只算 A.ext 的一个 hard read，却同时贡献 A 的少量 soft mass。
r10 原始链是 B 但这里分给 A.ext，说明 BAM 原始链归组不锁死后续候选归属。

## 5. Assignment 后：还没有离开区间

### 5.1 folded-TSS shadow recovery：先问短链是否真有自己的起点

生产 `folded_tss_recovery_mode="require"` 在第一次 assignment 后检查 dormant shadow；不是把所有 shadow 都自动救回。依据：`fin/pipeline/config.py:71`、`fin/analysis/tss_shadow_recovery.py:175`。
需要 parent 在第一次候选轴上且 soft abundance>0；短链是纯 5′ truncation，正链为 exact suffix、负链为 exact prefix；至少 3 条 distinct shadow reads，其中≥80% 的 TES 距 parent 3′ 端≤25 bp。依据：`fin/analysis/tss_shadow_recovery.py:185`。
本例 parent B mass=5>0，r6 的链确实是 suffix，TES 距 B 为 `|1045181−1045183|=2 bp`。
但只有 1 条 shadow read，`1<3`，因此直接记录 too-few-shadow-reads，不拟造一个“1043700 新 TSS”，也不重跑 assignment。

用同一几何做一个**独立的高覆盖分支演算**：假设有 20 条 shadow read，18 条末端在±25 bp，`18/20=90%≥80%`，start 中位数=1043700。
在 B spliced frame 中该 start offset 为第一 exon 的 `1043454−1043010=444` 加第二 exon 的 `1043700−1043600=100`，合计 544 bp。
25 bp 分箱得到 `544//25=21`；parent spliced length=1427，搜索校正箱数是 `1427//25=57`。
测试 read 来自第一次对 parent 的 `R≥0.05` 行，每条计一次、不按 R 加权；背景取 parent hazard，并用局部背景下限保护降解热点。依据：`fin/analysis/tss_shadow_recovery.py:210`。
假设该箱 `n_at_risk=50`、`n_peak=20`、上游=30、MAD=2 bp，pooled 与 local 取大后的 hazard=0.02。
则 peak fraction=`20/50=0.4`，背景期望 `50×0.02=1`；满足至少 10 at-risk、至少 3 peak、peak fraction≥0.10、MAD≤12 的条件。依据：`fin/analysis/tss_shadow_recovery.py:232`。
若 2000 次 null bootstrap 中达到 20 的次数为 0，则 pointwise p=`1/2001`，搜索校正后 p=`57/2001≈0.028486≤0.05`，该教学分支通过。依据：`fin/analysis/tss_evidence.py:497`。
测试还有“没有上游 read、起点高度集中”的结构性支持分支，不应把所有 supported 都解释为通过了 bootstrap；本例高覆盖分支有上游 read，不走那条。依据：`fin/analysis/tss_evidence.py:441`。

支持后发出 novel 短链，5′ 取 in-bin shadow reads 中位数、3′ 取 parent 的终点、继承 family，仅 in-bin reads 作为生成支持；新增候选才触发一次第二轮 assignment，不递归恢复。依据：`fin/analysis/tss_shadow_recovery.py:27`、`fin/pipeline/runner.py:2425`。
若 recovered identity 已经存在，只给已有模型合入支持，不因为这次合入就必然重跑；例如已存在模型原有 5 条支持，新找到 3 条不同 read，支持 union 变 8，而不是再发出一个相同模型。依据：`fin/analysis/tss_shadow_recovery.py:283`、`fin/pipeline/runner.py:2419`。
第二轮带 recovery-specific likelihood：只在 raw tie **恰好** `{parent,recovered}` 时作用；在支持起点箱中给 parent 加 `−ln(hazard)`，给 recovered 加 0。依据：`fin/pipeline/assignment.py:708`。
例如 frozen hazard=0.02，parent 代价加 `−ln(0.02)=3.912023`，其他项相同则 responsibility 为 parent `0.02/1.02≈0.019608`、short `1/1.02≈0.980392`。
在 short 起点箱上游的 read 会 mask short；下游 read 此专用项弃权；含第三候选的 tie 不应用它。
专用对比还写入 support gate 使用的 NLL 字典，防止“刚恢复又因没有 junction 对比而被支持门误删”。普通 ends 项没有这种相同的字典同步。依据：`fin/pipeline/assignment.py:727`、`fin/pipeline/assignment.py:617`。

### 5.2 区间选择：先 resolve extent，再做支持和 containment，最后 mono

实际顺序不是“量化完直接全样本筛选”：`select_m2_interval()` 首先 resolve A/A.ext，随后 structural/support/containment gates，末尾才做 mono resolution。依据：`fin/pipeline/selection.py:94`、`fin/pipeline/selection.py:468`。
A.ext 的 5.995972647 大于 A 的 1.004027353，所以留下 A.ext；两列相加得到 `7.000000000`，hard read union 也是 7。
如果两者恰好各 3.5，选择注释 span A；不是默认偏向更长的 extent。依据：`fin/pipeline/selection.py:518`。
resolve 同时重指向 R、ties、NLL 和 fold redirect；后面的门看到一个汇集证据的 extent 模型，而不是两个互相分走支持的候选。

cluster recheck 是额外的结构近邻复核，默认窗口 20 bp、fraction=0.15，不要与 discovery family 的 6 bp 搞混。依据：`fin/pipeline/config.py:611`。
本例最强模型质量 7，低量参考线是 `7×0.15=1.05`；B=5、C=2 都高于它，故本例不因低 fraction 被清除；GTF 还存在直接 junction read-support guard。依据：`fin/pipeline/selection.py:225`。
novel junction support 门要求每个 junction 至少 2 个直接观测 read，坐标容差默认 2 bp；GTF/fusion/mono 豁免。依据：`fin/pipeline/selection.py:259`。
B 第一 junction 6 条、C 第一 junction 2 条，都合格；若 C 只有 r5 一条，它会在 snap 前被删，不能用后面的 snap 阈值替它“预先救活”。
M2/M1 support 门要求 multi 至少赢得一条 read 的 sole-best AS，或一条 scored tie 的最低 NLL；默认最低 NLL 并列也算，GTF 不整体豁免。依据：`fin/scoring/m2_junction_nll.py:458`。
本例为 A.ext、B、C 各设有 sole-best 或最低 NLL read，因而通过；单有 0.2 soft mass 不等于满足这道证据门。

containment-cluster 门针对 novel 严格较短子链，用 abundance ratio 与 discovery read-count ratio **同时**判断，并排除已被删除的 parent。依据：`fin/pipeline/selection.py:379`。
例如另一个短链模型 abundance=2、parent=10，生成支持=2、parent=10：`2≤10×0.3=3` 两项都成立，且短链 read 数 2≤默认 cap 10，可被删除。
若短链支持=11，则即使两项比例低，仍超过该 cap；若 read ratio=4/10=0.4，也不满足双低条件。依据：`fin/scoring/m2_junction_nll.py:904`、`fin/pipeline/config.py:674`。
本例 B、C 都与 A.ext 有两个 intron，不是更短子链，因此这道门不把它们合并。

mono resolution 检查存活 mono 模型的 hard reads 是否落在存活 multi 的单个 exon 内，链向一致、带默认 10 bp terminal slop；多个可容纳 parent 中选 abundance 最大者。依据：`fin/pipeline/selection.py:468`。
独立例子：mono 的四条 read 中三条分别落在 `[1044600,1044700)`、`[1044650,1044800)`、`[1044700,1044900)`，被 A.ext 的末 exon 容纳，另有一条不被任何 multi 容纳。
三条 fragment 整条归 parent，不按 0.5/0.5 分；剩下 mono 只有 1 条，`1<2`，达不到区间 mono floor，被删除。
这与后面 global novel-mono 的 5-read floor 是两道不同规则；本例主线没有 mono，14 条主账不额外增加这四条。依据：`fin/pipeline/config.py:196`、`fin/pipeline/config.py:48`。

## 6. 全样本 finalize：身份变化以后，重新结清 read 质量

### 6.1 unify extent pairs，然后 aggregate

跨区间先 `unify_extent_pairs()`，再 `aggregate_across_intervals()`；不是先把两个 extent ID 永久当成两个转录本累计。依据：`fin/pipeline/runner.py:243`。
假设其他 locus 的同一注释在区间 I 保留 A mass=2，在区间 II 保留 A.ext mass=3，全样本比较 `3>2`，把 A 的结果身份重标为 A.ext，并保留 redirect；总量相等时保留 A。依据：`fin/pipeline/selection.py:639`。
本例只有一个相关区间，已留 A.ext=7，因此 unify 没有新变化，B=5、C=2 也不动。
aggregate 按 candidate_id 汇总，不按“同 gene”或“同 chain”合并；对于重复 hard read ID，取 union 并用 `unique/hard-observation-weight` 缩放累计 abundance。依据：`fin/analysis/quantification.py:533`、`fin/analysis/quantification.py:652`。
独立数字例子：I 有 `{u1,u2}`、abundance=1.8；II 有 `{u2,u3}`、abundance=1.6，hard observation 数为 4、unique 为 3。
aggregate 给出 `(1.8+1.6)×3/4=2.55`，不是 3.4，也不是把软责任度按 read 精确重算；精确 survivor 归一化要等 ledger refit。
没有 hard ID 的结果有回退分支，不能把这个 dedup 比例机械套到所有零-hard-read 候选。依据：`fin/analysis/quantification.py:659`。

### 6.2 global select：选择基于 refit 前的证据

生产 novel abundance floor 是 **严格 >1.0**，GTF 是 **≥1.0**，fusion 豁免；不启用人为 GTF abundance floor。依据：`fin/pipeline/config.py:38`、`fin/pipeline/selection.py:800`。
因此 novel 1.0000 删除，1.0001 保留；GTF 1.0000 保留，0.9999 删除。
A.ext 是 GTF、质量 7；B=5、C=2 是 novel，都通过。不要因 ID 带 `.ext` 就把 A.ext 用 novel 严格阈值处理。
novel multi 还需达到 family 内 dominant abundance 的 0.01；没有 family 时有 overlap fallback，GTF/fusion/mono 豁免这道 fraction gate。依据：`fin/pipeline/selection.py:845`。
本例 dominant=7，阈值=`7×0.01=0.07`，B=5、C=2 都远高于它；分母不是所有 isoform 的总 abundance=14。
若某 family dominant=200，一个 novel abundance=1.5 虽大于 1，也因 `1.5/200=0.0075<0.01` 被删。
novel mono 另要求 hard reads≥5，默认长度下限=0；所以 4 条失败、5 条可过，开关名 `drop_mono_exon_novel` 不表示无条件删除所有 novel mono。依据：`fin/pipeline/selection.py:932`、`fin/pipeline/config.py:686`。
本 profile 的 soft/hard ratio、full-length fraction、polyA+5′ retention gate 都关闭；严格正质量 GTF 特殊门也没启用，不能和普通 GTF≥1 混为一谈。依据：`fin/pipeline/config.py:38`、`fin/pipeline/config.py:214`。

### 6.3 junction snap：更强的观测 junction 才能替换 novel 坐标

全样本重新扫描 primary BAM 的精确 junction，按 chromosome/strand 计数；snap 仅处理 novel 多外显子，外侧 TSS/TES 不动。依据：`fin/pipeline/junction_snap.py:19`、`fin/pipeline/junction_snap.py:217`。
每个目标 junction 必须 donor 和 acceptor **分别**距当前≤6 bp，满足 canonical 要求；按支持数降序、总移动距离、坐标选最优。依据：`fin/pipeline/junction_snap.py:51`。
替换需要 `target_count≥max(min_support,floor(current_count×ratio)+1)`，默认 min_support=2、ratio=2.0；恰好两倍不够。依据：`fin/pipeline/config.py:50`、`fin/pipeline/junction_snap.py:82`。
本例 C 第一 junction `(1043451,1043600)` 有 2 条，B 的 `(1043454,1043600)` 有 6 条，A 的坐标有 5 条。
到 B donor 移动 `|1043454−1043451|=3`、acceptor 移动 0；目标 B 比目标 A 的 5 条更强，先选 B。
要求 `max(2,floor(2×2.0)+1)=max(2,5)=5`，而 `6≥5`，所以 C 第一 donor 从 1043451 改成 1043454。
B 自己当前 6 条，要更强替换需 `floor(6×2)+1=13`，没有这样的目标，因此保持；第二 junction 本来相同，也不动。
如果 current=1、target=6，则需要 `max(2,floor(1×2)+1)=3`，能通过；如果 current=3、target=6，则需要 7，不能通过。
即便计数通过，还要保持前 exon start < donor < acceptor < 后 exon end，防止产生无效 exon。依据：`fin/pipeline/junction_snap.py:89`。

snap 后 B exons 不变；C exons 变成 `(1043030,1043454),(1043600,1043900),(1044500,1045183)`。
二者现在 intron chain 相同，但 TSS 仍相差 20 bp；snap 没有把 C 的 start 改成 B 的 start。
A.ext donor 1043457 保持不动：虽然旁边 B 支持更强，它是 GTF，根本不进入 snap。依据：`fin/pipeline/junction_snap.py:217`。

### 6.4 第一种合并：完整 exon tuple 一样

snap 后先按 novel/source、chromosome、strand、完整 exon tuple 分组；GTF/fusion/mono 即便几何相同也保留独立身份。依据：`fin/pipeline/junction_snap.py:244`。
本例 B 第一 exon 起点=1043010，C=1043030，tuple 不一样，因此这第一轮**不合并** B/C。
若另一个教学 novel D 的 post-snap exons 与 B 完全相同，B mass=5、D mass=1.5，则这轮可得到 6.5；不是要求 TSS 差≤25，而是完整 tuple 相等。
代表由 hard reads 最大优先，再按 abundance、candidate_id 最大值决定；合并加 abundance、union hard IDs，保留代表几何。依据：`fin/pipeline/junction_snap.py:134`。
这里的 D 是独立算例，不加入本例 14 条 read 的主账。

### 6.5 第二种合并：snap-created chain collision 的八个条件

第二轮默认开启 `snap_merge_chain`；它允许有限 TSS 差，但必须证明关系是此次 snap 新产生的，而非本来就同链的 APA/alternative-TSS 模型。依据：`fin/pipeline/config.py:81`、`fin/pipeline/junction_snap.py:289`。
下面八项必须一起满足；每项都用 B/C 的通过情形和一个拒绝反例说明。

1. **所有原始成员都是 novel 多外显子。** B、C 各三 exon，合格；换成 A.ext 和 B 即失败，因为 A.ext 是 GTF，不因名字像“新模型”变 novel。依据：`fin/pipeline/junction_snap.py:244`、`fin/pipeline/junction_snap.py:314`。
2. **同 chromosome、同已知 strand。** B/C 都是 chr1、`+`，合格；C 若在 chr2、负链或 strand=`.`，不合格。依据：`fin/pipeline/junction_snap.py:314`。
3. **post-snap 完整 intron chain 相同。** B/C 都成为 `((1043454,1043600),(1043900,1044500))`；若 C 第二 acceptor 仍为 1044503，即使只差 3 bp，也不能按“近似同链”合并。依据：`fin/pipeline/junction_snap.py:317`。
4. **每个第一轮 survivor 的 pre-snap chain 集合非空且两两不相交。** B 的集合是 `{链B}`，C 是 `{链C}`，无交集；若两者在 snap 前就都叫链B，则失败，不能删掉已有端点差别。依据：`fin/pipeline/junction_snap.py:334`。
5. **至少一个原始成员确实被 snap 改过。** C donor 移动了 3 bp，合格；如果两者都原封不动，就不是此次 snap 创建的 collision。依据：`fin/pipeline/junction_snap.py:340`。
6. **所有原始成员只有一个非空 discovery family。** B/C 都来自 F，合格；若第一轮先合并了 F 和 G，不能只取代表显示的 F 来蒙混过关；任何原始成员 family 缺失也不合格。依据：`fin/pipeline/junction_snap.py:305`。
7. **链向 3′ 端完全相同。** B/C 都为 1045183，合格；若一个为 1045180，哪怕只差 3 bp 也拒绝，不是用 ±25 bp 合并 APA。负链的 3′ 端看 start。依据：`fin/pipeline/junction_snap.py:102`、`fin/pipeline/junction_snap.py:317`。
8. **5′ group diameter≤25 bp。** B/C 的 `1043030−1043010=20`，合格；若第三个起点=1043050，三者直径40，不能借相邻两两差20连锁合成一组。分组按排序后有界直径，而非 single-linkage。依据：`fin/pipeline/junction_snap.py:110`、`fin/pipeline/config.py:291`。

八项都满足后，本例 B hard=5、C hard=2，代表选 B，保留 B 自己的 post-snap 几何，而不是平均起点 1043020。
临时 abundance=`5+2=7`，hard union=7，redirect 为 C→B；chain-collision 合并把旧的 full-length fraction 标为 −1，不能假装旧端点指标仍适用于新 read 集合。依据：`fin/pipeline/junction_snap.py:134`、`fin/pipeline/junction_snap.py:350`。
因此当前两个输出身份是 A.ext=7 与 B=7；不是把 A.ext、B、C 三者都视为一个 transcript。

### 6.6 refit：重路由、归一化，而不是重新 mappy 比对

生产在 snap/merge 后用 responsibility ledger 重算：先沿 interval fold、extent、snap redirects 转移旧列权重，再在存活目标间逐 read 归一化。依据：`fin/pipeline/runner.py:438`、`fin/analysis/abundance_refit.py:143`。
本例 r2 原来 A=0.004027353、A.ext=0.995972647，A→A.ext 后相加为 1；归一化仍是 A.ext=1。
r5、r14 原来各给 C=1，C→B 后各给 B=1；没有删掉这两条 read，也没有凭模型合并再复制两条。
主账最终仍是 A.ext=7、B=7，总质量=14；本例所有原始质量都有存活去处，**refit 不一定改变总量**。

再看一个独立删列例子：某 read 原来 P=0.4、Q=0.2、被删 X=0.4，Q 又 merge 到 P。
路由后存活 P 的权重为 `0.4+0.2=0.6`，X 不在 survivor 集合，归一化 P=`0.6/0.6=1`；不是继续只报 0.6。
若两个 survivor 权重为 P=0.2、Q=0.3，另 0.5 被删，则 refit 后为 `0.2/0.5=0.4` 与 `0.3/0.5=0.6`。依据：`fin/analysis/abundance_refit.py:234`。
同一 read/candidate 的跨区间 ledger 权重先相加再归一化；不是把每次 interval observation 当一条新分子。依据：`fin/analysis/abundance_refit.py:185`。

mono-resolved read 不走普通软归一化，而跟随 forced target；target 存活就给它 1，所有 forced target 都消失时记 orphan，不随便转投其他模型。依据：`fin/analysis/abundance_refit.py:218`。
普通 read 若只剩被删除的目标，也记 selection-orphaned=1；这和“根本没进入有效 assignment 行”的 read 是不同诊断计数。依据：`fin/analysis/abundance_refit.py:250`、`fin/analysis/abundance_refit.py:338`。
例如 16 条 assignable reads 中 15 条有去处、1 条 orphan，质量校验是 `16=15+1`，不是强制给所有 survivor 分满16。
refit 重新算 abundance、hard IDs、confidence、max_R；forced reads 不参与 confidence 平均，不能把规则强制的 1 当成观测置信度。依据：`fin/analysis/abundance_refit.py:260`、`fin/analysis/abundance_refit.py:290`。
选择集合冻结于 refit 前，不再跑一轮普通 abundance/global gate；严格正质量 GTF 的后验清理函数在这里虽被调用，但本 profile 的 enabled 为 False。依据：`fin/analysis/abundance_refit.py:338`、`fin/pipeline/runner.py:460`、`fin/pipeline/config.py:214`。

### 6.7 写出：GTF 讲结构，TSV 讲分配结果

writer 前为 GTF source 按 `extent_of` 或候选 ID 查回 gene_id；没有 gene_id 的模型回退为自己的 candidate_id，不把“一个 family”强行当成“一个 gene”。依据：`fin/pipeline/finalize.py:94`。
假设 A 注释 gene_id 是 GENE1，A.ext 仍写 GENE1；B 没有注释 gene_id，则使用 B 的实际生成 ID。
GTF 按 gene/transcript/exon 层次写，内部 start 加 1，end 不变；transcript 的 strand 用模型自己的 strand。依据：`fin/io/io_gtf.py:563`、`fin/io/io_gtf.py:589`。
本例 B 的第一 exon `[1043010,1043454)` 写成 start=1043011、end=1043454，长度仍是 `1043454−1043011+1=444`。
GTF 九列依次是 chromosome、source、feature、start、end、score、strand、frame、attributes；本 writer source 默认 pyfin，score/frame 写 `.`。依据：`fin/io/io_gtf.py:550`、`fin/io/io_gtf.py:620`。
transcript attributes 包含 gene_id、transcript_id、abundance、confidence、num_reads、transcript_source 和三个旧 score；exon attributes 含 exon_number。依据：`fin/io/io_gtf.py:608`。
`transcript_source=gtf` 不保证端点等于原注释：A.ext 就是反例；GTF transcript attributes 不包含 TSV 的 TPM、max_R。依据：`fin/io/io_gtf.py:608`。

TSV 列的含义与本例数值如下，start/end 保留内部 0-based 表示。依据：`fin/io/io_tsv.py:19`、`fin/analysis/quantification.py:24`。

| 列 | A.ext 示例 | B 示例 | 含义 |
| --- | --- | --- | --- |
| candidate_id | A.ext | B | 实际运行中是完整模型 ID；这里用教学简称 |
| gene_id | GENE1 | B | 注释 gene 或缺失时的 ID 回退 |
| chrom | chr1 | chr1 | 基因组染色体 |
| strand | + | + | transcript 链向 |
| start | 1043000 | 1043010 | 0-based 左边界，不总是生物学 TSS |
| end | 1045183 | 1045183 | half-open 右边界，不总是生物学 TES |
| source | gtf | novel | 模型来源类别 |
| abundance | 7.0000 | 7.0000 | refit 后 responsibility 列和 |
| confidence | 1.0000 | 1.0000 | 非 forced hard-assigned reads 的平均 winning probability |
| coherence_score | 0.0000 | 0.0000 | 保留的旧 schema 字段 |
| discrimination_score | 0.0000 | 0.0000 | 保留的旧 schema 字段 |
| combined_score | 0.0000 | 0.0000 | 保留的旧 schema 字段 |
| num_reads | 7 | 7 | refit 后 hard-assigned read 数 |
| tpm | 497732.8218 | 502267.1782 | 下文假设样本只剩这两个模型时的 TPM |
| max_R | 1.0000 | 1.0000 | 任一 read 对此模型的最大 responsibility |

三个旧 score 在 assembly 路径不再被填充，0 不代表“信号特别差”；它们为输出 schema 稳定而保留。依据：`fin/analysis/quantification.py:42`。
A.ext spliced length=`457+300+683=1440 bp`，B=`444+300+683=1427 bp`；TPM 使用 spliced length，不用 2183/2173 bp genomic span。
仅有这两个模型时，A.ext RPK=`7/1.440=4.861111`、B RPK=`7/1.427≈4.905396`；各除以 RPK 总和再乘百万，得到表中 TPM。依据：`fin/analysis/quantification.py:500`。
真实全样本还有其他模型时，分母包括它们，不能把本 locus 的两个 TPM 自己归一到百万后当正式输出。
主线的 confidence 最终都为 1，是因为 resolve/merge 后每条 read 只剩一个有权重的身份；它不等于“这两个生物学 isoform 的存在概率都是 100%”。

## 7. 哪些默认关闭：不要误接到生产主线

以下一行一项；值为 0 的 gate 表示禁用，而不是阈值已通过。

- `abundance_feedback=False`：不迭代 theta；`abundance_length_norm=False` 也不参与分配反馈。依据：`fin/cli.py:261`。
- `junction_correction_mode="off"`：不在 discovery 前逐 read 改写 junction chain。依据：`fin/pipeline/config.py:208`。
- `junction_election_mode="off"`：不跑另一个候选级 junction election 决议层。依据：`fin/pipeline/config.py:204`。
- `containment_evidence_mode="off"`：额外的 TES/TSS containment evidence layer 关闭，不等于 production containment-cluster 关闭。依据：`fin/pipeline/config.py:207`。
- `complete_chain_rescue_reads=0`：不使用完整链 clean-depth 的特别救援豁免。依据：`fin/pipeline/config.py:203`。
- `require_gtf_positive_mass=False`：不启用额外严格正质量 GTF 门。依据：`fin/pipeline/config.py:214`。
- `max_soft_mass_ratio=0`：不按 soft abundance/hard reads 比值过滤。依据：`fin/pipeline/config.py:45`。
- `min_fulllen_fraction=0`：不按 full-length fraction 过滤。依据：`fin/pipeline/config.py:46`。
- `min_polya5p_reads=0`：候选 retention 的 polyA+5′ 门关闭，不等于逐 read 的 ends 通道关闭。依据：`fin/pipeline/config.py:47`。
- `floor_gtf_abundance=False`：不人为把 GTF 质量抬到 abundance floor。依据：`fin/pipeline/config.py:42`。
- `m2_support_gate_gtf_exempt=False`：没有额外 annotation-wide support-gate 豁免。依据：`fin/pipeline/config.py:368`。
- `polya_end_accept_adapter=False`：默认不把 ADAPTER 当 PASS 端点证据。依据：`fin/pipeline/config.py:342`。
- `endpoint_qc_conditional=False`：不扩大到按 QC 类别分别拟合/使用端点核。依据：`fin/pipeline/config.py:353`。
- `endpoint_refine=False`：不在 finalize 另跑 endpoint states 的拆分流程。依据：`fin/pipeline/config.py:292`。
- `candidate_existence_mode="off"`：统一 existence audit 不作为默认流程。依据：`fin/pipeline/config.py:319`。
- `ranking_mode="off"`、`candidate_evidence=False`：不启用冻结模型排名过滤及其可选候选诊断表。依据：`fin/pipeline/config.py:243`。
- `catalog_mode="off"`：不为 benchmark catalog contract 补齐未观察注释。依据：`fin/pipeline/config.py:389`。
- `containment_collapse` 默认关闭：不是 production 的 discovery exact-subchain shadow，也不是默认 containment-cluster gate。依据：`fin/pipeline/selection.py:411`。
- `guided_junction_min_reads` 默认关闭：不是 novel 每 junction 至少 2 reads 的生产门。依据：`fin/cli.py:232`。
- `junction_dominance_filter` 默认关闭：registry 记录不自动变成 pre-EM dominance gate。依据：`fin/cli.py:243`。

相反，folded-TSS recovery=require、canonical gate、junction snap、snap-chain merge、post-selection refit 都是此 profile 的默认 ON，不能误放进“仅供实验”的清单。依据：`fin/pipeline/config.py:38`。
生成时 mono fold 的原始配置虽为 True，但 production post-EM mono resolution 使它不在生成阶段生效；这是被另一阶段取代，不是两次都折。依据：`fin/pipeline/config.py:178`。

## 8. 最容易误解的五件事

1. **同链不等于同一个 transcript，family 更不等于同链。** r1/r2 精确链相同、端点不同；A/A.ext 是明确关联的两种 extent；B/C 是 snap 才造成的同链 collision。它们分别走精确分组、extent resolve、八条件 merge，不存在一个通用“同链全合并”操作。依据：`fin/candidates/discovery.py:337`、`fin/pipeline/selection.py:518`、`fin/pipeline/junction_snap.py:334`。

2. **名字叫 EM，不表示默认在迭代 abundance theta。** r2 从两个固定 cost 算出 0.004027/0.995973，所有 read 列和才是 abundance；covered-vote prior、条件性第二次 assignment、final ledger refit 也不是 theta feedback。依据：`fin/pipeline/assignment.py:1039`、`fin/pipeline/assignment.py:1086`、`fin/analysis/abundance_refit.py:143`。

3. **普通 5′ assignment 不需要 start 簇；discovery contradiction 与 shadow recovery 才有自己的群体要求。** 一个 upstream read 可以给晚起点候选加 penalty；一个下游开始的 r6 却不足以建立新 TSS。不能把三种 5′ 用法读成同一个阈值。依据：`fin/pipeline/assignment.py:684`、`fin/pipeline/runner.py:1956`、`fin/analysis/tss_shadow_recovery.py:197`。

4. **recall 不变是实测，不是算法保证。** 仓库记录 snap-chain merge 在评测网格上 recall 保留到两位小数，并有 SIRV4、HEYA8、人类 GENCODE、Sequin 的提升/不变结果；代码没有查询 reference match 来证明合并不损伤任何真实 isoform。八条件只是操作性准入规则，换样本仍需验证。依据：`fin/pipeline/config.py:72`、`fin/pipeline/junction_snap.py:197`。

5. **snap 只动 novel 的内部 junction，不动 GTF，也不负责端点重估。** 本例 C donor 改 3 bp，B/C 通过八条件后用代表 B 的 TSS；A.ext 的 TES 变长发生在 discovery 并由 extent posterior 决议。最后 abundance 的账则由 refit 结清，不能拿“snap 前 5+2”代替所有情况下的最终 abundance。依据：`fin/candidates/discovery.py:299`、`fin/pipeline/junction_snap.py:217`、`fin/analysis/abundance_refit.py:234`。
