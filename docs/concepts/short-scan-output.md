# `short scan` 的输出表怎么读

本页解释 `alnkit short scan` 产出的每一个文件、每一列的含义，以及**最容易读错的两张表**
（`Insert_Interference` 与 `site_homology`）。示例一律用占位名（`sampleA`、`chr1`、
`insertA`/`insertB`），和你自己的产物对照时请换成实际名字。

其他命令的表（长读长链路、绘图）逐列解释见[各命令的输出表](output-tables.md)。

一句话版本：

| 你想回答的问题 | 看哪张表 |
|---|---|
| 位点在哪、有多少证据 | `<sample>.Insert_Position.result.xls` 或 `<sample>.sites.canonical.tsv` |
| 机器接口 / 和下游 join | `<sample>.sites.canonical.tsv`（列名带基准后缀） |
| 每条支持 read 的比对长什么样 | `<sample>.reads_info.out.xls` |
| 证据是否与**别的插入序列**共享（多外源面板） | `<sample>.Insert_Interference.result.xls` |
| 这个位点在基因组里是否**与插入序列同源**（外源 ↔ 内源） | `<sample>.site_homology.tsv`（要 `--reference-fasta` 或 `--homology-bed`） |

## 1. 产物清单

| 文件 | 何时产出 | 格式 |
|---|---|---|
| `<sample>.Insert_Position.result.xls` | 总是 | 8 列，TAB 分隔（历史格式，表头见 §3） |
| `<sample>.reads_info.out.xls` | 总是 | 表头 11 列、**数据 12 列**（历史缺陷，见 §4） |
| `<sample>.sites.canonical.tsv` | 总是 | 注册表 `insert_sites`（11 列，机器接口） |
| `<sample>.Insert_Interference.result.xls` | 总是 | 15 列，TAB 分隔（历史格式） |
| `<sample>.site_homology.tsv` | 仅在给了 `--reference-fasta` 或 `--homology-bed` 时 | 注册表 `site_homology`（12 列） |

文件名里的 `.xls` 只是历史命名，内容是纯文本 TSV。`--sample-name` 决定 `<sample>` 前缀。

## 2. 数据流（一行一个阶段）

```
BAM ──过滤(read 级)：flag/mapQ/错配/soft-clip 段数/候选闸门 ──► 候选 read
     ──归属：把 read 判给 query / 别的插入 / 不明确 / 没有插入命中
     ──证据：split read（断点处的 soft clip）+ PESR（双端同断点）
     ──按 (chrom, 断点) 聚合 ──► BreakpointEvidence（含各种 read 名单 + 干扰源头）
     ──窗口合并（--merge-window-bp）──► 代表位点
     ──► 三张 legacy 表 + canonical 表（+ 可选 site_homology）
```

窗口合并只吸收 **clip 端点抖动（几个 bp）**，默认 0（不合并）。`--merge-window-bp 10`
意味着相距 ≤10 bp 的断点算一个位点；把它调到几百 bp 会把同一座位两侧的两个连接点并成一行，
**结构信息就没了**（两侧是不同方向的 clip，见 §7 的误读清单）。

## 3. `<sample>.Insert_Position.result.xls`（位点主表）

| 列 | 含义 |
|---|---|
| `Chromosome` | 参考序列名 |
| `Insert_Position` | 断点边界（0-based 数值 = 断点左侧碱基的 1-based 坐标，两者数值相同） |
| `split_read_support` | SR 支持 read 数（断点处有 soft clip 的 read） |
| `pair_split_read_support` | PESR 支持数（双端分别给出同一断点） |
| `dup_read_count` | 被 markdup 标记、但落在这个位点的 read 数（**默认不作为证据**，只报数） |
| `split_read_name` / `pair_split_read_name` / `dup_read_name` | 对应的 read 名单，`;` 分隔；空列表示 0 条 |

读法：**先看 `split_read_support` 的量级**（真位点通常比邻居高一个数量级），再看
`split_read_name` 里是否是同一批 read 在多个位点重复出现（`--merge-window-bp` 该调大还是调小，
就看这个）。`dup_read_count` 不是“证据数”，它是**质量信号**：某位点 dup 很高说明这段区域重复。

## 4. `<sample>.reads_info.out.xls`（逐比对明细）

历史缺陷（**刻意保留**，改了就破坏与旧产物的逐字节对账）：表头是 11 个字段，最后一项写成
`stand,mapping_quality`（一个字段里塞两列），**数据行实际是 12 列**。请按**列序**解析：

| 序 | 列 | 含义 |
|---|---|---|
| 1 | `read_name` | read 名 |
| 2 | `read_pair` | `read1` / `read2` |
| 3 | `cigarstring` | 原始 CIGAR |
| 4 | `is_supplementary` | 该记录是否 supplementary（`True`/`False`） |
| 5 | `chrom` | 这条记录比对到哪条参考 |
| 6–7 | `M_ref_pos_start` / `M_ref_pos_end` | 证据块（M 块）的参考起点/终点 |
| 8–9 | `M_read_pos_start` / `M_read_pos_end` | 同一块在 read 轴上的起止 |
| 10 | `M_var_len` | 该块长度 |
| 11 | `stand` | `+` / `-` |
| 12 | `mapping_quality` | MAPQ |

用途：排查“这条 read 为什么被算成证据 / 为什么被丢弃”。它只含**通过质量过滤**的比对，
所以被丢掉的 read 不会出现在这里——丢弃原因看日志里的 `statistics` 行。

## 5. `<sample>.sites.canonical.tsv`（规范位点表，机器接口）

列定义见[表契约](tables.md)的 `insert_sites`。与主表
（`<sample>.Insert_Position.result.xls`）**按 `(chrom, insert_pos_1based)` 一一对应**——同样的位点、
同样的合并结果。恒等式：

```
insert_pos_1based - 1 == insert_pos_0based
insert_boundary_0based == insert_pos_1based        # 边界值 = 左侧碱基的 1-based 坐标
```

表里**不带 `#` 注释行**（0.7.0 起），第一行就是表头，可以直接 `pandas.read_csv(sep="\t")`。

`nearest_other_insert` / `other_inserts` 来自干扰判定，但**明确把 query 自己排除在外**：
“query 出现”只说明这条 read 也命中了本次的 query，不算竞争者（见 §6.3）。

## 5.1 `<sample>.site_homology.tsv`（可选：位点是不是内源同源）

只有给了 `--reference-fasta` 或 `--homology-bed` 才生成。列定义见[表契约](tables.md)的
`site_homology`；这里说怎么读。表里同样**不带 `#` 注释行**。

| 列 | 读法 |
|---|---|
| `is_homology_site` | `yes` = 该位点落在一条“与插入序列同源的参考区域”里（两侧各放宽 50 bp）；`no` = 没落在任何同源区域 |
| `region_sources` | 命中来源，`;` 分隔：`insert:<插入名>` = 把**插入序列**回比该位点周围的参考窗口得到的区域；`bed:<文件名/名字>` = 你给的 BED 区域 |
| `n_regions` | 命中的区域个数（多个来源重叠时会 >1） |
| `best_region_*` | 重叠区域里“最可信”的那一段：`insert:` 来源按**一致率**取最高，BED 来源按**长度**取最长；坐标 1-based 闭区间 |
| `best_region_identity` | 回比得到的一致率（0-1，按 `cs` 标签的真实差异算）；BED 来源没有一致率，填 `-` |
| `best_region_query_coverage` | 插入序列被这段区域覆盖的比例（0-1）；cDNA 样插入隔着内含子时这个数会明显 <1；BED 来源填 `-` |

判读要点：

- `region_sources` 里出现 `insert:<query>` 是**正常且期望**的：说明这条外源序列在基因组里有内源拷贝，
  而该位点正好落在拷贝上 → 这个位点多半是内源干扰，不是整合位点；
- 反过来，真整合位点应该 `no`（插入序列在那附近比不中基因组）；
- 这条判据与 §6 的干扰表**互相独立**：干扰表只比较“插入之间”，看不到“外源 ↔ 基因组同源”；
  外源与内源同源时，**只有这张表能分开真假位点**。

## 6. `<sample>.Insert_Interference.result.xls`（干扰诊断，**最容易读错**）

### 6.1 列

| 列 | 含义 |
|---|---|
| `Chromosome` / `Insert_Position` | 与主表同口径的位点 |
| `merged_position_count` / `merged_positions` | 这一行由哪几个原始断点合并而来（窗口 0 时是它自己） |
| `query_insert` | 本次运行的 `--query-chr`（整表同值） |
| `is_homology_interference` | 该位点**是否存在**干扰源头（= 下面那个名单非空） |
| `interference_source_count` | 干扰源头的**插入序列个数** |
| `interference_source_support_reads` | 每个源头各有几条支持 read，格式 `源头:条数;源头:条数` |
| `support_reads_total` | 该位点支持 read 总数（下面五项之和） |
| `support_reads_owned_query` | 归属 query 的 read 数 |
| `support_reads_owned_other` | 归属**其它**插入序列的 read 数（0.7.0 起恒为 0，见 §6.5） |
| `support_reads_consensus_assigned` | 原本不明确、被同断点投票推定归属的 read 数（`consensus` 策略） |
| `support_reads_ambiguous` | 归属不明确的 read 数 |
| `support_reads_no_insert_hit` | 完全没有插入命中的 read 数（0.7.0 起恒为 0，见 §6.5） |
| `ambiguous_read_names` | 上面“不明确”那批 read 的名字，`;` 分隔 |

### 6.2 判定链（query 永远不算源头）

1. 每条比对贡献“插入命中”：**自身命中**（source=`alignment`/`secondary`）＋ **XA 标签**里的
   备选命中（source=`xa`）；
2. “共享信号”只取 **secondary 与 XA** 命中（`alignment` 不算）。注意：**supplementary 属于
   `alignment`**，所以“断点 read 的另一段落在插入序列上”（正常形态）不会触发干扰；
3. 汇总前先把 **query 从命中集合与候选集合里剔除**（0.7.0 起恒如此），然后：
   - 还有其它插入命中 → 这些就是源头；
   - 归属不明确（`ambiguous`）→ 该 read 命中过的**其它**插入都算源头；
   - 都没有 → 该 read 不贡献源头。

> **历史：0.6.0 及更早**有一个 `--interference-scope {all,other}`。`all`（当时的默认）不剔除
> query，于是“断点 read 在插入面板上有次要命中”这条**正常形态**也会被记成源头——
> 单外源面板下真位点与内源同源位点都变成 `True`，`interference_source_support_reads` 里还会
> 出现 query 自己的名字（`insertA:3`）。这个口径无法区分生物学差异，已删除；下面 §6.3 说的
> 都不再受它影响。

### 6.3 它的价值：**不是“基因组同源性”**

这个表名里的 homology 指的是**“证据与插入面板共享”**（同源 read 在多个插入序列间横跳），
**不是**“外源序列与基因组同源”。两者后果完全不同：

| 场景 | 该看什么 | 为什么 |
|---|---|---|
| 多外源面板（`insertA`/`insertB`/`insertC`），某些位点的 read 同时命中别的外源 | 本表 | 这就是它擅长的：挑出“跨外源共享证据”的位点 |
| 外源序列与**基因组内源**序列同源（例如 query 是某内源基因的 cDNA 拷贝） | `<sample>.site_homology.tsv` | 判据必须是“把插入序列回比基因组”，而不是“read 是否命中插入面板” |

**单外源面板下这张表恒为空**（`is_homology_interference` 全 `False`、源头计数全 0）：
面板里只有 query 一条插入，剔除 query 后没有可比较的对象。这不是 bug，也不是“没有干扰”——
这类数据的干扰在**外源 ↔ 基因组**之间，本表按定义看不到。

结论：**不要用 `is_homology_interference` 做同源过滤**。要排除“外源 ↔ 内源同源”造成的假位点，
请用 `--reference-fasta`（插入序列回比位点周围窗口）或 `--homology-bed`（已知同源区域），
看 `site_homology.tsv` 的 `is_homology_site`。

### 6.4 同一个名字，两张表里含义不同

- 本表 `interference_source_*` 里出现 `insertB`（**另一条**插入）＝ 证据与别的外源共享，
  该位点可能来自外源之间的同源/共用片段；
- `site_homology.tsv` 的 `region_sources` 里出现 `insert:insertA`（**query 自己**）＝
  **同源证据**：这段参考序列与插入序列一致（内源拷贝就在那里）。这是**正常且期望**的标记。

一张表只谈“插入之间”，另一张才谈“插入 ↔ 基因组”；结论不要互相印证。

### 6.5 其他注意

- `interference_source_count` 数的是**源头个数**，不是 read 数；read 数在
  `interference_source_support_reads` 里；
- `support_reads_owned_other` / `support_reads_no_insert_hit` 自 0.7.0 起**恒为 0**（这两类 read
  在过滤阶段就被丢弃，不会再进到位点的支持 read 里）；列保留只是为了表格式不变；
- 干扰判定与过滤策略无关：`support_reads_*` 各项恒等于该位点的归属分布，
  `--insert-origin-policy` 只影响哪些 read 会被当作证据（见 §7）。

## 7. 归属策略：`--insert-origin-policy`（哪些 read 才算证据）

一条 read 往往同时命中好几条序列（外源、内源同源、拼接旁系同源）。这个参数决定**哪些 read
被当成“这个位点的证据”**，因此它直接改的是 `split_read_support` / `pair_split_read_support`
（主表）与 `support_reads_*`（干扰表），注意它**不改**已保留 read 的坐标。

判定顺序（`shortread/origin.py`）：先按“比对块长度 / 错配 / clip 长度”给每个插入序列打分，
得分唯一 → 归属该插入；并列 → `ambiguous`；一个插入都没命中 → `no_insert_hit`。

| 策略 | 保留谁 | 丢掉谁 | 什么时候用 |
|---|---|---|---|
| `owned`（默认） | 归属 **query** 的 read | 归属其它插入、`ambiguous`、`no_insert_hit` | 绝大多数情况从这里开始 |
| `consensus` | `owned` + 用**同断点投票**推定 `ambiguous` 的归属（附近 ±`--merge-window-bp` 内归属明确的 read 一致指向 query 才保留） | 投票为空或指向不一致的 ambiguous | 归属不明的 read 很多、又不想白白丢掉证据时（**推断不是测量**） |

另有一个**可叠加的开关**：

| 开关 | 作用 | 说明 |
|---|---|---|
| `--drop-shared-reads` | 在 `owned`/`consensus` 之上，再剔除“归属 query、但在**其它**插入上还有次要/XA 命中”的共有 read | 外源之间高度相似（共享接头/载体）时，避免共有 read 把特异位点撑起来 |

### 7.1 推荐怎么选（按数据形态，不按个人偏好）

**三步走**，不要一上来就换策略：

```bash
# 第 1 步：先用默认策略跑一遍，只看日志里这一行（统计会只列非 0 的键）
alnkit short scan --bam sampleA.bam --query-chr insertA \
    --insert-chrs insertA,insertB --merge-window-bp 10 --out-dir out_owned
#   -> INFO 统计: {'reads_kept': ..., 'drop_other_insert': ..., 'drop_ambiguous_insert': ..., ...}
```

- `drop_ambiguous_insert` 很大（几十上百）→ 有一批 read 在**多个插入之间**分不清归属 → 第 2 步试 `consensus`；
- `drop_other_insert` 很大而 `drop_ambiguous_insert` ≈ 0 → 这是**正常的**（外源面板里该位点的 read 属于别条外源），不要动策略；
- 两个都小、位点支持数也够 → **就用 `owned`**，收工。

**场景 → 推荐**：

| 你的情形 | 推荐 | 为什么 / 注意 |
|---|---|---|
| 单外源；外源与基因组内源序列同源（cDNA 拷贝、含内源片段） | **`owned`** + `--reference-fasta`（或 `--homology-bed`） | 这类干扰**不在插入之间**，任何策略都解决不了；靠 `site_homology.tsv` 的 `is_homology_site` 判。策略保持默认即可 |
| 多外源，各外源序列**互不同源** | **`owned`** | 实测这种数据上 `owned`/`consensus` 结果逐字节相同，换策略没有收益 |
| 多外源，且存在**归属不明**的 read（日志 `drop_ambiguous_insert` 大） | **`consensus`** | 用同断点投票把“模糊但上下文明确”的 read 救回来；投票依据是 `--merge-window-bp` 内的邻居，所以窗口别设 0 |
| 多外源，且**共享接头/载体/骨架**，担心共有 read 撑起假位点 | `owned` 或 `consensus` + **`--drop-shared-reads`** | 只剔“归属 query 但另在**其它**插入上有次要/XA 命中”的 read；单外源面板下它是 no-op |
| 位点支持数太少（个位数）且 `drop_ambiguous_insert ≈ 0` | **先别换策略** | 先降门槛（`--min-mapq`、`--max-mismatch`）或查 `reads_info.out.xls`；换策略也救不回来 |
| 要复跑 0.6.0 及更早的整表口径做逐字节对账 | 用 **`alnkit==0.6.0`**（或更早）跑，别再在当前版本上找开关 | `legacy` / `--legacy-profile` 已在 0.7.0 删除：它会让位点 11 → 25、真位点 SR 125 → 87，不适合正式分析 |

**四条不推荐的用法**：

1. 指望用某个开关复现 0.6.0 的整表口径（`legacy` 已删除，需要就装旧版本）；
2. 用 `--drop-shared-reads` 解决“同源位点”——它只看插入之间，管不了外源 ↔ 基因组；
3. 把 `consensus` 当“更宽松的 owned”用在单外源面板上（没有 ambiguous read 时它等于 `owned`）；
4. `--merge-window-bp 0` 时直接换 `consensus`：投票只在**完全同位**的邻居之间发生，
   实际几乎救不回任何 read（投票窗口就是 `--merge-window-bp`，推荐 10）。

**判读时永远看这两张表，不要只看支持数**：

- `site_homology.tsv` → 这个位点是不是内源同源（要 `--reference-fasta`/`--homology-bed`）；
- `Insert_Interference.result.xls` → 证据是否与插入面板共享（§6，注意它不是同源判定）。

> 0.7.0 起策略从 5 个减到 2 个：`owned-strict` 变成开关 `--drop-shared-reads`（可与
> `consensus` 组合，旧的枚举值做不到），`lenient`（无条件保留归属不明的 read）删除——它把
> “不知道属于谁”的 read 当成证据，正是这类工具最该避免的方向；`legacy`（不做归属判定、
> 全部保留）也删除——它与 `--legacy-profile` 一起，是唯一能复现旧脚本整表口径的路径，
> 代价是位点数虚高、真位点支持数被稀释。旧命令写这些值会被 argparse 直接拒绝（退出码 2）。
>
> **对账影响**：0.7.0 起只有 `Insert_Position.result.xls` 与 `reads_info.out.xls` 还能与旧脚本
> 逐字节对上（黄金回归也只比这两张）；`Insert_Interference.result.xls` 的语义已改（§6.2），
> 不再与上游产物比较。

几点容易踩的：

- **`--drop-shared-reads` 只看“其它插入”**：query 自己的次要命中不算（那是正常形态），所以
  单外源面板下它与不开这个开关完全等价。
- **两种策略可能完全等价**：同一份数据（两条互不同源的外源）上，`owned` 与 `consensus`
  产出的三张 legacy 表**逐字节相同** —— 因为这份数据里没有 `ambiguous`、也没有跨外源的共有
  read。它们的分歧**只在**“存在 ambiguous 或共有 read”时才出现。
- 策略只影响“哪些 read 算证据”，不影响干扰表的判定口径（那张表见 §6）。

## 8. 常见误读清单

| 误读 | 实际 |
|---|---|
| `is_homology_interference=True` ⇒ 这个位点是内源同源 | 不是，见 §6.3；要用 `site_homology.tsv` |
| `dup_read_count` 是证据数 | 不是，是被 markdup 标记的 read 数；默认不当证据（`--no-duplicate-filter` 才保留） |
| `--merge-window-bp` 调大能降噪 | 只是把相邻位点并成一行，会抹掉“同一座位两侧两个连接点”的结构；降噪用 `--min-sr-support` |
| `--min-clip-length` 能降噪 | 会同步削掉真位点的支持 read（同一批 read 的 clip 长度不一），只能用来抬高“证据更硬”的门槛 |
| `reads_info` 可以按表头名取列 | 不能，表头 11 项 / 数据 12 列，按列序取（§4） |
| `Insert_Position` 的数字是 1-based 还是 0-based | 两者数值相同（“断点左侧碱基的 1-based 坐标”＝“0-based 边界”），要明确基准请看 `sites.canonical.tsv` |
| minimap2 的结果和 bwa 等价 | 不等价：minimap2 默认**不写 XA 标签**，归属判定只剩 secondary，支持数会明显下降；二代链路请用 `bwa mem -Y`（**不要** `-M`） |

## 9. 怎么定 `--merge-window-bp`（实测方法，别猜）

**默认值是 0**（不合并，历史产物逐字节不变）。要不要合并、合多少，用下面三行在自己的数据上扫一遍：

```bash
for w in 0 3 10 20 50; do
  alnkit short scan --bam sampleA.markdup.bam --query-chr insertA \
      --insert-chrs insertA,insertB --merge-window-bp $w --sample-name w$w --out-dir out_w$w
done
cut -f1,2,3,4 out_w*/w*.Insert_Position.result.xls    # 比位点数与支持数
```

判读标准只有两条：**峰值位点的支持数什么时候不再涨**、**同一座位两侧的两个连接点什么时候被并掉**。

一次实测（短读长 PE，插入位点两侧是方向相反的两个 clip）的典型形状：

| 簇 | 边界点落在多少个位置 | 跨度（bp） |
|---|---|---|
| 连接点 A（支持数低） | 5 | 19 |
| 连接点 B（峰值所在） | 7 | 9 |
| A 与 B 之间 | — | 52 |

对应扫描结果：

| 窗口 | 位点数 | 峰值连接点的 SR | 备注 |
|---|---|---|---|
| 0 | 23 | 113 | 每个边界各自一行 |
| 3 | 15 | 115 | **同一连接点被拆成 7 行**（抖动 > 3） |
| 10 | 11 | 125 | 峰值连接点印成一簇；低支持那一簇还分 2 行 |
| 20 | 10 | 125 | 两簇各自成一簇，两个连接点仍分开（52 > 20） |

两条结论：

1. **抖动不是测序错误**，是比对器对 clip 端点的判定（junction 附近有微 indel / 同源时尤其明显），
   实测可达 ~20 bp —— 所以 `3` 偏小，会把同一连接点拆散、让 `--min-sr-support` 误杀真证据；
2. **聚类是锚点式**（与簇内最小坐标比、簇宽不超过窗口，不会链式无限合并），
   所以“一簇跨度 19 bp”就需要窗口 ≥19 才能收成一行；
3. 反过来别用几十到几百 bp：同一位座两侧的两个连接点通常相距几十~几百 bp，并掉就只剩一行，
   **结构信息永久丢失**。降噪请用 `--min-sr-support`，不要用窗口。

## 10. 自查：我的表跑对了吗

```bash
alnkit short scan --bam sampleA.markdup.bam --query-chr insertA \
    --insert-chrs insertA,insertB --merge-window-bp 10 --out-dir out --check-markdup
```

- 日志里 `statistics` 行给出每个阶段的计数（读到多少比对、丢掉多少、留下多少、写出多少行）——
  支持数异常时先看这里；
- `--check-markdup` 会在日志里报告 markdup 标记是否可信（未标记、标记不一致、单条被标 dup）；
- 与历史产物的兼容性由测试守住：`Insert_Position.result.xls` 与 `reads_info.out.xls` 在相同
  输入/参数下与旧脚本**逐字节一致**；canonical 表除 `sample` 列外一致；
  `Insert_Interference.result.xls` 自 0.7.0 起语义已改（§6.2），不再参与逐字节对账。
