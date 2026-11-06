# `short scan` 的输出表怎么读

本页解释 `alnkit short scan` 产出的每一个文件、每一列的含义，以及**最容易读错的两张表**
（`Insert_Interference` 与 `site_homology`）。示例一律用占位名（`sampleA`、`chr1`、
`insertA`/`insertB`），和你自己的产物对照时请换成实际名字。

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

列定义见[表契约](tables.md)的 `insert_sites`。与主表**按 `(chrom, insert_pos_1based)` 一一对应**
（同样的位点、同样的合并结果），恒等式：

```
insert_pos_1based - 1 == insert_pos_0based
insert_boundary_0based == insert_pos_1based
```

`nearest_other_insert` / `other_inserts` 来自干扰判定，但**明确把 query 自己排除在外**：
“query 出现”只说明这条 read 也命中了本次的 query，不算竞争者（见 §6.3）。

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
| `support_reads_owned_other` | 归属**其它**插入序列的 read 数（`legacy` 策略下才会非 0） |
| `support_reads_consensus_assigned` | 原本不明确、被同断点投票推定归属的 read 数（`consensus` 策略） |
| `support_reads_ambiguous` | 归属不明确的 read 数 |
| `support_reads_no_insert_hit` | 完全没有插入命中的 read 数（`legacy` 策略下才会非 0） |
| `ambiguous_read_names` | 上面“不明确”那批 read 的名字，`;` 分隔 |

### 6.2 判定链（为什么源头里会出现 query 自己）

1. 每条比对贡献“插入命中”：**自身命中**（source=`alignment`/`secondary`）＋ **XA 标签**里的
   备选命中（source=`xa`）；
2. “同源信号”只取 **secondary 与 XA** 命中（`alignment` 不算）。注意：**supplementary 属于
   `alignment`**，所以“断点 read 的另一段落在插入序列上”（正常形态）不会触发干扰；
3. `--interference-scope` 决定怎么汇总：
   - `all`（默认，历史口径）：secondary/XA 命中**任何**插入序列（**含 query**）就算源头；
     归属不明确时把该 read 的全部候选插入序列都算源头；
   - `other`：先把 **query 从命中集合与候选集合里剔除**，再按同样规则汇总。

所以 `all` 下的 `interference_source_support_reads` 里出现 `insertA:3`（= query 自己的名字）
**是设计出来的行为，不是串列**：它表示“这个位点有 3 条支持 read 在插入面板上还有次要/XA 命中”。

### 6.3 它的价值：**不是“基因组同源性”**

这个表名里的 homology 指的是**“证据与插入面板共享”**（同源 read 在多个插入序列间横跳），
**不是**“外源序列与基因组同源”。两者后果完全不同：

| 场景 | 该看什么 | 为什么 |
|---|---|---|
| 多外源面板（`insertA`/`insertB`/`insertC`），某些位点的 read 同时命中别的外源 | 本表 + `--interference-scope other` | 这就是它擅长的：挑出“跨外源共享”的位点 |
| 外源序列与**基因组内源**序列同源（例如 query 是某内源基因的 cDNA 拷贝） | `<sample>.site_homology.tsv` | 判据必须是“把插入序列回比基因组”，而不是“read 是否命中插入面板” |

在“单外源 + 该外源与基因组同源”的数据上，本表的两个口径**都会失去区分力**：

- `all`：真位点也会被标 `True`，源头是 query 自己——因为断点 read 在插入面板上有次要命中是
  正常形态；同源位点同样 `True`，两者**长得一样**；
- `other`：query 被剔除后，单外源面板没有“别的插入”可比，于是**所有位点都变成 `False`**。

结论：**不要用 `is_homology_interference` 做同源过滤**。要排除“外源 ↔ 内源同源”造成的假位点，
请用 `--reference-fasta`（插入序列回比位点周围窗口）或 `--homology-bed`（已知同源区域），
看 `site_homology.tsv` 的 `is_homology_site`。

### 6.3.1 仍然有用的读法：看**比例**，别看布尔值

布尔列（`is_homology_interference`）没有区分力，但
`interference_source_support_reads / support_reads_total` 这个**比例**有：

| 位点类型 | 比例 | 解读 |
|---|---|---|
| 内源同源位点（read 真正是“两处都能比上”） | ≈ 100% | 几乎所有支持 read 都有插入面板上的次要命中 |
| 真整合位点 | 几个百分点 | 绝大多数 read 是干净的 chr↔插入 split read，只有个别 read 带次要命中 |
| 噪声位点（单条 read 落错） | 0 | 那一条 read 没有面板上的次要命中 |

在“单外源 + 外源与基因组内源同源”的数据上实测就是这个形状（100% / 2–3% / 0）。
两条使用前提：**同一台 aligner**（bwa 写 XA、minimap2 不写，跨比对器不可比）、
**同一次运行**（比例不是生物学常数，只是这批 read 的歧义度）。

### 6.4 同一个名字，两张表里含义不同

- 本表 `interference_source_*` 里出现 `insertA`（query）＝ **共享证据**的信号（很可能是噪声，
  见 §6.3）；
- `site_homology.tsv` 的 `region_sources` 里出现 `insert:insertA` ＝ **同源证据**：
  这段参考序列与插入序列一致（内源拷贝就在那里）。这是**正常且期望**的标记，不是错误。

两张表都“包含 query 的名字”，但一个是“嫌疑”，一个是“结论”，不要互相印证。

### 6.5 其他注意

- 表里**没有**记录本次用的 `--interference-scope`（历史格式不带配置列），看日志首行的参数或
  命令历史确认；
- `interference_source_count` 数的是**源头个数**，不是 read 数；read 数在
  `interference_source_support_reads` 里；
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
| `owned`（默认） | 归属 **query** 的 read | 归属其它插入、`ambiguous`、`no_insert_hit` | 多外源面板的常规选择：只信“确定属于本外源”的证据 |
| `owned-strict` | 同 `owned`，但再剔除“归属 query、却在**其它**插入上还有次要/XA 命中”的 read | 上面那些 + 共有 read | 外源之间高度相似（共享接头/载体）时，避免共有 read 把特异位点撑起来 |
| `lenient` | `owned` + `ambiguous` | 只丢归属其它插入 / 无插入命中 | 判定过严导致支持数太低时的兜底；**不投票，最宽松** |
| `consensus` | `owned` + 用**同断点投票**推定 `ambiguous` 的归属（附近 ±`--merge-window-bp` 内归属明确的 read 一致指向 query 才保留） | 投票为空或指向不一致的 ambiguous | 同源干扰严重、又不想丢证据时；比 `lenient` 保守（**推断不是测量**） |
| `legacy` | **全部保留**（不做任何归属判定），并自动关闭候选 read 闸门 | 只有质量过滤（MAPQ / 错配 / clip 段数） | 复现旧脚本产物做逐字节对账；`--legacy-profile` 会隐含它 |

几点容易踩的：

- **`legacy` 不是“更宽松的 owned”**：它连“没命中任何插入序列的 read”都算证据，位点数会明显变多、
  真位点的支持数反而可能变小（证据被摊到更多位点上）。实测一份 PE 数据：
  `owned` 11 个位点、真位点 SR=125；`legacy` 25 个位点、真位点 SR=87。
- **`owned-strict` 只看“其它插入”**：query 自己的次要命中不算（那是正常形态），所以单外源面板下
  它与 `owned` 完全等价。
- **四种非 legacy 策略可能完全一样**：同一份数据（两条互不同源的外源）上，
  `owned` / `owned-strict` / `lenient` / `consensus` 产出的三张 legacy 表**逐字节相同**，
  只有 `legacy` 不同 —— 因为这份数据里没有 `ambiguous`、也没有跨外源的共有 read。
  它们的分歧**只在**“存在 ambiguous 或共有 read”时才出现。
- 策略只影响“哪些 read 算证据”，不影响干扰表的判定口径（那张表见 §6；它的 `interference_scope`
  是另一个旋钮）。

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
- 与历史产物的兼容性由测试守住：legacy 三表在相同输入/参数下**逐字节一致**，
  canonical 表除 `sample` 列外一致。
