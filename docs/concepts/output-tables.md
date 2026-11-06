# 各命令的输出表：谁写什么、每列什么意思

本页是**表格产物的全景索引**：哪个命令写哪些文件、每列的基准与含义。三层分工：

| 文档 | 内容 |
|---|---|
| [表契约](tables.md) | **注册表**的规范列定义（由 `core/tables.py` 生成，改代码必须同步） |
| 本页 | 命令 → 产物清单；**未注册**（旧脚本表头 / 手写表头）的表逐列解释 |
| [`short scan` 输出表](short-scan-output.md) | `short scan` 的逐列**读法**（含干扰表与同源表的坑） |

通用约定（所有表都一样）：

- 坐标列名带 `_0based` / `_1based` 后缀；旧脚本表头（`ref_start` / `Insert_Position` 等）沿用历史
  口径，**不要按名字猜基准**——§1 的“基准”一列写明了；
- 浮点固定 4 位小数；空值写 `-`；多值用 `,` 或 `;` 分隔（表内注明）；
- 行序由写入方显式排序（表内注明排序键），同一输入两次运行逐字节一致。

## 1. 总览

| 命令 | 表格产物 | 列定义 / 基准 |
|---|---|---|
| `short scan` | `<sample>.Insert_Position.result.xls`<br>`<sample>.reads_info.out.xls`<br>`<sample>.Insert_Interference.result.xls`<br>`<sample>.sites.canonical.tsv`<br>`<sample>.site_homology.tsv`（可选） | 旧脚本表头 → [scan 输出表](short-scan-output.md) §3–§6<br>注册表 `insert_sites` / `site_homology` → [tables.md](tables.md) |
| `plot insertion` | `<name>.clips.tsv`、`<name>.consensus.tsv` | 注册表 `insertion_clip_reads` / `insertion_clip_consensus` |
| `plot coverage` | `<name>.region.tsv`、`<name>.bins.tsv`、`<name>.bases.tsv`（`--bases-table`） | 注册表 `coverage_region` / `coverage_bins` / `coverage_bases` |
| `plot tracks` | `<name>.segments.tsv` | 注册表 `segments` |
| `plot ribbon` | `<name>.blocks.tsv`、`<name>.placements.tsv`、`<name>.regions.tsv` | §4.1（前两张在文件里带 `#` 列说明） |
| `seg table` | `<out>`（规范列，或 `--columns legacy` 历史 14 列） | 注册表 `segments` / §4.3 |
| `long classify` | `alignment_segments.tsv`、`read_classification.tsv`、`junction_events.tsv`、`junction_sites.tsv`、`locus_regions.tsv`、`classification_summary.tsv` + `related_reads.fastq` | §3.1 |
| `long graph` | `overlap_edges.tsv`、`overlap_components.tsv` | §3.2 |
| `long breakpoints` | `read_class_assignment.tsv`、`breakpoint_sites.tsv`、`breakpoint_classes.tsv`、`breakpoint_shared_reads.tsv`、`segment_ends.tsv`、`class_read_pools.tsv` + `classes/` | §3.3 |
| `long assemble` | `<class>/overlap_extension/contigs.fasta`、`assembly_stats.tsv`、`overlap_assembly_summary.tsv`、`assembly_support.tsv` | §3.4 + 注册表 `assembly_support` |
| `long locate` | `queries.tsv`、`class_assembly_positions.tsv`、`class_position_summary.tsv`、`contig_map.tsv` + 丝带图那三张表 | §3.5 + §4.1 |

> 长读长链路里除 `assembly_support.tsv` 外，全部是**旧脚本表头**（保留逐字节兼容），列名不是
> 带基准后缀的规范名；§3 逐列注明了基准。

## 2. 短读长

`short scan` 的 5 个产物的逐列读法在 [scan 输出表](short-scan-output.md)：主表 §3、
逐比对明细 §4、规范表 §5、同源表 §5.1、干扰表 §6。这里只给一句话：

| 文件 | 一句话 |
|---|---|
| `Insert_Position.result.xls` | 位点 + SR/PESR/dup 支持数与 read 名单 |
| `reads_info.out.xls` | 参与统计的每条比对的明细（**表头 11 项 / 数据 12 列**，按列序解析） |
| `Insert_Interference.result.xls` | 该位点的证据是否与**别的**插入序列共享（判不了内源同源） |
| `sites.canonical.tsv` | 同一批位点的机器接口（`insert_pos_1based/_0based/boundary` 三列并列） |
| `site_homology.tsv` | 位点是否落在“与插入序列同源的参考区域”里（要 `--reference-fasta`/`--homology-bed`） |

`plot insertion` / `plot coverage` 产出的都是**注册表**，列定义见 [tables.md](tables.md)；
插入图的画法与表格的关系见 README「二代出图」一节。

## 3. 长读长

### 3.1 `long classify`（6 张表）

**`read_classification.tsv`** — 一行一条目标 read（只有“有合格外源比对”的 read 会进来）

| 列 | 含义 |
|---|---|
| `read_name` / `read_length` | read 名与长度（bp，FASTQ 原始长度） |
| `read_category` | `junction`（同时有外源与基因组证据）/ `vector_region`（只有外源证据）；`genome_only`/`other` 不写进本表 |
| `n_exog_segments` / `exog_names` / `exog_aln_bp` | 外源侧比对段数、命中的外源名（`,` 分隔）、这些段在 read 轴上的比对碱基总数 |
| `n_genome_segments` / `genome_aln_bp` / `genome_chroms` | 基因组侧段数、比对碱基数（**read 轴区间并集**，不是参考跨度）、涉及的参考名 |
| `n_genome_segments_dropped` / `genome_aln_bp_dropped` | 被丢弃的基因组段数与其碱基数（例如长度不足门槛） |
| `max_mapq` | 该 read 全部比对段的最高 MAPQ |
| `notes` | 备注（判定过程中的提示，`;` 分隔；无则 `-`） |

**`alignment_segments.tsv`** — 同上那批 read 的**逐段**比对（旧表头）

| 列 | 含义 |
|---|---|
| `read_name` / `segment_role` | read 名；段角色 `genome` / `exog`（`vector`/`target` 也出现在别的链路里） |
| `chrom` / `chrom_length` | 比对到的参考名与长度 |
| `ref_start` / `ref_end` | 参考区间，**0-based 半开**（旧表头，无后缀） |
| `read_start` / `read_end` | read 轴区间，已统一到**正向 read**（反向记录做过镜像），0-based 半开 |
| `strand` / `mapq` | 方向 `+`/`-`；MAPQ |
| `ref_aln_len` / `read_aln_len` | 参考跨度 / read 跨度（bp） |
| `simplified_cigar` | 简化 CIGAR：两个 M 之间总长 ≤5 bp 的 I/D 并进 M；**不是** BAM 原始 CIGAR（原始 CIGAR 不保存在表里） |
| `kept` / `drop_reason` | `1`=参与分类、`0`=不参与；被丢弃的原因（未丢弃填 `-`） |

**`junction_events.tsv`** — 每个“基因组末端紧邻外源”的接头事件一行

| 列 | 含义 |
|---|---|
| `read_name` / `chrom` / `genome_side` | read 名；参考名；基因组在边界的哪一侧（`genome_left` = 基因组在左、插入从边界开始） |
| `boundary` | 接头边界，**0-based**（`genome_side=genome_left` 时数值 = 最后一个基因组碱基的 1-based 坐标；`genome_right` 时 = 插入后第一个基因组碱基的 0-based 坐标） |
| `region_id` / `site_id` | 该事件归属的区域与位点（见下面两张表） |
| `vector_name` | 外源侧的序列名 |
| `gap_bp` | 基因组段与外源段在 **read 轴**上的间隔（0 = 紧邻或重叠；负数不会出现） |
| `genome_read_start` / `genome_read_end` | 基因组段在 read 轴上的区间（0-based 半开） |
| `vector_read_start` / `vector_read_end` | 外源段在 read 轴上的区间 |
| `genome_mapq` / `vector_mapq` | 两段的 MAPQ |
| `genome_identity` / `vector_identity` | 两段的一致率（0-1，按错配/比对长度算） |

**`junction_sites.tsv`** — 把事件按“同染色体、同 `genome_side`、边界差 ≤ `--site-gap-bp`（默认 20 bp）”聚成位点

| 列 | 含义 |
|---|---|
| `site_id` / `region_id` | 位点编号与所属区域 |
| `chrom` / `genome_side` | 参考名；基因组侧（`genome_left` = 基因组在边界左侧） |
| `boundary_min` / `boundary_max` / `boundary_median` | 该位点内事件边界的最小/最大/中位数（**0-based**） |
| `n_events` / `n_reads` | 事件数 / 支持 read 数 |
| `exog_names` / `read_names` | 命中的外源名与支持 read 名（`,` / `;` 分隔） |

**`locus_regions.tsv`** — 位点再按“同染色体、边界相距 ≤ `--region-gap-bp`（默认 200 kb）”聚成区域

| 列 | 含义 |
|---|---|
| `region_id` / `region_label` | 区域编号与可读标签（`chr:start-end` 形态） |
| `chrom` / `start` / `end` | 参考名与区域边界（**0-based**） |
| `n_sites` / `n_reads` | 区域内位点数 / 支持 read 数 |
| `genome_sides` / `site_ids` | 涉及的基因组侧（集合）/ 区域内的位点编号 |
| `exog_names` | 区域内命中的外源名 |

**`classification_summary.tsv`** — 一行汇总

| 列 | 含义 |
|---|---|
| `total_reads` | FASTQ 总 reads（一次扫描统计） |
| `target_reads` | 有合格外源比对、进入分类的 read 数 |
| `junction` / `vector_region` | 两类目标 read 的数量 |
| `non_target_reads` | `total_reads - target_reads`（这些 read 不加载比对段、不输出逐条记录） |
| `target_fraction` | `target_reads / total_reads`（0-1） |

`related_reads.fastq` 是上表那批 `junction` + `vector_region` read 的原始序列，供下游组装用。

### 3.2 `long graph`（2 张表）

**`overlap_edges.tsv`** — read-vs-read overlap 边（每对 read 只保留 `block_bp` 最大的一条）

| 列 | 含义 |
|---|---|
| `read_a` / `read_b` | 参与 overlap 的两条 read（顺序固定为名字字典序，保证可复现） |
| `block_bp` | overlap 块长度（bp）；门槛 `--min-overlap` |
| `identity` | 一致率（0-1）；门槛 `--min-identity` |
| `strand` | 两条 read 的相对方向（`+`/`-`） |
| `read_a_start` / `read_a_end` | 该块在 `read_a` 自身坐标上的区间（**0-based 半开**） |
| `read_b_start` / `read_b_end` | 同上，`read_b` 侧 |

**`overlap_components.tsv`** — 全局连通分量（只用于观察 overlap 结构）

| 列 | 含义 |
|---|---|
| `component_id` | 连通分量编号（按成员数降序、同数按首成员名字） |
| `n_reads` | 分量内 read 数 |
| `categories` | 分量内出现过的 read 类别（`,` 分隔） |
| `n_vector_reads` | 其中 `vector_region` 类 read 数 |
| `read_names` | 全部成员名（`;` 分隔） |

> 方案 A 的“类”不由这张表决定，而是由断点共享 reads 决定（见 §3.3）。

### 3.3 `long breakpoints`（6 张表 + `classes/`）

**`segment_ends.tsv`** — 每个基因组片段的左右末端一行（断点候选的底表）

| 列 | 含义 |
|---|---|
| `read_name` / `chrom` / `strand` | read 名；参考名；该片段方向 |
| `segment_start` / `segment_end` | 片段在参考上的区间（**0-based 半开**，旧表头） |
| `end_side` | 这是片段的哪一端：`left` / `right` |
| `end_coordinate` | 该末端的参考坐标（0-based） |
| `genome_side` | 该末端“朝外”的基因组侧（`genome_left` = 基因组在左，即插入点在右） |
| `adjacent_type` / `adjacent_name` | 在 **read 轴**上紧邻这个末端的片段类型与序列名（`vector` / `genome` / `none`） |
| `adjacent_coordinate` | 那个紧邻片段落在什么坐标（按 `adjacent_type` 解释：外源侧用外源坐标，基因组侧用参考坐标） |
| `junction_type` | `vector` = 基因组末端紧邻外源；`genome` = 紧邻另一段基因组；`none` = 没有紧邻段 |
| `junction_id` | 该末端对应的断点编号（`S##` 为 vector 断点；`G##a`/`G##b` 为同一个 genome–genome 断点的两个末端） |
| `gap_bp` | 末端与紧邻片段在 read 轴上的间隔（0 = 紧邻/重叠） |
| `window_50bp` | 该末端 50 bp 锚点窗口的判定信息（`ok` = 窗口完整可用） |

**`breakpoint_sites.tsv`** — 按断点聚合（`S##` / `G##a` / `G##b`）

| 列 | 含义 |
|---|---|
| `site_id` / `region_id` | 断点编号、所属区域（`region_id` 供 `long locate` 用） |
| `chrom` / `genome_side` | 参考名、基因组侧 |
| `boundary_min` / `boundary_max` / `boundary_median` | 该断点内边界的最小/最大/中位（**0-based**） |
| `n_events` / `n_reads` | 支持事件数 / 支持 read 数 |
| `exog_names` / `read_names` | 涉及的外源名 / 支持 read 名 |
| `junction_type` | `vector`（基因组紧邻外源）/ `genome`（基因组紧邻基因组） |
| `junction_pair` | 同一个 genome–genome 断点的另一端（`G##a` ↔ `G##b`）；vector 断点填 `-` |

**`breakpoint_shared_reads.tsv`** — 两两位点之间共享多少 reads（合并成类的依据）

| 列 | 含义 |
|---|---|
| `site_a` / `site_b` | 两个断点编号 |
| `n_shared_reads` | 同时支持两者的 read 数；`--min-shared`（默认 2）是合并门槛 |
| `support_a` / `support_b` | 各自的支持 read 数 |
| `cross_overlap_pairs` | A 独有 reads 与 B 独有 reads 之间的 overlap 对数（支持度参考，**不参与**合并） |
| `merge` | `1` = `n_shared_reads >= --min-shared`，即合并成同一类 |

**`breakpoint_classes.tsv`** — 方案 A 的类（共享 reads 够多就合并）

| 列 | 含义 |
|---|---|
| `class_id` | 类编号（= `classes/` 下的目录名） |
| `region_id` / `region_label` | 所属区域 |
| `chrom` / `start` / `end` | 类的参考范围（**0-based**） |
| `n_sites` / `n_reads` | 类内断点数 / junction read 数 |
| `genome_sides` / `site_ids` | 涉及的基因组侧 / 断点编号 |
| `exog_names` / `junction_types` | 涉及的外源名 / 断点类型集合 |
| `read_names` | 类内 junction reads（并集，`;` 分隔） |

**`read_class_assignment.tsv`** — read → 类 / 断点的对照

| 列 | 含义 |
|---|---|
| `read_name` / `read_category` | read 名、类别（同 §3.1 的 `read_category`） |
| `class_ids` | 该 read 支持的断点所属的类；载体 reads 没有断点，填 `-` |
| `site_ids` | 该 read 支持的断点（`S##` / `G##a` / `G##b`） |

**`class_read_pools.tsv`** — 每类的组装输入池

| 列 | 含义 |
|---|---|
| `class_id` | 类编号 |
| `n_junction_reads` / `n_vector_reads` / `n_total_reads` | 池内 junction reads / 载体区域 reads / 合计 |
| `pool_read_bp` | 池内所有 read 的长度之和（**不是**组装长度） |
| `n_vector_with_overlap_support` | 其中能在 `overlap_edges.tsv` 里找到支持边的载体 read 数 |
| `reads_fastq` | 该池写出的 FASTQ 路径 |

`classes/<class_id>/` 下另有该类的 `reads.fastq`、`segments.tsv` 等中间文件。

### 3.4 `long assemble`（2 张表 + 序列 + 支撑度表）

| 文件 | 列 |
|---|---|
| `assembly_stats.tsv` / `overlap_assembly_summary.tsv` | `class_id`、`n_contigs`、`assembly_bp`（该类的总组装碱基）、`longest_contig_bp`——两个文件列相同，前者是“按类”的当前结果，后者是历史脚本的汇总口径，**保留两份是为了逐字节对账** |
| `<class>/overlap_extension/contigs.fasta` | 组装序列（非表格） |
| `assembly_support.tsv` | **注册表** `assembly_support`：把类内 reads 用 minimap2 打回 contig 的回验证据（`covered_fraction` / `mean_depth` / `depth_cliffs` / `status` 等），列定义见 [tables.md](tables.md) |

### 3.5 `long locate`（4 张表 + 丝带图三张表）

**`queries.tsv`** — 每个类实际用的查询片段

| 列 | 含义 |
|---|---|
| `query_name` / `query_type` | 查询名；`flank`（断点旁 50 bp 锚点）或 `exog`（外源序列） |
| `source_name` / `source_region` | 取序列的来源名与区间（外源查询没有参考坐标） |
| `chrom` / `start` / `end` / `length` | 参考名与区间（**0-based 半开**）；外源查询此三列填 `-`/长度 |

**`class_assembly_positions.tsv`** — 每个查询片段落在组装 contig 的哪里、又还原到参考的哪里

| 列 | 含义 |
|---|---|
| `class_id` / `assembly_contig` | 类编号、所用的 contig |
| `query_name` / `query_type` | 查询片段名与类型（`flank` / `exog`） |
| `query_start` / `query_end` / `query_strand` | 查询片段自身坐标（0-based 半开）与相对 contig 的方向 |
| `contig_start` / `contig_end` | 在 contig 上的区间（0-based 半开） |
| `identity` | 一致率（0-1） |
| `reference_chrom` / `reference_start` / `reference_end` | **还原出的参考坐标**（0-based 半开）。`query_type=flank` 的行才是“从 contig 反推的基因组坐标”，`exog` 行填 `-` |
| `query_source` | 该查询片段的**提取坐标**（与 `reference_*` 不是一回事） |
| `distance_to_breakpoint_bp` | 还原区间到侧翼对应接头边界的距离（0 = 跨过边界） |

**`class_position_summary.tsv`** — 每类一行汇总

| 列 | 含义 |
|---|---|
| `class_id` / `n_contigs` | 类编号、contig 数 |
| `n_flank_queries` / `n_flank_queries_hit` | 侧翼 50 bp 锚点数 / 成功定位数 |
| `n_exog_queries` / `n_exog_queries_hit` | 外源查询数 / 成功落点数 |
| `flank_positions` | 侧翼锚点的最佳位置（多值时 `[总数]` 写在前面，最多列 5 处） |
| `exog_positions` | 外源序列的落点（同上，最多 5 处） |
| `figures` | 本次生成的丝带图文件名 |

**`contig_map.tsv`** — 类 → contig 的索引

| 列 | 含义 |
|---|---|
| `class_id` / `contig` | 类编号、contig 名 |
| `contig_length_bp` | contig 长度（bp） |
| `assembly_fasta` | 该 contig 所在的 fasta 路径 |

## 4. 绘图

### 4.1 `plot ribbon`（3 张表）

| 文件 | 列 |
|---|---|
| `<name>.blocks.tsv` | 图里用到的每个映射块：`query_fragment`、`query_start`、`query_end`、`strand`、`reference`、`ref_start`、`ref_end`、`block_bp`、`identity`。**旧表头**（列名不带 `_0based`，但坐标是 0-based 半开；`query_*` 是查询片段自身坐标）。注册表里的规范名是 `ribbon_blocks`（`query_start_0based` 等），迁移期有意不切换 |
| `<name>.placements.tsv`（`--placement-table`） | 一行一个落点：`query_name`、`query_length_bp`、`placement_index`、`n_placements`、`reference_name`、`reference_length_bp`、`ref_start`、`ref_end`、`ref_length_bp`、`query_start`、`query_end`、`query_length_in_hit_bp`、`strand`、`identity_percent`、`covers_full_query`。**文件头自带 `#` 逐列说明** |
| `<name>.regions.tsv`（`--region-table`） | 一行一个同源区域：`block_id`、`ref_name`、`ref_length_bp`、`ref_start`、`ref_end`、`ref_region_bp`、`gap_from_previous_ref_region_bp`、`ref_coverage`、`query_order`、`query_name`、`query_length_bp`、`query_start`、`query_end`、`query_region_bp`、`query_coverage`、`strand`、`identity_percent`、`ref_region`、`query_region`。**文件头自带 `#` 逐列说明** |

### 4.2 `plot tracks` / `seg table`

| 命令 | 文件 | 列 |
|---|---|---|
| `plot tracks` | `<name>.segments.tsv` | **注册表** `segments`（规范列，见 [tables.md](tables.md)） |
| `seg table` | `<out>`（`--out` 指定） | 默认**注册表** `segments`；`--columns legacy` 时写历史 14 列：<br>`read_name, read_pair, cigarstring, is_supplementary, chrom, M_ref_pos_start, M_ref_pos_end, M_read_pos_start, M_read_pos_end, M_var_len, stand, mapping_quality, ref_aln_len, read_aln_len`（旧表头，无基准后缀；`M_*` 为证据 M 块在参考/read 轴上的区间） |

## 5. 自查：表是不是被改坏了

```bash
python scripts/gen_table_docs.py --check      # 注册表的列定义与文档一致
python -m pytest -o addopts="" tests/test_tables.py   # 坐标列后缀规则
```

- 注册表的列增删必须同时改 `core/tables.py` 的 `COLUMNS` + `COLUMN_DOC`，再重生成
  `docs/concepts/tables.md`（CI 会跑 `--check`）；
- 旧脚本表头（§3 全部、`short scan` 的前三张、ribbon 的 `blocks.tsv`）**不要**改名——它们承担
  与历史产物逐字节对账的职责；要规范名请用对应的注册表。
