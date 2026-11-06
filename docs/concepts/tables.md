# 表契约

> 本页由 `alnkit.core.tables` 的 `COLUMNS` / `COLUMN_DOC` 生成，**不要手改**。
> 重新生成：`python scripts/gen_table_docs.py`。

规则：

- 列序即下表顺序，写表只能通过 `alnkit.core.tables.write_tsv`；
- 坐标列名必须带 `_0based` / `_1based` 后缀；长度列（`*_bp`、`*_len`、`*_length`、`*_count`）不带；
- 旧表（legacy `.xls` / 旧 TSV）列名保持不变，由 `--columns legacy` 复现；
- 浮点固定 4 位小数；空值写 `-`。

各命令**写哪些表**、以及未注册的旧脚本表头怎么读，见
[各命令的输出表](output-tables.md)。


## `segments`

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `read_name` | — | read 名（与 BAM/FASTQ 一致） |
| 2 | `read_length` | — | read 长度（bp）；未知时留空 |
| 3 | `segment_role` | — | 段角色：genome / exog / vector / target |
| 4 | `chrom` | — | 参考序列名 |
| 5 | `chrom_length` | — | 参考序列长度（bp） |
| 6 | `ref_start_0based` | 0-based | 参考起点，0-based 半开 |
| 7 | `ref_end_0based` | 0-based | 参考终点，0-based 半开（不含） |
| 8 | `read_start_0based` | 0-based | read 轴起点，已翻到正向 read，0-based 半开 |
| 9 | `read_end_0based` | 0-based | read 轴终点，已翻到正向 read，0-based 半开（不含） |
| 10 | `strand` | — | 比对方向：+ / - |
| 11 | `mapq` | — | MAPQ |
| 12 | `ref_aln_len` | — | 参考跨度（含首尾 D/N），与简化 CIGAR 的 M 之和可不等 |
| 13 | `read_aln_len` | — | read 跨度 |
| 14 | `cigar` | — | 简化 CIGAR（两个 M 之间总长 <=5 bp 的 I/D 并入 M） |
| 15 | `kept` | — | 1/0：是否参与分类；只有 `long classify` 的产物有真实含义，`seg table` 恒为 1 |
| 16 | `drop_reason` | — | 被丢弃原因；未丢弃填 - |
| 17 | `is_secondary` | — | 是否 secondary 比对（1/0） |
| 18 | `is_supplementary` | — | 是否 supplementary 比对（1/0） |
| 19 | `nm` | — | 编辑距离（NM 标签；缺失时为空） |
| 20 | `match_columns` | — | CIGAR 里 M/=/X 的碱基总数 |
| 21 | `inserted_bases` | — | CIGAR 里 I 的碱基数 |
| 22 | `deleted_bases` | — | CIGAR 里 D 的碱基数 |

## `insert_sites`

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `chrom` | — | 参考序列名 |
| 2 | `insert_pos_1based` | 1-based | 断点左侧碱基的 1-based 坐标（= legacy Insert_Position） |
| 3 | `insert_pos_0based` | 0-based | 同一碱基的 0-based 坐标（= insert_pos_1based - 1） |
| 4 | `insert_boundary_0based` | 0-based | 断点边界，0-based（数值等于 insert_pos_1based） |
| 5 | `n_sr` | — | SR 证据数 |
| 6 | `n_pesr` | — | PESR 证据数 |
| 7 | `n_dup` | — | 该位点上被 markdup 标记的 read 数 |
| 8 | `sample` | — | 样本名（BAM 文件名首个 . 之前） |
| 9 | `query_insert` | — | 本次运行针对的插入片段（= --query-chr） |
| 10 | `nearest_other_insert` | — | 除 query 外、支持该位点 read 命中最多的插入片段（并列取名字小的）；没有其它插入片段命中时填 - |
| 11 | `other_inserts` | — | 去掉 query 与最近者之后，其余命中过的插入片段，`,` 分隔；无则 - |

## `site_homology`

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `sample` | — | 输出前缀（--sample-name） |
| 2 | `chrom` | — | 位点所在参考序列 |
| 3 | `insert_pos_1based` | 1-based | 断点左侧碱基的 1-based 坐标（= Insert_Position 列） |
| 4 | `is_homology_site` | — | yes/no：是否落在一条同源区域里（两侧各放宽 50 bp） |
| 5 | `n_regions` | — | 命中的同源区域条数（同一段可能被 splice / asm20 两个预设各命中一次，所以可能 >1） |
| 6 | `region_sources` | — | 命中的来源：insert:<插入名> 或 bed:<名字>，多个用分号连接 |
| 7 | `best_region_chrom` | — | 最可信的同源区域所在参考序列（无命中写 -） |
| 8 | `best_region_start_1based` | 1-based | 最可信同源区域的起始（1-based，闭区间） |
| 9 | `best_region_end_1based` | 1-based | 最可信同源区域的结束（1-based，闭区间） |
| 10 | `best_region_length_bp` | — | 最可信同源区域的长度（bp） |
| 11 | `best_region_identity` | — | 回比得到的一致率（0-1，按 cs 标签的真实差异算）；BED 来源写 - |
| 12 | `best_region_query_coverage` | — | 插入序列被这段区域覆盖的比例（0-1）；BED 来源写 - |

## `assembly_support`

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `class_id` | — | 类 ID（= 目录名） |
| 2 | `contig` | — | contig 名（contigs.fasta 里的记录名） |
| 3 | `contig_length_bp` | — | contig 长度（bp） |
| 4 | `n_reads_total` | — | 该类 reads 总数（FASTQ 行数 / 4） |
| 5 | `n_reads_mapped` | — | 能比回该 contig 的 read 数（MAPQ >= 10） |
| 6 | `covered_bp` | — | 被至少一条 read 覆盖的碱基数（合并区间后） |
| 7 | `covered_fraction` | — | covered_bp / contig_length_bp，0-1 浮点 |
| 8 | `mean_depth` | — | 比对碱基总数 / contig 长度 |
| 9 | `depth_p10` | — | 非零窗口深度的 10% 分位（窗口 200 bp） |
| 10 | `depth_cliffs` | — | 相邻 200 bp 窗口深度降到前一个 1/5 以下的次数（chimera 信号） |
| 11 | `ends_covered` | — | 首尾各 100 bp 中是否至少有一半被 read 覆盖（yes/no） |
| 12 | `status` | — | ok / low_coverage / chimera_suspect / unsupported_ends / no_reads |

## `ribbon_blocks`

!!! warning "尚未实现"
    目前没有命令输出这张表（`plot ribbon` 仍按旧表头写 blocks.tsv）；
    列定义先冻结在这里，等 C1 的口径统一切换。

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `query_fragment` | — | 查询片段名 |
| 2 | `query_start_0based` | 0-based | 查询片段自身坐标起点，0-based 半开 |
| 3 | `query_end_0based` | 0-based | 查询片段自身坐标终点，0-based 半开（不含） |
| 4 | `strand` | — | 查询相对参考的方向：+ / - |
| 5 | `reference` | — | 参考序列名 |
| 6 | `ref_start_0based` | 0-based | 参考起点，0-based 半开 |
| 7 | `ref_end_0based` | 0-based | 参考终点，0-based 半开（不含） |
| 8 | `block_bp` | — | 块长度（= ref_end - ref_start） |
| 9 | `identity` | — | 一致率，0-1 浮点，固定 4 位小数 |

## `insertion_clip_reads`

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `sample` | — | 样本名（BAM 文件名首个 . 之前的部分） |
| 2 | `chrom` | — | 参考序列名 |
| 3 | `breakpoint_0based` | 0-based | 断点边界，0-based（数值等于断点左侧碱基的 1-based 坐标） |
| 4 | `clip_side` | — | clip 在比对块的哪一侧：left / right |
| 5 | `read_name` | — | read 名（与 BAM 一致） |
| 6 | `read_length` | — | read 长度（含 soft clip 的碱基，bp） |
| 7 | `strand` | — | 比对方向：+ / - |
| 8 | `mapq` | — | MAPQ |
| 9 | `ref_start_0based` | 0-based | 比对块起点，0-based 半开 |
| 10 | `ref_end_0based` | 0-based | 比对块终点，0-based 半开（不含） |
| 11 | `clip_length_bp` | — | soft clip 长度（bp） |
| 12 | `clip_sequence` | — | clip 碱基，与参考同向；left = 比对块左侧那一段（反向链 read 上是 read 的 3' 端） |

## `insertion_clip_consensus`

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `sample` | — | 样本名（BAM 文件名首个 . 之前的部分） |
| 2 | `chrom` | — | 参考序列名 |
| 3 | `breakpoint_0based` | 0-based | 断点边界，0-based |
| 4 | `breakpoint_1based` | 1-based | 同一个断点，按“左侧碱基的 1-based 坐标”写（数值与 0based 列相同） |
| 5 | `clip_side` | — | left / right |
| 6 | `n_reads` | — | 该侧参与投票的 clip read 数（候选数，不代表每个碱基都投了票） |
| 7 | `consensus_bp` | — | 一致性序列长度（bp） |
| 8 | `consensus_fraction` | — | 一致率 0-1（每列多数碱基数之和 / 全部条数之和） |
| 9 | `stop_reason` | — | 中断原因：all columns / coverage<2 (only N reads) / 2-read overlap mismatch / no clip |
| 10 | `consensus_sequence` | — | 一致性序列（该侧无序列写 -）；left 侧以接头端为第 0 列、向断点外侧延伸（参考方向） |

## `coverage_region`

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `sample` | — | 样本名（BAM 文件名首个 . 之前的部分） |
| 2 | `chrom` | — | 参考序列名 |
| 3 | `region_start_1based` | 1-based | 统计区间起点，1-based 闭区间 |
| 4 | `region_end_1based` | 1-based | 统计区间终点，1-based 闭区间 |
| 5 | `region_length_bp` | — | 区间长度（bp） |
| 6 | `n_bins` | — | 分箱个数（= coverage_bins 的行数） |
| 7 | `n_reads_total` | — | fetch 到的比对记录数（含被过滤掉的） |
| 8 | `n_reads_selected` | — | 通过 flag / mapQ / 长度过滤的 read 数 |
| 9 | `n_reads_filtered` | — | 被过滤掉的 read 数（= 总数 - 选中数） |
| 10 | `mean_mapq` | — | 选中 read 的平均比对质量 |
| 11 | `mean_baseq` | — | 被覆盖位置上的平均碱基质量（只出现在终端摘要，图中不显示） |
| 12 | `covered_bases` | — | 被覆盖的碱基数（深度 >= min_depth 的位置数） |
| 13 | `percent_covered` | — | 覆盖广度百分比 0-100（= covered_bases / region_length_bp） |
| 14 | `mean_depth` | — | 全区间平均深度（分子只累加被覆盖位置，分母是区间长度） |

## `coverage_bins`

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `bin_index` | — | 第几个 bin，从 1 开始 |
| 2 | `bin_start_1based` | 1-based | bin 起点，1-based 闭区间 |
| 3 | `bin_end_1based` | 1-based | bin 终点，1-based 闭区间（相邻 bin 首尾相接） |
| 4 | `bin_width_bp` | — | bin 宽度（bp）；最后一个 bin 吸收除不尽的余数 |
| 5 | `covered_bases` | — | bin 内被覆盖的碱基数 |
| 6 | `depth_sum` | — | bin 内逐碱基深度之和（除以宽度就是该 bin 的平均深度） |

## `coverage_bases`

| # | 列名 | 基准 | 说明 |
|---|---|---|---|
| 1 | `pos_1based` | 1-based | 碱基坐标，1-based |
| 2 | `depth` | — | 该位置的深度（已按 --max-depth 封顶） |
