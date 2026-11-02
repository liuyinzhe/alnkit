# 表契约

> 本页由 `alnkit.core.tables` 的 `COLUMNS` / `COLUMN_DOC` 生成，**不要手改**。
> 重新生成：`python scripts/gen_table_docs.py`。

规则：

- 列序即下表顺序，写表只能通过 `alnkit.core.tables.write_tsv`；
- 坐标列名必须带 `_0based` / `_1based` 后缀；长度列（`*_bp`、`*_len`、`*_length`、`*_count`）不带；
- 旧表（legacy `.xls` / 旧 TSV）列名保持不变，由 `--columns legacy` 复现；
- 浮点固定 4 位小数；空值写 `-`。


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
