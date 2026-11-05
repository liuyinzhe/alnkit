# 快速开始

## 1. 自检

```bash
alnkit doctor      # 环境与外部命令
alnkit selftest    # 内置断言，不需要输入数据
```

`selftest` 全绿说明坐标换算与表契约正常。

## 2. 检查一条 read 上都有什么

```bash
# 直接从 BAM 画
alnkit plot tracks --bam reads.bam --read "readA" \
    --exog insertA,insertB --highlight-ref chr1:1000001 --outdir out

# 也可以先导出段表，手改后再画
alnkit seg table --bam reads.bam --out out/segments.tsv
alnkit plot tracks --segments-tsv out/segments.tsv --outdir out
```

产物：`out/read_tracks.png|svg|pdf` + `out/read_tracks.segments.tsv`。

## 3. 看一段序列落在哪些组装的哪里

```bash
alnkit plot ribbon --reference "insert=exog.fa:insertA" \
    --query "C01=assembly.fa:oe_contig1" --query "C02=assembly.fa:oe_contig2" \
    --outdir out --name exog_vs_classes
```

产物：`out/exog_vs_classes.png|pdf` + `out/exog_vs_classes.blocks.tsv`。

## 4. 二代插入断点证据

```bash
alnkit short scan --bam-dir ./bams --query-chr target_ref \
    --insert-chrs decoy_ref,homolog_ref,target_ref \
    --merge-window-bp 10 --out-dir out
```

产物：三张 legacy `.xls`（列名与旧脚本一致，供人读/Excel）
+ `<sample>.sites.canonical.tsv`（机器接口：`insert_pos_1based` / `insert_pos_0based` /
`insert_boundary_0based` 三列并列，恒等式由测试守住）。

### 参数怎么选（推荐做法）

| 参数 | 建议 | 为什么 |
|---|---|---|
| `--merge-window-bp` | **保持小值（10–50）** | 它只该吸收 **clip 端点抖动（几个 bp）**。放大到几百 bp 会把同一座位两侧的两个连接点（例如 73,520,999 与 73,521,229，相距 230 bp、clip 方向相反）并成一行，丢掉"这是一个双连接点座位"的结构信息 |
| `--min-sr-support` | 测序深度够时用 **2** | 降噪应该靠"支持数"，不是靠合并窗口。SR=1 的位点多为旁系同源/比对噪声 |
| `--min-clip-length` | **一般保持 0** | 它会把真位点的支持 read 一起削掉（同一批 read 的 clip 可能较短），不适合降噪 |
| `--reference-fasta` 或 `--homology-bed` | 有同源外源时**务必给** | 外源与基因组内源序列同源时（如 GLP1R-FLAG vs 内源 GLP1R，98.8%），内源位点的 SR/PESR 同样很高，只能靠 `site_homology.tsv` 区分 |

```bash
alnkit short scan --bam sample.markdup.bam \
    --query-chr target_ref --insert-chrs decoy_ref,homolog_ref,target_ref \
    --merge-window-bp 10 --min-sr-support 2 \
    --reference-fasta hg38i.fa \
    --homology-bed known_homologs.bed \
    --out-dir out
```

### 怎么读结果（人工判断，不自动剔除）

1. 先看 `sites.canonical.tsv` 的 `n_sr` / `n_pesr`：真位点通常支撑数比邻居高一个量级；
2. 再看 `site_homology.tsv` 的 `is_homology_site`：`yes` 的位点先怀疑是内源同源（若该外源本就
   含内源片段）；`no` 的位点是候选整合位点；
3. 两个相距 ~100–1000 bp 且 clip 方向相反的高支持位点，往往是**同一座位的两个连接点**
   （中间的参考片段仍被 reads 覆盖，不是缺失），不要合并调用；
4. 低支持（1–2 条）又落在外源同源区附近的位点，先在 IGV 上看 clip 方向与断点是否齐平再定。

## 5. 二代出图（位点示意图 / 覆盖度）

```bash
# 位点示意图：1 个断点（单侧/同断点双侧）或 2-3 个断点（多断点），画法自动判断
alnkit short insertion-plot --bam reads.sorted.bam --region chr1:99101-100900 \
    --breakpoint 100000 --out-dir out --name both

# 覆盖度：上广度、下深度，可选 BED 基因轨道
alnkit short coverage --bam reads.sorted.bam --region chr1:1-20000 --bins 400 \
    --bed genes.bed --out-dir out
```

两个命令都遵循“表是接口”：先写 TSV（`insertion_clip_reads` / `insertion_clip_consensus`、
`coverage_region` / `coverage_bins`），再由表渲染 PNG。参数说明见
[`alnkit short insertion-plot`](cli/alnkit-short-insertion-plot.md) 与
[`alnkit short coverage`](cli/alnkit-short-coverage.md)。

## 6. 三代断点→类→组装→落点

```bash
alnkit long classify  --bam reads.bam --reads reads.fastq --merged-fai merged.fa.fai \
    --exog insertA --exog insertB --outdir 03_classify
minimap2 -x ava-pb -c -t 8 merged.fa merged.fa > 04_overlap/overlaps.paf   # all-vs-all
alnkit long graph     --paf 04_overlap/overlaps.paf \
    --classification 03_classify/read_classification.tsv --outdir 04_overlap
alnkit long breakpoints --segments 03_classify/alignment_segments.tsv \
    --classification 03_classify/read_classification.tsv \
    --overlaps 04_overlap/overlap_edges.tsv --reads 03_classify/related_reads.fastq \
    --merged-fai merged.fa.fai --outdir 05_breakpoints
alnkit long assemble --classes-dir 05_breakpoints/classes --paf 04_overlap/overlaps.paf \
    --outdir 06_assembly        # 额外写 assembly_support.tsv（需 minimap2）
alnkit long locate --regions 05_breakpoints/breakpoint_classes.tsv \
    --sites 05_breakpoints/breakpoint_sites.tsv --reference-fasta reference.fa \
    --exog exog/insertA.fa --exog exog/insertB.fa \
    --assembly-dir 06_assembly --outdir 07_ribbon_positions
```

## 约定速记

- 坐标 0-based 半开；列名带基准后缀；区间串 1-based 闭。详见[坐标约定](concepts/coordinates.md)。
- 退出码：`0` 成功 / `1` 运行期错误 / `2` 参数或前置校验失败。
