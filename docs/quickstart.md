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

## 5. 三代断点→类→组装→落点

```bash
alnkit long classify  --bam reads.bam --reads reads.fastq --fai merged.fa.fai \
    --exog insertA --exog insertB --outdir 03_classify
alnkit long breakpoints --segments 03_classify/alignment_segments.tsv \
    --classification 03_classify/read_classification.tsv \
    --overlaps 04_overlap/overlap_edges.tsv --reads 03_classify/related_reads.fastq \
    --merged-fai merged.fa.fai --outdir 05_breakpoints
alnkit long assemble --classes-dir 05_breakpoints/classes --paf 04_overlap/overlaps.paf \
    --outdir 06_assembly
alnkit long locate --regions 05_breakpoints/breakpoint_classes.tsv \
    --sites 05_breakpoints/breakpoint_sites.tsv --reference-fasta reference.fa \
    --exog exog/insertA.fa --exog exog/insertB.fa \
    --assembly-dir 06_assembly --outdir 07_ribbon_positions
```

## 约定速记

- 坐标 0-based 半开；列名带基准后缀；区间串 1-based 闭。详见[坐标约定](concepts/coordinates.md)。
- 退出码：`0` 成功 / `1` 运行期错误 / `2` 参数或前置校验失败。
