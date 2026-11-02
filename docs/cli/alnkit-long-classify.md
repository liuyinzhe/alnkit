# `alnkit long classify`

先按外源证据筛出目标 read，再逐段判定外源/基因组/可疑错配，最后把接头事件聚成位点与区域，并导出相关 reads 子集。

## 用法

```bash
alnkit long classify [全局选项] [--bam …] [--reads …] [--merged-fai …] [--exog …] [--outdir …] [--min-mapq …] [--min-exog-block-bp …] [--contig-edge-bp …] [--short-genome-block-bp …] [--junction-max-gap-bp …] [--site-gap-bp …] [--region-gap-bp …]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--bam` | 是 | — | reads -> 合并参考 的 BAM |
| `--reads` | 是 | — | 原始 fastq |
| `--merged-fai` | 是 | — | 合并参考 .fai（取序列长度） |
| `--exog` | 是 | — | 外源序列名，可重复 |
| `--outdir` | 是 | — | 输出目录 |
| `--min-mapq` | 否 | `20` | 比对质量下限（默认 20） |
| `--min-exog-block-bp` | 否 | `100` | 外源比对段最小长度（默认 100 bp） |
| `--contig-edge-bp` | 否 | `100` | contig 两端多近算边缘（默认 100 bp） |
| `--short-genome-block-bp` | 否 | `500` | 边缘处基因组段多短算可疑错配（默认 500 bp） |
| `--junction-max-gap-bp` | 否 | `10000` | 基因组段与外源段在 read 轴上多近算相邻（默认 10 kb） |
| `--site-gap-bp` | 否 | `20` | 同一位点的边界允许错位（默认 20 bp） |
| `--region-gap-bp` | 否 | `200000` | 同染色体位点归为同一区域的最大间距（默认 200 kb） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```text
上游 BAM 的生成要求（断点证据来自 CIGAR 的 S，所以这两条是硬要求）：
  minimap2 -ax map-hifi -Y -N 50 -t 8 merged.fa reads.fq \
      | samtools sort -o reads.sorted.bam -  &&  samtools index reads.sorted.bam
  -Y     补充比对用 soft clip（否则被 hard clip，取不到断点/归属）
  -N 50  保留足够的次要比对（共有 read 的归属判定要用）
```

```bash
alnkit long classify --bam reads_to_merged.sorted.bam 
          --reads input.fastq --merged-fai merged.fa.fai 
          --exog insertA --exog insertB --outdir 03_classify
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
