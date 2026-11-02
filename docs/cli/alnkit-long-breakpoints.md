# `alnkit long breakpoints`

在 read 轴上按“紧邻片段”定义断点（基因组–外源 = vector S##；基因组–基因组 = genome G##a/G##b），两个断点共享 reads >= --min-shared 就合并为一类。

## 用法

```bash
alnkit long breakpoints [全局选项] [--segments …] [--classification …] [--overlaps …] [--reads …] [--merged-fai …] [--outdir …] [--site-gap …] [--min-shared …]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--segments` | 是 | — | alignment_segments.tsv |
| `--classification` | 是 | — | read_classification.tsv |
| `--overlaps` | 是 | — | overlap_edges.tsv |
| `--reads` | 是 | — | related_reads.fastq（输出每类子集用） |
| `--merged-fai` | 是 | — | 合并参考 .fai（取序列长度） |
| `--outdir` | 是 | — | 输出目录（建议 05_breakpoints） |
| `--site-gap` | 否 | `20` | 同一 site 的边界允许错位（默认 20 bp） |
| `--min-shared` | 否 | `2` | 合并断点所需最少共享 reads（默认 2） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit long breakpoints \
          --segments 03_classify/alignment_segments.tsv \
          --classification 03_classify/read_classification.tsv \
          --overlaps 04_overlap/overlap_edges.tsv \
          --reads 03_classify/related_reads.fastq \
          --merged-fai merged.fa.fai --outdir 05_breakpoints
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
