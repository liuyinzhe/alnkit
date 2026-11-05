# `alnkit short coverage`

按区间统计覆盖度：逐碱基深度 -> 分箱 -> 上图广度 / 下图深度（可选 BED 基因轨道）。出图只读 TSV：先写 region/bins 两张表，再由表渲染。

## 用法

```bash
alnkit short coverage [全局选项] [--bam …] [--region …] [--out-dir …] [--name …] [--bins …] [--bed …] [--title …] [--width …] [--dpi …] [--bases-table …] [--min-mapq …] [--min-len …] [--min-depth …] [--max-depth …] [--ff …] [--rf …] [--list]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--bam` PATH | 是 | — | BAM/CRAM（按区间取数需要索引） |
| `--region` CHROM[:START-END] | 否 | — | 区间（1-based 闭区间）；只写 contig 名就是整条序列；除 --list 外必填 |
| `--out-dir` | 否 | `.` | 输出目录，默认当前目录 |
| `--name` | 否 | — | 输出文件名前缀，默认 <bam 主名>.<contig>[_start-end] |
| `--bins` INT | 否 | — | bin 个数；默认约 25 bp/bin |
| `--bed` BED | 否 | — | 基因轨道 BED（给了才画） |
| `--title` | 否 | — | 覆盖默认标题（英文，避免缺字体） |
| `--width` | 否 | `10.0` | 图宽（英寸，默认 10） |
| `--dpi` | 否 | `200` | PNG 分辨率（默认 200） |
| `--bases-table` TSV | 否 | — | 额外写逐碱基深度表 |
| `--list` | 否 | 关 | 只列出参考序列与比对上的 reads 数后退出 |
| `--min-mapq` | 否 | `0` | 最小 MAPQ（默认 0） |
| `--min-len` | 否 | `0` | 最小 query 长度（默认 0） |
| `--min-depth` | 否 | `1` | 低于该深度的位置不算覆盖（默认 1） |
| `--max-depth` | 否 | `1000000` | 每个位置的深度上限，0 = 不限（默认 1000000） |
| `--ff` | 否 | `UNMAP,SECONDARY,QCFAIL,DUP` | 带这些 flag 的 read 跳过（默认 UNMAP,SECONDARY,QCFAIL,DUP；空串=不过滤） |
| `--rf` | 否 | — | 缺这些 flag 的 read 跳过（默认空） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```text
指标:
  只有 M/=/X 给深度加 1；D/N 占参考位置但没有碱基，不算覆盖；S/I/H 不占参考位置。
  depth >= --min-depth 的位置才算“被覆盖”；平均深度的分母是区间长度。
```

```bash
alnkit short coverage --bam reads.bam --region insertA --out-dir out
  alnkit short coverage --bam reads.bam --region chr1:1-20000 --bins 400 \
      --bed genes.bed --out-dir out --name region
  alnkit short coverage --bam reads.bam --list
```

```text
产物:
  <name>.png
  + <name>.region.tsv（coverage_region）
  + <name>.bins.tsv（coverage_bins）
  另加 --bases-table 时再写逐碱基表（coverage_bases）。
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
