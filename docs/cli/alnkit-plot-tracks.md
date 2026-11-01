# `alnkit plot tracks`

一条 read 一行，横轴为 read 自身坐标；每个比对片段画成一个块。输入三选一：--bam / --paf / --segments-tsv。

## 用法

```bash
alnkit plot tracks [全局选项] [--bam …] [--paf …] [--segments-tsv …] [--read …] [--reads-file …] [--exog …] [--bam-mode …] [--min-mapq …] [--highlight-ref …] [--outdir …] [--name …] [--title …] [--dpi …] [--fig-width …] [--label-mode …] [--max-axis-bp …] [--include-secondary] [--full-read-names] [--exact-highlight]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--bam` | 否 | — | BAM（pysam 读取） |
| `--paf` | 否 | — | minimap2 PAF（需要 -c） |
| `--segments-tsv` | 否 | — | 已整理的片段表 |
| `--read` | 否 | — | 只画这些 read；可重复 |
| `--reads-file` | 否 | — | read 名列表文件（一行一个，空行/# 跳过） |
| `--exog` | 否 | — | 外源序列名；可重复或用逗号分隔 |
| `--bam-mode` | 否 | `auto` | BAM 取数方式：auto（默认）/scan/subset（可选：`auto`, `scan`, `subset`） |
| `--min-mapq` | 否 | `0` | 最小 MAPQ（默认 0） |
| `--include-secondary` | 否 | 关 | 包含 secondary 比对 |
| `--highlight-ref` | 否 | — | 要标注的参考断点，格式 chrom:pos；可重复 |
| `--outdir` | 是 | — | 输出目录 |
| `--name` | 否 | `read_tracks` | 输出文件名前缀 |
| `--title` | 否 | `Read-axis alignment segments` | 图标题（英文，避免缺字体） |
| `--full-read-names` | 否 | 关 | 左侧显示完整 read 名 |
| `--dpi` | 否 | `150` | PNG 分辨率（默认 150） |
| `--fig-width` | 否 | `14.0` | 图宽（英寸，默认 14） |
| `--exact-highlight` | 否 | 关 | 参考断点按 CIGAR 精确映射（BAM 输入） |
| `--label-mode` | 否 | `auto` | 片段标签：auto（默认）/inline/side（可选：`auto`, `inline`, `side`） |
| `--max-axis-bp` | 否 | — | 超过该长度折叠内容空隙（断轴） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit plot tracks --segments-tsv seg.tsv --outdir out
     alnkit plot tracks --bam reads.bam --read readA --exog insertA --highlight-ref chr1:1000001 --outdir out
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
