# `alnkit seg table`

输出规范列（默认）或历史 14 列（--columns legacy）。

## 用法

```bash
alnkit seg table [全局选项] [--bam …] [--out …] [--fai …] [--columns …] [--read-file …] [--skip-secondary]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--bam` | 是 | — | 比对 BAM（顺序扫描，不需要索引） |
| `--out` | 是 | — | 输出 TSV |
| `--fai` | 否 | — | 参考序列的 .fai（提供后带 chrom_length 列） |
| `--columns` | 否 | `canonical` | canonical（默认，带基准后缀）或 legacy（历史 14 列，逐字节兼容）（可选：`canonical`, `legacy`） |
| `--read-file` | 否 | — | 只导出这些 read（一行一个） |
| `--skip-secondary` | 否 | 关 | 跳过 secondary 比对 |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit seg table --bam reads.sorted.bam --fai ref.fa.fai --out seg.tsv
     alnkit seg table --bam reads.bam --out legacy.tsv --columns legacy
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
