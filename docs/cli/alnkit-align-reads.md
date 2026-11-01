# `alnkit align reads`

等价于旧 shell：minimap2 建 .mmi -> minimap2 -ax map-hifi -Y -N 50 -> samtools sort -> samtools index。

## 用法

```bash
alnkit align reads [全局选项] [--reference …] [--reads …] [--outdir …] [--name …] [--threads …]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--reference` | 是 | — | 合并参考 FASTA（需 .fai） |
| `--reads` | 是 | — | 原始 reads（FASTQ/FASTA） |
| `--outdir` | 是 | — | 输出目录 |
| `--name` | 否 | `reads_to_merged` | 输出 BAM 文件名前缀 |
| `--threads` | 否 | `12` | 线程数（默认 12） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit align reads --reference 01_merge_ref/merged_reference.fa 
          --reads input.fastq --outdir 02_align --threads 12
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
