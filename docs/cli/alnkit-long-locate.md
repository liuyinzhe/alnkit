# `alnkit long locate`

为每个类构造查询序列（断点旁 50 bp 锚点 + 全部外源），与类组装比对后逐条 contig 出图，并把锚点/外源在 contig 上的落点还原成参考基因组坐标。

## 用法

```bash
alnkit long locate [全局选项] [--regions …] [--sites …] [--reference-fasta …] [--exog …] [--assembly-dir …] [--outdir …] [--threads …] [--flank-bp …]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--regions` | 是 | — | breakpoint_classes.tsv |
| `--sites` | 是 | — | breakpoint_sites.tsv |
| `--reference-fasta` | 是 | — | 参考基因组 fasta（需 .fai） |
| `--exog` | 是 | — | 外源 fasta，可重复 |
| `--assembly-dir` | 是 | — | 06_assembly |
| `--outdir` | 是 | — | 输出目录 |
| `--threads` | 否 | `12` | minimap2 线程数（默认 12） |
| `--flank-bp` | 否 | `50` | 断点锚点长度（默认 50 bp） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit long locate 
          --regions 05_breakpoints/breakpoint_classes.tsv 
          --sites 05_breakpoints/breakpoint_sites.tsv 
          --reference-fasta reference.fa 
          --exog exog/insertA.fa --exog exog/insertB.fa 
          --assembly-dir 06_assembly --outdir 07_ribbon_positions
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
