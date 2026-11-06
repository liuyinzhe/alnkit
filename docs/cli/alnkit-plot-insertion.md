# `alnkit plot insertion`

适用于短读长（二代 / PE）比对：从 CIGAR 首尾的 soft clip 取断点证据，按“给了几个断点 + 每个断点有几侧 clip”自动选画法。长读长（三代）的比对片段请用 `plot tracks` / `plot ribbon`。出图只读 TSV：先写 clips/consensus 两张表，再由表渲染。

## 用法

```bash
alnkit plot insertion [全局选项] [--bam …] [--region …] [--breakpoint …] [--out-dir …] [--name …] [--dpi …] [--no-consensus-text] [--no-formula]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--bam` PATH | 是 | — | 带 soft clip 的 BAM（需 .bai 索引） |
| `--region` CHROM:START-END | 是 | — | 窗口（1-based 闭区间），只画窗口内的 read |
| `--breakpoint` POS | 是 | — | 断点（1-based，= 断点左侧碱基的坐标）；给 1 个=单断点图，2-3 个=多断点图；可重复 |
| `--out-dir` | 否 | `.` | 输出目录，默认当前目录 |
| `--name` | 否 | `insertion` | 输出文件名前缀，默认 insertion |
| `--dpi` | 否 | `150` | PNG 分辨率（默认 150） |
| `--no-consensus-text` | 否 | 关 | 不在图底部写一致性序列文字 |
| `--no-formula` | 否 | 关 | 不在图右下角写 Consensus fraction 公式 |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```text
画法:
  1 个断点、只有一侧 clip -> 单侧图（未知条 + 彩色块贴在断点一侧）
  1 个断点、两侧都有 clip -> 同断点双侧图（左块 + 合并未知区 + 右块，连线收到断点）
  2-3 个断点              -> 多断点图（每断点各自 [未知][块]，块间强制留间隔）
```

```text
配色: 左 clip 蓝 / 右 clip 紫 / 未知区灰斜线 / 断点红。
```

```bash
alnkit plot insertion --bam reads.bam --region chr1:99201-100800 \
      --breakpoint 100000 --out-dir out --name demo
  alnkit plot insertion --bam reads.bam --region chr1:699101-700900 \
      --breakpoint 700000 --breakpoint 700230 --out-dir out --name two
```

```text
产物:
  <name>.png（示意图）+ <name>.clips.tsv（insertion_clip_reads）
  + <name>.consensus.tsv（insertion_clip_consensus）
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
