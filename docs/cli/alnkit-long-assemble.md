# `alnkit long assemble`

每个类先拿一条最长的接头 read 作种子，再从 overlap 边里反复挑一条能往左/右延伸的 read，把非重叠悬垂序列接到 contig 两端。

产物：<outdir>/<类>/overlap_extension/contigs.fasta + assembly_stats.tsv，外加 <outdir>/assembly_support.tsv（默认开，需 minimap2；--no-support-report 关）。
接缝门槛默认关闭：--min-junction-depth-ratio 0.1 让“深度骤降处断开而不是硬拼”。
方案 C（--recruit-edges，默认关）：贪心之后再拿 contig 末端当探针招募 read 继续延长。

## 用法

```bash
alnkit long assemble [全局选项] [--classes-dir …] [--paf …] [--outdir …] [--min-overlap …] [--min-identity …] [--agreement-bp …] [--min-junction-depth-ratio …] [--bridge-window-bp …] [--recruit-min-support …] [--recruit-rounds …] [--no-support-report …] [--threads …] [--recruit-edges]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--classes-dir` | 是 | — | 05_breakpoints/classes |
| `--paf` | 是 | — | all-vs-all PAF |
| `--outdir` | 是 | — | 输出目录 |
| `--min-overlap` | 否 | `500` | 延伸所需最小 overlap（默认 500 bp） |
| `--min-identity` | 否 | `0.9` | 延伸所需最小一致率（默认 0.90） |
| `--agreement-bp` | 否 | `200` | 多个邻居给出同一位置时允许的起点差（默认 200 bp） |
| `--min-junction-depth-ratio` | 否 | `0.0` | 接缝局部深度 ÷ contig 平均深度 低于此比值就拒绝延伸（默认 0 = 关闭）。实测 0.1/0.2/0.3 会把 5 条 contig 切碎成 10/16/17 条而锚点回收不变（8/8），所以默认不启用；需要“宁可断开也不拼错”时再开 |
| `--bridge-window-bp` | 否 | `200` | 接缝判定窗口总宽（默认 200 bp，即接缝两侧各 100 bp）；只影响这道门槛，不改 assembly_support.tsv 的统计窗口（固定 200 bp） |
| `--recruit-edges` | 否 | 关 | 方案 C（借鉴 DEGAP）：贪心组装后，再用 contig 末端 500 bp 当探针招募 read，把悬垂移植过来继续延长（每轮必须延长 ≥100 bp 才接受；默认关闭，需要 minimap2） |
| `--recruit-min-support` | 否 | `2` | 一个悬垂至少要有多少条 read 支持才移植（默认 2） |
| `--recruit-rounds` | 否 | `6` | 端点招募最多几轮（默认 6） |
| `--no-support-report` | 否 | 关 | 不生成 assembly_support.tsv（默认会用 minimap2 回映 reads 做支撑度统计） |
| `--threads` | 否 | `4` | minimap2 回映线程数（默认 4） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit long assemble \
          --classes-dir 05_breakpoints/classes \
          --paf 04_overlap/overlaps.paf --outdir 06_assembly
```

```text
加了支撑度报告后，先用
  column -t -s$'\t' 06_assembly/assembly_support.tsv
看每类 contig 的 covered_fraction / mean_depth / depth_cliffs / status。
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
