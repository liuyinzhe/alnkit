# `alnkit long assemble`

每个类先拿一条最长的接头 read 作种子，再从 overlap 边里反复挑一条能往左/右延伸的 read，把非重叠悬垂序列接到 contig 两端。

## 用法

```bash
alnkit long assemble [全局选项] [--classes-dir …] [--paf …] [--outdir …] [--min-overlap …] [--min-identity …] [--agreement-bp …] [--min-junction-depth-ratio …] [--bridge-window-bp …] [--no-support-report …] [--threads …]
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
| `--no-support-report` | 否 | 开 | 不生成 assembly_support.tsv（默认会用 minimap2 回映 reads 做支撑度统计） |
| `--threads` | 否 | `4` | minimap2 回映线程数（默认 4） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit long assemble 
          --classes-dir 05_breakpoints/classes 
          --paf 04_overlap/overlaps.paf --outdir 06_assembly
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
