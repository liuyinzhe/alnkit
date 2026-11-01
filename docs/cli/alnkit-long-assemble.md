# `alnkit long assemble`

每个类先拿一条最长的接头 read 作种子，再从 overlap 边里反复挑一条能往左/右延伸的 read，把非重叠悬垂序列接到 contig 两端。

## 用法

```bash
alnkit long assemble [全局选项] [--classes-dir …] [--paf …] [--outdir …] [--min-overlap …] [--min-identity …] [--agreement-bp …]
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
