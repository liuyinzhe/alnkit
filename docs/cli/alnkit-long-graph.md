# `alnkit long graph`

读 all-vs-all PAF（minimap2 -x ava-pb -c），每对 read 只保留最好的一条 overlap，再求全局连通分量。方案 A 的类别由断点共享 reads 决定，这里只用于观察 overlap 结构。

## 用法

```bash
alnkit long graph [全局选项] [--paf …] [--classification …] [--outdir …] [--min-overlap …] [--min-identity …]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--paf` | 是 | — | all-vs-all PAF（minimap2 -x ava-pb -c） |
| `--classification` | 是 | — | read_classification.tsv |
| `--outdir` | 是 | — | 输出目录 |
| `--min-overlap` | 否 | `300` | 边所需最小 overlap（默认 300 bp） |
| `--min-identity` | 否 | `0.85` | 边所需最小一致率（默认 0.85） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit long graph --paf overlaps.paf \
          --classification 03_classify/read_classification.tsv --outdir 04_overlap
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
