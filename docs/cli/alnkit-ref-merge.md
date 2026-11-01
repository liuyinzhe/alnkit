# `alnkit ref merge`

外源按顺序转大写后拼到参考末尾；同时输出 reference.fa(.fai)、reference_headers.tsv、exogenous_uppercased.fa。

## 用法

```bash
alnkit ref merge [全局选项] [--reference …] [--exog …] [--outdir …]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--reference` | 是 | — | 参考基因组 FASTA（可 .gz） |
| `--exog` | 是 | — | 外源 FASTA，可重复（按给出顺序拼接） |
| `--outdir` | 是 | — | 输出目录 |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit ref merge --reference ref.fa.gz 
          --exog exog/insertA.fa --exog exog/insertB.fa --outdir 01_merge_ref
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
