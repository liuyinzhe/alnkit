# `alnkit plot ribbon`

方式 1：给 --reference + --query，内部跑 minimap2；方式 2：给 --paf（现成比对），直接画图。

## 用法

```bash
alnkit plot ribbon [全局选项] [--reference …] [--query …] [--query-file …] [--paf …] [--reference-name …] [--query-name …] [--outdir …] [--name …] [--title …] [--reference-position …] [--order …] [--order-file …] [--dpi …] [--min-label-bp …] [--label-coordinate …] [--max-labels-per-row …] [--placement-table …] [--region-table …] [--min-identity …] [--min-block-bp …] [--preset …] [--query-preset …] [--threads …]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--reference` | 否 | — | 参考侧，格式 标签=fasta路径[:片段名]（--paf 模式可省略） |
| `--query` | 否 | — | 查询侧，格式同上；可重复 |
| `--query-file` | 否 | — | 查询清单文件（每行一个 标签=路径[:片段名]） |
| `--paf` | 否 | — | 现成 PAF（方式 2；需带 -c 才有匹配列） |
| `--reference-name` | 否 | — | --paf 模式下指定哪条是参考行 |
| `--query-name` | 否 | — | 只画这些查询片段；可重复 |
| `--outdir` | 是 | — | 输出目录 |
| `--name` | 是 | — | 输出文件名前缀 |
| `--title` | 否 | — | 图标题（英文，避免缺字体） |
| `--reference-position` | 否 | `top` | 参考行放顶部还是底部（可选：`top`, `bottom`） |
| `--order` | 否 | — | 行顺序，用 > 分隔 |
| `--order-file` | 否 | — | 行顺序清单文件 |
| `--dpi` | 否 | `100` | PNG 分辨率（默认 100） |
| `--min-label-bp` | 否 | `300` | 只标长度不小于它的块（默认 300） |
| `--label-coordinate` | 否 | `query` | 标签用哪套坐标（可选：`query`, `reference`, `both`） |
| `--max-labels-per-row` | 否 | `8` | 每行最多标签数（默认 8） |
| `--placement-table` | 否 | — | 额外输出位置汇总表 |
| `--region-table` | 否 | — | 额外输出同源区域表 |
| `--min-identity` | 否 | `0.8` | 保留块的最小一致率（默认 0.8） |
| `--min-block-bp` | 否 | `200` | 保留块的最小长度（默认 200） |
| `--preset` | 否 | `asm20` | minimap2 预设（默认 asm20；auto=按长度分档） |
| `--query-preset` 标签=预设[ 额外参数] | 否 | — | 按查询标签覆盖比对参数；可重复 |
| `--threads` | 否 | `8` | 线程数（默认 8） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```bash
alnkit plot ribbon --paf aln.paf --reference-name assembly_a --outdir out --name ex4
     alnkit plot ribbon --reference "ins=exog.fa:insertA" --query "C01=asm.fa:oe_contig1" --outdir out --name exog_vs_class
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
