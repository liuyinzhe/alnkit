# `alnkit short scan`

一条跨越“参考基因组 <-> 外源插入序列”的嵌合比对可定位断点；本命令按断点聚合证据（SR / PESR / DUP）。

## 用法

```bash
alnkit short scan [全局选项] [--bam …] [--bam-dir …] [--bam-glob …] [--out-dir …] [--sample-name …] [--query-chr …] [--insert-chrs …] [--insert-origin-policy …] [--merge-window-bp …] [--min-mapq …] [--max-mismatch …] [--max-softclip-blocks …] [--min-clip-length …] [--mismatch-source …] [--on-missing-tag …] [--bam-threads …] [--no-recursive] [--keep-secondary] [--no-split-match] [--no-candidate-gate] [--no-duplicate-filter] [--require-dup-flags] [--legacy-coordinates] [--legacy-profile] [--self-test]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--bam` PATH | 否 | — | 直接指定 BAM（可重复）。与 --bam-dir 二选一 |
| `--bam-dir` DIR | 否 | `当前目录` | 扫描 BAM 的目录，默认当前目录 |
| `--bam-glob` PATTERN | 否 | `*.sorted.bam` | 在 --bam-dir 下匹配 BAM 的通配符，默认 *.sorted.bam |
| `--no-recursive` | 否 | 关 | 不递归子目录（默认递归） |
| `--out-dir` DIR | 否 | `当前目录` | 结果输出目录，默认当前目录 |
| `--sample-name` NAME | 否 | — | 样本名（输出文件名前缀），默认取 BAM 文件名首个 . 之前的部分 |
| `--query-chr` NAME | 否 | — | 本次要出结果的外源插入序列名（必填，除非 --self-test） |
| `--insert-chrs` NAME1,NAME2,... | 否 | — | 参考序列里全部外源/同源序列名，逗号分隔；默认只含 --query-chr |
| `--insert-origin-policy` | 否 | `owned` | 共有 read 的区分策略：owned（默认）/owned-strict/lenient/consensus/legacy（可选：`owned`, `owned-strict`, `lenient`, `consensus`, `legacy`） |
| `--merge-window-bp` INT | 否 | `0` | 位点窗口合并：相距 <= 该值的断点合并为一个位点，0=不合并（默认）；推荐 10 |
| `--min-mapq` INT | 否 | `20` | 比对质量 <= 该值则丢弃，默认 20 |
| `--max-mismatch` INT | 否 | `3` | 错配数 > 该值则丢弃，默认 3 |
| `--max-softclip-blocks` INT | 否 | `1` | soft clip 段数 > 该值则丢弃，默认 1 |
| `--min-clip-length` INT | 否 | `0` | 作为断点证据时要求的最小 soft clip 长度，默认 0（不限制） |
| `--keep-secondary` | 否 | 关 | 保留 secondary 比对（默认丢弃） |
| `--mismatch-source` | 否 | `nm` | 错配数来源：nm=NM 标签（默认）；md=由 MD 推导（历史行为）（可选：`nm`, `md`） |
| `--on-missing-tag` | 否 | `skip` | 缺少 NM/MD 标签时：skip=丢弃该比对（默认）；error=报错（可选：`skip`, `error`） |
| `--no-split-match` | 否 | 关 | 关闭「目标序列上要求存在等长 split 比对」的校验（默认开启） |
| `--no-candidate-gate` | 否 | 关 | 关闭候选 read 闸门（等价 2.4.0 全量逐条处理）；legacy 策略下自动关闭 |
| `--no-duplicate-filter` | 否 | 关 | 保留被 markdup 标记为重复的 read 作为证据（默认剔除，只写 dup 列） |
| `--require-dup-flags` | 否 | 关 | BAM 里没有任何 FLAG 0x400（未做 markdup）时报错退出（退出码 2） |
| `--bam-threads` INT | 否 | `1` | BAM 解压线程数，默认 1；大文件可调到 2-4 |
| `--legacy-coordinates` | 否 | 关 | 使用历史的 M 块选取规则（取最后一个 M 块） |
| `--legacy-profile` | 否 | 关 | 一键复现历史结果（= legacy 坐标 + legacy 归属 + 由 MD 推导错配） |
| `--self-test` | 否 | 关 | 只跑内置自检（不读 BAM） |

## 全局选项（每个子命令都有）

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--log-level` LEVEL | 否 | `INFO` | 日志级别：DEBUG / INFO / WARNING / ERROR |
| `--quiet` | 否 | 关 | 只输出错误（等价 --log-level ERROR） |

## 例子

```text
上游 BAM 的比对要求（详见原脚本 README）:
  必须 -Y（补充比对用 soft clip，只从 CIGAR 的 S 取断点）；
  禁止 -M（会把补充比对标成 secondary，证据变少）；
  建议 -a -h 20；比对后不要用 samtools -F 2048 / -F 256 / -q 过滤。
```

```text
示例:
  alnkit short scan --bam-dir ./bam --query-chr target_ref \
      --insert-chrs decoy_ref,homolog_ref,target_ref \
      --merge-window-bp 10 --out-dir out
  alnkit short scan --self-test
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
