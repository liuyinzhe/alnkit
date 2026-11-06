# `alnkit short scan`

一条跨越“参考基因组 <-> 外源插入序列”的嵌合比对可定位断点；本命令按断点聚合证据（SR / PESR / DUP）。

## 用法

```bash
alnkit short scan [全局选项] [--bam …] [--bam-dir …] [--bam-glob …] [--out-dir …] [--sample-name …] [--query-chr …] [--insert-chrs …] [--insert-origin-policy …] [--merge-window-bp …] [--min-mapq …] [--max-mismatch …] [--max-softclip-blocks …] [--min-clip-length …] [--mismatch-source …] [--on-missing-tag …] [--min-sr-support …] [--interference-scope …] [--reference-fasta …] [--homology-bed …] [--homology-window-bp …] [--homology-min-identity …] [--homology-min-length …] [--bam-threads …] [--no-recursive] [--keep-secondary] [--no-split-match] [--no-candidate-gate] [--no-duplicate-filter] [--check-markdup] [--require-dup-flags] [--legacy-coordinates] [--legacy-profile] [--self-test]
```

## 参数

| 参数 | 必填 | 默认 | 说明 |
|---|---|---|---|
| `--bam` PATH | 否 | — | 直接指定 BAM（可重复）。与 --bam-dir 二选一 |
| `--bam-dir` DIR | 否 | `当前目录` | 扫描 BAM 的目录，默认当前目录 |
| `--bam-glob` PATTERN | 否 | `*.sorted.bam` | 在 --bam-dir 下匹配 BAM 的通配符，默认 *.sorted.bam |
| `--no-recursive` | 否 | 关 | 不递归子目录（默认递归） |
| `--out-dir` DIR | 否 | `当前目录` | 结果输出目录，默认当前目录 |
| `--sample-name` NAME | 否 | — | 输出前缀（所有产物文件名的开头），默认取 BAM 文件名首个 . 之前的部分 |
| `--query-chr` NAME | 否 | — | 本次要出结果的外源插入序列名（必填，除非 --self-test） |
| `--insert-chrs` NAME1,NAME2,... | 否 | — | 参考序列里全部外源/同源序列名，逗号分隔；默认只含 --query-chr |
| `--insert-origin-policy` | 否 | `owned` | 共有 read 的区分策略：owned（默认）/owned-strict/lenient/consensus/legacy（可选：`owned`, `owned-strict`, `lenient`, `consensus`, `legacy`） |
| `--merge-window-bp` INT | 否 | `0` | 位点窗口合并：与簇内最小坐标相距 <= 该值的断点合并为一个位点（贪心、簇宽不超窗），0=不合并（默认，历史产物逐字节不变）。推荐 10–20：测序本身不产生这种抖动，抖动来自比对器对 clip 端点的判定，实测同一连接点的边界会散在 ~20 bp 内；取 3 会把同一连接点拆成多行、峰值支持数掉 10%%；反过来调到几十-几百 bp 会把同一座位两侧的两个连接点并成一行 |
| `--min-mapq` INT | 否 | `20` | 比对质量 <= 该值则丢弃，默认 20 |
| `--max-mismatch` INT | 否 | `3` | 错配数 > 该值则丢弃，默认 3 |
| `--max-softclip-blocks` INT | 否 | `1` | soft clip 段数 > 该值则丢弃，默认 1 |
| `--min-clip-length` INT | 否 | `0` | 作为断点证据时要求的最小 soft clip 长度，默认 0（不限制） |
| `--keep-secondary` | 否 | 关 | 保留 secondary 比对（默认丢弃） |
| `--mismatch-source` | 否 | `nm` | 错配数来源：nm=NM 标签（默认）；md=由 MD 推导（历史行为）。两者都不自动回退：选 nm 却没有 NM 时按 --on-missing-tag 丢弃或报错（可选：`nm`, `md`） |
| `--on-missing-tag` | 否 | `skip` | 缺少 NM/MD 标签时：skip=丢弃该比对（默认）；error=报错（可选：`skip`, `error`） |
| `--no-split-match` | 否 | 关 | 关闭「目标序列上要求存在等长 split 比对」的校验（默认开启） |
| `--no-candidate-gate` | 否 | 关 | 关闭候选 read 闸门（等价 2.4.0 全量逐条处理）；legacy 策略下自动关闭 |
| `--no-duplicate-filter` | 否 | 关 | 保留被 markdup 标记为重复的 read 作为证据（默认剔除，只写 dup 列） |
| `--min-sr-support` INT | 否 | `1` | 位点至少要有多少条支持 read（SR+PESR）才输出，默认 1（不过滤，保持历史行为）；设 2–3 可去掉单 read 噪声位点 |
| `--interference-scope` | 否 | `all` | 同源干扰的判据范围：all=支持 read 命中任何插入序列都算干扰（默认，历史口径，多外源面板下偏灵敏）；other=只算 query 之外的插入（多外源面板更精确）。注意：两者都只比较“插入之间”，无法识别“外源与基因组同源”造成的干扰（可选：`all`, `other`） |
| `--reference-fasta` | 否 | — | 含插入序列与基因组的参考 FASTA；给了它就把插入序列回比基因组，自动找出“外源同源区”并标注位点（需要 minimap2） |
| `--homology-bed` BED | 否 | — | 已知的内源同源基因/区域（BED，0-based 半开，可给多次）；与 --reference-fasta 可并用 |
| `--homology-window-bp` INT | 否 | `50000` | 局部回比时位点两侧各取多长（默认 50000 bp）；cDNA 样插入的外显子分散在几十 kb 里，窗口太小只能比到单个外显子 |
| `--homology-min-identity` FLOAT | 否 | `0.9` | 回比时认定的最小一致率，默认 0.90 |
| `--homology-min-length` INT | 否 | `200` | 回比时认定的最小比对长度（bp），默认 200（过滤偶然的短命中） |
| `--check-markdup` | 否 | 关 | 核对 markdup 标记是否可信：把 read 按「比对坐标 + 错配数」分组，报告全未标记组（markdup 没做/不完整）、标记不一致组、以及单条被标 dup 的 read |
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

```bash
alnkit short scan --bam-dir ./bam --query-chr target_ref \
      --insert-chrs decoy_ref,homolog_ref,target_ref \
      --merge-window-bp 10 --out-dir out
  alnkit short scan --self-test
```

```text
产物:
  三张 legacy .xls（列名与旧脚本一致）+ <sample>.sites.canonical.tsv（机器接口）。
  逐列含义与读法：docs/concepts/short-scan-output.md
  注意: reads_info 的 legacy 表头把 stand 与 mapping_quality 合并成一列（11 列），
  数据行其实是 12 列——按表头名解析会错位，按位置或用 canonical 表。
```

```text
错配来源与 markdup:
  --mismatch-source nm（默认）读 NM；md 用 MD 推导（历史口径）。都不自动回退，
  标签缺失时按 --on-missing-tag 丢弃或报错；--legacy-profile 隐含 md。
  markdup 的 FLAG 0x400 与错配统计无关，只决定 dup 两列；用 --check-markdup 核对
  标记是否可信（按两端 5' 坐标分组，再用 MD 判断被标 dup 的 read 是否真同序列）。
```

## 退出码

`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。
