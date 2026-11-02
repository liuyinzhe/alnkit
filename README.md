# alnkit

BAM/PAF 的插入证据工具箱：断点提取 + 独立绘图。
二代与三代各一条链路，绘图工具可独立使用。


## 上游准备（自己跑，两步）

`alnkit` 从"合并参考 + BAM"开始，这两步用现成命令做即可：

```bash
# 1) 合并参考：外源序列转大写后拼到参考末尾，再建索引
cat ref.fa <(awk '/^>/{print;next}{print toupper($0)}' exogA.fa exogB.fa) > merged.fa
samtools faidx merged.fa

# 2) 比对（长读长）。-Y 与 -N 50 不是可选项，理由见下
minimap2 -ax map-hifi -Y -N 50 -t 8 merged.fa reads.fq \
    | samtools sort -o reads.sorted.bam -
samtools index reads.sorted.bam
```

| 参数 | 为什么必须 | 不加会怎样 |
|---|---|---|
| `-Y` | 补充比对用 **soft** clip 表达；`long classify` 只从 CIGAR 的 `S` 取断点与归属 | 变成 hard clip，断点证据直接丢失 |
| `-N 50` | 保留足够的次要比对，用于"共有 read"归属判定 | 同源区段的 read 归属判不出来，断点被错分 |

短读长同理：`bwa mem -Y`（**不要** `-M`，否则补充比对被标成 secondary），比对后不要用
`samtools view -F 2048/-F 256/-q` 过滤掉它们。

## 三条链路

> 示例里的序列名（`target_ref`/`homolog_ref`/`insertA`…）是**占位名**，换成你自己 BAM/FASTA header 里真实存在的名字；`long classify` 会校验并列出可用名。

```bash
# 三代：参考 → 比对 → 分类 → overlap → 断点/类 → 组装 → 落点
# 上游两步见上面「上游准备」：得到 merged.fa(.fai) 与 reads.sorted.bam
alnkit long classify  --bam reads.sorted.bam --reads input.fastq \
                      --merged-fai merged.fa.fai --exog insertA --exog insertB \
                      --outdir 03_classify
alnkit long graph     --paf 04_overlap/overlaps.paf \
                      --classification 03_classify/read_classification.tsv --outdir 04_overlap
alnkit long breakpoints --segments 03_classify/alignment_segments.tsv \
                      --classification 03_classify/read_classification.tsv \
                      --overlaps 04_overlap/overlap_edges.tsv --reads 03_classify/related_reads.fastq \
                      --merged-fai merged.fa.fai --outdir 05_breakpoints
alnkit long assemble  --classes-dir 05_breakpoints/classes --paf 04_overlap/overlaps.paf \
                      --outdir 06_assembly        # 额外写 assembly_support.tsv（需 minimap2）
alnkit long locate    --regions 05_breakpoints/breakpoint_classes.tsv \
                      --sites 05_breakpoints/breakpoint_sites.tsv --reference-fasta ref.fa \
                      --exog exog/h.fa --exog exog/l.fa --assembly-dir 06_assembly \
                      --outdir 07_ribbon

# 二代：BAM → 三张证据表 + 规范位点表
alnkit short scan     --bam-dir ./bams --query-chr target_ref \
                      --insert-chrs decoy_ref,homolog_ref,target_ref --merge-window-bp 10 --out-dir out

# 绘图（零上游依赖）
alnkit plot tracks    --segments-tsv seg.tsv --outdir out --name tracks
alnkit plot ribbon    --paf aln.paf --reference-name assembly_a --outdir out --name ex4
```

完整的参数说明见文档站的 **CLI 参考**（由 `--help` 自动生成）。
>
> `alnkit long assemble` 修掉了旧实现的非确定性（旧版遍历 `set`，同一输入会随
> `PYTHONHASHSEED` 漂）。除 `layout.tsv`（本包不再输出，见「已知限制」）与
> `assembly_support.tsv`（本包新增的回验表）之外，其余产物与旧脚本逐字节一致。
> 计划与契约见仓库内 `docs/`。

## 安装

```bash
pip install alnkit            # 分析链路（核心依赖：pysam）
pip install "alnkit[plot]"    # 额外装绘图（matplotlib）
```

外部命令按需安装（不在 wheel 里）：`minimap2`（>= 2.24，长读长链路的 `long classify` /
`long locate` / `plot ribbon` 内部比对 / `long assemble` 默认的支撑度报告都要它）；
`samtools`（>= 1.16，**可选**，只有 `plot tracks --bam-mode subset` 会调用外部 `samtools`）。

## 快速开始

```bash
alnkit --help          # 命令总览
alnkit doctor          # 打印 Python / 依赖 / 外部命令的可用性与版本
alnkit selftest        # 内置断言（坐标/表契约 + 二代断点口径，不需要输入数据），全绿即环境正常
```

## 现在能做什么

| 命令 | 作用 | 状态 |
|---|---|---|
| `alnkit doctor` | 环境自检：pysam、matplotlib、minimap2、samtools | 可用 |
| `alnkit selftest` | 坐标换算与表契约断言 | 可用 |
| `alnkit seg table` | BAM → 比对段表（规范列 / `--columns legacy` 历史 14 列） | 可用 |
| `alnkit plot tracks` | 单 read 轨道图（BAM / PAF / 段表） | 可用 |
| `alnkit plot ribbon` | 共线性丝带图（内部比对 / 现成 PAF） | 可用 |
| `alnkit short scan` | 二代插入断点证据（SR / PESR / DUP）三张表 | 可用 |
| `alnkit long classify` | 三代目标 reads 分类 + 接头位点/区域 | 可用 |
| `alnkit long graph` | read-vs-read overlap 边与连通分量 | 可用 |
| `alnkit long breakpoints` | 三代断点（S##/G##）+ 共享 reads 合并成类 + 类池 | 可用 |
| `alnkit long assemble` | 三代每类 junction 局部延伸组装 + 组装支撑度报告 | 可用 |
| `alnkit long locate` | 三代 50 bp 锚点 + 外源 → 丝带图 + 参考坐标还原 | 可用 |

## 关键约定（写文档和写脚本前必读）

- **坐标**：内部一律 0-based 半开 `[start, end)`；输出列名必须带 `_0based` / `_1based` 后缀；
  区间展示串一律 1-based 闭区间 `chrom:start-end`。见 `docs/concepts/coordinates.md`。
- **表是接口**：所有分析先落成 TSV，图从 TSV 出；列的增删受契约测试约束。
- **退出码**：`0` 成功 / `1` 运行期错误 / `2` 参数或前置校验失败。

## 为什么自带组装器（而不是让你用 hifiasm / miniasm）

`long assemble` 不是通用组装器，它只做一件事：**从 junction read 出发，用 read-vs-read overlap
把断点两侧往外延伸**，目的是让下游 `long locate` 能把 50 bp 断点锚点比回 contig、还原参考坐标。

在一套多类、每类百余条 HiFi reads 的基准上实测三种方案，
判据是**下游真正要的东西**：全部断点锚点能否被还原出参考坐标。

| 方案 | 锚点落回 | contig 数 | 总长 | 备注 |
|---|---|---|---|---|
| **`alnkit long assemble`** | **全部** | 5 | 340 kb | junction 作种子，天然沿着断点延伸 |
| hifiasm（每类 reads 各自组装） | 只落回 3/8 | — | — | 通用组装器在重复/嵌合处会断开，锚点正好落在断点缝上 |
| miniasm（吃 `overlaps.paf`） | **0/8** | — | — | 只保留覆盖度足够的 tiling path，断点两侧被裁掉 |

结论：**两者都无法替代**，因此保留自带组装器。若你要的是"整段插入的完整序列"而不是"断点验证"，
hifiasm 仍然是更合适的选择——`long locate` 对任何来源的 `06_assembly/<C>/overlap_extension/contigs.fasta`
都能用，两者可以并存（先 hifiasm 看全貌，再 `long assemble` 看断点）。

**仍然存在的限制**：本组装器只保证断点附近正确；远离断点的延伸无法验证，
所以它**不输出逐 read 的布局表**（`layout.tsv` 已移除），只给 contig 序列与统计。

## 已知限制

- **`long assemble` 的组装结果需要复核**。2026-11-01 修掉了一个反链 overlap 落位差
  `aligned_bp-1` 的坐标 bug（这类数据里相当比例的 overlap 是反链）。`layout.tsv` 已不再输出：
  它记录的坐标属于“当时那一轮”的坐标系，与最终 contig 直接对比只有约 28% 逐位一致，
  这是**输出记账**问题（算法本身的延伸决策自洽），不是拼错——这也是改用
  `assembly_support.tsv`（比对回验）做外部证据的原因。
  **把 `long assemble` 的 contig 当作“参考级”证据前，先看支撑度状态**；`long locate` 的锚点
  命中与参考坐标还原不受影响（它只用 50 bp 锚点与全 contig 比对）。
  需要“宁可断开也不拼错”时加 `--min-junction-depth-ratio 0.1`（默认 `0` 关闭；
  实测 0.1 会把 5 条 contig 切成 10 条而锚点回收不变，见 CHANGELOG）。
- **`assembly_support.tsv` 需要 minimap2**（`-ax map-hifi -Y`）。缺失时只警告并跳过该报告，
  其余产物不受影响；`--no-support-report` 可静默跳过。状态里的 `chimera_suspect` 是按
  “相邻窗口深度骤降到 1/5 以下且至少 3 处”判的，属于提示而不是定论。
- **`short scan` 的 legacy 三表列名保持历史写法**（`Insert_Position` 等），不要按列名推断基准；
  要机器可读、基准明确的版本请用同一目录下的 `<sample>.sites.canonical.tsv`
  （带 `_1based`/`_0based` 后缀，与 legacy 主表按 `(chrom, insert_pos_1based)` 一一对应）。
- `long locate` 的 PDF 含生成时间，不能与其他运行逐字节比较（测试只比大小）；
  你自己生成的 BAM 与本项目历史产物的字节也不会相同（`@PG` 记录命令行），
  但比对记录本身可用于对账。

## 许可

**GPL-3.0-or-later**。见 `LICENSE`。
