# Changelog

格式遵循 [Keep a Changelog](https://keepachangelog.com/zh-CN/1.1.0/)，版本号遵循 SemVer。
条目按功能 ID（`A*` 三代 / `B*` 二代 / `C*` 绘图 / `D*` 公共层）分组，ID 定义见规划文档。

## [0.4.0] - 2026-11-02

### Added

- **A6 `long assemble` 组装支撑度报告**（`assembly_support.tsv`，默认开启）：把类内 reads 用
  `minimap2 -ax map-hifi -Y` 打回 contig，输出 `covered_fraction` / `mean_depth`（×）/
  `depth_p10` / `depth_cliffs` / `ends_covered` / `status`
  （`ok` / `low_coverage` / `chimera_suspect` / `unsupported_ends` / `no_reads`）。
  只解析 SAM stdout，不落 BAM、不需要 samtools。没有 minimap2 时只警告并跳过；
  `--no-support-report` 静默跳过。新增表 `assembly_support`（12 列，已进 `docs/concepts/tables.md`）。
- `long assemble` 新增 `--threads`（默认 4），供支撑度报告的 minimap2 回映使用。
### Removed

- `long assemble` 不再输出 `layout.tsv`（0.1.0 起就在做，本版补记）：它没有下游消费者，
  且其中记录的是“每轮当时的坐标系”，与最终 contig 逐位一致率只有约 28%——属**输出记账**问题，
  贪心延伸的决策本身自洽（把坐标随左延伸平移反而会让序列被悄悄改短）。
  排障需要该表时可直接跑上游 `scripts/overlap_assembly.py`。
- **评估后未采纳**（为避免留下无效开关，实现已整块删除）：接缝绝对支持数门槛
  （接缝处邻居数按构造恒为 1，门槛 2 永不触发、门槛 3 会把 5 个类拆成 642 条 contig）、
  接缝相对深度门槛（0.1/0.2/0.3 会把 5 条 contig 切成 10/16/17 条，断点锚点回收始终 8/8，
  无收益）、DEGAP 式边锚定招募（contig 长 +18–31% 但延长段常落重复区，同样不改善锚点回收）。

## [0.3.0] - 2026-11-02

### Removed

- **移除 `alnkit ref merge` 与 `alnkit align reads`**（以及与它们配套的测试与对账）：
  这两步用现成命令做即可（`cat` + `samtools faidx`、`minimap2 -Y -N 50 | samtools sort`），
  放在包里属于重复封装。README 新增「上游准备」一节，写明命令与两个关键参数的必要性
  （`-Y` 让补充比对用 soft clip，否则断点证据丢失；`-N 50` 保留次要比对用于归属判定）。
  `long classify --help` 也会直接打印这段要求。

### Fixed

- `doctor` 之前把 `bwa` 的报错当成版本号（`bwa --version` 不被支持，会输出
  "unrecognized command"）。现在版本探测在 `--version` 失败时退化成无参调用，
  从帮助文本里取 `Version:` 一行。
- 给 `long|plot|seg` 的子命令加了 `RawDescriptionHelpFormatter`：多行 `--help` 不再被重排成一行。

## [0.2.0] - 2026-11-02

### Security

- **全量脱敏**（0.1.0 发布后发现仓库与包内含内部信息，本版一并清理）：
  - 本机绝对路径 → 改由 `ALNKIT_UPSTREAM_*` 环境变量或 `--source/--bam-dir` 指定；
  - 样本名/细胞系名/项目号 → `sampleA` / `sampleB`；
  - 外源构建体与靶点名 → `target_ref` / `homolog_ref` / `decoy_ref` / `insertA` / `insertB`；
  - 真实参考序列号与坐标 → `chr1`/`chr2` 与 100 万级合成坐标；
  - 测试 fixture 的 read 名与坐标同步替换（黄金用上游工具按新 fixture 重新生成）；
  - 移除作者邮箱；真实序列名只保留在**本地** `manifest.json`（不入库），测试从 manifest 读。

### Changed

- **`samtools` 从必需变成可选**：`align reads` 的排序/建索引改用 `pysam` 内置的 samtools
  （`pysam.sort` / `pysam.index`，0.24 对应 samtools 1.23），比对记录与旧 shell 产物逐字节相同；
  只有 pysam 封装不可用时才回退到 PATH 里的 `samtools`。
- **项目描述收敛**为最短版本：「BAM/PAF 的插入证据工具箱：断点提取 + 独立绘图」
  （英文 `Insertion evidence toolkit for BAM/PAF: breakpoint extraction and standalone plots`）。
- `sites.canonical.tsv` 增加插入归属三列：`query_insert` / `nearest_other_insert` / `other_inserts`
  （一条 read 可能同时命中多个插入序列，"最近"按支持 read 数定义，query 自己不算竞争者）。

## [0.1.1] - 2026-11-01

### Added

- **`<sample>.sites.canonical.tsv`**（`short scan` 新增产物，0.1.0 缺失的那张表）：
  二代位点的机器接口，带 `insert_pos_1based` / `insert_pos_0based` / `insert_boundary_0based`
  三列，恒等式 `insert_pos_0based + 1 == insert_pos_1based == insert_boundary_0based`；
  与 `<sample>.Insert_Position.result.xls` 按 `(chrom, insert_pos_1based)` 逐行对应。

## [0.1.0] - 2026-11-01

### Fixed

- **基础安装可用**：`alnkit` 之前会在 import 阶段无条件导入绘图模块，未装 matplotlib 时
  连 `alnkit --version` 都崩（`pip install alnkit` 只装 pysam）。改为惰性导入，缺依赖时
  给 `pip install "alnkit[plot]"` 的可行动提示 + 退出码 2。
- **反链 overlap 落位差 `aligned_bp-1`**：`placement_from_block` 在 `strand="-"` 时把
  unplaced 起点与 placed 起点配对，反链应为头尾相对（起点对终点）。本数据集 38% 的 overlap
  是反链，所以 06 组装的 contig 与 07 的图都会变；同一修复已提交到上游（`82a50b7`）。
- **合并参考漏一条序列**：输入 FASTA 末尾缺换行时，`ref merge` 的拼接会把 `>下一条` 粘到
  上一条序列末尾，`faidx` 仍然成功、下游静默少一条参考。现在边界处强制补换行。
- `long classify` 的 `--exog` 打错名字不再静默输出 0 条结果，直接报可用序列名。
- `long locate` 的位点表与参考 `.fai` 不匹配时给 `UsageError`（原来是 `KeyError` traceback）。
- `short scan` 两个 BAM 推导出同名样本时提前报错（原来后一个会静默覆盖前一个的输出）。
- `short scan` 缺 `NM` 标签的插入侧比对不再被当成“零错配”（那是最优 tie-break，会让它抢到归属）。

### Added

- 契约层加固（独立复审后修）：strand 空值不再静默当 `+`；可选整数列坏值报 `UsageError`
  而不是裸 `ValueError`；坐标列后缀校验不再能被 `ref_end_bp` 这类名字绕过；
  原子写保留目标文件原权限（`mkstemp` 的 0600 不再泄漏）；`write_legacy_tsv` 也改成原子写
  （失败不覆盖旧文件）；零长区间/`left_base_0based(0)`/千分位位置错误一律报错。
- 黄金回归的**出身校验**：`tests/fixtures/build_private_*.py` 默认要求旧脚本存在
  （`--trust-current` 需显式给），每个黄金目录写 `PROVENANCE.txt`（generator/命令/脚本 sha256），
  `tests/test_golden_provenance.py` 拒绝 `generator=alnkit` 的黄金。
  builder 同时补齐了 A5 需要的三个输入（以前只能手工准备）。
- 私有数据已设 `ALNKIT_TEST_DATA` 时，缺输入从 `skip` 改为 `fail`。

- A4 `alnkit long graph`：read-vs-read overlap 边与连通分量，移植 `overlap_graph.py`；
  三张产物逐字节一致。
- A1 `alnkit ref merge`：把 `01_merge_ref.sh` 改写成 Python（不再依赖 awk）。
  `reference_headers.tsv`（含 awk 的单空格重组语义）、`exogenous_uppercased.fa`、两个 `.fai`
  逐字节一致；2.3 GB 的 `merged_reference.fa` md5 相同。
- A2 `alnkit align reads`：把 `02_align_merged.sh` 改写成外部命令层调用。
  **BAM 不能逐字节比**（minimap2 把命令行含输出路径写进 `@PG`），验收改为比对记录
  `samtools view -F 0x900` 的 md5 一致（实测相同）。
- A7 `alnkit long locate`：50 bp 断点锚点 + 外源 → 丝带图 + 参考坐标还原，移植
  `ribbon_positions.py`。非 PDF 产物（7 张 PNG、`class_assembly_positions.tsv`、
  `contig_map.tsv`、`work/*.paf`、`summary.md`、`tables/*`）与旧脚本逐字节一致。
  两处**有意改动**：① 不再 subprocess 调 vendored `tools/synteny_ribbon.py` 1.8.0，
  改为进程内调用本包 `plot.ribbon`（1.10.1）；audit 的版本影响评估已实测两者对当前
  参数产物相同，本次也复核了 PNG 逐字节一致。② PDF 不参与逐字节对账——matplotlib 在
  PDF 里写生成时间，旧脚本自己连跑两次也不一致，测试只比大小。
- A6 `alnkit long assemble`：三代每类 overlap 贪心延伸组装，移植 `overlap_assembly.py`。
  旧实现用 `for unplaced_name in list(unplaced)` 遍历 `set`，同一输入在不同
  `PYTHONHASHSEED` 下会产出不同 `layout.tsv`（实测三次三种结果）。该问题已**先在上游修**
  （`upstream` commit `20b4dda`，本包 `src/alnkit/longread/assemble.py` 为同构实现），
  修完两边 `diff -r` 整目录逐字节一致（含 `layout.tsv`）。
- core 补齐：`paf.parse_paf`、`fasta.{reverse_complement,write_fasta_records,read_fastq_sequences}`；
  新增 `python -m alnkit` / `python -m alnkit.cli` 入口。
- A5 `alnkit long breakpoints`：三代断点（S##/G##）+ 共享 reads 并查集合并成类 + 类池，
  移植 `breakpoints.py`；`diff -r` 整目录逐字节一致（含 `classes/<C>/*.fastq` 与 `graph_summary.txt`）。
  唯一非内容差异：`class_read_pools.tsv` 的 `reads_fastq` 列是绝对路径，比对时按 `<OUTDIR>` 归一化。
- core 补齐：`tables.read_legacy_tsv`（旧语义读取，行是字符串列表）。
- A3 `alnkit long classify`：三代目标 reads 分类 + 接头事件/位点/区域聚类，
  移植 `03_classify_reads.py`；**整目录**（4 张表 + 4 个 FASTQ 子集 + `regions/`）
  与旧脚本逐字节一致。
- core 补齐：`Segment` 增加 `nm` / `match_columns` / `inserted_bases` / `deleted_bases`；
  `tables.write_legacy_tsv`（旧语义写出）；`coords.merge_intervals`。
- `shortread` 按依赖分成 7 个模块（`model` / `bam` / `interference` / `evidence` /
  `origin` / `merge` / `report`）+ 编排层 `scan`；拆分后 owned / markdup / consensus
  三种口径的黄金回归仍然逐字节一致。
- B1 `alnkit short scan`：二代插入断点证据（SR / PESR / DUP），移植 `parse_bam.py` v3.0.0；
  三张表在 owned / markdup / consensus 三种口径下与旧脚本逐字节一致（黄金回归需私有数据，
  公开 CI 跳过）。`consensus` 策略下 `candidate_insert_chromosomes` 先用后赋值的潜伏
  `NameError` 已**先在上游修**（已在上游同步修复）。
  另有一处**有意差异**：规范表 `segments` 的 `nm` 列在 NM 标签缺失时记 `-`（旧 14 列表记 0，
  该列不在 legacy 表里，故不影响对账）。
- B2 `parse_bam --self-test` 的 30+ 条口径断言并入 `alnkit selftest`。
- C1 `alnkit plot ribbon`：共线性丝带图，移植 `synteny_ribbon.py` v1.10.1；`--paf` 模式下
  `blocks.tsv` 与 PNG 与旧脚本逐字节一致。相对旧版只加了一条输入校验：空标签直接报错
  （旧版会画出一行无名轨道）。
- D2 `alnkit seg table`：BAM → 比对段表；`--columns legacy` 复现历史 14 列（与旧脚本
  `alignment_table.py` 逐字节一致），默认写带基准后缀的规范列。含 CIGAR 简化视图与
  read 轴镜像归正（`core/bam_segments.py`）。
- C2 `alnkit plot tracks`：单 read 比对片段轨道图（BAM / PAF / 段表三选一），
  由 `read_segment_tracks.py`（git 51bf68f）忠实移植；`segments.tsv` 与 PNG 与旧脚本逐字节一致。

- D6 `alnkit doctor`：Python / pysam / matplotlib / minimap2 / samtools / bwa 的可用性与版本。
- D6 `alnkit selftest`：坐标换算、表后缀规则、三套历史段表头的归一断言。
- D1–D3 契约层：`coords`（唯一 ±1 入口）、`tables`（列注册表 + 原子写 + 后缀校验）、
  `segments`（`Segment` 模型 + 别名读取）、`external`（工具探测）。
- 文档骨架：`docs/index.md`、`docs/install.md`、`docs/quickstart.md`、`docs/concepts/coordinates.md`、
  `docs/concepts/tables.md`（由 `core.tables` 生成）。

[Unreleased]: https://github.com/liuyinzhe/alnkit/compare/v0.4.0...HEAD
[0.4.0]: https://github.com/liuyinzhe/alnkit/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/liuyinzhe/alnkit/releases/tag/v0.3.0
[0.2.0]: https://github.com/liuyinzhe/alnkit/releases/tag/v0.2.0
[0.1.1]: https://github.com/liuyinzhe/alnkit/releases/tag/v0.1.1
[0.1.0]: https://github.com/liuyinzhe/alnkit/releases/tag/v0.1.0
