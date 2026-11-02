# 安装

## 从 PyPI

```bash
pip install alnkit            # 分析链路
pip install "alnkit[plot]"    # + 绘图（matplotlib）
```

要求 Python >= 3.10。

- **必需**：`pysam`（读 BAM/FASTA；自带 samtools 能力，但 0.3.0 起包内不再做排序/建索引）；
- **长读长链路必需**：`minimap2`（`long classify` 的 overlap、`long locate` 的锚点比对、
  `plot ribbon` 的内部比对、`long assemble` 默认的 `assembly_support.tsv` 都要它；
  目前没有纯 Python 等价替代）；
- **可选**：`matplotlib`（`alnkit[plot]`，只有绘图子命令需要）。
- 不要装 `mappy` 也能用：它只是 minimap2 的 Python 绑定，只返回比对命中，不产出 BAM/PAF，
  无法替代本包对 `minimap2` 命令行输出的依赖。

## 从源码（开发）

```bash
git clone git@github.com:liuyinzhe/alnkit.git
cd alnkit
pip install -e ".[dev,docs]"
python -m pytest
```

## 外部命令（不在 wheel 里）

| 命令 | 用途 | 建议版本 | 装法 |
|---|---|---|---|
| `minimap2` | 比对（长读长 / overlap / 丝带图内部比对） | >= 2.24 | `conda install -c bioconda minimap2` |
| `samtools` | **可选**：只有 `plot tracks --bam-mode subset` 会调用外部 `samtools`；准备上游文件（`faidx`/`sort`/`index`）时你自己需要 | >= 1.16 | `conda install -c bioconda samtools` |
| `bwa` | 二代比对（可选，用户自己跑） | — | `conda install -c bioconda bwa` |

装完自检：

```bash
alnkit doctor
```

它会打印每个组件的版本；必需项缺失时退出码为 1，可选项缺失只提示。

## 不装外部命令会怎样

`long locate`、`long assemble`（默认开启的支撑度报告，可用 `--no-support-report` 关掉）、
`plot ribbon` 的「内部比对」模式都会调用 `minimap2`——它们要吃 `minimap2` 的 SAM/PAF 输出，
没有替代实现。缺 `minimap2` 时：`long assemble` 只警告并跳过该报告（其余产物照常），
`long locate` 与 `plot ribbon` 直接报错（退出码 1）。
其余命令只吃你自带的 BAM/PAF，不需要任何外部比对器；
`samtools` 只有 `plot tracks --bam-mode subset` 会用。
