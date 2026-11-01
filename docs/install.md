# 安装

## 从 PyPI

```bash
pip install alnkit            # 分析链路
pip install "alnkit[plot]"    # + 绘图（matplotlib）
```

要求 Python >= 3.10。

- **必需**：`pysam`（读 BAM/FASTA，并**内置 samtools 1.23**，参考索引/排序/建索引/
  `markdup` 都由它完成——所以 `samtools` 不是必需的外部命令）；
- **可选**：`matplotlib`（`alnkit[plot]`，只有绘图子命令需要）；`minimap2`（比对/overlap
  要用，目前没有纯 Python 等价替代，见下）。
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
| `samtools` | **可选**：参考索引、BAM 排序/去重。`pysam` 已内置同名功能（1.23），alnkit 默认用 pysam 调它们，所以不装也能跑 | >= 1.16 | `conda install -c bioconda samtools` |
| `bwa` | 二代比对（可选，用户自己跑） | — | `conda install -c bioconda bwa` |

装完自检：

```bash
alnkit doctor
```

它会打印每个组件的版本；必需项缺失时退出码为 1，可选项缺失只提示。

## 不装外部命令会怎样

分析链路里所有命令都接受**用户自带的 BAM / PAF**，此时不会调用 `minimap2`/`samtools`。
只有“把比对也交给 alnkit 跑”的便利命令才需要它们。
