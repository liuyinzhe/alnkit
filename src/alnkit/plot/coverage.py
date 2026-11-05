"""``short coverage`` 的绘图层：从两张表画覆盖广度 / 平均深度（可选基因轨道）。

输入是 :mod:`alnkit.shortread.coverage` 写出的 ``coverage_bins``（逐 bin）与
``coverage_region``（区域汇总）两张表；本模块只做“把表画出来”，不算深度。

图的结构（三块共用一个 1-based 基因组 x 轴）::

    上：Breadth (%)   每个 bin 里有覆盖的碱基占比
    下：Depth (×)     每个 bin 的平均深度，虚线标全区间平均深度
    最下：Gene track  只有给了 --bed 才画，重叠的特征自动分层
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  （backend 必须先设置）
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

from ..core import coords, tables  # noqa: E402
from ..core.errors import RunError, UsageError  # noqa: E402
from ..core.logging import get_logger  # noqa: E402
from ..shortread import coverage as coverage_stats  # noqa: E402

__all__ = [
    "CoverageSettings",
    "Gene",
    "default_bins",
    "parse_bed",
    "render",
    "run",
]

### 配色与字号：整张图共用这一套 ###

BREADTH_COLOR = "#3b7dd8"  # 上图（覆盖广度）
DEPTH_COLOR = "#17a89a"  # 下图（平均深度）
GRID_COLOR = "#e9edf2"  # 水平网格线
SPINE_COLOR = "#c9d1da"  # 坐标轴边框
TEXT_COLOR = "#1d2430"  # 标题与基因名
MUTED_COLOR = "#6b7787"  # 副标题、轴标签
#: BED 里没写颜色时按这个顺序取色（相邻色块必然不同色）
GENE_PALETTE = ["#5b8def", "#e8804f", "#4cb894", "#dfae3c", "#9a8be0", "#d884b4", "#7fb7c9"]
BREADTH_PANEL_HEIGHT = 0.40  # 上图相对下图的高度比
DEPTH_PANEL_HEIGHT = 1.0  # 下图为基准 1.0
GENE_PANEL_HEIGHT = 0.30  # 基因轨道第一行的高度比
GENE_ROW_GROWTH = 0.26  # 每多一行基因，轨道高度比 +0.26（画布同步放大）
GENE_LABEL_INSIDE_RATIO = 0.16  # 色块占区间超过这个比例，名字写在块内
LIGHT_BACKGROUND_LUMA = 0.62  # 色块比这个更亮时，块内文字改用深色

### 分箱默认值 ###


def default_bins(region_length: int, requested_bins: int | None) -> int:
    """没有显式给 ``--bins`` 时，按约 25 bp 一个 bin 决定分箱数（限制在 60..2000）。

    输入示例: default_bins(1416, None)   → 输出示例: 60（1.4 kb 区间取下限）
    输入示例: default_bins(100_000, None) → 输出示例: 2000（约 50 bp/bin 封顶）
    输入示例: default_bins(50, None)     → 输出示例: 50（bin 数不能超过碱基数）
    输入示例: default_bins(1000, 7)      → 输出示例: 7（显式指定就照办）
    """
    if requested_bins:
        return max(1, min(requested_bins, region_length))
    return max(1, min(region_length, max(60, min(2000, region_length // 25))))


### 数值格式化：图与终端摘要用同一套写法 ###


def format_human_length(base_pairs: float) -> str:
    """碱基数转成好读的字符串。

    输入示例: 1416  → 输出示例: "1.4 Kbp"
    输入示例: 87    → 输出示例: "87 bp"
    """
    for unit, scale in (("Mbp", 1_000_000), ("Kbp", 1_000)):
        if base_pairs >= scale:
            return f"{base_pairs / scale:.1f} {unit}"
    return f"{int(base_pairs)} bp"


def format_with_commas(number: float) -> str:
    """整数加千分位。输入示例: 2313 → 输出示例: "2,313"""
    return f"{int(round(number)):,}"


def format_decimal(number: float, digits: int = 1) -> str:
    """保留固定小数位并加千分位。

    输入示例: (231.942, 1) → 输出示例: "231.9"
    输入示例: (12345.6, 2) → 输出示例: "12,345.60"
    """
    return f"{number:,.{digits}f}"


### BED：基因区间 + 可选颜色 ###


@dataclass(frozen=True, slots=True)
class Gene:
    """BED 里的一条特征。

    示例: Gene("chr1", 1389, 1413, "tag", "#e8804f")
          （start/end 是 0-based 半开，color 来自 itemRgb；没写则为 None）
    """

    chrom: str
    start: int  # 0-based
    end: int  # 0-based，不含
    name: str
    color: str | None = None


#: UCSC itemRgb 写法：'255,128,0'
RGB_TRIPLE_PATTERN = re.compile(r"^(\d{1,3})\s*,\s*(\d{1,3})\s*,\s*(\d{1,3})$")
#: '#rgb' 与 '#rrggbb'
HEX_COLOR_PATTERN = re.compile(r"^#(?:[0-9a-fA-F]{3}|[0-9a-fA-F]{6})$")


def parse_bed_color(color_text: str | None) -> str | None:
    """把 BED 里的颜色写法转成 ``#rrggbb``；不是颜色则返回 None。

    输入示例: parse_bed_color("255,128,0") → 输出示例: "#ff8000"（UCSC itemRgb）
    输入示例: parse_bed_color("#1e90ff")  → 输出示例: "#1e90ff"
    输入示例: parse_bed_color("tomato")   → 输出示例: None（只认数字与 hex）
    """
    if not color_text:
        return None
    color_text = color_text.strip()
    if color_text in ("", ".", "0", "-1"):
        return None
    rgb_triple = RGB_TRIPLE_PATTERN.match(color_text)
    if rgb_triple:
        return "#" + "".join(f"{min(255, int(value)):02x}" for value in rgb_triple.groups())
    if HEX_COLOR_PATTERN.match(color_text):
        return color_text
    return None


def parse_bed(bed_path: str | Path) -> list[Gene]:
    """读 BED 文件，返回特征列表（按文件顺序）。

    输入示例（制表符分隔，第 9 列是 itemRgb）::

        # chrom  start  end   name   ...  itemRgb
        chr1     0      1389  geneA  ...  59,141,239
        chr1     1389   1413  tag    ...  232,128,79

    输出示例::

        [Gene("chr1", 0, 1389, "geneA", "#3b8def"),
         Gene("chr1", 1389, 1413, "tag", "#e8804f")]

    颜色优先级：第 9 列 itemRgb > 第 5 列（以 '#' 开头时）；都没有则为 None。
    注释行、track/browser 行、列数不足或坐标不是数字的行会被跳过。
    """
    genes: list[Gene] = []
    with open(bed_path, encoding="utf-8") as bed_file:
        for line in bed_file:
            line = line.rstrip("\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                fields = line.split()
            if len(fields) < 4:
                continue
            try:
                gene_start, gene_end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            if gene_end <= gene_start:  # 坐标非法（start >= end）的行直接跳过
                continue
            color = parse_bed_color(fields[8]) if len(fields) >= 9 else None
            if color is None and len(fields) >= 5 and fields[4].startswith("#"):
                color = parse_bed_color(fields[4])
            genes.append(Gene(fields[0], gene_start, gene_end, fields[3], color))
    return genes


def layout_genes(genes: Sequence[Gene]) -> tuple[list[tuple[Gene, int]], int]:
    """给每条基因分配一个行号，重叠的基因不同行。

    输入示例: [Gene(…, 0, 100, "A"), Gene(…, 50, 150, "B"), Gene(…, 200, 260, "C")]
    输出示例: ([(A, 0), (B, 1), (C, 0)], 2)（A 与 B 重叠，C 回到第 0 行）
    """
    placed_genes: list[tuple[Gene, int]] = []
    row_intervals: list[list[tuple[int, int]]] = []
    for gene in sorted(genes, key=lambda gene: (gene.start, gene.end)):
        for row_index, row in enumerate(row_intervals):
            if all(gene.start >= end or gene.end <= start for start, end in row):
                row.append((gene.start, gene.end))
                placed_genes.append((gene, row_index))
                break
        else:  # 和所有已有行都重叠 -> 新开一行
            row_intervals.append([(gene.start, gene.end)])
            placed_genes.append((gene, len(row_intervals) - 1))
    return placed_genes, len(row_intervals)


def resolve_gene_color(color_from_bed: str | None, gene_index: int) -> str:
    """决定一条基因用什么颜色：BED 给了合法颜色就用它，否则用调色板。

    输入示例: resolve_gene_color("#ff8000", 0) → 输出示例: "#ff8000"
    输入示例: resolve_gene_color("bogus", 2)   → 输出示例: "#4cb894"（回落到调色板）
    """
    if color_from_bed:
        import matplotlib.colors as matplotlib_colors

        if matplotlib_colors.is_color_like(color_from_bed):
            return color_from_bed
    return GENE_PALETTE[gene_index % len(GENE_PALETTE)]


def resolve_label_color(background_color: str) -> str:
    """色块里的文字用什么颜色：深色块用白字，浅色块用深灰字。

    输入示例: resolve_label_color("#3b8def") → 输出示例: "white"
    输入示例: resolve_label_color("#ffffff") → 输出示例: "#1d2430"
    """
    import matplotlib.colors as matplotlib_colors

    red, green, blue = matplotlib_colors.to_rgb(background_color)
    luma = 0.2126 * red + 0.7152 * green + 0.0722 * blue
    return TEXT_COLOR if luma > LIGHT_BACKGROUND_LUMA else "white"


### 主绘图函数 ###


def bin_arrays(bin_rows: Sequence[dict[str, Any]]) -> tuple[list[int], list[float], list[float]]:
    """从 ``coverage_bins`` 的行算出画图需要的两个序列与横轴边界。

    输入示例: 两个 bin：宽 25/覆盖 25/深度和 50，宽 25/覆盖 0/深度和 0，起点 1
    输出示例: ([1, 26], [100.0, 0.0], [2.0, 0.0])   # (边界, 广度%, 平均深度)
    """
    if not bin_rows:
        raise RunError("coverage_bins 表里没有任何 bin")
    boundaries: list[int] = []
    breadth_percent: list[float] = []
    mean_depth: list[float] = []
    for index, row in enumerate(bin_rows):
        bin_start = int(row["bin_start_1based"])
        bin_end = int(row["bin_end_1based"])
        bin_width = int(row["bin_width_bp"])
        if index == 0:
            boundaries.append(bin_start)
        elif bin_start != boundaries[-1] + 1:  # 分箱必须首尾相接，否则图画歪了也不报错
            raise RunError(
                f"coverage_bins 的第 {index + 1} 个 bin 与上一个不相接："
                f"{bin_start} 接不上 {boundaries[-1]}"
            )
        boundaries.append(bin_end)
        breadth_percent.append(100.0 * int(row["covered_bases"]) / bin_width if bin_width else 0.0)
        mean_depth.append(int(row["depth_sum"]) / bin_width if bin_width else 0.0)
    return boundaries, breadth_percent, mean_depth


def render(
    region_row: dict[str, Any],
    bin_rows: Sequence[dict[str, Any]],
    output_path: str | Path,
    genes: Sequence[Gene] | None = None,
    *,
    title: str | None = None,
    width: float = 10.0,
    dpi: int = 200,
) -> str:
    """把两张表画成一张 PNG，返回写出的文件路径。

    输入示例::
        region_row  = {"chrom": "chr1", "region_start_1based": 1, "region_end_1based": 1416,
                       "n_reads_selected": 2313, "percent_covered": 100.0, "mean_depth": 231.9,
                       "mean_mapq": 60.0, "region_length_bp": 1416, "n_bins": 60, …}
        bin_rows    = coverage_bins 的行（bin_arrays 能接受的那种）
        genes       = parse_bed("genes.bed")        # 不给则不画基因轨道

    输出示例: "out/coverage.png"（上图广度、下图深度，最下基因轨道）
    """
    boundaries, breadth_percent, mean_depth = bin_arrays(bin_rows)
    region_name = str(region_row["chrom"])
    region_start = int(region_row["region_start_1based"])
    region_end = int(region_row["region_end_1based"])
    region_length = int(region_row["region_length_bp"])
    region_mean_depth = float(region_row["mean_depth"])
    region_start_0based = coords.to_0based(region_start)
    selected_reads = int(region_row["n_reads_selected"])

    # 只画落在当前区间内的基因，坐标裁到区间边界（复制，不改调用方的 Gene 对象）
    visible_genes: list[tuple[Gene, int]] = []  # (裁剪后的基因, 在原始 BED 里的序号)
    for bed_index, gene in enumerate(genes or []):
        if gene.chrom != region_name or gene.end <= region_start_0based or gene.start >= region_end:
            continue
        clipped = replace(
            gene,
            start=max(gene.start, region_start_0based),
            end=min(gene.end, region_end),
        )
        visible_genes.append((clipped, bed_index))
    color_index_by_gene = {id(clipped): bed_index for clipped, bed_index in visible_genes}
    placed_genes, gene_row_count = layout_genes([clipped for clipped, _ in visible_genes])

    ### 画布高度：基因轨道占多少行就多留多少高度 ###
    panel_height_ratios = [BREADTH_PANEL_HEIGHT, DEPTH_PANEL_HEIGHT]
    extra_rows = max(0, gene_row_count - 1)
    gene_panel_height = GENE_PANEL_HEIGHT + GENE_ROW_GROWTH * extra_rows
    figure_height = 4.4
    if visible_genes:
        panel_height_ratios.append(gene_panel_height)
        # 画布按同一个比例放大，保证每多一行的高度与轨道里多出的高度一致
        figure_height += 0.5 * (gene_panel_height / GENE_PANEL_HEIGHT)

    plt.rcParams.update(
        {
            "font.size": 9,
            "axes.edgecolor": SPINE_COLOR,
            "axes.labelcolor": MUTED_COLOR,
            "text.color": TEXT_COLOR,
            "xtick.color": MUTED_COLOR,
            "ytick.color": MUTED_COLOR,
        }
    )
    figure = plt.figure(figsize=(width, figure_height), constrained_layout=True)
    gridspec = figure.add_gridspec(len(panel_height_ratios), 1, height_ratios=panel_height_ratios)

    breadth_axis = figure.add_subplot(gridspec[0])
    depth_axis = figure.add_subplot(gridspec[1], sharex=breadth_axis)
    gene_axis = figure.add_subplot(gridspec[2], sharex=breadth_axis) if visible_genes else None

    ### 上图：覆盖广度 %（每个 bin 里被覆盖的碱基占比） ###
    breadth_axis.stairs(
        breadth_percent,
        boundaries,
        baseline=0,
        fill=True,
        facecolor=matplotlib.colors.to_rgba(BREADTH_COLOR, 0.10),
        edgecolor=BREADTH_COLOR,
        linewidth=1.3,
    )
    breadth_axis.set_ylim(0, 104)
    breadth_axis.set_yticks([0, 100])
    breadth_axis.set_ylabel(
        "Breadth\n(%)", rotation=0, ha="right", va="center", fontsize=8, labelpad=6
    )

    ### 下图：平均深度 ×（每个 bin 的深度均值），虚线标全区间平均深度 ###
    depth_axis.stairs(
        mean_depth,
        boundaries,
        baseline=0,
        fill=True,
        facecolor=matplotlib.colors.to_rgba(DEPTH_COLOR, 0.26),
        edgecolor=DEPTH_COLOR,
        linewidth=1.0,
    )
    tallest_bin_depth = max(mean_depth) if mean_depth else 0.0
    depth_axis.set_ylim(0, max(1.0, tallest_bin_depth * 1.15))
    depth_axis.yaxis.set_major_locator(MaxNLocator(3, min_n_ticks=2))
    depth_axis.set_ylabel("Depth\n(×)", rotation=0, ha="right", va="center", fontsize=8, labelpad=6)
    if tallest_bin_depth > 0:
        depth_axis.axhline(
            region_mean_depth, color=MUTED_COLOR, linewidth=0.7, linestyle=(0, (4, 3))
        )
        depth_axis.annotate(
            f"mean {region_mean_depth:,.1f}×",
            xy=(region_end, region_mean_depth),
            xytext=(-4, 3),
            textcoords="offset points",
            ha="right",
            va="bottom",
            fontsize=7,
            color=MUTED_COLOR,
            zorder=5,
        )

    ### 基因轨道：一条基因一根色块，窄的写在块上方 ###
    if gene_axis is not None:
        for gene, row_index in placed_genes:
            gene_index = color_index_by_gene[id(gene)]
            row_y = gene_row_count - 1 - row_index
            feature_start = coords.to_1based(gene.start)
            feature_span = gene.end - gene.start
            bar_color = resolve_gene_color(gene.color, gene_index)
            label_color = resolve_label_color(bar_color)
            gene_axis.barh(
                row_y,
                feature_span,
                left=feature_start,
                height=0.46,
                color=bar_color,
                linewidth=0.6,
                zorder=3,
                # 浅色条带用深色描边，否则在白底上看不见
                edgecolor="white" if label_color == "white" else SPINE_COLOR,
            )
            if feature_span > region_length * GENE_LABEL_INSIDE_RATIO:
                gene_axis.text(
                    feature_start + feature_span / 2,
                    row_y,
                    gene.name,
                    ha="center",
                    va="center",
                    fontsize=7.5,
                    fontweight="bold",
                    color=label_color,
                    zorder=4,
                    clip_on=True,
                )
            else:  # 太窄：名字放到块上方，靠右的改为右对齐
                anchor_right = feature_start + feature_span / 2 > region_start + region_length * 0.6
                gene_axis.text(
                    feature_start + feature_span if anchor_right else feature_start,
                    row_y + 0.32,
                    gene.name,
                    ha="right" if anchor_right else "left",
                    va="bottom",
                    fontsize=7,
                    color=TEXT_COLOR,
                    zorder=4,
                    clip_on=False,
                )
        gene_axis.set_ylim(-0.42, gene_row_count - 0.30)
        gene_axis.set_yticks([])
        for spine_name in ("top", "right", "left", "bottom"):
            gene_axis.spines[spine_name].set_visible(False)
        gene_axis.tick_params(axis="x", length=0, labelsize=8.5)

    ### 坐标轴美化：只留左/下边框，浅色水平网格 ###
    for axis in (breadth_axis, depth_axis):
        axis.set_axisbelow(True)
        axis.grid(axis="y", color=GRID_COLOR, linewidth=0.7)
        axis.tick_params(axis="y", length=0)
        for spine_name in ("top", "right"):
            axis.spines[spine_name].set_visible(False)
    depth_axis.set_xlim(region_start, region_end)
    for axis in (breadth_axis, depth_axis):
        axis.tick_params(axis="x", length=0, labelbottom=False)

    # x 轴刻度只画在最下面那块上
    bottom_axis = gene_axis if gene_axis is not None else depth_axis
    bottom_axis.xaxis.set_major_locator(MaxNLocator(8, integer=True, steps=[1, 2, 2.5, 5, 10]))
    bottom_axis.xaxis.set_major_formatter(
        FuncFormatter(lambda value, _pos: format_with_commas(value))
    )
    bottom_axis.tick_params(axis="x", length=0, labelbottom=True, pad=4)
    bottom_axis.set_xlabel("Position (bp)", fontsize=8, labelpad=6)

    ### 标题与副标题（标题下的一行摘要） ###
    figure_title = title or (
        f"{region_name}  ·  {format_with_commas(region_start)}–{format_with_commas(region_end)} bp"
    )
    breadth_axis.set_title(
        figure_title, loc="left", fontsize=12, fontweight="bold", color=TEXT_COLOR, pad=20
    )
    subtitle = (
        f"{format_human_length(region_length)}"
        f"  ·  {format_with_commas(selected_reads)} reads"
        f"  ·  {float(region_row['percent_covered']):.1f}% covered"
        f"  ·  mean depth {region_mean_depth:,.1f}×"
        f"  ·  mapQ {float(region_row['mean_mapq']):.1f}"
        f"  ·  {len(breadth_percent)} bins"
    )
    breadth_axis.text(
        0,
        1.015,
        subtitle,
        transform=breadth_axis.transAxes,
        fontsize=8.5,
        color=MUTED_COLOR,
        va="bottom",
    )

    ### 区间内没有 read 时给一句提示，图仍然写出来 ###
    if selected_reads == 0:
        breadth_axis.text(
            0.5,
            0.5,
            "no reads in this region",
            transform=breadth_axis.transAxes,
            ha="center",
            va="center",
            fontsize=11,
            color=MUTED_COLOR,
        )

    figure.savefig(output_path, dpi=dpi, facecolor="white")
    plt.close(figure)
    return str(output_path)


### 区间解析与参考列表 ###


def resolve_region(alignment_file, region_text: str) -> coverage_stats.Region:
    """把 ``--region`` 字符串变成 :class:`~alnkit.shortread.coverage.Region`。

    输入示例: resolve_region(bam, "chr1:101-200") → 输出示例: Region("chr1", 100, 200)
    输入示例: resolve_region(bam, "chr1")         → 输出示例: Region("chr1", 0, <该序列长度>)

    异常: 参考名不在 BAM 里时抛 ``UsageError`` 并列出可用名（最多 5 个）。
    """
    region_text = region_text.strip()
    reference_lengths = dict(zip(alignment_file.references, alignment_file.lengths, strict=True))
    if not region_text:
        raise UsageError("--region 必填，形如 chr1 或 chr1:1-2000（1-based 闭区间）")
    if region_text in reference_lengths:  # 整个 contig（名字本身带 ':' 也走这条）
        return coverage_stats.Region(region_text, 0, reference_lengths[region_text])
    try:
        chrom, start, end = coords.parse_region_1based(region_text)
    except ValueError as error:
        raise UsageError(f"{error}（--region 写成 chr1 或 chr1:1-2000）") from None
    if chrom not in reference_lengths:
        names = ", ".join(list(reference_lengths)[:5])
        raise UsageError(f"{chrom!r} 不是该 BAM 的参考序列；可用：{names} …")
    end = min(end, reference_lengths[chrom])
    if end <= start:
        raise UsageError(f"--region {region_text} 是空区间（超出该序列长度？）")
    return coverage_stats.Region(chrom, start, end)


def default_name(bam_path: str | Path, region: coverage_stats.Region, reference_length: int) -> str:
    """没给 ``--name`` 时的默认前缀：``<bam 主名>.<contig>[_start-end]``。

    输入示例: default_name("sampleA.markdup.bam", Region("chr1", 0, 1416), 1416)
    输出示例: "sampleA.chr1"
    输入示例: default_name("sampleA.markdup.bam", Region("chr1", 99, 200), 170805979)
    输出示例: "sampleA.chr1_100-200"
    """
    bam_stem = Path(bam_path).name.split(".")[0] or "coverage"
    covers_whole_reference = region.start == 0 and region.end == reference_length
    label = (
        region.name
        if covers_whole_reference
        else f"{region.name}_{coords.to_1based(region.start)}-{region.end}"
    )
    return f"{bam_stem}.{re.sub(r'[^A-Za-z0-9._-]', '_', label)}"


### 编排：BAM → 两张表 → PNG ###


@dataclass(frozen=True, slots=True)
class CoverageSettings:
    """一次 ``short coverage`` 的全部输入。

    示例::

        CoverageSettings(bam="reads.bam", region="chr1:1-2000", out_dir="out")
    """

    bam: str
    region: str = ""
    out_dir: str = "."
    name: str = ""  # 空 = 按 default_name 生成
    bins: int | None = None
    bed: str | None = None
    title: str | None = None
    dpi: int = 200
    width: float = 10.0
    fail_flags: str = "UNMAP,SECONDARY,QCFAIL,DUP"
    require_flags: str = ""
    min_mapq: int = 0
    min_len: int = 0
    min_depth: int = 1
    max_depth: int = 1_000_000
    bases_table: str | None = None
    list_only: bool = False


def _open_alignment_file(bam_path: str | Path):
    """打开 BAM/CRAM；缺文件或没有索引时给可行动的报错，而不是 pysam 的 traceback。"""
    import pysam

    path = Path(bam_path)
    if not path.exists():
        raise UsageError(f"--bam 不存在：{path}")
    try:
        return pysam.AlignmentFile(str(path))
    except (OSError, ValueError) as error:
        raise RunError(f"打不开 BAM：{path}（{error}）") from error


def run(settings: CoverageSettings) -> dict[str, Any]:
    """打开 BAM → 写两张表 → 读回表 → 出图；``list_only`` 时只列参考序列。

    输入示例: CoverageSettings(bam="reads.bam", region="chr1:1-2000", bins=4)
    输出示例: {"png": "out/sampleA.chr1.png", "region_tsv": "…", "bins_tsv": "…",
               "region": "chr1:1-2000", "reads": 12, "percent_covered": 75.0,
               "mean_depth": 1.25, "mean_mapq": 60.0, "mean_baseq": 35.0}

    为什么要“写表再读表”：表是接口。图只允许看到表里的东西，读回来一次就能保证
    表本身足够画出这张图（少列、错列会在这里立刻报错，而不是画出一张空图）。
    """
    alignment_file = _open_alignment_file(settings.bam)
    with alignment_file:
        if settings.list_only:
            return {"references": _reference_rows(alignment_file)}

        if settings.width <= 0 or settings.dpi <= 0:
            raise UsageError("--width 与 --dpi 必须为正数")
        if settings.bins is not None and settings.bins < 1:
            raise UsageError("--bins 必须 >= 1")
        if "/" in settings.name or "\\" in settings.name:
            raise UsageError("--name 只能是文件名前缀，不能带目录（用 --out-dir）")

        region = resolve_region(alignment_file, settings.region)
        filters = coverage_stats.Filters(
            fail_flags=coverage_stats.parse_flags(settings.fail_flags),
            require_flags=coverage_stats.parse_flags(settings.require_flags),
            min_mapq=settings.min_mapq,
            min_len=settings.min_len,
            min_depth=settings.min_depth,
            max_depth=settings.max_depth,
        )
        bins = default_bins(region.length, settings.bins)
        try:
            result = coverage_stats.compute_coverage(alignment_file, region, filters, bins)
        except ValueError as error:  # 没有索引时 fetch 会报这个
            raise UsageError(
                f"{settings.bam} 无法按区间取数：{error}（先跑 samtools index）"
            ) from error
        reference_length = alignment_file.get_reference_length(region.name)

    genes = _load_genes(settings.bed, region)

    output_directory = Path(settings.out_dir)
    output_directory.mkdir(parents=True, exist_ok=True)
    name = settings.name or default_name(settings.bam, region, reference_length)
    sample = Path(settings.bam).name.split(".")[0] or "sample"
    region_path = output_directory / f"{name}.region.tsv"
    bins_path = output_directory / f"{name}.bins.tsv"
    png_path = output_directory / f"{name}.png"
    comments = (
        f"sample={sample}",
        f"region={coords.format_region_1based(region.name, region.start, region.end)}",
        "filters="
        f"fail:{settings.fail_flags or 'none'} "
        f"require:{settings.require_flags or 'none'} "
        f"min_mapq:{settings.min_mapq} min_len:{settings.min_len} "
        f"min_depth:{settings.min_depth} max_depth:{settings.max_depth}",
        "producer=alnkit short coverage",
    )
    written = coverage_stats.write_coverage_tables(
        result,
        sample,
        region_path=region_path,
        bins_path=bins_path,
        bases_path=settings.bases_table,
        comments=comments,
    )

    # 从这里开始只认表：图由读回来的行画出
    region_header, region_rows = tables.read_tsv(region_path)
    bins_header, bins_rows = tables.read_tsv(bins_path)
    for header, wanted, path in (
        (region_header, tables.columns_of("coverage_region"), region_path),
        (bins_header, tables.columns_of("coverage_bins"), bins_path),
    ):
        if tuple(header) != tuple(wanted):
            raise RunError(f"{path} 的列与注册表不一致：{header}")
    region_row = region_rows[0]
    render(
        region_row,
        bins_rows,
        png_path,
        genes,
        title=settings.title,
        width=settings.width,
        dpi=settings.dpi,
    )
    return {
        "png": str(png_path),
        "region_tsv": str(region_path),
        "bins_tsv": str(bins_path),
        "reference_length": reference_length,
        "region": f"{region.name}:{coords.to_1based(region.start)}-{region.end}",
        "region_name": region.name,
        "region_start_1based": coords.to_1based(region.start),
        "region_end": region.end,
        "n_bins": written["bins"],
        "reads": result.selected_read_count,
        "percent_covered": result.percent_covered,
        "mean_depth": result.mean_coverage,
        "mean_mapq": result.mean_mapq,
        "mean_baseq": result.mean_baseq,
    }


def _reference_rows(alignment_file) -> list[tuple[str, int, int]]:
    """``--list`` 的行：参考序列名、长度、比对上来的 reads 数（按 reads 数降序）。

    输入示例: 一个有两万条 reads 落在 chr1 的 BAM
    输出示例: [("chr1", 170805979, 3769), ("insertA", 1416, 2691)]
    """
    try:
        statistics = sorted(alignment_file.get_index_statistics(), key=lambda row: -row.mapped)
    except ValueError as error:  # 未建索引
        raise UsageError(f"--list 需要索引（{error}）；先跑 samtools index") from None
    return [
        (row.contig, alignment_file.get_reference_length(row.contig), row.mapped)
        for row in statistics
        if row.mapped
    ]


def _load_genes(bed_path: str | None, region: coverage_stats.Region) -> list[Gene] | None:
    """读 BED（给了才读）；没有落在该区间的特征时只警告一次。"""
    if not bed_path:
        return None
    path = Path(bed_path)
    if not path.exists():
        raise UsageError(f"--bed 不存在：{path}")
    try:
        genes = parse_bed(path)
    except (OSError, UnicodeDecodeError) as error:
        raise RunError(f"读不了 BED {path}：{error}") from error
    if genes and not any(gene.chrom == region.name for gene in genes):
        get_logger(__name__).warning("BED 里没有 %s 上的特征，基因轨道为空", region.name)
    return genes
