"""``short insertion-plot`` 的绘图层：从表渲染外源插入位点示意图。

输入只有两样东西：**已经算好的** :class:`~alnkit.shortread.insertion.ClipGroup` 列表
（来自 ``insertion_clip_reads`` 表）与一个查看窗口；本模块不碰 BAM、不算一致性序列。

图按"给了几个断点 + 每个断点有几侧 clip"自动分三种::

    单侧        1 个断点、只有一侧 clip  → [未知区][彩色块] 贴在断点一侧
    同断点双侧  1 个断点、两侧都有 clip  → [左块][未知（合并）][右块] + 两条连线收到断点
    多断点      2–3 个断点              → 每个断点各自 [未知区][块] / [块][未知区]

配色固定：左 clip = 蓝、右 clip = 紫、未知区 = 灰斜线、断点 = 红。
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402  （backend 必须先设置）
from matplotlib.patches import Rectangle  # noqa: E402
from matplotlib.ticker import FuncFormatter  # noqa: E402

from ..core import coords  # noqa: E402
from ..core.errors import RunError, UsageError  # noqa: E402
from ..shortread import insertion  # noqa: E402

__all__ = [
    "InsertionSettings",
    "draw_multiple_breakpoints_figure",
    "draw_single_breakpoint_figure",
    "plan_multi_breakpoint_blocks",
    "render",
    "run",
]

### 配色与尺寸（改这里就能整体换风格） ###

LEFT_CLIP_COLOR = "#2980b9"
LEFT_CLIP_EDGE = "#1a5276"
RIGHT_CLIP_COLOR = "#8e44ad"
RIGHT_CLIP_EDGE = "#4a235a"
#: 未知区（无证据）与断点的颜色
UNKNOWN_COLOR = "#7f8c8d"
BREAKPOINT_COLOR = "#c0392b"
REFERENCE_LINE_COLOR = "#566573"
#: read 上虚线框标出 clip 余段（不属于任何一方的中性色）
CLIP_OUTLINE_COLOR = "#e67e22"

#: 未知区条形的宽度 = max(该侧最大 clip 长度, 下限) × 该倍数（示意，不是真实 bp）
UNKNOWN_WIDTH_FACTOR = 2.6
UNKNOWN_WIDTH_MIN_CLIP = 100
#: 相邻的不同断点之间至少留出的空白（窗宽的 3%，下限 30 bp）
BLOCK_GAP_FRACTION = 0.03
BLOCK_GAP_MIN_BP = 30
#: 嵌入块内的标签放不下时，块宽低于该值就改成块外标注
NARROW_BLOCK_BP = 90


def side_colors(clip_side: str) -> tuple[str, str]:
    """该侧使用的（填充色, 边框色）。

    输入示例: clip_side="left"   →  输出示例: ("#2980b9", "#1a5276")
    输入示例: clip_side="right"  →  输出示例: ("#8e44ad", "#4a235a")
    """
    return (
        (LEFT_CLIP_COLOR, LEFT_CLIP_EDGE)
        if clip_side == "left"
        else (RIGHT_CLIP_COLOR, RIGHT_CLIP_EDGE)
    )


### 坐标轴与画布工具 ###


def format_position_axis(axes_object) -> None:
    """把 x 轴刻度改成"千分位整数"，避免 matplotlib 的 +7.35e7 偏移记号。

    输入示例: 轴范围 1000000–1000900
    输出示例: 刻度显示 "1,000,000"、"1,000,100" …
    """
    axes_object.xaxis.set_major_formatter(FuncFormatter(lambda value, _pos: f"{int(value):,}"))


def expand_axis_limits(
    span_start: int, span_end: int, block_edges: list[int]
) -> tuple[float, float]:
    """把 x 轴范围扩到能完整容纳所有块，避免块被轴边界裁掉。

    输入示例: span_start=299990, span_end=300010, block_edges=[299784, 300000]
    输出示例: (299734, 300060)（左侧扩到块起点外 50 bp，右侧扩到窗口终点外 50 bp）
    """
    padding = max(50, int((span_end - span_start) * 0.04))
    return min([span_start, *block_edges]) - padding, max([span_end, *block_edges]) + padding


def save_figure(figure, output_path: str | Path, dpi: int) -> None:
    """保存并关闭画布（tight 裁剪保证底部说明与公式不被切掉）。

    输入示例: figure=<Figure>, output_path=out/demo.png, dpi=150
    输出示例: 生成 out/demo.png（像素宽 = 图宽英寸 × 150）

    写不了时转成 ``RunError``（退出码 1），并保证画布无论如何都被关掉。
    """
    try:
        figure.savefig(output_path, dpi=dpi, bbox_inches="tight", pad_inches=0.25)
    except OSError as error:
        raise RunError(f"写不了 {output_path}：{error}") from error
    finally:
        plt.close(figure)


def build_figure(*, figure_width: float = 13.0, figure_height: float = 8.8):
    """建一张"上插入轨道 / 下 reads 轨道"的双面板画布。

    输入示例: 默认尺寸
    输出示例: (figure, insert_axis, reads_axis)，两面板上下排列、底部留出说明空间
    """
    figure = plt.figure(figsize=(figure_width, figure_height))
    grid = figure.add_gridspec(
        2,
        1,
        height_ratios=[1.0, 3.4],
        hspace=0.42,
        left=0.06,
        right=0.98,
        top=0.9,
        bottom=0.22,
    )
    insert_axis = figure.add_subplot(grid[0])
    reads_axis = figure.add_subplot(grid[1])
    insert_axis.set_ylim(-0.6, 8.8)
    insert_axis.set_yticks([])
    insert_axis.spines[["left", "right", "top"]].set_visible(False)
    format_position_axis(insert_axis)
    format_position_axis(reads_axis)
    insert_axis.tick_params(axis="x", labelsize=8)
    reads_axis.tick_params(axis="x", labelsize=8)
    return figure, insert_axis, reads_axis


### 绘图元素：插入块 / 未知区 / 断点线 ###


def draw_insert_block(
    axes_object,
    start_position: int,
    end_position: int,
    clip_group: insertion.ClipGroup,
    *,
    top: float = 3.0,
    height: float = 3.8,
) -> None:
    """在插入轨道上画一个彩色 clip 块，块内居中写长度（如 "60bp"）。

    输入示例: start_position=99940, end_position=100000, group.consensus 长 60 bp
    输出示例: 一块蓝色矩形，内部白字 "60bp"；块太窄时标签移到块外侧
    """
    fill_color, edge_color = side_colors(clip_group.side)
    axes_object.add_patch(
        Rectangle(
            (start_position, top),
            end_position - start_position,
            height,
            facecolor=fill_color,
            edgecolor=edge_color,
            alpha=0.95,
        )
    )
    label = f"{len(clip_group.consensus)}bp" if clip_group.consensus else "n/a"
    if end_position - start_position >= NARROW_BLOCK_BP:
        axes_object.text(
            (start_position + end_position) / 2,
            top + height / 2,
            label,
            ha="center",
            va="center",
            fontsize=8.5,
            color="white",
        )
        return
    # 块太窄放不下白字：改到块外侧用深色字，靠边的块不会压到中间的未知区
    axes_object.text(
        start_position - 6,
        top + height / 2,
        label,
        ha="right",
        va="center",
        fontsize=7.5,
        color=edge_color,
        bbox={
            "boxstyle": "round,pad=0.25",
            "fc": "white",
            "ec": edge_color,
            "linewidth": 0.5,
            "alpha": 0.95,
        },
    )


def draw_unknown_block(
    axes_object, start_position: int, end_position: int, *, top: float = 3.0, height: float = 3.8
) -> None:
    """画灰斜线矩形表示"没有 read 证据"的未知区。

    输入示例: start_position=99784, end_position=99940
    输出示例: 一条斜线填充的矩形（与 clip 块同一水平线）
    """
    axes_object.add_patch(
        Rectangle(
            (start_position, top),
            end_position - start_position,
            height,
            facecolor="none",
            edgecolor=UNKNOWN_COLOR,
            hatch="///",
            linewidth=0.9,
        )
    )


def draw_breakpoint_line(axes_object, position: int, *, top: float) -> None:
    """画一条从面板底部到给定高度的断点虚线（永远是红色）。

    输入示例: position=100000, top=3.0
    输出示例: 一条竖向红色虚线，贯穿插入轨道下方到块底边
    """
    axes_object.plot(
        [position, position], [-0.5, top], color=BREAKPOINT_COLOR, linestyle="--", linewidth=1.2
    )


### 说明文字 ###


def format_consensus_caption(group: insertion.ClipGroup, prefix: str = "") -> str:
    """一行英文说明：多少条 read、多长、一致率多少。

    输入示例: group.reads 20 条, group.consensus 60 bp, consensus_fraction=1.0
    输出示例: "left: soft-clip consensus sequence (20 reads support, 60bp),
               Consensus fraction: 100.00%"
    """
    lead = f"{prefix} " if prefix else ""
    return (
        f"{lead}soft-clip consensus sequence ({len(group.reads)} reads support, "
        f"{len(group.consensus)}bp), Consensus fraction: {group.consensus_fraction:.2%}"
    )


def draw_consensus_notes(
    figure, entries: list[tuple[str, str]], *, show_consensus: bool, show_formula: bool
) -> None:
    """在图底部写一致性序列说明（可选），右下角写一致率公式（可选）。

    输入示例: entries=[("left: soft-clip consensus…", "GCATTT…")]
    输出示例: 图底部出现该行文字与序列；右下角出现 Consensus fraction 公式
    """
    if show_consensus:
        line_y = 0.085
        for caption, sequence in entries:
            figure.text(0.055, line_y, caption, ha="left", va="top", fontsize=9)
            figure.text(
                0.055,
                line_y - 0.026,
                sequence[:130] + ("…" if len(sequence) > 130 else ""),
                ha="left",
                va="top",
                fontsize=7,
                family="monospace",
            )
            line_y -= 0.062

    if show_formula:
        figure.text(
            0.985,
            0.022,
            "$\\mathrm{Consensus\\ fraction} = \\dfrac{\\sum_{c} \\max_{b} n_{c,b}}"
            "{\\sum_{c} \\sum_{b} n_{c,b}}$   ($n_{c,b}$: reads with base $b$ at column $c$)",
            ha="right",
            va="bottom",
            fontsize=8,
            color=REFERENCE_LINE_COLOR,
        )


### reads 轨道 ###


def draw_reads_track(
    axes_object,
    clip_groups: list[insertion.ClipGroup],
    chrom: str,
    span_start: int,
    span_end: int,
) -> None:
    """画 reads 轨道：左 clip 排在上方、右 clip 排在下方，并在每组第一条 read 旁标注。

    输入示例: clip_groups=[ClipGroup(breakpoint=300000, side="left", reads=15 条)]
    输出示例: 15 行蓝色 read 方块（左端带橙色虚线 clip 框），最上面一行左侧
              出现标注 "breakpoint:300,000(15 soft-clip)"

    说明: 只有两侧都有 clip 时才画中间那条参考基线（它起分区作用）；
          单侧图不画，也不为另一侧预留空间。
    """
    has_left = any(group.side == "left" for group in clip_groups)
    has_right = any(group.side == "right" for group in clip_groups)
    if has_left and has_right:
        axes_object.axhline(0, color=REFERENCE_LINE_COLOR, linewidth=1.2)

    for breakpoint in sorted({group.breakpoint for group in clip_groups}):
        axes_object.axvline(breakpoint, color=BREAKPOINT_COLOR, linestyle="--", linewidth=1.2)

    left_rows = 0
    right_rows = 0
    for group in clip_groups:
        fill_color, edge_color = side_colors(group.side)
        label = f"breakpoint:{coords.left_base_1based(group.breakpoint):,}"
        for index, read in enumerate(group.reads):
            if group.side == "left":
                row_y = 1.4 + left_rows
                left_rows += 1
                clip_left_edge = read.ref_start - read.clip_length
                annotation_x = clip_left_edge
                annotation_offset = -8
                ha = "right"
            else:
                row_y = -1.4 - right_rows
                right_rows += 1
                clip_left_edge = read.ref_end
                annotation_x = clip_left_edge + read.clip_length
                annotation_offset = 8
                ha = "left"

            axes_object.add_patch(
                Rectangle(
                    (read.ref_start, row_y),
                    read.ref_end - read.ref_start,
                    0.7,
                    color=fill_color,
                    alpha=0.9,
                )
            )
            axes_object.add_patch(
                Rectangle(
                    (clip_left_edge, row_y),
                    read.clip_length,
                    0.7,
                    facecolor="none",
                    edgecolor=CLIP_OUTLINE_COLOR,
                    linestyle=":",
                    linewidth=0.9,
                )
            )
            if index == 0:
                axes_object.annotate(
                    f"{label}({len(group.reads)} soft-clip)",
                    xy=(annotation_x, row_y + 0.35),
                    xytext=(annotation_offset, 0),
                    textcoords="offset points",
                    ha=ha,
                    va="center",
                    fontsize=8,
                    color=edge_color,
                    bbox={
                        "boxstyle": "round,pad=0.25",
                        "fc": "white",
                        "ec": fill_color,
                        "linewidth": 0.6,
                        "alpha": 0.95,
                    },
                )

    upper_limit = (1.4 + left_rows + 0.8) if left_rows else 0.8
    lower_limit = (-1.4 - right_rows - 0.8) if right_rows else -0.8
    axes_object.set_ylim(lower_limit, upper_limit)
    axes_object.set_yticks([])
    axes_object.set_xlim(span_start, span_end)
    axes_object.set_xlabel(f"{chrom} coordinate (0-based)", fontsize=9)


### 图 1 / 图 2：单个断点 ###


def draw_single_breakpoint_figure(
    clip_groups: list[insertion.ClipGroup],
    chrom: str,
    span_start: int,
    span_end: int,
    breakpoint: int,
    output_path: str | Path,
    *,
    dpi: int = 150,
    show_consensus: bool = True,
    show_formula: bool = True,
) -> dict[str, object]:
    """画"一个断点"的图：只有一侧 clip 时是单侧图，两侧都有时是同断点双侧图。

    输入示例::
        clip_groups=[ClipGroup(100000,"left",20 条), ClipGroup(100000,"right",20 条)],
        chrom="chr1", span_start=99200, span_end=100800, breakpoint=100000
    输出示例::
        单侧 → [未知区][60bp 块] 贴着断点，reads 轨道只有该侧
        双侧 → [60bp][未知区（合并）][55bp] + 两条连线收到断点
        返回 {"xlim": (…), "junctions": [100000], "groups": [("left", 20), …]}
    """
    left_group = next((group for group in clip_groups if group.side == "left"), None)
    right_group = next((group for group in clip_groups if group.side == "right"), None)
    both_sides = left_group is not None and right_group is not None
    groups = [group for group in (left_group, right_group) if group is not None]
    if not groups:  # 先拦下来：不然会先建一张画布、再抛错
        raise UsageError(f"断点 {coords.left_base_1based(breakpoint):,} 在本区间没有任何 clip read")

    figure, insert_axis, reads_axis = build_figure()

    if both_sides:
        middle_half_width = (
            max(left_group.maximum_clip_length, right_group.maximum_clip_length) * 1.3
        )
        middle_start = breakpoint - middle_half_width
        middle_end = breakpoint + middle_half_width
        left_block_start = middle_start - left_group.maximum_clip_length
        right_block_end = middle_end + right_group.maximum_clip_length

        draw_insert_block(insert_axis, left_block_start, middle_start, left_group)
        draw_unknown_block(insert_axis, middle_start, middle_end)
        draw_insert_block(insert_axis, middle_end, right_block_end, right_group)

        apex_y = insert_axis.get_ylim()[0]
        for block_inner_edge, clip_side in ((middle_start, "left"), (middle_end, "right")):
            insert_axis.plot(
                [block_inner_edge, breakpoint],
                [3.0, apex_y],
                color=side_colors(clip_side)[0],
                linewidth=1.8,
                alpha=0.9,
                solid_capstyle="round",
            )
        insert_axis.scatter(
            [breakpoint], [apex_y], s=28, color=BREAKPOINT_COLOR, zorder=5, clip_on=False
        )
        block_edges = [left_block_start, right_block_end]
    else:
        group = groups[0]
        unknown_width = (
            max(group.maximum_clip_length, UNKNOWN_WIDTH_MIN_CLIP) * UNKNOWN_WIDTH_FACTOR
        )
        if group.side == "left":
            block_start = breakpoint - group.maximum_clip_length
            draw_unknown_block(insert_axis, block_start - unknown_width, block_start)
            draw_insert_block(insert_axis, block_start, breakpoint, group)
            block_edges = [block_start - unknown_width, breakpoint]
        else:
            block_end = breakpoint + group.maximum_clip_length
            draw_insert_block(insert_axis, breakpoint, block_end, group)
            draw_unknown_block(insert_axis, block_end, block_end + unknown_width)
            block_edges = [breakpoint, block_end + unknown_width]

    draw_breakpoint_line(insert_axis, breakpoint, top=3.0)
    x_min, x_max = expand_axis_limits(span_start, span_end, [int(edge) for edge in block_edges])
    insert_axis.set_xlim(x_min, x_max)

    draw_reads_track(reads_axis, groups, chrom, span_start, span_end)
    reads_axis.set_xlim(x_min, x_max)  # 与插入面板一致，保证上下对齐

    draw_consensus_notes(
        figure,
        [(format_consensus_caption(group), group.consensus) for group in groups],
        show_consensus=show_consensus,
        show_formula=show_formula,
    )
    save_figure(figure, output_path, dpi)
    return {
        "xlim": (x_min, x_max),
        "junctions": [breakpoint],
        "groups": [(group.side, len(group.reads)) for group in groups],
    }


### 图 3：多个断点 ###


def plan_multi_breakpoint_blocks(
    groups: list[insertion.ClipGroup], span_start: int, span_end: int
) -> tuple[list[tuple[int, int, str, insertion.ClipGroup | None]], tuple[float, float]]:
    """为多断点图规划每个块的横向范围，并裁掉会撞上下一块的未知区。

    输入示例::
        groups=[组(700000,left), 组(700230,right)]，窗口 699100–700900
    输出示例::

        ([(699680, 699940, "unknown", None), (699940, 700000, "clip", 组A),
          (700230, 700285, "clip", 组B), (700285, 700545, "unknown", None)],
         (699608, 700972))

    说明: 不同断点属于不同插入，未知区不得越过下一块的起点，并至少留 ``gap`` bp，
          避免两块贴在一起被误读成同一段插入序列。
    """
    widest_clip = max((group.maximum_clip_length for group in groups), default=100)
    unknown_width = int(max(widest_clip, UNKNOWN_WIDTH_MIN_CLIP) * UNKNOWN_WIDTH_FACTOR)

    blocks: list[tuple[int, int, str, insertion.ClipGroup | None]] = []
    for group in groups:
        if group.side == "left":
            block_start = group.breakpoint - group.maximum_clip_length
            blocks.append((block_start - unknown_width, block_start, "unknown", None))
            blocks.append((block_start, group.breakpoint, "clip", group))
        else:
            block_end = group.breakpoint + group.maximum_clip_length
            blocks.append((group.breakpoint, block_end, "clip", group))
            blocks.append((block_end, block_end + unknown_width, "unknown", None))

    raw_min = min([span_start, *(block[0] for block in blocks)])
    raw_max = max([span_end, *(block[1] for block in blocks)])
    gap = max(BLOCK_GAP_MIN_BP, int((raw_max - raw_min) * BLOCK_GAP_FRACTION))

    blocks.sort(key=lambda block: block[0])
    trimmed: list[tuple[int, int, str, insertion.ClipGroup | None]] = []
    for index, (block_start, block_end, kind, group) in enumerate(blocks):
        if kind == "unknown":
            for next_start, _next_end, _next_kind, _next_group in blocks[index + 1 :]:
                if next_start < block_end:  # 任何后继块都算阻挡（含相对而伸的未知区）
                    block_end = max(block_start, next_start - gap)
                    break
        if block_end - block_start > 20:  # 太窄的碎条不画，避免出现"针尖"状未知区
            trimmed.append((block_start, block_end, kind, group))

    x_min, x_max = expand_axis_limits(
        span_start, span_end, [block[0] for block in trimmed] + [block[1] for block in trimmed]
    )
    return trimmed, (x_min, x_max)


def draw_multiple_breakpoints_figure(
    clip_groups: list[insertion.ClipGroup],
    chrom: str,
    span_start: int,
    span_end: int,
    breakpoints: list[int],
    output_path: str | Path,
    *,
    dpi: int = 150,
    show_consensus: bool = True,
    show_formula: bool = True,
) -> dict[str, object]:
    """画 2–3 个不同断点的图：每个断点自己一组 [未知区][块] / [块][未知区]。

    输入示例: clip_groups=[组(700000,left), 组(700230,right)]，窗口 699100–700900
    输出示例: 两条互不重叠、块间留白的分组；reads 轨道上左断点在上、右断点在下
              返回 {"xlim": (…), "junctions": [700000, 700230], "groups": […]}
    """
    blocks, (x_min, x_max) = plan_multi_breakpoint_blocks(clip_groups, span_start, span_end)

    figure, insert_axis, reads_axis = build_figure()
    insert_axis.set_xlim(x_min, x_max)

    for block_start, block_end, kind, group in blocks:
        if kind == "unknown" or group is None:
            draw_unknown_block(insert_axis, block_start, block_end, top=4.2, height=3.4)
        else:
            draw_insert_block(insert_axis, block_start, block_end, group, top=4.2, height=3.4)

    for breakpoint in sorted({group.breakpoint for group in clip_groups}):
        draw_breakpoint_line(insert_axis, breakpoint, top=4.2)

    draw_reads_track(reads_axis, clip_groups, chrom, span_start, span_end)
    reads_axis.set_xlim(x_min, x_max)

    show_breakpoint_prefix = len(clip_groups) > 1
    entries = []
    for group in clip_groups:
        prefix = (
            f"{coords.left_base_1based(group.breakpoint):,} {group.side}:"
            if show_breakpoint_prefix
            else ""
        )
        entries.append((format_consensus_caption(group, prefix), group.consensus))
    draw_consensus_notes(figure, entries, show_consensus=show_consensus, show_formula=show_formula)

    save_figure(figure, output_path, dpi)
    return {
        "xlim": (x_min, x_max),
        "junctions": sorted({group.breakpoint for group in clip_groups}),
        "groups": [(group.side, len(group.reads)) for group in clip_groups],
    }


def render(
    clip_groups: list[insertion.ClipGroup],
    chrom: str,
    span_start: int,
    span_end: int,
    breakpoints: list[int],
    output_path: str | Path,
    *,
    dpi: int = 150,
    show_consensus: bool = True,
    show_formula: bool = True,
) -> dict[str, object]:
    """按断点个数选画法并出图（1 个 → 单侧/双侧图；2 个以上 → 多断点图）。

    输入示例: breakpoints=[300000], clip_groups=[组(300000,"left",15 条)]
    输出示例: out/demo.png + {"xlim": (299734, 300060), "junctions": [300000], …}
    """
    arguments = {
        "dpi": dpi,
        "show_consensus": show_consensus,
        "show_formula": show_formula,
    }
    if len(breakpoints) == 1:
        return draw_single_breakpoint_figure(
            clip_groups, chrom, span_start, span_end, breakpoints[0], output_path, **arguments
        )
    return draw_multiple_breakpoints_figure(
        clip_groups, chrom, span_start, span_end, breakpoints, output_path, **arguments
    )


### 编排：BAM → 两张表 → PNG ###


@dataclass(frozen=True, slots=True)
class InsertionSettings:
    """一次 ``short insertion-plot`` 的全部输入。

    示例::

        InsertionSettings(bam="reads.bam", chrom="chr1", window_start=99200,
                          window_end=100800, breakpoints=(100000,), out_dir="out")
    """

    bam: str
    chrom: str
    window_start: int  # 0-based 半开
    window_end: int
    breakpoints: tuple[int, ...]  # 0-based 断点边界
    out_dir: str = "."
    name: str = "insertion_plot"
    dpi: int = 150
    show_consensus: bool = True
    show_formula: bool = True


def run(settings: InsertionSettings) -> dict[str, object]:
    """读 BAM → 写 ``insertion_clip_reads`` → 读回该表 → 写一致性表 → 由表出图。

    输入示例: InsertionSettings(bam="reads.bam", chrom="chr1", window_start=99200,
                                window_end=100800, breakpoints=(100000,))
    输出示例: {"png": "out/insertion_plot.png", "clips_tsv": "…", "consensus_tsv": "…",
               "junctions": [100000], "n_reads": 40, "xlim": (99734, 100660)}

    为什么要“写表再读表”：表是接口。图只允许看到表里的东西，读回来一次就能保证
    表本身足够画出这张图（少列、错列会在这里立刻报错，而不是画出一张空图）。
    """
    if settings.window_start < 0 or settings.window_start >= settings.window_end:
        raise UsageError("--region 的起点必须 >= 1 且小于终点（1-based 闭区间）")
    if settings.dpi <= 0:
        raise UsageError("--dpi 必须为正数")
    if "/" in settings.name or "\\" in settings.name:
        raise UsageError("--name 只能是文件名前缀，不能带目录（用 --out-dir）")

    clip_reads = insertion.load_clip_reads(
        settings.bam, settings.chrom, settings.window_start, settings.window_end
    )
    if not clip_reads:
        raise RunError(
            f"{settings.chrom}:{settings.window_start + 1}-{settings.window_end} 里没有带 "
            "soft clip 的主比对（检查 --region、BAM 的 -Y 比对参数）"
        )
    breakpoints = insertion.validate_breakpoints(settings.breakpoints, clip_reads)

    output_directory = Path(settings.out_dir)
    sample = Path(settings.bam).name.split(".")[0] or "sample"
    clips_path = output_directory / f"{settings.name}.clips.tsv"
    consensus_path = output_directory / f"{settings.name}.consensus.tsv"
    png_path = output_directory / f"{settings.name}.png"
    try:
        output_directory.mkdir(parents=True, exist_ok=True)
    except OSError as error:  # cli.main 只捕 AlnkitError，不能把 OSError 漏出去
        raise RunError(f"无法创建输出目录 {output_directory}：{error}") from error

    comments = (
        f"sample={sample}",
        "window="
        + coords.format_region_1based(settings.chrom, settings.window_start, settings.window_end),
        "producer=alnkit short insertion-plot",
    )
    insertion.write_clip_table(clips_path, clip_reads, settings.chrom, sample, comments)

    # 从这里开始只认表：图与一致性表都由“表里读回来的证据”算出，输入换了就一起换
    table_reads, chrom = insertion.read_clip_table(clips_path)
    groups = insertion.build_groups(table_reads, breakpoints)
    insertion.write_consensus_table(consensus_path, groups, chrom, sample, comments)

    layout = render(
        groups,
        chrom,
        settings.window_start,
        settings.window_end,
        list(breakpoints),
        png_path,
        dpi=settings.dpi,
        show_consensus=settings.show_consensus,
        show_formula=settings.show_formula,
    )
    return {
        "png": str(png_path),
        "clips_tsv": str(clips_path),
        "consensus_tsv": str(consensus_path),
        "junctions": layout["junctions"],
        "n_reads": len(table_reads),
        "xlim": layout["xlim"],
        "groups": layout["groups"],
    }
