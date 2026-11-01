"""坐标换算与区间字符串：全包唯一的 1-based / 0-based 转换入口。

约定：

* 内部计算一律 **0-based 半开** ``[start, end)``；
* 区间展示串一律 **1-based 闭区间** ``chrom:start-end``（IGV/samtools/VCF 习惯）；
* 断点这类“单点”用 ``*_pos_1based`` 表示**边界左侧碱基**的 1-based 坐标，
  其 0-based 边界值与之数值相同（历史巧合，见 ``left_base_*``）。

任何模块都不得手写 ``+1`` / ``-1``，只能调用本模块。
"""

from __future__ import annotations

import re
from collections.abc import Iterable

__all__ = [
    "REGION_PATTERN",
    "merge_intervals",
    "flip_read_interval",
    "format_region_1based",
    "interval_0based_to_1based",
    "interval_1based_to_0based",
    "left_base_0based",
    "left_base_1based",
    "parse_region_1based",
    "to_0based",
    "to_1based",
]

# 数字允许千分位逗号，但只在正确位置：1,000 合法，1,5 不合法
REGION_PATTERN = re.compile(
    r"^(?P<chrom>[^:\s]+):(?P<start>\d{1,3}(?:,\d{3})*|\d+)-(?P<end>\d{1,3}(?:,\d{3})*|\d+)$"
)


def to_1based(position_0based: int) -> int:
    """碱基坐标：0-based -> 1-based。"""
    if position_0based < 0:
        raise ValueError(f"0-based 坐标不能为负：{position_0based}")
    return position_0based + 1


def to_0based(position_1based: int) -> int:
    """碱基坐标：1-based -> 0-based。"""
    if position_1based < 1:
        raise ValueError(f"1-based 坐标必须 >= 1：{position_1based}")
    return position_1based - 1


def interval_0based_to_1based(start_0based: int, end_0based: int) -> tuple[int, int]:
    """半开区间 ``[s, e)`` -> 1-based 闭区间 ``(s+1, e)``。"""
    if start_0based < 0:
        raise ValueError(f"0-based 起点不能为负：{start_0based}")
    if end_0based < start_0based:
        raise ValueError(f"区间终点小于起点：[{start_0based}, {end_0based})")
    if end_0based == start_0based:
        # 零长区间无法表示为 1-based 闭区间（start+1 > end），直接拒绝而不是产出解析不回去的串
        raise ValueError(f"零长区间无法表示：{start_0based}")
    return start_0based + 1, end_0based


def interval_1based_to_0based(start_1based: int, end_1based: int) -> tuple[int, int]:
    """1-based 闭区间 ``(s, e)`` -> 半开区间 ``[s-1, e)``。"""
    if start_1based < 1:
        raise ValueError(f"1-based 起点必须 >= 1：{start_1based}")
    if end_1based < start_1based:
        raise ValueError(f"区间终点小于起点：{start_1based}-{end_1based}")
    return start_1based - 1, end_1based


def left_base_1based(boundary_0based: int) -> int:
    """边界左侧碱基的 1-based 坐标（= legacy ``Insert_Position`` 的口径）。

    0-based 边界 B 表示“其左侧有 B 个碱基”，故左侧碱基的 1-based 编号恰为 B。
    """
    if boundary_0based < 0:
        raise ValueError(f"边界不能为负：{boundary_0based}")
    return boundary_0based


def left_base_0based(boundary_0based: int) -> int:
    """边界左侧碱基的 0-based 坐标（= 1-based 值 - 1）。

    边界为 0 表示插入发生在第一个碱基之前，没有“左侧碱基”，抛 ``ValueError``。
    """
    if boundary_0based < 1:
        raise ValueError(f"边界 {boundary_0based} 没有左侧碱基")
    return left_base_1based(boundary_0based) - 1


def flip_read_interval(read_length: int, start: int, end: int) -> tuple[int, int]:
    """把反向链比对段的 read 区间翻到正向 read 轴（0-based 半开）。"""
    if start < 0 or end > read_length or end < start:
        raise ValueError(f"read 区间越界：[{start}, {end}) / 长度 {read_length}")
    return read_length - end, read_length - start


def parse_region_1based(text: str) -> tuple[str, int, int]:
    """``"chr1:1,000-2,000"`` -> ``("chr1", 999, 2000)``（0-based 半开）。"""
    match = REGION_PATTERN.match(text.strip())
    if not match:
        raise ValueError(f"区间格式应为 chrom:start-end（1-based）：{text!r}")
    chrom = match.group("chrom")
    start_1based = int(match.group("start").replace(",", ""))
    end_1based = int(match.group("end").replace(",", ""))
    start_0based, end_0based = interval_1based_to_0based(start_1based, end_1based)
    return chrom, start_0based, end_0based


def format_region_1based(chrom: str, start_0based: int, end_0based: int) -> str:
    """``("chr1", 999, 2000)`` -> ``"chr1:1000-2000"``（0-based 半开输入）。"""
    start_1based, end_1based = interval_0based_to_1based(start_0based, end_0based)
    return f"{chrom}:{start_1based}-{end_1based}"


def merge_intervals(intervals: Iterable[tuple[int, int]]) -> list[tuple[int, int]]:
    """合并重叠或相邻的半开区间（0-based）。

    例：``[(0, 10), (5, 20), (30, 40)]`` -> ``[(0, 20), (30, 40)]``。
    """
    merged: list[list[int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]
