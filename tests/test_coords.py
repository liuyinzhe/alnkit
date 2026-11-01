"""坐标换算与区间串（core/coords.py）。"""

from __future__ import annotations

import pytest

from alnkit.core import coords


@pytest.mark.parametrize("position_0based", [0, 1, 9424, 1_000_000])
def test_position_roundtrip(position_0based: int) -> None:
    assert coords.to_0based(coords.to_1based(position_0based)) == position_0based


def test_position_rejects_negative_and_zero() -> None:
    with pytest.raises(ValueError):
        coords.to_1based(-1)
    with pytest.raises(ValueError):
        coords.to_0based(0)


def test_interval_conversions_are_exact() -> None:
    # 0-based 半开 [2_000_000, 2_000_050) 是 50 bp，对应 1-based 闭 2000001-2000050
    assert coords.interval_0based_to_1based(2_000_000, 2_000_050) == (2_000_001, 2_000_050)
    assert coords.interval_1based_to_0based(2_000_001, 2_000_050) == (2_000_000, 2_000_050)


def test_interval_rejects_inverted_or_negative() -> None:
    with pytest.raises(ValueError):
        coords.interval_0based_to_1based(10, 5)
    with pytest.raises(ValueError):
        coords.interval_0based_to_1based(-1, 5)
    with pytest.raises(ValueError):
        coords.interval_1based_to_0based(0, 5)


def test_region_parse_and_format_roundtrip() -> None:
    text = coords.format_region_1based("chr1", 2_000_000, 2_000_050)
    assert text == "chr1:2000001-2000050"
    assert coords.parse_region_1based(text) == ("chr1", 2_000_000, 2_000_050)


def test_region_parse_accepts_thousands_separator() -> None:
    assert coords.parse_region_1based("chr1:1,000-2,000") == ("chr1", 999, 2000)


@pytest.mark.parametrize("text", ["chr1", "chr1:100", "chr1:0-10", "chr1:20-10", "chr1:a-b"])
def test_region_parse_rejects_bad_input(text: str) -> None:
    with pytest.raises(ValueError):
        coords.parse_region_1based(text)


def test_left_base_identity_matches_insert_position() -> None:
    """legacy Insert_Position：数值 = 0-based 边界 = 左侧碱基的 1-based 坐标。"""
    boundary = 1_000_001
    assert coords.left_base_1based(boundary) == boundary
    assert coords.left_base_0based(boundary) == boundary - 1
    assert coords.left_base_0based(boundary) + 1 == coords.left_base_1based(boundary)


def test_flip_read_interval() -> None:
    # 1S31204M9424S 的反向链比对段：参考两端 -> 正向 read 轴
    assert coords.flip_read_interval(40_592, 9_424, 40_591) == (1, 31_168)
    assert coords.flip_read_interval(150, 0, 150) == (0, 150)


def test_flip_read_interval_rejects_out_of_range() -> None:
    with pytest.raises(ValueError):
        coords.flip_read_interval(100, 10, 200)


def test_merge_intervals() -> None:
    assert coords.merge_intervals([(0, 10), (5, 20), (30, 40)]) == [(0, 20), (30, 40)]
    assert coords.merge_intervals([(10, 20), (0, 10)]) == [(0, 20)]  # 相邻也合并
    assert coords.merge_intervals([]) == []


def test_zero_length_interval_cannot_be_formatted() -> None:
    # 零长区间无 1-based 闭表示（start+1 > end），必须报错而不是产出解析不回去的串
    with pytest.raises(ValueError):
        coords.format_region_1based("chr1", 5, 5)


def test_left_base_requires_a_left_base() -> None:
    with pytest.raises(ValueError):
        coords.left_base_0based(0)
    assert coords.left_base_0based(1) == 0


@pytest.mark.parametrize("text", ["chr1:1,5-20", "chr1:1234,567-2000", "chr1:1,000,0-5"])
def test_region_parse_rejects_misplaced_commas(text: str) -> None:
    with pytest.raises(ValueError):
        coords.parse_region_1based(text)


def test_region_parse_accepts_proper_thousands() -> None:
    assert coords.parse_region_1based("chr1:1,000,000-1,000,500") == ("chr1", 999_999, 1_000_500)
