"""``alnkit selftest``：不依赖任何输入数据的断言集。

M0 覆盖坐标换算与表/段契约；二代 B1 的 30+ 条口径断言在 P2 并入这里
。任何一条失败 → 退出码 1。
"""

from __future__ import annotations

import argparse
from collections.abc import Callable
from typing import Any

from ..core import coords, segments, tables
from ..core.errors import UsageError
from ._common import GLOBAL_OPTIONS

_CHECKS: list[tuple[str, Callable[[], None]]] = []


def check(name: str) -> Callable[[Callable[[], None]], Callable[[], None]]:
    """把一个无参函数注册成一条自检（失败即抛异常）。"""

    def register(func: Callable[[], None]) -> Callable[[], None]:
        _CHECKS.append((name, func))
        return func

    return register


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "selftest",
        parents=[GLOBAL_OPTIONS],
        help="运行内置自检（不需要输入数据）",
        description="校验坐标换算与表契约；任何一条失败则退出码为 1。",
        epilog="例子：\n  alnkit selftest",
    )
    parser.set_defaults(handler=run)


# ------------------------------------------------------------------ coords
@check("coords: 1-based <-> 0-based 往返")
def _coords_position_roundtrip() -> None:
    for position_0based in (0, 1, 9424, 1_000_000, 10**9):
        back = coords.to_0based(coords.to_1based(position_0based))
        assert back == position_0based, (position_0based, back)


@check("coords: 区间串往返（1-based 闭 <-> 0-based 半开）")
def _coords_region_roundtrip() -> None:
    text = coords.format_region_1based("chr1", 2_000_000, 2_000_050)
    assert text == "chr1:2000001-2000050", text
    assert coords.parse_region_1based(text) == ("chr1", 2_000_000, 2_000_050)
    assert coords.parse_region_1based("chr1:1,000-2,000") == ("chr1", 999, 2000)


@check("coords: 区间换算互逆")
def _coords_interval_roundtrip() -> None:
    for start_0based, end_0based in ((0, 1), (99, 200), (2_000_000, 2_000_050)):
        start_1based, end_1based = coords.interval_0based_to_1based(start_0based, end_0based)
        assert (start_1based, end_1based) == (start_0based + 1, end_0based)
        assert coords.interval_1based_to_0based(start_1based, end_1based) == (
            start_0based,
            end_0based,
        )


@check("coords: 断点边界 / 左碱基恒等式（Insert_Position 口径）")
def _coords_boundary_identity() -> None:
    boundary = 1_000_001
    assert coords.left_base_1based(boundary) == 1_000_001
    assert coords.left_base_0based(boundary) == 1_000_000
    assert coords.left_base_0based(boundary) + 1 == coords.left_base_1based(boundary)


@check("coords: 反向链 read 区间翻转")
def _coords_flip_read_interval() -> None:
    assert coords.flip_read_interval(40_592, 9_424, 40_591) == (1, 31_168)
    assert coords.flip_read_interval(150, 0, 150) == (0, 150)
    try:
        coords.flip_read_interval(100, 10, 200)
    except ValueError:
        pass
    else:  # pragma: no cover - 防御式断言
        raise AssertionError("越界区间应当报错")


# ------------------------------------------------------------------ tables
@check("tables: 所有注册表满足坐标后缀规则")
def _tables_registry_valid() -> None:
    for table in tables.COLUMNS:
        tables.validate_columns(table)


@check("tables: 坐标列缺少基准后缀会被拒绝")
def _tables_reject_bare_coordinate() -> None:
    try:
        tables.validate_column_names("demo", ("ref_start", "read_end"))
    except AssertionError:
        return
    raise AssertionError("裸坐标列名 ref_start/read_end 应当被拒绝")


@check("tables: 写入的列序与注册表一致")
def _tables_column_order() -> None:
    columns = tables.columns_of("insert_sites")
    assert columns[0] == "chrom"
    assert columns[1] == "insert_pos_1based"
    assert columns[3] == "insert_boundary_0based"


# ---------------------------------------------------------------- segments
_SAME_NUMBERS = {
    "read_name": "r1",
    "chrom": "chr1",
    "ref_start": 9_000_000,
    "ref_end": 9_031_165,
    "read_start": 9_424,
    "read_end": 40_591,
    "strand": "-",
    "mapq": 60,
}


@check("segments: 三套历史表头解析出同一组坐标")
def _segments_three_legacy_headers() -> None:
    alignment_table_header = [
        "chrom",
        "chrom_length",
        "ref_pos_start",
        "ref_pos_end",
        "read_name",
        "read_pos_start",
        "read_pos_end",
        "strand",
        "mapping_quality",
        "is_secondary",
        "is_supplementary",
        "ref_aln_len",
        "read_aln_len",
        "simplified_cigar",
    ]
    alignment_table_row = {
        "chrom": "chr1",
        "chrom_length": "100000000",
        "ref_pos_start": "9000000",
        "ref_pos_end": "9031165",
        "read_name": "r1",
        "read_pos_start": "9424",
        "read_pos_end": "40591",
        "strand": "-",
        "mapping_quality": "60",
        "is_secondary": "0",
        "is_supplementary": "0",
        "ref_aln_len": "31165",
        "read_aln_len": "31167",
        "simplified_cigar": "1S31204M9424S",
    }
    classify_header = [
        "read_name",
        "segment_role",
        "chrom",
        "chrom_length",
        "ref_start",
        "ref_end",
        "read_start",
        "read_end",
        "strand",
        "mapq",
        "ref_aln_len",
        "read_aln_len",
        "simplified_cigar",
        "kept",
        "drop_reason",
    ]
    classify_row = {
        "read_name": "r1",
        "segment_role": "genome",
        "chrom": "chr1",
        "chrom_length": "100000000",
        "ref_start": "9000000",
        "ref_end": "9031165",
        "read_start": "9424",
        "read_end": "40591",
        "strand": "-",
        "mapq": "60",
        "ref_aln_len": "31165",
        "read_aln_len": "31167",
        "simplified_cigar": "1S31204M9424S",
        "kept": "1",
        "drop_reason": "-",
    }
    tracks_header = [
        "read_name",
        "read_length",
        "segment_type",
        "target_name",
        "target_start",
        "target_end",
        "read_start",
        "read_end",
        "strand",
        "mapq",
    ]
    tracks_row = {
        "read_name": "r1",
        "read_length": "40592",
        "segment_type": "genome",
        "target_name": "chr1",
        "target_start": "9000000",
        "target_end": "9031165",
        "read_start": "9424",
        "read_end": "40591",
        "strand": "-",
        "mapq": "60",
    }

    parsed = [
        segments.segments_from_rows(header, [row])[0]
        for header, row in (
            (alignment_table_header, alignment_table_row),
            (classify_header, classify_row),
            (tracks_header, tracks_row),
        )
    ]
    expected = (
        _SAME_NUMBERS["ref_start"],
        _SAME_NUMBERS["ref_end"],
        _SAME_NUMBERS["read_start"],
        _SAME_NUMBERS["read_end"],
    )
    for segment in parsed:
        got = (
            segment.ref_start_0based,
            segment.ref_end_0based,
            segment.read_start_0based,
            segment.read_end_0based,
        )
        assert got == expected, (got, expected)
        assert segment.chrom == "chr1"
        assert segment.strand == "-"
        assert segment.mapq == 60


@check("segments: 缺列 / 倒置区间 / 坏 strand 都报 UsageError")
def _segments_error_paths() -> None:
    good_header = [
        "read_name",
        "chrom",
        "ref_start",
        "ref_end",
        "read_start",
        "read_end",
        "strand",
        "mapq",
    ]
    good_row = {
        "read_name": "r1",
        "chrom": "chr1",
        "ref_start": "1",
        "ref_end": "2",
        "read_start": "0",
        "read_end": "1",
        "strand": "+",
        "mapq": "60",
    }

    cases: list[tuple[list[str], dict[str, Any]]] = [
        (good_header, {**good_row, "ref_start": "10", "ref_end": "5"}),
        (good_header, {**good_row, "read_start": "9", "read_end": "1"}),
        (good_header, {**good_row, "strand": "?"}),
        (["read_name", "chrom"], {"read_name": "r1", "chrom": "chr1"}),
    ]
    for header, row in cases:
        try:
            segments.segments_from_rows(header, [row])
        except UsageError:
            continue
        raise AssertionError(f"应当报 UsageError：{header} {row}")


@check("short scan: 断点坐标 / 归属 / 窗口合并口径（原 parse_bam --self-test）")
def _short_scan_self_test() -> None:
    import contextlib
    import io

    from ..shortread import scan

    buffer = io.StringIO()
    with contextlib.redirect_stdout(buffer):
        scan.run_self_test()
    assert "self-test 通过" in buffer.getvalue()


def run(args: argparse.Namespace) -> int:
    del args
    failures: list[str] = []
    total = len(_CHECKS)
    for name, func in _CHECKS:
        try:
            func()
        except Exception as exc:  # noqa: BLE001 - 自检要汇总所有失败
            failures.append(name)
            print(f"FAIL  {name}\n        {type(exc).__name__}: {exc}")
        else:
            print(f"PASS  {name}")
    print(f"\n{total - len(failures)}/{total} 通过")
    return 1 if failures else 0
