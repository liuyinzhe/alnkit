"""比对段模型与表读写：三套历史表头的统一入口。

历史上有三套“比对段表”（14 列 / 15 列 / 10 列，列名各不相同），语义一致、
基准一致。本模块只做一件事：
**把任意一套读成同一个 ``Segment``**，写出时统一用规范列名。
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass, fields
from typing import Any

from .errors import UsageError
from .tables import read_tsv

__all__ = ["ALIASES", "REQUIRED_FIELDS", "Segment", "read_segments", "segment_rows"]

#: 规范字段 -> 可接受的输入列名（第一项是规范名）
ALIASES: dict[str, tuple[str, ...]] = {
    "read_name": ("read_name",),
    "read_length": ("read_length",),
    "segment_role": ("segment_role", "segment_type"),
    "chrom": ("chrom", "target_name", "Chromosome"),
    "chrom_length": ("chrom_length",),
    "ref_start_0based": (
        "ref_start_0based",
        "ref_start",
        "ref_pos_start",
        "M_ref_pos_start",
        "target_start",
        "segment_start",
    ),
    "ref_end_0based": (
        "ref_end_0based",
        "ref_end",
        "ref_pos_end",
        "M_ref_pos_end",
        "target_end",
        "segment_end",
    ),
    "read_start_0based": (
        "read_start_0based",
        "read_start",
        "read_pos_start",
        "M_read_pos_start",
    ),
    "read_end_0based": ("read_end_0based", "read_end", "read_pos_end", "M_read_pos_end"),
    "strand": ("strand", "stand"),  # stand：历史拼写错误
    "mapq": ("mapq", "mapping_quality"),
    "ref_aln_len": ("ref_aln_len",),
    "read_aln_len": ("read_aln_len",),
    "cigar": ("cigar", "simplified_cigar"),
    "kept": ("kept",),
    "drop_reason": ("drop_reason",),
    "is_secondary": ("is_secondary",),
    "is_supplementary": ("is_supplementary",),
    "nm": ("nm",),
    "match_columns": ("match_columns",),
    "inserted_bases": ("inserted_bases",),
    "deleted_bases": ("deleted_bases",),
}

REQUIRED_FIELDS: tuple[str, ...] = (
    "read_name",
    "chrom",
    "ref_start_0based",
    "ref_end_0based",
    "read_start_0based",
    "read_end_0based",
    "strand",
    "mapq",
)

#: 空值标记（表里统一用 "-" 表示空）
_EMPTY = {"", "-", "NA", "N/A", "None", "null"}
#: strand 列例外："-" 是合法值，不能当空
_EMPTY_TEXT = {"", "NA", "N/A", "None", "null"}


@dataclass(slots=True, frozen=True)
class Segment:
    """一条 read 的一段比对（参考 + read 双坐标，0-based 半开）。"""

    read_name: str
    chrom: str
    ref_start_0based: int
    ref_end_0based: int
    read_start_0based: int
    read_end_0based: int
    strand: str
    mapq: int
    read_length: int | None = None
    segment_role: str | None = None
    chrom_length: int | None = None
    ref_aln_len: int | None = None
    read_aln_len: int | None = None
    cigar: str | None = None
    kept: bool = True
    drop_reason: str | None = None
    is_secondary: bool = False
    is_supplementary: bool = False
    nm: int | None = None
    match_columns: int | None = None
    inserted_bases: int | None = None
    deleted_bases: int | None = None

    @property
    def ref_span(self) -> int:
        """参考跨度（由坐标算，不看 ``ref_aln_len``）。"""
        return self.ref_end_0based - self.ref_start_0based

    @property
    def read_span(self) -> int:
        """read 跨度（由坐标算）。"""
        return self.read_end_0based - self.read_start_0based


def _clean(value: str | None) -> str | None:
    if value is None:
        return None
    text = value.strip()
    return None if text in _EMPTY else text


def _clean_text(value: str | None) -> str | None:
    """只把空白/NA 当空，保留 "-"（strand 用）。"""
    if value is None:
        return None
    text = value.strip()
    return None if text in _EMPTY_TEXT else text


def _to_int(value: str | None, field: str, where: str) -> int:
    text = _clean(value)
    if text is None:
        raise UsageError(f"{where}: 列 {field} 为空，但它是必需的整数列")
    try:
        return int(text)
    except ValueError:
        raise UsageError(f"{where}: 列 {field} 不是整数：{text!r}") from None


def _to_opt_int(value: str | None, field: str, where: str) -> int | None:
    text = _clean(value)
    if text is None:
        return None
    try:
        return int(text)
    except ValueError:
        raise UsageError(f"{where}: 列 {field} 不是整数：{text!r}") from None


def _strand_or_fail(value: str | None, where: str) -> str:
    """strand 是必需语义列：空/NA 直接报错，绝不默认成 "+"（会静默算错断点侧）。"""
    text = _clean_text(value)
    if text is None:
        raise UsageError(f"{where}: strand 列为空；必须是 + 或 -")
    return text


def _to_bool(value: str | None, default: bool = True) -> bool:
    text = _clean(value)
    if text is None:
        return default
    return text.lower() in {"1", "true", "yes", "y"}


def _resolve_columns(header: Sequence[str], where: str) -> dict[str, str]:
    """规范字段 -> 该表实际列名（找不到就不放进结果）。"""
    present = set(header)
    resolved: dict[str, str] = {}
    for canonical, candidates in ALIASES.items():
        for candidate in candidates:
            if candidate in present:
                resolved[canonical] = candidate
                break
    missing = [name for name in REQUIRED_FIELDS if name not in resolved]
    if missing:
        raise UsageError(
            f"{where}: 缺少必需列 {missing}；表头是 {list(header)}"
            "（可用 --columns legacy:* 指定历史表头）"
        )
    return resolved


def segments_from_rows(
    header: Sequence[str], rows: Iterable[Mapping[str, str]], where: str = "<rows>"
) -> list[Segment]:
    """把 (表头, 行) 转成 ``Segment`` 列表。"""
    mapping = _resolve_columns(header, where)

    def value(row: Mapping[str, str], canonical: str) -> str | None:
        source = mapping.get(canonical)
        return None if source is None else row.get(source)

    segments: list[Segment] = []
    for index, row in enumerate(rows, start=2):  # 2 = 首个数据行
        here = f"{where}:{index}"
        segment = Segment(
            read_name=_clean(value(row, "read_name")) or "",
            chrom=_clean(value(row, "chrom")) or "",
            ref_start_0based=_to_int(value(row, "ref_start_0based"), "ref_start_0based", here),
            ref_end_0based=_to_int(value(row, "ref_end_0based"), "ref_end_0based", here),
            read_start_0based=_to_int(value(row, "read_start_0based"), "read_start_0based", here),
            read_end_0based=_to_int(value(row, "read_end_0based"), "read_end_0based", here),
            strand=_strand_or_fail(value(row, "strand"), here),
            mapq=_to_int(value(row, "mapq"), "mapq", here),
            read_length=_to_opt_int(value(row, "read_length"), "read_length", here),
            segment_role=_clean(value(row, "segment_role")),
            chrom_length=_to_opt_int(value(row, "chrom_length"), "chrom_length", here),
            ref_aln_len=_to_opt_int(value(row, "ref_aln_len"), "ref_aln_len", here),
            read_aln_len=_to_opt_int(value(row, "read_aln_len"), "read_aln_len", here),
            cigar=_clean(value(row, "cigar")),
            kept=_to_bool(value(row, "kept")),
            drop_reason=_clean(value(row, "drop_reason")),
            is_secondary=_to_bool(value(row, "is_secondary"), default=False),
            is_supplementary=_to_bool(value(row, "is_supplementary"), default=False),
            nm=_to_opt_int(value(row, "nm"), "nm", here),
            match_columns=_to_opt_int(value(row, "match_columns"), "match_columns", here),
            inserted_bases=_to_opt_int(value(row, "inserted_bases"), "inserted_bases", here),
            deleted_bases=_to_opt_int(value(row, "deleted_bases"), "deleted_bases", here),
        )
        if segment.ref_end_0based < segment.ref_start_0based:
            raise UsageError(
                f"{here}: 参考区间倒置 [{segment.ref_start_0based}, {segment.ref_end_0based})"
            )
        if segment.read_end_0based < segment.read_start_0based:
            raise UsageError(
                f"{here}: read 区间倒置 [{segment.read_start_0based}, {segment.read_end_0based})"
            )
        if segment.strand not in {"+", "-"}:
            raise UsageError(f"{here}: strand 只能是 + 或 -，得到 {segment.strand!r}")
        segments.append(segment)
    return segments


def read_segments(path: str | os.PathLike[str]) -> list[Segment]:
    """读一比对段表（自动识别三套历史表头）。"""
    header, rows = read_tsv(path)
    if not header:
        raise UsageError(f"{os.fspath(path)}: 空表或没有表头")
    return segments_from_rows(header, rows, os.fspath(path))


def segment_rows(segments: Iterable[Segment]) -> list[dict[str, Any]]:
    """``Segment`` 列表 -> 可直接交给 ``tables.write_tsv(..., "segments", ...)`` 的行。"""
    names = {field.name for field in fields(Segment)}
    return [
        {key: value for key, value in asdict(segment).items() if key in names}
        for segment in segments
    ]
