"""``plot insertion`` 的数据层：BAM → 逐 read 的 clip 证据表 + 一致性序列表。

流水（``→`` 是数据流，不是调用顺序的猜测）::

    identify_soft_clip → load_clip_reads ─┐
                                          ├→ build_clip_group → 一致性行
    build_clip_consensus ─────────────────┘
                                          └→ 写表 → plot/insertion.render（只读表）

本模块一行 matplotlib 都没有：图由 :mod:`alnkit.plot.insertion` 从写出的表渲染，
表由 :class:`~alnkit.core.tables.COLUMNS` 注册（``insertion_clip_reads`` /
``insertion_clip_consensus``）。

只认 CIGAR 首尾的 ``S``（soft clip）：``H`` 与 ``P`` 不产碱基，先剔除再判断首尾。
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pysam

from ..core import coords, tables
from ..core.errors import RunError, UsageError

__all__ = [
    "MAX_BREAKPOINTS",
    "MIN_SUPPORT_READS",
    "ClipGroup",
    "ClipRead",
    "build_clip_consensus",
    "build_clip_group",
    "build_groups",
    "breakpoint_from_1based",
    "identify_soft_clip",
    "load_clip_reads",
    "read_clip_table",
    "validate_breakpoints",
    "write_clip_table",
    "write_consensus_table",
]

#: 一致性序列每列至少要这么多条 read 支持，不足即中断
MIN_SUPPORT_READS = 2
#: 一张图最多画多少个断点（多断点图；再多会挤成一团）
MAX_BREAKPOINTS = 3

### 数据结构 ###


@dataclass(frozen=True, slots=True)
class ClipRead:
    """一条带 soft clip 的主比对（只保留画图与投票需要的字段）。

    输入示例（来自 BAM）::

        read_name="readA", ref_start=100000, ref_end=100100, read_length=160,
        strand="+", mapq=60, clip_side="left", clip_length=60,
        clip_sequence="GCATTT…"(60 bp), breakpoint=100000

    输出示例（``breakpoint`` 由 clip 侧决定）: left → 比对起点；right → 比对终点。
    """

    read_name: str
    ref_start: int
    ref_end: int
    read_length: int
    strand: str
    mapq: int
    clip_side: str
    clip_length: int
    clip_sequence: str
    breakpoint: int


@dataclass(frozen=True, slots=True)
class ClipGroup:
    """一个断点某一侧的全部 read 与它们的一致性序列。

    输入示例::

        breakpoint=100000, side="left", reads=(ClipRead, …)

    输出示例（由 :func:`build_clip_group` 生成）::

        consensus="GCATTT…"(60 bp), consensus_fraction=1.0, stop_reason="all columns"
    """

    breakpoint: int
    side: str
    reads: tuple[ClipRead, ...]
    consensus: str
    consensus_fraction: float
    stop_reason: str

    @property
    def maximum_clip_length(self) -> int:
        """该侧最长的 clip 长度（决定块宽与 read 行宽）。

        输入示例: 各 read 的 clip_length = [60, 58, 55]  →  输出示例: 60
        """
        return max((read.clip_length for read in self.reads), default=0)


### 读取 BAM：识别 soft clip ###


def identify_soft_clip(alignment: pysam.AlignedSegment) -> tuple[str, int, str]:
    """判断一条比对是"左 clip / 右 clip / 无 clip"，并取出 clip 碱基。

    只认 soft clip（CIGAR op=4）；硬剪贴（op=5）没有碱基，不参与一致性序列。

    输入示例::
        CIGAR "60S100M"  → ("left", 60, <序列前 60 bp>)
        CIGAR "100M55S"  → ("right", 55, <序列后 55 bp>)
        CIGAR "100M7H"   → ("none", 0, "")
    输出示例: ("left", 60, "GCATTT…")，序列已统一大写、与参考同向

    边界约定:
        * 首尾同时有 S 时只取较长的一侧（并列取 left），短侧证据被丢弃；
        * 硬剪贴 H 与 padding P 先被剔除，所以 "7H60S100M" 仍能识别出 60 bp 左 clip；
        * 整条 read 都是 clip（只剩 1 个 op）视为异常记录，返回 ("none", 0, "")。
    """
    if not alignment.cigartuples:
        return "none", 0, ""
    operations = [item for item in alignment.cigartuples if item[0] not in (5, 6)]
    if len(operations) < 2:  # 整条 read 都是 clip / 全是 H：异常记录，不当证据
        return "none", 0, ""

    first_operation, first_length = operations[0]
    last_operation, last_length = operations[-1]
    left_clip_length = first_length if first_operation == 4 else 0
    right_clip_length = last_length if last_operation == 4 else 0
    sequence = (alignment.query_sequence or "").upper()

    if left_clip_length and left_clip_length >= right_clip_length:
        return "left", left_clip_length, sequence[:left_clip_length]
    if right_clip_length:
        return "right", right_clip_length, sequence[len(sequence) - right_clip_length :]
    return "none", 0, ""


def _reference_names(alignment_file: pysam.AlignmentFile, limit: int = 5) -> str:
    """参考序列名预览，用于"名字打错"的报错。"""
    names = list(alignment_file.references or ())
    preview = ", ".join(names[:limit])
    return f"{preview} …（共 {len(names)} 条）" if len(names) > limit else preview


def load_clip_reads(
    bam_path: str | Path, chrom: str, span_start: int, span_end: int
) -> list[ClipRead]:
    """读取区间内带 soft clip 的 primary 比对，组装成 :class:`ClipRead` 列表。

    输入示例: bam_path=reads.bam, chrom="chr1", span_start=99200, span_end=100800
    输出示例::

        [ClipRead(read_name="readA", ref_start=99900, ref_end=100000, clip_side="right",
                  clip_length=55, clip_sequence="AGCC…", breakpoint=100000), …]

    说明:
        * 只取 primary 比对（secondary / supplementary 会把同一条 read 数两遍）；
        * 不按 MAPQ 与 duplicate 过滤（保持与上游脚本一致的口径）；
        * 窗口向外扩 1 bp，避免断点正好落在窗口边界的 read 被漏掉。
    """
    path = Path(bam_path)
    if not path.exists():
        raise UsageError(f"--bam 不存在：{path}")
    try:
        alignment_file = pysam.AlignmentFile(str(path))
    except (OSError, ValueError) as error:
        raise RunError(f"打不开 BAM：{path}（{error}）") from error

    with alignment_file:
        if not alignment_file.has_index():
            raise UsageError(f"{path} 没有索引：先跑 samtools index {path}")
        try:  # 先把窗口内的比对取出来：参考名打错时 pysam 在这里抛 ValueError
            alignments = list(alignment_file.fetch(chrom, max(0, span_start - 1), span_end + 1))
        except ValueError as error:
            raise UsageError(
                f"--chrom {chrom} 不在 {path} 里；可用：{_reference_names(alignment_file)}"
            ) from error

        clip_reads: list[ClipRead] = []
        for alignment in alignments:
            if alignment.is_unmapped or alignment.is_secondary or alignment.is_supplementary:
                continue
            if alignment.reference_start is None or alignment.reference_end is None:
                continue
            clip_side, clip_length, clip_sequence = identify_soft_clip(alignment)
            if clip_side == "none":
                continue
            clip_reads.append(
                ClipRead(
                    read_name=alignment.query_name,
                    ref_start=alignment.reference_start,
                    ref_end=alignment.reference_end,
                    read_length=alignment.query_length or 0,
                    strand="-" if alignment.is_reverse else "+",
                    mapq=alignment.mapping_quality,
                    clip_side=clip_side,
                    clip_length=clip_length,
                    clip_sequence=clip_sequence,
                    breakpoint=(
                        alignment.reference_start
                        if clip_side == "left"
                        else alignment.reference_end
                    ),
                )
            )
    return clip_reads


### 一致性序列：逐列投票 ###


def build_clip_consensus(
    clip_reads: Sequence[ClipRead], clip_side: str, min_support: int = MIN_SUPPORT_READS
) -> tuple[str, list[int], int, float, str]:
    """把同侧 soft clip 按"贴比对的那一端"对齐后逐列投票，得到一致性序列。

    规则:
        1. 每列至少 ``min_support`` 条 read 支持，不足则中断（``coverage<2``）；
        2. 3 条以上重叠时取多数碱基，少数差异被压过；
        3. 只有 2 条重叠且不一致时立刻中断（``2-read overlap mismatch``）；
        4. 多数并列时按碱基字母定序，保证同一输入永远得到同一结果。

    输入示例: 15 条左 clip（其中 14 条 60 bp、一条 10 bp），clip_side="left"
    输出示例: ("GCATTT…"(60 bp), [15, 15, …], 15, 1.0, "all columns")
              第 11 列起只剩 14 条 read，仍 >= min_support，所以一路走完；
              若长 read 只剩 1 条，则会在第 11 列给出 "coverage<2 (only 1 reads)"。

    返回: (一致性序列, 每列支持深度, 参与 read 数, 一致率, 中断原因)
    """
    oriented_sequences = [
        (read.clip_sequence[::-1] if clip_side == "left" else read.clip_sequence)
        for read in clip_reads
        if read.clip_side == clip_side and read.clip_sequence
    ]
    if not oriented_sequences:
        return "", [], 0, 0.0, "no clip"

    longest = max(len(sequence) for sequence in oriented_sequences)
    consensus_bases: list[str] = []
    column_depths: list[int] = []
    matched_bases = 0
    total_bases = 0
    stop_reason = "all columns"

    for column in range(longest):
        column_bases = [
            sequence[column]
            for sequence in oriented_sequences
            if column < len(sequence) and sequence[column] != "N"
        ]
        if len(column_bases) < min_support:
            stop_reason = f"coverage<{min_support} (only {len(column_bases)} reads)"
            break
        base_counts = Counter(column_bases)
        if len(column_bases) == 2 and len(base_counts) > 1:
            stop_reason = "2-read overlap mismatch"
            break
        winner, winner_count = sorted(base_counts.items(), key=lambda item: (-item[1], item[0]))[0]
        consensus_bases.append(winner)
        column_depths.append(sum(base_counts.values()))
        matched_bases += winner_count
        total_bases += sum(base_counts.values())

    consensus = "".join(consensus_bases)
    fraction = matched_bases / total_bases if total_bases else 0.0
    return consensus, column_depths, len(oriented_sequences), fraction, stop_reason


def build_clip_group(
    clip_reads: Iterable[ClipRead], breakpoint: int, clip_side: str
) -> ClipGroup | None:
    """取出某个断点某一侧的 read 并算好一致性序列；该侧没有 read 时返回 None。

    输入示例: clip_reads=[…], breakpoint=100000, clip_side="left"
    输出示例: ClipGroup(breakpoint=100000, side="left", reads=(…15 条…),
                        consensus="GCATTT…", consensus_fraction=1.0,
                        stop_reason="all columns")
    """
    side_reads = sorted(
        (
            read
            for read in clip_reads
            if read.breakpoint == breakpoint and read.clip_side == clip_side
        ),
        key=lambda read: (-read.clip_length, read.read_name),
    )
    if not side_reads:
        return None
    consensus, _depths, _count, fraction, stop_reason = build_clip_consensus(side_reads, clip_side)
    return ClipGroup(
        breakpoint=breakpoint,
        side=clip_side,
        reads=tuple(side_reads),
        consensus=consensus,
        consensus_fraction=fraction,
        stop_reason=stop_reason,
    )


def build_groups(clip_reads: Iterable[ClipRead], breakpoints: Sequence[int]) -> list[ClipGroup]:
    """把 read 按（断点, 侧）分组，顺序固定为断点升序、左侧在前。

    输入示例: clip_reads=[…], breakpoints=[700000, 700230]
    输出示例: [ClipGroup(700000,"left",…), ClipGroup(700230,"right",…)]
              （某侧没有 read 就不出现）
    """
    groups: list[ClipGroup] = []
    for breakpoint in sorted(breakpoints):
        for clip_side in ("left", "right"):
            group = build_clip_group(clip_reads, breakpoint, clip_side)
            if group is not None:
                groups.append(group)
    return groups


### 参数校验 ###


def breakpoint_from_1based(breakpoint_1based: int) -> int:
    """把 ``--breakpoint`` 的 1-based 输入换成内部的 0-based 断点边界。

    输入示例: breakpoint_from_1based(100000) → 输出示例: 100000

    说明: 单点坐标的 0-based 边界值与“断点左侧碱基的 1-based 坐标”相同
          （即 ``core.coords.left_base_1based``），所以这里是恒等映射。单独写成函数
          是为了让调用处一眼看出“没有做 ±1”，并拦住 <= 0 的输入。
    """
    if breakpoint_1based < 1:
        raise UsageError(f"--breakpoint 是 1-based 坐标，必须 >= 1：{breakpoint_1based}")
    return coords.left_base_1based(breakpoint_1based)


def validate_breakpoints(
    breakpoints: Sequence[int], clip_reads: Sequence[ClipRead]
) -> tuple[int, ...]:
    """校验断点参数：去重、限制数量，并要求每个断点在数据里确实有 clip read。

    输入示例: breakpoints=[100000, 100000], clip_reads=[…breakpoint=100000…]
    输出示例: (100000,)（重复值已合并）

    异常: 超过 :data:`MAX_BREAKPOINTS` 个、或某断点没有 read 时抛 ``UsageError``，
          后者会列出本区间可用的断点（1-based）。
    """
    unique_breakpoints = tuple(sorted(set(breakpoints)))
    if not unique_breakpoints:
        raise UsageError("--breakpoint 至少要给 1 个")
    if len(unique_breakpoints) > MAX_BREAKPOINTS:
        raise UsageError(f"--breakpoint 最多给 {MAX_BREAKPOINTS} 个（推荐 2 个）")

    available = sorted({read.breakpoint for read in clip_reads})
    missing = [bp for bp in unique_breakpoints if bp not in available]
    if missing:
        listed = "、".join(str(coords.left_base_1based(bp)) for bp in available[:10])
        suffix = f" …（共 {len(available)} 个）" if len(available) > 10 else ""
        raise UsageError(
            f"断点 {[coords.left_base_1based(bp) for bp in missing]} 在本区间没有任何 clip read；"
            f"可用的断点(1-based)：{listed}{suffix}"
        )
    return unique_breakpoints


### 写表 / 读表 ###


def clip_read_rows(clip_reads: Iterable[ClipRead], chrom: str, sample: str) -> list[dict[str, Any]]:
    """把 ``ClipRead`` 转成 ``insertion_clip_reads`` 的行（列序由表契约决定）。

    输入示例: clip_reads=[ClipRead(ref_start=99940, clip_side="left", …)], chrom="chr1"
    输出示例::

        [{"sample": "sampleA", "chrom": "chr1", "breakpoint_0based": 100000,
          "clip_side": "left", "read_name": "readA", "clip_length_bp": 60, …}]
    """
    return [
        {
            "sample": sample,
            "chrom": chrom,
            "breakpoint_0based": read.breakpoint,
            "clip_side": read.clip_side,
            "read_name": read.read_name,
            "read_length": read.read_length,
            "strand": read.strand,
            "mapq": read.mapq,
            "ref_start_0based": read.ref_start,
            "ref_end_0based": read.ref_end,
            "clip_length_bp": read.clip_length,
            "clip_sequence": read.clip_sequence or "-",
        }
        for read in clip_reads
    ]


def consensus_rows(groups: Iterable[ClipGroup], chrom: str, sample: str) -> list[dict[str, Any]]:
    """把 ``ClipGroup`` 转成 ``insertion_clip_consensus`` 的行。

    输入示例: groups=[ClipGroup(breakpoint=100000, side="left", consensus="GCAT…", …)]
    输出示例::

        [{"sample": "sampleA", "chrom": "chr1", "breakpoint_0based": 100000,
          "breakpoint_1based": 100000, "clip_side": "left", "n_reads": 15,
          "consensus_bp": 60, "consensus_fraction": 1.0, "stop_reason": "all columns",
          "consensus_sequence": "GCAT…"}]

    说明: ``breakpoint_1based`` 按 ``core/coords.py`` 的单点口径取“左侧碱基的 1-based
    坐标”，数值与 0-based 边界相同（不是边界 + 1）。
    """
    return [
        {
            "sample": sample,
            "chrom": chrom,
            "breakpoint_0based": group.breakpoint,
            "breakpoint_1based": coords.left_base_1based(group.breakpoint),
            "clip_side": group.side,
            "n_reads": len(group.reads),
            "consensus_bp": len(group.consensus),
            "consensus_fraction": group.consensus_fraction,
            "stop_reason": group.stop_reason,
            "consensus_sequence": group.consensus or "-",
        }
        for group in groups
    ]


def write_clip_table(
    path: str | Path,
    clip_reads: Iterable[ClipRead],
    chrom: str,
    sample: str,
    comments: Sequence[str] = (),
) -> int:
    """写 ``insertion_clip_reads``（原子写），返回行数。

    输入示例: clip_reads=30 条, path=out/demo.clips.tsv
    输出示例: 生成 12 列表格，返回 30
    """
    rows = clip_read_rows(clip_reads, chrom, sample)
    return tables.write_tsv(
        path,
        "insertion_clip_reads",
        rows,
        comments=comments,
        sort_key=lambda row: (row["breakpoint_0based"], row["clip_side"], row["read_name"]),
    )


def write_consensus_table(
    path: str | Path,
    groups: Iterable[ClipGroup],
    chrom: str,
    sample: str,
    comments: Sequence[str] = (),
) -> int:
    """写 ``insertion_clip_consensus``（原子写），返回行数。

    输入示例: groups=3 组, path=out/demo.consensus.tsv
    输出示例: 生成 10 列表格，返回 3
    """
    rows = consensus_rows(groups, chrom, sample)
    return tables.write_tsv(
        path,
        "insertion_clip_consensus",
        rows,
        comments=comments,
        sort_key=lambda row: (row["breakpoint_0based"], row["clip_side"]),
    )


def read_clip_table(path: str | Path) -> tuple[list[ClipRead], str]:
    """读 ``insertion_clip_reads`` 表，返回 (逐 read 证据, 参考序列名)。

    输入示例: path=out/demo.clips.tsv
    输出示例: ([ClipRead(ref_start=99940, clip_side="left", …), …], "chr1")

    异常: 表头不是注册列序、或表里没有行时抛 ``RunError``（退出码 1）。
    """
    header, rows = tables.read_tsv(path)
    if tuple(header) != tables.columns_of("insertion_clip_reads"):
        raise RunError(f"{path} 的列与注册表不一致：{header}")
    if not rows:
        raise RunError(f"{path} 里没有任何 clip read")

    clip_reads = [
        ClipRead(
            read_name=str(row["read_name"]),
            ref_start=int(row["ref_start_0based"]),
            ref_end=int(row["ref_end_0based"]),
            read_length=int(row["read_length"]),
            strand=str(row["strand"]),
            mapq=int(row["mapq"]),
            clip_side=str(row["clip_side"]),
            clip_length=int(row["clip_length_bp"]),
            clip_sequence="" if row["clip_sequence"] == "-" else str(row["clip_sequence"]),
            breakpoint=int(row["breakpoint_0based"]),
        )
        for row in rows
    ]
    return clip_reads, str(rows[0]["chrom"])
