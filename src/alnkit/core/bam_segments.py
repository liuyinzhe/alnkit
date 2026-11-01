"""BAM -> 规范 ``Segment``（D2）：CIGAR 解析、简化视图、read 轴归正。

来源：``upstream/scripts/alignment_table.py``（276 行）的移植。与移植前保持三件事不变：

* 参考坐标 0-based 半开，取 pysam 的 ``reference_start/reference_end``；
* read 坐标统一到**正向 read 轴**（反向记录做 ``read_length - stored`` 镜像），
  这样同一条 read 的基因组段与外源段才能直接比距离；
* ``cigar`` 是**简化视图**（两个 M 之间总长 <= 5 bp 的 I/D 并入 M），
  因此它的 M 之和可能小于 ``ref_aln_len`` 的真实跨度。

历史表头（``ref_pos_start`` 等 14 列）不在这里写，见 ``alnkit.cli.seg`` 的 legacy 写出器。
"""

from __future__ import annotations

import os

import pysam

from .coords import flip_read_interval
from .errors import UsageError
from .segments import Segment

__all__ = [
    "ABSORBABLE_GAP_MAX_LENGTH",
    "cigar_ops_from_pysam",
    "format_cigar_ops",
    "load_alignment_segments",
    "load_fai_lengths",
    "simplify_cigar_ops",
]

#: 被吸进 M 的小 I/D 总长上限
ABSORBABLE_GAP_MAX_LENGTH = 5

_OPERATION_BY_CODE = {0: "M", 1: "I", 2: "D", 3: "N", 4: "S", 5: "H", 6: "P", 7: "=", 8: "X"}


def cigar_ops_from_pysam(read: pysam.AlignedSegment) -> list[tuple[int, str]]:
    """pysam ``cigartuples`` -> ``[(长度, 操作符)]``。"""
    return [(length, _OPERATION_BY_CODE[code]) for code, length in (read.cigartuples or [])]


def simplify_cigar_ops(cigar_ops: list[tuple[int, str]]) -> list[tuple[int, str]]:
    """把“夹在两个 M 之间、总长 <= 5”的 I/D 并进 M。

    例：``[(39,'M'),(1,'D'),(66,'M'),(1,'D'),(92,'M'),(1,'I'),(445,'M')]`` -> ``[(645,'M')]``。
    S/H/N/P 与首尾的 I/D 一律不动；相邻同类型操作合并。
    """
    match_length = 0
    pending_ops: list[tuple[int, str]] = []
    pending_bases = 0
    simplified: list[tuple[int, str]] = []

    def close_match_run() -> None:
        nonlocal match_length, pending_ops, pending_bases
        if match_length:
            simplified.append((match_length, "M"))
            match_length = 0
        if pending_ops:
            simplified.extend(pending_ops)
            pending_ops = []
            pending_bases = 0

    for length, operation in cigar_ops:
        if operation in ("M", "=", "X"):
            if match_length:
                match_length += pending_bases + length
            else:
                close_match_run()
                match_length = length
            pending_ops = []
            pending_bases = 0
            continue
        if (
            operation in ("I", "D")
            and match_length
            and pending_bases + length <= ABSORBABLE_GAP_MAX_LENGTH
        ):
            pending_ops.append((length, operation))
            pending_bases += length
            continue
        close_match_run()
        simplified.append((length, operation))

    close_match_run()
    merged: list[tuple[int, str]] = []
    for length, operation in simplified:
        if merged and merged[-1][1] == operation:
            merged[-1] = (merged[-1][0] + length, operation)
        else:
            merged.append((length, operation))
    return merged


def format_cigar_ops(cigar_ops: list[tuple[int, str]]) -> str:
    """``[(418,'S'),(1358,'M'),(10151,'S')]`` -> ``"418S1358M10151S"``。"""
    return "".join(f"{length}{operation}" for length, operation in cigar_ops)


def load_fai_lengths(fai_path: str | os.PathLike[str]) -> dict[str, int]:
    """读 ``.fai``，返回 ``{序列名: 长度}``。"""
    lengths: dict[str, int] = {}
    try:
        with open(fai_path, encoding="utf-8") as handle:
            for line in handle:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    lengths[fields[0]] = int(fields[1])
    except FileNotFoundError:
        raise UsageError(f"找不到 .fai：{os.fspath(fai_path)}（用 samtools faidx 生成）") from None
    return lengths


def load_alignment_segments(
    bam_path: str | os.PathLike[str],
    chrom_lengths: dict[str, int] | None = None,
    selected_read_names: set[str] | None = None,
    skip_secondary: bool = False,
) -> list[Segment]:
    """读 BAM 的每条比对，转成规范 ``Segment``。

    ``selected_read_names`` 给定时只解析这些 read（“先找目标 read、再只处理目标”，
    大 BAM 上省内存）；``skip_secondary`` 在解析 CIGAR 之前就跳过 secondary 记录。
    """
    if not os.path.exists(bam_path):
        raise UsageError(f"找不到 BAM：{os.fspath(bam_path)}")

    segments: list[Segment] = []
    with pysam.AlignmentFile(os.fspath(bam_path), "rb") as bam_file:
        for read in bam_file:
            if read.is_unmapped:
                continue
            if skip_secondary and read.is_secondary:
                continue
            if selected_read_names is not None and read.query_name not in selected_read_names:
                continue
            # 优先 infer_read_length：含硬剪切，否则反向记录的镜像换算会偏
            read_length = read.infer_read_length() or read.query_length
            if not read_length:
                # 无 SEQ 的 BAM（samtools view --no-seq）算不出 read 轴坐标，宁可报错也不写 0
                raise UsageError(
                    f"{read.query_name}: 无法确定 read 长度（BAM 缺 SEQ？用原始 BAM 重跑）"
                )
            cigar_ops = cigar_ops_from_pysam(read)
            match_columns = sum(length for length, op in cigar_ops if op in ("M", "=", "X"))
            inserted_bases = sum(length for length, op in cigar_ops if op == "I")
            deleted_bases = sum(length for length, op in cigar_ops if op == "D")
            alignment_start = read.query_alignment_start or 0
            alignment_end = read.query_alignment_end or 0
            if read.is_reverse:
                read_start, read_end = flip_read_interval(
                    read_length, alignment_start, alignment_end
                )
            else:
                read_start, read_end = alignment_start, alignment_end
            nm_value = read.get_tag("NM") if read.has_tag("NM") else None
            segments.append(
                Segment(
                    read_name=read.query_name or "",
                    chrom=read.reference_name or "",
                    ref_start_0based=read.reference_start,
                    ref_end_0based=read.reference_end or read.reference_start,
                    read_start_0based=read_start,
                    read_end_0based=read_end,
                    strand="-" if read.is_reverse else "+",
                    mapq=read.mapping_quality,
                    read_length=read_length,
                    chrom_length=(chrom_lengths or {}).get(read.reference_name or "", None),
                    ref_aln_len=(read.reference_end or 0) - read.reference_start,
                    read_aln_len=alignment_end - alignment_start,
                    cigar=format_cigar_ops(simplify_cigar_ops(cigar_ops)),
                    is_secondary=bool(read.is_secondary),
                    is_supplementary=bool(read.is_supplementary),
                    nm=int(nm_value) if isinstance(nm_value, (int, float)) else None,
                    match_columns=match_columns,
                    inserted_bases=inserted_bases,
                    deleted_bases=deleted_bases,
                )
            )
    return segments
