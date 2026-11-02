"""组装支撑度（A 方案）：把类内 reads 用 minimap2 打回 contig，量化“这个 contig 能不能用”。

为什么要这一步：`long assemble` 只保证断点附近的延伸，远端无法自证；而这个组装器
不输出逐 read 的布局表（坐标记账不可靠）。改用**比对回验**：reads 能不能打回 contig、
覆盖多少、深度有没有断层，是外部可复核的证据。

只用 minimap2（`-ax map-hifi -Y`）+ 纯 Python 统计，不引入 samtools/BAM 中转。
"""

from __future__ import annotations

import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import pysam

from ..core.coords import merge_intervals, to_0based
from ..core.errors import RunError
from ..core.external import find_tool
from ..core.logging import get_logger

__all__ = [
    "SupportRow",
    "class_reads_path",
    "collect_all_support",
    "collect_support",
    "coverage_stats",
    "find_minimap2",
    "write_support_report",
]

logger = get_logger("assemble.support")

#: 统计覆盖深度时用的窗口大小（bp）
WINDOW_BP = 200
#: 两端各取多少 bp 判断“末端是否被 read 支撑”
END_WINDOW_BP = 100
#: 判定 chimera 的深度断层阈值：相邻窗口深度降到前一个的 1/5 以下才算断层。
#: 定标依据（本仓基准数据，5 个类）：2× 会把 3/5 条 contig 误标（插入/载体边界天然有深度落差），
#: 3× 标 1/5，5× 在本数据上 0 误标但仍能标出 C03 的 2 处深断层。
CLIFF_RATIO = 0.2
#: 至少出现几个断层才算可疑
CLIFF_MIN_COUNT = 3
#: 低于这个覆盖比例算覆盖不足
MIN_COVERED_FRACTION = 0.9
#: 低于这个平均深度（×）算覆盖不足
MIN_MEAN_DEPTH = 2.0


@dataclass(slots=True, frozen=True)
class SupportRow:
    """一个 contig 的支撑度。"""

    class_id: str
    contig: str
    contig_length_bp: int
    n_reads_total: int
    n_reads_mapped: int
    covered_bp: int
    covered_fraction: float
    mean_depth: float
    depth_p10: float
    depth_cliffs: int
    ends_covered: str
    status: str


def coverage_stats(
    contig_length: int, intervals: list[tuple[int, int]], window_bp: int = WINDOW_BP
) -> tuple[int, float, float, list[float]]:
    """由 read 的参考区间算 (覆盖碱基数, 覆盖比例, 平均深度, 每窗口深度)。

    ``mean_depth`` 用“比对碱基总数 / contig 长度”，与 `samtools depth` 的全局均值同义。
    返回的每窗口深度已按窗口实际长度归一化成**倍数（×）**，所以最后那个不满宽的
    窗口不会因为容量小而假报断层。
    """
    if contig_length <= 0:
        return 0, 0.0, 0.0, []
    clipped = [(max(0, start), min(contig_length, end)) for start, end in intervals]
    clipped = [(start, end) for start, end in clipped if end > start]
    covered_bp = sum(end - start for start, end in merge_intervals(clipped))
    aligned_bp = sum(end - start for start, end in clipped)
    mean_depth = aligned_bp / contig_length
    n_windows = (contig_length + window_bp - 1) // window_bp
    counts = [0.0] * n_windows
    for start, end in clipped:
        first, last = start // window_bp, (end - 1) // window_bp
        for index in range(first, last + 1):
            left = max(start, index * window_bp)
            right = min(end, (index + 1) * window_bp)
            if right > left:
                counts[index] += right - left
    depths = [
        counts[index] / min(window_bp, contig_length - index * window_bp)
        for index in range(n_windows)
    ]
    return covered_bp, covered_bp / contig_length, mean_depth, depths


def _count_cliffs(depths: list[float]) -> int:
    """深度断层：某窗口深度 < 前一个窗口的 ``CLIFF_RATIO`` 倍，且前一个窗口至少 2× 覆盖。"""
    cliffs = 0
    for previous, current in zip(depths, depths[1:]):
        if previous >= 2 and current < previous * CLIFF_RATIO:
            cliffs += 1
    return cliffs


def _percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    index = min(len(ordered) - 1, max(0, int(len(ordered) * fraction)))
    return ordered[index]


def _parse_contigs(fasta_path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    name = None
    for line in fasta_path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            lengths[name] = 0
        elif name:
            lengths[name] += len(line.strip())
    return lengths


def _parse_sam_line(line: str) -> tuple[str, str, int, int, int] | None:
    """从 SAM 行取 (read 名, 参考名, 起点, 参考跨度, MAPQ)。

    跳过未比对（0x4）与 secondary（0x100）记录；保留 supplementary（0x800），
    因为 `-Y` 下嵌合 read 的远段就是独立的 supplementary 记录。
    跨度只累加 M/D/N/=/X（软/硬剪贴与插入不占参考）。
    """
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 10 or fields[0].startswith("@"):
        return None
    flag = int(fields[1])
    if flag & 0x4 or flag & 0x100:
        return None
    mapq = int(fields[4])
    span = 0
    for length, operation in re.findall(r"(\d+)([MIDNSHP=X])", fields[5]):
        if operation in ("M", "D", "N", "=", "X"):
            span += int(length)
    if span <= 0:
        return None
    start = to_0based(int(fields[3]))
    return fields[0], fields[2], start, span, mapq


def _run_minimap2(command: list[str], class_id: str) -> str:
    logger.debug("运行：%s", " ".join(command))
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        tail = (completed.stderr or "").strip().splitlines()[-5:]
        detail = "；".join(tail) if tail else "无 stderr"
        raise RunError(
            f"minimap2 回映失败（{class_id}）退出码 {completed.returncode}：{detail}；"
            + " ".join(command)
        )
    return completed.stdout


def collect_support(
    class_id: str,
    contigs_fasta: Path,
    reads_fastq: Path,
    minimap2: str,
    *,
    min_mapq: int = 10,
    window_bp: int = WINDOW_BP,
    threads: int = 4,
) -> list[SupportRow]:
    """把类内 reads 打回该类的 contig，产出每 contig 一行支撑度。"""
    contig_lengths = _parse_contigs(contigs_fasta)
    if not contig_lengths:
        return []
    command = [
        minimap2,
        "-ax",
        "map-hifi",
        "-Y",
        "-t",
        str(threads),
        str(contigs_fasta),
        str(reads_fastq),
    ]
    intervals: dict[str, list[tuple[int, int]]] = {name: [] for name in contig_lengths}
    mapped_reads: dict[str, set[str]] = {name: set() for name in contig_lengths}
    for line in _run_minimap2(command, class_id).splitlines():
        parsed = _parse_sam_line(line)
        if parsed is None:
            continue
        read_name, reference, start, span, mapq = parsed
        if mapq < min_mapq or reference not in contig_lengths:
            continue
        intervals[reference].append((start, start + span))
        mapped_reads[reference].add(read_name)

    n_reads_total = sum(1 for _ in pysam.FastxFile(str(reads_fastq)))
    rows: list[SupportRow] = []
    for contig, length in sorted(contig_lengths.items()):
        covered_bp, covered_fraction, mean_depth, depths = coverage_stats(
            length, intervals[contig], window_bp
        )
        cliffs = _count_cliffs(depths)
        p10 = _percentile([depth for depth in depths if depth > 0], 0.10)
        ends_covered = (
            "yes"
            if intervals[contig]
            and _covered_bp(intervals[contig], 0, min(END_WINDOW_BP, length)) > 0
            and _covered_bp(intervals[contig], max(0, length - END_WINDOW_BP), length) > 0
            else "no"
        )
        rows.append(
            SupportRow(
                class_id=class_id,
                contig=contig,
                contig_length_bp=length,
                n_reads_total=n_reads_total,
                n_reads_mapped=len(mapped_reads[contig]),
                covered_bp=covered_bp,
                covered_fraction=covered_fraction,
                mean_depth=mean_depth,
                depth_p10=p10,
                depth_cliffs=cliffs,
                ends_covered=ends_covered,
                status=_status(covered_fraction, mean_depth, cliffs, ends_covered),
            )
        )
    return rows


def _covered_bp(intervals: list[tuple[int, int]], start: int, end: int) -> int:
    return sum(max(0, min(end, right) - max(start, left)) for left, right in intervals)


def _status(covered_fraction: float, mean_depth: float, cliffs: int, ends: str) -> str:
    """判状态。阈值 ``MIN_COVERED_FRACTION`` / ``MIN_MEAN_DEPTH`` 见模块常量。

    输入示例: 1.0, 35.1, 0, "yes"  →  输出示例: "ok"
    """
    if covered_fraction <= 0:
        return "no_reads"
    if covered_fraction < MIN_COVERED_FRACTION or mean_depth < MIN_MEAN_DEPTH:
        return "low_coverage"
    if cliffs >= CLIFF_MIN_COUNT:
        return "chimera_suspect"
    if ends == "no":
        return "unsupported_ends"
    return "ok"


def write_support_report(rows: list[SupportRow], output_path: Path) -> int:
    """写 assembly_support.tsv（列序由 core.tables 的 `assembly_support` 定义）。"""
    from ..core import tables

    payload = [
        {
            "class_id": row.class_id,
            "contig": row.contig,
            "contig_length_bp": row.contig_length_bp,
            "n_reads_total": row.n_reads_total,
            "n_reads_mapped": row.n_reads_mapped,
            "covered_bp": row.covered_bp,
            "covered_fraction": round(row.covered_fraction, 4),
            "mean_depth": round(row.mean_depth, 4),
            "depth_p10": round(row.depth_p10, 2),
            "depth_cliffs": row.depth_cliffs,
            "ends_covered": row.ends_covered,
            "status": row.status,
        }
        for row in rows
    ]
    return tables.write_tsv(
        output_path,
        "assembly_support",
        payload,
        comments=[
            "把该类 reads 用 minimap2 -ax map-hifi -Y 打回 contig 得到的支撑度。",
            "status: ok / low_coverage / chimera_suspect / unsupported_ends / no_reads。",
            "mean_depth / depth_p10 单位是倍数（×）；depth_cliffs = 相邻窗口深度骤降次数。",
        ],
        sort_key=lambda row: (row["class_id"], row["contig"]),
    )


def find_minimap2() -> str | None:
    """找到 minimap2 的路径；没有就返回 None（调用方降级为不写报告）。"""
    found = find_tool("minimap2")
    return found.path if found else None


def class_reads_path(classes_dir: Path, class_id: str) -> Path:
    return classes_dir / class_id / "reads.fastq"


def iter_assembly_dirs(assembly_dir: Path, class_names=None):
    """枚举 `<assembly>/<class>/overlap_extension/contigs.fasta`；``class_names`` 可限定子集。"""
    if not assembly_dir.exists():
        return
    allowed = set(class_names) if class_names is not None else None
    for class_dir in sorted(assembly_dir.iterdir()):
        if not class_dir.is_dir() or (allowed is not None and class_dir.name not in allowed):
            continue
        contigs = class_dir / "overlap_extension" / "contigs.fasta"
        if contigs.exists():
            yield class_dir.name, contigs


def collect_all_support(
    assembly_dir: Path,
    classes_dir: Path,
    minimap2: str,
    *,
    threads: int = 4,
    class_names=None,
) -> list[SupportRow]:
    rows: list[SupportRow] = []
    for class_id, contigs in iter_assembly_dirs(assembly_dir, class_names):
        reads = class_reads_path(classes_dir, class_id)
        if not reads.exists():
            logger.warning("%s: 找不到 %s，跳过支撑度统计", class_id, os.fspath(reads))
            continue
        rows.extend(collect_support(class_id, contigs, reads, minimap2, threads=threads))
    return rows
