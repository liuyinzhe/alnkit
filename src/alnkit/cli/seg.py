"""``alnkit seg table``：BAM -> 比对段表（D2 的 CLI 出口）。

默认写**规范列**（带基准后缀）；``--columns legacy`` 写历史上那张 14 列表，
用于与旧脚本逐字节对账。
"""

from __future__ import annotations

import argparse
import os

from ..core import tables
from ..core.bam_segments import load_alignment_segments, load_fai_lengths
from ..core.errors import UsageError
from ..core.segments import Segment, segment_rows
from ._common import GLOBAL_OPTIONS

#: 历史 14 列表：对外列名 -> Segment 字段（顺序不可动，黄金对账依赖它）
LEGACY_COLUMNS: tuple[tuple[str, str], ...] = (
    ("chrom", "chrom"),
    ("chrom_length", "chrom_length"),
    ("ref_pos_start", "ref_start_0based"),
    ("ref_pos_end", "ref_end_0based"),
    ("read_name", "read_name"),
    ("read_pos_start", "read_start_0based"),
    ("read_pos_end", "read_end_0based"),
    ("strand", "strand"),
    ("mapping_quality", "mapq"),
    ("is_secondary", "is_secondary"),
    ("is_supplementary", "is_supplementary"),
    ("ref_aln_len", "ref_aln_len"),
    ("read_aln_len", "read_aln_len"),
    ("simplified_cigar", "cigar"),
)


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    seg = subparsers.add_parser(
        "seg",
        help="比对段表（BAM -> TSV）",
        description="把 BAM 的每条比对导成比对段表；坐标 0-based 半开，read 轴已翻正。",
    )
    seg_subparsers = seg.add_subparsers(dest="seg_command", metavar="<subcommand>")
    seg_subparsers.required = True
    _add_table(seg_subparsers)


def _add_table(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "table",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="BAM -> 比对段表",
        description="输出规范列（默认）或历史 14 列（--columns legacy）。",
        epilog=(
            "例子：\n  alnkit seg table --bam reads.sorted.bam --fai ref.fa.fai --out seg.tsv\n"
            "     alnkit seg table --bam reads.bam --out legacy.tsv --columns legacy"
        ),
    )
    parser.add_argument("--bam", required=True, help="比对 BAM（顺序扫描，不需要索引）")
    parser.add_argument("--out", required=True, help="输出 TSV")
    parser.add_argument("--fai", default="", help="参考序列的 .fai（提供后带 chrom_length 列）")
    parser.add_argument(
        "--columns",
        choices=("canonical", "legacy"),
        default="canonical",
        help="canonical（默认，带基准后缀）或 legacy（历史 14 列，逐字节兼容）",
    )
    parser.add_argument("--read-file", dest="read_file", help="只导出这些 read（一行一个）")
    parser.add_argument(
        "--skip-secondary", dest="skip_secondary", action="store_true", help="跳过 secondary 比对"
    )
    parser.set_defaults(handler=run_table)


def _format_legacy(segment: Segment) -> list[str]:
    values: list[str] = []
    for _, field in LEGACY_COLUMNS:
        value = getattr(segment, field)
        if value is None:
            values.append("-")
        elif isinstance(value, bool):
            values.append("1" if value else "0")
        else:
            values.append(str(value))
    return values


def write_legacy_table(segments: list[Segment], path: str) -> int:
    """写历史 14 列表（列名/列序/排序/取值与旧脚本一致）。"""
    ordered = sorted(segments, key=lambda s: (s.read_name, s.chrom, s.ref_start_0based))
    directory = os.path.dirname(os.path.abspath(path))
    os.makedirs(directory, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        handle.write("\t".join(name for name, _ in LEGACY_COLUMNS) + "\n")
        for segment in ordered:
            handle.write("\t".join(_format_legacy(segment)) + "\n")
    return len(ordered)


def run_table(args: argparse.Namespace) -> int:
    chrom_lengths = load_fai_lengths(args.fai) if args.fai else None
    selected = None
    if args.read_file:
        with open(args.read_file, encoding="utf-8") as handle:
            selected = {
                line.strip() for line in handle if line.strip() and not line.startswith("#")
            }
        if not selected:
            raise UsageError(f"{args.read_file} 里没有 read 名")

    segments = load_alignment_segments(
        args.bam,
        chrom_lengths=chrom_lengths,
        selected_read_names=selected,
        skip_secondary=args.skip_secondary,
    )
    if args.columns == "legacy":
        count = write_legacy_table(segments, args.out)
    else:
        count = tables.write_tsv(
            args.out,
            "segments",
            segment_rows(segments),
            sort_key=lambda row: (
                row["read_name"],
                row["chrom"],
                row["ref_start_0based"],
            ),
        )
    print(f"比对表 {args.out}：{count} 行（columns={args.columns}）")
    return 0
