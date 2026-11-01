"""``alnkit align <reads>``：把 reads 比对到合并参考（A2）。"""

from __future__ import annotations

import argparse

from .. import prepare
from ._common import GLOBAL_OPTIONS


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    align = subparsers.add_parser(
        "align",
        help="比对（reads -> 合并参考）",
        description="调用 minimap2 + samtools 的薄封装；需要它们装在 PATH 里。",
    )
    sub = align.add_subparsers(dest="align_command", metavar="<subcommand>")
    sub.required = True

    parser = sub.add_parser(
        "reads",
        parents=[GLOBAL_OPTIONS],
        help="长读长 reads -> sorted BAM（-ax map-hifi -Y -N 50）",
        description=(
            "等价于旧 shell：minimap2 建 .mmi -> minimap2 -ax map-hifi -Y -N 50 -> "
            "samtools sort -> samtools index。"
        ),
        epilog=(
            "例子：\n  alnkit align reads --reference 01_merge_ref/merged_reference.fa \\\n"
            "          --reads input.fastq --outdir 02_align --threads 12"
        ),
    )
    parser.add_argument("--reference", required=True, help="合并参考 FASTA（需 .fai）")
    parser.add_argument("--reads", required=True, help="原始 reads（FASTQ/FASTA）")
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--name", default="reads_to_merged", help="输出 BAM 文件名前缀")
    parser.add_argument("--threads", type=int, default=12, help="线程数（默认 12）")
    parser.set_defaults(handler=run_align)


def run_align(args: argparse.Namespace) -> int:
    result = prepare.align_reads(
        prepare.AlignSettings(
            reference=args.reference,
            reads=args.reads,
            outdir=args.outdir,
            name=args.name,
            threads=args.threads,
        )
    )
    print(f"完成：{result['bam']}")
    return 0
