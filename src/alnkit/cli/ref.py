"""``alnkit ref <merge>``：参考准备（A1）。"""

from __future__ import annotations

import argparse

from .. import prepare
from ._common import GLOBAL_OPTIONS


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    ref = subparsers.add_parser(
        "ref",
        help="参考准备（合并外源序列）",
        description="把参考基因组与外源序列合成一条合并参考，供后续比对使用。",
    )
    sub = ref.add_subparsers(dest="ref_command", metavar="<subcommand>")
    sub.required = True

    parser = sub.add_parser(
        "merge",
        parents=[GLOBAL_OPTIONS],
        help="参考 + 外源 -> merged_reference.fa（含 .fai 与 header 表）",
        description=(
            "外源按顺序转大写后拼到参考末尾；同时输出 reference.fa(.fai)、"
            "reference_headers.tsv、exogenous_uppercased.fa。"
        ),
        epilog=(
            "例子：\n  alnkit ref merge --reference ref.fa.gz \\\n"
            "          --exog exog/insertA.fa --exog exog/insertB.fa --outdir 01_merge_ref"
        ),
    )
    parser.add_argument("--reference", required=True, help="参考基因组 FASTA（可 .gz）")
    parser.add_argument(
        "--exog",
        action="append",
        default=[],
        required=True,
        help="外源 FASTA，可重复（按给出顺序拼接）",
    )
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.set_defaults(handler=run_merge)


def run_merge(args: argparse.Namespace) -> int:
    result = prepare.merge_reference(
        prepare.MergeSettings(reference=args.reference, exog=tuple(args.exog), outdir=args.outdir)
    )
    print(f"完成：{result['sequences']} 条序列 -> {result['merged']}")
    return 0
