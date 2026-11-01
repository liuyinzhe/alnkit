"""绘图子命令：``alnkit plot <ribbon|tracks>``。

参数名与旧的单文件脚本保持一一对应，迁移期只需换命令名。
"""

from __future__ import annotations

import argparse

from ..core.errors import UsageError
from ._common import GLOBAL_OPTIONS


def _require_plot_dependency():
    """绘图子命令才需要 matplotlib；缺了要给可行动提示而不是 traceback。

    基础安装（``pip install alnkit``）不含 matplotlib，所以这里必须惰性导入，
    否则连 ``alnkit --version`` 都会在 import 阶段就崩（0.1.0 的真实缺陷）。
    """
    try:
        from ..plot import ribbon, tracks
    except ImportError as exc:  # pragma: no cover - 只在缺 matplotlib 时触发
        raise UsageError(
            f'绘图子命令需要 matplotlib（{exc}）。安装：pip install "alnkit[plot]"'
        ) from exc
    return ribbon, tracks


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    plot = subparsers.add_parser(
        "plot",
        help="绘图（丝带图 / 单 read 轨道图）",
        description="静态图输出：PNG + SVG + PDF；不依赖任何上游步骤产物。",
    )
    plot_subparsers = plot.add_subparsers(dest="plot_command", metavar="<subcommand>")
    plot_subparsers.required = True
    _add_tracks(plot_subparsers)
    _add_ribbon(plot_subparsers)


def _add_tracks(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "tracks",
        parents=[GLOBAL_OPTIONS],
        help="单条 read 的比对片段轨道图",
        description=(
            "一条 read 一行，横轴为 read 自身坐标；每个比对片段画成一个块。"
            "输入三选一：--bam / --paf / --segments-tsv。"
        ),
        epilog=(
            "例子：\n  alnkit plot tracks --segments-tsv seg.tsv --outdir out\n"
            "     alnkit plot tracks --bam reads.bam --read readA --exog insertA "
            "--highlight-ref chr1:1000001 --outdir out"
        ),
    )
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--bam", help="BAM（pysam 读取）")
    source.add_argument("--paf", help="minimap2 PAF（需要 -c）")
    source.add_argument("--segments-tsv", dest="segments_tsv", help="已整理的片段表")
    parser.add_argument("--read", action="append", default=[], help="只画这些 read；可重复")
    parser.add_argument(
        "--reads-file", dest="reads_file", help="read 名列表文件（一行一个，空行/# 跳过）"
    )
    parser.add_argument(
        "--exog", action="append", default=[], help="外源序列名；可重复或用逗号分隔"
    )
    parser.add_argument(
        "--bam-mode",
        dest="bam_mode",
        choices=("auto", "scan", "subset"),
        default="auto",
        help="BAM 取数方式：auto（默认）/scan/subset",
    )
    parser.add_argument(
        "--min-mapq", dest="min_mapq", type=int, default=0, help="最小 MAPQ（默认 0）"
    )
    parser.add_argument(
        "--include-secondary",
        dest="include_secondary",
        action="store_true",
        help="包含 secondary 比对",
    )
    parser.add_argument(
        "--highlight-ref",
        dest="highlight_ref",
        action="append",
        default=[],
        help="要标注的参考断点，格式 chrom:pos；可重复",
    )
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--name", default="read_tracks", help="输出文件名前缀")
    parser.add_argument(
        "--title", default="Read-axis alignment segments", help="图标题（英文，避免缺字体）"
    )
    parser.add_argument(
        "--full-read-names",
        dest="full_read_names",
        action="store_true",
        help="左侧显示完整 read 名",
    )
    parser.add_argument("--dpi", type=int, default=150, help="PNG 分辨率（默认 150）")
    parser.add_argument(
        "--fig-width", dest="fig_width", type=float, default=14.0, help="图宽（英寸，默认 14）"
    )
    parser.add_argument(
        "--exact-highlight",
        dest="exact_highlight",
        action="store_true",
        help="参考断点按 CIGAR 精确映射（BAM 输入）",
    )
    parser.add_argument(
        "--label-mode",
        dest="label_mode",
        choices=("auto", "inline", "side"),
        default="auto",
        help="片段标签：auto（默认）/inline/side",
    )
    parser.add_argument(
        "--max-axis-bp",
        dest="max_axis_bp",
        type=int,
        default=None,
        help="超过该长度折叠内容空隙（断轴）",
    )
    parser.set_defaults(handler=run_tracks)


def run_tracks(args: argparse.Namespace) -> int:
    _ribbon_unused, tracks = _require_plot_dependency()
    settings = tracks.TrackSettings(
        outdir=args.outdir,
        name=args.name,
        bam=args.bam,
        paf=args.paf,
        segments_tsv=args.segments_tsv,
        read=tuple(args.read),
        reads_file=args.reads_file,
        exog=tuple(args.exog),
        bam_mode=args.bam_mode,
        min_mapq=args.min_mapq,
        include_secondary=args.include_secondary,
        highlight_ref=tuple(args.highlight_ref),
        title=args.title,
        full_read_names=args.full_read_names,
        dpi=args.dpi,
        fig_width=args.fig_width,
        exact_highlight=args.exact_highlight,
        label_mode=args.label_mode,
        max_axis_bp=args.max_axis_bp,
    )
    result = tracks.run(settings)
    print(f"reads: {result['reads']}；segments: {result['segments']}")
    for key, label in (
        ("png", "图"),
        ("svg", "SVG"),
        ("pdf", "PDF"),
        ("segments_tsv", "片段表"),
    ):
        print(f"{label}: {result[key]}")
    return 0


def _add_ribbon(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "ribbon",
        parents=[GLOBAL_OPTIONS],
        help="共线性丝带图（参考 vs 多个查询片段）",
        description=(
            "方式 1：给 --reference + --query，内部跑 minimap2；"
            "方式 2：给 --paf（现成比对），直接画图。"
        ),
        epilog=(
            "例子：\n  alnkit plot ribbon --paf aln.paf --reference-name assembly_a "
            "--outdir out --name ex4\n"
            '     alnkit plot ribbon --reference "ins=exog.fa:insertA" '
            '--query "C01=asm.fa:oe_contig1" --outdir out --name exog_vs_class'
        ),
    )
    parser.add_argument(
        "--reference", help="参考侧，格式 标签=fasta路径[:片段名]（--paf 模式可省略）"
    )
    parser.add_argument("--query", action="append", default=[], help="查询侧，格式同上；可重复")
    parser.add_argument(
        "--query-file", dest="query_file", help="查询清单文件（每行一个 标签=路径[:片段名]）"
    )
    parser.add_argument("--paf", help="现成 PAF（方式 2；需带 -c 才有匹配列）")
    parser.add_argument(
        "--reference-name", dest="reference_name", help="--paf 模式下指定哪条是参考行"
    )
    parser.add_argument(
        "--query-name",
        dest="query_name",
        action="append",
        default=[],
        help="只画这些查询片段；可重复",
    )
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--name", required=True, help="输出文件名前缀")
    parser.add_argument("--title", default="", help="图标题（英文，避免缺字体）")
    parser.add_argument(
        "--reference-position",
        dest="reference_position",
        choices=("top", "bottom"),
        default="top",
        help="参考行放顶部还是底部",
    )
    order = parser.add_mutually_exclusive_group()
    order.add_argument("--order", help="行顺序，用 > 分隔")
    order.add_argument("--order-file", dest="order_file", help="行顺序清单文件")
    parser.add_argument("--dpi", type=int, default=100, help="PNG 分辨率（默认 100）")
    parser.add_argument(
        "--min-label-bp",
        dest="min_label_bp",
        type=int,
        default=300,
        help="只标长度不小于它的块（默认 300）",
    )
    parser.add_argument(
        "--label-coordinate",
        dest="label_coordinate",
        choices=("query", "reference", "both"),
        default="query",
        help="标签用哪套坐标",
    )
    parser.add_argument(
        "--max-labels-per-row",
        dest="max_labels_per_row",
        type=int,
        default=8,
        help="每行最多标签数（默认 8）",
    )
    parser.add_argument("--placement-table", dest="placement_table", help="额外输出位置汇总表")
    parser.add_argument("--region-table", dest="region_table", help="额外输出同源区域表")
    parser.add_argument(
        "--min-identity",
        dest="min_identity",
        type=float,
        default=0.8,
        help="保留块的最小一致率（默认 0.8）",
    )
    parser.add_argument(
        "--min-block-bp",
        dest="min_block_bp",
        type=int,
        default=200,
        help="保留块的最小长度（默认 200）",
    )
    parser.add_argument(
        "--preset", default="asm20", help="minimap2 预设（默认 asm20；auto=按长度分档）"
    )
    parser.add_argument(
        "--query-preset",
        dest="query_preset",
        action="append",
        default=[],
        metavar="标签=预设[ 额外参数]",
        help="按查询标签覆盖比对参数；可重复",
    )
    parser.add_argument("--threads", type=int, default=8, help="线程数（默认 8）")
    parser.set_defaults(handler=run_ribbon)


def run_ribbon(args: argparse.Namespace) -> int:
    ribbon, _tracks_unused = _require_plot_dependency()
    settings = ribbon.RibbonSettings(
        outdir=args.outdir,
        name=args.name,
        reference=args.reference,
        query=tuple(args.query),
        query_file=args.query_file,
        paf=args.paf,
        reference_name=args.reference_name,
        query_name=tuple(args.query_name),
        title=args.title,
        reference_position=args.reference_position,
        order=args.order,
        order_file=args.order_file,
        dpi=args.dpi,
        min_label_bp=args.min_label_bp,
        label_coordinate=args.label_coordinate,
        max_labels_per_row=args.max_labels_per_row,
        placement_table=args.placement_table,
        region_table=args.region_table,
        min_identity=args.min_identity,
        min_block_bp=args.min_block_bp,
        preset=args.preset,
        query_preset=tuple(args.query_preset),
        threads=args.threads,
    )
    result = ribbon.run(settings)
    print(f"图: {result['png']}")
    print(f"映射块表: {result['blocks_tsv']}")
    return 0
