"""绘图子命令：``alnkit plot <ribbon|tracks|insertion|coverage>``。

四条命令按**数据来源**分成两类：``ribbon`` / ``tracks`` 最初是为三代长读长写的
（一条 read 上多个比对片段、contig vs 参考），``insertion`` / ``coverage`` 是给二代
短读长比对用的（soft clip 断点、覆盖广度/深度）。参数名与原脚本保持一一对应。
"""

from __future__ import annotations

import argparse

from ..core import coords
from ..core.errors import UsageError
from ..shortread import insertion
from ._common import GLOBAL_OPTIONS, require_existing_files


def _require_plot_module(module_name: str):
    """取 ``alnkit.plot.<name>``；缺 matplotlib/numpy 时给可行动提示而不是 traceback。

    基础安装（``pip install alnkit``）不含 matplotlib，所以这里必须惰性导入，
    否则连 ``alnkit --version`` 都会在 import 阶段就崩（0.1.0 的真实缺陷）。
    """
    import importlib

    try:
        return importlib.import_module(f"alnkit.plot.{module_name}")
    except ImportError as exc:  # pragma: no cover - 只在缺 matplotlib 时触发
        raise UsageError(
            f'绘图子命令需要 matplotlib 与 numpy（{exc}）。安装：pip install "alnkit[plot]"'
        ) from exc


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    plot = subparsers.add_parser(
        "plot",
        help="绘图（丝带图 / read 轨道图 / 插入位点示意图 / 覆盖度）",
        description="静态图输出：PNG（tracks 另有 SVG/PDF）；不依赖任何上游步骤产物。",
    )
    plot_subparsers = plot.add_subparsers(dest="plot_command", metavar="<subcommand>")
    plot_subparsers.required = True
    _add_tracks(plot_subparsers)
    _add_ribbon(plot_subparsers)
    _add_insertion(plot_subparsers)
    _add_coverage(plot_subparsers)


def _add_tracks(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "tracks",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="long read 的比对片段轨道图",
        description=(
            "一条 read 一行，横轴为 read 自身坐标；每个比对片段画成一个块。"
            "输入三选一：--bam / --paf / --segments-tsv。\n"
            "设计初衷是长读长（三代）比对：一条 read 上有多个比对片段（外源 + 基因组 / "
            "contig 拼接）。二代 PE 也能画，但短读长的断点证据与覆盖度请看 "
            "`plot insertion` / `plot coverage`。"
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
        help="要标注的参考断点，格式 chrom:pos（pos 为 0-based，与表的坐标列一致）；可重复",
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
    require_existing_files(
        ("--bam", args.bam),
        ("--paf", args.paf),
        ("--segments-tsv", args.segments_tsv),
        ("--reads-file", args.reads_file),
    )
    tracks = _require_plot_module("tracks")
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
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="共线性丝带图（参考 vs 多个查询片段）",
        description=(
            "方式 1：给 --reference + --query，内部跑 minimap2；"
            "方式 2：给 --paf（现成比对），直接画图。\n"
            "设计初衷是长读长（三代）组装结果 vs 参考（contig 落在参考的哪些段、哪些是外源）；"
            "短读长的插入证据用 `plot insertion`、覆盖均匀性用 `plot coverage`。"
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
    order.add_argument("--order", help="行顺序，用逗号分隔（如 reference,insertA）")
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
    require_existing_files(
        ("--paf", args.paf), ("--order-file", args.order_file), ("--query-file", args.query_file)
    )
    ribbon = _require_plot_module("ribbon")
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


# ---------------------------------------------------------------- 二代的两种静图
# 数据/统计层在 `shortread/insertion.py` 与 `shortread/coverage.py`，这里只做 CLI 编排。
def _add_insertion(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "insertion",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="short reads 比对的插入位点示意图（单侧 / 同断点双侧 / 多断点）",
        description=(
            "适用于短读长（二代 / PE）比对：从 CIGAR 首尾的 soft clip 取断点证据，"
            "按“给了几个断点 + 每个断点有几侧 clip”自动选画法。长读长（三代）的比对片段请用 "
            "`plot tracks` / `plot ribbon`。出图只读 TSV：先写 clips/consensus 两张表，再由表渲染。"
        ),
        epilog=(
            "画法:\n"
            "  1 个断点、只有一侧 clip -> 单侧图（未知条 + 彩色块贴在断点一侧）\n"
            "  1 个断点、两侧都有 clip -> 同断点双侧图（左块 + 合并未知区 + 右块，连线收到断点）\n"
            "  2-3 个断点              -> 多断点图（每断点各自 [未知][块]，块间强制留间隔）\n"
            "\n配色: 左 clip 蓝 / 右 clip 紫 / 未知区灰斜线 / 断点红。\n"
            "\n例子:\n"
            "  alnkit plot insertion --bam reads.bam --region chr1:99201-100800 \\\n"
            "      --breakpoint 100000 --out-dir out --name demo\n"
            "  alnkit plot insertion --bam reads.bam --region chr1:699101-700900 \\\n"
            "      --breakpoint 700000 --breakpoint 700230 --out-dir out --name two\n"
            "\n产物:\n"
            "  <name>.png（示意图）+ <name>.clips.tsv（insertion_clip_reads）\n"
            "  + <name>.consensus.tsv（insertion_clip_consensus）"
        ),
    )
    parser.add_argument(
        "--bam", required=True, metavar="PATH", help="带 soft clip 的 BAM（需 .bai 索引）"
    )
    parser.add_argument(
        "--region",
        required=True,
        metavar="CHROM:START-END",
        help="窗口（1-based 闭区间），只画窗口内的 read",
    )
    parser.add_argument(
        "--breakpoint",
        action="append",
        required=True,
        default=[],
        type=int,
        metavar="POS",
        help=("断点（1-based，= 断点左侧碱基的坐标）；给 1 个=单断点图，2-3 个=多断点图；可重复"),
    )
    parser.add_argument("--out-dir", dest="out_dir", default=".", help="输出目录，默认当前目录")
    parser.add_argument("--name", default="insertion", help="输出文件名前缀，默认 insertion")
    parser.add_argument("--dpi", type=int, default=150, help="PNG 分辨率（默认 150）")
    parser.add_argument(
        "--no-consensus-text",
        dest="no_consensus_text",
        action="store_true",
        help="不在图底部写一致性序列文字",
    )
    parser.add_argument(
        "--no-formula",
        dest="no_formula",
        action="store_true",
        help="不在图右下角写 Consensus fraction 公式",
    )
    parser.set_defaults(handler=run_insertion_plot)


def run_insertion_plot(args: argparse.Namespace) -> int:
    insertion_plot = _require_plot_module("insertion")
    try:
        chrom, window_start, window_end = coords.parse_region_1based(args.region)
    except ValueError as error:
        raise UsageError(f"{error}（--region 要写成 chrom:start-end）") from None
    settings = insertion_plot.InsertionSettings(
        bam=args.bam,
        chrom=chrom,
        window_start=window_start,
        window_end=window_end,
        breakpoints=tuple(insertion.breakpoint_from_1based(value) for value in args.breakpoint),
        out_dir=args.out_dir,
        name=args.name,
        dpi=args.dpi,
        show_consensus=not args.no_consensus_text,
        show_formula=not args.no_formula,
    )
    result = insertion_plot.run(settings)
    print(f"图: {result['png']}")
    print(f"证据表: {result['clips_tsv']}")
    print(f"一致性序列表: {result['consensus_tsv']}")
    print(f"断点: {result['junctions']}  clip 证据: {result['n_reads']} 条 read")
    return 0


def _add_coverage(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "coverage",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="BAM 的区间覆盖广度与深度（PNG + 两张表）",
        description=(
            "按区间统计覆盖度：逐碱基深度 -> 分箱 -> 上图广度 / 下图深度（可选 BED 基因轨道）。"
            "设计上是给短读长（二代 / PE）数据看的：覆盖度均匀性、插入片段两侧的断层，"
            "长读长数据也能跑但通常用 `plot ribbon` 看结构。出图只读 TSV：先写 region/bins "
            "两张表，再由表渲染。"
        ),
        epilog=(
            "指标:\n"
            "  只有 M/=/X 给深度加 1；D/N 占参考位置但没有碱基，不算覆盖；S/I/H 不占参考位置。\n"
            "  depth >= --min-depth 的位置才算“被覆盖”；平均深度的分母是区间长度。\n"
            "\n例子:\n"
            "  alnkit plot coverage --bam reads.bam --region insertA --out-dir out\n"
            "  alnkit plot coverage --bam reads.bam --region chr1:1-20000 --bins 400 \\\n"
            "      --bed genes.bed --out-dir out --name region\n"
            "  alnkit plot coverage --bam reads.bam --list\n"
            "\n产物:\n"
            "  <name>.png\n"
            "  + <name>.region.tsv（coverage_region）\n"
            "  + <name>.bins.tsv（coverage_bins）\n"
            "  另加 --bases-table 时再写逐碱基表（coverage_bases）。"
        ),
    )
    parser.add_argument(
        "--bam", required=True, metavar="PATH", help="BAM/CRAM（按区间取数需要索引）"
    )
    parser.add_argument(
        "--region",
        default="",
        metavar="CHROM[:START-END]",
        help="区间（1-based 闭区间）；只写 contig 名就是整条序列；除 --list 外必填",
    )
    parser.add_argument("--out-dir", dest="out_dir", default=".", help="输出目录，默认当前目录")
    parser.add_argument(
        "--name",
        default="",
        help="输出文件名前缀，默认 <bam 主名>.<contig>[_start-end]",
    )
    parser.add_argument(
        "--bins", type=int, default=None, metavar="INT", help="bin 个数；默认约 25 bp/bin"
    )
    parser.add_argument("--bed", default=None, metavar="BED", help="基因轨道 BED（给了才画）")
    parser.add_argument("--title", default=None, help="覆盖默认标题（英文，避免缺字体）")
    parser.add_argument("--width", type=float, default=10.0, help="图宽（英寸，默认 10）")
    parser.add_argument("--dpi", type=int, default=200, help="PNG 分辨率（默认 200）")
    parser.add_argument(
        "--bases-table", dest="bases_table", default=None, metavar="TSV", help="额外写逐碱基深度表"
    )
    parser.add_argument(
        "--list", action="store_true", help="只列出参考序列与比对上的 reads 数后退出"
    )
    filters = parser.add_argument_group("过滤")
    filters.add_argument(
        "--min-mapq", dest="min_mapq", type=int, default=0, help="最小 MAPQ（默认 0）"
    )
    filters.add_argument(
        "--min-len", dest="min_len", type=int, default=0, help="最小 query 长度（默认 0）"
    )
    filters.add_argument(
        "--min-depth",
        dest="min_depth",
        type=int,
        default=1,
        help="低于该深度的位置不算覆盖（默认 1）",
    )
    filters.add_argument(
        "--max-depth",
        dest="max_depth",
        type=int,
        default=1_000_000,
        help="每个位置的深度上限，0 = 不限（默认 1000000）",
    )
    filters.add_argument(
        "--ff",
        dest="fail_flags",
        default="UNMAP,SECONDARY,QCFAIL,DUP",
        help="带这些 flag 的 read 跳过（默认 UNMAP,SECONDARY,QCFAIL,DUP；空串=不过滤）",
    )
    filters.add_argument(
        "--rf",
        dest="require_flags",
        default="",
        help="缺这些 flag 的 read 跳过（默认空）",
    )
    parser.set_defaults(handler=run_coverage)


def run_coverage(args: argparse.Namespace) -> int:
    coverage_plot = _require_plot_module("coverage")
    settings = coverage_plot.CoverageSettings(
        bam=args.bam,
        region=args.region,
        out_dir=args.out_dir,
        name=args.name,
        bins=args.bins,
        bed=args.bed,
        title=args.title,
        dpi=args.dpi,
        width=args.width,
        fail_flags=args.fail_flags,
        require_flags=args.require_flags,
        min_mapq=args.min_mapq,
        min_len=args.min_len,
        min_depth=args.min_depth,
        max_depth=args.max_depth,
        bases_table=args.bases_table,
        list_only=args.list,
    )
    result = coverage_plot.run(settings)
    if "references" in result:
        print("reference\tlength\tmapped")
        for contig, length, mapped in result["references"]:
            print(f"{contig}\t{length}\t{mapped}")
        return 0
    print(f"图: {result['png']}")
    print(f"区域表: {result['region_tsv']}")
    print(f"分箱表: {result['bins_tsv']}")
    print(
        f"{result['region']}  {result['reads']} reads  "
        f"{coverage_plot.format_decimal(result['percent_covered'])}% covered  "
        f"mean {coverage_plot.format_decimal(result['mean_depth'])}x  "
        f"mapQ {coverage_plot.format_decimal(result['mean_mapq'])}  "
        f"baseQ {coverage_plot.format_decimal(result['mean_baseq'])}"
    )
    return 0
