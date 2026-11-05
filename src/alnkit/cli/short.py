"""短读长子命令：``alnkit short scan`` / ``short insertion-plot`` / ``short coverage``。

``short scan`` 的参数名、默认值、帮助文本与旧脚本 ``parse_bam.py`` 一一对应；
``insertion-plot`` 与 ``coverage`` 是两个绘图命令，参数取自并入的两个独立脚本
（``insertion_schematic.py`` / ``coverage_viz``），但**出图改成从 TSV 出**。
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..core import coords
from ..core.errors import UsageError
from ..shortread import insertion, scan
from ._common import GLOBAL_OPTIONS


def _require_plot_module(module_name: str):
    """绘图子命令才需要 matplotlib / numpy；缺了要给可行动提示而不是 traceback。

    基础安装（``pip install alnkit``）只带 pysam，所以这两个命令必须惰性导入，
    否则连 ``alnkit --version`` 都会在 import 阶段就崩。
    """
    import importlib

    try:
        return importlib.import_module(f"alnkit.plot.{module_name}")
    except ImportError as exc:  # pragma: no cover - 只在缺 matplotlib 时触发
        raise UsageError(
            f'绘图子命令需要 matplotlib 与 numpy（{exc}）。安装：pip install "alnkit[plot]"'
        ) from exc


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    short = subparsers.add_parser(
        "short",
        help="短读长链路（二代插入断点证据）",
        description="从 BAM 里找插入断点的证据（SR / PESR / DUP）。",
    )
    short_subparsers = short.add_subparsers(dest="short_command", metavar="<subcommand>")
    short_subparsers.required = True
    _add_scan(short_subparsers)
    _add_insertion_plot(short_subparsers)
    _add_coverage(short_subparsers)


def _add_scan(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "scan",
        parents=[GLOBAL_OPTIONS],
        help="BAM -> 断点证据三张表",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "一条跨越“参考基因组 <-> 外源插入序列”的嵌合比对可定位断点；"
            "本命令按断点聚合证据（SR / PESR / DUP）。"
        ),
        epilog=(
            "上游 BAM 的比对要求（详见原脚本 README）:\n"
            "  必须 -Y（补充比对用 soft clip，只从 CIGAR 的 S 取断点）；\n"
            "  禁止 -M（会把补充比对标成 secondary，证据变少）；\n"
            "  建议 -a -h 20；比对后不要用 samtools -F 2048 / -F 256 / -q 过滤。\n"
            "\n例子:\n"
            "  alnkit short scan --bam-dir ./bam --query-chr target_ref \\\n"
            "      --insert-chrs decoy_ref,homolog_ref,target_ref \\\n"
            "      --merge-window-bp 10 --out-dir out\n"
            "  alnkit short scan --self-test\n"
            "\n产物:\n"
            "  三张 legacy .xls（列名与旧脚本一致）+ <sample>.sites.canonical.tsv（机器接口）。\n"
            "  注意: reads_info 的 legacy 表头把 stand 与 mapping_quality 合并成一列（11 列），\n"
            "  数据行其实是 12 列——按表头名解析会错位，按位置或用 canonical 表。\n"
            "\n错配来源与 markdup:\n"
            "  --mismatch-source nm（默认）读 NM；md 用 MD 推导（历史口径）。都不自动回退，\n"
            "  标签缺失时按 --on-missing-tag 丢弃或报错；--legacy-profile 隐含 md。\n"
            "  markdup 的 FLAG 0x400 与错配统计无关，只决定 dup 两列；用 --check-markdup 核对\n"
            "  标记是否可信（按两端 5' 坐标分组，再用 MD 判断被标 dup 的 read 是否真同序列）。"
        ),
    )
    parser.add_argument(
        "--bam",
        action="append",
        type=Path,
        default=[],
        metavar="PATH",
        help="直接指定 BAM（可重复）。与 --bam-dir 二选一",
    )
    parser.add_argument(
        "--bam-dir",
        action="store",
        type=Path,
        default=Path.cwd(),
        metavar="DIR",
        help="扫描 BAM 的目录，默认当前目录",
    )
    parser.add_argument(
        "--bam-glob",
        action="store",
        default="*.sorted.bam",
        metavar="PATTERN",
        help="在 --bam-dir 下匹配 BAM 的通配符，默认 *.sorted.bam",
    )
    parser.add_argument("--no-recursive", action="store_true", help="不递归子目录（默认递归）")
    parser.add_argument(
        "--out-dir",
        action="store",
        type=Path,
        default=Path.cwd(),
        metavar="DIR",
        help="结果输出目录，默认当前目录",
    )
    parser.add_argument(
        "--sample-name",
        action="store",
        default=None,
        metavar="NAME",
        help="输出前缀（所有产物文件名的开头），默认取 BAM 文件名首个 . 之前的部分",
    )
    parser.add_argument(
        "--query-chr",
        action="store",
        default=None,
        metavar="NAME",
        help="本次要出结果的外源插入序列名（必填，除非 --self-test）",
    )
    parser.add_argument(
        "--insert-chrs",
        action="store",
        default=None,
        metavar="NAME1,NAME2,...",
        help="参考序列里全部外源/同源序列名，逗号分隔；默认只含 --query-chr",
    )
    parser.add_argument(
        "--insert-origin-policy",
        action="store",
        choices=("owned", "owned-strict", "lenient", "consensus", "legacy"),
        default="owned",
        help="共有 read 的区分策略：owned（默认）/owned-strict/lenient/consensus/legacy",
    )
    parser.add_argument(
        "--merge-window-bp",
        action="store",
        type=int,
        default=0,
        metavar="INT",
        help="位点窗口合并：相距 <= 该值的断点合并为一个位点，0=不合并（默认）；推荐 10",
    )
    parser.add_argument(
        "--min-mapq",
        action="store",
        type=int,
        default=20,
        metavar="INT",
        help="比对质量 <= 该值则丢弃，默认 20",
    )
    parser.add_argument(
        "--max-mismatch",
        action="store",
        type=int,
        default=3,
        metavar="INT",
        help="错配数 > 该值则丢弃，默认 3",
    )
    parser.add_argument(
        "--max-softclip-blocks",
        action="store",
        type=int,
        default=1,
        metavar="INT",
        help="soft clip 段数 > 该值则丢弃，默认 1",
    )
    parser.add_argument(
        "--min-clip-length",
        action="store",
        type=int,
        default=0,
        metavar="INT",
        help="作为断点证据时要求的最小 soft clip 长度，默认 0（不限制）",
    )
    parser.add_argument(
        "--keep-secondary", action="store_true", help="保留 secondary 比对（默认丢弃）"
    )
    parser.add_argument(
        "--mismatch-source",
        action="store",
        choices=("nm", "md"),
        default="nm",
        help=(
            "错配数来源：nm=NM 标签（默认）；md=由 MD 推导（历史行为）。"
            "两者都不自动回退：选 nm 却没有 NM 时按 --on-missing-tag 丢弃或报错"
        ),
    )
    parser.add_argument(
        "--on-missing-tag",
        action="store",
        choices=("skip", "error"),
        default="skip",
        help="缺少 NM/MD 标签时：skip=丢弃该比对（默认）；error=报错",
    )
    parser.add_argument(
        "--no-split-match",
        action="store_true",
        help="关闭「目标序列上要求存在等长 split 比对」的校验（默认开启）",
    )
    parser.add_argument(
        "--no-candidate-gate",
        action="store_true",
        help="关闭候选 read 闸门（等价 2.4.0 全量逐条处理）；legacy 策略下自动关闭",
    )
    parser.add_argument(
        "--no-duplicate-filter",
        action="store_true",
        help="保留被 markdup 标记为重复的 read 作为证据（默认剔除，只写 dup 列）",
    )
    parser.add_argument(
        "--min-sr-support",
        dest="min_sr_support",
        type=int,
        default=1,
        metavar="INT",
        help=(
            "位点至少要有多少条支持 read（SR+PESR）才输出，默认 1（不过滤，保持历史行为）；"
            "设 2–3 可去掉单 read 噪声位点"
        ),
    )
    parser.add_argument(
        "--interference-scope",
        dest="interference_scope",
        choices=("all", "other"),
        default="all",
        help=(
            "同源干扰的判据范围：all=支持 read 命中任何插入序列都算干扰（默认，历史口径，"
            "多外源面板下偏灵敏）；other=只算 query 之外的插入（多外源面板更精确）。"
            "注意：两者都只比较“插入之间”，无法识别“外源与基因组同源”造成的干扰"
        ),
    )
    parser.add_argument(
        "--reference-fasta",
        dest="reference_fasta",
        help=(
            "含插入序列与基因组的参考 FASTA；给了它就把插入序列回比基因组，"
            "自动找出“外源同源区”并标注位点（需要 minimap2）"
        ),
    )
    parser.add_argument(
        "--homology-bed",
        dest="homology_bed",
        action="append",
        default=[],
        metavar="BED",
        help="已知的内源同源基因/区域（BED，0-based 半开，可给多次）；与 --reference-fasta 可并用",
    )
    parser.add_argument(
        "--homology-window-bp",
        dest="homology_window_bp",
        type=int,
        default=50000,
        metavar="INT",
        help=(
            "局部回比时位点两侧各取多长（默认 50000 bp）；cDNA 样插入的外显子分散在几十 kb 里，"
            "窗口太小只能比到单个外显子"
        ),
    )
    parser.add_argument(
        "--homology-min-identity",
        dest="homology_min_identity",
        type=float,
        default=0.90,
        metavar="FLOAT",
        help="回比时认定的最小一致率，默认 0.90",
    )
    parser.add_argument(
        "--homology-min-length",
        dest="homology_min_length",
        type=int,
        default=200,
        metavar="INT",
        help="回比时认定的最小比对长度（bp），默认 200（过滤偶然的短命中）",
    )
    parser.add_argument(
        "--check-markdup",
        action="store_true",
        help=(
            "核对 markdup 标记是否可信：把 read 按「比对坐标 + 错配数」分组，报告"
            "全未标记组（markdup 没做/不完整）、标记不一致组、以及单条被标 dup 的 read"
        ),
    )
    parser.add_argument(
        "--require-dup-flags",
        action="store_true",
        help="BAM 里没有任何 FLAG 0x400（未做 markdup）时报错退出（退出码 2）",
    )
    parser.add_argument(
        "--bam-threads",
        action="store",
        type=int,
        default=1,
        metavar="INT",
        help="BAM 解压线程数，默认 1；大文件可调到 2-4",
    )
    parser.add_argument(
        "--legacy-coordinates",
        action="store_true",
        help="使用历史的 M 块选取规则（取最后一个 M 块）",
    )
    parser.add_argument(
        "--legacy-profile",
        action="store_true",
        help="一键复现历史结果（= legacy 坐标 + legacy 归属 + 由 MD 推导错配）",
    )
    parser.add_argument("--self-test", action="store_true", help="只跑内置自检（不读 BAM）")
    parser.set_defaults(handler=run_scan)


def run_scan(args: argparse.Namespace) -> int:
    if args.self_test:
        scan.run_self_test()
        return 0
    settings = scan.ScanSettings(
        bam=tuple(args.bam),
        bam_dir=args.bam_dir,
        bam_glob=args.bam_glob,
        no_recursive=args.no_recursive,
        out_dir=args.out_dir,
        sample_name=args.sample_name,
        query_chr=args.query_chr,
        insert_chrs=args.insert_chrs,
        insert_origin_policy=args.insert_origin_policy,
        merge_window_bp=args.merge_window_bp,
        min_mapq=args.min_mapq,
        max_mismatch=args.max_mismatch,
        max_softclip_blocks=args.max_softclip_blocks,
        min_clip_length=args.min_clip_length,
        keep_secondary=args.keep_secondary,
        mismatch_source=args.mismatch_source,
        on_missing_tag=args.on_missing_tag,
        no_split_match=args.no_split_match,
        no_candidate_gate=args.no_candidate_gate,
        no_duplicate_filter=args.no_duplicate_filter,
        require_dup_flags=args.require_dup_flags,
        check_markdup=args.check_markdup,
        reference_fasta=Path(args.reference_fasta) if args.reference_fasta else None,
        homology_bed=tuple(Path(p) for p in args.homology_bed),
        homology_min_identity=args.homology_min_identity,
        homology_window_bp=args.homology_window_bp,
        homology_min_length=args.homology_min_length,
        minimum_site_support=args.min_sr_support,
        interference_scope=args.interference_scope,
        bam_threads=args.bam_threads,
        legacy_coordinates=args.legacy_coordinates,
        legacy_profile=args.legacy_profile,
        log_level=args.log_level,
        quiet=args.quiet,
    )
    result = scan.run(settings)
    print(
        f"完成：{len(result['bams'])} 个 BAM，"
        f"query={result['query_chr']}，输出目录 {result['out_dir']}"
    )
    return 0


def _add_insertion_plot(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "insertion-plot",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="外源插入位点示意图（单侧 / 同断点双侧 / 多断点）",
        description=(
            "窗口内按 soft clip 证据画插入位点示意图；脚本按“给了几个断点 + 每个断点有几侧 "
            "clip”自动选画法。出图只读 TSV：先写 clips/consensus 两张表，再由表渲染。"
        ),
        epilog=(
            "画法:\n"
            "  1 个断点、只有一侧 clip -> 单侧图（未知条 + 彩色块贴在断点一侧）\n"
            "  1 个断点、两侧都有 clip -> 同断点双侧图（左块 + 合并未知区 + 右块，连线收到断点）\n"
            "  2-3 个断点              -> 多断点图（每断点各自 [未知][块]，块间强制留间隔）\n"
            "\n配色: 左 clip 蓝 / 右 clip 紫 / 未知区灰斜线 / 断点红。\n"
            "\n例子:\n"
            "  alnkit short insertion-plot --bam reads.bam --region chr1:99201-100800 \\\n"
            "      --breakpoint 100000 --out-dir out --name demo\n"
            "  alnkit short insertion-plot --bam reads.bam --region chr1:699101-700900 \\\n"
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
    parser.add_argument(
        "--name", default="insertion_plot", help="输出文件名前缀，默认 insertion_plot"
    )
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
        help="区间覆盖广度与深度（PNG + 两张表）",
        description=(
            "按区间统计覆盖度：逐碱基深度 -> 分箱 -> 上图广度 / 下图深度（可选 BED 基因轨道）。"
            "出图只读 TSV：先写 region/bins 两张表，再由表渲染。"
        ),
        epilog=(
            "指标:\n"
            "  只有 M/=/X 给深度加 1；D/N 占参考位置但没有碱基，不算覆盖；S/I/H 不占参考位置。\n"
            "  depth >= --min-depth 的位置才算“被覆盖”；平均深度的分母是区间长度。\n"
            "\n例子:\n"
            "  alnkit short coverage --bam reads.bam --region insertA --out-dir out\n"
            "  alnkit short coverage --bam reads.bam --region chr1:1-20000 --bins 400 \\\n"
            "      --bed genes.bed --out-dir out --name region\n"
            "  alnkit short coverage --bam reads.bam --list\n"
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
