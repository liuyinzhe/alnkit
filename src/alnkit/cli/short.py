"""短读长子命令：``alnkit short scan``（BAM -> 断点证据三张表）。

参数名、默认值、帮助文本与旧脚本 ``parse_bam.py`` 一一对应。

出图命令不在这一支下面：插入位点示意图与覆盖度图都在 ``alnkit plot``
（``plot insertion`` / ``plot coverage``），与其他静图命令放在一起。
"""

from __future__ import annotations

import argparse
from pathlib import Path

from ..shortread import scan
from ._common import GLOBAL_OPTIONS


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    short = subparsers.add_parser(
        "short",
        help="短读长链路（二代插入断点证据）",
        description="从 BAM 里找插入断点的证据（SR / PESR / DUP）。",
    )
    short_subparsers = short.add_subparsers(dest="short_command", metavar="<subcommand>")
    short_subparsers.required = True
    _add_scan(short_subparsers)


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
