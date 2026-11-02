"""长读长子命令：``alnkit long <classify|graph|breakpoints|assemble|locate>``。

目前只有 ``classify``（A3）从旧流程迁入；其余随 P3 补齐。
"""

from __future__ import annotations

import argparse

from ..longread import assemble, breakpoints, classify, graph, locate
from ._common import GLOBAL_OPTIONS


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    long = subparsers.add_parser(
        "long",
        help="长读长链路（三代：断点 / 类 / 组装 / 落点）",
        description="PacBio/ONT 长读长上检测外源插入的结构：先分类，再建断点图与类。",
    )
    long_subparsers = long.add_subparsers(dest="long_command", metavar="<subcommand>")
    long_subparsers.required = True
    _add_graph(long_subparsers)
    _add_classify(long_subparsers)
    _add_breakpoints(long_subparsers)
    _add_assemble(long_subparsers)
    _add_locate(long_subparsers)


def _add_classify(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "classify",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="目标 reads 分类 + 接头事件/位点/区域聚类",
        description=(
            "先按外源证据筛出目标 read，再逐段判定外源/基因组/可疑错配，"
            "最后把接头事件聚成位点与区域，并导出相关 reads 子集。"
        ),
        epilog=(
            "上游 BAM 的生成要求（断点证据来自 CIGAR 的 S，所以这两条是硬要求）：\n"
            "  minimap2 -ax map-hifi -Y -N 50 -t 8 merged.fa reads.fq \\\n"
            "      | samtools sort -o reads.sorted.bam -  &&  samtools index reads.sorted.bam\n"
            "  -Y     补充比对用 soft clip（否则被 hard clip，取不到断点/归属）\n"
            "  -N 50  保留足够的次要比对（共有 read 的归属判定要用）\n"
            "\n例子：\n  alnkit long classify --bam reads_to_merged.sorted.bam \\\n"
            "          --reads input.fastq --merged-fai merged.fa.fai \\\n"
            "          --exog insertA --exog insertB --outdir 03_classify"
        ),
    )
    parser.add_argument("--bam", required=True, help="reads -> 合并参考 的 BAM")
    parser.add_argument("--reads", required=True, help="原始 fastq")
    parser.add_argument(
        "--merged-fai", dest="merged_fai", required=True, help="合并参考 .fai（取序列长度）"
    )
    parser.add_argument(
        "--exog", action="append", default=[], required=True, help="外源序列名，可重复"
    )
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument(
        "--min-mapq", dest="min_mapq", type=int, default=20, help="比对质量下限（默认 20）"
    )
    parser.add_argument(
        "--min-exog-block-bp",
        dest="min_exog_block_bp",
        type=int,
        default=100,
        help="外源比对段最小长度（默认 100 bp）",
    )
    parser.add_argument(
        "--contig-edge-bp",
        dest="contig_edge_bp",
        type=int,
        default=100,
        help="contig 两端多近算边缘（默认 100 bp）",
    )
    parser.add_argument(
        "--short-genome-block-bp",
        dest="short_genome_block_bp",
        type=int,
        default=500,
        help="边缘处基因组段多短算可疑错配（默认 500 bp）",
    )
    parser.add_argument(
        "--junction-max-gap-bp",
        dest="junction_max_gap_bp",
        type=int,
        default=10_000,
        help="基因组段与外源段在 read 轴上多近算相邻（默认 10 kb）",
    )
    parser.add_argument(
        "--site-gap-bp",
        dest="site_gap_bp",
        type=int,
        default=20,
        help="同一位点的边界允许错位（默认 20 bp）",
    )
    parser.add_argument(
        "--region-gap-bp",
        dest="region_gap_bp",
        type=int,
        default=200_000,
        help="同染色体位点归为同一区域的最大间距（默认 200 kb）",
    )
    parser.set_defaults(handler=run_classify)


def run_classify(args: argparse.Namespace) -> int:
    settings = classify.ClassifySettings(
        bam=args.bam,
        reads=args.reads,
        merged_fai=args.merged_fai,
        exog=tuple(args.exog),
        outdir=args.outdir,
        min_mapq=args.min_mapq,
        min_exog_block_bp=args.min_exog_block_bp,
        contig_edge_bp=args.contig_edge_bp,
        short_genome_block_bp=args.short_genome_block_bp,
        junction_max_gap_bp=args.junction_max_gap_bp,
        site_gap_bp=args.site_gap_bp,
        region_gap_bp=args.region_gap_bp,
    )
    result = classify.run(settings)
    print(
        f"完成：总 reads {result['total_reads']}，目标 {result['target_reads']}，"
        f"位点 {result['sites']}，区域 {result['regions']}；输出 {result['outdir']}"
    )
    return 0


def _add_breakpoints(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "breakpoints",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="断点图方案 A：紧邻片段断点 + 共享 reads 合并成类",
        description=(
            "在 read 轴上按“紧邻片段”定义断点（基因组–外源 = vector S##；"
            "基因组–基因组 = genome G##a/G##b），两个断点共享 reads >= --min-shared 就合并为一类。"
        ),
        epilog=(
            "例子：\n  alnkit long breakpoints \\\n"
            "          --segments 03_classify/alignment_segments.tsv \\\n"
            "          --classification 03_classify/read_classification.tsv \\\n"
            "          --overlaps 04_overlap/overlap_edges.tsv \\\n"
            "          --reads 03_classify/related_reads.fastq \\\n"
            "          --merged-fai merged.fa.fai --outdir 05_breakpoints"
        ),
    )
    parser.add_argument("--segments", required=True, help="alignment_segments.tsv")
    parser.add_argument("--classification", required=True, help="read_classification.tsv")
    parser.add_argument("--overlaps", required=True, help="overlap_edges.tsv")
    parser.add_argument("--reads", required=True, help="related_reads.fastq（输出每类子集用）")
    parser.add_argument(
        "--merged-fai", dest="merged_fai", required=True, help="合并参考 .fai（取序列长度）"
    )
    parser.add_argument("--outdir", required=True, help="输出目录（建议 05_breakpoints）")
    parser.add_argument(
        "--site-gap",
        dest="site_gap",
        type=int,
        default=20,
        help="同一 site 的边界允许错位（默认 20 bp）",
    )
    parser.add_argument(
        "--min-shared",
        dest="min_shared",
        type=int,
        default=2,
        help="合并断点所需最少共享 reads（默认 2）",
    )
    parser.set_defaults(handler=run_breakpoints)


def run_breakpoints(args: argparse.Namespace) -> int:
    settings = breakpoints.BreakpointSettings(
        segments=args.segments,
        classification=args.classification,
        overlaps=args.overlaps,
        reads=args.reads,
        merged_fai=args.merged_fai,
        outdir=args.outdir,
        site_gap=args.site_gap,
        min_shared=args.min_shared,
    )
    result = breakpoints.run(settings)
    print(
        f"完成：断点 {result['sites']} 个，类 {result['classes']} 个，"
        f"segment_ends {result['segment_ends']} 行；输出 {result['outdir']}"
    )
    return 0


def _add_assemble(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "assemble",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="按 overlap 关系贪心延伸组装每个类",
        description=(
            "每个类先拿一条最长的接头 read 作种子，再从 overlap 边里反复挑一条能往"
            "左/右延伸的 read，把非重叠悬垂序列接到 contig 两端。"
        ),
        epilog=(
            "例子：\n  alnkit long assemble \\\n"
            "          --classes-dir 05_breakpoints/classes \\\n"
            "          --paf 04_overlap/overlaps.paf --outdir 06_assembly"
        ),
    )
    parser.add_argument(
        "--classes-dir", dest="classes_dir", required=True, help="05_breakpoints/classes"
    )
    parser.add_argument("--paf", required=True, help="all-vs-all PAF")
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument(
        "--min-overlap",
        dest="min_overlap",
        type=int,
        default=500,
        help="延伸所需最小 overlap（默认 500 bp）",
    )
    parser.add_argument(
        "--min-identity",
        dest="min_identity",
        type=float,
        default=0.90,
        help="延伸所需最小一致率（默认 0.90）",
    )
    parser.add_argument(
        "--agreement-bp",
        dest="agreement_bp",
        type=int,
        default=200,
        help="多个邻居给出同一位置时允许的起点差（默认 200 bp）",
    )
    parser.set_defaults(handler=run_assemble)


def run_assemble(args: argparse.Namespace) -> int:
    settings = assemble.AssembleSettings(
        classes_dir=args.classes_dir,
        paf=args.paf,
        outdir=args.outdir,
        min_overlap=args.min_overlap,
        min_identity=args.min_identity,
        agreement_bp=args.agreement_bp,
    )
    result = assemble.run(settings)
    print(f"完成：{result['classes']} 个类，{result['contigs']} 条 contig；输出 {result['outdir']}")
    return 0


def _add_locate(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "locate",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="50 bp 断点锚点 + 外源查询 -> 丝带图 + 参考坐标还原",
        description=(
            "为每个类构造查询序列（断点旁 50 bp 锚点 + 全部外源），与类组装比对后逐条 contig 出图，"
            "并把锚点/外源在 contig 上的落点还原成参考基因组坐标。"
        ),
        epilog=(
            "例子：\n  alnkit long locate \\\n"
            "          --regions 05_breakpoints/breakpoint_classes.tsv \\\n"
            "          --sites 05_breakpoints/breakpoint_sites.tsv \\\n"
            "          --reference-fasta reference.fa \\\n"
            "          --exog exog/insertA.fa --exog exog/insertB.fa \\\n"
            "          --assembly-dir 06_assembly --outdir 07_ribbon_positions"
        ),
    )
    parser.add_argument("--regions", required=True, help="breakpoint_classes.tsv")
    parser.add_argument("--sites", required=True, help="breakpoint_sites.tsv")
    parser.add_argument(
        "--reference-fasta",
        dest="reference_fasta",
        required=True,
        help="参考基因组 fasta（需 .fai）",
    )
    parser.add_argument(
        "--exog", action="append", default=[], required=True, help="外源 fasta，可重复"
    )
    parser.add_argument("--assembly-dir", dest="assembly_dir", required=True, help="06_assembly")
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument("--threads", type=int, default=12, help="minimap2 线程数（默认 12）")
    parser.add_argument(
        "--flank-bp", dest="flank_bp", type=int, default=50, help="断点锚点长度（默认 50 bp）"
    )
    parser.set_defaults(handler=run_locate)


def run_locate(args: argparse.Namespace) -> int:
    settings = locate.LocateSettings(
        regions=args.regions,
        sites=args.sites,
        reference_fasta=args.reference_fasta,
        exog=tuple(args.exog),
        assembly_dir=args.assembly_dir,
        outdir=args.outdir,
        threads=args.threads,
        flank_bp=args.flank_bp,
    )
    result = locate.run(settings)
    print(
        f"完成：{result['classes']} 个类，{result['figures']} 张图，"
        f"位置行 {result['position_rows']} 行；输出 {result['outdir']}"
    )
    return 0


def _add_graph(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "graph",
        parents=[GLOBAL_OPTIONS],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        help="read-vs-read overlap 边与连通分量",
        description=(
            "读 all-vs-all PAF（minimap2 -x ava-pb -c），每对 read 只保留最好的一条 overlap，"
            "再求全局连通分量。方案 A 的类别由断点共享 reads 决定，这里只用于观察 overlap 结构。"
        ),
        epilog=(
            "例子：minimap2 -x ava-pb -c -t 8 related.fastq related.fastq > overlaps.paf\n"
            "     alnkit long graph --paf overlaps.paf \\\n"
            "          --classification 03_classify/read_classification.tsv --outdir 04_overlap"
        ),
    )
    parser.add_argument("--paf", required=True, help="all-vs-all PAF（minimap2 -x ava-pb -c）")
    parser.add_argument("--classification", required=True, help="read_classification.tsv")
    parser.add_argument("--outdir", required=True, help="输出目录")
    parser.add_argument(
        "--min-overlap",
        dest="min_overlap",
        type=int,
        default=300,
        help="边所需最小 overlap（默认 300 bp）",
    )
    parser.add_argument(
        "--min-identity",
        dest="min_identity",
        type=float,
        default=0.85,
        help="边所需最小一致率（默认 0.85）",
    )
    parser.set_defaults(handler=run_graph)


def run_graph(args: argparse.Namespace) -> int:
    settings = graph.GraphSettings(
        paf=args.paf,
        classification=args.classification,
        outdir=args.outdir,
        min_overlap=args.min_overlap,
        min_identity=args.min_identity,
    )
    result = graph.run(settings)
    print(f"完成：overlap 边 {result['edges']} 条，连通分量 {result['components']} 个")
    return 0
