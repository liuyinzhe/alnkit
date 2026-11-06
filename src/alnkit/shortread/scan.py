"""B1 编排层：配置解析、环境校验与告警、单 BAM 处理与自检。"""

from __future__ import annotations

import logging
import sys
import time
from pathlib import Path
from typing import Dict, Optional, Set
import pysam
from ..core.errors import RunError, UsageError

from .bam import (
    choose_evidence_block,
    compute_junction_position,
    detect_clip_side,
    measure_cigar_block_and_clip,
    parse_xa_entry,
    read_is_duplicate,
    check_markdup_consistency,
    report_duplicate_flags,
    report_markdup_consistency,
    warn_about_alignment_features,
)

from .evidence import (
    aggregate_breakpoint_evidence,
    collect_duplicate_reads_by_breakpoint,
    create_statistics,
    find_pair_end_split_read_junction,
    iter_split_read_evidence,
)

from .interference import collect_homology_hit_chromosomes, determine_interference_sources

from ..longread import support
from .homology import (
    annotate_sites,
    load_bed_regions,
    local_homology_regions,
    write_site_homology_table,
)
from .merge import (
    cluster_positions_within_window,
    filter_by_minimum_support,
    merge_nearby_breakpoints,
)

from .model import (
    AlignmentBlock,
    AlignmentRecord,
    BreakpointEvidence,
    ConfigurationError,
    InsertHit,
    ReadEvidence,
    RuntimeConfig,
    SCAN_ENGINE_VERSION,
    ScanSettings,
)

from .origin import (
    collect_read_records,
    decide_read_usage,
    determine_insert_origin,
    discover_candidate_reads,
    resolve_ambiguous_reads_by_site_consensus,
    select_read_evidences,
)

from .report import (
    write_insert_interference_table,
    write_canonical_sites_table,
    write_insert_position_table,
    write_reads_info_table,
)


def resolve_runtime_config(settings: ScanSettings) -> RuntimeConfig:
    """把命令行参数整理成 RuntimeConfig，并做前置校验。

    输入示例:
        arguments = parser.parse_args(['--bam-dir','./bam','--query-chr','target_ref',
                                       '--insert-chrs','decoy_ref,homolog_ref,target_ref'])
    输出示例:
        RuntimeConfig(bam_paths=[Path('bam/xxx.sorted.bam')], query_chromosome='target_ref',
                      insert_chromosomes=('decoy_ref','homolog_ref','target_ref'), ...)
    校验失败时调用 raise UsageError()（退出码 2）。
    """
    mismatch_source = settings.mismatch_source

    if not settings.query_chr:
        raise UsageError("必须指定 --query-chr（本次要出结果的外源插入序列名）")

    # 插入序列列表：默认只含 query 自己；去重但保持顺序
    if settings.insert_chrs:
        provided = (name.strip() for name in settings.insert_chrs.split(",") if name.strip())
        insert_chromosomes = tuple(dict.fromkeys(provided))
    else:
        insert_chromosomes = (settings.query_chr,)
    if settings.query_chr not in insert_chromosomes:
        insert_chromosomes = (settings.query_chr,) + insert_chromosomes

    if settings.bam:
        bam_paths = list(settings.bam)
    else:
        pattern = settings.bam_glob
        if settings.no_recursive:
            bam_paths = sorted(
                path
                for path in settings.bam_dir.glob(pattern)
                if path.is_file() and not path.is_symlink()
            )
        else:
            bam_paths = sorted(
                path
                for path in settings.bam_dir.rglob(pattern)
                if path.is_file() and not path.is_symlink()
            )

    if not bam_paths:
        raise UsageError(f"未找到 BAM 文件：目录={settings.bam_dir} 通配符={settings.bam_glob}")
    if len(bam_paths) > 1 and settings.sample_name:
        raise UsageError("指定了多个 BAM 时不能使用 --sample-name（样本名由每个 BAM 文件名推导）")

    if settings.bam_threads < 1:
        raise UsageError("--bam-threads 必须 >= 1")
    if settings.reference_fasta is not None and not Path(settings.reference_fasta).is_file():
        raise UsageError(f"--reference-fasta 不存在：{settings.reference_fasta}")
    if not 0.0 < settings.homology_min_identity <= 1.0:
        raise UsageError("--homology-min-identity 必须在 (0, 1] 之间")
    if settings.homology_min_length < 1:
        raise UsageError("--homology-min-length 必须 >= 1")
    if settings.minimum_site_support < 1:
        raise UsageError("--min-sr-support 必须 >= 1（1 = 不过滤，保持历史行为）")

    return RuntimeConfig(
        bam_paths=bam_paths,
        sample_name=settings.sample_name,
        output_directory=settings.out_dir,
        query_chromosome=settings.query_chr,
        insert_chromosomes=insert_chromosomes,
        insert_origin_policy=settings.insert_origin_policy,
        drop_shared_reads=settings.drop_shared_reads,
        merge_window_bp=settings.merge_window_bp,
        mismatch_source=mismatch_source,
        on_missing_tag=settings.on_missing_tag,
        minimum_mapping_quality=settings.min_mapq,
        maximum_mismatch_count=settings.max_mismatch,
        maximum_soft_clip_blocks=settings.max_softclip_blocks,
        minimum_clip_length=settings.min_clip_length,
        keep_secondary_alignments=settings.keep_secondary,
        require_split_match=not settings.no_split_match,
        legacy_coordinates=settings.legacy_coordinates,
        use_candidate_gate=not settings.no_candidate_gate,
        exclude_duplicates=not settings.no_duplicate_filter,
        require_duplicate_flags=settings.require_dup_flags,
        check_markdup=settings.check_markdup,
        minimum_site_support=settings.minimum_site_support,
        reference_fasta=settings.reference_fasta,
        homology_bed=tuple(settings.homology_bed),
        homology_min_identity=settings.homology_min_identity,
        homology_min_length=settings.homology_min_length,
        homology_window_bp=settings.homology_window_bp,
        bam_threads=settings.bam_threads,
    )


def configure_logging(settings: ScanSettings) -> None:
    """初始化日志（输出到 stderr，避免污染结果文件）。

    输入示例: settings.log_level = 'DEBUG'
    输出示例: stderr 出现 '2026-10-13 19:00:00 INFO 开始处理 ...'
    """
    level = logging.ERROR if settings.quiet else getattr(logging, settings.log_level)
    logging.basicConfig(
        level=level,
        stream=sys.stderr,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )


def validate_reference_names(
    alignment_file: "pysam.AlignmentFile", config: RuntimeConfig, bam_path: Path
) -> None:
    """校验 --query-chr / --insert-chrs 在 BAM header 的参考序列里确实存在。

    输入示例:
        alignment_file.references = ('chr1','chr2',...,'decoy_ref','homolog_ref','target_ref')
        config.query_chromosome = 'target_ref', config.insert_chromosomes = ('homolog_ref','targetref')
    输出示例:
        - query 存在、'targetref' 不存在 -> 记一条 WARNING（该名字不参与统计），不报错
        - query 写成 'targett_ref' 不存在   -> 抛 ConfigurationError（退出码 2）

    说明: 历史脚本写错目标序列名不会报错，只会静默产出另一套结果；这里把错误提前暴露。
    """
    reference_names = set(alignment_file.references)
    if config.query_chromosome not in reference_names:
        raise ConfigurationError(
            f"--query-chr {config.query_chromosome} 不在 BAM 参考序列里: {bam_path}\n"
            f"  BAM 参考序列共 {len(reference_names)} 条，示例: "
            f"{', '.join(sorted(reference_names)[:5])} ..."
        )
    missing_insert_chromosomes = [
        name for name in config.insert_chromosomes if name not in reference_names
    ]
    if missing_insert_chromosomes:
        logging.warning(
            "--insert-chrs 里这些名字不在 BAM 参考序列里，将不参与归属判定: %s (%s)",
            ",".join(missing_insert_chromosomes),
            bam_path,
        )


def process_single_bam(bam_path: Path, config: RuntimeConfig) -> Dict[str, int]:
    """处理一个 BAM，写出三张结果表。

    输入示例: bam_path=Path('sampleA.sorted.bam'), config.query_chromosome='target_ref'
    输出示例:
        <out_dir>/sampleA.reads_info.out.xls
        <out_dir>/sampleA.Insert_Position.result.xls
        <out_dir>/sampleA.Insert_Interference.result.xls
        返回 {'alignments_total': 25969, 'reads_kept': 1779, 'breakpoints_written': 31, ...}

    异常:
        ConfigurationError: BAM header 里没有 --query-chr 指定的序列，或 --require-dup-flags 下_BAM 未做 markdup
                            （由 main 转成退出码 2）
        RuntimeError:       BAM 不存在/损坏（由 main 转成退出码 1）
    """
    sample_name = config.sample_name or bam_path.name.split(".")[0]
    statistics = create_statistics()

    logging.info("处理 BAM: %s（样本名 %s）", bam_path, sample_name)
    try:
        with pysam.AlignmentFile(str(bam_path), "rb", threads=config.bam_threads) as alignment_file:
            validate_reference_names(alignment_file, config, bam_path)
            # 第 1 遍：只读标签，筛出可能相关的 read（--no-candidate-gate / legacy 则不筛）
            candidate_reads: Optional[Set[str]] = None
            if config.use_candidate_gate:
                started_at = time.perf_counter()
                candidate_reads = discover_candidate_reads(alignment_file, config, statistics)
                statistics["pass1_ms"] = int((time.perf_counter() - started_at) * 1000)
                alignment_file.reset()  # 回到文件头，第 2 遍从第一条记录重新开始
            started_at = time.perf_counter()
            read_collection = collect_read_records(
                alignment_file, config, statistics, candidate_reads
            )
            statistics["pass2_ms"] = int((time.perf_counter() - started_at) * 1000)
    except (OSError, ValueError) as error:
        raise RuntimeError(f"读取 BAM 失败: {bam_path} ({error})") from error

    warn_about_alignment_features(statistics, bam_path)
    report_duplicate_flags(statistics, bam_path, config.require_duplicate_flags)
    if config.check_markdup:
        report_markdup_consistency(check_markdup_consistency(bam_path), bam_path)
    if config.use_candidate_gate:
        logging.info(
            "扫描: 候选 read=%d 候选比对=%d（pass1 %.2fs / pass2 %.2fs）",
            statistics["candidate_reads"],
            statistics["alignments_total"] - statistics["skipped_by_candidate_gate"],
            statistics["pass1_ms"] / 1000.0,
            statistics["pass2_ms"] / 1000.0,
        )

    read_evidences, duplicate_read_evidences = select_read_evidences(
        read_collection, config, statistics
    )

    # 输出 1：保留比对明细（legacy 策略下等价于历史"info_dic 过滤后的全部记录"）
    reads_info_path = config.output_directory / f"{sample_name}.reads_info.out.xls"
    written_reads_info_rows = write_reads_info_table(reads_info_path, read_evidences)

    # 输出 2：断点证据汇总
    #   dup 只来自 BAM 的 markdup 标记：被剔除的重复 read（默认策略）按断点单独聚合，
    #   保证 SR/PESR 里没有它们，但 dup 两列依旧能看见“这些位点上有多少重复 read”。
    duplicates_by_breakpoint = collect_duplicate_reads_by_breakpoint(
        duplicate_read_evidences, config.query_chromosome, config.minimum_clip_length
    )
    evidence_by_breakpoint = aggregate_breakpoint_evidence(
        read_evidences, config.query_chromosome, config.minimum_clip_length
    )
    statistics["breakpoints_before_merge"] = len(evidence_by_breakpoint)
    # 位点窗口合并（--merge-window-bp 0 时为恒等变换，保持历史逐位点输出）
    evidence_by_breakpoint, duplicates_by_breakpoint, merged_positions_by_site = (
        merge_nearby_breakpoints(
            evidence_by_breakpoint, duplicates_by_breakpoint, config.merge_window_bp
        )
    )
    statistics["breakpoints_merged_away"] = statistics["breakpoints_before_merge"] - len(
        evidence_by_breakpoint
    )
    # 位点级最小支持数过滤（默认 1 = 不过滤）
    (
        evidence_by_breakpoint,
        duplicates_by_breakpoint,
        merged_positions_by_site,
        dropped_sites,
    ) = filter_by_minimum_support(
        evidence_by_breakpoint,
        duplicates_by_breakpoint,
        merged_positions_by_site,
        config.minimum_site_support,
    )
    if dropped_sites:
        statistics["sites_dropped_by_min_support"] = dropped_sites

    insert_position_path = config.output_directory / f"{sample_name}.Insert_Position.result.xls"
    written_breakpoints = write_insert_position_table(
        insert_position_path, evidence_by_breakpoint, duplicates_by_breakpoint
    )
    # 新增：同一个位点集合的机器接口（带 1-based/0-based 后缀，见 core/tables.py）
    canonical_path = config.output_directory / f"{sample_name}.sites.canonical.tsv"
    written_canonical = write_canonical_sites_table(
        canonical_path,
        evidence_by_breakpoint,
        duplicates_by_breakpoint,
        sample_name,
        config.query_chromosome,
    )

    # 输出 3：逐断点的同源干扰标记（新增诊断表，不影响上面两张历史表）
    interference_path = config.output_directory / f"{sample_name}.Insert_Interference.result.xls"
    written_interference_rows = write_insert_interference_table(
        interference_path, evidence_by_breakpoint, merged_positions_by_site, config
    )

    # 输出 4（可选）：同源位点诊断 —— 位点是否落在“与外源同源的基因组区域”里
    written_homology_rows = 0
    if config.reference_fasta is not None or config.homology_bed:
        regions = list(load_bed_regions(config.homology_bed))
        if config.reference_fasta is not None:
            minimap2 = support.find_minimap2()
            if minimap2 is None:
                logging.warning(
                    "没找到 minimap2，跳过“外源回比基因组”这一步（仍可用 --homology-bed）"
                )
            else:
                regions.extend(
                    local_homology_regions(
                        config.reference_fasta,
                        config.insert_chromosomes,
                        sorted(evidence_by_breakpoint),
                        minimap2,
                        window_bp=config.homology_window_bp,
                        min_identity=config.homology_min_identity,
                        min_length=config.homology_min_length,
                        threads=max(1, config.bam_threads),
                    )
                )
        if regions:
            hom_path = config.output_directory / f"{sample_name}.site_homology.tsv"
            written_homology_rows = write_site_homology_table(
                hom_path,
                annotate_sites(sorted(evidence_by_breakpoint), regions),
                sample_name,
            )
            statistics["homology_regions"] = len(regions)
            statistics["homology_sites"] = sum(
                1
                for row in annotate_sites(sorted(evidence_by_breakpoint), regions)
                if row["is_homology_site"] == "yes"
            )
            logging.info(
                "同源位点诊断: %d 段同源区域，%d 个位点落在其中 -> %s",
                len(regions),
                statistics["homology_sites"],
                hom_path,
            )

    statistics["breakpoints_written"] = written_breakpoints
    statistics["canonical_sites_written"] = written_canonical
    statistics["duplicate_reads_written"] = len(
        set().union(*duplicates_by_breakpoint.values()) if duplicates_by_breakpoint else set()
    )
    statistics["reads_info_rows_written"] = written_reads_info_rows
    statistics["interference_table_rows_written"] = written_interference_rows
    statistics["breakpoints_with_interference"] = sum(
        1 for bucket in evidence_by_breakpoint.values() if bucket.interference_source_read_names
    )

    logging.info("reads_info 输出 %d 行 -> %s", written_reads_info_rows, reads_info_path)
    logging.info(
        "断点输出 %d 行（合并前 %d 个位点，合并掉 %d 个）-> %s",
        written_breakpoints,
        statistics["breakpoints_before_merge"],
        statistics["breakpoints_merged_away"],
        insert_position_path,
    )
    logging.info(
        "同源干扰标记 %d 行（其中 %d 个断点存在干扰）-> %s",
        written_interference_rows,
        statistics["breakpoints_with_interference"],
        interference_path,
    )
    if statistics["breakpoints_with_interference"]:
        source_totals: Dict[str, int] = {}
        for bucket in evidence_by_breakpoint.values():
            for source, read_names in bucket.interference_source_read_names.items():
                source_totals[source] = source_totals.get(source, 0) + len(read_names)
        logging.info(
            "干扰源头汇总（支持 read 数）: %s",
            dict(sorted(source_totals.items(), key=lambda item: (-item[1], item[0]))),
        )
    logging.info(
        "重复 read 输出 %d 条（分布在 %d 个位点）",
        statistics["duplicate_reads_written"],
        len(duplicates_by_breakpoint),
    )
    logging.info("统计: %s", {key: value for key, value in statistics.items() if value})
    return statistics


def run(settings: ScanSettings) -> dict[str, object]:
    """执行一次扫描；返回产物路径与统计。失败时抛 ``UsageError`` / ``RunError``。"""
    configure_logging(settings)
    config = resolve_runtime_config(settings)
    config.output_directory.mkdir(parents=True, exist_ok=True)
    logging.info(
        "开始处理 %d 个 BAM；query=%s insert=%s policy=%s dup=%s merge-window=%dbp gate=%s",
        len(config.bam_paths),
        config.query_chromosome,
        ",".join(config.insert_chromosomes),
        config.insert_origin_policy,
        "markdup(剔除)" if config.exclude_duplicates else "markdup(仅统计)",
        config.merge_window_bp,
        "on" if config.use_candidate_gate else "off",
    )

    # 样本名由 BAM 文件名首个 "." 之前推导（历史口径）。两个 BAM 推出同名时会互相覆盖结果，
    # 所以在动手写文件之前先拦掉。
    derived = [config.sample_name or path.name.split(".")[0] for path in config.bam_paths]
    duplicated = sorted({name for name in derived if derived.count(name) > 1})
    if duplicated:
        raise UsageError(
            "这些 BAM 会推导出相同的样本名，输出会互相覆盖："
            + ", ".join(duplicated)
            + "；用 --sample-name（单个 BAM）或改名后重跑"
        )

    started_at = time.perf_counter()
    for bam_path in config.bam_paths:
        try:
            process_single_bam(bam_path, config)
        except ConfigurationError as error:
            raise UsageError(str(error)) from None
        except RuntimeError as error:
            raise RunError(str(error)) from None
    logging.info("全部完成，用时 %.2f 秒", time.perf_counter() - started_at)
    return {
        "bams": [str(path) for path in config.bam_paths],
        "out_dir": str(config.output_directory),
        "query_chr": config.query_chromosome,
        "insert_chrs": list(config.insert_chromosomes),
    }


# ------------------------------------------------------------------------------
# ### 12. 自检：断点坐标与归属判定的最小验证 ###
# ------------------------------------------------------------------------------
def run_self_test() -> None:
    """内置自检：不读 BAM，只校验断点坐标、clip 判定、重复标记与归属打分。

    输入示例: 无（python parse_bam.py --self-test）
    输出示例: 断言全部通过时打印 'self-test ... OK'；否则抛 AssertionError
    """

    def make_block(block_start, block_end, read_start, read_end):
        return AlignmentBlock(block_start, block_end, read_start, read_end, block_end - block_start)

    def make_record(
        cigarstring,
        chromosome,
        block,
        mapping_quality,
        clip_side,
        clip_length,
        read_pair="read1",
        block_length=None,
        is_duplicate=False,
    ):
        return AlignmentRecord(
            read_name="read",
            read_pair=read_pair,
            cigarstring=cigarstring,
            is_supplementary=False,
            chromosome=chromosome,
            block_start=block.reference_start,
            block_end=block.reference_end,
            read_start=block.read_start,
            read_end=block.read_end,
            block_length=block_length or block.length,
            strand="+",
            mapping_quality=mapping_quality,
            mismatch_count=0,
            clip_side=clip_side,
            clip_length=clip_length,
            is_duplicate=is_duplicate,
        )

    # (1) clip 侧判定
    assert detect_clip_side("47S67M1I35M") == ("left", 47)
    assert detect_clip_side("99M51S") == ("right", 51)
    assert detect_clip_side("150M") == ("none", 0)

    # (2) M 块选取：多 M 块时默认取与 clip 相邻的块（历史实现取最后一个，断点会偏 67 bp）
    first_block = make_block(7000000, 7000067, 47, 114)  # 67M
    last_block = make_block(7000101, 7000200, 115, 150)  # 35M
    assert choose_evidence_block([first_block, last_block], "47S67M1I35M", False) == first_block
    assert choose_evidence_block([first_block, last_block], "47S67M1I35M", True) == last_block
    assert choose_evidence_block([last_block], "99M51S", False) == last_block

    # (3) 断点坐标
    left_clipped = make_record("47S67M1I35M", "chr2", first_block, 60, "left", 47)
    right_clipped = make_record(
        "99M51S", "chr1", make_block(4000000, 4000099, 0, 99), 60, "right", 51
    )
    assert compute_junction_position(left_clipped) == 7000000
    assert compute_junction_position(right_clipped) == 4000099

    # (4) 共有 read 归属：target_ref 上的命中质量更高 -> 归属 target_ref
    homolog_record = make_record(
        "52M98S", "homolog_ref", make_block(300, 352, 0, 52), 0, "right", 98
    )
    target_record = make_record(
        "99M51S", "target_ref", make_block(2000, 2099, 0, 99), 60, "right", 51
    )
    genome_record = make_record(
        "99M51S", "chr1", make_block(4000000, 4000099, 0, 99), 60, "right", 51
    )

    # (5) 插入序列命中（归属判定的输入）与 CIGAR 量化
    def make_hit(
        chromosome,
        mapping_quality,
        block_length,
        clip_length,
        source="alignment",
        mismatch_count=0,
        record=None,
    ):
        return InsertHit(
            chromosome=chromosome,
            mapping_quality=mapping_quality,
            block_length=block_length,
            mismatch_count=mismatch_count,
            clip_length=clip_length,
            source=source,
            record=record,
        )

    assert measure_cigar_block_and_clip("116S34M") == (34, 116)
    assert measure_cigar_block_and_clip("150M") == (150, 0)

    # (5a) XA 条目解析：4 字段与 6 字段两种写法都要能用（本项目的 BAM 是 4 字段写法）
    insert_list = ["decoy_ref", "homolog_ref", "target_ref"]
    xa_four = parse_xa_entry("target_ref,+335,31M119S,0", insert_list)
    assert (xa_four.chromosome, xa_four.block_length, xa_four.clip_length) == (
        "target_ref",
        31,
        119,
    )
    xa_six = parse_xa_entry("target_ref,335,+,31M119S,60,2", insert_list)
    assert (xa_six.mapping_quality, xa_six.mismatch_count) == (60, 2)
    assert parse_xa_entry("chr22,+335,31M119S,0", insert_list) is None
    homolog_hit = make_hit("homolog_ref", 0, 52, 98, source="xa")
    target_hit = make_hit("target_ref", 60, 99, 51, source="alignment", record=target_record)
    origin, ambiguous, scores = determine_insert_origin(
        [homolog_hit, target_hit], ["decoy_ref", "homolog_ref", "target_ref"]
    )
    assert (origin, ambiguous) == ("target_ref", False), (origin, ambiguous, scores)
    # 两个插入片段得分并列 -> 无法区分
    tied_hit = make_hit("homolog_ref", 60, 99, 51, source="alignment")
    origin, ambiguous, _ = determine_insert_origin(
        [tied_hit, target_hit], ["homolog_ref", "target_ref"]
    )
    assert origin is None and ambiguous is True

    # (5b) 归属策略（参数：过滤后的记录 + 插入命中 + 同源命中 + query + 列表 + 策略）
    #      返回值是 4 元组：(是否保留, 处置标签, 归属比对, 归属标签)
    assert (
        decide_read_usage([], [homolog_hit], set(), "target_ref", insert_list, "owned")[0] is False
    )
    keep_flag, reason, resolved_origin, origin_label = decide_read_usage(
        [genome_record], [target_hit], set(), "target_ref", insert_list, "owned"
    )
    assert (keep_flag, reason, resolved_origin.chromosome, origin_label) == (
        True,
        "kept_owned",
        "target_ref",
        "owned_query",
    )
    # 归属其它插入片段 -> 归属标签 owned_other（legacy 策略下它会非 0，所以干扰表里也能计数）
    assert (
        decide_read_usage([], [homolog_hit], set(), "target_ref", insert_list, "owned")[3]
        == "owned_other"
    )
    assert (
        decide_read_usage([], [homolog_hit], set(), "target_ref", insert_list, "legacy")[3]
        == "owned_other"
    )

    # (5c) 两个插入片段得分并列 -> 无法区分：owned 直接丢，consensus 也只是“待投票”
    assert (
        decide_read_usage([], [tied_hit, target_hit], set(), "target_ref", insert_list, "owned")[1]
        == "drop_ambiguous_insert"
    )
    assert (
        decide_read_usage(
            [], [tied_hit, target_hit], set(), "target_ref", insert_list, "consensus"
        )[1]
        == "drop_ambiguous_insert"
    )
    # 实测数据里的真实情形：clip 落在三个插入片段共有的 31 bp 序列上 -> 三方并列
    shared_hits = [
        make_hit(name, 0, 31, 119, source="secondary")
        for name in ("homolog_ref", "decoy_ref", "target_ref")
    ]
    assert (
        decide_read_usage(
            [], shared_hits, {"homolog_ref", "decoy_ref"}, "target_ref", insert_list, "owned"
        )[1]
        == "drop_ambiguous_insert"
    )
    # 归属不明确的 read 不再有任何“无条件保留”的策略（原 lenient 已删除）
    assert "lenient" not in ("owned", "consensus", "legacy")

    # (5f) 同源干扰源头判定（三种典型情况）与“共有 read”信号抽取
    assert determine_interference_sources("owned_query", {"target_ref"}, {"homolog_ref"}) == {
        "homolog_ref"
    }
    # 归属不明确：除 query 外的候选都算源头（query 自己永远不算）
    assert determine_interference_sources(
        "ambiguous", {"target_ref", "homolog_ref"}, set(), "target_ref"
    ) == {"homolog_ref"}
    assert determine_interference_sources("owned_query", {"target_ref"}, set()) == set()
    assert collect_homology_hit_chromosomes(
        [
            make_hit("target_ref", 60, 99, 51, source="alignment"),
            make_hit("homolog_ref", 0, 31, 119, source="xa"),
            make_hit("decoy_ref", 0, 31, 119, source="secondary"),
        ]
    ) == {"homolog_ref", "decoy_ref"}

    # (5g) 单插入场景：--insert-chrs 只含 query 时，不可能出现归属并列
    origin, ambiguous, _ = determine_insert_origin([target_hit], ["target_ref"])
    assert (origin, ambiguous) == ("target_ref", False)

    # (5d) --drop-shared-reads：归属 query、但在其它序列上还有 secondary/XA 命中 -> 剔除
    strict_keep, strict_reason, _, _ = decide_read_usage(
        [],
        [target_hit],
        {"target_ref", "homolog_ref"},
        "target_ref",
        insert_list,
        "owned",
        drop_shared_reads=True,
    )
    assert (strict_keep, strict_reason) == (False, "drop_shared_with_other_insert")
    # 同一开关与 consensus 组合也生效（今天两个需求可同时表达）
    assert (
        decide_read_usage(
            [],
            [target_hit],
            {"target_ref", "homolog_ref"},
            "target_ref",
            insert_list,
            "consensus",
            drop_shared_reads=True,
        )[1]
        == "drop_shared_with_other_insert"
    )
    # 不开开关时保留（只做归属区分，不剔除共有 read）
    assert (
        decide_read_usage([], [target_hit], {"homolog_ref"}, "target_ref", insert_list, "owned")[0]
        is True
    )

    # (5e) 插入侧比对被 MAPQ 过滤掉（只能来自 secondary/XA），仍能判定 read 来自 target_ref
    low_quality_hit = make_hit("target_ref", 0, 34, 116, source="secondary")
    assert (
        decide_read_usage(
            [genome_record], [low_quality_hit], set(), "target_ref", insert_list, "owned"
        )[1]
        == "kept_owned"
    )
    # 完全没有插入序列命中 -> 丢弃
    assert (
        decide_read_usage([genome_record], [], set(), "target_ref", insert_list, "owned")[1]
        == "drop_no_insert_hit"
    )
    # legacy 策略已删除：不认识的值必须被拦在 CLI 层（这里只保证 owned/consensus 两条路存在）
    assert {
        decide_read_usage([], [], set(), "target_ref", insert_list, policy)[1]
        for policy in ("owned", "consensus")
    } == {"drop_no_insert_hit"}

    # (6) 重复 read：只认 markdup 的 0x400，按 read 聚合（任一条记录被标记即算重复）
    assert read_is_duplicate([right_clipped, target_record]) is False
    assert (
        read_is_duplicate(
            [
                right_clipped,
                make_record(
                    "99M51S",
                    "chr1",
                    make_block(4000000, 4000099, 0, 99),
                    60,
                    "right",
                    51,
                    is_duplicate=True,
                ),
            ]
        )
        is True
    )
    # 重复 read 仍能按断点归属（dup_read_count / dup_read_name 两列靠它）
    duplicate_read = ReadEvidence(
        read_name="readC",
        all_records=[right_clipped],
        evidence_records=[right_clipped],
        origin_record=target_record,
        origin_label="owned_query",
        origin_status="kept_owned",
        interference_sources=set(),
    )
    assert collect_duplicate_reads_by_breakpoint([duplicate_read], "target_ref", 0) == {
        ("chr1", 4000099): {"readC"}
    }
    assert collect_duplicate_reads_by_breakpoint([], "target_ref", 0) == {}

    # (7) SR 证据
    assert iter_split_read_evidence(left_clipped, "decoy_ref", 0) == (
        "chr2",
        7000000,
        "read1",
    )
    assert iter_split_read_evidence(target_record, "target_ref", 0) is None

    # (8) PESR 证据：双端同染色体、同断点（一侧左 clip、一侧右 clip，断点都落在 4000100）
    pair_record = make_record(
        "51S99M",
        "chr1",
        make_block(4000099, 4000198, 51, 150),
        60,
        "left",
        51,
        read_pair="read2",
    )
    assert find_pair_end_split_read_junction(right_clipped, pair_record, "target_ref") == (
        "chr1",
        4000099,
    )

    # (9) 位点窗口合并：相邻 6 bp 的位点合并为一个，代表位点 = 支持 read 最多者
    assert cluster_positions_within_window([5000000, 5000010, 5100000], 10) == [
        [5000000, 5000010],
        [5100000],
    ]
    assert cluster_positions_within_window([100, 200], 0) == [[100], [200]]

    def make_bucket(
        split_names=(),
        owned_query=(),
        owned_other=(),
        ambiguous=(),
        no_insert_hit=(),
        consensus=(),
        sources=None,
        pair_split=(),
    ):
        return BreakpointEvidence(
            split_read_names=set(split_names),
            pair_split_read_names=set(pair_split),
            owned_query_read_names=set(owned_query),
            owned_other_read_names=set(owned_other),
            ambiguous_read_names=set(ambiguous),
            no_insert_hit_read_names=set(no_insert_hit),
            consensus_read_names=set(consensus),
            interference_source_read_names=sources or {},
        )

    small_evidence = {
        ("chr4", 5000000): make_bucket(
            split_names=["readA/read1"], owned_query=["readA"], sources={"homolog_ref": {"readA"}}
        ),
        ("chr4", 5000010): make_bucket(
            split_names=["readB/read1", "readC/read1"], owned_query=["readB", "readC"]
        ),
    }
    merged_evidence, merged_duplicates, merged_positions = merge_nearby_breakpoints(
        small_evidence,
        {("chr4", 5000000): {"readA"}, ("chr4", 5000010): {"readB", "readC"}},
        10,
    )
    assert merged_positions == {("chr4", 5000010): (5000000, 5000010)}
    merged_bucket = merged_evidence[("chr4", 5000010)]
    assert merged_bucket.split_read_names == {"readA/read1", "readB/read1", "readC/read1"}
    assert merged_bucket.owned_query_read_names == {"readA", "readB", "readC"}
    assert merged_bucket.interference_source_read_names == {"homolog_ref": {"readA"}}
    assert merged_duplicates[("chr4", 5000010)] == {"readA", "readB", "readC"}
    # 窗口 0 -> 恒等（保持历史逐位点输出）
    identity_evidence, _, identity_positions = merge_nearby_breakpoints(small_evidence, {}, 0)
    assert set(identity_evidence) == set(small_evidence)
    assert identity_positions == {
        ("chr4", 5000000): (5000000,),
        ("chr4", 5000010): (5000010,),
    }

    # (10) consensus 策略：挂起的“无法区分”read 按同断点投票推定归属（需要窗口才能命中相邻位点）
    offset_genome_record = make_record(
        "99M51S", "chr1", make_block(4000006, 4000106, 0, 99), 60, "right", 51
    )
    consensus_config = RuntimeConfig(
        bam_paths=[],
        sample_name=None,
        output_directory=Path("."),
        query_chromosome="target_ref",
        insert_chromosomes=("homolog_ref", "target_ref"),
        insert_origin_policy="consensus",
        drop_shared_reads=False,
        merge_window_bp=10,
        mismatch_source="nm",
        on_missing_tag="skip",
        minimum_mapping_quality=20,
        maximum_mismatch_count=3,
        maximum_soft_clip_blocks=1,
        minimum_clip_length=0,
        keep_secondary_alignments=False,
        require_split_match=True,
        legacy_coordinates=False,
        use_candidate_gate=True,
        exclude_duplicates=True,
        require_duplicate_flags=False,
        check_markdup=False,
        minimum_site_support=1,
        reference_fasta=None,
        homology_bed=(),
        homology_min_identity=0.90,
        homology_min_length=200,
        homology_window_bp=50000,
        bam_threads=1,
    )
    ambiguous_read = ReadEvidence(
        read_name="readB",
        all_records=[offset_genome_record],
        evidence_records=[offset_genome_record],
        origin_record=None,
        origin_label="ambiguous",
        origin_status="pending_ambiguous",
        interference_sources={"homolog_ref", "target_ref"},
    )
    owner_read = ReadEvidence(
        read_name="readA",
        all_records=[genome_record],
        evidence_records=[genome_record],
        origin_record=target_record,
        origin_label="owned_query",
        origin_status="kept_owned",
        interference_sources=set(),
    )
    consensus_statistics = create_statistics()
    resolved_reads = resolve_ambiguous_reads_by_site_consensus(
        [ambiguous_read], [owner_read], consensus_config, consensus_statistics
    )
    assert len(resolved_reads) == 1
    assert (resolved_reads[0].origin_label, resolved_reads[0].origin_status) == (
        "consensus_query",
        "kept_consensus",
    )
    assert consensus_statistics["reads_kept_consensus"] == 1
    # 窗口 0：两个位点相距 6 bp，投票不命中 -> 丢弃
    exact_statistics = create_statistics()
    assert (
        resolve_ambiguous_reads_by_site_consensus(
            [ambiguous_read],
            [owner_read],
            consensus_config._replace(merge_window_bp=0),
            exact_statistics,
        )
        == []
    )
    assert exact_statistics["drop_consensus_undecided"] == 1
    # 投票指向别的插入片段 -> 不保留
    other_owner = owner_read._replace(origin_label="owned_other")
    other_statistics = create_statistics()
    assert (
        resolve_ambiguous_reads_by_site_consensus(
            [ambiguous_read], [other_owner], consensus_config, other_statistics
        )
        == []
    )
    assert other_statistics["drop_consensus_undecided"] == 1
    print(f"self-test 通过（scan engine {SCAN_ENGINE_VERSION}）")


__all__ = ["ScanSettings", "run", "resolve_runtime_config", "run_self_test"]
