"""短读长链路的共享数据类型与常量（对齐记录、命中、证据包、配置）。"""

from __future__ import annotations

import re
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Set, Tuple


SCAN_ENGINE_VERSION = "3.0.0"


# 热路径上的正则预先编译（底层 re 有 512 条缓存，仍会走一次 cache 查找；这里的调用量是“每条比对若干次”）
LEFT_CLIP_RE = re.compile(r"(\d+)S\d+M")


RIGHT_CLIP_RE = re.compile(r"(\d+)M(\d+)S$")


CIGAR_TOKEN_RE = re.compile(r"(\d+)([MIDSHX=])")


CIGAR_FULL_RE = re.compile(r"(?:\d+[MIDSHX=])+")


# 三类证据的日志/统计键名，集中定义，避免散落的字符串
STAT_KEYS = (
    "alignments_total",
    "alignments_with_sa",
    "alignments_with_xa",
    "supplementary_alignments",
    "secondary_alignments",
    "hard_clipped_alignments",
    "drop_low_mapping_quality",
    "drop_secondary",
    "drop_too_many_soft_clips",
    "drop_missing_tag",
    "drop_too_many_mismatches",
    "drop_failed_split_match",
    "drop_single_alignment",
    "drop_plain_pair",
    "drop_no_insert_hit",
    "drop_other_insert",
    "drop_ambiguous_insert",
    "drop_shared_with_other_insert",
    "reads_shared_with_other_insert",
    "reads_kept",
    "reads_kept_ambiguous",
    "reads_kept_consensus",
    "drop_consensus_undecided",
    "breakpoints_with_interference",
    "interference_table_rows_written",
    "breakpoints_before_merge",
    "breakpoints_merged_away",
    "breakpoints_written",
    "duplicate_alignments",
    "duplicate_reads",
    "duplicate_reads_written",
    "optical_duplicate_alignments",
    "drop_duplicate_read",
    "candidate_reads",
    "skipped_by_candidate_gate",
    "pass1_alignments",
    "pass1_ms",
    "pass2_ms",
    "reads_info_rows_written",
)


# ------------------------------------------------------------------------------
# ### 2. 运行配置：命令行参数 ###
# ------------------------------------------------------------------------------
class RuntimeConfig(NamedTuple):
    """一次运行的全部配置（由命令行参数解析而来）。

    输入示例（命令行）:
        --bam-dir ./bam --query-chr target_ref --insert-chrs decoy_ref,homolog_ref,target_ref
    本结构示例:
        RuntimeConfig(bam_paths=[...], query_chromosome='target_ref',
                      insert_chromosomes=('decoy_ref','homolog_ref','target_ref'), ...)
    """

    bam_paths: List[Path]
    sample_name: Optional[str]
    output_directory: Path
    query_chromosome: str
    insert_chromosomes: Tuple[str, ...]
    insert_origin_policy: str
    merge_window_bp: int
    mismatch_source: str
    on_missing_tag: str
    minimum_mapping_quality: int
    maximum_mismatch_count: int
    maximum_soft_clip_blocks: int
    minimum_clip_length: int
    keep_secondary_alignments: bool
    require_split_match: bool
    legacy_coordinates: bool
    use_candidate_gate: bool
    exclude_duplicates: bool
    require_duplicate_flags: bool
    bam_threads: int


@dataclass(slots=True)
class ScanSettings:
    """``alnkit short scan`` 的全部入参（与旧脚本参数一一对应）。"""

    bam: tuple[Path, ...] = ()
    bam_dir: Path = Path(".")
    bam_glob: str = "*.sorted.bam"
    no_recursive: bool = False
    out_dir: Path = Path(".")
    sample_name: str | None = None
    query_chr: str | None = None
    insert_chrs: str | None = None
    insert_origin_policy: str = "owned"
    merge_window_bp: int = 0
    min_mapq: int = 20
    max_mismatch: int = 3
    max_softclip_blocks: int = 1
    min_clip_length: int = 0
    keep_secondary: bool = False
    mismatch_source: str = "nm"
    on_missing_tag: str = "skip"
    no_split_match: bool = False
    no_candidate_gate: bool = False
    no_duplicate_filter: bool = False
    require_dup_flags: bool = False
    bam_threads: int = 1
    legacy_coordinates: bool = False
    legacy_profile: bool = False
    log_level: str = "INFO"
    quiet: bool = False


# ------------------------------------------------------------------------------
# ### 3. 数据结构：一条比对记录与一条 read 的证据包 ###
# ------------------------------------------------------------------------------
class AlignmentBlock(NamedTuple):
    """CIGAR 中一个 M/=/X 比对块（参考坐标 0-based 半开区间）。

    输入示例: cigarstring='47S67M1I35M'
    输出示例:
        AlignmentBlock(reference_start=7000000, reference_end=7000100, read_start=47, read_end=114, length=67)
        AlignmentBlock(reference_start=7000101, reference_end=7000200, read_start=115, read_end=150, length=35)
    """

    reference_start: int
    reference_end: int
    read_start: int
    read_end: int
    length: int


class AlignmentRecord(NamedTuple):
    """一条保留下来、要参与后续统计的比对（对应 reads_info 的一行）。

    输入示例: 一条 pysam 比对（read_name='LH...:1152:49645:20099', read_pair='read2',
              cigarstring='47S67M1I35M', chromosome='chr2', mapping_quality=60）
    输出示例:
        AlignmentRecord(read_name='LH...:1152:49645:20099', read_pair='read2',
                        cigarstring='47S67M1I35M', is_supplementary=False, chromosome='chr2',
                        block_start=7000000, block_end=7000100, read_start=47, read_end=114,
                        block_length=67, strand='+', mapping_quality=60, mismatch_count=0,
                        clip_side='left', clip_length=47)
    """

    read_name: str
    read_pair: str
    cigarstring: str
    is_supplementary: bool
    chromosome: str
    block_start: int
    block_end: int
    read_start: int
    read_end: int
    block_length: int
    strand: str
    mapping_quality: int
    mismatch_count: int
    clip_side: str
    clip_length: int
    is_duplicate: bool = (
        False  # 上游 markdup 的 FLAG 0x400；只用于重复 read 统计，不进 reads_info 表
    )


class ReadEvidence(NamedTuple):
    """一条 read（= 一个 read_name，可能含 read1/read2 多条比对）的汇总信息。

    输入示例: read_name='LH...:1152:49645:20099'，两条比对分别落在 chr2 与 target_ref 上，
              另在 homolog_ref 上有 secondary 命中
    输出示例:
        ReadEvidence(read_name='LH...:1152:49645:20099',
                     all_records=[chr2 上的比对, target_ref 上的比对],
                     evidence_records=[chr2 上的比对],
                     origin_record=target_ref 上的比对,
                     origin_label='owned_query',
                     origin_status='kept_owned',
                     interference_sources={'homolog_ref'})

    字段说明:
        origin_label         归属标签（**与过滤策略无关**，供统计/干扰表使用）：
                             owned_query / owned_other / ambiguous / no_insert_hit
        origin_status        本次运行的处置标签（kept_owned / kept_ambiguous / kept_legacy）
        interference_sources 与该 read 同源/共有的其它序列名（判据见 determine_interference_sources）
    """

    read_name: str
    all_records: List[AlignmentRecord]
    evidence_records: List[AlignmentRecord]
    origin_record: Optional[AlignmentRecord]
    origin_label: str
    origin_status: str
    interference_sources: Set[str]


class BreakpointEvidence(NamedTuple):
    """某个断点上的全部证据与同源干扰标记。

    输入示例: 断点 ('chr1', 4000100) 受 3 条不同归属状态的支持 read 支持
    输出示例:
        BreakpointEvidence(split_read_names={'readA/read2'},
                           pair_split_read_names=set(),
                           owned_query_read_names={'readA'},
                           owned_other_read_names={'readB'},
                           ambiguous_read_names={'readC'},
                           no_insert_hit_read_names=set(),
                           interference_source_counts={'homolog_ref': 1, 'target_ref': 1})

    字段说明:
        split_read_names / pair_split_read_names  历史表格用的名字（read_name/read_pair 与 read_name）
        owned_query_read_names                    归属明确 = query 的支持 read
        owned_other_read_names                    归属明确 = 其它插入片段的支持 read（legacy 策略下会非 0）
        ambiguous_read_names                      归属无法区分的支持 read（lenient 保留，本字段用于标记）
        no_insert_hit_read_names                  未命中任何插入序列的支持 read（legacy 策略下会非 0）
        consensus_read_names                      由“同断点归属投票”推定归属的支持 read（consensus 策略）
        interference_source_read_names            干扰源头 -> 与之同源/共有的支持 read 名单
                                                  （存名单而不是计数，合并位点时可精确取并集）
    """

    split_read_names: Set[str]
    pair_split_read_names: Set[str]
    owned_query_read_names: Set[str]
    owned_other_read_names: Set[str]
    ambiguous_read_names: Set[str]
    no_insert_hit_read_names: Set[str]
    consensus_read_names: Set[str]
    interference_source_read_names: Dict[str, Set[str]]


class InsertHit(NamedTuple):
    """一条“落在插入序列上的命中”，仅用于判定 read 归属（不要求是完整比对记录）。

    输入示例（三类来源）:
        1) 一条 mapq=37 的 target_ref 150M 比对（primary）
        2) 一条 mapq=0  的 homolog_ref 31M119S 比对（secondary）
        3) XA 条目 'target_ref,+335,31M119S,0'
    输出示例:
        InsertHit(chromosome='target_ref', mapping_quality=37, block_length=150,
                  mismatch_count=0, clip_length=0, source='alignment', record=AlignmentRecord(...))
        InsertHit(chromosome='homolog_ref', mapping_quality=0, block_length=31,
                  mismatch_count=0, clip_length=119, source='secondary', record=AlignmentRecord(...))
        InsertHit(chromosome='target_ref', mapping_quality=0, block_length=31,
                  mismatch_count=0, clip_length=119, source='xa', record=None)   # 只存在于 XA 中
    """

    chromosome: str
    mapping_quality: int
    block_length: int
    mismatch_count: int
    clip_length: int
    source: str
    record: Optional[AlignmentRecord]


class BamReadCollection(NamedTuple):
    """一个 BAM 的两份收集结果。

    输入示例: 遍历一个 BAM 后得到下列两份信息
    输出示例:
        BamReadCollection(
            records_by_read_name={'readA': [chr1 上的比对, target_ref 上的比对]},
            insert_hits_by_read_name={'readA': [InsertHit(target_ref, 37, 150, 0, 0, 'alignment', record),
                                                   InsertHit(homolog_ref, 0, 31, 0, 119, 'xa', None)]})

    字段说明:
        records_by_read_name     通过质量过滤的比对（作为断点证据的候选）
        insert_hits_by_read_name 所有插入序列上的命中（含被质量过滤掉的、secondary 与 XA）；
                                 “共有 read（同源命中）”集合由它推导，不另外存一份
    """

    records_by_read_name: "OrderedDict[str, List[AlignmentRecord]]"
    insert_hits_by_read_name: Dict[str, List[InsertHit]]


class ConfigurationError(Exception):
    """配置与 BAM 不匹配（如 --query-chr 不在 BAM header 里）：由 main 转成退出码 2。"""
