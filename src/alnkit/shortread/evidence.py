"""证据层：SR / PESR / DUP 三类证据的抽取与按位点聚合。"""

from __future__ import annotations

from typing import Dict, Iterator, Optional, Sequence, Set, Tuple

from .bam import compute_junction_position

from .model import AlignmentRecord, BreakpointEvidence, ReadEvidence, STAT_KEYS


def collect_duplicate_reads_by_breakpoint(
    duplicate_read_evidences: Sequence[ReadEvidence],
    query_chromosome: str,
    minimum_clip_length: int,
) -> Dict[Tuple[str, int], Set[str]]:
    """把被剔除的重复 read 按断点归属，供 dup_read_count / dup_read_name 两列使用。

    输入示例:
        duplicate_read_evidences = [readC（被 markdup 标记，支持 chr1:4000100）]
        query_chromosome = 'target_ref'
    输出示例: {('chr1', 4000100): {'readC'}}

    说明: 重复 read 默认被排除在证据之外（SR/PESR 不计它们），但“这些位点上有多少重复 read”
          仍然要能看见，所以单独聚合一次；断点定义与证据聚合完全一致（iter_read_site_evidence）。
    """
    duplicates_by_breakpoint: Dict[Tuple[str, int], Set[str]] = {}
    for read_evidence in duplicate_read_evidences:
        for breakpoint_key, _evidence_kind, _read_pair in iter_read_site_evidence(
            read_evidence.evidence_records, query_chromosome, minimum_clip_length
        ):
            duplicates_by_breakpoint.setdefault(breakpoint_key, set()).add(read_evidence.read_name)
    return duplicates_by_breakpoint


# ------------------------------------------------------------------------------
# ### 8. 断点证据统计：SR / PESR ###
# ------------------------------------------------------------------------------
def iter_split_read_evidence(
    record: AlignmentRecord, query_chromosome: str, minimum_clip_length: int
) -> Optional[Tuple[str, int, str]]:
    """判断单条比对是否构成 SR（split read）证据。

    输入示例:
        record(chromosome='chr2', cigarstring='47S67M1I35M', clip_side='left', clip_length=47)
        query_chromosome='decoy_ref'
    输出示例: ('chr2', 7000000, 'read2')
              落在 query_chromosome 上、没有 soft clip、或 clip 太短 -> None

    返回: (染色体, 断点位置, read1/read2) 或 None
    """
    if record.chromosome == query_chromosome:
        return None
    if record.clip_side == "none" or record.clip_length < minimum_clip_length:
        return None
    return record.chromosome, compute_junction_position(record), record.read_pair


def find_pair_end_split_read_junction(
    first_record: AlignmentRecord, second_record: AlignmentRecord, query_chromosome: str
) -> Optional[Tuple[str, int]]:
    """判断两条比对是否构成 PESR（双端同一断点）证据。

    输入示例:
        first_record  = 落在 chr2 上的比对（read1, '99M51S'）
        second_record = 落在 chr2 上的比对（read2, '51S99M'）
    输出示例: ('chr2', 4000100)      # 两侧断点相同
              read_pair 相同 / 染色体不同 / 落在 query 上 / 无 soft clip -> None

    说明: 与历史实现保持一致 —— 只有"双端落在同一条非目标染色体、且断点相同"才算。
    """
    if first_record.read_pair == second_record.read_pair:
        return None
    if first_record.chromosome != second_record.chromosome:
        return None
    if first_record.chromosome == query_chromosome:
        return None
    if first_record.clip_side == "none" or second_record.clip_side == "none":
        return None
    first_junction = compute_junction_position(first_record)
    second_junction = compute_junction_position(second_record)
    if first_junction != second_junction:
        return None
    return first_record.chromosome, first_junction


def iter_read_site_evidence(
    records: Sequence[AlignmentRecord], query_chromosome: str, minimum_clip_length: int
) -> Iterator[Tuple[Tuple[str, int], str, str]]:
    """产出一条 read 支持的全部“断点证据”（SR 与 PESR 的唯一实现入口）。

    输入示例:
        records = [chr1 上的 99M51S 比对（read2）, chr1 上的 51S99M 比对（read1）]
        query_chromosome = 'target_ref'
    输出示例:
        [(('chr1', 4000100), 'split', 'readA/read2'),
         (('chr1', 4000100), 'split', 'readA/read1'),
         (('chr1', 4000100), 'pair_split', 'readA')]
        -- 前两项是 SR（名字是 read_name/read_pair），第三项是 PESR（名字是 read_name）

    说明: 聚合统计与“归属投票”都调用它，保证两处用同一套断点定义。
    """
    # (1) SR：每条记录单独判定
    for record in records:
        split_evidence = iter_split_read_evidence(record, query_chromosome, minimum_clip_length)
        if split_evidence is None:
            continue
        chromosome, junction_position, read_pair = split_evidence
        yield (chromosome, junction_position), "split", read_pair

    # (2) PESR：两两组合判定（正常情况下每条 read 只有 2 条记录）
    for first_index in range(len(records)):
        for second_index in range(first_index + 1, len(records)):
            pair_evidence = find_pair_end_split_read_junction(
                records[first_index], records[second_index], query_chromosome
            )
            if pair_evidence is not None:
                yield pair_evidence, "pair_split", ""


def aggregate_breakpoint_evidence(
    read_evidences: Sequence[ReadEvidence], query_chromosome: str, minimum_clip_length: int
) -> Dict[Tuple[str, int], BreakpointEvidence]:
    """把全部 read 的证据聚合到断点上（含“同源干扰”标记与每个干扰源头的支持 read 数）。

    输入示例:
        read_evidences 含 2 条 read：readA（归属明确 = query）与 readB（无法区分，
        命中 target_ref + homolog_ref），两者都在 chr2:7000000 上产生 SR 证据
    输出示例:
        ({('chr2', 7000000): BreakpointEvidence(
             split_read_names={'readA/read2', 'readB/read2'},
             pair_split_read_names=set(),
             owned_query_read_names={'readA'},
             ambiguous_read_names={'readB'},
             interference_source_read_names={'target_ref': {'readA'}, 'homolog_ref': {'readB'}})},
         {('chr2', 7000000)})

    返回: 断点 -> BreakpointEvidence（断点集合就是它的键）
    """
    evidence_by_breakpoint: Dict[Tuple[str, int], BreakpointEvidence] = {}

    def bucket_for(breakpoint_key: Tuple[str, int]) -> BreakpointEvidence:
        if breakpoint_key not in evidence_by_breakpoint:
            evidence_by_breakpoint[breakpoint_key] = BreakpointEvidence(
                split_read_names=set(),
                pair_split_read_names=set(),
                owned_query_read_names=set(),
                owned_other_read_names=set(),
                ambiguous_read_names=set(),
                no_insert_hit_read_names=set(),
                consensus_read_names=set(),
                interference_source_read_names={},
            )
        return evidence_by_breakpoint[breakpoint_key]

    for read_evidence in read_evidences:
        sites_of_this_read: Set[Tuple[str, int]] = set()

        # (1)(2) SR / PESR 证据
        for breakpoint_key, evidence_kind, read_pair in iter_read_site_evidence(
            read_evidence.evidence_records, query_chromosome, minimum_clip_length
        ):
            sites_of_this_read.add(breakpoint_key)
            bucket = bucket_for(breakpoint_key)
            if evidence_kind == "split":
                bucket.split_read_names.add(f"{read_evidence.read_name}/{read_pair}")
            else:
                bucket.pair_split_read_names.add(read_evidence.read_name)

        # (3) 归属分布与同源干扰标记：每个被该 read 支持的断点都记一次
        for breakpoint_key in sites_of_this_read:
            bucket = bucket_for(breakpoint_key)
            if read_evidence.origin_label == "owned_query":
                bucket.owned_query_read_names.add(read_evidence.read_name)
            elif read_evidence.origin_label == "owned_other":
                bucket.owned_other_read_names.add(read_evidence.read_name)
            elif read_evidence.origin_label == "consensus_query":
                bucket.consensus_read_names.add(read_evidence.read_name)
            elif read_evidence.origin_label == "ambiguous":
                bucket.ambiguous_read_names.add(read_evidence.read_name)
            else:
                bucket.no_insert_hit_read_names.add(read_evidence.read_name)
            for source in read_evidence.interference_sources:
                bucket.interference_source_read_names.setdefault(source, set()).add(
                    read_evidence.read_name
                )

    return evidence_by_breakpoint


def collect_read_sites(
    records: Sequence[AlignmentRecord], query_chromosome: str, minimum_clip_length: int
) -> Set[Tuple[str, int]]:
    """给出一条 read 支持的全部断点（用于“归属投票”）。

    输入示例: 同 iter_read_site_evidence
    输出示例: {('chr1', 4000100)}
    """
    return {
        site
        for site, _, _ in iter_read_site_evidence(records, query_chromosome, minimum_clip_length)
    }


# ------------------------------------------------------------------------------
# ### 11. 主流程：逐 BAM 处理 ###
# ------------------------------------------------------------------------------
def create_statistics() -> Dict[str, int]:
    """初始化统计计数器。

    输入示例: 无
    输出示例: {'alignments_total': 0, 'drop_low_mapping_quality': 0, ..., 'reads_kept': 0}
    """
    return {key: 0 for key in STAT_KEYS}
