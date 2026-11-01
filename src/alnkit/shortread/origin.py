"""归属判定：一条 read 到底来自哪个插入片段（含 consensus 投票与两遍扫描）。"""

from __future__ import annotations

from collections import OrderedDict
from typing import Dict, List, Optional, Sequence, Set, Tuple
import pysam

from .bam import (
    build_alignment_record,
    choose_evidence_block,
    collect_insert_hits,
    count_mismatches,
    has_matching_split_alignment,
    mentions_any_contig,
    parse_alignment_blocks,
    read_is_duplicate,
    score_insert_hit,
)

from .evidence import collect_read_sites

from .interference import collect_homology_hit_chromosomes, determine_interference_sources

from .model import AlignmentRecord, BamReadCollection, InsertHit, ReadEvidence, RuntimeConfig


def determine_insert_origin(
    hits: Sequence[InsertHit], insert_chromosomes: Sequence[str]
) -> Tuple[Optional[str], bool, Dict[str, Tuple[int, int, int, int]]]:
    """判定一条 read 归属哪个插入片段（共有 read 的区分依据）。

    输入示例:
        hits = [InsertHit('chr2', 60, 99, 0, 51, r1),        # 基因组侧，不算候选
                InsertHit('homolog_ref', 0, 31, 0, 119, None),
                InsertHit('target_ref', 37, 150, 0, 0, r2)]
        insert_chromosomes = ['decoy_ref', 'homolog_ref', 'target_ref']
    输出示例:
        ('target_ref', False, {'homolog_ref': (0, 31, 0, -119), 'target_ref': (37, 150, 0, 0)})
        若多个插入片段得分并列 -> (None, True, {...})    # 无法区分（典型情形：clip 落在各插入片段共有的接头序列上）

    返回: (归属的插入序列名或 None, 是否无法区分, 每个候选插入序列的最高分)
    """
    best_score_by_chromosome: Dict[str, Tuple[int, int, int, int]] = {}
    for hit in hits:
        if hit.chromosome not in insert_chromosomes:
            continue
        score = score_insert_hit(hit)
        previous = best_score_by_chromosome.get(hit.chromosome)
        if previous is None or score > previous:
            best_score_by_chromosome[hit.chromosome] = score

    if not best_score_by_chromosome:
        return None, False, best_score_by_chromosome

    highest_score = max(best_score_by_chromosome.values())
    winners = [name for name, score in best_score_by_chromosome.items() if score == highest_score]
    if len(winners) > 1:
        return None, True, best_score_by_chromosome
    return winners[0], False, best_score_by_chromosome


def decide_read_usage(
    records: Sequence[AlignmentRecord],
    insert_hits: Sequence[InsertHit],
    homology_hit_chromosomes: Set[str],
    query_chromosome: str,
    insert_chromosomes: Sequence[str],
    policy: str,
) -> Tuple[bool, str, Optional[AlignmentRecord], str]:
    """按配置决定一条 read 是否参与本次（针对 --query-chr 的）断点统计。

    输入示例:
        records                 = 通过质量过滤的比对（用于后续断点统计）
        insert_hits             = 插入序列上的全部命中（含被质量过滤掉的、secondary、XA）
        homology_hit_chromosomes= {'homolog_ref'}   # 该 read 在其它序列上的 secondary/XA 命中
        query_chromosome        = 'target_ref'
    输出示例:
        policy='owned'        : 归属 target_ref -> (True,  'kept_owned', target_ref 上的比对, 'owned_query')
        policy='owned'        : 归属 homolog_ref -> (False, 'drop_other_insert', None, 'owned_other')
        policy='owned-strict' : 归属 target_ref 但另有同源命中
                                -> (False, 'drop_shared_with_other_insert', None, 'owned_query')
        policy='lenient'      : 归属不明确 -> (True,  'kept_ambiguous', None, 'ambiguous')
        policy='legacy'       : -> (True, 'kept_legacy', 归属比对或 None, 归属标签)
        （完全没有插入命中的 read 归属标签为 'no_insert_hit'）

    返回: (是否保留, 处置标签, 归属插入序列上的比对记录或 None, 归属标签)
          归属标签与过滤策略无关，供统计/干扰表使用。
    """
    origin_chromosome, ambiguous, _ = determine_insert_origin(insert_hits, insert_chromosomes)
    origin_record = None
    if origin_chromosome is not None:
        candidates_with_record = [
            hit
            for hit in insert_hits
            if hit.chromosome == origin_chromosome and hit.record is not None
        ]
        if candidates_with_record:
            origin_record = max(candidates_with_record, key=score_insert_hit).record

    if ambiguous:
        origin_label = "ambiguous"
    elif origin_chromosome is None:
        origin_label = "no_insert_hit"
    elif origin_chromosome == query_chromosome:
        origin_label = "owned_query"
    else:
        origin_label = "owned_other"

    if policy == "legacy":
        return True, "kept_legacy", origin_record, origin_label
    if ambiguous:
        if policy == "lenient":
            return True, "kept_ambiguous", origin_record, origin_label
        return False, "drop_ambiguous_insert", None, origin_label
    if origin_chromosome is None:
        return False, "drop_no_insert_hit", None, origin_label
    if origin_chromosome != query_chromosome:
        return False, "drop_other_insert", None, origin_label
    # 严格模式：即使归属 query，只要在其它插入序列上还有 secondary/XA 命中，也当作共有 read 剔除
    if policy == "owned-strict" and homology_hit_chromosomes - {query_chromosome}:
        return False, "drop_shared_with_other_insert", None, origin_label
    return True, "kept_owned", origin_record, origin_label


def resolve_ambiguous_reads_by_site_consensus(
    pending_ambiguous_reads: Sequence[ReadEvidence],
    read_evidences: Sequence[ReadEvidence],
    config: RuntimeConfig,
    statistics: Dict[str, int],
) -> List[ReadEvidence]:
    """用“同一断点上归属明确的 read”投票，推定无法区分归属的 read（consensus 策略）。

    输入示例:
        pending_ambiguous_reads = [readB（命中 target_ref 与 homolog_ref，无法区分）, ...]
        read_evidences          = 已保留的 read（其中归属明确者的断点参与投票）
        config.merge_window_bp  = 10
    输出示例:
        - readB 支持的断点 chr4:5000000 在 ±10 bp 内有 8 条 owned_query 投票、无相反投票
          -> 保留，origin_label='consensus_query'、origin_status='kept_consensus'
        - 某 read 的断点附近投票为空、或同时有 owned_query 与 owned_other -> 丢弃
          （统计 drop_consensus_undecided）

    说明: 这是**推断**而不是测量：只当该 read 支持的全部断点附近投票一致指向 query 时才保留，
          并将其干扰源头标记原样保留（方便人工复核）。窗口 = --merge-window-bp（0 = 精确匹配）。
    """
    if not pending_ambiguous_reads:
        return []

    # (1) 用归属明确的 read 建投票表：染色体 -> [(断点, 归属标签), ...]
    votes_by_chromosome: Dict[str, List[Tuple[int, str]]] = {}
    for read_evidence in read_evidences:
        if read_evidence.origin_label not in ("owned_query", "owned_other"):
            continue
        for chromosome, position in collect_read_sites(
            read_evidence.evidence_records, config.query_chromosome, config.minimum_clip_length
        ):
            votes_by_chromosome.setdefault(chromosome, []).append(
                (position, read_evidence.origin_label)
            )

    # (2) 逐条推定
    resolved_reads: List[ReadEvidence] = []
    for read_evidence in pending_ambiguous_reads:
        owner_labels: Set[str] = set()
        for chromosome, position in collect_read_sites(
            read_evidence.evidence_records, config.query_chromosome, config.minimum_clip_length
        ):
            for vote_position, vote_label in votes_by_chromosome.get(chromosome, ()):
                if abs(vote_position - position) <= config.merge_window_bp:
                    owner_labels.add(vote_label)
        if owner_labels == {"owned_query"}:
            statistics["reads_kept"] += 1
            statistics["reads_kept_consensus"] += 1
            resolved_reads.append(
                read_evidence._replace(
                    origin_label="consensus_query", origin_status="kept_consensus"
                )
            )
        else:
            statistics["drop_consensus_undecided"] += 1
    return resolved_reads


def select_evidence_records(
    records: Sequence[AlignmentRecord],
    query_chromosome: str,
    insert_chromosomes: Sequence[str],
    policy: str,
) -> List[AlignmentRecord]:
    """挑出参与断点统计的比对记录。

    输入示例:
        records = [chr2 上的比对, target_ref 上的比对]   # query_chromosome='target_ref'
    输出示例:
        policy != 'legacy' -> [chr2 上的比对]          # target_ref 自身坐标不算断点，另一个插入片段被剔除
        policy == 'legacy' -> [chr2 上的比对, target_ref 上的比对]

    说明: 非 legacy 策略下，落在"其它插入序列"上的比对不再计入断点，
          避免把 A 插入片段的 read 记成 B 插入片段的断点证据。
    """
    if policy == "legacy":
        return list(records)
    other_insert_chromosomes = set(insert_chromosomes) - {query_chromosome}
    return [record for record in records if record.chromosome not in other_insert_chromosomes]


def discover_candidate_reads(
    alignment_file: "pysam.AlignmentFile", config: RuntimeConfig, statistics: Dict[str, int]
) -> Set[str]:
    """第 1 遍：只读标签，筛出“可能相关”的 read 名（大 BAM 提速的关键）。

    输入示例:
        BAM 中两条比对:
          A: query_name='r1', reference_name='chr4', SA缺失, XA='target_ref,+335,31M119S,0;'
          B: query_name='r2', reference_name='chr7', SA缺失, XA缺失
    输出示例: {'r1'}      # r2 既不碰插入序列、也没有 SA/XA，永远不可能产生证据

    候选判据:
        ① 比对落在 --insert-chrs 的 contig 上
        ② 某条比对的 SA 里出现任一插入 contig 名
        ③ 某条比对的 XA 里出现任一插入 contig 名

    说明:
        - 必须按 **read 名**收候选，不能“逐条比对”决定丢不丢：一条 read 的基因组侧比对
          （带 clip）本身可能既无 SA 也无 XA，只因为它的 mate 落在插入 contig 上才相关。
          实测（本仓库 BAM）逐记录闸门丢 7/30 个位点，逐 read 闸门 0 丢——那 7 个位点的
          共同形态是「基因组侧 `3S147M`/`23S127M` 等无 SA/XA 的 clip 比对 + mate 落在插入 contig 上」。
        - **不把“有 soft clip”当成候选条件**：证据的判定在 `decide_read_usage()`，
          没有任何插入命中的 read 在 owned/lenient/consensus/owned-strict 下一律丢弃
          （legacy 另外关掉闸门），所以“只有 clip、没有插入线索”的 read 永远不可能成为证据，
          加进来只多读、不会多出结果（实测：11 万条 12S138M 无关 read 会让候选从 1.1 万涨到 15.5 万，结果不变）。
        - 顺手统计 is_duplicate（0x400）与光学重复（dt/DT），供“BAM 是否做过 markdup”判断。
        - 子串判定可能因 contig 名互为子串而多收几条 read —— 只会多看，不会漏。

    返回: 候选 read 名集合
    """
    candidate_reads: Set[str] = set()
    insert_names = set(config.insert_chromosomes)
    # until_eof=True：顺序遍历全文件，**不需要 .bai 索引**（pysam 的无参 fetch() 在无索引 BAM 上会报错）
    for alignment in alignment_file.fetch(until_eof=True):
        if alignment.is_unmapped or alignment.reference_name is None:
            continue
        statistics["pass1_alignments"] += 1
        if alignment.is_duplicate:
            statistics["duplicate_alignments"] += 1
            if alignment.has_tag("dt") or alignment.has_tag("DT"):
                statistics["optical_duplicate_alignments"] += 1
        read_name = alignment.query_name
        if alignment.reference_name in insert_names:
            candidate_reads.add(read_name)
            continue
        if alignment.has_tag("SA") and mentions_any_contig(alignment.get_tag("SA"), insert_names):
            candidate_reads.add(read_name)
            continue
        if alignment.has_tag("XA") and mentions_any_contig(alignment.get_tag("XA"), insert_names):
            candidate_reads.add(read_name)
    statistics["candidate_reads"] = len(candidate_reads)
    return candidate_reads


def collect_read_records(
    alignment_file: "pysam.AlignmentFile",
    config: RuntimeConfig,
    statistics: Dict[str, int],
    candidate_reads: Optional[Set[str]] = None,
) -> BamReadCollection:
    """遍历 BAM，按 read 名收集“通过过滤的比对记录”、“归属候选”与“共有 read 信号”。

    输入示例:
        BAM 中一条比对: query_name='LH...:2416:17423:10304', cigarstring='116S34M',
                        reference_name='target_ref', mapping_quality=0, is_secondary=True
        config.minimum_mapping_quality = 20
    输出示例:
        BamReadCollection(
            records_by_read_name={'LH...:1152:49645:20099': [chr2 上的比对]},
            insert_hits_by_read_name={'LH...:2416:17423:10304':
                [InsertHit('target_ref', 0, 34, 0, 116, 'secondary', record)]})
        解读：mapq=0 且 secondary 的 target_ref 比对不能作为断点证据（不进第一个字段），
              但它是归属判定与“共有 read”信号的依据。

    返回: BamReadCollection（两份字典，read 顺序与 BAM 读取顺序一致）

    参数 candidate_reads: 第 1 遍得到的候选 read 名；不为 None 时只处理这些 read 的比对
        （不仅不会丢证据，而且**必须**把它们的全部比对都收进来，含 mate 与次要比对）。
    """
    records_by_read_name: OrderedDict[str, List[AlignmentRecord]] = OrderedDict()
    insert_hits_by_read_name: Dict[str, List[InsertHit]] = {}

    for alignment in alignment_file.fetch(until_eof=True):
        statistics["alignments_total"] += 1
        if alignment.is_unmapped or alignment.reference_name is None:
            continue
        read_name = alignment.query_name
        if candidate_reads is not None and read_name not in candidate_reads:
            statistics["skipped_by_candidate_gate"] += 1
            continue

        # (0) 顺便统计比对特征，用于事后检查这个 BAM 的比对参数是否满足脚本要求
        #     （SA / soft clip / 补充比对 / 次要比对 / hard clip，详见 warn_about_alignment_features）
        if alignment.has_tag("SA"):
            statistics["alignments_with_sa"] += 1
        if alignment.has_tag("XA"):
            statistics["alignments_with_xa"] += 1
        if alignment.is_supplementary:
            statistics["supplementary_alignments"] += 1
        if alignment.is_secondary:
            statistics["secondary_alignments"] += 1
        if "H" in (alignment.cigarstring or ""):
            statistics["hard_clipped_alignments"] += 1

        # (1) 插入序列命中（归属判定与“共有 read”信号用，不受质量过滤影响）：
        #     落在插入序列上的比对（含 secondary）+ 该比对的 XA 备选命中
        own_record_for_hits = None
        if alignment.reference_name in config.insert_chromosomes:
            own_blocks = parse_alignment_blocks(alignment)
            if own_blocks:
                own_block = choose_evidence_block(
                    own_blocks, alignment.cigarstring or "", config.legacy_coordinates
                )
                own_mismatch = count_mismatches(
                    alignment, config.mismatch_source, config.on_missing_tag, statistics=None
                )
                # 缺 NM/MD 时 count_mismatches 返回 None。绝不能当 0 用——0 是"最完美"，
                # 会让这条比对在归属打分里胜出。用一个大数表示"没有信息"，只参与记录不争取归属。
                own_record_for_hits = build_alignment_record(
                    alignment, own_block, own_mismatch if own_mismatch is not None else 10**9
                )
        alignment_insert_hits = collect_insert_hits(
            alignment, config.insert_chromosomes, own_record_for_hits
        )
        if alignment_insert_hits:
            insert_hits_by_read_name.setdefault(read_name, []).extend(alignment_insert_hits)

        if alignment.mapping_quality <= config.minimum_mapping_quality:
            statistics["drop_low_mapping_quality"] += 1
            continue
        if alignment.is_secondary and not config.keep_secondary_alignments:
            statistics["drop_secondary"] += 1
            continue

        cigarstring = alignment.cigarstring or ""
        if cigarstring.count("S") > config.maximum_soft_clip_blocks:
            statistics["drop_too_many_soft_clips"] += 1
            continue

        blocks = parse_alignment_blocks(alignment)
        if not blocks:
            statistics["drop_too_many_soft_clips"] += 1
            continue

        mismatch_count = count_mismatches(
            alignment, config.mismatch_source, config.on_missing_tag, statistics
        )
        if mismatch_count is None:
            continue
        if mismatch_count > config.maximum_mismatch_count:
            statistics["drop_too_many_mismatches"] += 1
            continue

        if config.require_split_match and alignment.reference_name == config.query_chromosome:
            if not has_matching_split_alignment(alignment):
                statistics["drop_failed_split_match"] += 1
                continue

        block = choose_evidence_block(blocks, cigarstring, config.legacy_coordinates)
        record = build_alignment_record(alignment, block, mismatch_count)
        records_by_read_name.setdefault(record.read_name, []).append(record)

    return BamReadCollection(
        records_by_read_name=records_by_read_name,
        insert_hits_by_read_name=insert_hits_by_read_name,
    )


def select_read_evidences(
    read_collection: BamReadCollection, config: RuntimeConfig, statistics: Dict[str, int]
) -> Tuple[List[ReadEvidence], List[ReadEvidence]]:
    """对每条 read 做归属判定，产出参与统计的 ReadEvidence 列表。

    输入示例:
        records_by_read_name = {'LH...:1152:49645:20099': [chr2 记录, target_ref 记录],
                                'LH...:2416:17423:10304': [target_ref 记录]}
        config = query_chromosome='target_ref', insert_origin_policy='owned'
    输出示例:
        ([ReadEvidence(read_name='LH...:1152:49645:20099', origin_label='owned_query',
                       evidence_records=[chr2 记录], ...)],
         [ReadEvidence(read_name='LH...:2416:17423:10304', ...)])   # 被 markdup 标记的重复 read

    说明（单记录 read 的处理）：
        - --insert-origin-policy legacy 保留全部历史规则：比对记录数 <2 的 read 丢弃；
        - 新策略下不适用该规则：插入侧比对常常因为 MAPQ 被过滤，只剩基因组侧那一条带 clip 的
          比对，但它仍是有效的断点证据（归属由“插入命中”判定，不依赖该比对是否存在）。
    说明（重复 read 的处理）：
        - 只有“已通过归属判定”的 read 才可能被当作重复 read 剔除（否则 dup 列会把本该因归属
          被丢的 read 也算进去）；consensus 推定归属的 read 同样参与这个判断。
        - 剔除的 read 单独返回：它们不计入 SR/PESR，但仍用于填充 dup_read_count / dup_read_name。

    返回: (参与统计的 read, 被剔除的重复 read)，两者顺序与 BAM 一致
    """
    read_evidences: List[ReadEvidence] = []
    pending_ambiguous_reads: List[ReadEvidence] = []
    for read_name, records in read_collection.records_by_read_name.items():
        # (a) 历史规则（仅 --insert-origin-policy legacy）：比对记录数 <2 的 read 丢弃
        if config.insert_origin_policy == "legacy" and len(records) < 2:
            statistics["drop_single_alignment"] += 1
            continue
        # (b) 历史规则（仅 --insert-origin-policy legacy）：恰好两条比对、且 is_supplementary
        #     与 M 块起始坐标都相同，视为普通 read 丢弃
        #     （注意：历史实现的索引是 read_info[2]=is_supplementary 与 read_info[4]=M_ref_pos_start，
        #      注释里写的 "cigarstring 与 chrom" 与实际不符）
        if config.insert_origin_policy == "legacy" and len(records) == 2:
            first, second = records
            if (
                first.is_supplementary == second.is_supplementary
                and first.block_start == second.block_start
            ):
                statistics["drop_plain_pair"] += 1
                continue
        # (c) 插入片段归属判定（共有 read 的区分）
        #     归属判定使用 insert_hits（落在插入序列上的 primary/supplementary/secondary 比对 + XA 命中，
        #     包含被质量过滤掉的），避免“插入侧比对因 MAPQ/错配被过滤 -> 被误判为不涉及任何插入片段”。
        insert_hits = read_collection.insert_hits_by_read_name.get(read_name, [])
        homology_hit_chromosomes = collect_homology_hit_chromosomes(insert_hits)
        if homology_hit_chromosomes:
            statistics["reads_shared_with_other_insert"] += 1
        keep_read, reason, origin_record, origin_label = decide_read_usage(
            records,
            insert_hits,
            homology_hit_chromosomes,
            config.query_chromosome,
            config.insert_chromosomes,
            config.insert_origin_policy,
        )
        candidate_insert_chromosomes = {
            hit.chromosome for hit in insert_hits if hit.chromosome in config.insert_chromosomes
        }
        if not keep_read:
            if config.insert_origin_policy == "consensus" and reason == "drop_ambiguous_insert":
                # 先挂起，等全部 read 分类完再用“同断点投票”推定归属
                pending_ambiguous_reads.append(
                    ReadEvidence(
                        read_name=read_name,
                        all_records=list(records),
                        evidence_records=select_evidence_records(
                            records,
                            config.query_chromosome,
                            config.insert_chromosomes,
                            config.insert_origin_policy,
                        ),
                        origin_record=origin_record,
                        origin_label=origin_label,
                        origin_status="pending_ambiguous",
                        interference_sources=determine_interference_sources(
                            origin_label, candidate_insert_chromosomes, homology_hit_chromosomes
                        ),
                    )
                )
                continue
            statistics[reason] = statistics.get(reason, 0) + 1
            continue
        statistics["reads_kept"] += 1
        if reason == "kept_ambiguous":
            statistics["reads_kept_ambiguous"] += 1
        read_evidences.append(
            ReadEvidence(
                read_name=read_name,
                all_records=list(records),
                evidence_records=select_evidence_records(
                    records,
                    config.query_chromosome,
                    config.insert_chromosomes,
                    config.insert_origin_policy,
                ),
                origin_record=origin_record,
                origin_label=origin_label,
                origin_status=reason,
                interference_sources=determine_interference_sources(
                    origin_label, candidate_insert_chromosomes, homology_hit_chromosomes
                ),
            )
        )

    # consensus 策略：挂起的“归属不明确”read 按同断点归属投票推定（推断，不是测量）
    if pending_ambiguous_reads:
        resolved_reads = resolve_ambiguous_reads_by_site_consensus(
            pending_ambiguous_reads, read_evidences, config, statistics
        )
        read_evidences.extend(resolved_reads)

    # 重复 read：只认 BAM 的 markdup 标记
    #   不论是否剔除，都单独返回一份 —— dup_read_count / dup_read_name 两列始终有意义
    kept_read_evidences: List[ReadEvidence] = []
    duplicate_read_evidences: List[ReadEvidence] = []
    for read_evidence in read_evidences:
        target = (
            duplicate_read_evidences
            if read_is_duplicate(read_evidence.all_records)
            else kept_read_evidences
        )
        target.append(read_evidence)
    statistics["duplicate_reads"] = len(duplicate_read_evidences)
    if config.exclude_duplicates:
        statistics["drop_duplicate_read"] = len(duplicate_read_evidences)
        statistics["reads_kept"] -= len(duplicate_read_evidences)
        return kept_read_evidences, duplicate_read_evidences
    return read_evidences, duplicate_read_evidences
