"""位点窗口合并：相邻断点归并、代表位点选择。"""

from __future__ import annotations

from typing import Dict, List, Sequence, Set, Tuple

from .model import BreakpointEvidence


# ------------------------------------------------------------------------------
# ### 9. 位点窗口合并（同一个断点的 clip 边界摆动）###
# ------------------------------------------------------------------------------
def cluster_positions_within_window(positions: Sequence[int], window_bp: int) -> List[List[int]]:
    """把同一条染色体上的位点按窗口聚成若干簇。

    输入示例: positions=[5000000, 5000010, 5100000], window_bp=10
    输出示例: [[5000000, 5000010], [5100000]]

    说明: 贪心聚类，**锚点 = 簇的第一个位点**，因此簇宽 <= window_bp，
          不会出现“链式扩展”把一整片位点连成一簇。window_bp=0 时每个位点各自成簇。
    """
    clusters: List[List[int]] = []
    for position in sorted(positions):
        if clusters and position - clusters[-1][0] <= window_bp:
            clusters[-1].append(position)
        else:
            clusters.append([position])
    return clusters


def choose_representative_site(
    raw_positions: Sequence[int],
    chromosome: str,
    evidence_by_breakpoint: Dict[Tuple[str, int], BreakpointEvidence],
) -> int:
    """在合并簇里选代表位点：支持 read 最多者，并列则取坐标较小者（保证可复现）。

    输入示例: 簇内 chr4:5000000（1 条支持）与 chr4:5000010（8 条支持）
    输出示例: 5000010
    """

    def support_read_count(position: int) -> int:
        bucket = evidence_by_breakpoint.get((chromosome, position))
        if bucket is None:
            return 0
        return len(bucket.split_read_names) + len(bucket.pair_split_read_names)

    return min(raw_positions, key=lambda position: (-support_read_count(position), position))


def merge_nearby_breakpoints(
    evidence_by_breakpoint: Dict[Tuple[str, int], BreakpointEvidence],
    duplicates_by_breakpoint: Dict[Tuple[str, int], Set[str]],
    window_bp: int,
) -> Tuple[
    Dict[Tuple[str, int], BreakpointEvidence],
    Dict[Tuple[str, int], Set[str]],
    Dict[Tuple[str, int], Tuple[int, ...]],
]:
    """把同一染色体上相距 <= window_bp 的断点合并为一个位点（window_bp=0 时保持原样）。

    输入示例:
        evidence_by_breakpoint = {('chr4', 5000000): 1 条支持, ('chr4', 5000010): 8 条支持}
        duplicates_by_breakpoint = {('chr4', 5000010): {'readA', 'readB'}}
        window_bp = 10
    输出示例:
        ({('chr4', 5000010): 合并后的 BreakpointEvidence},
         {('chr4', 5000010): {'readA', 'readB'}},
         {('chr4', 5000010): (5000000, 5000010)})

    说明: read 名单取并集、干扰源头名单取并集、DUP 集合取并集，因此不会重复计数；
          第三个返回值记录“代表位点 <- 哪些原始位点”，用于输出合并审计信息。
    """
    positions_by_chromosome: Dict[str, List[int]] = {}
    for chromosome, position in evidence_by_breakpoint:
        positions_by_chromosome.setdefault(chromosome, []).append(position)

    merged_evidence: Dict[Tuple[str, int], BreakpointEvidence] = {}
    merged_duplicates: Dict[Tuple[str, int], Set[str]] = {}
    merged_positions_by_site: Dict[Tuple[str, int], Tuple[int, ...]] = {}

    for chromosome in sorted(positions_by_chromosome):
        clusters = cluster_positions_within_window(positions_by_chromosome[chromosome], window_bp)
        for cluster in clusters:
            representative = choose_representative_site(cluster, chromosome, evidence_by_breakpoint)
            site = (chromosome, representative)
            accumulator = BreakpointEvidence(
                split_read_names=set(),
                pair_split_read_names=set(),
                owned_query_read_names=set(),
                owned_other_read_names=set(),
                ambiguous_read_names=set(),
                no_insert_hit_read_names=set(),
                consensus_read_names=set(),
                interference_source_read_names={},
            )
            for position in cluster:
                bucket = evidence_by_breakpoint.get((chromosome, position))
                if bucket is None:
                    continue
                for field_name in (
                    "split_read_names",
                    "pair_split_read_names",
                    "owned_query_read_names",
                    "owned_other_read_names",
                    "ambiguous_read_names",
                    "no_insert_hit_read_names",
                    "consensus_read_names",
                ):
                    getattr(accumulator, field_name).update(getattr(bucket, field_name))
                for source, read_names in bucket.interference_source_read_names.items():
                    accumulator.interference_source_read_names.setdefault(source, set()).update(
                        read_names
                    )
            merged_evidence[site] = accumulator
            merged_positions_by_site[site] = tuple(cluster)

            duplicate_read_names: Set[str] = set()
            for position in cluster:
                duplicate_read_names.update(
                    duplicates_by_breakpoint.get((chromosome, position), set())
                )
            if duplicate_read_names:
                merged_duplicates[site] = duplicate_read_names

    return merged_evidence, merged_duplicates, merged_positions_by_site


def filter_by_minimum_support(
    evidence_by_breakpoint: Dict[Tuple[str, int], BreakpointEvidence],
    duplicates_by_breakpoint: Dict[Tuple[str, int], Set[str]],
    merged_positions_by_site: Dict[Tuple[str, int], Tuple[int, ...]],
    minimum_support: int,
) -> Tuple[
    Dict[Tuple[str, int], BreakpointEvidence],
    Dict[Tuple[str, int], Set[str]],
    Dict[Tuple[str, int], Tuple[int, ...]],
    int,
]:
    """按位点支持数（SR + PESR）过滤，丢弃支持不足的位点。

    输入示例:
        evidence = {('chr1', 100): 2 条 SR 的桶, ('chr1', 900): 1 条 SR 的桶}
        minimum_support = 2
    输出示例: ({('chr1', 100): 桶}, {}, {('chr1', 100): (100,)}, 1)   # 末位是被丢弃的位点数

    说明: minimum_support <= 1 时原样返回（默认值，保持历史输出逐字节不变）。
          三份结构同步过滤，保证四张表口径一致。
    """
    if minimum_support <= 1:
        return (
            evidence_by_breakpoint,
            duplicates_by_breakpoint,
            merged_positions_by_site,
            0,
        )
    kept = {
        key
        for key, bucket in evidence_by_breakpoint.items()
        if len(bucket.split_read_names) + len(bucket.pair_split_read_names) >= minimum_support
    }
    dropped = len(evidence_by_breakpoint) - len(kept)
    return (
        {key: bucket for key, bucket in evidence_by_breakpoint.items() if key in kept},
        {key: names for key, names in duplicates_by_breakpoint.items() if key in kept},
        {key: positions for key, positions in merged_positions_by_site.items() if key in kept},
        dropped,
    )
