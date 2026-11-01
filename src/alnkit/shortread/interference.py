"""同源干扰：哪些别的插入序列污染了当前位点的证据。"""

from __future__ import annotations

from typing import Sequence, Set

from .model import InsertHit


def collect_homology_hit_chromosomes(insert_hits: Sequence[InsertHit]) -> Set[str]:
    """从插入命中里挑出“同源/共有 read”信号：来自 secondary 比对或 XA 标签的命中。

    输入示例:
        insert_hits = [InsertHit('target_ref', 37, 150, 0, 0, 'alignment', record),
                       InsertHit('homolog_ref', 0, 31, 0, 119, 'xa', None)]
    输出示例: {'homolog_ref'}

    说明: 这是“该 read 与别的序列共有”的唯一判据；只落在 query 上的 primary/supplementary
          锚定（source='alignment'）不算同源干扰。
    """
    return {hit.chromosome for hit in insert_hits if hit.source != "alignment"}


def determine_interference_sources(
    origin_label: str, candidate_insert_chromosomes: Set[str], homology_hit_chromosomes: Set[str]
) -> Set[str]:
    """汇总一条 read 的“同源干扰源头”。

    输入示例（三种典型情况）:
        1) origin_label='owned_query',   candidate={'target_ref'}, homology={'homolog_ref'}
        2) origin_label='ambiguous',      candidate={'target_ref', 'homolog_ref', 'decoy_ref'}, homology=同集
        3) origin_label='owned_query',   candidate={'target_ref'}, homology=set()
    输出示例:
        1) {'homolog_ref'}                       # 明确归属 target，但另与同源序列一致
        2) {'target_ref', 'homolog_ref', 'decoy_ref'}   # 无法判定，全部列为干扰源头
        3) set()                               # 干净的唯一插入证据

    说明: 这是“是否同源干扰/干扰来自哪里”的唯一判据；表格里的干扰源头计数由它汇总而来。
    """
    if homology_hit_chromosomes:
        return set(homology_hit_chromosomes)
    if origin_label == "ambiguous":
        return set(candidate_insert_chromosomes)
    return set()
