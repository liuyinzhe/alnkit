"""三张结果表的写出（列名/列序/编码与历史产物一致）。

已知历史缺陷（**刻意保留**，改了就破坏与旧产物的逐字节对账）：
``reads_info`` 表的表头是 11 列，最后一项写成 ``"stand,mapping_quality"``（一个字段名里
塞了两列），所以数据行实际是 12 列。按**列序**解析，不要按表头名取列。
`short scan` 的 README/文档里必须提醒这一点。
"""

from __future__ import annotations

from pathlib import Path
from typing import Dict, List, Sequence, Set, Tuple

from .model import BreakpointEvidence, ReadEvidence, RuntimeConfig


# ------------------------------------------------------------------------------
# ### 10. 结果输出（表头与列顺序保持历史格式）###
# ------------------------------------------------------------------------------
READS_INFO_HEADER = "\t".join(
    [
        "read_name",
        "read_pair",
        "cigarstring",
        "is_supplementary",
        "chrom",
        "M_ref_pos_start",
        "M_ref_pos_end",
        "M_read_pos_start",
        "M_read_pos_end",
        "M_var_len",
        "stand,mapping_quality",
    ]
)


INSERT_POSITION_HEADER = "\t".join(
    [
        "Chromosome",
        "Insert_Position",
        "split_read_support",
        "pair_split_read_support",
        "dup_read_count",
        "split_read_name",
        "pair_split_read_name",
        "dup_read_name",
    ]
)


# 新增的第三张表（可选的诊断表，**不影响**上面两张历史表的结构）：
# 逐断点标出“是否同源干扰 / 干扰源头 / 各源头支持 read 数”。
INSERT_INTERFERENCE_HEADER = "\t".join(
    [
        "Chromosome",
        "Insert_Position",
        "merged_position_count",
        "merged_positions",
        "query_insert",
        "is_homology_interference",
        "interference_source_count",
        "interference_source_support_reads",
        "support_reads_total",
        "support_reads_owned_query",
        "support_reads_owned_other",
        "support_reads_consensus_assigned",
        "support_reads_ambiguous",
        "support_reads_no_insert_hit",
        "ambiguous_read_names",
    ]
)


def write_reads_info_table(output_path: Path, read_evidences: Sequence[ReadEvidence]) -> int:
    """写出 <sample>.reads_info.out.xls（保留比对的明细）。

    输入示例:
        read_evidences = [ReadEvidence(read_name='LH...:1350:9163:17577',
                                       all_records=[record(read_pair='read1', chromosome='chr1',
                                                           block_start=6000050, block_end=6000000, ...)])]
    输出示例（TAB 分隔，表头为历史保留的 11 列写法）:
        read_name\tread_pair\tcigarstring\tis_supplementary\tchrom\tM_ref_pos_start\t...
        LH...:1350:9163:17577\tread1\t68M82S\tFalse\tchr1\t6000050\t6000000\t0\t68\t68\t-\t60

    返回: 写出的数据行数
    """
    rows: List[List[str]] = []
    for read_evidence in read_evidences:
        for record in read_evidence.all_records:
            rows.append(
                [
                    record.read_name,
                    record.read_pair,
                    record.cigarstring,
                    str(record.is_supplementary),
                    record.chromosome,
                    str(record.block_start),
                    str(record.block_end),
                    str(record.read_start),
                    str(record.read_end),
                    str(record.block_length),
                    record.strand,
                    str(record.mapping_quality),
                ]
            )
    # 稳定排序：让输出可复现（历史实现依赖集合迭代顺序，同一输入两次运行结果不同）
    rows.sort(key=lambda row: (row[0], row[1], row[4], int(row[5]), row[2]))

    with output_path.open("w", encoding="utf-8") as output_handle:
        output_handle.write(READS_INFO_HEADER + "\n")
        for row in rows:
            output_handle.write("\t".join(row) + "\n")
    return len(rows)


def write_insert_position_table(
    output_path: Path,
    evidence_by_breakpoint: Dict[Tuple[str, int], BreakpointEvidence],
    duplicates_by_breakpoint: Dict[Tuple[str, int], Set[str]],
) -> int:
    """写出 <sample>.Insert_Position.result.xls（断点证据汇总，表头与列顺序保持历史格式）。

    输入示例:
        evidence_by_breakpoint = {('chr1', 6000000): BreakpointEvidence(
            split_read_names={'LH.../read1'}, pair_split_read_names=set(), ...)}
        duplicates_by_breakpoint = {('chr1', 6000000): {'LH...:1350:9163:17577'}}
    输出示例（TAB 分隔）:
        Chromosome\tInsert_Position\tsplit_read_support\tpair_split_read_support\tdup_read_count\t...
        chr1\t6000000\t1\t0\t1\tLH.../read1\t\tLH...:1350:9163:17577

    返回: 写出的数据行数
    """
    rows: List[List[str]] = []
    for breakpoint_key in sorted(evidence_by_breakpoint, key=lambda key: (key[0], key[1])):
        chromosome, junction_position = breakpoint_key
        bucket = evidence_by_breakpoint[breakpoint_key]
        split_read_names = sorted(bucket.split_read_names)
        pair_split_read_names = sorted(bucket.pair_split_read_names)
        duplicate_read_names = sorted(duplicates_by_breakpoint.get(breakpoint_key, set()))
        rows.append(
            [
                chromosome,
                str(junction_position),
                str(len(split_read_names)),
                str(len(pair_split_read_names)),
                str(len(duplicate_read_names)),
                ";".join(split_read_names),
                ";".join(pair_split_read_names),
                ";".join(duplicate_read_names),
            ]
        )

    with output_path.open("w", encoding="utf-8") as output_handle:
        output_handle.write(INSERT_POSITION_HEADER + "\n")
        for row in rows:
            output_handle.write("\t".join(row) + "\n")
    return len(rows)


def write_insert_interference_table(
    output_path: Path,
    evidence_by_breakpoint: Dict[Tuple[str, int], BreakpointEvidence],
    merged_positions_by_site: Dict[Tuple[str, int], Tuple[int, ...]],
    config: RuntimeConfig,
) -> int:
    """写出 <sample>.Insert_Interference.result.xls（逐断点的同源干扰标记与源头计数）。

    输入示例:
        断点 ('chr2', 7000000)：支持 read 3 条，其中 1 条归属 query、2 条无法区分；
        2 条无法区分的 read 同时命中 target_ref 与 homolog_ref
    输出示例（TAB 分隔，共 15 列）:
        Chromosome  Insert_Position  merged_position_count  merged_positions  query_insert
        is_homology_interference  interference_source_count  interference_source_support_reads
        support_reads_total  support_reads_owned_query  support_reads_owned_other
        support_reads_consensus_assigned  support_reads_ambiguous  support_reads_no_insert_hit
        ambiguous_read_names
        chr2  7000000  2  7000000;7000000  target_ref  True  2  target_ref:2;homolog_ref:2
        3  1  0  0  2  0  readB;readC

    说明:
        - merged_position_count / merged_positions：该行由哪几个原始位点合并而来
          （--merge-window-bp 0 时就是它自己）
        - is_homology_interference：该断点是否存在同源干扰（支持 read 里有同源命中或无法区分的）
        - interference_source_count：干扰源头个数
        - interference_source_support_reads：每个源头有多少条支持 read，格式 源头:条数;源头:条数
          （按条数降序、同数按名字排序，保证输出可复现）
        - 各支持 read 计数列与过滤策略无关；legacy 策略下 owned_other / no_insert_hit 会非 0，
          consensus 策略下 consensus_assigned 会非 0

    返回: 写出的数据行数
    """
    rows: List[List[str]] = []
    for breakpoint_key in sorted(evidence_by_breakpoint, key=lambda key: (key[0], key[1])):
        chromosome, junction_position = breakpoint_key
        bucket = evidence_by_breakpoint[breakpoint_key]
        source_reads = bucket.interference_source_read_names
        # 按（支持 read 数降序, 源头名）输出，保证字节级可复现
        source_text = ";".join(
            f"{source}:{len(read_names)}"
            for source, read_names in sorted(
                source_reads.items(), key=lambda item: (-len(item[1]), item[0])
            )
        )
        owned_query_names = sorted(bucket.owned_query_read_names)
        owned_other_names = sorted(bucket.owned_other_read_names)
        consensus_names = sorted(bucket.consensus_read_names)
        ambiguous_names = sorted(bucket.ambiguous_read_names)
        no_insert_hit_names = sorted(bucket.no_insert_hit_read_names)
        merged_positions = merged_positions_by_site.get(breakpoint_key, (junction_position,))
        rows.append(
            [
                chromosome,
                str(junction_position),
                str(len(merged_positions)),
                ";".join(str(position) for position in merged_positions),
                config.query_chromosome,
                str(bool(source_reads)),
                str(len(source_reads)),
                source_text,
                str(
                    len(owned_query_names)
                    + len(owned_other_names)
                    + len(consensus_names)
                    + len(ambiguous_names)
                    + len(no_insert_hit_names)
                ),
                str(len(owned_query_names)),
                str(len(owned_other_names)),
                str(len(consensus_names)),
                str(len(ambiguous_names)),
                str(len(no_insert_hit_names)),
                ";".join(ambiguous_names),
            ]
        )

    with output_path.open("w", encoding="utf-8") as output_handle:
        output_handle.write(INSERT_INTERFERENCE_HEADER + "\n")
        for row in rows:
            output_handle.write("\t".join(row) + "\n")
    return len(rows)


def _split_other_inserts(
    source_read_names: Dict[str, Set[str]], query_insert: str
) -> Tuple[str, List[str]]:
    """把一个位点的命中源头拆成（最近的其它插入片段, 其余其它插入片段）。

    一条 read 可能同时命中多个插入序列（同源/共有 read），所以一个位点会有多个源头。
    * **query 自己不算竞争者**——它出现只说明"这条 read 也命中了本次的 query"；
    * “最近”按**支持该插入片段的 read 数**定义（同源有多近没法用坐标衡量，用证据量更实在），
      并列时取名字小的，保证输出可复现。
    """
    others = {name: reads for name, reads in source_read_names.items() if name != query_insert}
    if not others:
        return "-", []
    ranked = sorted(others.items(), key=lambda item: (-len(item[1]), item[0]))
    return ranked[0][0], [name for name, _ in ranked[1:]]


def write_canonical_sites_table(
    output_path: Path,
    evidence_by_breakpoint: Dict[Tuple[str, int], BreakpointEvidence],
    duplicates_by_breakpoint: Dict[Tuple[str, int], Set[str]],
    sample_name: str,
    query_insert: str,
) -> int:
    """写出 <sample>.sites.canonical.tsv —— 二代位点的**机器接口**（新增，不替代 legacy 三表）。

    历史表里 ``Insert_Position`` 是"断点左侧碱基的 1-based 坐标"，但数值上又恰好等于该断点的
    0-based 边界，光看数字分不出来。这里把三种基准一次写清，恒等式由测试守住：

        insert_pos_1based - 1 == insert_pos_0based
        insert_boundary_0based == insert_pos_1based

    行与 ``<sample>.Insert_Position.result.xls`` 一一对应（同样的位点、同样的合并结果），
    可用 ``(chrom, insert_pos_1based)`` 与原表 join。

    插入归属三列：``query_insert`` 是本次运行的 query；``nearest_other_insert`` 是**除它之外**
    支持 read 最多的插入片段（一条 read 可能同时命中多个插入序列，所以"最近"按 read 支持数
    定义，而不是按坐标距离；query 自己不算竞争者）；``other_inserts`` 是其余被命中过的插入片段。
    """
    from ..core import tables

    rows = []
    for chromosome, junction_position in sorted(
        evidence_by_breakpoint, key=lambda key: (key[0], key[1])
    ):
        bucket = evidence_by_breakpoint[(chromosome, junction_position)]
        nearest, others = _split_other_inserts(bucket.interference_source_read_names, query_insert)
        rows.append(
            {
                "chrom": chromosome,
                # junction_position 是 0-based 边界；左侧碱基的 1-based 坐标在数值上等于它
                "insert_pos_1based": junction_position,
                "insert_pos_0based": junction_position - 1,
                "insert_boundary_0based": junction_position,
                "n_sr": len(bucket.split_read_names),
                "n_pesr": len(bucket.pair_split_read_names),
                "n_dup": len(duplicates_by_breakpoint.get((chromosome, junction_position), set())),
                "sample": sample_name,
                "query_insert": query_insert,
                "nearest_other_insert": nearest,
                "other_inserts": ",".join(others) if others else "-",
            }
        )
    return tables.write_tsv(
        output_path,
        "insert_sites",
        rows,
        comments=[
            "二代位点规范表（机器接口）；与 <sample>.Insert_Position.result.xls 按 "
            "(chrom, insert_pos_1based) 一一对应。",
            "insert_pos_1based = 断点左侧碱基的 1-based 坐标；insert_boundary_0based = 断点边界。",
        ],
        sort_key=lambda row: (row["chrom"], row["insert_pos_1based"]),
    )
