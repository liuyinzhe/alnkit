#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# A3：把 reads 分成「载体区域 / 接头区域」，并把接头事件聚成位点与区域。

# 来源：upstream/scripts/03_classify_reads.py（700 行）的忠实移植。逻辑未改；
# 相对原脚本只有三处包装：参数解析搬到 alnkit.cli.long、库层收 ClassifySettings、
# main() 改为 run()。输出表 / 子集 FASTQ 与旧脚本逐字节一致。
#
# 输入：
#   02_align_merged/reads_to_merged.sorted.bam   reads -> 合并参考（参考 + 外源）
#   01_merge_ref/merged_reference.fa.fai          取序列长度
#   data/rawdata/input.fastq                      原始 reads（输出子集用）
# 输出：
#   03_classify_reads/alignment_segments.tsv      只列有外源证据的 read 的非 secondary 比对段
#   03_classify_reads/read_classification.tsv     只列有外源证据的 read（junction / vector_region）
#   03_classify_reads/classification_summary.tsv   总数/目标数/非目标数（不输出非目标逐条记录）
#   以下旧 200 kb 中间表仅脚本内部兼容输出，shell 随后删除，不参与方案 A：
#   03_classify_reads/junction_events.tsv / junction_sites.tsv / locus_regions.tsv
#   03_classify_reads/related_reads.fastq         所有相关 reads（载体区域 + 接头）
#   03_classify_reads/vector_region_reads.fastq   只比对外源、连不上基因组的 reads
#   03_classify_reads/junction_reads.fastq        同时比对外源与基因组的 reads
#
# 关键口径：
#   * 外源证据：外源比对段 >= --min-exog-block-bp（默认 100 bp）且 MAPQ >= --min-mapq（20）。
#   * 基因组证据：参考基因组比对段 MAPQ 达标；若同时满足
#       “落在 contig 两端 --contig-edge-bp（100 bp）内”且“长度 <= --short-genome-block-bp（500 bp）”
#     则判为可疑错配，丢弃。
#   * 接头事件：基因组段与外源段在 read 轴上相邻（间隔 <= --junction-max-gap-bp，默认 10 kb；
#     实测接头两侧常有几 kb 未比对上的间隔，所以这里不要求严格紧邻）。
#     边界和 side 都按**参考基因组方向**给出（已考虑比对链）：
#       genome_left  = 基因组在边界左侧，插入从边界开始；
#       genome_right = 基因组在边界右侧，插入在边界结束。
#       （07 画丝带图时取边界旁 50 bp 作锚点。）
#   * 位点：同染色体、同 side 的边界差 <= --site-gap-bp（20 bp）聚成一个位点。
#   * 区域/类别：同染色体上位点边界相距 <= --region-gap-bp（200 kb）的位点归为一个区域；
#     区域是后续“按接头位点分类”的类别单位（一个插入可能有多个断点/两侧边界）。
#   * 内存：先按外源序列找出目标 read（通常远少于全部 reads），再只加载这些 read 的
#     非 secondary 比对段；genome_only/other 不进入本表，也不加载它们的比对段。
#
# 输入示例：
#   --bam 02_align_merged/reads_to_merged.sorted.bam
#   --reads data/rawdata/input.fastq --merged-fai 01_merge_ref/merged_reference.fa.fai
#   --exog insertA --exog insertB --outdir 03_classify_reads
# 输出示例（read_classification.tsv）：
#   read_name  read_length  read_category  n_exog_segments  exog_names  ...
#   sampleA-run/ccs  21569  junction  1  insertB ...
# ---------------------------------------------------------------------------
from __future__ import annotations

import collections
import os
from dataclasses import dataclass

import pysam

from ..core.bam_segments import load_alignment_segments as _load_alignment_segments
from ..core.bam_segments import load_fai_lengths
from ..core.coords import merge_intervals
from ..core.errors import UsageError
from ..core.tables import write_legacy_tsv as write_tsv

# 旧脚本的比对段是 dict；core 用 Segment 模型。这里做一层边界适配，
# 字段名与空值口径（chrom_length 缺失记 "-"）保持与 alignment_table 完全一致。
_SEGMENT_FIELDS = (
    ("chrom", "chrom"),
    ("read_name", "read_name"),
    ("ref_start", "ref_start_0based"),
    ("ref_end", "ref_end_0based"),
    ("read_start", "read_start_0based"),
    ("read_end", "read_end_0based"),
    ("read_length", "read_length"),
    ("strand", "strand"),
    ("mapq", "mapq"),
    ("nm", "nm"),
    ("match_columns", "match_columns"),
    ("inserted_bases", "inserted_bases"),
    ("deleted_bases", "deleted_bases"),
    ("ref_aln_len", "ref_aln_len"),
    ("read_aln_len", "read_aln_len"),
    ("cigar", "cigar"),
)


def _segment_to_legacy_dict(segment):
    record = {legacy: getattr(segment, field) for legacy, field in _SEGMENT_FIELDS}
    record["chrom_length"] = "-" if segment.chrom_length is None else segment.chrom_length
    record["is_secondary"] = int(segment.is_secondary)
    record["is_supplementary"] = int(segment.is_supplementary)
    return record


def load_alignment_segments(
    bam_path, chrom_lengths=None, selected_read_names=None, skip_secondary=False
):
    """与旧 ``alignment_table.load_alignment_segments`` 等价的 dict 版加载。"""
    return [
        _segment_to_legacy_dict(segment)
        for segment in _load_alignment_segments(
            bam_path,
            chrom_lengths=chrom_lengths,
            selected_read_names=selected_read_names,
            skip_secondary=skip_secondary,
        )
    ]


### 读取比对段 ###


def find_target_read_names(bam_path, exog_names, min_mapq, min_exog_block_bp):
    """先只找出“有合格外源比对”的 read 名（不把全部比对段读进内存）。

    输入示例: find_target_read_names("reads_to_merged.sorted.bam", ("insertA", "insertB"), 20, 100)
    输出示例: {"m.../ccs", "m.../ccs2", ...}
    说明：
      * 优先用 BAM 索引按外源序列 `fetch`，只扫描外源上的记录（很小）；
      * 外源名不在 BAM 里、或 BAM 没有索引时，退化成一次全量扫描；
      * 过滤条件与 classify_segment 的外源证据一致：非 secondary、MAPQ 达标、
        外源比对跨度 >= min_exog_block_bp。
    """
    if not os.path.exists(bam_path):
        raise UsageError(f"找不到 BAM：{bam_path}")

    target_read_names = set()

    def keep_exog_record(read):
        if read.is_unmapped or read.is_secondary:
            return False
        if read.reference_name not in exog_names:
            return False
        if read.mapping_quality < min_mapq:
            return False
        if (read.reference_end - read.reference_start) < min_exog_block_bp:
            return False
        return True

    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            available_references = set(bam_file.references)
            missing_exog = [name for name in exog_names if name not in available_references]
            if missing_exog:
                raise UsageError(
                    "--exog 里这些名字不在 BAM header 中："
                    + ", ".join(missing_exog)
                    + "；可用序列名示例："
                    + ", ".join(sorted(available_references)[:5])
                )
            for exog_name in exog_names:
                if exog_name not in available_references:
                    continue
                for read in bam_file.fetch(exog_name):
                    if keep_exog_record(read):
                        target_read_names.add(read.query_name)
    except (ValueError, OSError):
        # 没有 .bai 或 fetch 失败：退化成一次全量扫描，只保留目标 read 名（仍然不存比对段）
        with pysam.AlignmentFile(bam_path, "rb") as bam_file:
            for read in bam_file:
                if keep_exog_record(read):
                    target_read_names.add(read.query_name)
    return target_read_names


def load_segments_by_read(bam_path, chrom_lengths, selected_read_names):
    """读取 BAM，只保留指定 read 的非 secondary 比对段。

    输入示例: load_segments_by_read("reads_to_merged.sorted.bam", {"insertA": 1353, ...},
                                    {"m.../ccs", "m.../ccs2"})
    输出示例: {"m.../ccs": [{"chrom": "insertA", "ref_start": 0, "ref_end": 1353, ...}, ...]}
    说明：supplementary 保留（split alignment 是接头证据）；secondary 在解析 CIGAR 前就丢弃；
          selected_read_names 由 find_target_read_names() 得到，避免加载非目标 reads。
    """
    segments_by_read = collections.defaultdict(list)
    for segment in load_alignment_segments(
        bam_path,
        chrom_lengths=chrom_lengths,
        selected_read_names=selected_read_names,
        skip_secondary=True,
    ):
        segments_by_read[segment["read_name"]].append(segment)
    return segments_by_read


def classify_segment(segment, exog_names, settings):
    """判断一个比对段是“外源证据 / 基因组证据 / 可疑错配”，返回 (role, kept, reason)。

    输入示例: classify_segment({"chrom": "insertA", "ref_aln_len": 1353, "mapq": 60, ...},
                              ("insertA", "insertB"), args)
    输出示例: ("exog", True, "-")
    说明：
      * role = exog / genome / other；
      * 外源段太短或 MAPQ 太低 -> kept=False；
      * 基因组段落在 contig 两端 100 bp 内且 <= 500 bp -> kept=False（可疑错配）。
    """
    chrom = segment["chrom"]
    if chrom in exog_names:
        if segment["ref_aln_len"] < settings.min_exog_block_bp:
            return "exog", False, "exog_block_too_short"
        if segment["mapq"] < settings.min_mapq:
            return "exog", False, "low_mapq"
        return "exog", True, "-"

    if segment["mapq"] < settings.min_mapq:
        return "genome", False, "low_mapq"
    chrom_length = segment["chrom_length"]
    if isinstance(chrom_length, int):
        near_contig_edge = (
            segment["ref_start"] < settings.contig_edge_bp
            or chrom_length - segment["ref_end"] < settings.contig_edge_bp
        )
        if near_contig_edge and segment["ref_aln_len"] <= settings.short_genome_block_bp:
            return "genome", False, "contig_edge_short_block"
    return "genome", True, "-"


### 接头事件与位点聚类 ###


def find_junction_events(read_name, exog_segments, genome_segments, max_gap_bp):
    """找出“基因组段紧邻外源段”的接头事件（每个基因组段左右两侧各最多一个）。

    输入示例：
      exog_segments   = [{"chrom": "insertA", "read_start": 25000, "read_end": 26353, ...}]
      genome_segments = [{"chrom": "chr1", "read_start": 0, "read_end": 25000,
                          "ref_start": 2002000, "ref_end": 2004500, ...}]
      max_gap_bp = 100
    输出示例：
      [{"read_name": ..., "chrom": "chr1", "genome_side": "genome_left",
        "boundary": 2004500, "vector_name": "insertA", "gap_bp": 0, ...}]
    说明：
      * 先看外源段在 read 轴的哪一侧，再结合基因组段的比对链，换算成参考基因组方向：
          - 外源在 read 右侧、基因组链 + ：genome_left，边界 = ref_end
          - 外源在 read 右侧、基因组链 - ：genome_right，边界 = ref_start
          - 外源在 read 左侧、基因组链 + ：genome_right，边界 = ref_start
          - 外源在 read 左侧、基因组链 - ：genome_left，边界 = ref_end
      * boundary 是 0-based 边界；genome_left 时等于最后一个基因组碱基的 1-based 坐标。
    """
    events = []
    for genome_segment in genome_segments:
        # 右侧最近的外源段（外源段起点 >= 基因组段终点）
        right_vectors = [
            segment
            for segment in exog_segments
            if segment["read_start"] >= genome_segment["read_end"]
        ]
        if right_vectors:
            vector_segment = min(right_vectors, key=lambda segment: segment["read_start"])
            gap_bp = vector_segment["read_start"] - genome_segment["read_end"]
            if gap_bp <= max_gap_bp:
                # 外源在 read 右侧：链 + 时基因组在边界左侧，链 - 时在右侧
                if genome_segment["strand"] == "+":
                    genome_side, boundary = "genome_left", genome_segment["ref_end"]
                else:
                    genome_side, boundary = "genome_right", genome_segment["ref_start"]
                events.append(
                    {
                        "read_name": read_name,
                        "chrom": genome_segment["chrom"],
                        "genome_side": genome_side,
                        "boundary": boundary,
                        "vector_name": vector_segment["chrom"],
                        "gap_bp": gap_bp,
                        "genome_read_start": genome_segment["read_start"],
                        "genome_read_end": genome_segment["read_end"],
                        "vector_read_start": vector_segment["read_start"],
                        "vector_read_end": vector_segment["read_end"],
                        "genome_mapq": genome_segment["mapq"],
                        "vector_mapq": vector_segment["mapq"],
                        "genome_identity": estimate_identity(genome_segment),
                        "vector_identity": estimate_identity(vector_segment),
                    }
                )

        # 左侧最近的外源段（外源段终点 <= 基因组段起点）
        left_vectors = [
            segment
            for segment in exog_segments
            if segment["read_end"] <= genome_segment["read_start"]
        ]
        if left_vectors:
            vector_segment = max(left_vectors, key=lambda segment: segment["read_end"])
            gap_bp = genome_segment["read_start"] - vector_segment["read_end"]
            if gap_bp <= max_gap_bp:
                # 外源在 read 左侧：链 + 时基因组在边界右侧，链 - 时在左侧
                if genome_segment["strand"] == "+":
                    genome_side, boundary = "genome_right", genome_segment["ref_start"]
                else:
                    genome_side, boundary = "genome_left", genome_segment["ref_end"]
                events.append(
                    {
                        "read_name": read_name,
                        "chrom": genome_segment["chrom"],
                        "genome_side": genome_side,
                        "boundary": boundary,
                        "vector_name": vector_segment["chrom"],
                        "gap_bp": gap_bp,
                        "genome_read_start": genome_segment["read_start"],
                        "genome_read_end": genome_segment["read_end"],
                        "vector_read_start": vector_segment["read_start"],
                        "vector_read_end": vector_segment["read_end"],
                        "genome_mapq": genome_segment["mapq"],
                        "vector_mapq": vector_segment["mapq"],
                        "genome_identity": estimate_identity(genome_segment),
                        "vector_identity": estimate_identity(vector_segment),
                    }
                )
    return events


def estimate_identity(segment):
    """用简化 CIGAR 的 M 长度与 NM 估算一致率（没有 NM 时按 1.0）。

    输入示例: segment={"match_columns": 1353, "nm": 3, "inserted_bases": 0, "deleted_bases": 0}
    输出示例: 0.9978
    """
    match_columns = segment.get("match_columns") or 0
    if match_columns <= 0:
        return 1.0
    edits = (
        (segment.get("nm") or 0)
        - (segment.get("inserted_bases") or 0)
        - (segment.get("deleted_bases") or 0)
    )
    return max(0.0, 1.0 - max(0, edits) / match_columns)


def cluster_events_into_sites(events, site_gap_bp):
    """按 (染色体, genome_side) 聚位点：边界差 <= site_gap_bp 的事件归为同一位点。

    输入示例: events=[{"chrom": "chr4", "genome_side": "genome_right",
                       "boundary": 3000000}, ...]
    输出示例: [{"chrom": "chr4", "genome_side": "genome_right",
                "boundary_min": 2999995, "boundary_max": 3000005,
                "boundary_median": 3000000, "events": [...]}]
    说明：单链接聚类（与上一个事件的边界差 <= 阈值就并入），不做两两全比较。
    """
    grouped = collections.defaultdict(list)
    for event in events:
        grouped[(event["chrom"], event["genome_side"])].append(event)

    sites = []
    for (chrom, genome_side), site_events in grouped.items():
        site_events.sort(key=lambda event: event["boundary"])
        current = []
        for event in site_events:
            if current and event["boundary"] - current[-1]["boundary"] > site_gap_bp:
                sites.append(build_site(current))
                current = []
            current.append(event)
        if current:
            sites.append(build_site(current))

    sites.sort(key=lambda site: (site["chrom"], site["boundary_median"], site["genome_side"]))
    for index, site in enumerate(sites, 1):
        site["site_id"] = "S{0:02d}".format(index)
    return sites


def build_site(site_events):
    """把一个边界簇整理成位点记录（边界 min/max/中位数、reads、外源名）。

    输入示例: site_events=[{... "boundary": 3000000, "vector_name": "insertA"}, ...]
    输出示例: {"chrom": "chr1", "genome_side": "genome_right",
                "boundary_min": 3000000, "boundary_max": 3000000,
                "boundary_median": 3000000, "events": [...],
                "read_names": {...}, "vector_names": {...}}
    """
    boundaries = sorted(event["boundary"] for event in site_events)
    return {
        "chrom": site_events[0]["chrom"],
        "genome_side": site_events[0]["genome_side"],
        "boundary_min": boundaries[0],
        "boundary_max": boundaries[-1],
        "boundary_median": boundaries[len(boundaries) // 2],
        "events": site_events,
        "read_names": {event["read_name"] for event in site_events},
        "vector_names": {event["vector_name"] for event in site_events},
    }


def group_sites_into_regions(sites, region_gap_bp):
    """把同染色体上相距 <= region_gap_bp 的位点归为一个区域（类别单位）。

    输入示例: sites=[{"chrom": "chr4", "boundary_min": 3000000, ...},
                     {"chrom": "chr4", "boundary_min": 3000500, ...}]
    输出示例: [{"region_id": "R01", "chrom": "chr4", "start": 3000000,
                "end": 3000500, "sites": [...], "read_names": {...}}]
    说明：左侧和右侧边界也会被归进同一区域（一个插入的两端）；区域只是分类单位，
         写表时仍保留每个位点的 genome_side 与精确边界。
    """
    grouped = collections.defaultdict(list)
    for site in sites:
        grouped[site["chrom"]].append(site)

    regions = []
    for chrom, chrom_sites in grouped.items():
        chrom_sites.sort(key=lambda site: site["boundary_min"])
        current = []
        for site in chrom_sites:
            if (
                current
                and site["boundary_min"] - max(item["boundary_max"] for item in current)
                > region_gap_bp
            ):
                regions.append(build_region(chrom, current))
                current = []
            current.append(site)
        if current:
            regions.append(build_region(chrom, current))

    regions.sort(key=lambda region: (region["chrom"], region["start"]))
    for index, region in enumerate(regions, 1):
        region["region_id"] = "R{0:02d}".format(index)
        region["region_label"] = "{0}|{1}:{2}-{3}".format(
            region["region_id"], region["chrom"], region["start"], region["end"]
        )
        for site in region["sites"]:
            site["region_id"] = region["region_id"]
    return regions


def build_region(chrom, region_sites):
    """把一组位点整理成区域记录。

    输入示例: region_sites=[site1, site2]
    输出示例: {"chrom": chrom, "start": 3000000, "end": 3000500, "sites": [...],
                "read_names": {...}, "vector_names": {...}}
    """
    read_names = set()
    vector_names = set()
    for site in region_sites:
        read_names |= site["read_names"]
        vector_names |= site["vector_names"]
    return {
        "chrom": chrom,
        "start": min(site["boundary_min"] for site in region_sites),
        "end": max(site["boundary_max"] for site in region_sites),
        "sites": region_sites,
        "read_names": read_names,
        "vector_names": vector_names,
    }


### 输出 ###


def write_read_classification(rows, output_path):
    """写有外源证据的 read 的分类表（一行一条 read）。

    输入示例: rows=[{"read_name": "m.../ccs", "read_category": "junction", ...}]
    输出示例: read_name  read_length  read_category  n_exog_segments  ...
    """
    columns = [
        "read_name",
        "read_length",
        "read_category",
        "n_exog_segments",
        "exog_names",
        "exog_aln_bp",
        "n_genome_segments",
        "genome_aln_bp",
        "genome_chroms",
        "n_genome_segments_dropped",
        "genome_aln_bp_dropped",
        "max_mapq",
        "notes",
    ]
    table_rows = []
    for row in rows:
        table_rows.append(
            [
                row["read_name"],
                row["read_length"],
                row["read_category"],
                row["n_exog_segments"],
                ",".join(row["exog_names"]) or "-",
                row["exog_aln_bp"],
                row["n_genome_segments"],
                row["genome_aln_bp"],
                ",".join(row["genome_chroms"]) or "-",
                row["n_genome_segments_dropped"],
                row["genome_aln_bp_dropped"],
                row["max_mapq"],
                ";".join(row["notes"]) or "-",
            ]
        )
    write_tsv(
        output_path,
        columns,
        table_rows,
        [
            "只列有外源证据的 read；类别：junction = 同时有外源与基因组证据；"
            "vector_region = 只有外源证据。genome_only/other 不进入本表（不加载非目标 reads）。",
            "坐标：边界为 0-based；genome_aln_bp 是 read 轴上基因组比对区间的并集。",
        ],
    )


def write_junction_events(events, output_path):
    """写接头事件表（每个基因组段紧邻外源段的事件一行）。

    输入示例: events=[{"read_name": "m.../ccs", "chrom": "chr1", ...}]
    输出示例: read_name  chrom  side  boundary  region_id  site_id  vector_name  gap_bp ...
    """
    columns = [
        "read_name",
        "chrom",
        "genome_side",
        "boundary",
        "region_id",
        "site_id",
        "vector_name",
        "gap_bp",
        "genome_read_start",
        "genome_read_end",
        "vector_read_start",
        "vector_read_end",
        "genome_mapq",
        "vector_mapq",
        "genome_identity",
        "vector_identity",
    ]
    rows = []
    for event in events:
        rows.append(
            [
                event["read_name"],
                event["chrom"],
                event["genome_side"],
                event["boundary"],
                event["region_id"],
                event["site_id"],
                event["vector_name"],
                event["gap_bp"],
                event["genome_read_start"],
                event["genome_read_end"],
                event["vector_read_start"],
                event["vector_read_end"],
                event["genome_mapq"],
                event["vector_mapq"],
                "{0:.4f}".format(event["genome_identity"]),
                "{0:.4f}".format(event["vector_identity"]),
            ]
        )
    write_tsv(
        output_path,
        columns,
        rows,
        [
            "边界 boundary 为 0-based；genome_side=genome_left 时等于最后一个基因组碱基的 1-based 坐标，"
            "genome_side=genome_right 时等于插入后第一个基因组碱基的 0-based 坐标。",
            "gap_bp：基因组段与外源段在 read 轴上的间隔（0 = 紧邻/重叠）。",
        ],
    )


def write_junction_sites(sites, output_path):
    """写接头位点表（同 side、边界差 <= 20 bp 的事件聚成一行）。

    输入示例: sites=[{"site_id": "S01", "region_id": "R01", ...}]
    输出示例: site_id  region_id  chrom  side  boundary_min  boundary_max  boundary_median ...
    """
    columns = [
        "site_id",
        "region_id",
        "chrom",
        "genome_side",
        "boundary_min",
        "boundary_max",
        "boundary_median",
        "n_events",
        "n_reads",
        "exog_names",
        "read_names",
    ]
    rows = []
    for site in sites:
        rows.append(
            [
                site["site_id"],
                site["region_id"],
                site["chrom"],
                site["genome_side"],
                site["boundary_min"],
                site["boundary_max"],
                site["boundary_median"],
                len(site["events"]),
                len(site["read_names"]),
                ",".join(sorted(site["vector_names"])),
                ",".join(sorted(site["read_names"])),
            ]
        )
    write_tsv(
        output_path,
        columns,
        rows,
        [
            "位点 = 同染色体、同 genome_side、边界差 <= --site-gap-bp（默认 20 bp）的接头事件。",
            "boundary_* 为 0-based 边界；genome_side=genome_left 表示基因组在边界左侧（插入从边界开始）。",
        ],
    )


def write_locus_regions(regions, output_path):
    """写位点区域表（区域的类别单位，供后续按类组装）。

    输入示例: regions=[{"region_id": "R01", "chrom": "chr1", ...}]
    输出示例: region_id  region_label  chrom  start  end  n_sites  n_reads
              genome_sides  site_ids  exog_names
    """
    columns = [
        "region_id",
        "region_label",
        "chrom",
        "start",
        "end",
        "n_sites",
        "n_reads",
        "genome_sides",
        "site_ids",
        "exog_names",
    ]
    rows = []
    for region in regions:
        rows.append(
            [
                region["region_id"],
                region["region_label"],
                region["chrom"],
                region["start"],
                region["end"],
                len(region["sites"]),
                len(region["read_names"]),
                ",".join(sorted({site["genome_side"] for site in region["sites"]})),
                ",".join(site["site_id"] for site in region["sites"]),
                ",".join(sorted(region["vector_names"])),
            ]
        )
    write_tsv(
        output_path,
        columns,
        rows,
        [
            "区域 = 同染色体上位点边界相距 <= --region-gap-bp（默认 200 kb）的位点集合；"
            "区域是后续“按接头位点分类”的类别单位。",
            "坐标 start/end 为 0-based 边界范围。",
        ],
    )


def write_fastq_subsets(reads_path, read_names_by_output_path):
    """一次扫描 FASTQ，把所有子集写到各自文件（保持原始 read 顺序）。

    输入示例: write_fastq_subsets("input.fastq", {
                 "out/related.fastq": {"r1", "r2"},
                 "out/junction.fastq": {"r1"},
             })
    输出示例: out/related.fastq 有 r1/r2，out/junction.fastq 有 r1
    说明：旧实现每个子集都重新扫一遍 FASTQ（related/vector/junction/每个区域各一次）；
          这里只扫一遍，内存只保存目标 read 名集合，不保存序列。
    """
    handles = {}
    total_reads = 0
    try:
        for output_path in read_names_by_output_path:
            os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
            handles[output_path] = open(output_path, "w")
        with pysam.FastxFile(reads_path) as fastq_file:
            for record in fastq_file:
                total_reads += 1
                for output_path, selected_read_names in read_names_by_output_path.items():
                    if record.name in selected_read_names:
                        handles[output_path].write(
                            "@{0}\n{1}\n+\n{2}\n".format(
                                record.name, record.sequence, record.quality
                            )
                        )
    finally:
        for handle in handles.values():
            handle.close()
    return total_reads


### 主流程 ###


@dataclass(slots=True)
class ClassifySettings:
    """``alnkit long classify`` 的全部入参（与旧脚本参数一一对应）。"""

    bam: str = ""
    reads: str = ""
    merged_fai: str = ""
    exog: tuple[str, ...] = ()
    outdir: str = "."
    min_mapq: int = 20
    min_exog_block_bp: int = 100
    contig_edge_bp: int = 100
    short_genome_block_bp: int = 500
    junction_max_gap_bp: int = 10_000
    site_gap_bp: int = 20
    region_gap_bp: int = 200_000


def run(settings: ClassifySettings) -> dict[str, object]:
    os.makedirs(settings.outdir, exist_ok=True)
    exog_names = tuple(name.strip() for name in settings.exog if name and name.strip())
    if not exog_names:
        raise UsageError("--exog 至少要给一个非空的外源序列名")
    chrom_lengths = load_fai_lengths(settings.merged_fai)

    # 第一步：只找有合格外源比对的 read（通常远少于全部 reads）
    target_read_names = find_target_read_names(
        settings.bam, exog_names, settings.min_mapq, settings.min_exog_block_bp
    )
    # 第二步：只加载这些 read 的非 secondary 比对段
    segments_by_read = load_segments_by_read(settings.bam, chrom_lengths, target_read_names)
    read_lengths = {
        read_name: segments[0]["read_length"] for read_name, segments in segments_by_read.items()
    }

    ### 1. 逐段判定：外源证据 / 基因组证据 / 可疑错配 ###

    segment_rows = []
    read_records = []
    for read_name in read_lengths:
        segments = segments_by_read.get(read_name, [])
        exog_segments = []
        genome_segments = []
        dropped_genome_segments = []
        for segment in segments:
            role, kept, reason = classify_segment(segment, exog_names, settings)
            segment_rows.append(
                [
                    read_name,
                    role,
                    segment["chrom"],
                    segment["chrom_length"],
                    segment["ref_start"],
                    segment["ref_end"],
                    segment["read_start"],
                    segment["read_end"],
                    segment["strand"],
                    segment["mapq"],
                    segment["ref_aln_len"],
                    segment["read_aln_len"],
                    segment["cigar"],
                    int(kept),
                    reason,
                ]
            )
            if not kept:
                if role == "genome":
                    dropped_genome_segments.append(segment)
                continue
            if role == "exog":
                exog_segments.append(segment)
            elif role == "genome":
                genome_segments.append(segment)

        ### 2. 接头事件 ###

        events = (
            find_junction_events(
                read_name, exog_segments, genome_segments, settings.junction_max_gap_bp
            )
            if exog_segments and genome_segments
            else []
        )

        if exog_segments and genome_segments:
            read_category = "junction"
        elif exog_segments:
            read_category = "vector_region"
        elif genome_segments:
            read_category = "genome_only"
        else:
            read_category = "other"

        notes = []
        if dropped_genome_segments:
            notes.append("genome_blocks_dropped_as_contig_edge_artifact")
        if read_category == "junction" and not events:
            notes.append("genome_and_vector_not_adjacent")

        genome_intervals = merge_intervals(
            [(segment["read_start"], segment["read_end"]) for segment in genome_segments]
        )
        genome_aln_bp = sum(end - start for start, end in genome_intervals)
        exog_aln_bp = sum(segment["ref_aln_len"] for segment in exog_segments)
        all_mapqs = [segment["mapq"] for segment in segments]

        read_records.append(
            {
                "read_name": read_name,
                "read_length": read_lengths[read_name],
                "read_category": read_category,
                "events": events,
                "n_exog_segments": len(exog_segments),
                "exog_names": sorted({segment["chrom"] for segment in exog_segments}),
                "exog_aln_bp": exog_aln_bp,
                "n_genome_segments": len(genome_segments),
                "genome_aln_bp": genome_aln_bp,
                "genome_chroms": sorted({segment["chrom"] for segment in genome_segments}),
                "n_genome_segments_dropped": len(dropped_genome_segments),
                "genome_aln_bp_dropped": sum(
                    segment["ref_aln_len"] for segment in dropped_genome_segments
                ),
                "max_mapq": max(all_mapqs) if all_mapqs else 0,
                "notes": notes,
            }
        )

    ### 3. 位点、区域聚类，并把事件写回 read ###

    all_events = [event for record in read_records for event in record["events"]]
    sites = cluster_events_into_sites(all_events, settings.site_gap_bp)
    regions = group_sites_into_regions(sites, settings.region_gap_bp)

    for site in sites:
        for event in site["events"]:
            event["site_id"] = site["site_id"]
            event["region_id"] = site["region_id"]

    ### 4. 输出 ###

    write_tsv(
        os.path.join(settings.outdir, "alignment_segments.tsv"),
        [
            "read_name",
            "segment_role",
            "chrom",
            "chrom_length",
            "ref_start",
            "ref_end",
            "read_start",
            "read_end",
            "strand",
            "mapq",
            "ref_aln_len",
            "read_aln_len",
            "simplified_cigar",
            "kept",
            "drop_reason",
        ],
        segment_rows,
        [
            "只列有外源证据的 read 的非 secondary 比对段；kept=0 的段不参与分类。",
            "read_start/read_end 已统一到正向 read 轴（见 scripts/alignment_table.py）。",
            "simplified_cigar：把两个 M 之间总长 <=5 bp 的 I/D 并进 M 后的简化视图，"
            "不是 BAM 原始 CIGAR；原始 CIGAR 未保存在本表。",
        ],
    )
    write_read_classification(
        read_records, os.path.join(settings.outdir, "read_classification.tsv")
    )
    write_junction_events(all_events, os.path.join(settings.outdir, "junction_events.tsv"))
    write_junction_sites(sites, os.path.join(settings.outdir, "junction_sites.tsv"))
    write_locus_regions(regions, os.path.join(settings.outdir, "locus_regions.tsv"))

    related_read_names = {
        record["read_name"]
        for record in read_records
        if record["read_category"] in ("junction", "vector_region")
    }
    vector_read_names = {
        record["read_name"] for record in read_records if record["read_category"] == "vector_region"
    }
    junction_read_names = {
        record["read_name"] for record in read_records if record["read_category"] == "junction"
    }

    regions_dir = os.path.join(settings.outdir, "regions")
    os.makedirs(regions_dir, exist_ok=True)
    read_names_by_output_path = {
        os.path.join(settings.outdir, "related_reads.fastq"): related_read_names,
        os.path.join(settings.outdir, "vector_region_reads.fastq"): vector_read_names,
        os.path.join(settings.outdir, "junction_reads.fastq"): junction_read_names,
    }
    for region in regions:
        region_dir = os.path.join(regions_dir, region["region_id"])
        os.makedirs(region_dir, exist_ok=True)
        read_names_by_output_path[os.path.join(region_dir, "junction_reads.fastq")] = region[
            "read_names"
        ]
    total_reads = write_fastq_subsets(settings.reads, read_names_by_output_path)

    ### 5. 输出汇总（保留“非目标 reads 被跳过”的数量，不输出逐条记录） ###

    target_read_count = len(read_records)
    junction_count = len(junction_read_names)
    vector_count = len(vector_read_names)
    non_target_count = max(0, total_reads - target_read_count)
    target_fraction = (target_read_count / total_reads) if total_reads else 0.0
    write_tsv(
        os.path.join(settings.outdir, "classification_summary.tsv"),
        [
            "total_reads",
            "target_reads",
            "junction",
            "vector_region",
            "non_target_reads",
            "target_fraction",
        ],
        [
            [
                total_reads,
                target_read_count,
                junction_count,
                vector_count,
                non_target_count,
                "{0:.4f}".format(target_fraction),
            ]
        ],
        [
            "total_reads：FASTQ 总 reads（一次扫描统计）；target_reads：有合格外源比对、进入分类的 reads；",
            "non_target_reads = total_reads - target_reads；这些 reads 不加载比对段、不输出逐条记录。",
        ],
    )

    ### 6. 控制台摘要 ###

    category_counts = collections.Counter(record["read_category"] for record in read_records)
    print("有外源证据的 read 分类：{0}".format(dict(category_counts)))
    print(
        "FASTQ 总 reads：{0}；目标 reads：{1}；非目标跳过：{2}（不加载、不输出）".format(
            total_reads, target_read_count, non_target_count
        )
    )
    print("接头位点 {0} 个、区域 {1} 个".format(len(sites), len(regions)))
    for region in regions:
        print(
            "  {0}  {1}  sites={2}  reads={3}  exog={4}".format(
                region["region_id"],
                region["region_label"],
                len(region["sites"]),
                len(region["read_names"]),
                ",".join(sorted(region["vector_names"])),
            )
        )
    print(
        "相关 reads：{0} 条（载体区域 {1}、接头 {2}）".format(
            len(related_read_names), len(vector_read_names), len(junction_read_names)
        )
    )

    return {
        "outdir": settings.outdir,
        "total_reads": total_reads,
        "target_reads": target_read_count,
        "junction_reads": len(junction_read_names),
        "vector_reads": len(vector_read_names),
        "sites": len(sites),
        "regions": len(regions),
    }
