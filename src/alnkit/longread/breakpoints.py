#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# A5：按 read 轴上“紧邻片段”定义断点，并按共享 reads 合并成类
#
# 来源：upstream/scripts/breakpoints.py（711 行）的忠实移植。逻辑未改；
# 相对原脚本只有三处包装：参数解析搬到 alnkit.cli.long、库层收 BreakpointSettings、
# main() 改为 run()。输出表 / classes 子集 FASTQ / graph_summary.txt 与旧脚本逐字节一致。
# 已知（保留不改）：class_read_pools.tsv 的 reads_fastq 列是绝对路径，换 outdir 就会变；
# 下游要可移植的路径，得等这几个三代表统一改契约（单独立项）。
#
# 与旧 200 kb 区域分类的区别：
#   * 末端类型看 read 轴上紧邻的片段：
#       - 基因组–外源 -> vector 断点（S##）
#       - 基因组–基因组 -> genome 断点（G##a/G##b）
#       - 其他/无 -> 不进入断点图
#   * 两个断点之间 >= --min-shared 条直接共享 reads -> 合并为同一类（方案 A）。
#   * 3000000（genome_right）会被识别为 G01a：
#     它紧邻的是 genome B，不是外源；B 的另一端才是 S04。
#
# 输入：
#   03_classify_reads/alignment_segments.tsv   目标 read 的比对段（kept=1）
#   03_classify_reads/read_classification.tsv  read 类别与长度
#   04_overlap_graph/overlap_edges.tsv         read-read overlap 边
#   03_classify_reads/related_reads.fastq      输出每类 reads 子集
#   01_merge_ref/merged_reference.fa.fai       取序列长度（50bp 窗口截断）
# 输出（到 05_breakpoints/）：
#   breakpoint_sites.tsv        vector/genome 断点节点
#   breakpoint_classes.tsv      方案 A 的类（供 07 丝带图使用）
#   breakpoint_shared_reads.tsv 断点对的共享 reads / unique overlap 支持
#   segment_ends.tsv            每个 genome 片段末端的紧邻信息（含非朝向端）
#   class_read_pools.tsv        每类 read 池（供 06 组装）
#   classes/<class>/*.fastq     每类 reads / junction / vector
#   graph_summary.txt           控制台摘要
# ---------------------------------------------------------------------------
from __future__ import annotations

import collections
import os
from dataclasses import dataclass

import pysam

from ..core.bam_segments import load_fai_lengths
from ..core.tables import read_legacy_tsv as read_tsv
from ..core.tables import write_legacy_tsv as write_tsv

### 读取输入 ###


def load_segments_by_read(segments_path):
    """读 alignment_segments.tsv，返回 {read_name: [segment, ...]}（kept=1）。

    输入示例: load_segments_by_read("alignment_segments.tsv")
    输出示例: {"r1": [{"role": "genome", "chrom": "chr1", ...}, ...]}
    """
    header, rows = read_tsv(segments_path)
    index = {name: i for i, name in enumerate(header)}
    segments_by_read = collections.defaultdict(list)
    for row in rows:
        if row[index["kept"]] != "1":
            continue
        role = row[index["segment_role"]]
        if role not in ("genome", "exog"):
            continue
        segments_by_read[row[index["read_name"]]].append(
            {
                "read_name": row[index["read_name"]],
                "role": role,
                "chrom": row[index["chrom"]],
                "ref_start": int(row[index["ref_start"]]),
                "ref_end": int(row[index["ref_end"]]),
                "read_start": int(row[index["read_start"]]),
                "read_end": int(row[index["read_end"]]),
                "strand": row[index["strand"]],
                "mapq": int(row[index["mapq"]]),
            }
        )
    for segments in segments_by_read.values():
        segments.sort(key=lambda segment: (segment["read_start"], segment["read_end"]))
    return segments_by_read


def load_classification(classification_path):
    """读 read 分类表，返回 {read_name: {"length", "category"}}。

    输入示例: load_classification("read_classification.tsv")
    输出示例: {"r1": {"length": 19000, "category": "vector_region"}}
    """
    header, rows = read_tsv(classification_path)
    index = {name: i for i, name in enumerate(header)}
    classification = {}
    for row in rows:
        classification[row[index["read_name"]]] = {
            "length": int(row[index["read_length"]]),
            "category": row[index["read_category"]],
        }
    return classification


def load_overlap_edges(edges_path):
    """读 overlap_edges.tsv，返回无向邻接表。

    输入示例: load_overlap_edges("overlap_edges.tsv")
    输出示例: {"r1": {"r2"}, "r2": {"r1"}}
    """
    header, rows = read_tsv(edges_path)
    index = {name: i for i, name in enumerate(header)}
    adjacency = collections.defaultdict(set)
    for row in rows:
        read_a = row[index["read_a"]]
        read_b = row[index["read_b"]]
        adjacency[read_a].add(read_b)
        adjacency[read_b].add(read_a)
    return adjacency


### 片段末端与 read 轴紧邻关系 ###


def segment_end_info(segment, side):
    """给 genome 片段的左/右端算出参考边界与 genome_side。

    输入示例: segment=strand "-" 的 3100000-3101000, side="left"
    输出示例: {"side": "left", "boundary": 3101000,
               "genome_side": "genome_left"}
    说明：
      * 左端（read 轴左侧）：
          strand + -> ref_start，genome_right
          strand - -> ref_end，  genome_left
      * 右端（read 轴右侧）：
          strand + -> ref_end，  genome_left
          strand - -> ref_start，genome_right
    """
    if side == "left":
        if segment["strand"] == "+":
            return {"side": side, "boundary": segment["ref_start"], "genome_side": "genome_right"}
        return {"side": side, "boundary": segment["ref_end"], "genome_side": "genome_left"}
    if segment["strand"] == "+":
        return {"side": side, "boundary": segment["ref_end"], "genome_side": "genome_left"}
    return {"side": side, "boundary": segment["ref_start"], "genome_side": "genome_right"}


def find_neighbor(segment, segments, side):
    """找 genome 片段在 read 轴左/右的紧邻片段（不要求无间隔）。

    输入示例: segment=A（read 0-5334），side="right"
    输出示例: B（read 5331-7550，genome）
    说明：左邻 = 起点更小且终点更小的片段中 read_end 最大者；
          右邻 = 起点更大且终点更大的片段中 read_start 最小者。
          这样会跳过完全包含在当前片段内部的小片段。
    """
    if side == "left":
        candidates = [
            other
            for other in segments
            if other is not segment
            and other["read_start"] < segment["read_start"]
            and other["read_end"] < segment["read_end"]
        ]
        return max(candidates, key=lambda other: other["read_end"], default=None)
    candidates = [
        other
        for other in segments
        if other is not segment
        and other["read_start"] > segment["read_start"]
        and other["read_end"] > segment["read_end"]
    ]
    return min(candidates, key=lambda other: other["read_start"], default=None)


def neighbor_gap(segment, neighbor, side):
    """算 genome 片段与紧邻片段的 read 轴间隔（重叠记 0）。

    输入示例: segment=B（5331-7550），neighbor=exog（10876-12229），side="right"
    输出示例: 3326
    """
    if neighbor is None:
        return "-"
    if side == "left":
        return max(0, segment["read_start"] - neighbor["read_end"])
    return max(0, neighbor["read_start"] - segment["read_end"])


def window_50bp(genome_side, boundary, chrom_length):
    """按 genome_side 取边界旁 50 bp 窗口（1-based 展示）。

    输入示例: window_50bp("genome_left", 3101000, 100000000)
    输出示例: "chr1:...-..." 不拼 chrom；返回 (start, end)
    """
    if genome_side == "genome_left":
        start = max(0, boundary - 50)
        end = min(chrom_length, boundary)
    else:
        start = max(0, boundary)
        end = min(chrom_length, boundary + 50)
    return start, end


def build_read_axis_events(segments_by_read, chrom_lengths):
    """扫描每条 read 的片段，生成断点模板、segment_ends 行。

    输入示例: segments_by_read={"r1": [genomeA, genomeB, exog1]}
    输出示例: (templates, genome_pair_events, segment_end_rows)
      templates: [{"junction_type": "genome"/"vector", "chrom": ..., "genome_side": ...}]
      genome_pair_events: [(read_name, template_right, template_left), ...]
      segment_end_rows: 每个 genome 片段左右端一行
    """
    templates = []
    genome_pair_events = []
    segment_end_rows = []

    for read_name, segments in segments_by_read.items():
        for segment in segments:
            if segment["role"] != "genome":
                continue
            left_neighbor = find_neighbor(segment, segments, "left")
            right_neighbor = find_neighbor(segment, segments, "right")

            # genome 片段的左端
            left_info = segment_end_info(segment, "left")
            left_adjacent_coordinate = "-"
            if left_neighbor is not None and left_neighbor["role"] == "genome":
                left_adjacent_coordinate = segment_end_info(left_neighbor, "right")["boundary"]
            left_junction_type = (
                "genome"
                if left_neighbor is not None and left_neighbor["role"] == "genome"
                else "vector"
                if left_neighbor is not None and left_neighbor["role"] == "exog"
                else "none"
            )
            left_window = window_50bp(
                left_info["genome_side"],
                left_info["boundary"],
                chrom_lengths.get(segment["chrom"], left_info["boundary"] + 50),
            )
            segment_end_rows.append(
                {
                    "read_name": read_name,
                    "chrom": segment["chrom"],
                    "segment_start": segment["ref_start"],
                    "segment_end": segment["ref_end"],
                    "strand": segment["strand"],
                    "end_side": "left",
                    "end_coordinate": left_info["boundary"],
                    "genome_side": left_info["genome_side"],
                    "adjacent_type": left_neighbor["role"] if left_neighbor else "none",
                    "adjacent_name": left_neighbor["chrom"] if left_neighbor else "-",
                    "adjacent_coordinate": left_adjacent_coordinate,
                    "junction_type": left_junction_type,
                    "gap_bp": neighbor_gap(segment, left_neighbor, "left"),
                    "window_50bp": "{0}:{1}-{2}".format(
                        segment["chrom"], left_window[0] + 1, left_window[1]
                    ),
                    "junction_id": "-",
                }
            )
            if left_neighbor is not None and left_neighbor["role"] == "exog":
                templates.append(
                    {
                        "read_name": read_name,
                        "junction_type": "vector",
                        "chrom": segment["chrom"],
                        "genome_side": left_info["genome_side"],
                        "boundary": left_info["boundary"],
                        "exog_name": left_neighbor["chrom"],
                    }
                )

            # genome 片段的右端
            right_info = segment_end_info(segment, "right")
            right_adjacent_coordinate = "-"
            if right_neighbor is not None and right_neighbor["role"] == "genome":
                right_adjacent_coordinate = segment_end_info(right_neighbor, "left")["boundary"]
            right_junction_type = (
                "genome"
                if right_neighbor is not None and right_neighbor["role"] == "genome"
                else "vector"
                if right_neighbor is not None and right_neighbor["role"] == "exog"
                else "none"
            )
            right_window = window_50bp(
                right_info["genome_side"],
                right_info["boundary"],
                chrom_lengths.get(segment["chrom"], right_info["boundary"] + 50),
            )
            segment_end_rows.append(
                {
                    "read_name": read_name,
                    "chrom": segment["chrom"],
                    "segment_start": segment["ref_start"],
                    "segment_end": segment["ref_end"],
                    "strand": segment["strand"],
                    "end_side": "right",
                    "end_coordinate": right_info["boundary"],
                    "genome_side": right_info["genome_side"],
                    "adjacent_type": right_neighbor["role"] if right_neighbor else "none",
                    "adjacent_name": right_neighbor["chrom"] if right_neighbor else "-",
                    "adjacent_coordinate": right_adjacent_coordinate,
                    "junction_type": right_junction_type,
                    "gap_bp": neighbor_gap(segment, right_neighbor, "right"),
                    "window_50bp": "{0}:{1}-{2}".format(
                        segment["chrom"], right_window[0] + 1, right_window[1]
                    ),
                    "junction_id": "-",
                }
            )
            if right_neighbor is not None and right_neighbor["role"] == "exog":
                templates.append(
                    {
                        "read_name": read_name,
                        "junction_type": "vector",
                        "chrom": segment["chrom"],
                        "genome_side": right_info["genome_side"],
                        "boundary": right_info["boundary"],
                        "exog_name": right_neighbor["chrom"],
                    }
                )
            elif right_neighbor is not None and right_neighbor["role"] == "genome":
                right_template = {
                    "read_name": read_name,
                    "junction_type": "genome",
                    "chrom": segment["chrom"],
                    "genome_side": right_info["genome_side"],
                    "boundary": right_info["boundary"],
                    "exog_name": "-",
                }
                left_template = {
                    "read_name": read_name,
                    "junction_type": "genome",
                    "chrom": right_neighbor["chrom"],
                    "genome_side": segment_end_info(right_neighbor, "left")["genome_side"],
                    "boundary": segment_end_info(right_neighbor, "left")["boundary"],
                    "exog_name": "-",
                }
                templates.extend([right_template, left_template])
                genome_pair_events.append((read_name, right_template, left_template))

    return templates, genome_pair_events, segment_end_rows


### 断点聚类与合并 ###


def make_cluster(templates):
    """把一簇模板整理成 site 记录。

    输入示例: templates=[同一 side/边界的多个 read 模板]
    输出示例: {"junction_type": "vector", "chrom": ..., "genome_side": ...,
               "boundary_min": ..., "boundary_max": ..., "boundary_median": ...,
               "read_names": {...}, "exog_names": {...}, "templates": [...]}
    """
    boundaries = sorted(template["boundary"] for template in templates)
    return {
        "junction_type": templates[0]["junction_type"],
        "chrom": templates[0]["chrom"],
        "genome_side": templates[0]["genome_side"],
        "boundary_min": boundaries[0],
        "boundary_max": boundaries[-1],
        "boundary_median": boundaries[len(boundaries) // 2],
        "read_names": {template["read_name"] for template in templates},
        "exog_names": {
            template["exog_name"]
            for template in templates
            if template.get("exog_name") not in (None, "-")
        },
        "templates": list(templates),
        "site_id": None,
        "junction_pair": "-",
    }


def cluster_templates(templates, gap_bp):
    """按 (junction_type, chrom, genome_side) 聚断点；边界差 <= gap 为同一 site。

    输入示例: cluster_templates(templates, 20)
    输出示例: [cluster, ...]（按 chrom/boundary 排序）
    """
    grouped = collections.defaultdict(list)
    for template in templates:
        grouped[(template["junction_type"], template["chrom"], template["genome_side"])].append(
            template
        )
    clusters = []
    for key in sorted(grouped):
        group = sorted(grouped[key], key=lambda template: template["boundary"])
        current = []
        for template in group:
            if current and template["boundary"] - current[-1]["boundary"] > gap_bp:
                clusters.append(make_cluster(current))
                current = []
            current.append(template)
        if current:
            clusters.append(make_cluster(current))
    clusters.sort(
        key=lambda cluster: (cluster["chrom"], cluster["boundary_median"], cluster["junction_type"])
    )
    return clusters


def assign_vector_site_ids(clusters):
    """按染色体+边界给 vector 断点分配本地 S01/S02… 编号。

    输入示例: clusters=[vector cluster...]
    输出示例: 修改 cluster["site_id"]（S01、S02…）
    说明：独立项目不复用旧流程的 S## 编号，全部按本数据坐标重新编号。
    """
    vector_clusters = [cluster for cluster in clusters if cluster["junction_type"] == "vector"]
    vector_clusters.sort(key=lambda cluster: (cluster["chrom"], cluster["boundary_median"]))
    for index, cluster in enumerate(vector_clusters, 1):
        cluster["site_id"] = "S{0:02d}".format(index)


def assign_genome_site_ids(clusters, genome_pair_events, cluster_by_template):
    """按 genome–genome 配对给 genome cluster 编号 G01a/G01b。

    输入示例: genome_pair_events=[(read, template_right, template_left), ...]
    输出示例: 修改 cluster["site_id"]；返回 genome pair 列表
    """
    cluster_index = {id(cluster): index for index, cluster in enumerate(clusters)}
    pair_counts = collections.Counter()
    for _, right_template, left_template in genome_pair_events:
        right_cluster = cluster_by_template[id(right_template)]
        left_cluster = cluster_by_template[id(left_template)]
        pair_counts[
            tuple(sorted((cluster_index[id(right_cluster)], cluster_index[id(left_cluster)])))
        ] += 1

    pairs = []
    for pair_key in pair_counts:
        pair_clusters = [clusters[index] for index in pair_key]
        pair_clusters.sort(key=lambda cluster: cluster["boundary_median"])
        pairs.append((pair_clusters, pair_counts[pair_key]))
    pairs.sort(key=lambda item: (item[0][0]["chrom"], item[0][0]["boundary_median"]))

    for index, (pair_clusters, _) in enumerate(pairs, 1):
        pair_name = "G{0:02d}".format(index)
        for suffix, cluster in zip(("a", "b"), pair_clusters):
            if cluster["site_id"] is None:
                cluster["site_id"] = pair_name + suffix
                cluster["junction_pair"] = pair_name
        pairs[index - 1] = (pair_clusters, pair_name)
    return pairs


def assign_region_ids(clusters, shared_counts, min_shared):
    """按共享 reads >= min_shared 合并断点，返回 classes 与 site->class 映射。

    输入示例: shared_counts={(S05,G01a): 5, ...}, min_shared=2
    输出示例: classes=[{"class_id": "C01", "sites": [...]}, ...]
    """
    site_ids = [cluster["site_id"] for cluster in clusters]
    parent = {site_id: site_id for site_id in site_ids}

    def find(site_id):
        while parent[site_id] != site_id:
            parent[site_id] = parent[parent[site_id]]
            site_id = parent[site_id]
        return site_id

    def union(site_a, site_b):
        root_a, root_b = find(site_a), find(site_b)
        if root_a != root_b:
            parent[root_a] = root_b

    for (site_a, site_b), count in shared_counts.items():
        if count >= min_shared:
            union(site_a, site_b)

    groups = collections.defaultdict(list)
    for cluster in clusters:
        groups[find(cluster["site_id"])].append(cluster)

    class_sites = list(groups.values())
    class_sites.sort(
        key=lambda sites: (
            min(site["chrom"] for site in sites),
            min(site["boundary_min"] for site in sites),
        )
    )
    classes = []
    site_to_class = {}
    for index, sites in enumerate(class_sites, 1):
        class_id = "C{0:02d}".format(index)
        site_ids_sorted = sorted(site["site_id"] for site in sites)
        for site in sites:
            site_to_class[site["site_id"]] = class_id
        classes.append(
            {
                "class_id": class_id,
                "sites": sites,
                "site_ids": site_ids_sorted,
                "junction_read_names": set().union(*(site["read_names"] for site in sites)),
                "exog_names": set().union(*(site["exog_names"] for site in sites)),
            }
        )
    return classes, site_to_class


### 输出 ###


def fill_segment_end_junction_ids(segment_end_rows, clusters, gap_bp=20):
    """给 segment_ends 行填 junction_id（按 type/chrom/side/boundary 匹配 site）。

    输入示例: fill_segment_end_junction_ids(rows, clusters, gap_bp=20)
    输出示例: 行里的 junction_id 被填成最近的 site_id（超出 gap_bp 保持 "-"）

    说明: gap_bp 必须与位点聚类用的 --site-gap 一致，否则两张表会自相矛盾
          （断点存在、但末端找不到它）。
    """
    for row in segment_end_rows:
        candidates = [
            cluster
            for cluster in clusters
            if cluster["junction_type"] == row["junction_type"]
            and cluster["chrom"] == row["chrom"]
            and cluster["genome_side"] == row["genome_side"]
            and abs(cluster["boundary_median"] - row["end_coordinate"]) <= gap_bp
        ]
        if candidates:
            row["junction_id"] = min(
                candidates,
                key=lambda cluster: abs(cluster["boundary_median"] - row["end_coordinate"]),
            )["site_id"]


def write_fastq_subsets(reads_path, read_names_by_output_path):
    """一次扫描 FASTQ，写出所有类别子集。"""
    handles = {}
    try:
        for output_path in read_names_by_output_path:
            os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
            handles[output_path] = open(output_path, "w")
        with pysam.FastxFile(reads_path) as fastq:
            for record in fastq:
                for output_path, selected_names in read_names_by_output_path.items():
                    if record.name in selected_names:
                        handles[output_path].write(
                            "@{0}\n{1}\n+\n{2}\n".format(
                                record.name, record.sequence, record.quality
                            )
                        )
    finally:
        for handle in handles.values():
            handle.close()


### 主流程 ###


@dataclass(slots=True)
class BreakpointSettings:
    """``alnkit long breakpoints`` 的全部入参（与旧脚本参数一一对应）。"""

    segments: str = ""
    classification: str = ""
    overlaps: str = ""
    reads: str = ""
    merged_fai: str = ""
    outdir: str = "."
    site_gap: int = 20
    min_shared: int = 2


def run(settings: BreakpointSettings) -> dict[str, object]:
    os.makedirs(settings.outdir, exist_ok=True)
    chrom_lengths = load_fai_lengths(settings.merged_fai)
    segments_by_read = load_segments_by_read(settings.segments)
    classification = load_classification(settings.classification)
    overlap_adjacency = load_overlap_edges(settings.overlaps)

    ### 1. read 轴邻接 -> 模板 / segment ends ###

    templates, genome_pair_events, segment_end_rows = build_read_axis_events(
        segments_by_read, chrom_lengths
    )

    ### 2. 断点聚类与编号 ###

    clusters = cluster_templates(templates, settings.site_gap)
    assign_vector_site_ids(clusters)
    cluster_by_template = {}
    for cluster in clusters:
        for template in cluster["templates"]:
            cluster_by_template[id(template)] = cluster
    genome_pairs = assign_genome_site_ids(clusters, genome_pair_events, cluster_by_template)

    # 没有配对的 genome cluster（极端情况）：给一个独立 G 编号
    next_genome_index = len(genome_pairs) + 1
    for cluster in clusters:
        if cluster["junction_type"] == "genome" and cluster["site_id"] is None:
            cluster["site_id"] = "G{0:02d}u".format(next_genome_index)
            cluster["junction_pair"] = "G{0:02d}".format(next_genome_index)
            next_genome_index += 1

    ### 3. 共享 reads 与方案 A 类别 ###

    supported_sites_by_read = collections.defaultdict(set)
    for cluster in clusters:
        for template in cluster["templates"]:
            supported_sites_by_read[template["read_name"]].add(cluster["site_id"])

    shared_counts = collections.Counter()
    for read_name, site_ids in supported_sites_by_read.items():
        sorted_sites = sorted(site_ids)
        for left_index in range(len(sorted_sites)):
            for right_index in range(left_index + 1, len(sorted_sites)):
                shared_counts[(sorted_sites[left_index], sorted_sites[right_index])] += 1

    support_by_site = {cluster["site_id"]: set(cluster["read_names"]) for cluster in clusters}
    shared_rows = []
    for (site_a, site_b), shared in sorted(
        shared_counts.items(), key=lambda item: (-item[1], item[0])
    ):
        unique_a = support_by_site[site_a] - support_by_site[site_b]
        unique_b = support_by_site[site_b] - support_by_site[site_a]
        smaller, larger = (
            (unique_a, unique_b) if len(unique_a) <= len(unique_b) else (unique_b, unique_a)
        )
        cross_overlap = sum(
            1 for read_name in smaller if overlap_adjacency.get(read_name, set()) & larger
        )
        shared_rows.append(
            [
                site_a,
                site_b,
                shared,
                len(support_by_site[site_a]),
                len(support_by_site[site_b]),
                cross_overlap,
                1 if shared >= settings.min_shared else 0,
            ]
        )

    classes, site_to_class = assign_region_ids(clusters, shared_counts, settings.min_shared)

    # 给 cluster 填 region_id
    for cluster in clusters:
        cluster["region_id"] = site_to_class[cluster["site_id"]]

    fill_segment_end_junction_ids(segment_end_rows, clusters, settings.site_gap)

    ### 4. 写断点表 / 类表 / read 类别归属 / segment ends ###

    assignment_rows = []
    for read_name in sorted(classification):
        site_ids = sorted(supported_sites_by_read.get(read_name, set()))
        class_ids = sorted({site_to_class[site_id] for site_id in site_ids})
        assignment_rows.append(
            [
                read_name,
                classification[read_name]["category"],
                "+".join(class_ids) or "-",
                ",".join(site_ids) or "-",
            ]
        )
    write_tsv(
        os.path.join(settings.outdir, "read_class_assignment.tsv"),
        [
            "read_name",
            "read_category",
            "class_ids",
            "site_ids",
        ],
        assignment_rows,
        [
            "class_ids：该 read 支持的断点所属的方案 A 类；载体 reads 没有断点，填 -。",
            "site_ids：该 read 支持的断点（S##/G##a/G##b）。",
        ],
    )

    site_rows = []
    for cluster in sorted(
        clusters, key=lambda item: (item["chrom"], item["boundary_median"], item["site_id"])
    ):
        site_rows.append(
            [
                cluster["site_id"],
                cluster["region_id"],
                cluster["chrom"],
                cluster["genome_side"],
                cluster["boundary_min"],
                cluster["boundary_max"],
                cluster["boundary_median"],
                len(cluster["templates"]),
                len(cluster["read_names"]),
                ",".join(sorted(cluster["exog_names"])) or "-",
                ",".join(sorted(cluster["read_names"])),
                cluster["junction_type"],
                cluster["junction_pair"],
            ]
        )
    write_tsv(
        os.path.join(settings.outdir, "breakpoint_sites.tsv"),
        [
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
            "junction_type",
            "junction_pair",
        ],
        site_rows,
        [
            "junction_type=vector：基因组末端紧邻外源；genome：基因组末端紧邻另一段基因组。",
            "G##a/G##b 是同一个 genome–genome 断点的两个末端；S## 为 vector 断点。",
        ],
    )

    class_rows = []
    for class_info in classes:
        sites = class_info["sites"]
        chroms = sorted({site["chrom"] for site in sites})
        class_rows.append(
            [
                class_info["class_id"],
                "{0}|{1}".format(class_info["class_id"], ",".join(class_info["site_ids"])),
                ",".join(chroms),
                min(site["boundary_min"] for site in sites),
                max(site["boundary_max"] for site in sites),
                len(sites),
                len(class_info["junction_read_names"]),
                ",".join(sorted({site["genome_side"] for site in sites})),
                ",".join(class_info["site_ids"]),
                ",".join(sorted(class_info["exog_names"])) or "-",
                ",".join(sorted({site["junction_type"] for site in sites})),
                ",".join(sorted(class_info["junction_read_names"])),
            ]
        )
    write_tsv(
        os.path.join(settings.outdir, "breakpoint_classes.tsv"),
        [
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
            "junction_types",
            "read_names",
        ],
        class_rows,
        [
            "方案 A：site 之间共享 reads >= --min-shared（默认 2）就合并为一类。",
            "region_id 供 07 丝带图使用；read_names 是类内 junction reads（并集）。",
        ],
    )

    write_tsv(
        os.path.join(settings.outdir, "breakpoint_shared_reads.tsv"),
        [
            "site_a",
            "site_b",
            "n_shared_reads",
            "support_a",
            "support_b",
            "cross_overlap_pairs",
            "merge",
        ],
        shared_rows,
        [
            "n_shared_reads：同时支持两个断点的 reads 数。",
            "cross_overlap_pairs：A 独有 reads 与 B 独有 reads 之间的 overlap 对数（支持度参考，不参与方案 A 合并）。",
            "merge=1 表示 n_shared_reads >= --min-shared。",
        ],
    )

    write_tsv(
        os.path.join(settings.outdir, "segment_ends.tsv"),
        [
            "read_name",
            "chrom",
            "segment_start",
            "segment_end",
            "strand",
            "end_side",
            "end_coordinate",
            "genome_side",
            "adjacent_type",
            "adjacent_name",
            "adjacent_coordinate",
            "junction_type",
            "junction_id",
            "gap_bp",
            "window_50bp",
        ],
        [
            [
                row[column]
                for column in (
                    "read_name",
                    "chrom",
                    "segment_start",
                    "segment_end",
                    "strand",
                    "end_side",
                    "end_coordinate",
                    "genome_side",
                    "adjacent_type",
                    "adjacent_name",
                    "adjacent_coordinate",
                    "junction_type",
                    "junction_id",
                    "gap_bp",
                    "window_50bp",
                )
            ]
            for row in segment_end_rows
        ],
        [
            "每个 genome 片段的左右末端一行；adjacent_type 看的是 read 轴上紧邻的片段。",
            "非朝向载体的末端也会列出（junction_type=genome/none），不会被误当 vector 断点。",
        ],
    )

    ### 5. 每类 read 池与 FASTQ ###

    vector_read_names = {
        read_name
        for read_name, info in classification.items()
        if info["category"] == "vector_region"
    }
    class_pool_rows = []
    read_names_by_output = {}
    for class_info in classes:
        class_id = class_info["class_id"]
        class_dir = os.path.join(settings.outdir, "classes", class_id)
        os.makedirs(class_dir, exist_ok=True)
        junction_reads = class_info["junction_read_names"]
        all_reads = junction_reads | vector_read_names
        supported_vector = sum(
            1
            for vector_read in vector_read_names
            if overlap_adjacency.get(vector_read, set()) & junction_reads
        )
        pool_bp = sum(classification.get(read_name, {}).get("length", 0) for read_name in all_reads)
        reads_fastq = os.path.join(class_dir, "reads.fastq")
        class_pool_rows.append(
            [
                class_id,
                len(junction_reads),
                len(vector_read_names),
                len(all_reads),
                pool_bp,
                supported_vector,
                reads_fastq,
            ]
        )
        read_names_by_output[reads_fastq] = all_reads
        read_names_by_output[os.path.join(class_dir, "junction_reads.fastq")] = junction_reads
        read_names_by_output[os.path.join(class_dir, "vector_reads.fastq")] = vector_read_names

    write_tsv(
        os.path.join(settings.outdir, "class_read_pools.tsv"),
        [
            "class_id",
            "n_junction_reads",
            "n_vector_reads",
            "n_total_reads",
            "pool_read_bp",
            "n_vector_with_overlap_support",
            "reads_fastq",
        ],
        class_pool_rows,
        [
            "类池 = 该类断点的 junction reads（并集） + 全部载体区域 reads。",
            "pool_read_bp：类池所有 reads 长度之和（不是组装长度）。",
        ],
    )
    write_fastq_subsets(settings.reads, read_names_by_output)

    ### 6. 摘要 ###

    summary_lines = []
    summary_lines.append(
        "断点图方案 A：site 合并阈值 >= {0} 条共享 reads".format(settings.min_shared)
    )
    summary_lines.append("sites: {0}".format(len(clusters)))
    summary_lines.append("classes: {0}".format(len(classes)))
    summary_lines.append(
        "genome–genome 配对: {0}".format(
            "; ".join(
                "{0}({1})".format(pair_name, ",".join(site["site_id"] for site in pair_clusters))
                for pair_clusters, pair_name in genome_pairs
            )
            or "-"
        )
    )
    for class_info in classes:
        summary_lines.append(
            "  {0}: {1}  (junction reads {2})".format(
                class_info["class_id"],
                ",".join(class_info["site_ids"]),
                len(class_info["junction_read_names"]),
            )
        )
    with open(
        os.path.join(settings.outdir, "graph_summary.txt"), "w", encoding="utf-8"
    ) as summary_file:
        summary_file.write("\n".join(summary_lines) + "\n")
    print("\n".join(summary_lines))

    return {
        "outdir": settings.outdir,
        "sites": len(clusters),
        "classes": len(classes),
        "segment_ends": len(segment_end_rows),
    }
