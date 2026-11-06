#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# A4：read-vs-read overlap 边与连通分量
#
# 来源：upstream/scripts/overlap_graph.py（124 行）的移植。逻辑未改；
# 参数解析搬到 alnkit.cli.long、库层收 GraphSettings、main() 改为 run()。
# 三张产物与旧脚本逐字节一致。
#
# 输入：03_classify_reads/related_reads.fastq 的 all-vs-all PAF
#       03_classify_reads/read_classification.tsv
# 输出：04_overlap_graph/overlap_edges.tsv
#       04_overlap_graph/overlap_components.tsv
#       04_overlap_graph/overlap_summary.txt
# ---------------------------------------------------------------------------
from __future__ import annotations

import collections
import os
from dataclasses import dataclass

from ..core.paf import parse_paf
from ..core.tables import read_legacy_tsv as read_tsv
from ..core.tables import write_legacy_tsv as write_tsv


def load_read_categories(classification_path):
    """读 read 类别，返回 {read_name: read_category}。"""
    header, rows = read_tsv(classification_path)
    index = {name: i for i, name in enumerate(header)}
    return {row[index["read_name"]]: row[index["read_category"]] for row in rows}


def build_edges(blocks, min_overlap_bp, min_identity):
    """每对 read 只保留 block_bp 最大的一条 overlap；返回边与邻接表。"""
    best_by_pair = {}
    for block in blocks:
        if block["query"] == block["target"]:
            continue
        if block["block_bp"] < min_overlap_bp or block["identity"] < min_identity:
            continue
        pair = tuple(sorted((block["query"], block["target"])))
        previous = best_by_pair.get(pair)
        if previous is None or (block["block_bp"], block["identity"]) > (
            previous["block_bp"],
            previous["identity"],
        ):
            best_by_pair[pair] = block

    edges = []
    adjacency = collections.defaultdict(set)
    for pair, block in sorted(best_by_pair.items()):
        edges.append(block)
        adjacency[block["query"]].add(block["target"])
        adjacency[block["target"]].add(block["query"])
    return edges, adjacency


def connected_components(read_names, adjacency):
    """并查集求连通分量，按大小降序。"""
    parent = {name: name for name in read_names}

    def find(name):
        while parent[name] != name:
            parent[name] = parent[parent[name]]
            name = parent[name]
        return name

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[left_root] = right_root

    for read_name, neighbors in adjacency.items():
        for neighbor in neighbors:
            union(read_name, neighbor)

    groups = collections.defaultdict(set)
    for read_name in read_names:
        groups[find(read_name)].add(read_name)
    return sorted(groups.values(), key=lambda members: (-len(members), sorted(members)[0]))


@dataclass(slots=True)
class GraphSettings:
    """``alnkit long graph`` 的全部入参（与旧脚本参数一一对应）。"""

    paf: str = ""
    classification: str = ""
    outdir: str = "."
    min_overlap: int = 300
    min_identity: float = 0.85


def run(settings: GraphSettings) -> dict[str, object]:
    os.makedirs(settings.outdir, exist_ok=True)

    categories = load_read_categories(settings.classification)
    blocks = parse_paf(settings.paf)
    edges, adjacency = build_edges(blocks, settings.min_overlap, settings.min_identity)
    components = connected_components(set(categories), adjacency)

    write_tsv(
        os.path.join(settings.outdir, "overlap_edges.tsv"),
        [
            "read_a",
            "read_b",
            "block_bp",
            "identity",
            "strand",
            "read_a_start",
            "read_a_end",
            "read_b_start",
            "read_b_end",
        ],
        [
            [
                edge["query"],
                edge["target"],
                edge["block_bp"],
                "{0:.4f}".format(edge["identity"]),
                edge["strand"],
                edge["query_start"],
                edge["query_end"],
                edge["target_start"],
                edge["target_end"],
            ]
            for edge in edges
        ],
        [
            "每对 read 只保留 block_bp 最大的一条 overlap；坐标 0-based 半开。",
            "阈值：block_bp >= --min-overlap，identity >= --min-identity。",
        ],
    )

    component_rows = []
    for index, members in enumerate(components, 1):
        component_rows.append(
            [
                "C{0:02d}".format(index),
                len(members),
                ",".join(sorted({categories.get(name, "-") for name in members})),
                sum(1 for name in members if categories.get(name) == "vector_region"),
                ";".join(sorted(members)),
            ]
        )
    write_tsv(
        os.path.join(settings.outdir, "overlap_components.tsv"),
        [
            "component_id",
            "n_reads",
            "categories",
            "n_vector_reads",
            "read_names",
        ],
        component_rows,
        [
            "全局连通分量，只用于观察 overlap 结构；方案 A 的类别由断点共享 reads 决定。",
        ],
    )

    summary = [
        "overlap 边: {0}".format(len(edges)),
        "连通分量: {0}".format(len(components)),
        "最大分量 reads: {0}".format(len(components[0]) if components else 0),
    ]
    with open(
        os.path.join(settings.outdir, "overlap_summary.txt"), "w", encoding="utf-8"
    ) as output:
        output.write("\n".join(summary) + "\n")
    print("\n".join(summary))

    return {
        "outdir": settings.outdir,
        "edges": len(edges),
        "components": len(components),
    }
