#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# A6：用 read-vs-read overlap 关系做贪心延伸组装（overlap_extension）
#
# 来源：upstream/scripts/overlap_assembly.py（370 行）的忠实移植。逻辑未改；
# 参数解析搬到 alnkit.cli.long、库层收 AssembleSettings、main() 改为 run()。
# 产物（contigs.fasta / assembly_stats.tsv / 汇总表）与旧脚本逐字节一致。
#
# **不输出 layout.tsv**：它没有任何下游消费者（`long locate` 只读 contigs.fasta），
# 而独立校验发现它记录的 contig 内坐标与最终 contig 序列逐位一致率只有约 28%
# （贪心延伸把 read 拼上去时的坐标记账问题，非本次移植引入；上游同样如此）。
# 与其发布一个“看不出来错在哪”的表，不如不发布；需要排障时可直接跑上游脚本
# `upstream/scripts/overlap_assembly.py` 复现该表。
#
# 思路：
#   每个类先拿一条最长的“接头 read”作种子（它带着参考基因组侧翼），
#   然后反复从 overlap 边里挑一条能把手上的 contig 往左/往右延长的 read，
#   把它的非重叠悬垂序列接到 contig 两端。没有 overlap 支持的 read 不会被硬塞进来。
#
# 输入：
#   05_breakpoints/classes/<class_id>/reads.fastq        该类 read 池
#   05_breakpoints/classes/<class_id>/junction_reads.fastq 接头 reads（优先作种子）
#   04_overlap_graph/overlaps.paf                          全局 all-vs-all PAF
# 输出：
#   06_assembly/<class_id>/overlap_extension/contigs.fasta
#   06_assembly/<class_id>/overlap_extension/assembly_stats.tsv
#   06_assembly/overlap_assembly_summary.tsv
#   06_assembly/all_overlap_contigs.fasta
#
# 输入示例：
#   --classes-dir 05_breakpoints/classes --paf 04_overlap_graph/overlaps.paf
#   --outdir 06_assembly --min-overlap 500 --min-identity 0.90
# 输出示例（contigs.fasta）：
#   >oe_contig1
#   ACGT...（contig1 长度 38291 bp）
# ---------------------------------------------------------------------------
from __future__ import annotations

import collections
import os
from dataclasses import dataclass
from pathlib import Path

import pysam

from ..core.fasta import reverse_complement, write_fasta_records
from ..core.logging import get_logger
from ..core.paf import parse_paf
from ..core.tables import write_legacy_tsv as write_tsv
from . import support

### 读取输入 ###


def read_fastq_sequences(fastq_path):
    """读 fastq，返回 {read_name: 大写序列}。

    输入示例: read_fastq_sequences("classes/C03/reads.fastq")
    输出示例: {"m.../ccs": "ACGT...", ...}
    """
    return {record.name: record.sequence.upper() for record in pysam.FastxFile(fastq_path)}


def load_best_blocks(paf_path, read_names, min_overlap_bp, min_identity):
    """读 PAF，只保留类内 read 对，并且每对只留 block_bp 最大的一条。

    输入示例: load_best_blocks("overlaps.paf", {"r1","r2"}, 500, 0.90)
    输出示例: {("r1","r2"): {"query": "r1", "target": "r2", "block_bp": 12000, ...}}
    说明：重复序列可能让同一对 read 有多条比对；这里只留最好的一条，
          拿不准的多重位置会在贪心延伸时被判为 ambiguous 并跳过。
    """
    best_blocks = {}
    for block in parse_paf(paf_path):
        if block["query"] == block["target"]:
            continue
        if block["query"] not in read_names or block["target"] not in read_names:
            continue
        if block["block_bp"] < min_overlap_bp or block["identity"] < min_identity:
            continue
        pair = tuple(sorted((block["query"], block["target"])))
        previous = best_blocks.get(pair)
        if previous is None or (block["block_bp"], block["identity"]) > (
            previous["block_bp"],
            previous["identity"],
        ):
            best_blocks[pair] = block
    return best_blocks


def read_junction_read_names(junction_fastq_path):
    """读接头 reads 的名字，用于优先选种子。

    输入示例: read_junction_read_names("classes/C03/junction_reads.fastq")
    输出示例: {"m.../ccs", ...}
    """
    return {record.name for record in pysam.FastxFile(junction_fastq_path)}


### 由 overlap 推算 read 在 contig 上的位置 ###


def placement_from_block(
    placed_name,
    placed_start,
    placed_direction,
    placed_length,
    unplaced_name,
    unplaced_length,
    block,
):
    """根据一条 overlap，推算未放置 read 在 contig 坐标系里的起点与方向。

    输入示例：
      placed_name="p", placed_start=0, placed_direction=+1, placed_length=100
      unplaced_name="u", unplaced_length=100
      block={"query": "u", "target": "p", "strand": "+",
             "query_start": 80, "query_end": 100, "target_start": 0, "target_end": 20}
    输出示例: (u_start=-80, u_direction=+1)  # u 的前 80 bp 落在 contig 左边
    说明：
      * 先把 overlap 归一成“placed 坐标 -> unplaced 坐标”的线性映射；
      * 再结合 placed 在 contig 上的方向，得到 unplaced 的 contig 起点；
      * unplaced_direction=+1 表示它的序列方向与 contig 一致，-1 表示需要反向互补。
    """
    if block["query"] == unplaced_name and block["target"] == placed_name:
        # query = 未放置 read，target = 已放置 read
        placed_position_at_start = block["target_start"]
        placed_position_at_end = block["target_end"]
        unplaced_coordinate_at_start = block["query_start"]
        increasing = block["strand"] == "+"
    elif block["query"] == placed_name and block["target"] == unplaced_name:
        # query = 已放置 read，target = 未放置 read
        placed_position_at_start = block["query_start"]
        placed_position_at_end = block["query_end"]
        unplaced_coordinate_at_start = block["target_start"]
        increasing = block["strand"] == "+"
    else:
        return None

    # 与 unplaced read 起点配对的是 placed read 的哪一端？
    #   +：起点对起点；-：起点对**终点**（反链是头尾相对，起点对起点会差 aligned_bp-1）
    placed_position = placed_position_at_start if increasing else placed_position_at_end - 1

    # placed read 的那个坐标在 contig 上的位置
    if placed_direction == 1:
        contig_position_at_start = placed_start + placed_position
    else:
        contig_position_at_start = placed_start + (placed_length - 1 - placed_position)

    unplaced_direction = placed_direction * (1 if increasing else -1)
    if unplaced_direction == 1:
        unplaced_start = contig_position_at_start - unplaced_coordinate_at_start
    else:
        unplaced_start = (
            contig_position_at_start + unplaced_coordinate_at_start - (unplaced_length - 1)
        )
    return unplaced_start, unplaced_direction


def extension_lengths(unplaced_start, unplaced_length, contig_start, contig_end):
    """算这条 read 能把 contig 往左/右各延长多少 bp。

    输入示例: extension_lengths(-80, 100, 0, 100) -> 输出示例: (80, 0)
    """
    left_extension = max(0, contig_start - unplaced_start)
    right_extension = max(0, (unplaced_start + unplaced_length) - contig_end)
    return left_extension, right_extension


def choose_placement(candidate_placements, agreement_bp):
    """从多个邻居给出的位置里挑一个；位置互相矛盾就返回 None（判为 ambiguous）。

    输入示例: candidate_placements=[{"unplaced_start": -80, "unplaced_direction": 1,
                                     "left_extension": 80, "right_extension": 0,
                                     "overlap_score": 9000}, ...]
    输出示例: {"unplaced_start": -80, "unplaced_direction": 1, ...} 或 None
    说明：先把方向一致、起点相差 <= agreement_bp 的候选归为一组；
          取票数最多、延长量最大的一组；如果最大票数组和次大组打平，判为 ambiguous。
    """
    if not candidate_placements:
        return None

    groups = []
    for candidate in sorted(candidate_placements, key=lambda item: item["unplaced_start"]):
        for group in groups:
            representative = group[0]
            if (
                candidate["unplaced_direction"] == representative["unplaced_direction"]
                and abs(candidate["unplaced_start"] - representative["unplaced_start"])
                <= agreement_bp
            ):
                group.append(candidate)
                break
        else:
            groups.append([candidate])

    def group_score(group):
        representative = group[0]
        return (
            len(group),
            representative["left_extension"] + representative["right_extension"],
            representative["overlap_score"],
        )

    groups.sort(key=group_score, reverse=True)
    best_group = groups[0]
    if len(groups) > 1 and len(best_group) == len(groups[1]):
        return None
    return max(
        best_group,
        key=lambda item: (item["left_extension"] + item["right_extension"], item["overlap_score"]),
    )


### 贪心延伸主逻辑 ###


def extend_contig(
    contig_sequence, contig_start, contig_end, read_sequence, read_start, read_direction
):
    """把一条 read 的悬垂部分接到 contig 两端，返回新的 (序列, 起点, 终点)。

    输入示例: contig="ACGT"（坐标 0-4），read 放在 -2，方向 +1，read="TTACGT"
             -> 左悬垂 "TT" + "ACGT" = "TTACGT"，新起点 -2，终点 4
    输出示例: ("TTACGT", -2, 4)
    """
    oriented_sequence = read_sequence if read_direction == 1 else reverse_complement(read_sequence)
    left_extension = max(0, contig_start - read_start)
    right_extension = max(0, (read_start + len(read_sequence)) - contig_end)
    left_bases = oriented_sequence[:left_extension] if left_extension else ""
    right_bases = (
        oriented_sequence[len(oriented_sequence) - right_extension :] if right_extension else ""
    )
    return (
        left_bases + contig_sequence + right_bases,
        contig_start - left_extension,
        contig_end + right_extension,
    )


def count_junction_support(placed, read_sequences, junction, half_window_bp):
    """数一数接缝处有多少条 read 支撑（B 方案，TGS-GapCloser 的 `--min_nread` 思路）。

    为什么不是"跨接缝的 read 数"：延伸时新加的悬垂还没被任何已放置 read 覆盖，
    那样数出来永远是 0。这里用**接缝处的覆盖深度**做代理：

        support = 1（本次要放进去的 read）+ 已放置 read 中覆盖
                  [junction ± half_window] 的条数

    含义：这次延伸不是"一条 read 的孤独末端"，而是接在已经被若干 read 覆盖的序列上。
    深度骤降（重复区、嵌合）时这个数会掉下来，从而拦住盲目的延伸。

    输入示例: placed={"r1": (0, 1)}, read_sequences={"r1": "ACGT"*500}, junction=1500
    输出示例: 2（r1 覆盖接缝 + 新 read 自己）
    """
    left, right = junction - half_window_bp, junction + half_window_bp
    overlapping = 0
    for name, (start, _direction) in placed.items():
        # placed 里的 start 恒为“定向序列首碱基”的 contig 坐标，与 +/− 无关
        # （extend_contig 与 placed_reads 都按 [start, start+len) 切片，别按方向翻）。
        read_left = start
        read_right = start + len(read_sequences[name]) - 1
        if read_left <= right and read_right >= left:
            overlapping += 1
    return 1 + overlapping


def junction_depth_ratio(placed, read_sequences, junction, half_window_bp):
    """接缝局部深度 ÷ 该 contig 的平均深度（相对深度，抗测序深度差异）。

    绝对邻居数在接缝处按构造恒为 1（造出边界的那条 read 必然覆盖接缝，更深的 read
    覆盖不到还没拼进来的悬垂），所以**只能用相对量**：本函数把“接缝处有几条 read 覆盖”
    和“这条 contig 平均有几条 read 叠着”相比。实测（本仓基准数据 38392 次候选评估）：
    比值 min=0.057 / p10=0.53 / p50=1.71 / p90=2.62，门槛取 0.3 只拒掉 5.3%。

    输入示例: placed={"r1": (0, 1), "r2": (500, 1)}, 两条各 1000 bp, junction=900
    输出示例: 2 / mean_depth(≈2.0) ≈ 1.0
    """
    neighbours = count_junction_support(placed, read_sequences, junction, half_window_bp) - 1
    left = min(start for start, _ in placed.values())
    right = max(start + len(read_sequences[name]) - 1 for name, (start, _) in placed.items())
    span = max(1, right - left + 1)
    mean_depth = sum(len(read_sequences[name]) for name in placed) / span
    return neighbours / mean_depth if mean_depth > 0 else 0.0


def assemble_one_class(
    read_sequences,
    junction_read_names,
    best_blocks,
    agreement_bp,
    bridge_window_bp=200,
    min_junction_depth_ratio=0.0,
    statistics=None,
):
    """对一个类做贪心延伸，返回 contig 列表与布局信息。

    输入示例: read_sequences={"j1": "...", "v1": "..."}, junction_read_names={"j1"},
             best_blocks={("j1","v1"): {block}}
    输出示例: [{"name": "oe_contig1", "sequence": "...",
                "n_reads": 12}]
    说明：
      * 种子优先选最长的接头 read（带基因组锚点）；接头 read 用完后才用最长载体 read；
      * 每一轮挑“延长量最大”的候选；位置互相矛盾（重复序列）的 read 先跳过；
      * 没有 overlap 连到当前 contig 的 read 会另起一条 contig（不硬拼）。
    """
    if statistics is None:
        statistics = collections.Counter()
    pair_blocks = collections.defaultdict(list)
    for pair, block in best_blocks.items():
        pair_blocks[pair].append(block)

    unplaced = set(read_sequences)
    contigs = []

    while unplaced:
        seed_candidates = unplaced & junction_read_names
        seed_pool = seed_candidates if seed_candidates else unplaced
        # 排序后再取最大：并列时按名字定序，避免 set 迭代顺序让结果随 PYTHONHASHSEED 漂
        seed_name = max(sorted(seed_pool), key=lambda name: len(read_sequences[name]))
        unplaced.remove(seed_name)

        contig_sequence = read_sequences[seed_name]
        contig_start, contig_end = 0, len(contig_sequence)
        placed = {seed_name: (0, 1)}
        placed_reads = [[seed_name, 0, len(contig_sequence), "+", 1]]

        progress = True
        while progress:
            progress = False
            candidates = []
            for unplaced_name in sorted(unplaced):
                candidate_placements = []
                for placed_name, (placed_start, placed_direction) in placed.items():
                    pair = tuple(sorted((unplaced_name, placed_name)))
                    for block in pair_blocks.get(pair, []):
                        placement = placement_from_block(
                            placed_name,
                            placed_start,
                            placed_direction,
                            len(read_sequences[placed_name]),
                            unplaced_name,
                            len(read_sequences[unplaced_name]),
                            block,
                        )
                        if placement is None:
                            continue
                        unplaced_start, unplaced_direction = placement
                        left_extension, right_extension = extension_lengths(
                            unplaced_start,
                            len(read_sequences[unplaced_name]),
                            contig_start,
                            contig_end,
                        )
                        # 接缝 = 本次重叠在 contig 轴上的边界（往左延伸时取左端，否则右端）
                        if left_extension:
                            junction = max(unplaced_start, contig_start)
                        else:
                            junction = min(
                                unplaced_start + len(read_sequences[unplaced_name]) - 1,
                                contig_end - 1,
                            )
                        candidate_placements.append(
                            {
                                "unplaced_start": unplaced_start,
                                "unplaced_direction": unplaced_direction,
                                "left_extension": left_extension,
                                "right_extension": right_extension,
                                "overlap_score": block["block_bp"] * block["identity"],
                                "junction": junction,
                            }
                        )
                chosen = choose_placement(candidate_placements, agreement_bp)
                if chosen is None:
                    continue
                if min_junction_depth_ratio > 0:
                    ratio = junction_depth_ratio(
                        placed, read_sequences, chosen["junction"], bridge_window_bp // 2
                    )
                    if ratio < min_junction_depth_ratio:
                        statistics["rejected_by_depth_ratio"] += 1
                        continue
                candidates.append((chosen, unplaced_name))

            if not candidates:
                break

            chosen, unplaced_name = max(
                candidates,
                key=lambda item: (
                    item[0]["left_extension"] + item[0]["right_extension"],
                    item[0]["overlap_score"],
                ),
            )
            contig_sequence, contig_start, contig_end = extend_contig(
                contig_sequence,
                contig_start,
                contig_end,
                read_sequences[unplaced_name],
                chosen["unplaced_start"],
                chosen["unplaced_direction"],
            )
            # placed 里的坐标始终是“以种子为原点”的坐标系（contig_start 可为负），
            # 不能随 contig_start 平移，否则会和 extension_lengths 的基准错位。
            # 平移到最终 contig 坐标只在收尾时做一次。
            placed[unplaced_name] = (chosen["unplaced_start"], chosen["unplaced_direction"])
            unplaced.remove(unplaced_name)
            placed_reads.append(
                [
                    unplaced_name,
                    chosen["unplaced_start"],
                    chosen["unplaced_start"] + len(read_sequences[unplaced_name]),
                    "+" if chosen["unplaced_direction"] == 1 else "-",
                    len(placed_reads) + 1,
                ]
            )
            progress = True

        # 记账坐标是“以种子为原点”的全局坐标，左侧延伸后可能出现负值；
        # 这里统一平移到最终 contig 的 0-based 坐标（只用于统计，不再写表）。
        coordinate_shift = -contig_start
        if coordinate_shift:
            for row in placed_reads:
                row[1] += coordinate_shift
                row[2] += coordinate_shift
        for row in placed_reads:
            if row[1] < 0 or row[2] > len(contig_sequence):
                raise RuntimeError(
                    "layout coordinate out of contig: {0} not in [0, {1}]".format(
                        (row[1], row[2]), len(contig_sequence)
                    )
                )

        contigs.append(
            {
                "name": "oe_contig{0}".format(len(contigs) + 1),
                "sequence": contig_sequence,
                "placed": placed_reads,
                "n_reads": len(placed),
            }
        )

    return contigs


### 主流程 ###


logger = get_logger("assemble")


@dataclass(slots=True)
class AssembleSettings:
    """``alnkit long assemble`` 的全部入参（与旧脚本参数一一对应）。"""

    classes_dir: str = ""
    paf: str = ""
    outdir: str = "."
    min_overlap: int = 500
    min_identity: float = 0.90
    agreement_bp: int = 200
    min_junction_depth_ratio: float = 0.0
    bridge_window_bp: int = 200
    support_report: bool = True
    threads: int = 4


def run(settings: AssembleSettings) -> dict[str, object]:
    os.makedirs(settings.outdir, exist_ok=True)

    statistics = collections.Counter()
    class_names = sorted(
        name
        for name in os.listdir(settings.classes_dir)
        if os.path.isdir(os.path.join(settings.classes_dir, name))
    )
    all_contigs = {}
    summary_rows = []

    for class_name in class_names:
        class_dir = os.path.join(settings.classes_dir, class_name)
        read_sequences = read_fastq_sequences(os.path.join(class_dir, "reads.fastq"))
        junction_read_names = read_junction_read_names(
            os.path.join(class_dir, "junction_reads.fastq")
        )
        best_blocks = load_best_blocks(
            settings.paf, set(read_sequences), settings.min_overlap, settings.min_identity
        )
        contigs = assemble_one_class(
            read_sequences,
            junction_read_names,
            best_blocks,
            settings.agreement_bp,
            min_junction_depth_ratio=settings.min_junction_depth_ratio,
            bridge_window_bp=settings.bridge_window_bp,
            statistics=statistics,
        )

        output_dir = os.path.join(settings.outdir, class_name, "overlap_extension")
        os.makedirs(output_dir, exist_ok=True)

        contig_records = {}
        for contig in contigs:
            contig_records[contig["name"]] = contig["sequence"]
            all_contigs["{0}|{1}".format(class_name, contig["name"])] = contig["sequence"]

        write_fasta_records(contig_records, os.path.join(output_dir, "contigs.fasta"))

        placed_read_count = sum(contig["n_reads"] for contig in contigs)
        total_bp = sum(len(contig["sequence"]) for contig in contigs)
        longest_bp = max((len(contig["sequence"]) for contig in contigs), default=0)
        write_tsv(
            os.path.join(output_dir, "assembly_stats.tsv"),
            [
                "class_id",
                "n_contigs",
                "assembly_bp",
                "longest_contig_bp",
            ],
            [[class_name, len(contigs), total_bp, longest_bp]],
        )
        summary_rows.append([class_name, len(contigs), total_bp, longest_bp])

        print(
            "{0}: contigs={1} total={2} bp longest={3} bp reads={4}".format(
                class_name, len(contigs), total_bp, longest_bp, placed_read_count
            )
        )

    write_fasta_records(all_contigs, os.path.join(settings.outdir, "all_overlap_contigs.fasta"))
    write_tsv(
        os.path.join(settings.outdir, "overlap_assembly_summary.tsv"),
        [
            "class_id",
            "n_contigs",
            "assembly_bp",
            "longest_contig_bp",
        ],
        summary_rows,
    )

    support_rows = []
    if settings.support_report:
        minimap2 = support.find_minimap2()
        if minimap2 is None:
            logger.warning("没找到 minimap2，跳过组装支撑度报告（--no-support-report 可静默跳过）")
        else:
            support_rows = support.collect_all_support(
                Path(settings.outdir),
                Path(settings.classes_dir),
                minimap2,
                threads=max(1, settings.threads),
                class_names=class_names,
            )
            support_path = Path(settings.outdir) / "assembly_support.tsv"
            support.write_support_report(support_rows, support_path)
            logger.info("组装支撑度 -> %s（%d 条 contig）", support_path, len(support_rows))

    return {
        "outdir": settings.outdir,
        "classes": len(class_names),
        "contigs": len(all_contigs),
        "support_rows": len(support_rows),
        "rejected_by_depth_ratio": statistics.get("rejected_by_depth_ratio", 0),
        "support_status": collections.Counter(row.status for row in support_rows),
    }
