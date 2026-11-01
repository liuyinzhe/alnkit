#!/usr/bin/env python3
# ---------------------------------------------------------------------------
# A7：丝带图 + 断点锚点/外源位置表（分类优先主线的最后一步）
#
# 来源：upstream/scripts/ribbon_positions.py（714 行）的移植。算法与输出不变；
# 两处包装：参数解析搬到 alnkit.cli.long、main() 改为 run()。
# 另有一处**有意改动**：旧脚本用 subprocess 调 vendored tools/synteny_ribbon.py（1.8.0），
# 这里直接调用本包的 alnkit.plot.ribbon（1.10.1 移植版）。audit 的版本影响评估已实测
# 两者对当前参数产出相同的 PNG 与 placements/regions，因此产物对账仍然成立。
#
# 前提：步骤 6 已经用 overlap_extension 给每个类拼出组装片段。
# 本步骤只做三件事：
#   1. 为每个类构造查询序列：
#      - 50 bp 断点锚点（flank）：位点边界旁 50 bp，用来肉眼确认断点；
#      - 全部外源序列（exog）。
#   2. 把类组装与查询序列比对，逐条 contig 画丝带图（组装 contig 放在最下方）。
#   3. 输出位置表：50 bp 锚点/外源在组装 contig 上的坐标、方向、一致率；
#      锚点再通过“组装 contig -> 参考基因组”的比对，还原成真实参考基因组坐标。
#
# 输入：
#   05_breakpoints/breakpoint_classes.tsv、breakpoint_sites.tsv
#   01_merge_ref/reference.fa(.fai)                     取断点旁 50 bp
#   data/exog/*.fa                                      外源序列
#   06_assembly/<class>/overlap_extension/contigs.fasta
#   tools/synteny_ribbon.py
# 输出：
#   07_ribbon_positions/<class>/queries.fasta / queries.tsv
#   07_ribbon_positions/<class>/assembly.fasta
#   07_ribbon_positions/<class>/assembly_queries.paf
#   07_ribbon_positions/figures/<class>_<contig>.png
#   07_ribbon_positions/tables/<class>_<contig>.blocks.tsv / .placements.tsv / .regions.tsv
#   07_ribbon_positions/class_assembly_positions.tsv
#   07_ribbon_positions/class_position_summary.tsv
#   07_ribbon_positions/contig_map.tsv
#   07_ribbon_positions/summary.md
#
# 输入示例：
#   --regions 05_breakpoints/breakpoint_classes.tsv --sites 05_breakpoints/breakpoint_sites.tsv
#   --reference-fasta 01_merge_ref/reference.fa --assembly-dir 06_assembly
#   --outdir 07_ribbon_positions
# 输出示例（class_assembly_positions.tsv）：
#   class_id  assembly_contig  query_name  query_type  contig_start  contig_end
#   C03       oe_contig1       flank|G01b... flank       32022        32072
#   ...       ...              ...         ...         ...           ...
#   reference_chrom  reference_start  reference_end  query_source                     distance_to_breakpoint_bp
#   chr1      2003010         2003060       chr1:2003010-2003050   10
# ---------------------------------------------------------------------------
from __future__ import annotations

import collections
import os
import shutil
import subprocess
from dataclasses import dataclass

import pysam

from ..core.errors import UsageError
from ..core.fasta import read_fasta_records, write_fasta_records
from ..core.paf import parse_paf
from ..core.tables import read_legacy_tsv as read_tsv
from ..core.tables import write_legacy_tsv as write_tsv
from ..plot import ribbon


def run_checked(command, stdout_file=None, context="command"):
    """运行外部命令；非 0 退出码直接抛错，避免空 PAF/缺图静默通过。

    输入示例: run_checked(["minimap2", "-cx", "sr", ...], paf_file, "minimap2 flank")
    输出示例: 命令成功时返回 CompletedProcess；失败时抛 RuntimeError。
    说明：stdout 可以指向已打开的文件；stderr 捕获后只保留末尾几行，避免日志爆炸。
    """
    completed = subprocess.run(
        command,
        stdout=stdout_file if stdout_file is not None else subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        check=False,
    )
    if completed.returncode != 0:
        stderr_tail = " | ".join((completed.stderr or "").strip().splitlines()[-5:])
        raise RuntimeError(
            "{0} failed with exit code {1}: {2}".format(
                context, completed.returncode, stderr_tail or "(no stderr)"
            )
        )
    return completed


### 1. 构造查询序列 ###


def load_sites_by_region(sites_path):
    """读接头位点表，按 region_id 分组。

    输入示例: load_sites_by_region("05_breakpoints/breakpoint_classes.tsv")
    输出示例: {"C03": [{"site_id": "G01b", "chrom": "chr1",
                        "genome_side": "genome_left", "boundary_median": 2000050,
                        "read_names": {"m.../ccs", ...}}, ...]}
    """
    header, rows = read_tsv(sites_path)
    index = {name: position for position, name in enumerate(header)}
    sites_by_region = collections.defaultdict(list)
    for row in rows:
        site = {
            "site_id": row[index["site_id"]],
            "region_id": row[index["region_id"]],
            "chrom": row[index["chrom"]],
            "genome_side": row[index["genome_side"]],
            "boundary_median": int(row[index["boundary_median"]]),
            "read_names": {name for name in row[index["read_names"]].split(",") if name},
        }
        sites_by_region[site["region_id"]].append(site)
    for region_sites in sites_by_region.values():
        region_sites.sort(key=lambda site: (site["chrom"], site["boundary_median"]))
    return sites_by_region


def flank_interval(site, chrom_length, flank_bp):
    """算一个接头位点的 50 bp 断点锚点区间（0-based 半开）与 1-based 展示区间。

    输入示例: site={"genome_side": "genome_right", "boundary_median": 2003000},
             chrom_length=100000000, flank_bp=50
    输出示例: (start=2003000, end=2003050, "chr1:2003010-2003050")
    说明：
      * genome_left ：基因组在边界左侧，取 [boundary - 50, boundary)；
      * genome_right：基因组在边界右侧，取 [boundary, boundary + 50)。
        这是断点右侧的基因组序列，不是“断点外”的序列；它能和组装比对上，
        因为组装里也有这段基因组；若组装从边界内侧若干 bp 才开始，偏移会体现在
        query_start 或 distance_to_breakpoint_bp 上。
    """
    boundary = site["boundary_median"]
    if site["genome_side"] == "genome_left":
        start = max(0, boundary - flank_bp)
        end = min(chrom_length, boundary)
    else:
        start = max(0, boundary)
        end = min(chrom_length, boundary + flank_bp)
    display_region = "{0}:{1}-{2}".format(site["chrom"], start + 1, end)
    return start, end, display_region


def build_queries_for_class(region_sites, reference_fasta, chrom_lengths, exog_sequences, flank_bp):
    """为一个类构造查询序列（50 bp 锚点 + 外源），并返回 fasta 与说明行。

    输入示例: region_sites=[site...],
             exog_sequences={"insertA": "ACGT...", "insertB": "ACGT..."}, flank_bp=50
    输出示例: ({"flank|S01|genome_left|chr1:2000001-2000050": "ACGT...",
                "exog|insertA": "ACGT..."},
               [{"query_name": ..., "query_type": "flank", "source_name": "S01", ...},
                {"query_name": ..., "query_type": "exog", ...}])
    """
    query_records = {}
    query_rows = []
    for site in region_sites:
        chrom_length = chrom_lengths.get(site["chrom"])
        if chrom_length is None:
            raise UsageError(f"位点表里的序列 {site['chrom']} 不在参考 .fai 里（参考用错了？）")

        # 50 bp 断点锚点：用来肉眼确认断点
        start, end, display_region = flank_interval(site, chrom_length, flank_bp)
        query_name = "flank|{0}|{1}|{2}".format(
            site["site_id"], site["genome_side"], display_region
        )
        sequence = reference_fasta.fetch(site["chrom"], start, end).upper()
        query_records[query_name] = sequence
        query_rows.append(
            {
                "query_name": query_name,
                "query_type": "flank",
                "source_name": site["site_id"],
                "source_region": display_region,
                "chrom": site["chrom"],
                "start": start,
                "end": end,
                "length": len(sequence),
            }
        )

    for exog_name, sequence in exog_sequences.items():
        query_name = "exog|{0}".format(exog_name)
        query_records[query_name] = sequence
        query_rows.append(
            {
                "query_name": query_name,
                "query_type": "exog",
                "source_name": exog_name,
                "source_region": exog_name,
                "chrom": "-",
                "start": "-",
                "end": "-",
                "length": len(sequence),
            }
        )
    return query_records, query_rows


### 2. 组装与查询比对 ###


def load_assembly_contigs(assembly_path):
    """读类组装片段；文件不存在或为空时返回 {}。

    输入示例: load_assembly_contigs("06_assembly/C03/overlap_extension/contigs.fasta")
    输出示例: {"oe_contig1": "ACGT...", "oe_contig2": "ACGT..."}
    """
    if not os.path.exists(assembly_path):
        return {}
    return read_fasta_records(assembly_path)


def align_assembly_to_queries(assembly_fasta, queries_fasta, paf_path, threads):
    """把组装片段当参考、查询序列当 query 跑 minimap2，输出合并后的 PAF。

    输入示例: align_assembly_to_queries("assembly.fasta", "queries.fasta",
                                        "assembly_queries.paf", 12)
    输出示例: PAF 文件；target = 组装片段名，query = 50 bp 锚点/参考段/外源名
    说明：
      * `exog` 是长查询，用 `-x asm20`；
      * `flank` 只有 50 bp，`asm20` 不会稳定输出，单独用 `-x sr`；
      * 两批 PAF 直接拼接，丝带图/位置表按 query 名区分，互不影响。
    """
    if not os.path.exists(assembly_fasta) or not os.path.exists(queries_fasta):
        open(paf_path, "w").close()
        return

    query_records = read_fasta_records(queries_fasta)
    flank_records = {
        name: sequence for name, sequence in query_records.items() if name.startswith("flank|")
    }
    long_records = {
        name: sequence for name, sequence in query_records.items() if not name.startswith("flank|")
    }

    output_dir = os.path.dirname(paf_path)
    paf_parts = []
    for records, preset, part_name in (
        (long_records, "asm20", "long"),
        (flank_records, "sr", "short"),
    ):
        if not records:
            continue
        part_fasta = os.path.join(output_dir, "queries_{0}.fasta".format(part_name))
        part_paf = os.path.join(output_dir, "assembly_queries_{0}.paf".format(part_name))
        write_fasta_records(records, part_fasta)
        with open(part_paf, "w") as paf_file:
            run_checked(
                [
                    "minimap2",
                    "-cx",
                    preset,
                    "-N",
                    "500",
                    "-t",
                    str(threads),
                    assembly_fasta,
                    part_fasta,
                ],
                stdout_file=paf_file,
                context="minimap2 assembly-vs-{0}".format(part_name),
            )
        paf_parts.append(part_paf)

    with open(paf_path, "w") as merged_paf:
        for part_paf in paf_parts:
            if not os.path.exists(part_paf):
                continue
            with open(part_paf) as part_file:
                merged_paf.write(part_file.read())


### 3. 组装 contig -> 参考基因组的坐标还原 ###


def load_site_boundaries(sites_path):
    """读位点边界，返回 {site_id: {"boundary": int, "genome_side": str, "chrom": str}}。

    输入示例: load_site_boundaries("05_breakpoints/breakpoint_sites.tsv")
    输出示例: {"S02": {"boundary": 2003000, "genome_side": "genome_right",
                       "chrom": "chr1"}}
    """
    header, rows = read_tsv(sites_path)
    index = {name: position for position, name in enumerate(header)}
    boundaries = {}
    for row in rows:
        boundaries[row[index["site_id"]]] = {
            "boundary": int(row[index["boundary_median"]]),
            "genome_side": row[index["genome_side"]],
            "chrom": row[index["chrom"]],
        }
    return boundaries


def align_contigs_to_reference(contigs, reference_fasta, output_dir, threads):
    """把所有类组装 contig 一次比对到参考基因组，按 contig 分组返回 PAF 块。

    输入示例: align_contigs_to_reference({"C03|oe_contig1": "ACGT..."},
                                        "reference.fa", "work", 12)
    输出示例: {"C03|oe_contig1": [{"target": "chr1", "query_start": 41000,
                                  "query_end": 42000, "target_start": 2002000, ...}]}
    说明：contig 名加 `<class>|` 前缀，避免不同类的 `oe_contig1` 重名。
    """
    os.makedirs(output_dir, exist_ok=True)
    all_fasta = os.path.join(output_dir, "all_assembly_contigs.fasta")
    with open(all_fasta, "w") as fasta_file:
        for contig_name, sequence in contigs.items():
            fasta_file.write(">{0}\n".format(contig_name))
            for start in range(0, len(sequence), 60):
                fasta_file.write(sequence[start : start + 60] + "\n")

    paf_path = os.path.join(output_dir, "all_assembly_contigs_vs_reference.paf")
    with open(paf_path, "w") as paf_file:
        run_checked(
            [
                "minimap2",
                "-cx",
                "asm20",
                "-N",
                "50",
                "-t",
                str(threads),
                reference_fasta,
                all_fasta,
            ],
            stdout_file=paf_file,
            context="minimap2 contig-vs-reference",
        )

    blocks_by_contig = collections.defaultdict(list)
    for block in parse_paf(paf_path):
        blocks_by_contig[block["query"]].append(block)
    return blocks_by_contig


def restore_reference_interval(assembly_start, assembly_end, reference_blocks):
    """把组装 contig 上的一个区间，通过 contig->参考 的比对还原成参考坐标。

    输入示例: assembly_start=32022, assembly_end=32072,
             reference_blocks=[{"target": "chr1", "query_start": 7219,
                                "query_end": 32072, "target_start": 2003000,
                                "target_end": 2005000, "strand": "-", ...}]
    输出示例: {"chrom": "chr1", "start": 2003010, "end": 2003060,
               "overlap_bp": 50, "identity": 1.0}
    说明：只取与目标区间重叠最多、且一致率最高的 contig->参考 块；
          没找到时返回 None（表里填 -）。
    """
    best = None
    for block in reference_blocks:
        overlap_start = max(assembly_start, block["query_start"])
        overlap_end = min(assembly_end, block["query_end"])
        if overlap_end <= overlap_start:
            continue
        if block["strand"] == "+":
            ref_start = block["target_start"] + (overlap_start - block["query_start"])
            ref_end = block["target_start"] + (overlap_end - block["query_start"])
        else:
            ref_start = block["target_end"] - (overlap_end - block["query_start"])
            ref_end = block["target_end"] - (overlap_start - block["query_start"])
        candidate = {
            "chrom": block["target"],
            "start": ref_start,
            "end": ref_end,
            "overlap_bp": overlap_end - overlap_start,
            "identity": block["identity"],
        }
        if best is None or (candidate["overlap_bp"], candidate["identity"]) > (
            best["overlap_bp"],
            best["identity"],
        ):
            best = candidate
    return best


def distance_to_boundary(reference_interval, site_boundary):
    """算还原区间到预期接头边界的距离（0 = 区间跨过边界）。

    输入示例: reference_interval={"start": 2003010, "end": 2003060},
             site_boundary=2003000
    输出示例: 10
    """
    if reference_interval is None:
        return "-"
    start, end = reference_interval["start"], reference_interval["end"]
    if start <= site_boundary <= end:
        return 0
    return min(abs(start - site_boundary), abs(end - site_boundary))


### 4. 丝带图与位置表 ###


def write_contig_paf(blocks, output_path):
    """把某个 contig 的比对块写成单独 PAF（供丝带图 --paf 模式使用）。

    输入示例: write_contig_paf([block...], "C03_oe_contig1.paf")
    输出示例: PAF 文件（列顺序与 minimap2 一致）
    """
    with open(output_path, "w") as paf_file:
        for block in blocks:
            paf_file.write(
                "\t".join(
                    [
                        block["query"],
                        str(block["query_length"]),
                        str(block["query_start"]),
                        str(block["query_end"]),
                        block["strand"],
                        block["target"],
                        str(block["target_length"]),
                        str(block["target_start"]),
                        str(block["target_end"]),
                        str(block["matches"]),
                        str(block["block_bp"]),
                        str(block["mapq"]),
                    ]
                )
                + "\n"
            )


def run_ribbon_tool(paf_path, contig_name, class_id, output_dirs, title):
    """调用独立的丝带图工具，生成 png 与位置表（组装 contig 放在最下方）。

    输入示例: run_ribbon_tool("C03_oe_contig1.paf", "oe_contig1", "C03", {...},
                              "C03 assembly")
    输出示例: figures/C03_oe_contig1.png、tables/C03_oe_contig1.*.tsv
    """
    figure_name = "{0}_{1}".format(class_id, contig_name)
    ribbon.run(
        ribbon.RibbonSettings(
            outdir=output_dirs["figures"],
            name=figure_name,
            paf=paf_path,
            reference_name=contig_name,
            reference_position="bottom",
            title=title,
            placement_table=os.path.join(output_dirs["tables"], figure_name + ".placements.tsv"),
            region_table=os.path.join(output_dirs["tables"], figure_name + ".regions.tsv"),
            min_block_bp=20,
            min_label_bp=20,
            label_coordinate="reference",
            min_identity=0.8,
            max_labels_per_row=12,
            dpi=150,
        )
    )
    figure_path = os.path.join(output_dirs["figures"], figure_name + ".png")
    if not os.path.exists(figure_path):
        raise RuntimeError("synteny_ribbon did not create {0}".format(figure_path))
    return figure_path


def build_position_rows(
    class_id, blocks, query_metadata, site_boundaries, reference_blocks_by_contig
):
    """把一个类的全部比对块整理成位置表行。

    输入示例: blocks=[{"target": "oe_contig2", "query": "flank|S02|...", ...}]
    输出示例: 一行 [class_id, contig, query_name, "flank", query_start, query_end,
                    query_strand, contig_start, contig_end, identity,
                    reference_chrom, reference_start, reference_end,
                    query_source, distance_to_breakpoint_bp]
    """
    rows = []
    for block in sorted(
        blocks, key=lambda item: (item["target"], item["query"], item["target_start"])
    ):
        metadata = query_metadata.get(block["query"], {})
        query_type = metadata.get("query_type", "-")
        reference_interval = None
        offset = "-"
        if query_type == "flank":
            reference_interval = restore_reference_interval(
                block["target_start"],
                block["target_end"],
                reference_blocks_by_contig.get("{0}|{1}".format(class_id, block["target"]), []),
            )
            site = site_boundaries.get(metadata.get("source_name", ""))
            if site is not None:
                offset = distance_to_boundary(reference_interval, site["boundary"])
        rows.append(
            [
                class_id,
                block["target"],
                block["query"],
                query_type,
                block["query_start"],
                block["query_end"],
                block["strand"],
                block["target_start"],
                block["target_end"],
                "{0:.4f}".format(block["identity"]),
                reference_interval["chrom"] if reference_interval else "-",
                reference_interval["start"] if reference_interval else "-",
                reference_interval["end"] if reference_interval else "-",
                metadata.get("source_region", "-"),
                offset,
            ]
        )
    return rows


def describe_blocks(query_names, query_type, blocks_by_query, position_rows, max_blocks=None):
    """把一个类型（flank/exog）的查询位置压缩成短字符串。

    输入示例: query_names=["flank|S01|..."], query_type="flank",
             blocks_by_query={...}, position_rows=[...]
    输出示例: ["S01@oe_contig1:41755-41805(-)->chr1:2000000-2000050"]
    说明：flank 取最佳的一块（block_bp 最大）并带上还原参考坐标；
          外源最多列前 max_blocks 处，其余用 `,...` 省略。
    """
    descriptions = []
    for query_name in query_names:
        query_blocks = sorted(
            blocks_by_query.get(query_name, []),
            key=lambda block: (-block["block_bp"], block["target"]),
        )
        if not query_blocks:
            descriptions.append(
                "{0}:no_hit".format(query_name.split("|")[1] if "|" in query_name else query_name)
            )
            continue
        if query_type == "exog":
            shown = query_blocks[: max_blocks or 5]
            positions = ",".join(
                "{0}:{1}-{2}({3})".format(
                    block["target"], block["target_start"], block["target_end"], block["strand"]
                )
                for block in shown
            )
            descriptions.append(
                "{0}:{1} hits [{2}{3}]".format(
                    query_name,
                    len(query_blocks),
                    positions,
                    ",..." if len(query_blocks) > len(shown) else "",
                )
            )
            continue

        best_block = query_blocks[0]
        matching_rows = [
            row
            for row in position_rows
            if row[2] == query_name
            and row[7] == best_block["target_start"]
            and row[1] == best_block["target"]
        ]
        reference_text = "ref=-"
        if matching_rows and matching_rows[0][10] != "-":
            reference_text = "{0}:{1}-{2}".format(
                matching_rows[0][10], matching_rows[0][11], matching_rows[0][12]
            )
        descriptions.append(
            "{0}@{1}:{2}-{3}({4})->{5}".format(
                query_name.split("|")[1] if "|" in query_name else query_name,
                best_block["target"],
                best_block["target_start"],
                best_block["target_end"],
                best_block["strand"],
                reference_text,
            )
        )
    return descriptions


def summarize_class(
    class_id, blocks, query_metadata, position_rows, figure_paths, n_assembly_contigs
):
    """生成一个类的汇总表行（50 bp 锚点/外源各命中多少、落在哪）。

    输入示例: blocks=[...], query_metadata={...}, position_rows=[...]
    输出示例: [class_id, n_contigs,
               n_flank_queries, n_flank_queries_hit,
               n_exog_queries, n_exog_queries_hit,
               flank_positions, exog_positions, figures]
    """
    blocks_by_query = collections.defaultdict(list)
    for block in blocks:
        blocks_by_query[block["query"]].append(block)

    flank_queries = [name for name, meta in query_metadata.items() if meta["query_type"] == "flank"]
    exog_queries = [name for name, meta in query_metadata.items() if meta["query_type"] == "exog"]

    flank_descriptions = describe_blocks(flank_queries, "flank", blocks_by_query, position_rows)
    exog_descriptions = describe_blocks(
        exog_queries, "exog", blocks_by_query, position_rows, max_blocks=5
    )

    return [
        class_id,
        n_assembly_contigs,
        len(flank_queries),
        sum(1 for name in flank_queries if blocks_by_query.get(name)),
        len(exog_queries),
        sum(1 for name in exog_queries if blocks_by_query.get(name)),
        "; ".join(flank_descriptions),
        "; ".join(exog_descriptions),
        "; ".join(figure_paths),
    ]


def write_summary_markdown(summary_rows, output_path):
    """写人可读的 summary.md（每类一小节，含图与关键位置）。

    输入示例: write_summary_markdown(summary_rows, "summary.md")
    输出示例: Markdown 文件（含每类 50 bp 锚点/外源位置与图链接）
    """
    with open(output_path, "w") as markdown_file:
        markdown_file.write("# 07 丝带图与位置汇总\n\n")
        markdown_file.write(
            "坐标说明：组装坐标为 contig 自身 0-based 半开；"
            "参考坐标由“组装 contig -> 参考基因组”的比对还原，"
            "不是直接把 50 bp 锚点的提取坐标抄过来。\n\n"
        )
        for row in summary_rows:
            (
                class_id,
                n_contigs,
                n_flank,
                n_flank_hit,
                n_exog,
                n_exog_hit,
                flank_text,
                exog_text,
                figure_text,
            ) = row
            markdown_file.write("## {0}\n\n".format(class_id))
            markdown_file.write("- 组装：`overlap_extension`，{0} 条 contig\n".format(n_contigs))
            markdown_file.write("- 50 bp 断点锚点命中：{0}/{1}\n".format(n_flank_hit, n_flank))
            markdown_file.write("- 外源命中：{0}/{1}\n\n".format(n_exog_hit, n_exog))
            markdown_file.write("50 bp 断点锚点位置：\n\n")
            for item in flank_text.split("; "):
                markdown_file.write("- `{0}`\n".format(item))
            markdown_file.write("\n外源位置：\n\n")
            for item in exog_text.split("; "):
                markdown_file.write("- `{0}`\n".format(item))
            markdown_file.write("\n图：\n\n")
            for figure_path in figure_text.split("; "):
                markdown_file.write(
                    "![]({0})\n\n".format(
                        os.path.relpath(figure_path, os.path.dirname(output_path))
                    )
                )


### 主流程 ###


@dataclass(slots=True)
class LocateSettings:
    """``alnkit long locate`` 的全部入参（与旧脚本参数一一对应）。"""

    regions: str = ""
    sites: str = ""
    reference_fasta: str = ""
    exog: tuple[str, ...] = ()
    assembly_dir: str = ""
    outdir: str = "."
    threads: int = 12
    flank_bp: int = 50


def run(settings: LocateSettings) -> dict[str, object]:
    os.makedirs(settings.outdir, exist_ok=True)
    output_dirs = {
        "figures": os.path.join(settings.outdir, "figures"),
        "tables": os.path.join(settings.outdir, "tables"),
        "work": os.path.join(settings.outdir, "work"),
    }
    # 清掉上一次运行的图和中间 PAF，避免留下旧图
    for directory in output_dirs.values():
        shutil.rmtree(directory, ignore_errors=True)
        os.makedirs(directory, exist_ok=True)

    header, region_rows = read_tsv(settings.regions)
    region_index = {name: position for position, name in enumerate(header)}
    regions = [
        {
            "region_id": row[region_index["region_id"]],
            "region_label": row[region_index["region_label"]],
        }
        for row in region_rows
    ]
    sites_by_region = load_sites_by_region(settings.sites)
    site_boundaries = load_site_boundaries(settings.sites)
    exog_sequences = {}
    for exog_path in settings.exog:
        exog_sequences.update(read_fasta_records(exog_path))

    class_data = {}
    all_contigs = {}

    ### 1. 每类：查询序列 + 组装 vs 查询比对 ###

    with pysam.FastaFile(settings.reference_fasta) as reference_fasta:
        chrom_lengths = dict(zip(reference_fasta.references, reference_fasta.lengths))
        for region in regions:
            class_id = region["region_id"]
            class_dir = os.path.join(settings.outdir, class_id)
            os.makedirs(class_dir, exist_ok=True)

            assembly_path = os.path.join(
                settings.assembly_dir, class_id, "overlap_extension", "contigs.fasta"
            )
            contigs = load_assembly_contigs(assembly_path)
            if not contigs:
                print("{0}: 没有 overlap_extension 组装，跳过".format(class_id))
                continue

            query_records, query_rows = build_queries_for_class(
                sites_by_region.get(class_id, []),
                reference_fasta,
                chrom_lengths,
                exog_sequences,
                settings.flank_bp,
            )
            queries_fasta = os.path.join(class_dir, "queries.fasta")
            write_fasta_records(query_records, queries_fasta)
            write_tsv(
                os.path.join(class_dir, "queries.tsv"),
                [
                    "query_name",
                    "query_type",
                    "source_name",
                    "source_region",
                    "chrom",
                    "start",
                    "end",
                    "length",
                ],
                [
                    [
                        row["query_name"],
                        row["query_type"],
                        row["source_name"],
                        row["source_region"],
                        row["chrom"],
                        row["start"],
                        row["end"],
                        row["length"],
                    ]
                    for row in query_rows
                ],
                [
                    "flank = 断点旁 50 bp 锚点；start/end 为 0-based 半开；exog 查询没有参考坐标。",
                ],
            )

            assembly_fasta = os.path.join(class_dir, "assembly.fasta")
            write_fasta_records(contigs, assembly_fasta)
            paf_path = os.path.join(class_dir, "assembly_queries.paf")
            align_assembly_to_queries(assembly_fasta, queries_fasta, paf_path, settings.threads)

            blocks = parse_paf(paf_path)
            query_metadata = {row["query_name"]: row for row in query_rows}
            class_data[class_id] = {
                "region": region,
                "contigs": contigs,
                "blocks": blocks,
                "query_metadata": query_metadata,
            }
            for contig_name, sequence in contigs.items():
                all_contigs["{0}|{1}".format(class_id, contig_name)] = sequence

    ### 2. 所有 contig 一次比对到参考基因组（用于还原坐标） ###

    reference_blocks_by_contig = align_contigs_to_reference(
        all_contigs, settings.reference_fasta, output_dirs["work"], settings.threads
    )

    ### 3. 逐类画图、整理位置表 ###

    all_position_rows = []
    summary_rows = []
    contig_map_rows = []
    for class_id in sorted(class_data):
        data = class_data[class_id]
        contigs = data["contigs"]
        blocks = data["blocks"]
        query_metadata = data["query_metadata"]
        blocks_by_contig = collections.defaultdict(list)
        for block in blocks:
            blocks_by_contig[block["target"]].append(block)

        figure_paths = []
        for contig_name in sorted(contigs):
            contig_blocks = blocks_by_contig.get(contig_name, [])
            contig_map_rows.append(
                [
                    class_id,
                    contig_name,
                    len(contigs[contig_name]),
                    os.path.join(
                        settings.assembly_dir, class_id, "overlap_extension", "contigs.fasta"
                    ),
                ]
            )
            if not contig_blocks:
                print(
                    "  {0}: contig {1} 没有侧翼/外源命中，跳过丝带图".format(class_id, contig_name)
                )
                continue
            contig_paf_path = os.path.join(
                output_dirs["work"], "{0}_{1}.paf".format(class_id, contig_name)
            )
            write_contig_paf(contig_blocks, contig_paf_path)
            title = "{0} overlap_extension {1}: 50bp anchors and exog".format(class_id, contig_name)
            figure_path = run_ribbon_tool(
                contig_paf_path, contig_name, class_id, output_dirs, title
            )
            if os.path.exists(figure_path):
                figure_paths.append(figure_path)

        position_rows = build_position_rows(
            class_id, blocks, query_metadata, site_boundaries, reference_blocks_by_contig
        )
        all_position_rows.extend(position_rows)
        summary_rows.append(
            summarize_class(
                class_id, blocks, query_metadata, position_rows, figure_paths, len(contigs)
            )
        )
        print(
            "{0}: contigs={1}，图 {2} 张，位置行 {3} 行".format(
                class_id, len(contigs), len(figure_paths), len(position_rows)
            )
        )

    ### 4. 写表与 Markdown ###

    write_tsv(
        os.path.join(settings.outdir, "class_assembly_positions.tsv"),
        [
            "class_id",
            "assembly_contig",
            "query_name",
            "query_type",
            "query_start",
            "query_end",
            "query_strand",
            "contig_start",
            "contig_end",
            "identity",
            "reference_chrom",
            "reference_start",
            "reference_end",
            "query_source",
            "distance_to_breakpoint_bp",
        ],
        all_position_rows,
        [
            "query_type=flank 的行：reference_* 是通过组装 contig 与参考基因组比对还原的坐标，"
            "不是 50 bp 锚点的提取坐标（query_source 才是提取坐标）。",
            "query_type=exog 的行：reference_* 填 -（外源没有参考基因组坐标）。",
            "distance_to_breakpoint_bp：还原区间到侧翼对应接头边界的距离，0 表示跨过边界。",
            "坐标：query/contig/reference 均为 0-based 半开；query_strand 是查询相对组装 contig 的方向。",
        ],
    )

    write_tsv(
        os.path.join(settings.outdir, "class_position_summary.tsv"),
        [
            "class_id",
            "n_contigs",
            "n_flank_queries",
            "n_flank_queries_hit",
            "n_exog_queries",
            "n_exog_queries_hit",
            "flank_positions",
            "exog_positions",
            "figures",
        ],
        summary_rows,
        [
            "flank_positions：50 bp 断点锚点的最佳位置；外源位置最多列出前 5 处（总数写在方括号前）。",
        ],
    )

    write_tsv(
        os.path.join(settings.outdir, "contig_map.tsv"),
        [
            "class_id",
            "contig",
            "contig_length_bp",
            "assembly_fasta",
        ],
        contig_map_rows,
    )

    write_summary_markdown(summary_rows, os.path.join(settings.outdir, "summary.md"))

    return {
        "outdir": settings.outdir,
        "classes": len(class_data),
        "figures": len(
            [
                name
                for name in os.listdir(os.path.join(settings.outdir, "figures"))
                if name.endswith(".png")
            ]
        ),
        "position_rows": len(all_position_rows),
    }
