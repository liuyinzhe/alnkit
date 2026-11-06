"""单条 read 的比对片段轨道图（C2）。

一条 read 一行，横轴是 read 自身坐标；每个比对片段画成一个块：

* 颜色表示正负链：``+`` 蓝、``-`` 橙；外源片段用斜线填充，基因组片段实心；
* 块上方/下方交替标注片段类型、参考序列名、参考坐标、方向；
* ``highlight_ref`` 可标注某个参考断点在每条 read 上的位置。

输入三选一：BAM（``pysam``，反向链换算到正向 read 轴）、minimap2 PAF（需 ``-c``）、
已整理好的片段表（``segments_tsv``）。

产物：``<outdir>/<name>.png|.svg|.pdf`` 与 ``<outdir>/<name>.segments.tsv``。

来源：``read_segment_tracks`` v0.3（git 51bf68f）单文件实现的忠实移植，
产物与原脚本逐字节一致（见 ``tests/test_plot_tracks.py`` 的黄金回归）。
"""

from __future__ import annotations

import collections
import os
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass

import matplotlib
import pysam

from ..core.errors import UsageError
from ..core.segments import ALIASES

matplotlib.use("Agg")  # 无 DISPLAY 环境也必须能出图
import matplotlib.pyplot as plt  # noqa: E402  （backend 必须先设置）
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

STRAND_COLORS = {"+": "#1f77b4", "-": "#ff7f0e"}  # 蓝 = +，橙 = -
EXOG_HATCH = "///"
HIGHLIGHT_COLORS = ["#9467bd", "#17becf", "#e377c2", "#8c564b", "#7f7f7f"]


### 基本工具 ###


def parse_region(text):
    """把 `chrom:pos` 解析成 (chrom, pos)；pos 允许逗号。

    输入示例: parse_region("chr1:3,101,000")
    输出示例: ("chr1", 3101000)
    说明：只支持单点；边界点也可以（会按片段端点规则映射到 read 坐标）。
    """
    if ":" not in text:
        raise UsageError(f"--highlight-ref 需要 chrom:pos 格式，收到: {text}")
    chrom, position_text = text.rsplit(":", 1)
    if not chrom:
        raise UsageError(f"--highlight-ref 缺少序列名，收到: {text}")
    try:
        position = int(position_text.replace(",", ""))
    except ValueError:
        raise UsageError(f"--highlight-ref 的位置不是整数，收到: {text}") from None
    if position < 0:
        raise UsageError(f"--highlight-ref 的位置不能为负，收到: {text}")
    return chrom, position


def infer_segment_type(target_name, exog_names):
    """按目标序列名判断片段类型：外源 / 基因组。

    输入示例: infer_segment_type("exog_a", {"exog_a", "exog_b"})
    输出示例: "exog"
    """
    return "exog" if target_name in exog_names else "genome"


def short_read_label(read_name):
    """给图上左侧行标签取一个短名字。

    输入示例: short_read_label("sample/1001/ccs")
    输出示例: "1001"
    说明：优先取倒数第二段；取不到时返回完整 read 名。
    """
    parts = read_name.split("/")
    if len(parts) >= 2 and parts[-2]:
        return parts[-2]
    return read_name


def read_labels(read_names, full_read_names):
    """生成左侧行标签，撞名时自动加 #1/#2 区分。

    输入示例: read_labels(["a/7/ccs", "b/7/ccs"], False)
    输出示例: ["7#1", "7#2"]
    """
    raw_labels = [
        read_name if full_read_names else short_read_label(read_name) for read_name in read_names
    ]
    totals = collections.Counter(raw_labels)
    seen = collections.Counter()
    labels = []
    for label in raw_labels:
        if totals[label] > 1:
            seen[label] += 1
            labels.append(f"{label}#{seen[label]}")
        else:
            labels.append(label)
    return labels


def ref_pos_to_read(segment, chrom, position):
    """把参考坐标/边界换算成该片段在 read 轴上的坐标。

    输入示例:
      segment={"target_name": "chr1", "target_start": 100, "target_end": 200,
               "read_start": 10, "read_end": 110, "strand": "-"}
      ref_pos_to_read(segment, "chr1", 200) -> 10   # 片段右端边界，反向链
    输出示例: 10 / 59 / 110 / None
    说明：
      * 片段内部：+ 链 read_start + (pos - target_start)；
        - 链 read_end - 1 - (pos - target_start)；
      * pos == target_start：+ 链取 read_start，- 链取 read_end（左边界）；
      * pos == target_end  ：+ 链取 read_end，  - 链取 read_start（右边界）；
      * 不属于该片段返回 None。
    """
    if segment["target_name"] != chrom:
        return None
    target_start = segment["target_start"]
    target_end = segment["target_end"]
    read_start = segment["read_start"]
    read_end = segment["read_end"]
    if not (target_start <= position <= target_end):
        return None
    if position == target_start:
        return read_start if segment["strand"] == "+" else read_end
    if position == target_end:
        return read_end if segment["strand"] == "+" else read_start
    if segment["strand"] == "+":
        return read_start + (position - target_start)
    return read_end - 1 - (position - target_start)


### 读取三种输入 ###


def cigar_to_read_span(cigartuples, read_length, is_reverse):
    """把 CIGAR 的比对段换算成正向 read 轴上的 (read_start, read_end)。

    输入示例: cigar_to_read_span(((5, 5), (0, 100)), 105, True)   # 5H100M，反向
    输出示例: (0, 100)；无比对段或坐标越界时返回 None
    说明：
      * cigartuples 是 pysam 的 (operation, length) 列表，H/S 在两端；
      * infer_read_length() 含硬剪切，而 CIGAR 里的偏移是相对存储 SEQ（不含 H）的，
        所以先补前导硬剪切再按方向镜像；纯软剪切（无 H）时与旧公式逐字节一致。
    """
    if not cigartuples:
        return None
    query_ops = (pysam.CMATCH, pysam.CINS, pysam.CSOFT_CLIP, pysam.CEQUAL, pysam.CDIFF)
    index = 0
    leading_hard_clip = 0
    while index < len(cigartuples) and cigartuples[index][0] == pysam.CHARD_CLIP:
        leading_hard_clip += cigartuples[index][1]
        index += 1
    leading_soft_clip = 0
    while index < len(cigartuples) and cigartuples[index][0] == pysam.CSOFT_CLIP:
        leading_soft_clip += cigartuples[index][1]
        index += 1
    index = len(cigartuples) - 1
    while index >= 0 and cigartuples[index][0] == pysam.CHARD_CLIP:
        index -= 1
    trailing_soft_clip = 0
    while index >= 0 and cigartuples[index][0] == pysam.CSOFT_CLIP:
        trailing_soft_clip += cigartuples[index][1]
        index -= 1
    stored_length = sum(length for operation, length in cigartuples if operation in query_ops)
    aligned_length = stored_length - leading_soft_clip - trailing_soft_clip
    if aligned_length <= 0:
        return None
    query_start = leading_soft_clip
    query_end = leading_soft_clip + aligned_length
    if is_reverse:
        read_start = read_length - leading_hard_clip - query_end
        read_end = read_length - leading_hard_clip - query_start
    else:
        read_start = leading_hard_clip + query_start
        read_end = leading_hard_clip + query_end
    if read_start < 0 or read_end <= read_start:
        return None
    return read_start, read_end


def cigar_ref_to_read(cigartuples, target_start, read_length, is_reverse, position):
    """参考坐标 -> 正向 read 坐标（沿 CIGAR 走，O(操作数) 时间、O(1) 内存）。

    输入示例: cigar_ref_to_read(((0, 100), (2, 100), (0, 100)), 1000, 200, False, 1200)
    输出示例: 100   # 前 100 bp 对在 read 0-99，中间 100 bp 缺失不占 read
    说明：
      * 只映射落在 M/=/X 碱基上的位置；落在 D/N 缺失里或片段外返回 None；
      * 与 get_aligned_pairs() 结果一致，但不物化任何列表，逐条记录调用也没有内存压力；
      * 边界（pos == target_start / target_end）仍走 ref_pos_to_read()，规则不变。
    """
    if not cigartuples:
        return None
    reference_position = target_start
    query_position = 0
    index = 0
    while index < len(cigartuples) and cigartuples[index][0] == pysam.CHARD_CLIP:
        index += 1
    # 硬剪切不在存储 SEQ 里：只影响正向轴偏移，不推进 SEQ 下标
    leading_hard_clip = sum(length for operation, length in cigartuples[:index])
    for operation, length in cigartuples[index:]:
        if operation in (pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF):
            if reference_position <= position < reference_position + length:
                query_position += position - reference_position
                if is_reverse:
                    return read_length - 1 - leading_hard_clip - query_position
                return leading_hard_clip + query_position
            reference_position += length
            query_position += length
        elif operation == pysam.CINS:
            query_position += length
        elif operation in (pysam.CDEL, pysam.CREF_SKIP):
            if reference_position <= position < reference_position + length:
                return None
            reference_position += length
        elif operation == pysam.CSOFT_CLIP:
            query_position += length
    return None


def add_segment(reads, segment, read_length):
    """把一个片段加进 {read_name: {read_length, segments}}，并保持 read 长度一致。

    输入示例: add_segment({}, {"read_name": "r1", ...}, 19000)
    输出示例: {"r1": {"read_length": 19000, "segments": [{...}]}}
    """
    record = reads.setdefault(segment["read_name"], {"read_length": read_length, "segments": []})
    record["read_length"] = max(record["read_length"], read_length, segment["read_end"])
    record["segments"].append(segment)


#: 片段表必需列（注册表规范名；历史表头由 core.segments.ALIASES 映射过来）
_TRACK_COLUMNS = (
    "read_name",
    "read_length",
    "chrom",
    "ref_start_0based",
    "ref_end_0based",
    "read_start_0based",
    "read_end_0based",
    "strand",
)


def load_segments_tsv(path, exog_names):
    """读片段表；列见模块文档 / README。

    输入示例: load_segments_tsv("example_segments.tsv", {"exog_a"})
    输出示例: {"r1": {"read_length": 19000, "segments": [...]}}
    说明：用 utf-8-sig 打开，Excel 导出的 BOM 不会污染列名；列名容忍历史写法，
          所以 `alnkit seg table` 写出的规范表（chrom / ref_start_0based …）也能直接画。
    """
    reads = {}
    header = None
    with open(path, encoding="utf-8-sig") as table:
        for line_number, line in enumerate(table, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if header is None:
                header = fields
                index = {name: i for i, name in enumerate(header)}
                # 历史表头（target_* / read_start）与注册表列名都能读，见 core.segments.ALIASES
                resolved = {}
                missing = []
                for canonical in (*_TRACK_COLUMNS, "segment_role", "mapq"):
                    for alias in ALIASES[canonical]:
                        if alias in index:
                            resolved[canonical] = index[alias]
                            break
                    else:
                        if canonical in _TRACK_COLUMNS:
                            missing.append(canonical)
                if missing:
                    raise UsageError("片段表缺少列: {0}".format(", ".join(sorted(missing))))
                continue
            if len(fields) != len(header):
                raise UsageError(
                    f"片段表第 {line_number} 行有 {len(fields)} 列，应为 {len(header)} 列"
                )
            row = {canonical: fields[position] for canonical, position in resolved.items()}
            read_name = row["read_name"]
            read_length = int(row["read_length"])
            segment_type = row.get("segment_role") or ""
            if segment_type not in ("genome", "exog"):
                segment_type = infer_segment_type(row["chrom"], exog_names)
            strand = row["strand"]
            if strand not in ("+", "-"):
                raise UsageError(f"片段表第 {line_number} 行 strand 只能是 + 或 -，收到: {strand}")
            segment = {
                "read_name": read_name,
                "segment_type": segment_type,
                "target_name": row["chrom"],
                "target_start": int(row["ref_start_0based"]),
                "target_end": int(row["ref_end_0based"]),
                "read_start": int(row["read_start_0based"]),
                "read_end": int(row["read_end_0based"]),
                "strand": strand,
                "mapq": int(row.get("mapq", 0) or 0),
            }
            if segment["read_end"] > read_length:
                print(
                    "警告：片段表第 {0} 行 read_end({1}) > read_length({2})，已按 read_end 调整".format(
                        line_number, segment["read_end"], read_length
                    ),
                    file=sys.stderr,
                )
                read_length = segment["read_end"]
            add_segment(reads, segment, read_length)
    return reads


def load_segments_bam(bam_path, exog_names, read_names, min_mapq, include_secondary):
    """从 BAM 提取片段；只保留 read_names（None = 全部）。

    输入示例: load_segments_bam("x.bam", {"exog_a"}, {"r1"}, 0, False)
    输出示例: {"r1": {"read_length": 19000, "segments": [...]}}
    说明：
      * 跳过 unmapped；默认跳过 secondary，保留 supplementary；
      * read 坐标统一到正向 read 轴，换算交给 cigar_to_read_span()（补前导硬剪切再镜像）。
    """
    reads = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam:
            if read.is_unmapped:
                continue
            if read.is_secondary and not include_secondary:
                continue
            if read_names is not None and read.query_name not in read_names:
                continue
            if read.mapping_quality < min_mapq:
                continue
            if read.reference_name is None or read.reference_end is None:
                continue
            read_length = read.infer_read_length() or read.query_length or 0
            if read_length <= 0:
                continue
            span = cigar_to_read_span(read.cigartuples, read_length, read.is_reverse)
            if span is None:
                continue
            read_start, read_end = span
            segment = {
                "read_name": read.query_name,
                "segment_type": infer_segment_type(read.reference_name, exog_names),
                "target_name": read.reference_name,
                "target_start": read.reference_start,
                "target_end": read.reference_end,
                "read_start": read_start,
                "read_end": read_end,
                "strand": "-" if read.is_reverse else "+",
                "mapq": read.mapping_quality,
                "read_length": read_length,
                "cigar": read.cigartuples,
            }
            add_segment(reads, segment, read_length)
    return reads


def bam_has_index(bam_path):
    """判断 BAM 是否带索引（samtools 按名取记录才快）。"""
    base = bam_path[:-4] if bam_path.lower().endswith(".bam") else bam_path
    return any(
        os.path.exists(path)
        for path in (bam_path + ".bai", bam_path + ".csi", base + ".bai", base + ".csi")
    )


def load_segments_bam_subset(bam_path, exog_names, read_names, min_mapq, include_secondary):
    """用 `samtools view -N` 只抽出目标 read，再交给 load_segments_bam()。

    输入示例: load_segments_bam_subset("x.bam", {"exog_a"}, {"r1"}, 0, False)
    输出示例: {"r1": {...}}；samtools 不可用或调用失败时返回 None（调用方回退 scan）
    说明：结果与全量 scan 完全一致（只影响速度）；没索引也能跑，但有索引才真正省时间。
    """
    if not read_names:
        return None
    samtools = shutil.which("samtools")
    if samtools is None:
        return None
    with tempfile.TemporaryDirectory() as directory:
        names_path = os.path.join(directory, "reads.txt")
        with open(names_path, "w") as names:
            names.write("".join(f"{name}\n" for name in sorted(read_names)))
        subset_path = os.path.join(directory, "subset.bam")
        completed = subprocess.run(
            [samtools, "view", "-b", "-N", names_path, "-o", subset_path, bam_path],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        if completed.returncode != 0 or not os.path.exists(subset_path):
            return None
        return load_segments_bam(subset_path, exog_names, None, min_mapq, include_secondary)


def load_bam_input(settings, exog_names, read_names):
    """按 --bam-mode 选 scan / subset；auto 时有索引且给了 read 名单就 subset。"""
    mode = settings.bam_mode
    if mode == "auto":
        mode = "subset" if read_names and bam_has_index(settings.bam) else "scan"
    if mode == "subset":
        if not read_names:  # 先分清“没给名单”和“没有 samtools”，否则报错会指错方向
            raise UsageError(
                "--bam-mode subset 需要 --read 或 --reads-file（subset 只取指定 read）"
            )
        reads = load_segments_bam_subset(
            settings.bam, exog_names, read_names, settings.min_mapq, settings.include_secondary
        )
        if reads is not None:
            return reads
        if settings.bam_mode == "subset":
            raise UsageError("--bam-mode subset 需要 samtools（PATH 里找不到或调用失败）")
    return load_segments_bam(
        settings.bam, exog_names, read_names, settings.min_mapq, settings.include_secondary
    )


def load_segments_paf(paf_path, exog_names, read_names, min_mapq, include_secondary):
    """从 minimap2 PAF 提取片段（需要 -c）。

    输入示例: load_segments_paf("aln.paf", {"exog_a"}, {"r1"}, 0, False)
    输出示例: {"r1": {"read_length": 19000, "segments": [...]}}
    说明：
      * PAF 的 query 坐标已经是正向 read 轴，不需要像 BAM 那样镜像；
      * PAF 的 secondary 写在可选标签 `tp:A:S` 里，默认跳过（与 BAM 一致）；
        `--include-secondary` 时保留。
    """
    reads = {}
    with open(paf_path) as paf:
        for line in paf:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12 or fields[5] == "*":
                continue
            if not include_secondary and "tp:A:S" in fields[12:]:
                continue
            read_name = fields[0]
            if read_names is not None and read_name not in read_names:
                continue
            mapq = int(fields[11])
            if mapq < min_mapq:
                continue
            read_length = int(fields[1])
            read_start = int(fields[2])
            read_end = int(fields[3])
            target_start = int(fields[7])
            target_end = int(fields[8])
            if read_end <= read_start or target_end <= target_start:
                continue
            segment = {
                "read_name": read_name,
                "segment_type": infer_segment_type(fields[5], exog_names),
                "target_name": fields[5],
                "target_start": target_start,
                "target_end": target_end,
                "read_start": read_start,
                "read_end": read_end,
                "strand": fields[4],
                "mapq": mapq,
            }
            add_segment(reads, segment, read_length)
    return reads


### 输出 ###


def write_segments_tsv(reads, output_path):
    """把实际用于绘图的片段写成 tsv，便于核对。

    输入示例: write_segments_tsv(reads, "out/tracks.segments.tsv")
    输出示例: 一行一个片段，列为 read_name/read_length/segment_type/...
    """
    columns = [
        "read_name",
        "read_length",
        "segment_type",
        "target_name",
        "target_start",
        "target_end",
        "read_start",
        "read_end",
        "strand",
        "mapq",
    ]
    with open(output_path, "w") as output:
        output.write("\t".join(columns) + "\n")
        for read_name in sorted(reads):
            record = reads[read_name]
            for segment in sorted(record["segments"], key=lambda s: s["read_start"]):
                output.write(
                    "\t".join(
                        [
                            read_name,
                            str(record["read_length"]),
                            segment["segment_type"],
                            segment["target_name"],
                            str(segment["target_start"]),
                            str(segment["target_end"]),
                            str(segment["read_start"]),
                            str(segment["read_end"]),
                            segment["strand"],
                            str(segment["mapq"]),
                        ]
                    )
                    + "\n"
                )


def resolve_highlights(record, highlights, exact_highlights=False):
    """算出某条 read 上每组断点的 read 坐标。

    输入示例: resolve_highlights(record, [("chr1", 1250, "#9467bd")], exact=True)
    输出示例: [("#9467bd", [])]   # 该位置落在 500 bp 缺失里，没有对应 read 碱基
    说明：
      * 默认按片段两端线性插值（ref_pos_to_read）；exact=True 时片段内部按 CIGAR 精确映射；
      * pos == target_start / target_end 的边界始终按 ref_pos_to_read() 的约定取端点；
      * 返回的坐标已排序去重。
    """
    resolved = []
    for chrom, position, color in highlights:
        positions = set()
        for segment in record["segments"]:
            if segment["target_name"] != chrom:
                continue
            # 边界优先：pos == target_start / target_end 永远按文档规则取端点（与方向无关）
            if position in (segment["target_start"], segment["target_end"]):
                read_position = ref_pos_to_read(segment, chrom, position)
            elif exact_highlights and segment.get("cigar"):
                read_position = cigar_ref_to_read(
                    segment["cigar"],
                    segment["target_start"],
                    segment.get("read_length", record["read_length"]),
                    segment["strand"] == "-",
                    position,
                )
            else:
                read_position = ref_pos_to_read(segment, chrom, position)
            # 片段 ref 跨度远大于 read 跨度（大 N/D）时线性插值会算出轴外坐标，直接丢掉
            if read_position is not None and 0 <= read_position <= record["read_length"]:
                positions.add(read_position)
        resolved.append((color, sorted(positions)))
    return resolved


def _overlaps(first, second, margin=0.5):
    """两个显示坐标 bbox 是否重叠（留一点空隙，单位像素）。"""
    return not (
        first.x1 + margin <= second.x0
        or second.x1 + margin <= first.x0
        or first.y1 + margin <= second.y0
        or second.y1 + margin <= first.y0
    )


def _add_segment_label(axes, x, y, text, offset):
    """在 (x, y+offset) 放一个片段标签，返回 Text 对象。"""
    return axes.text(
        x,
        y + offset,
        text,
        ha="center",
        va="bottom" if offset > 0 else "top",
        fontsize=5.8,
        color="#333333",
        bbox=dict(facecolor="white", alpha=0.72, edgecolor="none", pad=0.3),
        zorder=5,
    )


def place_segment_labels(axes, entries):
    """按候选高度逐条放标签，返回“四个位置都撞”的条数。

    entries: [(x, y, 文本, 首选偏移)]，首选偏移为 +0.24 / -0.24（上下交替）。
    说明：用渲染器的文字 bbox 做碰撞检测；首选位置空着就保持老样子，
          撞了就换另一侧、再撞就换更远的两个 lane；四层都撞时放回首选位置并计数，
          调用方（--label-mode auto）看到计数 > 0 会改用侧边注释表。
    """
    renderer = axes.figure.canvas.get_renderer()
    placed = []
    unplaced = 0
    for x, y, text, preferred in entries:
        offsets = (preferred, -preferred, 0.55, -0.55)
        for offset in offsets:
            label = _add_segment_label(axes, x, y, text, offset)
            box = label.get_window_extent(renderer)
            if not any(_overlaps(box, other) for other in placed):
                placed.append(box)
                break
            label.remove()
        else:
            label = _add_segment_label(axes, x, y, text, preferred)
            placed.append(label.get_window_extent(renderer))
            unplaced += 1
    return unplaced


def side_table_lines(reads, full_read_names):
    """生成侧边注释表的文本行（纯函数，便于测试）。

    输入示例: side_table_lines(reads, False)
    输出示例: ["1002  (read 21,544 bp)", "  1. genome chr1 3,000,000-40,000,000 (-)  read 0-13,698", ...]
    """
    read_names = sorted(reads)
    labels = read_labels(read_names, full_read_names)
    lines = []
    for label, read_name in zip(labels, read_names):
        record = reads[read_name]
        shown = label if len(label) <= 40 else label[:40] + "..."
        lines.append("{0}  (read {1:,} bp)".format(shown, record["read_length"]))
        segments = sorted(record["segments"], key=lambda s: (s["read_start"], s["target_name"]))
        for index, segment in enumerate(segments, 1):
            lines.append(
                "  {0}. {1} {2} {3:,}-{4:,} ({5})  read {6:,}-{7:,}".format(
                    index,
                    segment["segment_type"],
                    segment["target_name"],
                    segment["target_start"],
                    segment["target_end"],
                    segment["strand"],
                    segment["read_start"],
                    segment["read_end"],
                )
            )
    return lines


def draw_side_table(figure, lines, font_size=5.5):
    """在图的右侧画注释表（必要时分两列，仍放不下就截断并提示看 .segments.tsv）。

    输入示例: draw_side_table(figure, side_table_lines(reads, False))
    输出示例: 图上多出若干行等宽文字，返回实际画出的行数
    """
    if not lines:
        return 0
    width_in, height_in = figure.get_size_inches()
    left, right, top, bottom = 0.55, 0.99, 0.96, 0.02
    line_height = font_size * 1.35 / 72 / height_in
    rows_per_column = max(1, int((top - bottom) / line_height))
    columns = 1 if len(lines) <= rows_per_column else 2
    column_width_in = (right - left) * width_in / columns
    max_chars = max(20, int(column_width_in / (font_size * 0.6 / 72)))
    capacity = rows_per_column * columns
    if len(lines) > capacity:
        lines = lines[: capacity - 1] + [
            f"... {len(lines) - capacity + 1} more rows: see .segments.tsv"
        ]
    for index, line in enumerate(lines):
        column = index // rows_per_column
        row = index % rows_per_column
        figure.text(
            left + column * (right - left) / columns,
            top - row * line_height,
            line[:max_chars],
            fontsize=font_size,
            family="monospace",
            va="top",
            ha="left",
            color="#333333",
        )
    return len(lines)


def build_x_map(reads, max_axis_bp):
    """超长 read 的断轴映射：把内容之间的大空隙折叠，内容区间内保持真实比例。

    输入示例: build_x_map(reads, 50000)   # 超过 50 kb 才折叠，空隙阈值 = max(200, 50000//25)
    输出示例: (xmap, breaks, pieces, total)；不触发时 xmap 是恒等函数、breaks 为空
    说明：xmap 是分段线性函数，只让显示坐标变短，read 坐标本身不变。
    """
    max_read_length = max(record["read_length"] for record in reads.values())
    if max_axis_bp is None or max_read_length <= max_axis_bp:
        return (lambda value: value), [], [], max_read_length
    gap_threshold = max(200, max_axis_bp // 25)
    intervals = []
    for record in reads.values():
        length = record["read_length"]
        intervals.append((0, min(1, length)))
        intervals.append((max(0, length - 1), length))
        for segment in record["segments"]:
            intervals.append((segment["read_start"], segment["read_end"]))
    intervals.sort()
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1] + gap_threshold:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    collapse = 400.0
    pieces = []
    breaks = []
    offset = 0.0
    for start, end in merged:
        if pieces:
            breaks.append(offset)
            offset += collapse
        pieces.append((start, end, offset))
        offset += end - start
    total = offset

    def xmap(value):
        for start, end, display in pieces:
            if value <= end:
                return display + max(0.0, value - start)
        return total

    return xmap, breaks, pieces, total


def _render_tracks(
    reads, title, highlights, full_read_names, fig_width, exact_highlights, label_mode, max_axis_bp
):
    """搭一张轨道图，返回 (figure, 放不下的标签数)。"""
    read_names = sorted(reads)
    labels = read_labels(read_names, full_read_names)
    xmap, breaks, pieces, total_x = build_x_map(reads, max_axis_bp)
    row_step = 1.6
    figure_height = max(3.5, 0.95 * len(read_names) + 2.0)
    figure, axes = plt.subplots(figsize=(fig_width, figure_height))
    axes.set_xlim(0, total_x * 1.02)
    axes.set_ylim(-0.9, (len(read_names) - 1) * row_step + 0.9)
    axes.set_title(title, fontsize=12, fontweight="bold")
    if pieces:
        axes.set_xlabel("Read coordinate (0-based, bp; gaps collapsed)")
        positions, tick_labels = [], []
        for piece_index, (start, end, _display) in enumerate(pieces):
            for value_index, value in enumerate((start, end)):
                position = xmap(value)
                if positions and abs(position - positions[-1]) < 0.5:
                    continue
                positions.append(position)
                # 折叠点两侧的刻度挤在一起：第二列换行显示，避免标签互相压住
                prefix = "\n" if piece_index and value_index == 0 else ""
                tick_labels.append(f"{prefix}{int(value):,}")
        axes.set_xticks(positions)
        axes.set_xticklabels(tick_labels, fontsize=7)
        for break_position in breaks:
            axes.axvline(break_position, color="#999999", linestyle=":", linewidth=0.8, zorder=2)
            axes.text(
                break_position,
                axes.get_ylim()[0] + 0.05,
                "//",
                ha="center",
                va="bottom",
                fontsize=7,
                color="#555555",
                zorder=6,
            )
    else:
        axes.set_xlabel("Read coordinate (0-based, bp)")

    y_positions = []
    y_labels = []
    label_entries = []
    for row_index, read_name in enumerate(read_names):
        y = (len(read_names) - 1 - row_index) * row_step
        y_positions.append(y)
        y_labels.append(labels[row_index])
        record = reads[read_name]
        axes.plot(
            [xmap(0), xmap(record["read_length"])], [y, y], color="#BBBBBB", linewidth=1.2, zorder=1
        )
        segments = sorted(record["segments"], key=lambda s: (s["read_start"], s["target_name"]))
        for segment_index, segment in enumerate(segments):
            color = STRAND_COLORS.get(segment["strand"], "#7f7f7f")
            hatch = EXOG_HATCH if segment["segment_type"] == "exog" else None
            start_x = xmap(segment["read_start"])
            width = xmap(segment["read_end"]) - start_x
            axes.add_patch(
                Rectangle(
                    (start_x, y - 0.16),
                    width,
                    0.32,
                    facecolor=color,
                    alpha=0.85,
                    edgecolor="#333333",
                    linewidth=0.6,
                    hatch=hatch,
                    zorder=3,
                )
            )
            label_entries.append(
                (
                    start_x + width / 2,
                    y,
                    "{0} {1}\n{2:,}-{3:,} ({4})".format(
                        segment["segment_type"],
                        segment["target_name"],
                        segment["target_start"],
                        segment["target_end"],
                        segment["strand"],
                    ),
                    0.24 if segment_index % 2 == 0 else -0.24,
                )
            )
        # 参考断点在每条 read 上的位置
        for color, positions in resolve_highlights(record, highlights, exact_highlights):
            for read_position in positions:
                x = xmap(read_position)
                axes.plot(
                    [x, x],
                    [y - 0.24, y + 0.24],
                    color=color,
                    linestyle="--",
                    linewidth=1.0,
                    zorder=4,
                )

    axes.set_yticks(y_positions)
    axes.set_yticklabels(y_labels, fontsize=9)
    axes.tick_params(axis="y", length=0)
    axes.spines[["top", "right", "left"]].set_visible(False)

    handles = [
        Patch(facecolor=STRAND_COLORS["+"], alpha=0.85, label="strand +"),
        Patch(facecolor=STRAND_COLORS["-"], alpha=0.85, label="strand -"),
    ]
    if any(
        segment["segment_type"] == "exog"
        for record in reads.values()
        for segment in record["segments"]
    ):
        handles.append(
            Patch(
                facecolor="white",
                edgecolor="#333333",
                hatch=EXOG_HATCH,
                label="exogenous (hatched)",
            )
        )
    for chrom, position, color in highlights:
        handles.append(Line2D([0], [0], color=color, linestyle="--", label=f"{chrom}:{position:,}"))
    # 图例固定在坐标轴下方约 0.6 in：若按轴高比例，行数一多就会被甩出画布
    legend_gap = 0.6 / figure_height
    axes.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -legend_gap),
        ncol=min(5, len(handles)),
        fontsize=7,
        frameon=False,
    )

    unplaced = 0
    if label_mode == "side":
        # 右侧注释表：不占坐标轴位置，图例仍留在轴下方（tight_layout 会预留空间）
        figure.tight_layout(rect=(0, 0, 0.52, 1))
        draw_side_table(figure, side_table_lines(reads, full_read_names))
    else:
        figure.tight_layout()
        unplaced = place_segment_labels(axes, label_entries)
    return figure, unplaced


def plot_tracks(
    reads,
    output_path,
    title,
    highlights,
    full_read_names,
    dpi=150,
    fig_width=14.0,
    exact_highlights=False,
    label_mode="auto",
    max_axis_bp=None,
):
    """画“一条 read 一行”的片段轨道图。

    输入示例: plot_tracks(reads, "out/tracks.png", "Read tracks", [("chr1", 3101000)],
                        fig_width=6.5, label_mode="auto", max_axis_bp=50000)
    输出示例: (png_path, svg_path, pdf_path)
    说明：
      * 颜色只表示链：+ 蓝、- 橙；外源片段斜线填充，基因组实心；
      * highlights 用来在每条 read 上标参考断点的位置（短竖线 + 图例）；
      * exact_highlights=True 时片段内部位置按 CIGAR 精确映射（块本身不变）；
      * label_mode: inline=块上下标签；side=底部注释表；auto=先 inline，放不下自动切 side；
      * max_axis_bp 给定时，超过该长度的 read 会把内容之间的大空隙折叠（画 // 标记）；
      * fig_width 同时作用于三种输出格式（PDF 排版时用 6.5 / 3.5 这类物理宽）。
    """
    read_names = sorted(reads)
    if not read_names:
        raise UsageError("没有可画的 reads")
    if not (fig_width > 0 and dpi > 0):  # NaN 也走这里
        raise UsageError("--fig-width / --dpi 必须为正数")
    if label_mode not in ("auto", "inline", "side"):
        raise UsageError("--label-mode 只能是 auto / inline / side")
    if max_axis_bp is not None and max_axis_bp <= 0:
        raise UsageError("--max-axis-bp 必须为正数")
    normalized_highlights = []
    for index, item in enumerate(highlights):
        chrom, position = item[0], item[1]
        color = item[2] if len(item) > 2 else HIGHLIGHT_COLORS[index % len(HIGHLIGHT_COLORS)]
        normalized_highlights.append((chrom, position, color))
    highlights = normalized_highlights

    mode = "side" if label_mode == "side" else "inline"
    figure, unplaced = _render_tracks(
        reads, title, highlights, full_read_names, fig_width, exact_highlights, mode, max_axis_bp
    )
    if label_mode == "auto" and unplaced:
        plt.close(figure)
        figure, unplaced = _render_tracks(
            reads,
            title,
            highlights,
            full_read_names,
            fig_width,
            exact_highlights,
            "side",
            max_axis_bp,
        )
    elif label_mode == "inline" and unplaced:
        print(
            f"提示：{unplaced} 个片段标签放不下（inline 模式），可改用 --label-mode auto/side",
            file=sys.stderr,
        )

    # 输出可复现：SVG 的 id 哈希、PDF 的时间戳默认都随每次运行变化
    matplotlib.rcParams["svg.hashsalt"] = "read_segment_tracks"
    figure.savefig(output_path, dpi=dpi)
    svg_path = output_path.rsplit(".", 1)[0] + ".svg"
    pdf_path = output_path.rsplit(".", 1)[0] + ".pdf"
    figure.savefig(svg_path, metadata={"Date": None})
    figure.savefig(pdf_path, metadata={"CreationDate": None, "ModDate": None})
    plt.close(figure)
    return output_path, svg_path, pdf_path


### 主流程 ###


@dataclass(slots=True)
class TrackSettings:
    """``alnkit plot tracks`` 的全部入参（与旧脚本参数一一对应）。"""

    outdir: str = "."
    name: str = "read_tracks"
    bam: str | None = None
    paf: str | None = None
    segments_tsv: str | None = None
    read: tuple[str, ...] = ()
    reads_file: str | None = None
    exog: tuple[str, ...] = ()
    bam_mode: str = "auto"
    min_mapq: int = 0
    include_secondary: bool = False
    highlight_ref: tuple[str, ...] = ()
    title: str = "Read-axis alignment segments"
    full_read_names: bool = False
    dpi: int = 150
    fig_width: float = 14.0
    exact_highlight: bool = False
    label_mode: str = "auto"
    max_axis_bp: int | None = None


def flatten_names(values):
    names = set()
    for value in values:
        for part in value.split(","):
            if part.strip():
                names.add(part.strip())
    return names


def load_read_names(settings):
    read_names = set(settings.read)
    if settings.reads_file:
        with open(settings.reads_file) as read_file:
            for line in read_file:
                name = line.strip()
                if name and not name.startswith("#"):
                    read_names.add(name)
    return read_names or None


def run(settings: TrackSettings) -> dict[str, object]:
    sources = [settings.bam, settings.paf, settings.segments_tsv]
    if sum(1 for source in sources if source) != 1:
        raise UsageError("必须且只能指定 bam / paf / segments_tsv 之一")
    exog_names = flatten_names(settings.exog)
    read_names = load_read_names(settings)
    if settings.exact_highlight and not settings.bam:
        print(
            "提示：--exact-highlight 只对 --bam 输入生效（PAF/TSV 没有 CIGAR，仍按线性插值）",
            file=sys.stderr,
        )
    if settings.bam:
        reads = load_bam_input(settings, exog_names, read_names)
    elif settings.paf:
        reads = load_segments_paf(
            settings.paf, exog_names, read_names, settings.min_mapq, settings.include_secondary
        )
    else:
        reads = load_segments_tsv(settings.segments_tsv, exog_names)

    if settings.min_mapq > 0:
        for record in reads.values():
            record["segments"] = [
                segment for segment in record["segments"] if segment["mapq"] >= settings.min_mapq
            ]
        reads = {name: record for name, record in reads.items() if record["segments"]}

    if read_names is not None:
        reads = {name: record for name, record in reads.items() if name in read_names}
        missing = sorted(read_names - set(reads))
        if missing:
            print(
                "警告：以下 read 没有可画的比对段：{0}".format(", ".join(missing[:10])),
                file=sys.stderr,
            )
    if not reads:
        raise UsageError("没有可画的 reads（检查 --read/--reads-file、--exog、--min-mapq）")

    highlights = []
    for index, text in enumerate(settings.highlight_ref):
        chrom, position = parse_region(text)
        highlights.append((chrom, position, HIGHLIGHT_COLORS[index % len(HIGHLIGHT_COLORS)]))

    os.makedirs(settings.outdir, exist_ok=True)
    output_base = os.path.join(settings.outdir, settings.name)
    segments_path = output_base + ".segments.tsv"
    write_segments_tsv(reads, segments_path)
    png_path, svg_path, pdf_path = plot_tracks(
        reads,
        output_base + ".png",
        settings.title,
        highlights,
        settings.full_read_names,
        settings.dpi,
        settings.fig_width,
        settings.exact_highlight,
        settings.label_mode,
        settings.max_axis_bp,
    )
    segment_count = sum(len(record["segments"]) for record in reads.values())
    return {
        "reads": len(reads),
        "segments": segment_count,
        "png": png_path,
        "svg": svg_path,
        "pdf": pdf_path,
        "segments_tsv": segments_path,
    }
