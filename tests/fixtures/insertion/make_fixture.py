#!/usr/bin/env python3
"""生成 ``short insertion-plot`` 的合成 fixture（确定性，无时间戳）。

全部是编造的序列与坐标（chr1 取 200 kb，断点都在 200 kb 以内），不含任何真实数据。

用法::

    python tests/fixtures/insertion/make_fixture.py

产物：``reads.bam``（+ ``.bai``）；位点布局见 ``SITES``。

``expected.consensus.legacy.tsv`` 由**独立脚本**（并入前的
``insertion_schematic.py --site-table``）对同一个 BAM 生成，是黄金对账的来源；
回归命令（在 alnkit 仓库根目录执行）::

    for site in "" "--junction 100000" "--junction 150000" "--junction 180000"; do ...
    python ../insertion_schematic/insertion_schematic.py --bam tests/fixtures/insertion/reads.bam \\
        --chrom chr1 --start 0 --end 200000 $site --outdir /tmp/legacy \\
        --site-table /tmp/legacy/sites.tsv

四个位点各跑一次后按 ``breakpoint_1based`` 排序拼接（表头只留一行）。
注意旧脚本的 ``breakpoint_1based`` 是“边界 + 1”，alnkit 的口径是“左侧碱基的
1-based 坐标”（数值等于 0-based 边界），所以对账时差 1。
"""

from __future__ import annotations

import random
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
REFERENCE_NAME = "chr1"
REFERENCE_LENGTH = 200_000
INSERT_NAME = "insertA"
INSERT_LENGTH = 1200
ALIGN_LENGTH = 100

#: 位点布局：0-based 断点 -> (该侧, read 数)；L4 是两个断点的多断点图
SITES = {
    100_000: {"left": 20, "right": 20},  # 同断点双侧
    120_000: {"left": 12},  # 多断点（左）
    120_230: {"right": 12},  # 多断点（右）
    150_000: {"left": 15},  # 单侧（左）
    180_000: {"right": 15},  # 单侧（右）
}


def _random_sequence(seed: int, length: int) -> str:
    """造一条固定随机种子的序列（保证每次生成的 fixture 逐字节一致）。

    输入示例: seed=11, length=8  →  输出示例: "GCATTTAC"（只含 ACGT）
    """
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def _make_read(
    name: str,
    reference_start: int,
    clip_side: str,
    clip_sequence: str,
    reference: str,
    *,
    reverse: bool = False,
    secondary: bool = False,
    supplementary: bool = False,
) -> pysam.AlignedSegment:
    """造一条带 soft clip 的单端比对（100M + clip）。"""
    aligned_bases = reference[reference_start : reference_start + ALIGN_LENGTH]
    if clip_side == "left":
        read = pysam.AlignedSegment()
        read.query_sequence = clip_sequence + aligned_bases
        read.cigartuples = [(4, len(clip_sequence)), (0, ALIGN_LENGTH)]
    else:
        read = pysam.AlignedSegment()
        read.query_sequence = aligned_bases + clip_sequence
        read.cigartuples = [(0, ALIGN_LENGTH), (4, len(clip_sequence))]
    read.query_name = name
    read.flag = (16 if reverse else 0) | (256 if secondary else 0) | (2048 if supplementary else 0)
    read.reference_id = 0
    read.reference_start = reference_start
    read.mapping_quality = 60
    read.query_qualities = pysam.qualitystring_to_array("I" * len(read.query_sequence))
    read.set_tag("NM", 0)
    read.set_tag("MD", str(ALIGN_LENGTH))
    return read


def build_bam(path: Path) -> None:
    """写一个含五个位点的合成 BAM（排好序 + 建索引）。

    输入示例: path=tests/fixtures/insertion/reads.bam
    输出示例: 生成 reads.bam(.bai)：双侧 / 多断点 ×2 / 单侧左 / 单侧右，
              另含硬剪贴、负链、secondary、supplementary 各一条用于覆盖边界分支。
    """
    reference = _random_sequence(11, REFERENCE_LENGTH)
    insert_sequence = _random_sequence(29, INSERT_LENGTH)
    left_clip = insert_sequence[:60]  # 左 clip 取自外源 5' 端
    right_clip = insert_sequence[-55:]  # 右 clip 取自外源 3' 端

    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [
            {"SN": REFERENCE_NAME, "LN": REFERENCE_LENGTH},
            {"SN": INSERT_NAME, "LN": INSERT_LENGTH},
        ],
    }
    raw_path = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(raw_path), "wb", header=header) as output:
        for breakpoint, sides in SITES.items():
            for index in range(sides.get("left", 0)):
                output.write(
                    _make_read(
                        f"bp{breakpoint}_left_{index:02d}",
                        breakpoint,
                        "left",
                        left_clip,
                        reference,
                        reverse=(index == 0),  # 负链也要计入
                    )
                )
            for index in range(sides.get("right", 0)):
                output.write(
                    _make_read(
                        f"bp{breakpoint}_right_{index:02d}",
                        breakpoint - ALIGN_LENGTH,
                        "right",
                        right_clip,
                        reference,
                    )
                )
        # 边界分支：硬剪贴不带碱基（应被忽略）、secondary/supplementary 不该重复计数
        hard_clip = pysam.AlignedSegment()
        hard_clip.query_name = "hardclip_read"
        hard_clip.query_sequence = reference[150_000 : 150_000 + ALIGN_LENGTH]
        hard_clip.cigartuples = [(0, ALIGN_LENGTH), (5, 7)]
        hard_clip.flag = 0
        hard_clip.reference_id = 0
        hard_clip.reference_start = 150_000
        hard_clip.mapping_quality = 60
        output.write(hard_clip)
        output.write(
            _make_read("secondary_read", 150_000, "left", left_clip, reference, secondary=True)
        )
        output.write(
            _make_read(
                "supplementary_read", 150_000, "left", left_clip, reference, supplementary=True
            )
        )
    pysam.sort("-o", str(path), str(raw_path))
    pysam.index(str(path))
    raw_path.unlink()


def main() -> int:
    """生成 fixture 并打印产物大小。"""
    build_bam(HERE / "reads.bam")
    bam = HERE / "reads.bam"
    print(f"{bam}  ({bam.stat().st_size // 1024} KB)")
    print(f"{bam}.bai")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
