#!/usr/bin/env python3
"""生成 ``short coverage`` 的合成 fixture（确定性，无时间戳）。

全部是编造的序列与坐标（chr1 取 2 kb），不含任何真实数据。

用法::

    python tests/fixtures/coverage/make_fixture.py

产物：``reads.bam``（+ ``.bai``），布局见 ``READS``。

``expected.bins.legacy.tsv`` 由**独立脚本**（并入前的 ``coverage-viz --dump``）
对同一个 BAM 生成，是黄金对账的来源::

    python -m coverage_viz tests/fixtures/coverage/reads.bam -r chr1:1-2000 --bins 8 \\
        --dump tests/fixtures/coverage/expected.bins.legacy.tsv \\
        -o /tmp/legacy.png

列是旧工具的 ``start end covered_bases breadth_pct mean_depth``（start/end 为 1-based）。
"""

from __future__ import annotations

import random
from pathlib import Path

import pysam

HERE = Path(__file__).resolve().parent
REFERENCE_NAME = "chr1"
REFERENCE_LENGTH = 2200  # 比统计区间（1-2000）长，末端 read 才不会被截断

#: (read 名, 0-based 起点, CIGAR, MAPQ, flag)；flag 里 256=secondary、1024=duplicate
READS: tuple[tuple[str, int, tuple[tuple[int, int], ...], int, int], ...] = (
    # 平铺覆盖：0..1500 每 60 bp 起一条 100M，深度约 2（跨 bin 边界）
    *((f"tile_{index:02d}", index * 60, ((0, 100),), 60, 0) for index in range(25)),
    # 30D：缺失段不携带碱基，不该计入深度（1651..1680 应保持 0）
    ("deletion_read", 1600, ((0, 50), (2, 30), (0, 50)), 60, 0),
    # 30S70M：soft clip 不占参考位置（1710..1740 应保持 0）
    ("softclip_read", 1740, ((4, 30), (0, 70)), 60, 0),
    # 20I：插入不占参考位置，两段 M 之间的参考位置仍然连续
    ("insertion_read", 1820, ((0, 40), (1, 20), (0, 40)), 60, 0),
    # 低 MAPQ：默认不过滤（--min-mapq 0），计入深度
    ("lowmapq_read", 1920, ((0, 100),), 5, 0),
    # secondary 与 duplicate：默认过滤掉
    ("secondary_read", 1900, ((0, 100),), 60, 256),
    ("duplicate_read", 1905, ((0, 100),), 60, 1024),
)


def _random_sequence(seed: int, length: int) -> str:
    """造一条固定随机种子的序列（保证 fixture 逐字节一致）。"""
    generator = random.Random(seed)
    return "".join(generator.choice("ACGT") for _ in range(length))


def build_bam(path: Path) -> None:
    """写一个含平铺覆盖与四种 CIGAR 边界的合成 BAM（排好序 + 建索引）。

    输入示例: path=tests/fixtures/coverage/reads.bam
    输出示例: 生成 reads.bam(.bai)，覆盖 ``READS`` 里列出的全部记录
    """
    reference = _random_sequence(7, REFERENCE_LENGTH)
    header = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": REFERENCE_NAME, "LN": REFERENCE_LENGTH}],
    }
    raw_path = path.with_suffix(".unsorted.bam")
    with pysam.AlignmentFile(str(raw_path), "wb", header=header) as output:
        for name, reference_start, cigar, mapq, flag in READS:
            read = pysam.AlignedSegment()
            read.query_name = name
            query_length = sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
            read.query_sequence = reference[reference_start : reference_start + query_length]
            read.flag = flag
            read.reference_id = 0
            read.reference_start = reference_start
            read.mapping_quality = mapq
            read.cigartuples = list(cigar)
            read.query_qualities = pysam.qualitystring_to_array("I" * query_length)
            output.write(read)
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
