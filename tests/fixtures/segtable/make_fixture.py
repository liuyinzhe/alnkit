#!/usr/bin/env python3
"""生成 seg table 的合成 fixture（确定性，无时间戳）。

用法：python tests/fixtures/segtable/make_fixture.py
产物：ref.fa、ref.fa.fai、reads.bam（+ .bai）；
     expected.legacy.tsv 由旧脚本 alignment_table.py 生成，见 tests/test_seg_table.py 注释。
"""
import pysam

HERE = __import__("pathlib").Path(__file__).resolve().parent
REF = {"chr1": "ACGT" * 500, "exog_a": "TTTT" * 100}

(HERE / "ref.fa").write_text("".join(f">{name}\n{seq}\n" for name, seq in REF.items()))
pysam.faidx(str(HERE / "ref.fa"))

header = {
    "HD": {"VN": "1.6", "SO": "coordinate"},
    "SQ": [{"SN": n, "LN": len(s)} for n, s in REF.items()],
}
records = [
    # read, ref, start, cigar, reverse, secondary, supplementary
    ("readA", "chr1", 100, [(0, 100), (4, 50)], False, False, False),
    ("readA", "exog_a", 0, [(4, 100), (0, 50)], False, False, False),
    ("readB", "chr1", 400, [(4, 60), (0, 140)], True, False, False),
    ("readB", "exog_a", 100, [(4, 60), (0, 140)], True, False, False),
    ("readC", "chr1", 800, [(0, 120), (3, 30), (0, 50), (4, 20)], False, False, False),
    ("readC", "exog_a", 0, [(4, 20), (0, 180)], False, False, False),
    ("readC", "chr1", 900, [(4, 20), (0, 180)], False, True, False),
    ("readD", "chr1", 1200, [(0, 200)], False, False, True),
]
with pysam.AlignmentFile(str(HERE / "reads.bam"), "wb", header=header) as bam:
    for name, ref, start, cigar, reverse, secondary, supplementary in records:
        read = pysam.AlignedSegment()
        read.query_name = name
        read.query_sequence = "A" * sum(length for op, length in cigar if op in (0, 1, 4, 7, 8))
        read.flag = (16 if reverse else 0) | (256 if secondary else 0) | (2048 if supplementary else 0)
        read.reference_id = bam.get_tid(ref)
        read.reference_start = start
        read.mapping_quality = 60
        read.cigartuples = cigar
        bam.write(read)
print("fixture 已生成：", HERE)
