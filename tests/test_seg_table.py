"""``alnkit seg table``（D2）：CIGAR 简化定点断言 + 与旧脚本的逐字节对账。

``expected.legacy.tsv`` 由旧脚本 ``upstream/scripts/alignment_table.py`` 在
``fixtures/segtable/reads.bam`` + ``ref.fa.fai`` 上生成，用于钉住历史 14 列表的
列名、列序、排序与取值。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from alnkit.cli import main
from alnkit.core import tables
from alnkit.core.bam_segments import (
    cigar_ops_from_pysam,
    format_cigar_ops,
    load_alignment_segments,
    load_fai_lengths,
    simplify_cigar_ops,
)
from alnkit.core.errors import UsageError

FIXTURES = Path(__file__).parent / "fixtures" / "segtable"
BAM = FIXTURES / "reads.bam"
FAI = FIXTURES / "ref.fa.fai"
EXPECTED_LEGACY = FIXTURES / "expected.legacy.tsv"


# ------------------------------------------------------------ CIGAR 简化
def test_simplify_absorbs_small_indels_between_matches() -> None:
    ops = [(39, "M"), (1, "D"), (66, "M"), (1, "D"), (92, "M"), (1, "I"), (445, "M")]
    assert simplify_cigar_ops(ops) == [(645, "M")]
    assert format_cigar_ops(simplify_cigar_ops(ops)) == "645M"


def test_simplify_keeps_big_indels_and_clips() -> None:
    ops = [(10, "S"), (100, "M"), (10, "D"), (100, "M"), (20, "S")]
    assert simplify_cigar_ops(ops) == [(10, "S"), (100, "M"), (10, "D"), (100, "M"), (20, "S")]


def test_simplify_does_not_touch_skips_or_edge_indels() -> None:
    # N（skip）不算可吸收差异；开头/结尾的 I 也不吸收
    ops = [(5, "I"), (100, "M"), (30, "N"), (100, "M"), (5, "D")]
    assert simplify_cigar_ops(ops) == [(5, "I"), (100, "M"), (30, "N"), (100, "M"), (5, "D")]


def test_simplify_merges_adjacent_same_ops() -> None:
    assert simplify_cigar_ops([(10, "M"), (5, "M")]) == [(15, "M")]


# --------------------------------------------------------------- 段抽取
def test_load_segments_reads_axis_and_spans() -> None:
    segments = load_alignment_segments(BAM, chrom_lengths=load_fai_lengths(FAI))
    assert len(segments) == 8
    by_key = {(s.read_name, s.chrom, s.ref_start_0based): s for s in segments}

    forward = by_key[("readA", "chr1", 100)]
    assert (forward.read_start_0based, forward.read_end_0based) == (0, 100)
    assert forward.ref_span == 100 == forward.ref_aln_len
    assert forward.cigar == "100M50S"

    # readB 是反向：60S140M，正向 read 轴上应仍是 0..140
    reverse = by_key[("readB", "chr1", 400)]
    assert reverse.strand == "-"
    assert (reverse.read_start_0based, reverse.read_end_0based) == (0, 140)

    # readC 的 120M30N50M20S：ref_aln_len 含 N 跨度，read 轴只算 170
    skipped = by_key[("readC", "chr1", 800)]
    assert skipped.ref_span == 200
    assert skipped.read_span == 170
    assert skipped.read_aln_len == 170


def test_load_segments_filters() -> None:
    only = load_alignment_segments(BAM, selected_read_names={"readA"})
    assert {s.read_name for s in only} == {"readA"}
    no_secondary = load_alignment_segments(BAM, skip_secondary=True)
    assert not any(s.is_secondary for s in no_secondary)
    assert len(no_secondary) == 7


def test_missing_bam_or_fai_raise_usage_error(tmp_path: Path) -> None:
    with pytest.raises(UsageError, match="找不到 BAM"):
        load_alignment_segments(tmp_path / "nope.bam")
    with pytest.raises(UsageError, match="找不到 .fai"):
        load_fai_lengths(tmp_path / "nope.fai")


# ------------------------------------------------------------ 黄金对账
def test_legacy_columns_match_old_script_byte_for_byte(tmp_path: Path) -> None:
    target = tmp_path / "legacy.tsv"
    assert (
        main(
            [
                "seg",
                "table",
                "--bam",
                str(BAM),
                "--fai",
                str(FAI),
                "--out",
                str(target),
                "--columns",
                "legacy",
            ]
        )
        == 0
    )
    assert target.read_bytes() == EXPECTED_LEGACY.read_bytes(), "历史 14 列表与旧脚本不再一致"


def test_canonical_columns_use_registered_order(tmp_path: Path) -> None:
    target = tmp_path / "canonical.tsv"
    assert main(["seg", "table", "--bam", str(BAM), "--out", str(target)]) == 0
    header, rows = tables.read_tsv(target)
    assert header == list(tables.columns_of("segments"))
    assert len(rows) == 8
    first = rows[0]
    assert first["ref_start_0based"] == "100"
    assert first["is_secondary"] == "0"


LEGACY_TO_CANONICAL = {
    "chrom": "chrom",
    "chrom_length": "chrom_length",
    "ref_pos_start": "ref_start_0based",
    "ref_pos_end": "ref_end_0based",
    "read_name": "read_name",
    "read_pos_start": "read_start_0based",
    "read_pos_end": "read_end_0based",
    "strand": "strand",
    "mapping_quality": "mapq",
    "is_secondary": "is_secondary",
    "is_supplementary": "is_supplementary",
    "ref_aln_len": "ref_aln_len",
    "read_aln_len": "read_aln_len",
    "simplified_cigar": "cigar",
}


def test_canonical_values_match_legacy_golden(tmp_path: Path) -> None:
    """规范列的值必须与历史 14 列黄金逐字段一致（只差列名）。

    这条替代了「只看首行两个字段」的弱断言：14 列 × 8 行全部比对，
    列错位（例如 read_start/read_end 互换、cigar 落错列）会立刻暴露。
    """
    target = tmp_path / "canonical.tsv"
    assert main(["seg", "table", "--bam", str(BAM), "--fai", str(FAI), "--out", str(target)]) == 0

    _, canonical = tables.read_tsv(target)
    _, legacy = tables.read_tsv(EXPECTED_LEGACY)
    assert len(canonical) == len(legacy) == 8

    for index, (new_row, old_row) in enumerate(zip(canonical, legacy, strict=True)):
        for legacy_name, canonical_name in LEGACY_TO_CANONICAL.items():
            assert new_row[canonical_name] == old_row[legacy_name], (
                f"第 {index + 1} 行 {canonical_name} != legacy {legacy_name}"
            )


def test_read_file_filter(tmp_path: Path) -> None:
    names = tmp_path / "reads.txt"
    names.write_text("# 注释\nreadD\n", encoding="utf-8")
    target = tmp_path / "only_d.tsv"
    assert (
        main(["seg", "table", "--bam", str(BAM), "--out", str(target), "--read-file", str(names)])
        == 0
    )
    _, rows = tables.read_tsv(target)
    assert [row["read_name"] for row in rows] == ["readD"]


def test_missing_bam_is_usage_error(tmp_path: Path) -> None:
    assert (
        main(
            ["seg", "table", "--bam", str(tmp_path / "nope.bam"), "--out", str(tmp_path / "x.tsv")]
        )
        == 2
    )


def test_cigar_ops_from_pysam_fixture() -> None:
    import pysam

    with pysam.AlignmentFile(str(BAM), "rb") as handle:
        read = next(iter(handle))
    assert cigar_ops_from_pysam(read) == [(100, "M"), (50, "S")]
