"""比对段表归一（core/segments.py）：三套历史表头读到同一个 Segment。"""

from __future__ import annotations

import pytest

from alnkit.core import segments
from alnkit.core.errors import UsageError

# 三套真实表头（取自坐标审计：alignment_table / 03_classify_reads / read_segment_tracks）
HEADER_ALIGNMENT_TABLE = [
    "chrom",
    "chrom_length",
    "ref_pos_start",
    "ref_pos_end",
    "read_name",
    "read_pos_start",
    "read_pos_end",
    "strand",
    "mapping_quality",
    "is_secondary",
    "is_supplementary",
    "ref_aln_len",
    "read_aln_len",
    "simplified_cigar",
]
HEADER_CLASSIFY = [
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
]
HEADER_TRACKS = [
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

ROW_ALIGNMENT_TABLE = {
    "chrom": "chr1",
    "chrom_length": "100000000",
    "ref_pos_start": "9000000",
    "ref_pos_end": "9031165",
    "read_name": "sampleA/read1/ccs",
    "read_pos_start": "9424",
    "read_pos_end": "40591",
    "strand": "-",
    "mapping_quality": "60",
    "is_secondary": "0",
    "is_supplementary": "0",
    "ref_aln_len": "31165",
    "read_aln_len": "31167",
    "simplified_cigar": "1S31204M9424S",
}
ROW_CLASSIFY = {
    "read_name": "sampleA/read1/ccs",
    "segment_role": "genome",
    "chrom": "chr1",
    "chrom_length": "100000000",
    "ref_start": "9000000",
    "ref_end": "9031165",
    "read_start": "9424",
    "read_end": "40591",
    "strand": "-",
    "mapq": "60",
    "ref_aln_len": "31165",
    "read_aln_len": "31167",
    "simplified_cigar": "1S31204M9424S",
    "kept": "1",
    "drop_reason": "-",
}
ROW_TRACKS = {
    "read_name": "sampleA/read1/ccs",
    "read_length": "40592",
    "segment_type": "genome",
    "target_name": "chr1",
    "target_start": "9000000",
    "target_end": "9031165",
    "read_start": "9424",
    "read_end": "40591",
    "strand": "-",
    "mapq": "60",
}


@pytest.mark.parametrize(
    ("header", "row"),
    [
        (HEADER_ALIGNMENT_TABLE, ROW_ALIGNMENT_TABLE),
        (HEADER_CLASSIFY, ROW_CLASSIFY),
        (HEADER_TRACKS, ROW_TRACKS),
    ],
    ids=["alignment_table", "classify", "tracks"],
)
def test_three_legacy_headers_yield_same_segment(header, row) -> None:
    segment = segments.segments_from_rows(header, [row])[0]
    assert (segment.ref_start_0based, segment.ref_end_0based) == (9000000, 9031165)
    assert (segment.read_start_0based, segment.read_end_0based) == (9424, 40591)
    assert segment.chrom == "chr1"
    assert segment.strand == "-"
    assert segment.mapq == 60


def test_spans_are_derived_from_coordinates() -> None:
    segment = segments.segments_from_rows(HEADER_CLASSIFY, [ROW_CLASSIFY])[0]
    assert segment.ref_span == 31165 == segment.ref_aln_len
    assert segment.read_span == 31167 == segment.read_aln_len


def test_optional_columns_become_none() -> None:
    header = [
        "read_name",
        "chrom",
        "ref_start",
        "ref_end",
        "read_start",
        "read_end",
        "strand",
        "mapq",
        "read_length",
        "cigar",
    ]
    row = {
        "read_name": "r",
        "chrom": "chr1",
        "ref_start": "0",
        "ref_end": "10",
        "read_start": "0",
        "read_end": "10",
        "strand": "+",
        "mapq": "60",
        "read_length": "-",
        "cigar": "10M",
    }
    segment = segments.segments_from_rows(header, [row])[0]
    assert segment.read_length is None
    assert segment.cigar == "10M"
    assert segment.kept is True  # 缺列时默认保留
    assert segment.segment_role is None


def test_canonical_names_take_precedence_over_aliases() -> None:
    header = [
        "read_name",
        "chrom",
        "ref_start_0based",
        "ref_start",
        "ref_end",
        "read_start",
        "read_end",
        "strand",
        "mapq",
    ]
    row = {
        "read_name": "r",
        "chrom": "chr1",
        "ref_start_0based": "100",
        "ref_start": "999",
        "ref_end": "200",
        "read_start": "0",
        "read_end": "100",
        "strand": "+",
        "mapq": "60",
    }
    assert segments.segments_from_rows(header, [row])[0].ref_start_0based == 100


def test_stand_typo_is_accepted() -> None:
    header = [
        "read_name",
        "chrom",
        "ref_start",
        "ref_end",
        "read_start",
        "read_end",
        "stand",
        "mapq",
    ]
    row = {
        "read_name": "r",
        "chrom": "chr1",
        "ref_start": "0",
        "ref_end": "10",
        "read_start": "0",
        "read_end": "10",
        "stand": "-",
        "mapq": "60",
    }
    assert segments.segments_from_rows(header, [row])[0].strand == "-"


def test_missing_required_column_raises_usage_error() -> None:
    with pytest.raises(UsageError, match="缺少必需列"):
        segments.segments_from_rows(["read_name", "chrom"], [{"read_name": "r", "chrom": "chr1"}])


@pytest.mark.parametrize(
    "override",
    [
        {"ref_start": "10", "ref_end": "5"},
        {"read_start": "9", "read_end": "1"},
        {"strand": "?"},
        {"mapq": ""},
        {"ref_start": "abc"},
    ],
)
def test_bad_values_raise_usage_error(override) -> None:
    row = {**ROW_CLASSIFY, **override}
    with pytest.raises(UsageError):
        segments.segments_from_rows(HEADER_CLASSIFY, [row])


def test_read_segments_from_file(tmp_path) -> None:
    path = tmp_path / "alignment_segments.tsv"
    path.write_text(
        "# 注释行\n"
        + "\t".join(HEADER_CLASSIFY)
        + "\n"
        + "\t".join(ROW_CLASSIFY[column] for column in HEADER_CLASSIFY)
        + "\n",
        encoding="utf-8",
    )
    parsed = segments.read_segments(path)
    assert len(parsed) == 1
    assert parsed[0].ref_start_0based == 9000000


def test_segment_rows_matches_registered_columns() -> None:
    from alnkit.core import tables

    parsed = segments.segments_from_rows(HEADER_CLASSIFY, [ROW_CLASSIFY])
    row = segments.segment_rows(parsed)[0]
    assert set(row) == set(tables.columns_of("segments"))


def test_empty_strand_is_an_error_not_a_default() -> None:
    """strand 空值绝不能静默当 "+"：它决定断点取 ref_start 还是 ref_end。"""
    row = {**ROW_CLASSIFY, "strand": ""}
    with pytest.raises(UsageError, match="strand"):
        segments.segments_from_rows(HEADER_CLASSIFY, [row])


def test_bad_optional_integer_is_usage_error() -> None:
    row = {**ROW_CLASSIFY, "read_length": "abc"}
    header = [*HEADER_CLASSIFY, "read_length"]
    with pytest.raises(UsageError, match="read_length"):
        segments.segments_from_rows(header, [row])
