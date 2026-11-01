"""``alnkit plot tracks``（C2）：纯逻辑定点断言 + 与旧脚本产物的黄金回归。

黄金文件 ``fixtures/tracks/expected.segments.tsv`` 由旧的单文件实现
``read_segment_tracks.py`` 生成（同一份 ``example_segments.tsv`` 输入），
用于证明移植后的输出逐字节一致。
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from alnkit.cli import main
from alnkit.core.errors import UsageError
from alnkit.plot import tracks

FIXTURES = Path(__file__).parent / "fixtures" / "tracks"
INPUT_TSV = FIXTURES / "example_segments.tsv"
EXPECTED_TSV = FIXTURES / "expected.segments.tsv"


# ---------------------------------------------------------------- 纯逻辑
def test_parse_region_accepts_commas() -> None:
    assert tracks.parse_region("chr1:3,101,000") == ("chr1", 3101000)


@pytest.mark.parametrize("text", ["chr1", ":100", "chr1:abc", "chr1:-5"])
def test_parse_region_rejects_bad_input(text: str) -> None:
    with pytest.raises(UsageError):
        tracks.parse_region(text)


@pytest.mark.parametrize(
    ("cigar", "read_length", "is_reverse", "expected"),
    [
        (((5, 5), (0, 100)), 105, True, (0, 100)),  # 5H100M 反向：硬剪切不计入正向轴
        (((5, 5), (0, 100)), 105, False, (5, 105)),
        (((0, 100),), 100, False, (0, 100)),
        (((4, 10), (0, 100)), 110, True, (0, 100)),  # 10S100M 反向：软剪切落在 3' 端
        (((4, 10), (0, 100), (4, 10)), 120, False, (10, 110)),
    ],
)
def test_cigar_to_read_span(cigar, read_length, is_reverse, expected) -> None:
    assert tracks.cigar_to_read_span(cigar, read_length, is_reverse) == expected


def test_cigar_to_read_span_returns_none_without_match() -> None:
    assert tracks.cigar_to_read_span(None, 100, False) is None
    assert tracks.cigar_to_read_span(((4, 100),), 100, False) is None


def test_cigar_ref_to_read_walks_cigar() -> None:
    cigar = ((0, 100), (2, 100), (0, 100))  # 100M 100D 100M
    assert tracks.cigar_ref_to_read(cigar, 1000, 200, False, 1200) == 100
    assert tracks.cigar_ref_to_read(cigar, 1000, 200, True, 1200) == 99
    assert tracks.cigar_ref_to_read(cigar, 1000, 200, False, 1100) is None  # 落在缺失里


def test_infer_segment_type_and_labels() -> None:
    assert tracks.infer_segment_type("exog_a", {"exog_a", "exog_b"}) == "exog"
    assert tracks.infer_segment_type("chr1", {"exog_a"}) == "genome"
    assert tracks.short_read_label("sample/1002/ccs") == "1002"
    assert tracks.short_read_label("plain_read") == "plain_read"


def test_requires_exactly_one_source() -> None:
    with pytest.raises(UsageError, match="必须且只能指定"):
        tracks.run(tracks.TrackSettings(outdir="."))
    with pytest.raises(UsageError, match="必须且只能指定"):
        tracks.run(tracks.TrackSettings(outdir=".", bam="a.bam", paf="b.paf"))


# ------------------------------------------------------------ 黄金回归
def _run_fixture(tmp_path: Path) -> dict[str, object]:
    return tracks.run(
        tracks.TrackSettings(
            segments_tsv=str(INPUT_TSV),
            exog=("exog_a", "exog_b"),
            highlight_ref=("chr1:3101000",),
            outdir=str(tmp_path),
            name="golden",
        )
    )


def test_golden_segments_tsv_is_byte_identical(tmp_path: Path) -> None:
    result = _run_fixture(tmp_path)
    assert (result["reads"], result["segments"]) == (5, 18)
    produced = Path(str(result["segments_tsv"])).read_bytes()
    assert produced == EXPECTED_TSV.read_bytes(), "片段表与旧脚本产物不再一致"


def test_outputs_are_written(tmp_path: Path) -> None:
    result = _run_fixture(tmp_path)
    for key in ("png", "svg", "pdf"):
        path = Path(str(result[key]))
        assert path.exists(), key
        assert path.stat().st_size > 5_000, key
    svg = Path(str(result["svg"])).read_text(encoding="utf-8")
    assert "read1" in svg  # read 短标签
    assert "exog" in svg


def test_cli_plot_tracks_end_to_end(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "plot",
            "tracks",
            "--segments-tsv",
            str(INPUT_TSV),
            "--exog",
            "exog_a,exog_b",
            "--highlight-ref",
            "chr1:3101000",
            "--outdir",
            str(tmp_path),
            "--name",
            "cli",
        ]
    )
    assert code == 0
    out = capsys.readouterr().out
    assert "reads: 5" in out
    assert (tmp_path / "cli.segments.tsv").exists()


def test_cli_plot_requires_subcommand() -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(["plot"])
    assert excinfo.value.code == 2


def test_bam_input_smoke(tmp_path: Path) -> None:
    """BAM 路径冒烟：造一条 2 段比对 + 1 条外源比对的 read。"""
    reference = tmp_path / "ref.fa"
    reference.write_text(">chr1\n" + "ACGT" * 250 + "\n>exog_a\n" + "TTTT" * 50 + "\n")
    pysam.faidx(str(reference))

    header = {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 1000}, {"SN": "exog_a", "LN": 200}]}
    bam_path = tmp_path / "reads.bam"
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as bam:
        for name, ref, start, cigar in (
            ("r1", "chr1", 100, [(0, 100), (4, 200)]),
            ("r1", "exog_a", 0, [(4, 200), (0, 100)]),
        ):
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * 300
            read.flag = 16  # 反向
            read.reference_id = bam.get_tid(ref)
            read.reference_start = start
            read.mapping_quality = 60
            read.cigartuples = cigar
            bam.write(read)

    settings = tracks.TrackSettings(
        bam=str(bam_path),
        read=("r1",),
        exog=("exog_a",),
        outdir=str(tmp_path),
        name="bam",
    )
    result = tracks.run(settings)
    assert result["reads"] == 1
    lines = Path(str(result["segments_tsv"])).read_text(encoding="utf-8").splitlines()
    assert len(lines) == 3  # 表头 + 2 段
