"""``plot insertion``：一致性序列规则、表契约、黄金对账、画法与错误路径。

全部用合成 fixture（``fixtures/insertion/reads.bam``，chr1 200 kb，坐标全是编造的）；
黄金 ``fixtures/insertion/expected.consensus.legacy.tsv`` 由并入前的独立脚本生成
（见 ``fixtures/insertion/make_fixture.py`` 的 docstring）。
"""

from __future__ import annotations

import shutil
import struct
from dataclasses import replace
from pathlib import Path

import pysam
import pytest

from alnkit.cli import main
from alnkit.core import tables
from alnkit.core.errors import RunError, UsageError
from alnkit.plot import insertion as insertion_plot
from alnkit.shortread import insertion

FIXTURES = Path(__file__).parent / "fixtures" / "insertion"
BAM = FIXTURES / "reads.bam"
LEGACY_CONSENSUS = FIXTURES / "expected.consensus.legacy.tsv"

# fixture 的五个位点（0-based 断点边界）
BOTH_SIDES = 100_000
MULTI_LEFT = 120_000
MULTI_RIGHT = 120_230
LEFT_ONLY = 150_000
RIGHT_ONLY = 180_000


def png_size(path: Path) -> tuple[int, int]:
    """PNG 的 (宽, 高)：直接读 IHDR，不引入 Pillow 依赖。"""
    header = Path(path).read_bytes()[:24]
    return struct.unpack(">II", header[16:24])


@pytest.fixture(scope="module")
def clip_reads() -> list[insertion.ClipRead]:
    """整个 fixture BAM 的 clip 证据（窗口覆盖全部五个位点）。"""
    return insertion.load_clip_reads(BAM, "chr1", 0, 200_000)


# ------------------------------------------------------------------ soft clip
@pytest.mark.parametrize(
    ("cigar", "sequence_length", "expected"),
    [
        (((4, 60), (0, 100)), 160, ("left", 60)),  # 60S100M
        (((0, 100), (4, 55)), 155, ("right", 55)),  # 100M55S
        (((0, 100), (5, 7)), 100, ("none", 0)),  # 100M7H：硬剪贴没有碱基
        (((5, 7), (4, 60), (0, 100)), 160, ("left", 60)),  # 7H60S100M：H 不能遮住 S
        (((4, 100),), 100, ("none", 0)),  # 整条都是 clip：异常记录
        (((4, 10), (0, 80), (4, 20)), 110, ("right", 20)),  # 两侧都有：取较长一侧
        (((4, 20), (0, 80), (4, 10)), 110, ("left", 20)),
        (((4, 15), (0, 80), (4, 15)), 110, ("left", 15)),  # 并列取 left
    ],
)
def test_identify_soft_clip(cigar, sequence_length, expected) -> None:
    alignment = pysam.AlignedSegment()
    alignment.query_name = "read"
    alignment.query_sequence = "ACGT" * (sequence_length // 4) + "A" * (sequence_length % 4)
    alignment.query_qualities = pysam.qualitystring_to_array("I" * sequence_length)
    alignment.reference_id = 0
    alignment.reference_start = 0
    alignment.mapping_quality = 60
    alignment.cigartuples = cigar

    side, length, sequence = insertion.identify_soft_clip(alignment)
    assert (side, length) == expected
    assert len(sequence) == expected[1]


def test_identify_soft_clip_uppercases_sequence() -> None:
    """小写碱基要统一成大写（BAM 里通常是大写，但不能假设）。"""
    alignment = pysam.AlignedSegment()
    alignment.query_name = "read"
    alignment.query_sequence = "acgtACGT" + "ACGT" * 2
    alignment.query_qualities = pysam.qualitystring_to_array("I" * 16)
    alignment.reference_id = 0
    alignment.reference_start = 0
    alignment.mapping_quality = 60
    alignment.cigartuples = [(4, 8), (0, 8)]

    side, length, sequence = insertion.identify_soft_clip(alignment)
    assert (side, length, sequence) == ("left", 8, "ACGTACGT")


# ------------------------------------------------------- 一致性序列（投票规则）
def _read(clip_sequence: str, *, clip_side: str = "left") -> insertion.ClipRead:
    """造一条只用于投票的最小 ClipRead。"""
    return insertion.ClipRead(
        read_name="r",
        ref_start=0,
        ref_end=100,
        read_length=100 + len(clip_sequence),
        strand="+",
        mapq=60,
        clip_side=clip_side,
        clip_length=len(clip_sequence),
        clip_sequence=clip_sequence,
        breakpoint=0,
    )


def test_consensus_requires_two_reads() -> None:
    sequence, _depths, count, _fraction, reason = insertion.build_clip_consensus(
        [_read("A" * 10)], "left"
    )
    assert (count, sequence) == (1, "")
    assert reason.startswith("coverage<2")


def test_consensus_majority_wins() -> None:
    reads = [_read("AAAAAAAAAC"), _read("AAAAAAAAAA"), _read("AAAAAAAAAA")]
    sequence, _depths, _count, fraction, reason = insertion.build_clip_consensus(reads, "left")
    # 左 clip 反向对齐：第 0 列是 1 个 C + 2 个 A；“取第一条 read”会错选成 C
    assert sequence[0] == "A"
    assert 0.0 < fraction < 1.0
    assert reason == "all columns"


def test_consensus_stops_on_two_read_mismatch() -> None:
    reads = [_read("AAAAAAAAAA"), _read("CCCCCCCCCC")]
    sequence, _depths, _count, _fraction, reason = insertion.build_clip_consensus(reads, "left")
    assert sequence == ""
    assert reason == "2-read overlap mismatch"


def test_consensus_aligns_at_junction_end() -> None:
    reads = [_read("ACGTACGTAC"), _read("ACGTACGTAC")]
    sequence, _depths, _count, fraction, _reason = insertion.build_clip_consensus(reads, "left")
    assert sequence == "CATGCATGCA"  # 左 clip 反向，接头端在第 0 列
    assert fraction == pytest.approx(1.0)


def test_consensus_short_read_truncates() -> None:
    reads = [_read("A" * 10), _read("AAAAAA")]
    sequence, _depths, _count, _fraction, reason = insertion.build_clip_consensus(reads, "left")
    assert sequence == "A" * 6
    assert reason == "coverage<2 (only 1 reads)"


def test_consensus_without_reads() -> None:
    assert insertion.build_clip_consensus([], "left") == ("", [], 0, 0.0, "no clip")


# ------------------------------------------------------- 读 BAM（fixture 口径）
def test_fixture_counts_only_primary_soft_clips(clip_reads) -> None:
    """硬剪贴 / secondary / supplementary 都不该进入证据。"""
    left_only = [read for read in clip_reads if read.breakpoint == LEFT_ONLY]
    assert len(left_only) == 15
    assert {read.clip_length for read in left_only} == {60}
    assert {read.clip_side for read in left_only} == {"left"}
    both = [read for read in clip_reads if read.breakpoint == BOTH_SIDES]
    assert len(both) == 40


def test_fixture_keeps_reverse_strand_reads(clip_reads) -> None:
    left_only = [read for read in clip_reads if read.breakpoint == LEFT_ONLY]
    assert {read.strand for read in left_only} == {"+", "-"}
    _sequence, _depths, count, fraction, reason = insertion.build_clip_consensus(left_only, "left")
    assert (count, fraction, reason) == (15, pytest.approx(1.0), "all columns")


def test_validate_breakpoints(clip_reads) -> None:
    assert insertion.validate_breakpoints([BOTH_SIDES, BOTH_SIDES], clip_reads) == (BOTH_SIDES,)
    with pytest.raises(UsageError, match="最多给 3 个"):
        insertion.validate_breakpoints([1, 2, 3, 4], clip_reads)
    with pytest.raises(UsageError, match="没有任何 clip read"):
        insertion.validate_breakpoints([123], clip_reads)
    with pytest.raises(UsageError, match="至少要给 1 个"):
        insertion.validate_breakpoints([], clip_reads)


def test_breakpoint_from_1based_is_identity() -> None:
    assert insertion.breakpoint_from_1based(100_000) == 100_000
    with pytest.raises(UsageError):
        insertion.breakpoint_from_1based(0)


# ------------------------------------------------------------------- 表契约
def test_clip_table_round_trip(tmp_path: Path, clip_reads) -> None:
    path = tmp_path / "clips.tsv"
    written = insertion.write_clip_table(path, clip_reads, "chr1", "sampleA", ("producer=test",))
    assert written == len(clip_reads)

    lines = path.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "# producer=test"
    assert lines[1].split("\t") == list(tables.columns_of("insertion_clip_reads"))

    table_reads, chrom = insertion.read_clip_table(path)
    assert chrom == "chr1"
    assert table_reads == sorted(
        clip_reads, key=lambda read: (read.breakpoint, read.clip_side, read.read_name)
    )


def test_consensus_table_uses_boundary_convention(tmp_path: Path, clip_reads) -> None:
    groups = insertion.build_groups(clip_reads, [LEFT_ONLY])
    path = tmp_path / "consensus.tsv"
    insertion.write_consensus_table(path, groups, "chr1", "sampleA")

    _header, rows = tables.read_tsv(path)
    assert rows[0]["breakpoint_0based"] == str(LEFT_ONLY)
    # 单点坐标：1-based 列取“左侧碱基”，数值与 0-based 边界相同（不是边界 + 1）
    assert rows[0]["breakpoint_1based"] == str(LEFT_ONLY)
    assert rows[0]["consensus_bp"] == "60"
    assert rows[0]["consensus_fraction"] == "1.0000"
    assert rows[0]["stop_reason"] == "all columns"


def test_read_clip_table_rejects_foreign_table(tmp_path: Path) -> None:
    path = tmp_path / "wrong.tsv"
    path.write_text("chrom\tpos_1based\tdepth\nchr1\t1\t2\n", encoding="utf-8")
    with pytest.raises(RunError, match="列与注册表不一致"):
        insertion.read_clip_table(path)


def test_consensus_golden_matches_standalone_script(tmp_path: Path, clip_reads) -> None:
    """与并入前的独立脚本逐列对账（口径差 1 已在 fixture 注释里说明）。"""
    breakpoints = [BOTH_SIDES, MULTI_LEFT, MULTI_RIGHT, LEFT_ONLY, RIGHT_ONLY]
    groups = insertion.build_groups(clip_reads, breakpoints)
    path = tmp_path / "consensus.tsv"
    insertion.write_consensus_table(path, groups, "chr1", "sampleA")
    _header, produced = tables.read_tsv(path)
    produced_by_key = {(int(row["breakpoint_0based"]), row["clip_side"]): row for row in produced}

    legacy_header, legacy_rows = tables.read_legacy_tsv(LEGACY_CONSENSUS)
    assert legacy_header == [
        "breakpoint_1based",
        "clip_side",
        "n_reads",
        "consensus_bp",
        "consensus_fraction",
        "stop_reason",
        "consensus_sequence",
    ]
    assert len(legacy_rows) == len(produced)
    for breakpoint_1based, side, n_reads, consensus_bp, fraction, reason, sequence in legacy_rows:
        key = (int(breakpoint_1based) - 1, side)
        row = produced_by_key[key]
        assert row["n_reads"] == n_reads, key
        assert row["consensus_bp"] == consensus_bp, key
        assert row["consensus_fraction"] == fraction, key
        assert row["stop_reason"] == reason, key
        assert row["consensus_sequence"] == sequence, key


def test_fixture_coordinates_stay_below_one_megabase() -> None:
    """fixture 的坐标必须都是编造的小坐标（< 1 Mb），避免真实基因组坐标入仓。"""
    _header, legacy_rows = tables.read_legacy_tsv(LEGACY_CONSENSUS)
    assert legacy_rows
    for row in legacy_rows:
        assert int(row[0]) < 1_000_000


# ------------------------------------------------------------------ 画法选择
def test_single_side_figure(tmp_path: Path, clip_reads) -> None:
    groups = insertion.build_groups(clip_reads, [RIGHT_ONLY])
    output = tmp_path / "single.png"
    layout = insertion_plot.draw_single_breakpoint_figure(
        groups, "chr1", 179_000, 181_000, RIGHT_ONLY, output
    )
    assert output.exists() and output.stat().st_size > 20_000
    assert layout["junctions"] == [RIGHT_ONLY]
    assert layout["groups"] == [("right", 15)]


def test_both_sides_figure_is_wider_than_single(tmp_path: Path, clip_reads) -> None:
    """双侧图多画了合并未知区与两条连线，文件应比单侧大。"""
    single = tmp_path / "single.png"
    both = tmp_path / "both.png"
    insertion_plot.draw_single_breakpoint_figure(
        insertion.build_groups(clip_reads, [LEFT_ONLY]), "chr1", 149_000, 151_000, LEFT_ONLY, single
    )
    insertion_plot.draw_single_breakpoint_figure(
        insertion.build_groups(clip_reads, [BOTH_SIDES]),
        "chr1",
        99_000,
        101_000,
        BOTH_SIDES,
        both,
    )
    assert both.stat().st_size > single.stat().st_size


def test_multi_breakpoint_figure(tmp_path: Path, clip_reads) -> None:
    groups = insertion.build_groups(clip_reads, [MULTI_LEFT, MULTI_RIGHT])
    output = tmp_path / "multi.png"
    layout = insertion_plot.draw_multiple_breakpoints_figure(
        groups, "chr1", 119_000, 121_000, output
    )
    assert output.exists() and output.stat().st_size > 20_000
    assert layout["junctions"] == [MULTI_LEFT, MULTI_RIGHT]


def test_render_without_notes_is_smaller(tmp_path: Path, clip_reads) -> None:
    groups = insertion.build_groups(clip_reads, [BOTH_SIDES])
    full = tmp_path / "full.png"
    plain = tmp_path / "plain.png"
    style = {"show_consensus": False, "show_formula": False}
    insertion_plot.render(groups, "chr1", 99_000, 101_000, [BOTH_SIDES], full)
    insertion_plot.render(groups, "chr1", 99_000, 101_000, [BOTH_SIDES], plain, **style)
    assert plain.stat().st_size < full.stat().st_size


def test_dpi_changes_pixel_size(tmp_path: Path, clip_reads) -> None:
    groups = insertion.build_groups(clip_reads, [LEFT_ONLY])
    low = tmp_path / "low.png"
    high = tmp_path / "high.png"
    insertion_plot.render(groups, "chr1", 149_000, 151_000, [LEFT_ONLY], low, dpi=80)
    insertion_plot.render(groups, "chr1", 149_000, 151_000, [LEFT_ONLY], high, dpi=160)
    assert png_size(high)[0] > png_size(low)[0]


def test_narrow_window_expands_xlim(tmp_path: Path, clip_reads) -> None:
    """窗口很窄时也要把 x 轴扩到容纳整个块，否则未知条会被裁掉。"""
    groups = insertion.build_groups(clip_reads, [LEFT_ONLY])
    layout = insertion_plot.render(
        groups, "chr1", 149_990, 150_010, [LEFT_ONLY], tmp_path / "narrow.png"
    )
    x_min, x_max = layout["xlim"]
    assert x_min < 149_990 and x_max > 150_010
    assert x_min <= LEFT_ONLY - 60  # 块起点在范围内


def test_plan_multi_breakpoint_blocks_keeps_gap() -> None:
    """相距 100 bp 的两个同侧断点：前一块的未知区必须被裁短、留出 gap。"""
    read = _read("A" * 60)
    groups = [
        insertion.ClipGroup(
            breakpoint=100_000,
            side="left",
            reads=(replace(read, breakpoint=100_000),),
            consensus="A" * 60,
            consensus_fraction=1.0,
            stop_reason="all columns",
        ),
        insertion.ClipGroup(
            breakpoint=100_100,
            side="left",
            reads=(replace(read, breakpoint=100_100),),
            consensus="A" * 60,
            consensus_fraction=1.0,
            stop_reason="all columns",
        ),
    ]
    blocks, (x_min, x_max) = insertion_plot.plan_multi_breakpoint_blocks(groups, 99_000, 101_000)
    ordered = sorted(blocks, key=lambda block: block[0])
    for left_block, right_block in zip(ordered, ordered[1:], strict=False):
        assert left_block[1] <= right_block[0], f"块重叠：{left_block} vs {right_block}"
    assert x_min <= ordered[0][0] and x_max >= ordered[-1][1]


# --------------------------------------------------------------- CLI 端到端
def test_cli_writes_png_and_both_tables(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "plot",
            "insertion",
            "--bam",
            str(BAM),
            "--region",
            "chr1:119001-121000",
            "--breakpoint",
            str(MULTI_LEFT),
            "--breakpoint",
            str(MULTI_RIGHT),
            "--out-dir",
            str(tmp_path),
            "--name",
            "demo",
            "--dpi",
            "100",
        ]
    )
    assert code == 0
    output = capsys.readouterr().out
    assert "demo.png" in output

    assert (tmp_path / "demo.png").stat().st_size > 20_000
    assert (tmp_path / "demo.clips.tsv").exists()
    _header, rows = tables.read_tsv(tmp_path / "demo.consensus.tsv")
    assert sorted(row["clip_side"] for row in rows) == ["left", "right"]


def test_cli_no_formula_makes_smaller_figure(tmp_path: Path) -> None:
    base = [
        "plot",
        "insertion",
        "--bam",
        str(BAM),
        "--region",
        "chr1:99101-100900",
        "--breakpoint",
        str(BOTH_SIDES),
        "--out-dir",
        str(tmp_path),
        "--dpi",
        "100",
    ]
    assert main(base + ["--name", "with_formula"]) == 0
    assert main(base + ["--name", "no_formula", "--no-formula", "--no-consensus-text"]) == 0
    full = (tmp_path / "with_formula.png").stat().st_size
    plain = (tmp_path / "no_formula.png").stat().st_size
    assert plain < full  # 烟测：说明文字确实进了图


def test_cli_rejects_name_with_directory(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "plot",
            "insertion",
            "--bam",
            str(BAM),
            "--region",
            "chr1:149001-151000",
            "--breakpoint",
            str(LEFT_ONLY),
            "--out-dir",
            str(tmp_path),
            "--name",
            "sub/demo",
        ]
    )
    assert code == 2
    assert "--name" in capsys.readouterr().err


def test_cli_rejects_unknown_breakpoint(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "plot",
            "insertion",
            "--bam",
            str(BAM),
            "--region",
            "chr1:149001-151000",
            "--breakpoint",
            "123",
            "--out-dir",
            str(tmp_path),
        ]
    )
    assert code == 2
    assert "没有任何 clip read" in capsys.readouterr().err


def test_cli_rejects_bad_region(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "plot",
            "insertion",
            "--bam",
            str(BAM),
            "--region",
            "chr1:300-100",
            "--breakpoint",
            "200",
            "--out-dir",
            str(tmp_path),
        ]
    )
    assert code == 2
    assert "--region" in capsys.readouterr().err


def test_cli_rejects_missing_index(tmp_path: Path, capsys) -> None:
    """删掉 .bai 后要给可行动的提示（退出码 2），而不是 pysam 的 traceback。"""
    copy = tmp_path / "noindex.bam"
    shutil.copy(BAM, copy)
    code = main(
        [
            "plot",
            "insertion",
            "--bam",
            str(copy),
            "--region",
            "chr1:149001-151000",
            "--breakpoint",
            str(LEFT_ONLY),
            "--out-dir",
            str(tmp_path),
        ]
    )
    assert code == 2
    assert "索引" in capsys.readouterr().err


def test_cli_rejects_unknown_reference(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "plot",
            "insertion",
            "--bam",
            str(BAM),
            "--region",
            "chrZ:1-1000",
            "--breakpoint",
            "500",
            "--out-dir",
            str(tmp_path),
        ]
    )
    assert code == 2
    assert "chr1" in capsys.readouterr().err
