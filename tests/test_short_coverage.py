"""``short coverage``：深度口径、分箱、表契约、黄金对账、绘图与错误路径。

fixture（``fixtures/coverage/reads.bam``，chr1 2.2 kb，坐标全是编造的）里有四种 CIGAR
边界，黄金 ``fixtures/coverage/expected.bins.legacy.tsv`` 由并入前的独立脚本
（``coverage-viz --dump``）生成，见 ``fixtures/coverage/make_fixture.py``。
"""

from __future__ import annotations

import struct
from pathlib import Path

import pysam
import pytest

from alnkit.cli import main
from alnkit.core import tables
from alnkit.core.errors import RunError, UsageError
from alnkit.plot import coverage as coverage_plot
from alnkit.shortread import coverage as coverage_stats

FIXTURES = Path(__file__).parent / "fixtures" / "coverage"
BAM = FIXTURES / "reads.bam"
LEGACY_BINS = FIXTURES / "expected.bins.legacy.tsv"

#: fixture 里几条关键 read 的 0-based 起点与它们造成的“空洞”
DELETION_GAP = (1650, 1680)  # 30D：没有碱基 -> 深度 0
SOFTCLIP_GAP = (1730, 1740)  # 30S70M 之前的参考位置 -> 深度 0
INSERTION_SPAN = (1860, 1900)  # 20I 不占参考位置 -> 深度连续


def png_size(path: Path) -> tuple[int, int]:
    """PNG 的 (宽, 高)：直接读 IHDR，不引入 Pillow 依赖。"""
    header = Path(path).read_bytes()[:24]
    return struct.unpack(">II", header[16:24])


@pytest.fixture(scope="module")
def alignment_file():
    with pysam.AlignmentFile(str(BAM)) as handle:
        yield handle


def depths(alignment_file, region, filters=None) -> list[int]:
    """某区间的逐碱基深度（测试里只关心整数深度）。"""
    result = coverage_stats.compute_coverage(
        alignment_file, region, filters or coverage_stats.Filters(), requested_bins=1
    )
    return [int(value) for value in result.depth]


# ------------------------------------------------------------------ 参数解析
@pytest.mark.parametrize(
    ("text", "expected"),
    [("", 0), ("UNMAP,DUP", 0x4 | 0x400), ("0x504", 0x504), ("1284", 1284), ("qcfail", 0x200)],
)
def test_parse_flags(text: str, expected: int) -> None:
    assert coverage_stats.parse_flags(text) == expected


def test_parse_flags_rejects_unknown_name() -> None:
    with pytest.raises(UsageError, match="未知的 flag 名"):
        coverage_stats.parse_flags("NOT_A_FLAG")


# ---------------------------------------------------------------------- 分箱
@pytest.mark.parametrize(
    ("region_length", "bins", "expected"),
    [
        (100, 4, [0, 25, 50, 75, 100]),
        (103, 4, [0, 25, 50, 75, 103]),  # 最后一个 bin 收余数
        (5, 50, [0, 1, 2, 3, 4, 5]),  # bin 数不能超过碱基数
        (2000, 8, [0, 250, 500, 750, 1000, 1250, 1500, 1750, 2000]),
    ],
)
def test_make_bin_edges(region_length: int, bins: int, expected: list[int]) -> None:
    assert coverage_stats.make_bin_edges(region_length, bins).tolist() == expected


def test_default_bins() -> None:
    assert coverage_plot.default_bins(1416, None) == 60  # 下限
    assert coverage_plot.default_bins(100_000, None) == 2000  # 上限
    assert coverage_plot.default_bins(50, None) == 50  # 不超过碱基数
    assert coverage_plot.default_bins(1000, 7) == 7  # 显式指定就照办


# ------------------------------------------------------------------ 深度口径
def test_deletion_carries_no_bases(alignment_file) -> None:
    region = coverage_stats.Region("chr1", DELETION_GAP[0], DELETION_GAP[1])
    assert set(depths(alignment_file, region)) == {0}
    # 缺失两侧的 M 块各贡献 1
    assert depths(alignment_file, coverage_stats.Region("chr1", 1600, 1601)) == [1]
    assert depths(alignment_file, coverage_stats.Region("chr1", 1680, 1681)) == [1]


def test_soft_clip_carries_no_bases(alignment_file) -> None:
    region = coverage_stats.Region("chr1", SOFTCLIP_GAP[0], SOFTCLIP_GAP[1])
    assert set(depths(alignment_file, region)) == {0}
    assert depths(alignment_file, coverage_stats.Region("chr1", 1740, 1741)) == [1]


def test_insertion_keeps_reference_contiguous(alignment_file) -> None:
    region = coverage_stats.Region("chr1", INSERTION_SPAN[0], INSERTION_SPAN[1])
    assert set(depths(alignment_file, region)) == {1}


def test_primary_filter_skips_secondary_and_duplicate(alignment_file) -> None:
    region = coverage_stats.Region("chr1", 1890, 1930)
    default = depths(alignment_file, region)
    assert default[20] == 0  # 1910 起只有 secondary/duplicate，默认被过滤
    unfiltered = depths(alignment_file, region, coverage_stats.Filters(fail_flags=0, min_depth=1))
    assert unfiltered[20] == 2  # 两条都被算进来


def test_min_mapq_and_min_depth(alignment_file) -> None:
    region = coverage_stats.Region("chr1", 1920, 1930)  # 只有那条 mapQ=5 的 read
    assert depths(alignment_file, region) == [1] * 10
    filtered = coverage_stats.Filters(min_mapq=20)
    assert depths(alignment_file, region, filtered) == [0] * 10

    tiled = coverage_stats.Region("chr1", 485, 515)  # 两条平铺 read 互相重叠，深度 2
    assert set(depths(alignment_file, tiled)) == {2}
    capped = coverage_stats.Filters(max_depth=1)
    assert set(depths(alignment_file, tiled, capped)) == {1}


def test_region_summary_metrics(alignment_file) -> None:
    result = coverage_stats.compute_coverage(
        alignment_file,
        coverage_stats.Region("chr1", 0, 2000),
        coverage_stats.Filters(),
        requested_bins=8,
    )
    assert (
        result.read_count,
        result.selected_read_count,
        result.filtered_read_count,
    ) == (31, 29, 2)
    assert result.covered_bases == 1870
    assert result.percent_covered == pytest.approx(93.5)
    assert result.mean_coverage == pytest.approx(1.415)
    assert result.mean_mapq == pytest.approx(58.1034, abs=1e-4)
    assert result.mean_baseq == pytest.approx(40.0)  # fixture 的碱基质量统一是 Q40


def test_oversized_region_is_rejected(alignment_file) -> None:
    with pytest.raises(UsageError, match="缩小 --region"):
        coverage_stats.compute_coverage(
            alignment_file, coverage_stats.Region("chr1", 0, coverage_stats.MAX_REGION_LENGTH + 1)
        )


# --------------------------------------------------------------- 表契约/黄金
def test_bins_golden_matches_standalone_script(tmp_path: Path, alignment_file) -> None:
    result = coverage_stats.compute_coverage(
        alignment_file,
        coverage_stats.Region("chr1", 0, 2000),
        coverage_stats.Filters(),
        requested_bins=8,
    )
    bin_path = tmp_path / "bins.tsv"
    region_path = tmp_path / "region.tsv"
    coverage_stats.write_coverage_tables(
        result, "sampleA", region_path=region_path, bins_path=bin_path, comments=("producer=test",)
    )

    header, produced = tables.read_tsv(bin_path)
    assert list(header) == list(tables.columns_of("coverage_bins"))
    legacy_header, legacy_rows = tables.read_legacy_tsv(LEGACY_BINS)
    assert legacy_header == ["start", "end", "covered_bases", "breadth_pct", "mean_depth"]
    assert len(produced) == len(legacy_rows)

    for row, (start, end, covered, breadth, mean_depth) in zip(produced, legacy_rows, strict=True):
        assert row["bin_start_1based"] == start
        assert row["bin_end_1based"] == end
        assert row["covered_bases"] == covered
        width = int(row["bin_width_bp"])
        assert f"{100 * int(row['covered_bases']) / width:.4f}" == breadth
        assert f"{int(row['depth_sum']) / width:.4f}" == mean_depth

    region_header, region_rows = tables.read_tsv(region_path)
    assert list(region_header) == list(tables.columns_of("coverage_region"))
    assert region_rows[0]["n_reads_selected"] == "29"
    assert region_rows[0]["percent_covered"] == "93.5000"
    assert region_rows[0]["mean_depth"] == "1.4150"


def test_bases_table_covers_every_position(tmp_path: Path, alignment_file) -> None:
    result = coverage_stats.compute_coverage(
        alignment_file, coverage_stats.Region("chr1", 0, 10), requested_bins=2
    )
    path = tmp_path / "bases.tsv"
    coverage_stats.write_coverage_tables(
        result,
        "sampleA",
        region_path=tmp_path / "r.tsv",
        bins_path=tmp_path / "b.tsv",
        bases_path=path,
    )
    _header, rows = tables.read_tsv(path)
    assert [row["pos_1based"] for row in rows] == [str(position) for position in range(1, 11)]


# -------------------------------------------------------------------- 绘图
def test_bin_arrays_rejects_non_contiguous_bins() -> None:
    rows = [
        {
            "bin_start_1based": "1",
            "bin_end_1based": "10",
            "bin_width_bp": "10",
            "covered_bases": "10",
            "depth_sum": "20",
        },
        {
            "bin_start_1based": "20",
            "bin_end_1based": "30",
            "bin_width_bp": "10",
            "covered_bases": "0",
            "depth_sum": "0",
        },
    ]
    with pytest.raises(RunError, match="不相接"):
        coverage_plot.bin_arrays(rows)


def test_render_writes_png(tmp_path: Path, alignment_file) -> None:
    result = coverage_stats.compute_coverage(
        alignment_file, coverage_stats.Region("chr1", 0, 2000), requested_bins=8
    )
    region_path = tmp_path / "region.tsv"
    bins_path = tmp_path / "bins.tsv"
    coverage_stats.write_coverage_tables(
        result, "sampleA", region_path=region_path, bins_path=bins_path
    )
    _region_header, region_rows = tables.read_tsv(region_path)
    _bins_header, bins_rows = tables.read_tsv(bins_path)

    output = tmp_path / "coverage.png"
    coverage_plot.render(region_rows[0], bins_rows, output, width=8.0, dpi=100)
    width, height = png_size(output)
    assert width >= 700 and height >= 400
    assert output.stat().st_size > 20_000


def test_render_with_gene_track(tmp_path: Path, alignment_file) -> None:
    result = coverage_stats.compute_coverage(
        alignment_file, coverage_stats.Region("chr1", 0, 2000), requested_bins=8
    )
    region_path = tmp_path / "region.tsv"
    bins_path = tmp_path / "bins.tsv"
    coverage_stats.write_coverage_tables(
        result, "sampleA", region_path=region_path, bins_path=bins_path
    )
    _region_header, region_rows = tables.read_tsv(region_path)
    _bins_header, bins_rows = tables.read_tsv(bins_path)

    bed = tmp_path / "genes.bed"
    bed.write_text(
        "chr1\t0\t1200\tgeneA\t0\t+\t0\t1200\t59,141,239\n"
        "chr1\t1100\t1400\tgeneB\t0\t+\t1100\t1400\t.\n"
        "chr1\t1900\t2000\ttag\t0\t+\t1900\t2000\t255,0,0\n",
        encoding="utf-8",
    )
    genes = coverage_plot.parse_bed(bed)
    assert [gene.name for gene in genes] == ["geneA", "geneB", "tag"]
    assert genes[0].color == "#3b8def"
    assert genes[2].color == "#ff0000"

    plain = tmp_path / "plain.png"
    with_genes = tmp_path / "genes.png"
    coverage_plot.render(region_rows[0], bins_rows, plain, dpi=100)
    coverage_plot.render(region_rows[0], bins_rows, with_genes, genes, dpi=100)
    assert with_genes.stat().st_size > plain.stat().st_size
    assert png_size(with_genes)[1] > png_size(plain)[1]  # 多一块轨道 -> 更高


def test_render_handles_empty_region(tmp_path: Path) -> None:
    region_row = {
        "chrom": "chr1",
        "region_start_1based": "1",
        "region_end_1based": "100",
        "region_length_bp": "100",
        "n_bins": "1",
        "n_reads_selected": "0",
        "percent_covered": "0.0000",
        "mean_depth": "0.0000",
        "mean_mapq": "0.0000",
    }
    bins = [
        {
            "bin_index": "1",
            "bin_start_1based": "1",
            "bin_end_1based": "100",
            "bin_width_bp": "100",
            "covered_bases": "0",
            "depth_sum": "0",
        }
    ]
    output = tmp_path / "empty.png"
    coverage_plot.render(region_row, bins, output, dpi=80)
    assert output.stat().st_size > 5_000


# -------------------------------------------------------------------- BED
@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("255,128,0", "#ff8000"),
        ("12, 34, 56", "#0c2238"),
        ("#1e90ff", "#1e90ff"),
        ("300,0,0", "#ff0000"),  # 超过 255 会被截断
        ("0", None),
        ("tomato", None),
        (None, None),
    ],
)
def test_parse_bed_color(text, expected) -> None:
    assert coverage_plot.parse_bed_color(text) == expected


def test_layout_genes_puts_overlaps_on_new_rows() -> None:
    genes = [
        coverage_plot.Gene("chr1", 0, 100, "A"),
        coverage_plot.Gene("chr1", 50, 150, "B"),
        coverage_plot.Gene("chr1", 200, 260, "C"),
    ]
    placed, row_count = coverage_plot.layout_genes(genes)
    rows = {gene.name: row for gene, row in placed}
    assert row_count == 2
    assert rows["A"] != rows["B"]
    assert rows["C"] == rows["A"]


def test_resolve_gene_color_falls_back_to_palette() -> None:
    assert coverage_plot.resolve_gene_color("#ff8000", 0) == "#ff8000"
    assert coverage_plot.resolve_gene_color(None, 1) == coverage_plot.GENE_PALETTE[1]
    assert coverage_plot.resolve_gene_color("bogus", 2) == coverage_plot.GENE_PALETTE[2]
    assert coverage_plot.resolve_label_color("#3b8def") == "white"
    assert coverage_plot.resolve_label_color("#ffffff") == coverage_plot.TEXT_COLOR


# --------------------------------------------------------------- CLI 端到端
def test_cli_end_to_end(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "short",
            "coverage",
            "--bam",
            str(BAM),
            "--region",
            "chr1:1-2000",
            "--bins",
            "8",
            "--out-dir",
            str(tmp_path),
            "--name",
            "cov",
        ]
    )
    assert code == 0
    output = capsys.readouterr().out
    assert "93.5% covered" in output
    assert "29 reads" in output
    for name in ("cov.png", "cov.region.tsv", "cov.bins.tsv"):
        assert (tmp_path / name).exists(), name


def test_cli_default_name_comes_from_bam_and_region(tmp_path: Path) -> None:
    code = main(
        [
            "short",
            "coverage",
            "--bam",
            str(BAM),
            "--region",
            "chr1:1-2000",
            "--bins",
            "8",
            "--out-dir",
            str(tmp_path),
        ]
    )
    assert code == 0
    assert (tmp_path / "reads.chr1_1-2000.png").exists()
    assert (tmp_path / "reads.chr1_1-2000.bins.tsv").exists()


def test_cli_whole_contig_and_bases_table(tmp_path: Path) -> None:
    code = main(
        [
            "short",
            "coverage",
            "--bam",
            str(BAM),
            "--region",
            "chr1",
            "--out-dir",
            str(tmp_path),
            "--name",
            "whole",
            "--bases-table",
            str(tmp_path / "bases.tsv"),
        ]
    )
    assert code == 0
    _header, rows = tables.read_tsv(tmp_path / "bases.tsv")
    assert len(rows) == 2200  # 整条参考序列
    _header, region_rows = tables.read_tsv(tmp_path / "whole.region.tsv")
    assert region_rows[0]["region_end_1based"] == "2200"


def test_cli_list_references(capsys) -> None:
    assert main(["short", "coverage", "--bam", str(BAM), "--list"]) == 0
    out = capsys.readouterr().out
    assert out.splitlines()[0] == "reference\tlength\tmapped"
    assert out.splitlines()[1].startswith("chr1\t2200\t")


def test_cli_rejects_missing_bam(tmp_path: Path, capsys) -> None:
    assert main(["short", "coverage", "--bam", str(tmp_path / "nope.bam"), "--region", "chr1"]) == 2
    assert "不存在" in capsys.readouterr().err


def test_cli_rejects_unknown_reference(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "short",
            "coverage",
            "--bam",
            str(BAM),
            "--region",
            "chrZ:1-100",
            "--out-dir",
            str(tmp_path),
        ]
    )
    assert code == 2
    assert "不是该 BAM 的参考序列" in capsys.readouterr().err


def test_cli_rejects_empty_region(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "short",
            "coverage",
            "--bam",
            str(BAM),
            "--region",
            "chr1:2201-2300",
            "--out-dir",
            str(tmp_path),
        ]
    )
    assert code == 2
    assert "空区间" in capsys.readouterr().err


def test_cli_rejects_missing_bed(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "short",
            "coverage",
            "--bam",
            str(BAM),
            "--region",
            "chr1:1-2000",
            "--bed",
            str(tmp_path / "nope.bed"),
            "--out-dir",
            str(tmp_path),
        ]
    )
    assert code == 2
    assert "--bed 不存在" in capsys.readouterr().err
