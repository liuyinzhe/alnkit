"""同源位点诊断：BED 解析、窗口回比、位点标注、写表。"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from alnkit.core import tables
from alnkit.shortread import homology


def test_load_bed_regions_basic(tmp_path: Path) -> None:
    bed = tmp_path / "known.bed"
    bed.write_text(
        "# comment\ntrack name=known\nchr6\t39048000\t39087000\tGLP1R\nchr1\t100\t200\n",
        encoding="utf-8",
    )
    regions = homology.load_bed_regions([bed])
    assert len(regions) == 2
    assert regions[0].chrom == "chr6"
    assert (regions[0].start_0based, regions[0].end_0based) == (39048000, 39087000)
    assert regions[0].source == "bed:GLP1R"
    assert regions[0].identity is None
    assert regions[1].source == "bed:known"  # 没有第 4 列时用文件名


def test_load_bed_regions_rejects_bad_input(tmp_path: Path) -> None:
    bad = tmp_path / "bad.bed"
    bad.write_text("chr1\t100\n", encoding="utf-8")
    from alnkit.core.errors import UsageError

    with pytest.raises(UsageError):
        homology.load_bed_regions([bad])
    bad.write_text("chr1\t200\t100\n", encoding="utf-8")
    with pytest.raises(UsageError):
        homology.load_bed_regions([bad])


def test_count_cs_differences() -> None:
    """cs 标签里数替换与 indel；内含子（~）不计入。"""
    # `+tt` 是 2 个插入碱基；`~gt…ag` 是内含子，不计入
    assert homology._count_cs_differences("=ACGT*ag+tt=GG~gt1234ag=TT") == (1, 2)
    assert homology._count_cs_differences("=AA*ca=+(1)") == (1, 0)  # 括号形式不计（非 cs=long）
    assert homology._count_cs_differences("=ACGTACGT") == (0, 0)


def test_annotate_sites_marks_inside_and_outside() -> None:
    """位点落在区域内（含 50 bp 放宽）标 yes，否则 no。"""
    regions = [homology.HomologyRegion("chr6", 1000, 2000, "insert:X", 0.99)]
    rows = homology.annotate_sites([("chr6", 1500), ("chr6", 3000), ("chr1", 1500)], regions)
    assert [row["is_homology_site"] for row in rows] == ["yes", "no", "no"]
    assert rows[0]["best_region_start_1based"] == 1001  # 1-based 闭区间
    assert rows[0]["best_region_end_1based"] == 2000
    assert rows[0]["best_region_identity"] == 0.99
    assert rows[1]["region_sources"] == "-"


def test_annotate_sites_padding_boundary() -> None:
    """边界外 50 bp 内仍算同源（断点常落在同源区边界外）。"""
    regions = [homology.HomologyRegion("chr6", 1000, 2000, "bed:X", None)]
    assert homology.annotate_sites([("chr6", 960)], regions)[0]["is_homology_site"] == "yes"
    assert homology.annotate_sites([("chr6", 940)], regions)[0]["is_homology_site"] == "no"


def _write_fasta(path: Path, records: dict[str, str]) -> Path:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n{sequence}\n")
    return path


@pytest.fixture()
def stub_minimap2(tmp_path: Path):
    """假 minimap2：从 stdin 参数里取窗口与插入 FASTA，直接吐一行 PAF。"""

    def build(qspan: int = 400) -> Path:
        script = tmp_path / "minimap2"
        script.write_text(
            "#!/bin/sh\n"
            "window=$(eval echo \\${$#})\n"
            "cat <<'EOF'\n"
            f"GLP1R-FLAG\t400\t0\t{qspan}\t+\twindow\t10000\t1000\t2000\t{qspan}\t{qspan}\t60\tcs:Z:=ACGT\n"
            "EOF\n",
            encoding="utf-8",
        )
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
        return script

    return build


def test_local_homology_regions_with_stub(tmp_path: Path, stub_minimap2) -> None:
    """窗口回比：把 PAF 坐标换算回基因组坐标，并算出插入覆盖率。"""
    reference = _write_fasta(tmp_path / "ref.fa", {"chr6": "A" * 10000, "GLP1R-FLAG": "A" * 400})
    minimap2 = stub_minimap2(qspan=400)
    regions = homology.local_homology_regions(
        reference,
        ["GLP1R-FLAG"],
        [("chr6", 5000)],
        str(minimap2),
        window_bp=5000,
        min_identity=0.90,
        min_length=200,
        threads=1,
    )
    assert len(regions) == 1
    region = regions[0]
    assert region.source == "insert:GLP1R-FLAG"
    assert (region.start_0based, region.end_0based) == (1000, 2000)  # 窗口内偏移 + 窗口起点
    assert region.query_coverage == 1.0


def test_write_site_homology_table(tmp_path: Path) -> None:
    rows = homology.annotate_sites(
        [("chr6", 1500)], [homology.HomologyRegion("chr6", 1000, 2000, "insert:X", 0.99, 0.95)]
    )
    output = tmp_path / "s.site_homology.tsv"
    written = homology.write_site_homology_table(output, rows, "s")
    assert written == 1
    body = [
        line for line in output.read_text(encoding="utf-8").splitlines() if not line.startswith("#")
    ]
    assert len(body[0].split("\t")) == len(tables.COLUMNS["site_homology"])
    assert body[-1].split("\t")[0] == "s"
