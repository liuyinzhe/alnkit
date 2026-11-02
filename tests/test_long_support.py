"""`long assemble` 的组装支撑度报告（A 方案）。"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from alnkit.longread import support

### 支撑度统计 ###


def test_coverage_stats_normalises_last_partial_window() -> None:
    """末窗口按实际长度归一化，满覆盖的短 contig 不应报断层。"""
    contig_length = 250
    intervals = [(0, contig_length)] * 2  # 2 条 read 满覆盖
    covered_bp, fraction, mean_depth, depths = support.coverage_stats(contig_length, intervals)
    assert covered_bp == contig_length
    assert fraction == pytest.approx(1.0)
    assert mean_depth == pytest.approx(2.0)
    assert depths == pytest.approx([2.0, 2.0])  # 200 bp 与 50 bp 两个窗口都是 2×
    assert support._count_cliffs(depths) == 0


def test_coverage_stats_merges_overlaps_for_covered_bp() -> None:
    """covered_bp 去重、mean_depth 按比对碱基累加（与 samtools depth 同义）。"""
    covered_bp, fraction, mean_depth, depths = support.coverage_stats(100, [(0, 60), (50, 100)])
    assert covered_bp == 100
    assert fraction == pytest.approx(1.0)
    assert mean_depth == pytest.approx(1.1)
    assert depths == pytest.approx([1.1])


def test_count_cliffs_positive_case() -> None:
    """深度骤降要能被数出来（此前只测了“不断层”）。"""
    assert support._count_cliffs([4.0, 0.5]) == 1
    assert support._count_cliffs([4.0, 3.0, 0.5]) == 1
    assert support._count_cliffs([1.5, 0.1]) == 0  # 前窗口 <2× 不算（噪声）


def test_status_thresholds() -> None:
    assert support._status(0.0, 0.0, 0, "no") == "no_reads"
    assert support._status(1.0, 1.0, 0, "yes") == "low_coverage"
    assert support._status(0.5, 50.0, 0, "yes") == "low_coverage"
    assert support._status(1.0, 50.0, support.CLIFF_MIN_COUNT, "yes") == "chimera_suspect"
    assert support._status(1.0, 50.0, 0, "no") == "unsupported_ends"
    assert support._status(1.0, 50.0, 0, "yes") == "ok"


### SAM 解析 ###


def _sam_line(flag: int, name: str = "r1", pos: int = 101, cigar: str = "50M100D50M") -> str:
    fields = [name, str(flag), "ctg1", str(pos), "60", cigar, "*", "0", "0", "*", "*"]
    return "\t".join(fields)


def test_parse_sam_line_filters_secondary_keeps_supplementary() -> None:
    assert support._parse_sam_line(_sam_line(0)) == ("r1", "ctg1", 100, 200, 60)
    assert support._parse_sam_line(_sam_line(0x100)) is None  # secondary
    assert support._parse_sam_line(_sam_line(0x4)) is None  # unmapped
    assert support._parse_sam_line("@SQ\tSN:ctg1\tLN:1000") is None
    supplementary = support._parse_sam_line(_sam_line(0x800, cigar="100S50M"))
    assert supplementary == ("r1", "ctg1", 100, 50, 60)


### 端到端：假 minimap2 打桩，验证 collect_support 落表 ###


@pytest.fixture()
def stub_minimap2(tmp_path: Path) -> Path:
    script = tmp_path / "minimap2"
    script.write_text(
        "#!/bin/sh\n"
        "cat <<'EOF'\n"
        "@SQ\tSN:ctg1\tLN:400\n"
        + _sam_line(0, "r1", 1, "400M")
        + "\n"
        + _sam_line(0, "r2", 1, "400M")
        + "\n"
        + _sam_line(0x100, "r3", 1, "400M")
        + "\n"
        + "EOF\n",
        encoding="utf-8",
    )
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    return script


def test_collect_support_ignores_secondary(tmp_path: Path, stub_minimap2: Path) -> None:
    contigs = tmp_path / "contigs.fasta"
    contigs.write_text(">ctg1\n" + "ACGT" * 100 + "\n", encoding="utf-8")
    reads = tmp_path / "reads.fastq"
    reads.write_text("".join(f"@r{index}\nACGT\n+\nIIII\n" for index in range(3)), encoding="utf-8")

    rows = support.collect_support("C01", contigs, reads, str(stub_minimap2), window_bp=200)

    assert len(rows) == 1
    row = rows[0]
    assert row.n_reads_total == 3
    assert row.n_reads_mapped == 2  # r3 是 secondary，被排除
    assert row.covered_fraction == pytest.approx(1.0)
    assert row.mean_depth == pytest.approx(2.0)
    assert row.depth_p10 == pytest.approx(2.0)
    assert row.status == "ok"

    output = tmp_path / "assembly_support.tsv"
    support.write_support_report(rows, output)
    body = [
        line for line in output.read_text(encoding="utf-8").splitlines() if not line.startswith("#")
    ]
    assert len(body[0].split("\t")) == 12  # 表头必须与 assembly_support 的 12 列一致
    assert body[-1].split("\t")[0] == "C01"
    assert body[-1].split("\t")[-1] == "ok"


def test_collect_support_reports_stderr(tmp_path: Path) -> None:
    """minimap2 失败时把 stderr 一起带进错误消息（否则用户看不到原因）。"""
    failing = tmp_path / "minimap2"
    failing.write_text("#!/bin/sh\necho 'failed to open file' >&2\nexit 1\n", encoding="utf-8")
    failing.chmod(failing.stat().st_mode | stat.S_IEXEC)
    contigs = tmp_path / "contigs.fasta"
    contigs.write_text(">ctg1\nACGT\n", encoding="utf-8")
    reads = tmp_path / "reads.fastq"
    reads.write_text("@r1\nACGT\n+\nIIII\n", encoding="utf-8")

    with pytest.raises(Exception) as excinfo:
        support.collect_support("C01", contigs, reads, str(failing))
    assert "failed to open file" in str(excinfo.value)
    assert os.fspath(contigs) in str(excinfo.value)


def test_parse_contigs_drops_duplicate_names(tmp_path: Path) -> None:
    """重名 contig 只保留第一条（重名会让 minimap2 把两条序列混进同一个桶）。"""
    fasta = tmp_path / "contigs.fasta"
    fasta.write_text(">ctg1\nACGT\n>ctg1\nACGTACGT\n>\nACGT\n", encoding="utf-8")
    lengths = support._parse_contigs(fasta)
    assert lengths["ctg1"] == 4  # 只保留第一条，不把两条长度相加
    assert lengths.get("", 0) == 0  # 空 header 的记录直接忽略，不抛异常


def test_covered_bp_merges_overlaps() -> None:
    """窗口覆盖率要按“合并后的覆盖”算：3 条 read 各盖同一段 20 bp 仍是 20 bp。"""
    assert support._covered_bp([(0, 20)] * 3, 0, 100) == 20
    assert support._covered_bp([(0, 60), (40, 100)], 0, 100) == 100


def test_ends_covered_requires_real_coverage() -> None:
    """末端只被重复覆盖一小段时报 no，而不是被 3 条 read 的重复计数灌成 yes。"""
    intervals = [(0, 20)] * 3
    rows = support.coverage_stats(2000, intervals)
    assert rows[0] == 20
    # 直接验证 status 判据用的窗口覆盖量
    assert support._covered_bp(intervals, 0, 100) == 20  # < 50 → 不算撑住


def test_class_names_accepts_single_string(tmp_path: Path, stub_minimap2: Path) -> None:
    """库调用方传字符串时按“一个类”处理，不会把 "C01" 拆成字符集。"""
    assembly = tmp_path / "assembly"
    (assembly / "C01" / "overlap_extension").mkdir(parents=True)
    (assembly / "C01" / "overlap_extension/contigs.fasta").write_text(
        ">ctg1\n" + "ACGT" * 100 + "\n", encoding="utf-8"
    )
    classes = tmp_path / "classes" / "C01"
    classes.mkdir(parents=True)
    reads = classes / "reads.fastq"
    reads.write_text("".join(f"@r{i}\nACGT\n+\nIIII\n" for i in range(3)), encoding="utf-8")

    rows = support.collect_all_support(
        assembly, tmp_path / "classes", str(stub_minimap2), class_names="C01"
    )
    assert [row.class_id for row in rows] == ["C01"]


def test_parse_sam_line_skips_zero_position() -> None:
    """POS=0 的坏行直接跳过，不抛 ValueError。"""
    assert support._parse_sam_line(_sam_line(0, pos=0)) is None
