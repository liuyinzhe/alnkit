"""`long assemble` 的 A/B 两项增量：支撑度报告（A）与接缝支持门槛（B）。"""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from alnkit.longread import assemble, support

### 接缝支持数（反链坐标回归） ###


@pytest.mark.parametrize("direction", [1, -1])
def test_junction_support_counts_both_strands(direction: int) -> None:
    """placed 里的 start 恒为定向序列首碱基坐标，与 +/− 无关（曾经的 bug）。"""
    sequences = {"r1": "A" * 1000}
    placed = {"r1": (500, direction)}
    # 接缝落在 r1 的 [500, 1500) 内 → r1 覆盖 + 本次 read 自己 = 2
    assert assemble.count_junction_support(placed, sequences, 900, 100) == 2
    # 接缝远在右侧（r1 只到 1499，窗口左侧 1350 仍相交 → 仍算覆盖）
    assert assemble.count_junction_support(placed, sequences, 1450, 100) == 2
    # 接缝在左侧 100 bp 外 → 只有本次 read 自己
    assert assemble.count_junction_support(placed, sequences, 300, 100) == 1


def test_junction_support_sums_neighbours() -> None:
    """邻居越多支持数越大；窗口变大能圈进更远的 read。"""
    sequences = {"r1": "A" * 1000, "r2": "A" * 1000, "r3": "A" * 2000}
    placed = {"r1": (500, 1), "r2": (600, -1), "r3": (0, 1)}
    # 窗口 ±100 覆盖到 r1/r2/r3 三者 → 3 邻居 + 本次 = 4
    assert assemble.count_junction_support(placed, sequences, 900, 100) == 4
    # 接缝在 r3 右端外（r3 只到 1999）→ 没有邻居，只剩本次 read
    assert assemble.count_junction_support(placed, sequences, 2500, 100) == 1
    # 窗口放宽到 ±800 才圈到 r3 → 1 邻居 + 本次 = 2
    assert assemble.count_junction_support(placed, sequences, 2500, 800) == 2


def test_junction_depth_ratio_gate() -> None:
    """相对深度门槛：调高到不可能满足时不再延伸，类里只剩种子自己一条 contig。"""
    sequences = {"j1": "ACGT" * 500, "v1": "TGCA" * 500}
    blocks = {
        ("j1", "v1"): {
            "query": "v1",
            "target": "j1",
            "strand": "+",
            "block_bp": 500,
            "identity": 1.0,
            "query_start": 0,
            "query_end": 500,
            "target_start": 1500,
            "target_end": 2000,
        }
    }
    loose = assemble.assemble_one_class(
        sequences, {"j1"}, blocks, agreement_bp=200, min_junction_depth_ratio=0.0
    )
    assert len(loose) == 1
    assert len(loose[0]["placed"]) == 2

    statistics = assemble.collections.Counter()
    strict = assemble.assemble_one_class(
        sequences,
        {"j1"},
        blocks,
        agreement_bp=200,
        min_junction_depth_ratio=99.0,
        statistics=statistics,
    )
    assert statistics["rejected_by_depth_ratio"] > 0
    assert [len(contig["placed"]) for contig in strict] == [1, 1]


def test_junction_depth_ratio_is_relative() -> None:
    """比值 = 接缝邻居数 ÷ 已放置 read 的平均深度，而不是绝对条数。"""
    sequences = {"r1": "A" * 1000, "r2": "A" * 1000}
    # 两条 read 完全重叠 → 平均深度 2×，接缝处邻居 2 → 比值约 1.0
    placed = {"r1": (0, 1), "r2": (0, 1)}
    ratio = assemble.junction_depth_ratio(placed, sequences, 500, 100)
    assert ratio == pytest.approx(1.0)
    # 只有种子时平均深度 1×、邻居 1 → 仍是 1.0（尺度无关）
    assert assemble.junction_depth_ratio({"r1": (0, 1)}, sequences, 500, 100) == 1.0


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
    body = output.read_text(encoding="utf-8").splitlines()
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


def test_parse_contigs_survives_bad_fasta(tmp_path: Path) -> None:
    """重复 contig 名不清零、空 header 不抛 IndexError（防御性）。"""
    fasta = tmp_path / "contigs.fasta"
    fasta.write_text(">ctg1\nACGT\n>ctg1\nACGTACGT\n>\nACGT\n", encoding="utf-8")
    lengths = support._parse_contigs(fasta)
    assert lengths["ctg1"] == 12  # 同名累加，不是被第二条覆盖成 8
    assert lengths.get("", 0) == 0  # 空 header 的记录直接忽略，不抛异常


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
