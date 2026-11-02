"""方案 C（借鉴 DEGAP 的边锚定招募）：悬垂解析、簇投票、重复守卫、端到端延伸。"""

from __future__ import annotations

import stat
from pathlib import Path

import pytest

from alnkit.longread import recruit

### CIGAR 走位 ###


def test_walk_cigar_offsets() -> None:
    assert recruit._walk_cigar("500M", 0) == 0
    assert recruit._walk_cigar("500M", 500) == 500
    assert recruit._walk_cigar("100S400M200M", 400) == 500  # 前导 soft clip 占 read 偏移
    assert recruit._walk_cigar("200M100D300M", 100) == 100  # 删除占参考不占 read
    assert recruit._walk_cigar("200M100I300M", 200) == 200  # 插入在 ref 200 之后
    assert recruit._walk_cigar("200M100I300M", 250) == 350  # 走过插入后 read 偏移多 100


### SAM 解析 ###


def _sam(flag: int, cigar: str, nm: int = 5, mapq: int = 60, pos: int = 1, name: str = "r1") -> str:
    fields = [
        name,
        str(flag),
        "edge",
        str(pos),
        str(mapq),
        cigar,
        "*",
        "0",
        "0",
        "*",
        "*",
        f"NM:i:{nm}",
    ]
    return "\t".join(fields)


def test_parse_alignment_basics() -> None:
    alignment = recruit.parse_alignment(_sam(0, "500M"), 500, 20, 0.99)
    assert alignment is not None
    assert (alignment.reference_start, alignment.reference_end) == (0, 500)
    assert alignment.identity == pytest.approx(0.99)


def test_parse_alignment_filters() -> None:
    assert recruit.parse_alignment(_sam(0x4, "500M"), 500, 20, 0.99) is None  # 未比对
    assert recruit.parse_alignment(_sam(0x100, "500M"), 500, 20, 0.99) is None  # secondary
    assert recruit.parse_alignment(_sam(0x800, "500M"), 500, 20, 0.99) is None  # supplementary
    assert recruit.parse_alignment(_sam(0, "500M", mapq=5), 500, 20, 0.99) is None  # 低 MAPQ
    assert recruit.parse_alignment(_sam(0, "500M", nm=50), 500, 20, 0.99) is None  # 一致率不够
    assert recruit.parse_alignment("@SQ\tSN:edge\tLN:500", 500, 20, 0.99) is None


def test_parse_alignment_soft_clip_hangover_coordinates() -> None:
    """右端：前导 0、尾部 soft clip 就是悬垂；左端：前导 soft clip 才是悬垂。"""
    right = recruit.parse_alignment(_sam(0, "500M300S"), 500, 20, 0.99)
    assert right is not None and right.query_end == 500
    left = recruit.parse_alignment(_sam(0, "300S500M"), 500, 20, 0.99)
    assert left is not None and left.query_start == 300


### 簇投票 ###


def test_pick_hangover_prefers_biggest_support() -> None:
    candidates = [
        ("a", "ACGT" * 100 + "A" * 100),  # 支持 2 条
        ("b", "ACGT" * 100 + "C" * 100),
        ("c", "TTTT" * 100 + "A" * 100),  # 支持 1 条
    ]
    chosen = recruit._pick_hangover(candidates, min_support=2, min_progress_bp=100)
    assert chosen is not None and chosen.startswith("ACGT")


def test_pick_hangover_rejects_thin_or_short() -> None:
    thin = [("a", "ACGT" * 100), ("b", "TTTT" * 100)]
    assert recruit._pick_hangover(thin, min_support=2, min_progress_bp=100) is None
    short = [("a", "ACGT" * 5), ("b", "ACGT" * 5)]
    assert recruit._pick_hangover(short, min_support=2, min_progress_bp=100) is None


### 重复守卫 ###


def test_repeat_guard_uses_anchor_depth_as_baseline() -> None:
    """本底用接缝锚点的出现条数；悬垂前缀远超本底才判重复。"""
    anchor = "GGCTTAGCAAGGCTTAGCAAGGCTTAGCAAGGCTTAGCAAGGCTTAGCAA"  # 50 bp，非回文
    hangover_motif = "TACCGATTGCTACCGATTGCTACCGATTGCTACCGATTGCTACCGATTGC"  # 50 bp
    hangover = hangover_motif + "CATTGCCATT" * 30

    baseline = {f"b{index}": "GAATTC" * 200 + anchor for index in range(10)}
    for index in range(10):
        baseline[f"h{index}"] = hangover_motif + "GCTAAG" * 150
    assert recruit._looks_repetitive(hangover, anchor, baseline) is False

    repeated = dict(baseline)
    for index in range(100):
        repeated[f"x{index}"] = hangover_motif + "TTGCAA" * 150
    assert recruit._looks_repetitive(hangover, anchor, repeated) is True


### 端到端：假 minimap2 ###


@pytest.fixture()
def stub_minimap2(tmp_path: Path):
    def build(cigar: str, names: tuple[str, ...] = ("r1", "r2")) -> Path:
        lines = "".join(_sam(0, cigar, nm=0, name=name) + "\n" for name in names)
        script = tmp_path / "minimap2"
        script.write_text("#!/bin/sh\n" + "cat <<'EOF'\n" + lines + "EOF\n", encoding="utf-8")
        script.chmod(script.stat().st_mode | stat.S_IEXEC)
        return script

    return build


@pytest.fixture()
def reads_setup(tmp_path: Path):
    probe = "ACGT" * 125  # 500 bp 探针序列
    reads = {"r1": probe + "TTGG" * 75, "r2": probe + "TTGG" * 75}
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "".join(f"@{name}\n{seq}\n+\n{'I' * len(seq)}\n" for name, seq in reads.items())
    )
    return probe, reads, fastq


def test_extend_right_end(stub_minimap2, reads_setup) -> None:
    probe, reads, fastq = reads_setup
    minimap2 = stub_minimap2("500M300S")
    extended, right, left = recruit.extend_contig_ends(
        probe, fastq, reads, str(minimap2), recruit.RecruitSettings(max_rounds=1)
    )
    assert right == 1 and left == 0
    assert extended == probe + "TTGG" * 75  # 300 bp 悬垂拼在右端


def test_extend_left_end(stub_minimap2, reads_setup) -> None:
    probe, reads, fastq = reads_setup
    # 左侧悬垂是“探针之前”的 300 bp：CIGAR 300S500M
    reads_left = {"r1": "TTGG" * 75 + probe, "r2": "TTGG" * 75 + probe}
    minimap2 = stub_minimap2("300S500M")
    extended, right, left = recruit.extend_contig_ends(
        probe, fastq, reads_left, str(minimap2), recruit.RecruitSettings(max_rounds=1)
    )
    assert right == 0 and left == 1
    assert extended == "TTGG" * 75 + probe


def test_extend_stops_when_progress_too_small(stub_minimap2, reads_setup) -> None:
    probe, _reads, fastq = reads_setup
    reads = {"r1": probe + "TT" * 10, "r2": probe + "TT" * 10}  # 悬垂只有 20 bp
    minimap2 = stub_minimap2("500M20S")
    extended, right, left = recruit.extend_contig_ends(
        probe, fastq, reads, str(minimap2), recruit.RecruitSettings(max_rounds=3)
    )
    assert (right, left) == (0, 0)
    assert extended == probe


def test_extend_stops_when_support_too_thin(stub_minimap2, reads_setup) -> None:
    probe, _reads, fastq = reads_setup
    minimap2 = stub_minimap2("500M300S", names=("r1",))  # 只有一条 read 支持
    single = {"r1": probe + "TTGG" * 75}
    extended, right, _left = recruit.extend_contig_ends(
        probe, fastq, single, str(minimap2), recruit.RecruitSettings(max_rounds=1, min_support=2)
    )
    assert right == 0
    assert extended == probe
