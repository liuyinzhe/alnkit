"""``alnkit long breakpoints``（A5）：黄金回归（断点 / 类 / 类池 / segment_ends）。

黄金回归需要真实测序数据（不入公开仓库）：

    python tests/fixtures/build_private_longread.py
    export ALNKIT_TEST_DATA=~/.local/share/alnkit-testdata
    python -m pytest tests/test_long_breakpoints.py -m slow

``class_read_pools.tsv`` 里有一列是**类池 fastq 的绝对路径**（旧实现把 outdir 拼进去），
黄金文件在生成时把该前缀归一化成 ``<OUTDIR>``，测试同样归一化后再比，
其它文件一律逐字节比。
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from _fs import golden_files

from alnkit.cli import main

PRIVATE = (
    Path(os.environ["ALNKIT_TEST_DATA"]) / "longread"
    if os.environ.get("ALNKIT_TEST_DATA")
    else None
)
requires_data = pytest.mark.skipif(
    PRIVATE is None or not PRIVATE.exists(),
    reason="需要 ALNKIT_TEST_DATA（含真实 BAM/FASTQ 与黄金产物，不入公开仓库）",
)

REQUIRED_INPUTS = (
    "alignment_segments.tsv",
    "read_classification.tsv",
    "overlap_edges.tsv",
    "related_reads.fastq",
    "merged_reference.fa.fai",
)


def _normalize_pools(text: str, outdir: str) -> str:
    return text.replace(outdir, "<OUTDIR>")


@requires_data
@pytest.mark.slow
def test_golden_outputs_match(tmp_path: Path) -> None:
    assert PRIVATE is not None
    inputs = PRIVATE / "input"
    golden = PRIVATE / "golden" / "breakpoints"
    missing = [name for name in REQUIRED_INPUTS if not (inputs / name).exists()]
    if missing:
        pytest.fail(f"缺少输入 {missing}（运行 tests/fixtures/build_private_longread.py）")

    code = main(
        [
            "long",
            "breakpoints",
            "--segments",
            str(inputs / "alignment_segments.tsv"),
            "--classification",
            str(inputs / "read_classification.tsv"),
            "--overlaps",
            str(inputs / "overlap_edges.tsv"),
            "--reads",
            str(inputs / "related_reads.fastq"),
            "--merged-fai",
            str(inputs / "merged_reference.fa.fai"),
            "--outdir",
            str(tmp_path),
        ]
    )
    assert code == 0

    produced = golden_files(tmp_path)
    expected = golden_files(golden)
    assert produced == expected, "产物集合与旧脚本不同"

    for name in expected:
        got = (tmp_path / name).read_text(encoding="utf-8")
        want = (golden / name).read_text(encoding="utf-8")
        if name == "class_read_pools.tsv":
            got, want = _normalize_pools(got, str(tmp_path)), _normalize_pools(want, "<OUTDIR>")
        assert got == want, f"{name} 与旧脚本产物不再一致"


@requires_data
@pytest.mark.slow
def test_class_pools_point_into_outdir(tmp_path: Path) -> None:
    """类池表里的 fastq 路径必须指向本次 outdir（旧实现就是绝对路径）。"""
    assert PRIVATE is not None
    inputs = PRIVATE / "input"
    if not (inputs / "alignment_segments.tsv").exists():
        pytest.fail("缺输入（运行 tests/fixtures/build_private_longread.py）")
    code = main(
        [
            "long",
            "breakpoints",
            "--segments",
            str(inputs / "alignment_segments.tsv"),
            "--classification",
            str(inputs / "read_classification.tsv"),
            "--overlaps",
            str(inputs / "overlap_edges.tsv"),
            "--reads",
            str(inputs / "related_reads.fastq"),
            "--merged-fai",
            str(inputs / "merged_reference.fa.fai"),
            "--outdir",
            str(tmp_path),
        ]
    )
    assert code == 0
    rows = [
        line.split("\t")
        for line in (tmp_path / "class_read_pools.tsv").read_text(encoding="utf-8").splitlines()
        if not line.startswith("#")
    ][1:]
    assert rows, "类池表没有数据行"
    for row in rows:
        fastq = Path(row[-1])
        assert fastq.exists(), fastq
        assert str(tmp_path) in str(fastq)
