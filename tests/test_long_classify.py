"""``alnkit long classify``（A3）：纯逻辑定点断言 + 私有数据黄金回归。

黄金回归需要真实测序数据（不入公开仓库）：

    python tests/fixtures/build_private_longread.py
    export ALNKIT_TEST_DATA=~/.local/share/alnkit-testdata
    python -m pytest tests/test_long_classify.py -m slow

没有数据时相关用例 skip（带原因），CI 默认走这条路径。
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
from _fs import golden_files

from alnkit.cli import main
from alnkit.longread import classify

PRIVATE = (
    Path(os.environ["ALNKIT_TEST_DATA"]) / "longread"
    if os.environ.get("ALNKIT_TEST_DATA")
    else None
)
requires_data = pytest.mark.skipif(
    PRIVATE is None or not PRIVATE.exists(),
    reason="需要 ALNKIT_TEST_DATA（含真实 BAM/FASTQ 与黄金产物，不入公开仓库）",
)


def _segment(**overrides) -> dict:
    base = {
        "read_name": "r1",
        "chrom": "chr1",
        "chrom_length": 1000,
        "ref_start": 200,
        "ref_end": 300,
        "read_start": 0,
        "read_end": 100,
        "strand": "+",
        "mapq": 60,
        "ref_aln_len": 100,
        "read_aln_len": 100,
        "cigar": "100M",
    }
    return {**base, **overrides}


def _manifest() -> dict:
    """私有数据里的真实序列名只存在本机 manifest 中（仓库里不出现）。"""
    assert PRIVATE is not None
    path = PRIVATE / "input/manifest.json"
    if not path.exists():
        pytest.fail(f"缺 {path}：跑 tests/fixtures/build_private_longread.py")
    return json.loads(path.read_text(encoding="utf-8"))


def _settings(**overrides) -> classify.ClassifySettings:
    base = {
        "exog": ("insertA",),
        "min_mapq": 20,
        "min_exog_block_bp": 100,
        "contig_edge_bp": 100,
        "short_genome_block_bp": 500,
    }
    return classify.ClassifySettings(**{**base, **overrides})


# ------------------------------------------------------------ 逐段判定
def test_exog_evidence_needs_length_and_mapq() -> None:
    assert classify.classify_segment(
        _segment(chrom="insertA", ref_aln_len=200), ("insertA",), _settings()
    ) == ("exog", True, "-")
    assert classify.classify_segment(
        _segment(chrom="insertA", ref_aln_len=50), ("insertA",), _settings()
    ) == ("exog", False, "exog_block_too_short")
    assert classify.classify_segment(
        _segment(chrom="insertA", ref_aln_len=200, mapq=10), ("insertA",), _settings()
    ) == ("exog", False, "low_mapq")


def test_genome_evidence() -> None:
    assert classify.classify_segment(_segment(), ("insertA",), _settings()) == ("genome", True, "-")
    assert classify.classify_segment(_segment(mapq=10), ("insertA",), _settings()) == (
        "genome",
        False,
        "low_mapq",
    )


def test_short_block_near_contig_edge_is_suspicious() -> None:
    # 落在 contig 两端 100 bp 内、且长度 <= 500 bp 的基因组段判为可疑错配
    edge = _segment(ref_start=5, ref_end=105, ref_aln_len=100)
    assert classify.classify_segment(edge, ("insertA",), _settings()) == (
        "genome",
        False,
        "contig_edge_short_block",
    )
    # 同样的短块，但离两端足够远 -> 正常保留
    middle = _segment(ref_start=200, ref_end=300, ref_aln_len=100)
    assert classify.classify_segment(middle, ("insertA",), _settings())[1] is True


def test_estimate_identity() -> None:
    segment = {"match_columns": 100, "nm": 2, "inserted_bases": 1, "deleted_bases": 0}
    assert classify.estimate_identity(segment) == pytest.approx(0.99)
    assert classify.estimate_identity({"match_columns": 0}) == 1.0


# ------------------------------------------------------------------- CLI
def test_cli_requires_inputs(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(["long", "classify", "--outdir", str(tmp_path)])
    assert excinfo.value.code == 2


# -------------------------------------------------------- 私有数据黄金回归
@requires_data
@pytest.mark.slow
def test_golden_outputs_byte_identical(tmp_path: Path) -> None:
    assert PRIVATE is not None
    inputs, golden = PRIVATE / "input", PRIVATE / "golden" / "classify"
    missing = [
        name
        for name in ("reads_to_merged.sorted.bam", "input.fastq", "merged_reference.fa.fai")
        if not (inputs / name).exists()
    ]
    if missing:
        pytest.fail(f"缺输入 {missing}：跑 tests/fixtures/build_private_longread.py")

    code = main(
        [
            "long",
            "classify",
            "--bam",
            str(inputs / "reads_to_merged.sorted.bam"),
            "--reads",
            str(inputs / "input.fastq"),
            "--merged-fai",
            str(inputs / "merged_reference.fa.fai"),
            "--exog",
            _manifest()["exog_names"][0],
            "--exog",
            _manifest()["exog_names"][1],
            "--outdir",
            str(tmp_path),
        ]
    )
    assert code == 0

    produced = golden_files(tmp_path)
    expected = golden_files(golden)
    assert produced == expected, "产物集合与旧脚本不同"

    for name in expected:
        assert (tmp_path / name).read_bytes() == (golden / name).read_bytes(), (
            f"{name} 与旧脚本产物不再一致"
        )
