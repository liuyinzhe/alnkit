"""A1 ``ref merge`` / A2 ``align reads``：与旧 shell 产物等价。

这两步原来是 shell，验收口径与三代分析步骤不同：

* **A1 可逐字节对账**：合并是纯文本操作。黄金收编了旧 shell 的产物
  （`reference_headers.tsv` / `exogenous_uppercased.fa` / 两个 `.fai`），
  2.3 GB 的 `merged_reference.fa` 只记 md5。
* **A2 不能逐字节对账**：minimap2 把命令行（含输出目录路径）写进 BAM 的 `@PG` 行，
  换个 outdir 就不可能同字节。因此比的是**比对记录**的 md5
  （`samtools view -F 0x900`，排除 secondary/supplementary 以外的 header 噪音）——
  记录一致就说明比对结果一致。

需要私有数据与外部命令（minimap2 / samtools）：

    python tests/fixtures/build_private_longread.py
    export ALNKIT_TEST_DATA=~/.local/share/alnkit-testdata
    python -m pytest tests/test_prepare.py -m slow
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
from pathlib import Path

import pytest

from alnkit import prepare

PRIVATE = (
    Path(os.environ["ALNKIT_TEST_DATA"]) / "longread"
    if os.environ.get("ALNKIT_TEST_DATA")
    else None
)
requires_data = pytest.mark.skipif(
    PRIVATE is None or not PRIVATE.exists(),
    reason="需要 ALNKIT_TEST_DATA（含真实数据与旧 shell 黄金，不入公开仓库）",
)


def _source_dir() -> Path:
    return Path(os.environ.get("ALNKIT_UPSTREAM_LONGREAD", "."))


@pytest.mark.slow
@requires_data
def test_ref_merge_matches_shell_outputs(tmp_path: Path) -> None:
    assert PRIVATE is not None
    inputs = PRIVATE / "input"
    # 用私有 fixture 里的参考与外源（文件名已脱敏；contig 名来自 fasta header）
    reference = inputs / "reference.fa"
    exog_h, exog_l = inputs / "exog/insertA.fa", inputs / "exog/insertB.fa"
    assert PRIVATE is not None
    golden = PRIVATE / "golden/refmerge"
    missing = [str(path) for path in (reference, exog_h, exog_l) if not path.exists()]
    if missing:
        pytest.fail(f"缺输入 {missing}")

    result = prepare.merge_reference(
        prepare.MergeSettings(
            reference=str(reference), exog=(str(exog_h), str(exog_l)), outdir=str(tmp_path)
        )
    )
    assert result["sequences"] == 649

    for name in (
        "reference_headers.tsv",
        "exogenous_uppercased.fa",
        "reference.fa.fai",
        "merged_reference.fa.fai",
    ):
        assert (tmp_path / name).read_bytes() == (golden / name).read_bytes(), (
            f"{name} 与 shell 不同"
        )

    expected_md5 = (golden / "merged_reference.md5").read_text(encoding="utf-8").strip()
    actual_md5 = hashlib.md5((tmp_path / "merged_reference.fa").read_bytes()).hexdigest()
    assert actual_md5 == expected_md5, "merged_reference.fa 与 shell 产物不同"


@pytest.mark.slow
@requires_data
def test_align_reads_records_match(tmp_path: Path) -> None:
    """比对记录必须与旧 shell 的结果一致（BAM 字节因 @PG 命令行含路径而不可比）。"""
    if shutil.which("minimap2") is None or shutil.which("samtools") is None:
        pytest.fail("缺 minimap2/samtools：装好后重跑")
    assert PRIVATE is not None
    golden = PRIVATE / "golden/align/records.md5"
    if not golden.exists():
        pytest.fail("缺 A2 黄金：跑 tests/fixtures/build_private_longread.py")

    merged = tmp_path / "merged_reference.fa"
    # 注意必须用**合并参考**（参考 + 外源），不是纯参考：比对目标不同，记录就不同
    source_merged = PRIVATE / "input/merged_reference.fa"
    if not source_merged.exists():
        pytest.fail("缺 input/merged_reference.fa（跑 builder）")
    shutil.copy2(Path(os.path.realpath(source_merged)), merged)
    shutil.copy2(Path(str(os.path.realpath(source_merged)) + ".fai"), str(merged) + ".fai")

    reads = PRIVATE / "input/input.fastq"
    result = prepare.align_reads(
        prepare.AlignSettings(
            reference=str(merged), reads=str(reads), outdir=str(tmp_path), threads=8
        )
    )
    bam = Path(str(result["bam"]))
    assert bam.exists() and Path(str(result["index"])).exists()

    records = subprocess.run(
        ["samtools", "view", "-F", "0x900", str(bam)], capture_output=True, check=True
    )
    actual = hashlib.md5(records.stdout).hexdigest()
    expected = golden.read_text(encoding="utf-8").strip()
    assert actual == expected, "比对记录与旧 shell 结果不同"


@requires_data
def test_prepare_rejects_missing_inputs(tmp_path: Path) -> None:
    from alnkit.core.errors import UsageError

    with pytest.raises(UsageError, match="找不到参考 FASTA"):
        prepare.merge_reference(
            prepare.MergeSettings(
                reference=str(tmp_path / "nope.fa"), exog=("x.fa",), outdir=str(tmp_path)
            )
        )
    dummy = tmp_path / "ref.fa"
    dummy.write_text(">c\nACGT\n", encoding="utf-8")
    with pytest.raises(UsageError, match="至少给一个"):
        prepare.merge_reference(prepare.MergeSettings(reference=str(dummy), outdir=str(tmp_path)))
