"""黄金文件的出身校验（防止「自己证明自己」）。

背景：私有黄金由 ``tests/fixtures/build_private_*.py`` 生成。这些脚本支持
``--trust-current``（用当前实现生成黄金），一旦被用来重建，A3/A5/A6/B1 的黄金回归
就退化成"当前实现和自己比"，而其余测试照样全绿。2026-11-01 的独立复审把这条列为
最严重的问题，本文件就是它的守卫。

没有私有数据时（``ALNKIT_TEST_DATA`` 未设置）整个文件跳过。
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

DATA = Path(os.environ["ALNKIT_TEST_DATA"]) if os.environ.get("ALNKIT_TEST_DATA") else None
pytestmark = pytest.mark.skipif(
    DATA is None or not DATA.exists(),
    reason="需要 ALNKIT_TEST_DATA（私有黄金目录，不入公开仓库）",
)

GOLDEN_DIRS = (
    "shortread/golden/owned",
    "shortread/golden/markdup",
    "shortread/golden/consensus",
    "longread/golden/classify",
    "longread/golden/breakpoints",
    "longread/golden/assemble",
)


def _provenance(directory: Path) -> dict[str, str]:
    path = directory / "PROVENANCE.txt"
    if not path.exists():
        return {}
    fields: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        key, _, value = line.partition("=")
        fields[key.strip()] = value.strip()
    return fields


@pytest.mark.parametrize("relative", GOLDEN_DIRS)
def test_golden_has_provenance(relative: str) -> None:
    assert DATA is not None
    directory = DATA / relative
    if not directory.exists():
        pytest.fail(f"{directory} 不存在：先跑 tests/fixtures/build_private_*.py")
    fields = _provenance(directory)
    assert fields, f"{directory} 缺 PROVENANCE.txt（无法证明黄金来自旧脚本）"
    for key in ("generator", "generated_at", "command"):
        assert key in fields, f"{relative}/PROVENANCE.txt 缺字段 {key}"


@pytest.mark.parametrize("relative", GOLDEN_DIRS)
def test_golden_was_generated_by_upstream(relative: str) -> None:
    """黄金必须由旧脚本生成；`generator=alnkit` 说明回归已经变成自证。"""
    assert DATA is not None
    directory = DATA / relative
    if not directory.exists():
        pytest.fail(f"{directory} 不存在：先跑 tests/fixtures/build_private_*.py")
    fields = _provenance(directory)
    if not fields:
        pytest.fail(f"{relative} 缺 PROVENANCE.txt")
    assert fields["generator"] == "upstream", (
        f"{relative} 的黄金是 {fields['generator']} 生成的（--trust-current？），"
        "黄金回归已失效；请用旧脚本重建"
    )
    assert fields.get("script_sha256"), f"{relative} 的 PROVENANCE 缺旧脚本 sha256，无法核对来源"


def test_stage_inputs_exist() -> None:
    """黄金能用还不够，输入也得齐（早先 A5 的输入压根没被 builder 创建）。"""
    assert DATA is not None
    required = [
        "longread/input/alignment_segments.tsv",
        "longread/input/read_classification.tsv",
        "longread/input/related_reads.fastq",
        "longread/input/overlap_edges.tsv",
        "longread/input/overlaps.paf",
        "longread/input/classes",
        "shortread/bam",
    ]
    missing = [name for name in required if not (DATA / name).exists()]
    assert not missing, f"缺输入 {missing}：跑 tests/fixtures/build_private_*.py"
