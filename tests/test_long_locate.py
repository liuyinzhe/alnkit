"""``alnkit long locate``（A7）：落点/坐标还原的黄金回归。

与 A3/A5/A6 相同的口径，但有两条现实约束：

1. **PDF 不可逐字节比**：matplotlib 会在 PDF 里写入生成时间，旧脚本自己连跑两次也不一致
   （实测）。因此 PDF 只断言「存在 + 字节数与黄金相同」；
2. 文本产物里有 outdir 的绝对路径（`class_position_summary.tsv`、`summary.md`），
   黄金在生成时归一化成 `<OUTDIR>`，测试同样归一化后比。

其余产物（含 7 张 PNG、`class_assembly_positions.tsv`、`contig_map.tsv`、`work/*.paf`）逐字节比。

需要私有数据：``export ALNKIT_TEST_DATA=~/.local/share/alnkit-testdata``，
并先跑 ``python tests/fixtures/build_private_longread.py``（内置 2.3 GB 参考的软链接）。
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest
from _fs import golden_files, output_files

from alnkit.cli import main
from alnkit.core import tables

PRIVATE = (
    Path(os.environ["ALNKIT_TEST_DATA"]) / "longread"
    if os.environ.get("ALNKIT_TEST_DATA")
    else None
)
pytestmark = pytest.mark.skipif(
    PRIVATE is None or not PRIVATE.exists(),
    reason="需要 ALNKIT_TEST_DATA（含真实数据与旧脚本黄金，不入公开仓库）",
)


def _inputs() -> dict[str, Path]:
    assert PRIVATE is not None
    inputs = PRIVATE / "input"
    return {
        "regions": PRIVATE / "golden/breakpoints/breakpoint_classes.tsv",
        "sites": PRIVATE / "golden/breakpoints/breakpoint_sites.tsv",
        "reference": inputs / "reference.fa",
        "exog_h": inputs / "exog/insertA.fa",
        "exog_l": inputs / "exog/insertB.fa",
        "assembly": inputs / "assembly",
    }


def _require_inputs(paths: dict[str, Path]) -> None:
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        pytest.fail(f"缺输入 {missing}：跑 tests/fixtures/build_private_longread.py")


def _run(outdir: Path) -> int:
    paths = _inputs()
    _require_inputs(paths)
    return main(
        [
            "long",
            "locate",
            "--regions",
            str(paths["regions"]),
            "--sites",
            str(paths["sites"]),
            "--reference-fasta",
            str(paths["reference"]),
            "--exog",
            str(paths["exog_h"]),
            "--exog",
            str(paths["exog_l"]),
            "--assembly-dir",
            str(paths["assembly"]),
            "--outdir",
            str(outdir),
            "--threads",
            "8",
        ]
    )


@pytest.mark.slow
def test_outputs_match_upstream(tmp_path: Path) -> None:
    assert PRIVATE is not None
    golden = PRIVATE / "golden/locate"
    if not golden.exists():
        pytest.fail(f"缺黄金 {golden}：跑 tests/fixtures/build_private_longread.py")
    assert _run(tmp_path) == 0

    produced = output_files(tmp_path)
    expected = golden_files(golden)
    assert produced == expected, "产物集合与旧脚本不同"

    text_suffixes = {".tsv", ".md"}
    for name in expected:
        got, want = tmp_path / name, golden / name
        suffix = Path(name).suffix
        if suffix == ".pdf":
            # 见模块 docstring：PDF 自带生成时间，只比大小
            assert got.stat().st_size == want.stat().st_size, f"{name} 大小不同"
            continue
        if suffix in text_suffixes:
            got_text = got.read_text(encoding="utf-8").replace(str(tmp_path), "<OUTDIR>")
            want_text = want.read_text(encoding="utf-8")
            assert got_text == want_text, f"{name} 与旧脚本产物不再一致"
            continue
        assert got.read_bytes() == want.read_bytes(), f"{name} 与旧脚本产物不再一致"


@pytest.mark.slow
def test_position_tables_are_self_consistent(tmp_path: Path) -> None:
    """不依赖黄金的自洽检查：还原的参考区间必须落在染色体长度内。"""
    assert _run(tmp_path) == 0
    rows = tables.read_tsv(tmp_path / "class_assembly_positions.tsv")[1]
    # 行数随组装结构变化（2026-11-01 修了反链放置后是 26 行，之前 35 行），
    # 所以这里不断言绝对行数，只断言结构不变式。
    assert rows, "位置表为空"

    flank = [row for row in rows if row["query_type"] == "flank"]
    assert flank, "没有 flank 行"
    for row in flank:
        start, end = int(row["reference_start"]), int(row["reference_end"])
        assert 0 <= start < end
        # 50 bp 锚点：还原区间长度应接近 50（允许比对边界差几个 bp）
        assert 45 <= end - start <= 55, (row["query_name"], start, end)

    exog = [row for row in rows if row["query_type"] == "exog"]
    for row in exog:
        assert row["reference_start"] == "-", "外源没有参考坐标，应填 -"

    # contig_map.tsv 的列名是 contig（不是 assembly_contig）
    contig_lengths = {
        (row["class_id"], row["contig"]): int(row["contig_length_bp"])
        for row in tables.read_tsv(tmp_path / "contig_map.tsv")[1]
    }
    for row in rows:
        key = (row["class_id"], row["assembly_contig"])
        assert key in contig_lengths, key
        assert int(row["contig_start"]) >= 0
        assert int(row["contig_end"]) <= contig_lengths[key], key
