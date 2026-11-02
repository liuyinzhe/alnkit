"""``alnkit long assemble``（A6）：整目录逐字节对账 + 确定性。

背景：旧实现用 ``for unplaced_name in list(unplaced)`` 遍历 ``set``，结果随
``PYTHONHASHSEED`` 漂（同输入连跑三次出三种 layout）。该 bug 已在上游
（commit ``20b4dda``）与 alnkit 同步修掉：遍历改为 ``sorted(...)``。

验收：

1. ``prop 产物``（contigs.fasta / assembly_stats.tsv / 汇总表）与上游脚本**逐字节一致**；
   ``layout.tsv`` 不再输出（无下游消费者，且与 contig 的一致性只有约 28%，见模块文档）；
2. 不同 ``PYTHONHASHSEED`` 下两次运行产物一致（确定性回归）；
3. 上游写 ``layout.tsv``、alnkit 不写，所以黄金里也删掉该文件（见 builder）。

需要真实数据（不入公开仓库）：

    python tests/fixtures/build_private_longread.py
    export ALNKIT_TEST_DATA=~/.local/share/alnkit-testdata
    python -m pytest tests/test_long_assemble.py -m slow
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest
from _fs import golden_files

from alnkit.cli import main
from alnkit.longread.assemble import placement_from_block

PRIVATE = (
    Path(os.environ["ALNKIT_TEST_DATA"]) / "longread"
    if os.environ.get("ALNKIT_TEST_DATA")
    else None
)
requires_data = pytest.mark.skipif(
    PRIVATE is None or not PRIVATE.exists(),
    reason="需要 ALNKIT_TEST_DATA（含真实数据与旧脚本黄金，不入公开仓库）",
)


def _inputs() -> tuple[Path, Path]:
    assert PRIVATE is not None
    return PRIVATE / "input/classes", PRIVATE / "input/overlaps.paf"


def _run(outdir: Path, env: dict[str, str] | None = None) -> None:
    classes, paf = _inputs()
    if env is None:
        code = main(
            [
                "long",
                "assemble",
                "--classes-dir",
                str(classes),
                "--paf",
                str(paf),
                "--outdir",
                str(outdir),
            ]
        )
        assert code == 0
        return
    subprocess.run(
        [
            sys.executable,
            "-m",
            "alnkit.cli",
            "long",
            "assemble",
            "--classes-dir",
            str(classes),
            "--paf",
            str(paf),
            "--outdir",
            str(outdir),
        ],
        check=True,
        env={**os.environ, **env},
        capture_output=True,
    )


@requires_data
@pytest.mark.slow
def test_outputs_match_upstream_script(tmp_path: Path) -> None:
    """整目录必须与上游脚本（已修确定性）逐字节一致。

    A 方案要在对账前关掉：`--no-support-report` 不额外写 `assembly_support.tsv`。
    支撑度报告默认开启，属于黄金基线之外的**增量产物**；其余产物不受影响。
    """
    assert PRIVATE is not None
    classes, paf = _inputs()
    if not classes.exists() or not paf.exists():
        pytest.fail("缺少输入（运行 tests/fixtures/build_private_longread.py）")
    golden = PRIVATE / "golden/assemble"

    assert (
        main(
            [
                "long",
                "assemble",
                "--classes-dir",
                str(classes),
                "--paf",
                str(paf),
                "--outdir",
                str(tmp_path),
                "--no-support-report",
            ]
        )
        == 0
    )

    produced = golden_files(tmp_path)
    expected = golden_files(golden)
    assert produced == expected, "产物集合与上游脚本不同"
    for name in expected:
        assert (tmp_path / name).read_bytes() == (golden / name).read_bytes(), (
            f"{name} 与上游脚本产物不再一致"
        )


@requires_data
@pytest.mark.slow
def test_assembly_is_deterministic_across_hash_seeds(tmp_path: Path) -> None:
    """同一输入、不同 PYTHONHASHSEED 下，产物必须完全一致。"""
    assert PRIVATE is not None
    classes, paf = _inputs()
    if not classes.exists() or not paf.exists():
        pytest.fail("缺少输入（运行 tests/fixtures/build_private_longread.py）")
    first, second = tmp_path / "seed1", tmp_path / "seed2"
    _run(first, {"PYTHONHASHSEED": "1"})
    _run(second, {"PYTHONHASHSEED": "2"})

    files_first = golden_files(first)
    files_second = golden_files(second)
    assert files_first == files_second
    for name in files_first:
        assert (first / name).read_bytes() == (second / name).read_bytes(), name


def test_placement_from_block_reverse_strand_partner() -> None:
    """反链 overlap：unplaced 的起点要与 placed 的**终点**配对。

    回归 2026-11-01 复审发现的 bug：旧实现起点对起点，导致反链邻居被放偏
    (aligned_bp - 1) bp，contig 连接处出现重复序列。
    """
    block = {
        "query": "u",
        "target": "p",
        "strand": "-",
        "query_start": 80,
        "query_end": 100,
        "target_start": 0,
        "target_end": 20,
    }
    # p 占 contig [0,100)；u[99]↔p[0]（contig 0），u[0] 落在 contig 99 → start=0, direction=-1
    assert placement_from_block("p", 0, 1, 100, "u", 100, block) == (0, -1)

    # 镜像情形：block 里 query/target 互换，结论相同
    mirror = {
        "query": "p",
        "target": "u",
        "strand": "-",
        "query_start": 0,
        "query_end": 20,
        "target_start": 80,
        "target_end": 100,
    }
    assert placement_from_block("p", 0, 1, 100, "u", 100, mirror) == (0, -1)

    # 正链不受影响（原 docstring 的例子）
    plus = dict(block, strand="+")
    assert placement_from_block("p", 0, 1, 100, "u", 100, plus) == (-80, 1)
