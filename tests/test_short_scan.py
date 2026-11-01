"""``alnkit short scan``（B1/B2）：配置校验、纯逻辑定点断言、私有数据黄金回归。

黄金回归需要真实测序数据（不入公开仓库）：

    export ALNKIT_TEST_DATA=~/.local/share/alnkit-testdata
    python tests/fixtures/build_private_shortread.py     # 首次准备（见该脚本）
    python -m pytest tests/test_short_scan.py -m slow

没有数据时这些用例会 skip（带原因），CI 默认走这条路径。
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest
from _fs import golden_files, output_files

from alnkit.cli import main
from alnkit.core import tables
from alnkit.core.errors import UsageError
from alnkit.shortread import scan

FIXTURES = Path(__file__).parent / "fixtures"


def _manifest() -> dict:
    """私有数据里的真实序列名只存在本机 manifest 中（仓库里不出现）。"""
    assert PRIVATE is not None
    path = PRIVATE / "manifest.json"
    if not path.exists():
        pytest.fail(f"缺 {path}：跑 tests/fixtures/build_private_shortread.py")
    return json.loads(path.read_text(encoding="utf-8"))


PRIVATE = (
    Path(os.environ["ALNKIT_TEST_DATA"]) / "shortread"
    if os.environ.get("ALNKIT_TEST_DATA")
    else None
)
requires_data = pytest.mark.skipif(
    PRIVATE is None or not PRIVATE.exists(),
    reason="需要 ALNKIT_TEST_DATA（含真实 BAM 与黄金表，不入公开仓库）",
)


def _settings(**overrides) -> scan.ScanSettings:
    base = {"query_chr": "target_ref", "bam": (FIXTURES / "segtable" / "reads.bam",)}
    return scan.ScanSettings(**{**base, **overrides})


# ------------------------------------------------------------ 配置校验
def test_query_chr_is_required() -> None:
    with pytest.raises(UsageError, match="--query-chr"):
        scan.resolve_runtime_config(scan.ScanSettings(bam=(Path("x.bam"),)))


def test_insert_chrs_defaults_to_query_and_keeps_order() -> None:
    config = scan.resolve_runtime_config(_settings())
    assert config.insert_chromosomes == ("target_ref",)

    config = scan.resolve_runtime_config(_settings(insert_chrs="decoy_ref, homolog_ref,target_ref"))
    assert config.insert_chromosomes == ("decoy_ref", "homolog_ref", "target_ref")

    # query 不在列表里时自动补到最前面，而不是报错
    config = scan.resolve_runtime_config(_settings(insert_chrs="decoy_ref,homolog_ref"))
    assert config.insert_chromosomes[0] == "target_ref"


def test_legacy_profile_expands_to_legacy_switches() -> None:
    config = scan.resolve_runtime_config(_settings(legacy_profile=True))
    assert config.insert_origin_policy == "legacy"
    assert config.mismatch_source == "md"
    assert config.legacy_coordinates is True
    # legacy 不做归属判定，候选闸门必须自动关闭
    assert config.use_candidate_gate is False


def test_candidate_gate_defaults_on_and_can_be_disabled() -> None:
    assert scan.resolve_runtime_config(_settings()).use_candidate_gate is True
    assert (
        scan.resolve_runtime_config(_settings(no_candidate_gate=True)).use_candidate_gate is False
    )


def test_duplicate_filter_defaults() -> None:
    config = scan.resolve_runtime_config(_settings())
    assert config.exclude_duplicates is True
    assert config.require_duplicate_flags is False
    assert (
        scan.resolve_runtime_config(_settings(no_duplicate_filter=True)).exclude_duplicates is False
    )


def test_empty_bam_dir_is_usage_error(tmp_path: Path) -> None:
    with pytest.raises(UsageError, match="未找到 BAM"):
        scan.resolve_runtime_config(_settings(bam=(), bam_dir=tmp_path))


def test_multiple_bams_with_sample_name_is_usage_error(tmp_path: Path) -> None:
    first = tmp_path / "a.sorted.bam"
    second = tmp_path / "b.sorted.bam"
    for path in (first, second):
        path.write_bytes(b"")
    with pytest.raises(UsageError, match="--sample-name"):
        scan.resolve_runtime_config(_settings(bam=(first, second), sample_name="S1"))


def test_bam_threads_must_be_positive() -> None:
    with pytest.raises(UsageError, match="--bam-threads"):
        scan.resolve_runtime_config(_settings(bam_threads=0))


# ------------------------------------------------------------ 纯逻辑
def test_cluster_positions_within_window() -> None:
    assert scan.cluster_positions_within_window([5000000, 5000010, 5100000], 10) == [
        [5000000, 5000010],
        [5100000],
    ]
    assert scan.cluster_positions_within_window([100, 200], 0) == [[100], [200]]


def test_parse_xa_entry_keeps_only_insert_chromosomes() -> None:
    hit = scan.parse_xa_entry("homolog_ref,+6009,69M81S,3,0;", ["homolog_ref", "target_ref"])
    assert hit is not None
    assert (hit.chromosome, hit.source, hit.mismatch_count, hit.clip_length) == (
        "homolog_ref",
        "xa",
        3,
        81,
    )
    assert scan.parse_xa_entry("chr1,-100,50M,2,0;", ["homolog_ref"]) is None


def test_determine_interference_sources() -> None:
    assert scan.determine_interference_sources("owned_query", {"target_ref"}, {"homolog_ref"}) == {
        "homolog_ref"
    }
    assert scan.determine_interference_sources("owned_query", {"target_ref"}, set()) == set()


def test_builtin_self_test_passes(capsys) -> None:
    """B2：旧脚本的 30+ 条口径断言，抽出来直接可跑。"""
    scan.run_self_test()
    assert "self-test 通过" in capsys.readouterr().out


def test_cli_self_test_flag() -> None:
    assert main(["short", "scan", "--self-test"]) == 0


def test_cli_missing_query_chr_exits_2(tmp_path: Path) -> None:
    assert (
        main(["short", "scan", "--bam", str(tmp_path / "nope.bam"), "--out-dir", str(tmp_path)])
        == 2
    )


# -------------------------------------------------------- 私有数据黄金回归
@requires_data
@pytest.mark.slow
@pytest.mark.parametrize(
    ("case", "bam_name", "extra"),
    [
        (
            "owned",
            "sampleA.sorted.bam",
            ["--merge-window-bp", "10"],
        ),
        (
            "markdup",
            "sampleA.markdup.bam",
            [],
        ),
        (
            "consensus",
            "sampleA.sorted.bam",
            ["--insert-origin-policy", "consensus", "--merge-window-bp", "10"],
        ),
    ],
)
def test_golden_tables_byte_identical(
    tmp_path: Path, case: str, bam_name: str, extra: list[str]
) -> None:
    assert PRIVATE is not None
    bam = PRIVATE / "bam" / bam_name
    if not bam.exists():
        pytest.fail(f"缺少 {bam}（运行 tests/fixtures/build_private_shortread.py）")
    golden = PRIVATE / "golden" / case

    code = main(
        [
            "short",
            "scan",
            "--bam",
            str(bam),
            "--query-chr",
            _manifest()["query_chr"],
            "--insert-chrs",
            ",".join(_manifest()["insert_chrs"]),
            "--out-dir",
            str(tmp_path),
            "--log-level",
            "ERROR",
            *extra,
        ]
    )
    assert code == 0

    produced = set(output_files(tmp_path))
    expected = set(golden_files(golden))
    # 黄金来自上游脚本（只有三张 legacy 表）；alnkit 额外写 sites.canonical.tsv，
    # 所以断言"上游产物全部存在且逐字节一致 + 多出来的只有规范表"。
    assert expected <= produced, f"{case}: 缺产物 {sorted(expected - produced)}"
    assert produced - expected == {"sampleA.sites.canonical.tsv"}, (
        f"{case}: 多出未预期的产物 {sorted(produced - expected)}"
    )

    for name in sorted(expected):
        assert (tmp_path / name).read_bytes() == (golden / name).read_bytes(), (
            f"{case}/{name} 与旧脚本产物不再一致"
        )


@requires_data
@pytest.mark.slow
def test_canonical_sites_table_matches_legacy(tmp_path: Path) -> None:
    """规范表必须与 legacy 主表逐行对齐，并满足 1-based/0-based 恒等式。"""
    assert PRIVATE is not None
    bam = PRIVATE / "bam" / "sampleA.markdup.bam"
    if not bam.exists():
        pytest.fail(f"缺 {bam}（跑 tests/fixtures/build_private_shortread.py）")

    assert (
        main(
            [
                "short",
                "scan",
                "--bam",
                str(bam),
                "--query-chr",
                _manifest()["query_chr"],
                "--insert-chrs",
                ",".join(_manifest()["insert_chrs"]),
                "--out-dir",
                str(tmp_path),
                "--log-level",
                "ERROR",
            ]
        )
        == 0
    )

    sample = "sampleA"
    _, legacy = tables.read_tsv(tmp_path / f"{sample}.Insert_Position.result.xls")
    _, canonical = tables.read_tsv(tmp_path / f"{sample}.sites.canonical.tsv")
    assert len(canonical) == len(legacy) > 0

    for old, new in zip(legacy, canonical, strict=True):
        assert new["chrom"] == old["Chromosome"]
        assert new["insert_boundary_0based"] == old["Insert_Position"]
        assert new["n_sr"] == old["split_read_support"]
        assert new["n_pesr"] == old["pair_split_read_support"]
        assert new["n_dup"] == old["dup_read_count"]
        # 恒等式：0-based 碱基 + 1 == 1-based 碱基 == 0-based 边界
        assert (
            int(new["insert_pos_0based"]) + 1
            == int(new["insert_pos_1based"])
            == int(new["insert_boundary_0based"])
        )
        assert new["sample"] == sample
        # 插入归属：query 自己不算竞争者；其余按支持 read 数排序
        query_chr = _manifest()["query_chr"]
        assert new["query_insert"] == query_chr
        assert new["nearest_other_insert"] != query_chr
        assert query_chr not in new["other_inserts"].split(",")


def test_split_other_inserts_ignores_query_and_ranks_by_support() -> None:
    from alnkit.shortread.report import _split_other_inserts

    sources = {
        "target_ref": {"r1", "r2", "r3"},  # query 自己，不算竞争者
        "homolog_ref": {"r1", "r2"},  # 2 条支持 → 最近
        "decoy_ref": {"r3"},  # 1 条 → 其余
    }
    nearest, others = _split_other_inserts(sources, "target_ref")
    assert nearest == "homolog_ref"
    assert others == ["decoy_ref"]


def test_split_other_inserts_handles_ties_and_empty() -> None:
    from alnkit.shortread.report import _split_other_inserts

    # 并列时取名字小的，保证可复现
    nearest, others = _split_other_inserts({"b_ref": {"r1"}, "a_ref": {"r2"}}, "q_ref")
    assert (nearest, others) == ("a_ref", ["b_ref"])
    # 只有 query 命中 → 没有竞争者
    assert _split_other_inserts({"q_ref": {"r1"}}, "q_ref") == ("-", [])
    assert _split_other_inserts({}, "q_ref") == ("-", [])
