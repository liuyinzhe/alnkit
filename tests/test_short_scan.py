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


def test_removed_legacy_profile_and_scope_flags_are_rejected() -> None:
    """0.7.0 删掉的 --legacy-profile / --interference-scope 必须由 argparse 拒绝。"""
    for flags in (
        ["--legacy-profile"],
        ["--interference-scope", "other"],
        ["--insert-origin-policy", "legacy"],
    ):
        with pytest.raises(SystemExit) as error:
            main(["short", "scan", "--bam-dir", ".", *flags])
        assert error.value.code == 2


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
    assert scan.determine_interference_sources(
        "owned_query", {"target_ref"}, {"homolog_ref"}, "target_ref"
    ) == {"homolog_ref"}
    assert (
        scan.determine_interference_sources("owned_query", {"target_ref"}, set(), "target_ref")
        == set()
    )


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

    # 干扰表的语义在 0.7.0 变了（query 不再算源头），不能再用上游产物逐字节对账；
    # 上游的两张表（Insert_Position / reads_info）仍然逐字节比。
    not_comparable = {"sampleA.Insert_Interference.result.xls"}
    produced = set(output_files(tmp_path)) - not_comparable
    expected = set(golden_files(golden)) - not_comparable
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

    # 干扰表不再与上游逐字节比（语义已改），但不能完全没覆盖：结构必须仍然成立
    interference_path = tmp_path / "sampleA.Insert_Interference.result.xls"
    assert interference_path.exists(), "干扰表没有产出"
    interference_lines = [
        line for line in interference_path.read_text(encoding="utf-8").splitlines() if line
    ]
    assert len(interference_lines[0].split("\t")) == 15, "干扰表列数变了"
    position_lines = [
        line
        for line in (tmp_path / "sampleA.Insert_Position.result.xls")
        .read_text(encoding="utf-8")
        .splitlines()
        if line
    ]
    assert len(interference_lines) == len(position_lines), "干扰表与主表行数不一致"


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


### markdup 核对（坐标 + MD） ###


def _write_bam(path: Path, reads, contig: str = "ctg1", length: int = 10_000) -> Path:
    """写一个最小 BAM：reads = [(name, start, is_read1, is_duplicate, MD), ...]。"""
    import pysam

    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": contig, "LN": length}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, start, is_read1, is_duplicate, md in sorted(reads, key=lambda item: item[1]):
            record = pysam.AlignedSegment()
            record.query_name = name
            record.query_sequence = "A" * 100
            record.flag = (64 if is_read1 else 128) | (1024 if is_duplicate else 0)
            record.reference_id = 0
            record.reference_start = start
            record.mapping_quality = 60
            record.cigartuples = [(0, 100)]
            record.query_qualities = pysam.qualitystring_to_array("I" * 100)
            record.set_tag("MD", md)
            out.write(record)
    return path


def test_duplicate_signature_uses_both_ends() -> None:
    """双端签名必须把两端 5' 坐标都放进去，否则不同分子会被误判成一组。"""
    from alnkit.shortread.bam import duplicate_signature

    bam = _write_bam(
        Path(os.environ.get("PYTEST_TMPDIR", "/tmp")) / "sig.bam",
        [("r1", 100, True, False, "100"), ("r1", 900, False, False, "100")],
    )
    import pysam

    with pysam.AlignmentFile(str(bam)) as handle:
        records = list(handle.fetch(until_eof=True))
    first = next(record for record in records if record.is_read1)
    second = next(record for record in records if record.is_read2)
    signature = duplicate_signature(first, second)
    assert signature[0] == "pair"
    assert signature[1:] == ("ctg1", "+", 100, "ctg1", "+", 900)


def test_check_markdup_consistency_healthy(tmp_path: Path) -> None:
    """正常 markdup：每组恰好 1 条未标记；MD 一致时计入 seq_identical。"""
    from alnkit.shortread.bam import check_markdup_consistency

    bam = _write_bam(
        tmp_path / "ok.bam",
        [
            ("a", 100, True, False, "100"),
            ("a", 900, False, False, "100"),
            ("b", 100, True, True, "100"),
            ("b", 900, False, True, "100"),
            ("c", 100, True, True, "100"),
            ("c", 900, False, True, "100"),
        ],
    )
    statistics = check_markdup_consistency(bam)
    assert statistics["groups_total"] == 1
    assert statistics["groups_ok"] == 1
    assert statistics["groups_multi_unflagged"] == 0
    assert statistics["duplicates_checked"] == 2
    assert statistics["duplicates_seq_identical"] == 2


def test_check_markdup_consistency_detects_missing_markdup(tmp_path: Path) -> None:
    """没有 FLAG 0x400 时，整组都未标记 → 应报 groups_multi_unflagged。"""
    from alnkit.shortread.bam import check_markdup_consistency

    bam = _write_bam(
        tmp_path / "nodupe.bam",
        [
            ("a", 100, True, False, "100"),
            ("a", 900, False, False, "100"),
            ("b", 100, True, False, "100"),
            ("b", 900, False, False, "100"),
        ],
    )
    statistics = check_markdup_consistency(bam)
    assert statistics["groups_total"] == 1
    assert statistics["groups_ok"] == 0
    assert statistics["groups_multi_unflagged"] == 1
    assert statistics["duplicates_checked"] == 0


def test_check_markdup_reports_sequence_difference(tmp_path: Path) -> None:
    """MD 模式不同 → 位置巧合而非同一分子，计入 seq_different 但不算错误。"""
    from alnkit.shortread.bam import check_markdup_consistency

    bam = _write_bam(
        tmp_path / "md.bam",
        [
            ("a", 100, True, False, "100"),
            ("a", 900, False, False, "100"),
            ("b", 100, True, True, "40A59"),
            ("b", 900, False, True, "40A59"),
        ],
    )
    statistics = check_markdup_consistency(bam)
    assert statistics["groups_ok"] == 1
    assert statistics["duplicates_checked"] == 1
    assert statistics["duplicates_seq_different"] == 1


def test_scan_cli_accepts_check_markdup(tmp_path: Path) -> None:
    """--check-markdup 是纯日志功能，不改变任何产物集合。"""
    bam = _write_bam(
        tmp_path / "sampleA.bam",
        [
            ("a", 100, True, False, "60S40M"),
            ("a", 900, False, False, "60S40M"),
            ("b", 100, True, True, "60S40M"),
            ("b", 900, False, True, "60S40M"),
        ],
    )
    out_plain, out_check = tmp_path / "plain", tmp_path / "check"
    for out, extra in ((out_plain, []), (out_check, ["--check-markdup"])):
        code = main(
            [
                "short",
                "scan",
                "--bam",
                str(bam),
                "--query-chr",
                "ctg1",
                "--insert-chrs",
                "ctg1",
                "--out-dir",
                str(out),
                *extra,
            ]
        )
        assert code == 0
    assert output_files(out_plain) == output_files(out_check)


### 位点级最小支持数 / 干扰口径 ###


def test_interference_sources_always_exclude_the_query() -> None:
    """query 自己永远不算干扰源头（断点 read 的另一段落在 query 上是正常形态）。"""
    from alnkit.shortread.interference import determine_interference_sources

    assert determine_interference_sources("owned_query", {"q"}, {"q"}, "q") == set()
    assert determine_interference_sources("owned_query", {"q"}, {"q", "h"}, "q") == {"h"}
    assert determine_interference_sources("ambiguous", {"q", "h"}, set(), "q") == {"h"}
    # ambiguous 且只命中 query 本身 → 没有其它插入，不构成干扰
    assert determine_interference_sources("ambiguous", {"q"}, set(), "q") == set()


def test_filter_by_minimum_support() -> None:
    """位点级过滤：按 SR+PESR 总支持数丢位点，结构三份同步；默认 1 原样返回。"""
    from alnkit.shortread.merge import filter_by_minimum_support
    from alnkit.shortread.model import BreakpointEvidence

    def bucket(sr: int, pesr: int) -> BreakpointEvidence:
        return BreakpointEvidence(
            split_read_names={f"s{i}" for i in range(sr)},
            pair_split_read_names={f"p{i}" for i in range(pesr)},
            owned_query_read_names=set(),
            owned_other_read_names=set(),
            ambiguous_read_names=set(),
            no_insert_hit_read_names=set(),
            consensus_read_names=set(),
            interference_source_read_names={},
        )

    evidence = {("chr1", 100): bucket(2, 0), ("chr1", 900): bucket(1, 0)}
    duplicates = {("chr1", 100): {"d1"}, ("chr1", 900): {"d2"}}
    merged = {("chr1", 100): (100,), ("chr1", 900): (900,)}

    same = filter_by_minimum_support(evidence, duplicates, merged, 1)
    assert same[0] is evidence and same[3] == 0  # 默认不过滤，原对象返回

    kept_evidence, kept_duplicates, kept_merged, dropped = filter_by_minimum_support(
        evidence, duplicates, merged, 2
    )
    assert set(kept_evidence) == {("chr1", 100)}
    assert set(kept_duplicates) == {("chr1", 100)}
    assert set(kept_merged) == {("chr1", 100)}
    assert dropped == 1


def test_min_sr_support_must_be_positive() -> None:
    from alnkit.shortread.model import ScanSettings

    with pytest.raises(UsageError, match="--min-sr-support"):
        scan.resolve_runtime_config(
            ScanSettings(
                bam=(Path("x.bam"),), query_chr="ctg1", insert_chrs="ctg1", minimum_site_support=0
            )
        )


def test_removed_origin_policies_are_rejected() -> None:
    """0.7.0 删掉的 owned-strict / lenient 必须由 argparse 拒绝（不许静默改行为）。"""
    for policy in ("owned-strict", "lenient"):
        with pytest.raises(SystemExit) as error:
            main(["short", "scan", "--bam-dir", ".", "--insert-origin-policy", policy])
        assert error.value.code == 2


def test_drop_shared_reads_flag_is_accepted_with_owned_and_consensus(tmp_path) -> None:
    """新开关与 owned/consensus 都能组合（旧的枚举值做不到这一点）。"""
    for policy in ("owned", "consensus"):
        code = main(
            [
                "short",
                "scan",
                "--bam",
                str(tmp_path / "missing.bam"),
                "--query-chr",
                "insertA",
                "--insert-origin-policy",
                policy,
                "--drop-shared-reads",
                "--out-dir",
                str(tmp_path),
            ]
        )
        # 参数被接受（否则 argparse 会以 2 退出并说 invalid choice），失败发生在读 BAM（1）
        assert code == 1
