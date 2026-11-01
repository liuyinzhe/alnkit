"""表契约（core/tables.py）：列序、后缀规则、写入行为。"""

from __future__ import annotations

import pytest

from alnkit.core import tables


def test_every_registered_table_passes_suffix_rule() -> None:
    for table in tables.COLUMNS:
        tables.validate_columns(table)


def test_bare_coordinate_column_is_rejected() -> None:
    with pytest.raises(AssertionError):
        tables.validate_column_names("demo", ("ref_start",))


def test_length_like_columns_are_not_coordinates() -> None:
    for name in ("ref_aln_len", "chrom_length", "block_bp", "n_reads", "end_side"):
        assert not tables.is_coordinate_column(name), name
    for name in ("ref_start_0based", "insert_pos_1based", "boundary_median_0based"):
        assert tables.is_coordinate_column(name), name


def test_columns_of_unknown_table_raises() -> None:
    with pytest.raises(ValueError):
        tables.columns_of("nope")


def test_write_tsv_uses_registered_order_and_formats(tmp_path) -> None:
    target = tmp_path / "sites.tsv"
    written = tables.write_tsv(
        target,
        "insert_sites",
        [
            {
                "chrom": "chr2",
                "insert_pos_1based": 200,
                "insert_pos_0based": 199,
                "insert_boundary_0based": 200,
                "n_sr": 1,
                "n_pesr": 0,
                "n_dup": 0,
                "sample": "s2",
            },
            {
                "chrom": "chr1",
                "insert_pos_1based": 100,
                "insert_pos_0based": 99,
                "insert_boundary_0based": 100,
                "n_sr": 3,
                "n_pesr": 1,
                "n_dup": None,
                "sample": "s1",
            },
        ],
        comments=["坐标见 docs/concepts/coordinates.md"],
        sort_key=lambda row: (row["chrom"], row["insert_pos_1based"]),
    )
    assert written == 2
    header, rows = tables.read_tsv(target)
    assert header == list(tables.columns_of("insert_sites"))
    assert [row["chrom"] for row in rows] == ["chr1", "chr2"]  # 排序生效
    assert rows[0]["n_dup"] == "-"  # None -> "-"
    lines = target.read_text(encoding="utf-8").splitlines()
    assert lines[0] == "# 坐标见 docs/concepts/coordinates.md"
    assert lines[1].startswith("chrom\tinsert_pos_1based\t")


def test_write_tsv_formats_float_and_bool(tmp_path) -> None:
    target = tmp_path / "blocks.tsv"
    tables.write_tsv(
        target,
        "ribbon_blocks",
        [
            {
                "query_fragment": "a",
                "query_start_0based": 0,
                "query_end_0based": 1200,
                "strand": "+",
                "reference": "r",
                "ref_start_0based": 0,
                "ref_end_0based": 1200,
                "block_bp": 1200,
                "identity": 1.0,
            }
        ],
    )
    row = tables.read_tsv(target)[1][0]
    assert row["identity"] == "1.0000"


def test_write_tsv_rejects_unknown_column(tmp_path) -> None:
    with pytest.raises(ValueError, match="未注册的列"):
        tables.write_tsv(tmp_path / "x.tsv", "insert_sites", [{"chrom": "chr1", "oops": 1}])


def test_write_tsv_is_atomic_and_leaves_no_temp(tmp_path) -> None:
    target = tmp_path / "sites.tsv"
    tables.write_tsv(target, "insert_sites", [])
    assert target.exists()
    assert [p.name for p in tmp_path.iterdir()] == ["sites.tsv"]


def test_read_tsv_skips_comments_and_blanks(tmp_path) -> None:
    path = tmp_path / "t.tsv"
    path.write_text("# 说明\n\nchrom\tinsert_pos_1based\nchr1\t100\n", encoding="utf-8")
    header, rows = tables.read_tsv(path)
    assert header == ["chrom", "insert_pos_1based"]
    assert rows == [{"chrom": "chr1", "insert_pos_1based": "100"}]


def test_write_legacy_tsv_keeps_old_semantics(tmp_path) -> None:
    target = tmp_path / "legacy.tsv"
    written = tables.write_legacy_tsv(
        target, ["a", "b"], [["x", 1], ["y", 2.5]], comments=["说明：旧格式"]
    )
    assert written == 2
    assert target.read_text(encoding="utf-8") == "# 说明：旧格式\na\tb\nx\t1\ny\t2.5\n"


def test_read_legacy_tsv_returns_rows_as_lists(tmp_path) -> None:
    path = tmp_path / "t.tsv"
    path.write_text("# 注释\n\na\tb\n1\t2\n", encoding="utf-8")
    header, rows = tables.read_legacy_tsv(path)
    assert header == ["a", "b"]
    assert rows == [["1", "2"]]


def test_coordinate_word_plus_length_suffix_is_still_a_coordinate() -> None:
    """`ref_end_bp` 这种混搭必须被拒（否则契约校验能被长度后缀整个绕过）。"""
    with pytest.raises(AssertionError):
        tables.validate_column_names("demo", ("ref_end_bp",))
    # 真正的长度列不受影响
    for name in ("chrom_length", "block_bp", "n_reads", "end_side", "junction_type"):
        assert not tables.is_coordinate_column(name), name


def test_write_tsv_keeps_existing_permissions(tmp_path) -> None:
    target = tmp_path / "sites.tsv"
    tables.write_tsv(target, "insert_sites", [])
    assert target.stat().st_mode & 0o777 != 0o600, "mkstemp 的 0600 不该泄漏成产物权限"
    target.chmod(0o640)
    tables.write_tsv(target, "insert_sites", [])
    assert target.stat().st_mode & 0o777 == 0o640


def test_write_tsv_failure_preserves_existing_file(tmp_path) -> None:
    """失败路径：目标文件内容不变，且不留临时文件。"""
    target = tmp_path / "sites.tsv"
    tables.write_tsv(target, "insert_sites", [{"chrom": "chr1", "n_sr": 1}])
    before = target.read_bytes()
    with pytest.raises(ValueError, match="未注册的列"):
        tables.write_tsv(target, "insert_sites", [{"chrom": "chr1", "oops": 1}])
    assert target.read_bytes() == before
    assert [p.name for p in tmp_path.iterdir()] == ["sites.tsv"]


def test_write_legacy_tsv_is_atomic_too(tmp_path) -> None:
    class Boom:
        def __str__(self) -> str:
            raise RuntimeError("boom")

    target = tmp_path / "legacy.tsv"
    tables.write_legacy_tsv(target, ["a"], [["1"]])
    before = target.read_bytes()
    with pytest.raises(RuntimeError):
        tables.write_legacy_tsv(target, ["a"], [["write ok"], [Boom()]])  # 写到一半抛错
    assert target.read_bytes() == before, "半截表不该覆盖旧文件"
    assert [p.name for p in tmp_path.iterdir()] == ["legacy.tsv"]
