"""``alnkit plot ribbon``（C1）：纯逻辑定点断言 + 与旧脚本的黄金回归。

``fixtures/ribbon/expected.blocks.tsv`` 由旧脚本 ``synteny_ribbon.py`` v1.10.1
在 ``example4_alignment.paf`` 上生成；同一输入下新实现必须逐字节一致。
"""

from __future__ import annotations

from pathlib import Path

import pytest

from alnkit.cli import main
from alnkit.core import tables
from alnkit.core.errors import UsageError
from alnkit.plot import ribbon

FIXTURES = Path(__file__).parent / "fixtures" / "ribbon"
PAF = FIXTURES / "example4_alignment.paf"
EXPECTED_BLOCKS = FIXTURES / "expected.blocks.tsv"


def _render(tmp_path: Path, **overrides) -> dict[str, object]:
    settings = ribbon.RibbonSettings(
        outdir=str(tmp_path), name="ex4", paf=str(PAF), reference_name="assembly_a", **overrides
    )
    return ribbon.run(settings)


# ---------------------------------------------------------------- 纯逻辑
def test_parse_spec_with_and_without_contig() -> None:
    assert ribbon.parse_spec("ins=/x/y.fa:contig1") == ("ins", "/x/y.fa", "contig1")
    assert ribbon.parse_spec("ins=/x/y.fa") == ("ins", "/x/y.fa", None)


@pytest.mark.parametrize("spec", ["no-equals-sign", "=path.fa", "name="])
def test_parse_spec_rejects_bad_input(spec: str) -> None:
    with pytest.raises(UsageError):
        ribbon.parse_spec(spec)


def test_parse_paf_text_columns() -> None:
    first_line = PAF.read_text(encoding="utf-8").splitlines()[0]
    block = ribbon.parse_paf_text(first_line)[0]
    assert block["record"] == "insert_alpha"
    assert (block["q_start"], block["q_end"]) == (0, 1200)
    assert (block["r_start"], block["r_end"]) == (25500, 26700)
    assert block["block_bp"] == 1200
    assert block["identity"] == 1.0


def test_resolve_row_order_variants() -> None:
    default = ribbon.resolve_row_order(["b", "a"], None, "top")
    assert default == (["reference", "b", "a"], [])
    assert ribbon.resolve_row_order(["b", "a"], None, "bottom")[0] == ["b", "a", "reference"]
    # 逗号分隔，reference 可插在中间
    assert ribbon.resolve_row_order(["b", "a"], "a,reference,b")[0] == ["a", "reference", "b"]
    # 名字写错要有告警而不是静默丢行
    order, warnings = ribbon.resolve_row_order(["b"], "b,typo")
    assert order == ["reference", "b"]
    assert warnings and "typo" in warnings[0]


def test_order_and_order_file_are_mutually_exclusive(tmp_path: Path) -> None:
    with pytest.raises(UsageError, match="只能用一个"):
        ribbon.resolve_row_order(["b"], "b", "top", ("b",))


def test_preset_tiers() -> None:
    assert ribbon.tier_of_length(1000) == 0  # >=1000 -> asm20
    assert ribbon.tier_of_length(999) == 1  # 40-999 -> map-ont 小种子
    assert ribbon.tier_of_length(40) == 1
    assert ribbon.tier_of_length(39) is None


def test_query_preset_parsing() -> None:
    label, extra = ribbon.parse_query_preset("x=map-ont -k11 -w5")
    assert label == "x"
    assert extra == ["-x", "map-ont", "-k11", "-w5"]
    assert ribbon.parse_query_preset("x=sr")[1] == ["-x", "sr"]


def test_unknown_query_preset_label_is_usage_error() -> None:
    with pytest.raises(UsageError, match="标签"):
        ribbon.collect_preset_overrides(["nope=sr"], {"a"})


def test_bezier_endpoints() -> None:
    points = ribbon.bezier_points((0, 0), (0, 1), (10, 1), (10, 0))
    assert points[0] == (0.0, 0.0)
    assert points[-1] == (10.0, 0.0)


# ------------------------------------------------------------ 黄金回归
def test_blocks_tsv_is_byte_identical(tmp_path: Path) -> None:
    result = _render(tmp_path)
    assert result["reference"] == "assembly_a"
    assert result["queries"] == 1
    assert result["blocks"] == 3
    produced = Path(str(result["blocks_tsv"])).read_bytes()
    assert produced == EXPECTED_BLOCKS.read_bytes(), "blocks.tsv 与旧脚本产物不再一致"


def test_png_and_pdf_written(tmp_path: Path) -> None:
    result = _render(tmp_path)
    png = Path(str(result["png"]))
    pdf = png.with_suffix(".pdf")
    assert png.stat().st_size > 5_000
    assert pdf.stat().st_size > 1_000


def test_labels_and_filter_thresholds(tmp_path: Path) -> None:
    # 阈值调高到 2000 bp：1200 bp 的块全被滤掉，只写表头（旧行为：告警但不报错）
    result = _render(tmp_path, min_block_bp=2_000)
    assert result["blocks"] == 0
    header = Path(str(result["blocks_tsv"])).read_text(encoding="utf-8").splitlines()
    assert len(header) == 1
    assert header[0].startswith("query_fragment\tquery_start\tquery_end")


def test_placement_and_region_tables(tmp_path: Path) -> None:
    placement = tmp_path / "placements.tsv"
    regions = tmp_path / "regions.tsv"
    result = _render(tmp_path, placement_table=str(placement), region_table=str(regions))
    assert result["blocks"] == 3

    # placements.tsv：一行一处落点（该片段在参考上有 3 处）
    header, rows = tables.read_tsv(placement)
    assert header == [name for name, _ in ribbon.PLACEMENT_COLUMNS]
    assert len(rows) == 3
    assert {(row["ref_start"], row["ref_end"], row["strand"]) for row in rows} == {
        ("5000", "6200", "+"),
        ("6400", "7600", "-"),
        ("7800", "9000", "+"),
    }
    for row in rows:
        assert row["query_name"] == "insert_alpha"
        assert row["n_placements"] == "3"
        assert row["identity_percent"] == "100.00%"
        assert row["covers_full_query"] == "yes"
    assert [row["placement_index"] for row in rows] == ["1", "2", "3"]

    # regions.tsv：以参考为坐标系，三个互不重叠的区域，按坐标升序编号
    header, rows = tables.read_tsv(regions)
    assert header == [name for name, _ in ribbon.REGION_TABLE_COLUMNS]
    assert [row["block_id"] for row in rows] == ["1", "2", "3"]
    assert [(row["ref_start"], row["ref_end"]) for row in rows] == [
        ("5000", "6200"),
        ("6400", "7600"),
        ("7800", "9000"),
    ]
    assert [row["strand"] for row in rows] == ["+", "-", "+"]
    # 第二行的“与上一区域间隔”= 6400-6200
    assert rows[1]["gap_from_previous_ref_region_bp"] == "200"


def test_cli_ribbon_end_to_end(tmp_path: Path, capsys) -> None:
    code = main(
        [
            "plot",
            "ribbon",
            "--paf",
            str(PAF),
            "--reference-name",
            "assembly_a",
            "--outdir",
            str(tmp_path),
            "--name",
            "cli",
        ]
    )
    assert code == 0
    assert (tmp_path / "cli.blocks.tsv").exists()
    assert "映射块表" in capsys.readouterr().out


def test_cli_ribbon_requires_name(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(["plot", "ribbon", "--paf", str(PAF), "--outdir", str(tmp_path)])
    assert excinfo.value.code == 2


def test_fasta_mode_requires_working_minimap2(tmp_path: Path, monkeypatch) -> None:
    """方式 1 需要 minimap2：缺失时应给出可读错误而不是栈。"""
    (tmp_path / "ref.fa").write_text(">ref\n" + "ACGT" * 500 + "\n", encoding="utf-8")
    (tmp_path / "query.fa").write_text(">q1\n" + "ACGT" * 250 + "\n", encoding="utf-8")
    monkeypatch.setattr(ribbon, "align_file_to_reference", _raise_missing_minimap2)
    settings = ribbon.RibbonSettings(
        outdir=str(tmp_path),
        name="x",
        reference=f"ref={tmp_path / 'ref.fa'}",
        query=(f"q={tmp_path / 'query.fa'}",),
    )
    with pytest.raises(UsageError, match="minimap2"):
        ribbon.run(settings)


def _raise_missing_minimap2(*args, **kwargs):
    raise UsageError("找不到 minimap2")
