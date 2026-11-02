"""文档一致性：生成物必须与代码同步。"""

from __future__ import annotations

from pathlib import Path

from alnkit._docs import render_table_docs

ROOT = Path(__file__).resolve().parent.parent


def test_tables_doc_is_up_to_date() -> None:
    target = ROOT / "docs" / "concepts" / "tables.md"
    assert target.exists(), "缺少 docs/concepts/tables.md，运行 python scripts/gen_table_docs.py"
    assert target.read_text(encoding="utf-8") == render_table_docs(), (
        "docs/concepts/tables.md 与 core.tables 不一致：运行 python scripts/gen_table_docs.py"
    )


def test_docs_pages_exist() -> None:
    for name in ("index.md", "install.md", "quickstart.md", "concepts/coordinates.md"):
        assert (ROOT / "docs" / name).exists(), f"缺少文档页 docs/{name}"


def test_readme_mentions_exit_codes() -> None:
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "0` 成功 / `1` 运行期错误 / `2`" in readme


def test_cli_docs_are_up_to_date() -> None:
    """docs/cli/*.md 必须与 argparse 定义同步（生成脚本 --check 的等价断言）。"""
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "gen_cli_docs", ROOT / "scripts" / "gen_cli_docs.py"
    )
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    pages = module.render_all()
    assert len(pages) >= 12, f"命令页只有 {len(pages)} 个，是否漏了子命令？"
    for name, content in pages.items():
        target = ROOT / "docs" / "cli" / name
        assert target.exists(), f"缺 docs/cli/{name}：运行 python scripts/gen_cli_docs.py"
        assert target.read_text(encoding="utf-8") == content, (
            f"docs/cli/{name} 与 argparse 定义不一致：运行 python scripts/gen_cli_docs.py"
        )


def test_cli_docs_are_in_nav() -> None:
    """mkdocs strict 模式要求每页都在 nav 里，这里提前抓漏。"""
    nav_text = (ROOT / "mkdocs.yml").read_text(encoding="utf-8")
    pages = sorted(p.name for p in (ROOT / "docs" / "cli").glob("*.md"))
    missing = [name for name in pages if f"cli/{name}" not in nav_text]
    assert not missing, f"mkdocs.yml 的 nav 缺：{missing}"
