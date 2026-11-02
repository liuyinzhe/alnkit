#!/usr/bin/env python3
"""从 argparse 命令树生成 ``docs/cli/*.md``（CLI 参考页）。

参数表不手写：唯一来源是 ``alnkit.cli.build_parser()`` 的定义，
所以文档永远与 `--help` 一致。CI 用 ``--check`` 校验"文档是否最新"。

用法：
    python scripts/gen_cli_docs.py           # 写入
    python scripts/gen_cli_docs.py --check   # 只校验，过期则退出码 1
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DOCS = ROOT / "docs" / "cli"
MKDOCS = ROOT / "mkdocs.yml"


def _subparsers(
    parser: argparse.ArgumentParser,
) -> dict[str, tuple[argparse.ArgumentParser, str]]:
    """{子命令名: (子 parser, 一行 help)}。

    argparse 把 `add_parser(..., help=...)` 的文本存在 `_ChoicesPseudoAction` 里，
    并不挂在子 parser 上，所以要从 subparsers action 里取。
    """
    for action in parser._actions:  # noqa: SLF001
        if isinstance(action, argparse._SubParsersAction):  # noqa: SLF001
            helps = {pa.dest: pa.help for pa in getattr(action, "_choices_actions", [])}
            return {name: (sub, helps.get(name) or "") for name, sub in action.choices.items()}
    return {}


def _leaf_commands(
    parser: argparse.ArgumentParser, prefix: tuple[str, ...] = (), summary: str = ""
) -> list[tuple[tuple[str, ...], argparse.ArgumentParser, str]]:
    children = _subparsers(parser)
    if not children:
        return [(prefix, parser, summary)]
    leaves: list[tuple[tuple[str, ...], argparse.ArgumentParser, str]] = []
    for name, (child, helper) in children.items():
        leaves.extend(_leaf_commands(child, (*prefix, name), helper))
    return leaves


def _format_default(action: argparse.Action) -> str:
    default = action.default
    if isinstance(default, Path):
        # 默认值是 Path.cwd() 的参数不能让文档跟着"在哪台机器上生成"变
        return "`当前目录`" if default == Path.cwd() else f"`{default}`"
    if isinstance(action, argparse._StoreFalseAction):  # noqa: SLF001
        # --no-xxx 这类开关：默认 True 表示“功能默认开着”，但开关本身默认是关闭的
        return "关" if default else "开"
    if action.nargs == 0 or isinstance(action, argparse._StoreTrueAction):  # noqa: SLF001
        return "关" if not default else "开"
    if default in (None, "", []):
        return "—"
    if default is argparse.SUPPRESS or default is sys.modules.get("argparse").SUPPRESS:
        return "—"
    return f"`{default}`"


#: 通过 parents=[GLOBAL_OPTIONS] 继承来的选项，单独列在“全局选项”一节
_GLOBAL_DESTS = frozenset({"log_level", "quiet"})


def _argument_row(action: argparse.Action) -> str | None:
    if action.dest in _GLOBAL_DESTS or action.dest == "help" or not action.option_strings:
        return None
    names = ", ".join(f"`{name}`" for name in action.option_strings)
    metavar = getattr(action, "metavar", None)
    if metavar and not isinstance(action, argparse._StoreTrueAction):  # noqa: SLF001
        names += f" {metavar}"
    required = "是" if action.required else "否"
    help_text = (action.help or "").replace("|", "\\|").replace("\n", " ")
    if action.choices:
        help_text += "（可选：" + ", ".join(f"`{c}`" for c in action.choices) + "）"
    return f"| {names} | {required} | {_format_default(action)} | {help_text} |"


def _usage_line(prog: str, parser: argparse.ArgumentParser) -> str:
    parts = [prog, "[全局选项]"]
    positionals = [a for a in parser._actions if not a.option_strings and a.dest != "help"]  # noqa: SLF001
    for action in positionals:
        parts.append(f"<{action.dest}>" if action.required else f"[<{action.dest}>]")
    options = [
        a
        for a in parser._actions  # noqa: SLF001
        if a.option_strings
        and a.dest not in _GLOBAL_DESTS
        and a.dest != "help"
        and not isinstance(a, argparse._StoreTrueAction)  # noqa: SLF001
    ]
    for action in options:
        parts.append(f"[{action.option_strings[0]} …]")
    flags = [
        a
        for a in parser._actions  # noqa: SLF001
        if isinstance(a, argparse._StoreTrueAction) and a.dest not in _GLOBAL_DESTS
    ]
    for action in flags:
        parts.append(f"[{action.option_strings[0]}]")
    return " ".join(parts)


def render_command(path: tuple[str, ...], parser: argparse.ArgumentParser) -> str:
    """渲染一个叶子命令页。"""
    prog = " ".join(("alnkit", *path))
    lines = [f"# `{prog}`", "", (parser.description or "").strip(), ""]

    lines += ["## 用法", "", "```bash", _usage_line(prog, parser), "```", ""]

    rows = [row for action in parser._actions if (row := _argument_row(action))]  # noqa: SLF001
    if rows:
        lines += ["## 参数", "", "| 参数 | 必填 | 默认 | 说明 |", "|---|---|---|---|", *rows, ""]

    global_rows = [
        ("`--log-level` LEVEL", "否", "`INFO`", "日志级别：DEBUG / INFO / WARNING / ERROR"),
        ("`--quiet`", "否", "关", "只输出错误（等价 --log-level ERROR）"),
    ]
    lines += [
        "## 全局选项（每个子命令都有）",
        "",
        "| 参数 | 必填 | 默认 | 说明 |",
        "|---|---|---|---|",
        *(f"| {name} | {req} | {default} | {text} |" for name, req, default, text in global_rows),
        "",
    ]

    epilog = (parser.epilog or "").strip()
    if epilog:
        examples = epilog.split("\n\n")
        lines += ["## 例子", ""]
        for block in examples:
            if "例子" in block.split("\n")[0] and block.startswith("例子"):
                body = block.split("\n", 1)[1] if "\n" in block else ""
                # 保留行尾续行符：去掉的话读者复制多行命令会直接失败
                lines += ["```bash", body.strip(), "```", ""]
            else:
                lines += ["```text", block, "```", ""]

    lines += [
        "## 退出码",
        "",
        "`0` 成功（可有警告）/ `1` 运行期错误 / `2` 参数或前置校验失败。",
        "",
    ]
    return "\n".join(lines).rstrip() + "\n"


def render_index(
    leaves: list[tuple[tuple[str, ...], argparse.ArgumentParser, str]],
) -> str:
    lines = [
        "# CLI 参考",
        "",
        "> 本页与 `cli/*.md` 由 `scripts/gen_cli_docs.py` 从 argparse 定义生成，**不要手改**。",
        "> 重新生成：`python scripts/gen_cli_docs.py`。",
        "",
        "| 命令 | 说明 |",
        "|---|---|",
    ]
    for path, parser, helper in leaves:
        summary = (helper or parser.description or "").strip().split("\n")[0]
        page = "-".join(("alnkit", *path))
        lines.append(f"| [`alnkit {' '.join(path)}`]({page}.md) | {summary} |")
    return "\n".join(lines) + "\n"


def render_all() -> dict[str, str]:
    """返回 {相对 docs/cli 的文件名: 内容}。"""
    parser = __import__("alnkit.cli", fromlist=["build_parser"]).build_parser()
    leaves = _leaf_commands(parser, ())
    pages = {"index.md": render_index(leaves)}
    for path, leaf_parser, _helper in leaves:
        pages["-".join(("alnkit", *path)) + ".md"] = render_command(path, leaf_parser)
    return pages


def _nav_check(pages: dict[str, str]) -> list[str]:
    """mkdocs strict 模式要求每页都在 nav 里，这里提前抓漏。"""
    if not MKDOCS.exists():
        return []
    nav_text = MKDOCS.read_text(encoding="utf-8")
    missing = [name for name in pages if f"cli/{name}" not in nav_text]
    return missing


def main(argv: list[str]) -> int:
    pages = render_all()
    problems: list[str] = []

    if "--check" in argv:
        for name, content in pages.items():
            target = DOCS / name
            if not target.exists() or target.read_text(encoding="utf-8") != content:
                problems.append(f"docs/cli/{name} 已过期或缺失")
    else:
        DOCS.mkdir(parents=True, exist_ok=True)
        for name, content in pages.items():
            (DOCS / name).write_text(content, encoding="utf-8")
        stale = {p.name for p in DOCS.glob("*.md")} - set(pages)
        for name in sorted(stale):
            (DOCS / name).unlink()
            print(f"删除不再存在的命令页：docs/cli/{name}")
        print(f"已写入 {len(pages)} 页到 {DOCS}")

    for name in _nav_check(pages):
        problems.append(f"mkdocs.yml 的 nav 缺 docs/cli/{name}")

    if problems:
        for problem in problems:
            print(f"✗ {problem}", file=sys.stderr)
        print("修复：python scripts/gen_cli_docs.py", file=sys.stderr)
        return 1
    print(f"✓ CLI 文档最新（{len(pages)} 页）")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
