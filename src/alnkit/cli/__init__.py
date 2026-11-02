"""命令行入口：``alnkit <command>``。

6 个顶层命令 / 11 个叶子命令均已可用。
"""

from __future__ import annotations

import argparse
import sys
from collections.abc import Sequence

from .. import __version__
from ..core import errors
from . import doctor, long, plot, seg, selftest, short
from ._common import apply_global_options

_COMMANDS = (doctor, selftest, seg, plot, short, long)


def build_parser() -> argparse.ArgumentParser:
    """构建完整命令树（供 CLI 与文档生成共用）。"""
    parser = argparse.ArgumentParser(
        prog="alnkit",
        description="BAM/PAF 的插入证据工具箱：断点提取 + 独立绘图",
        epilog="例子：\n  alnkit doctor   |   alnkit selftest   |   alnkit <command> --help",
    )
    parser.add_argument("--version", action="version", version=f"alnkit {__version__}")
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")
    for module in _COMMANDS:
        module.add_parser(subparsers)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """CLI 主入口；返回进程退出码。"""
    parser = build_parser()
    args = parser.parse_args(argv)
    if getattr(args, "command", None) is None:
        parser.print_help()
        return 0
    apply_global_options(args)
    try:
        return int(args.handler(args))
    except errors.AlnkitError as exc:
        print(f"alnkit {args.command}: {exc}", file=sys.stderr)
        return exc.exit_code
    except KeyboardInterrupt:
        print("已中断", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
