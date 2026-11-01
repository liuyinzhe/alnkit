"""``alnkit doctor``：打印运行环境，缺什么就说什么。"""

from __future__ import annotations

import argparse
import platform
from importlib.metadata import PackageNotFoundError, version

from .. import __version__
from ..core import errors
from ..core.external import ALIGNERS, REFERENCE_TOOLS, find_tool
from ._common import GLOBAL_OPTIONS

_REQUIRED_PACKAGES = ("pysam",)
_OPTIONAL_PACKAGES = (("matplotlib", "绘图子命令 alnkit plot *"),)


def add_parser(subparsers: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    parser = subparsers.add_parser(
        "doctor",
        parents=[GLOBAL_OPTIONS],
        help="检查运行环境（Python / 依赖 / 外部命令）",
        description="打印版本与可用性；必需项缺失时以退出码 1 结束。",
        epilog="例子：\n  alnkit doctor",
    )
    parser.set_defaults(handler=run)


def _package_version(name: str) -> str | None:
    try:
        return version(name)
    except PackageNotFoundError:
        return None


def run(args: argparse.Namespace) -> int:
    del args  # doctor 不使用全局参数以外的入参
    rows: list[tuple[str, str, str]] = []
    rows.append(("alnkit", __version__, "必需"))
    rows.append(("python", platform.python_version(), ">=3.10"))

    missing_required: list[str] = []
    for name in _REQUIRED_PACKAGES:
        found = _package_version(name)
        rows.append((name, found or "缺失", "必需"))
        if found is None:
            missing_required.append(name)

    missing_optional: list[str] = []
    for name, purpose in _OPTIONAL_PACKAGES:
        found = _package_version(name)
        rows.append((name, found or "缺失", f"可选：{purpose}"))
        if found is None:
            missing_optional.append(name)

    for name in ALIGNERS:
        tool = find_tool(name)
        rows.append((name, tool.describe(), "可选：比对（至少装一个）"))
    for name in REFERENCE_TOOLS:
        tool = find_tool(name)
        rows.append((name, tool.describe(), "可选：参考索引 / BAM 处理"))

    width = max(len(name) for name, _, _ in rows)
    print(f"alnkit {__version__}  平台 {platform.platform()}")
    print(f"{'组件':<{width}}  版本{'':<38}说明")
    for name, value, note in rows:
        print(f"{name:<{width}}  {value:<40}{note}")

    aligners = [name for name in ALIGNERS if find_tool(name).available]
    if missing_required:
        raise errors.RunError(
            "缺少必需依赖：" + ", ".join(missing_required) + "。安装：pip install alnkit"
        )
    if missing_optional:
        print(
            f"\n提示：未安装 {', '.join(missing_optional)}；"
            '绘图子命令需要 `pip install "alnkit[plot]"`。'
        )
    if not aligners:
        print(f"\n警告：{'/'.join(ALIGNERS)} 都未找到；分析链路需要用户自行提供 BAM/PAF。")
    return 0
