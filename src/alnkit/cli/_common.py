"""CLI 公共件：解析器骨架、全局参数、统一错误出口。

参数命名规则：kebab-case 长参数，尽量与旧脚本同名（`--min-mapq`、`--merge-window-bp`…），
迁移期用户只需换命令名。退出码：0 成功 / 1 运行期错误 / 2 参数或前置校验失败。
"""

from __future__ import annotations

import argparse

from ..core.logging import configure

GLOBAL_OPTIONS = argparse.ArgumentParser(add_help=False)
GLOBAL_OPTIONS.add_argument(
    "--log-level",
    default="INFO",
    choices=("DEBUG", "INFO", "WARNING", "ERROR"),
    help="日志级别（默认 INFO）",
)
GLOBAL_OPTIONS.add_argument("--quiet", action="store_true", help="只输出错误")


def apply_global_options(args: argparse.Namespace) -> None:
    """把全局参数应用到日志系统。"""
    configure(level=getattr(args, "log_level", "INFO"), quiet=getattr(args, "quiet", False))
