"""CLI 公共件：解析器骨架、全局参数、统一错误出口。

参数命名规则：kebab-case 长参数，尽量与旧脚本同名（`--min-mapq`、`--merge-window-bp`…），
迁移期用户只需换命令名。退出码：0 成功 / 1 运行期错误 / 2 参数或前置校验失败。
"""

from __future__ import annotations

import argparse
import os

from ..core.errors import UsageError
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


def require_existing_files(*pairs) -> None:
    """检查「标签, 路径」列表里的输入文件是否存在；缺失时抛 ``UsageError``（退出码 2）。

    输入示例: require_existing_files(("--paf", "a.paf"), ("--classification", "b.tsv"))
    输出示例: 文件都在 -> 正常返回；缺 --paf -> UsageError("--paf 不存在：a.paf")

    说明: 不检查就会在更深处抛 FileNotFoundError/OSError，而 ``cli.main`` 只捕 ``AlnkitError``，
          用户看到的是一整段 traceback。
    """
    missing = [
        (label, path)
        for label, path in pairs
        if path is not None and not os.path.isfile(os.fspath(path))
    ]
    if missing:
        label, path = missing[0]
        raise UsageError(f"{label} 不存在或不是文件：{path}")
