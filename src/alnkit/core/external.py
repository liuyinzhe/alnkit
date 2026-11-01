"""外部命令探测（minimap2 / samtools / bwa）。

只负责“找得到吗、什么版本”，命令拼装与执行留给各自的分析模块，
避免核心层耦合具体流程。
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass

__all__ = ["Tool", "ALIGNERS", "REFERENCE_TOOLS", "find_tool"]

#: 比对器（至少装一个）
ALIGNERS = ("minimap2", "bwa")
#: 参考/比对文件处理
REFERENCE_TOOLS = ("samtools",)

_TIMEOUT_SECONDS = 10


@dataclass(slots=True, frozen=True)
class Tool:
    """一个外部命令的探测结果。"""

    name: str
    path: str | None
    version: str | None

    @property
    def available(self) -> bool:
        return self.path is not None

    def describe(self) -> str:
        if not self.available:
            return "未安装"
        return self.version or "已安装（版本未知）"


def _read_version(name: str, path: str) -> str | None:
    try:
        completed = subprocess.run(
            [path, "--version"],
            capture_output=True,
            text=True,
            timeout=_TIMEOUT_SECONDS,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    for line in (completed.stdout + completed.stderr).splitlines():
        if line.strip():
            return line.strip()[:120]
    return None


def find_tool(name: str) -> Tool:
    """探测一个命令；找不到时返回 ``available=False`` 而不是抛异常。"""
    path = shutil.which(name)
    if path is None:
        return Tool(name=name, path=None, version=None)
    return Tool(name=name, path=path, version=_read_version(name, path))
