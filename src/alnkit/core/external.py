"""外部命令探测（minimap2 / samtools / bwa）。

只负责“找得到吗、什么版本”，命令拼装与执行留给各自的分析模块，
避免核心层耦合具体流程。
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections.abc import Sequence
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


#: 版本号形态（``0.7.19-r1273`` / ``2.31-r1302`` / ``1.24``）；用于从帮助文本里抠版本
_VERSION_PATTERN = re.compile(r"\b(\d+\.\d+(?:[.\-+][0-9A-Za-z]+)*)\b")


def _extract_version(text: str) -> str | None:
    """从命令行输出里取版本号。

    ``minimap2 --version`` / ``samtools --version`` 直接给一行版本；
    ``bwa --version`` 会报 "unrecognized command" 并把帮助打到 stderr，
    因此要退化成"跑一次无参调用，从帮助里找 Version 行"。
    """
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    for line in lines:
        if "version" in line.lower():
            match = _VERSION_PATTERN.search(line)
            if match:
                return match.group(1)
    # 没有 "Version:" 字样的（minimap2 / samtools）就取第一个像版本号的 token
    for line in lines:
        match = _VERSION_PATTERN.search(line)
        if match:
            return match.group(1)
    return None


def _run_probe(path: str, args: Sequence[str]) -> str | None:
    try:
        completed = subprocess.run(
            [path, *args],
            capture_output=True,
            text=True,
            timeout=_TIMEOUT_SECONDS,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout + completed.stderr


def _read_version(name: str, path: str) -> str | None:
    """探测版本：先 ``--version``；不支持时退化成无参调用（bwa 属于后者）。"""
    for args in (("--version",), ()):
        output = _run_probe(path, args)
        if not output:
            continue
        version = _extract_version(output)
        if version:
            return version
    return None


def find_tool(name: str) -> Tool:
    """探测一个命令；找不到时返回 ``available=False`` 而不是抛异常。"""
    path = shutil.which(name)
    if path is None:
        return Tool(name=name, path=None, version=None)
    return Tool(name=name, path=path, version=_read_version(name, path))
