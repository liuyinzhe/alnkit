"""测试共用的小工具。"""

from __future__ import annotations

from pathlib import Path

#: 黄金目录里的元数据文件，不属于产物，比对时要排除
METADATA_NAMES = frozenset({"PROVENANCE.txt"})


def output_files(directory: Path, *, exclude_metadata: bool = False) -> list[str]:
    """列出目录下所有文件的相对路径（排序后），便于做"产物集合"断言。"""
    return sorted(
        str(path.relative_to(directory))
        for path in directory.rglob("*")
        if path.is_file() and not (exclude_metadata and path.name in METADATA_NAMES)
    )


def golden_files(directory: Path) -> list[str]:
    """黄金目录的产物清单（排除 PROVENANCE.txt）。"""
    return output_files(directory, exclude_metadata=True)
