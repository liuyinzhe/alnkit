"""脱敏守卫：仓库里不得出现真实样本名、构建体名或本机绝对路径。

真实名字只存在本机私有目录的 `manifest.json`（`exog_names` / `query_chr` / `insert_chrs`），
仓库统一用占位符（`sampleA/sampleB`、`target_ref/homolog_ref/decoy_ref`、`insertA/insertB`）。
这条测试把"名单"从私有目录读出来再扫仓库，所以**名单本身不会进仓库**。

没有私有数据时（`ALNKIT_TEST_DATA` 未设置）跳过。
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pytest

DATA = Path(os.environ["ALNKIT_TEST_DATA"]) if os.environ.get("ALNKIT_TEST_DATA") else None
pytestmark = pytest.mark.skipif(
    DATA is None or not DATA.exists(),
    reason="需要 ALNKIT_TEST_DATA（私有目录里的 manifest 提供真实名单）",
)

REPO = Path(__file__).resolve().parents[1]

#: 扫描的文本文件（只扫会随包/仓库分发的内容）；.tsv 覆盖 fixtures 里的黄金表
SCAN_SUFFIXES = {".py", ".md", ".toml", ".yml", ".yaml", ".txt", ".cfg", ".cff", ".tsv"}
SKIP_DIRS = {
    ".git",
    "site",
    "dist",
    "build",
    ".pytest_cache",
    ".mypy_cache",
    ".ruff_cache",
    ".venv",
}

#: 本机路径前缀（在文档里写死会让别人照抄失败，也泄漏目录结构）
LOCAL_PATH_PREFIXES = ("/home/askr", "C:\\Users\\")


def _private_names() -> set[str]:
    names: set[str] = set()
    for manifest in DATA.rglob("manifest.json"):
        payload = json.loads(manifest.read_text(encoding="utf-8"))
        for key in (
            "exog_names",
            "insert_chrs",
            "other_inserts",
            "query_chr",
            "sample",
            "sample_name",
        ):
            value = payload.get(key)
            if isinstance(value, str):
                names.add(value)
            elif isinstance(value, list):
                names.update(item for item in value if isinstance(item, str))
    return {name for name in names if len(name) >= 4}


def _iter_files():
    for path in REPO.rglob("*"):
        if path.is_dir() or path.suffix.lower() not in SCAN_SUFFIXES:
            continue
        if any(part in SKIP_DIRS for part in path.relative_to(REPO).parts):
            continue
        if path.resolve() == Path(__file__).resolve():
            continue  # 本文件自己定义禁词模式，跳过
        yield path


def test_repository_has_no_private_names() -> None:
    """真实样本/构建体名不得出现在仓库任何文本文件里。"""
    names = _private_names()
    assert names, "私有 manifest 里没读到任何真实名字，检查 ALNKIT_TEST_DATA 是否完整"
    leaks = []
    for path in _iter_files():
        text = path.read_text(encoding="utf-8", errors="ignore")
        for name in names:
            if name in text:
                leaks.append(f"{path.relative_to(REPO)}: {name}")
    assert not leaks, "发现私有名字泄漏：\n" + "\n".join(leaks)


def test_repository_has_no_local_absolute_paths() -> None:
    """文档/代码里不得写死本机绝对路径（别人照抄会失败）。"""
    leaks = []
    for path in _iter_files():
        text = path.read_text(encoding="utf-8", errors="ignore")
        for prefix in LOCAL_PATH_PREFIXES:
            if prefix in text:
                leaks.append(f"{path.relative_to(REPO)}: {prefix}")
    assert not leaks, "发现本机绝对路径：\n" + "\n".join(leaks)
