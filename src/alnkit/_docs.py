"""文档生成：表契约 → markdown（唯一来源是 ``core.tables``）。

放成私有模块（而不是 ``scripts/`` 里的独立脚本），是为了让测试能直接 import 并校验
“文档与代码一致”，不必做 sys.path 手脚。
"""

from __future__ import annotations

from .core import tables

_HEADER = """# 表契约

> 本页由 `alnkit.core.tables` 的 `COLUMNS` / `COLUMN_DOC` 生成，**不要手改**。
> 重新生成：`python scripts/gen_table_docs.py`。

规则：

- 列序即下表顺序，写表只能通过 `alnkit.core.tables.write_tsv`；
- 坐标列名必须带 `_0based` / `_1based` 后缀；长度列（`*_bp`、`*_len`、`*_length`、`*_count`）不带；
- 旧表（legacy `.xls` / 旧 TSV）列名保持不变，由 `--columns legacy:*` 复现；
- 浮点固定 4 位小数；空值写 `-`。
"""


#: 注册了但还没有生产者的表（文档里要写明，避免用户去找不存在的文件）
NOT_PRODUCED: frozenset[str] = frozenset({"ribbon_blocks"})


def render_table_docs() -> str:
    """渲染 docs/concepts/tables.md 的完整内容。"""
    chunks = [_HEADER]
    for table, columns in tables.COLUMNS.items():
        docs = tables.COLUMN_DOC.get(table, {})
        chunks.append(f"\n## `{table}`\n")
        if table in NOT_PRODUCED:
            chunks.append(
                '!!! warning "尚未实现"\n'
                "    目前没有命令输出这张表（`plot ribbon` 仍按旧表头写 blocks.tsv）；\n"
                "    列定义先冻结在这里，等 C1 的口径统一切换。\n"
            )
        chunks.append("| # | 列名 | 基准 | 说明 |")
        chunks.append("|---|---|---|---|")
        for index, column in enumerate(columns, start=1):
            if column.endswith("_0based"):
                basis = "0-based"
            elif column.endswith("_1based"):
                basis = "1-based"
            elif tables.is_coordinate_column(column):  # pragma: no cover - 注册时已拦住
                basis = "**缺失**"
            else:
                basis = "—"
            chunks.append(f"| {index} | `{column}` | {basis} | {docs.get(column, '')} |")
    return "\n".join(chunks) + "\n"
