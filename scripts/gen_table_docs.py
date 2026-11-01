#!/usr/bin/env python3
"""从 ``alnkit.core.tables`` 生成 ``docs/concepts/tables.md``（或校验是否最新）。

用法：
    python scripts/gen_table_docs.py           # 写入
    python scripts/gen_table_docs.py --check   # 只校验，过期则退出码 1
"""

from __future__ import annotations

import sys
from pathlib import Path

from alnkit._docs import render_table_docs

TARGET = Path(__file__).resolve().parent.parent / "docs" / "concepts" / "tables.md"


def main(argv: list[str]) -> int:
    rendered = render_table_docs()
    if "--check" in argv:
        current = TARGET.read_text(encoding="utf-8") if TARGET.exists() else ""
        if current != rendered:
            print(f"{TARGET} 已过期，请运行 python scripts/gen_table_docs.py", file=sys.stderr)
            return 1
        print(f"{TARGET} 最新")
        return 0
    TARGET.parent.mkdir(parents=True, exist_ok=True)
    TARGET.write_text(rendered, encoding="utf-8")
    print(f"已写入 {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
