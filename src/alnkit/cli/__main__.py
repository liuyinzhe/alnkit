"""``python -m alnkit.cli`` 等价于 ``alnkit``（便于在子进程里调用，不必依赖 PATH）。"""

from __future__ import annotations

import sys

from . import main

if __name__ == "__main__":
    sys.exit(main())
