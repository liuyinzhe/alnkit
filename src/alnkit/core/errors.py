"""异常类型与退出码约定。

退出码沿用被重构脚本的既有约定：
0 成功（可有警告）/ 1 运行期错误 / 2 参数或前置校验失败。
"""

from __future__ import annotations


class AlnkitError(Exception):
    """所有可预期错误的基类；``exit_code`` 决定进程退出码。"""

    exit_code = 1


class UsageError(AlnkitError):
    """参数错误或前置条件不满足（退出码 2）。"""

    exit_code = 2


class RunError(AlnkitError):
    """运行期失败：读文件、外部命令、写盘（退出码 1）。"""

    exit_code = 1
