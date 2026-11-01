"""统一日志：单点配置，避免各模块自建 handler。"""

from __future__ import annotations

import logging
import sys

LOGGER_NAME = "alnkit"
_FORMAT = "%(levelname)s %(name)s: %(message)s"


def configure(level: str = "INFO", quiet: bool = False) -> None:
    """配置 alnkit 根 logger。重复调用只更新级别，不叠加 handler。"""
    logger = logging.getLogger(LOGGER_NAME)
    logger.setLevel(logging.ERROR if quiet else level.upper())
    if not logger.handlers:
        handler = logging.StreamHandler(sys.stderr)
        handler.setFormatter(logging.Formatter(_FORMAT))
        logger.addHandler(handler)
    logger.propagate = False


def get_logger(name: str) -> logging.Logger:
    """取 ``alnkit.<name>`` 子 logger。"""
    return logging.getLogger(f"{LOGGER_NAME}.{name}")
