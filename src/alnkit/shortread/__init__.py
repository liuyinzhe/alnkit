"""短读长链路（B1/B2）：BAM 里的插入断点证据（SR / PESR / DUP）。

模块分层（**依赖是 DAG，不是一条链**；低层不依赖高层）::

    model          共享数据类型与常量（对齐记录、命中、证据包、配置）
      ├── bam      CIGAR/比对解析、断点坐标、插入命中抽取、比对特征统计
      │     ├── interference   同源干扰源头判定
      │     └── evidence       SR / PESR / DUP 证据抽取与按位点聚合
      ├── merge    位点窗口合并（只用 model 的类型）
      ├── report   三张结果表写出（只用 model 的类型）
      └── origin   归属判定（用到 bam / evidence / interference）
                    ↓
                  scan   编排层（导入以上全部），并 re-export 它们的公开名

``model`` 不导入任何同包模块；``scan`` 是唯一的总装点。
"""

from __future__ import annotations

__all__ = ["scan"]
