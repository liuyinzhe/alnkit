# alnkit

BAM/PAF 的插入证据工具箱：断点提取 + 独立绘图。

三个支柱，互不依赖：

| 支柱 | 回答的问题 | 入口 |
|---|---|---|
| 长读长（三代） | 外源插在基因组的哪里、结构如何、能不能拼出来 | `alnkit long *` |
| 短读长（二代） | 断点在哪、有多少 SR/PESR/DUP 证据、有没有同源干扰 | `alnkit short scan` |
| 绘图 | 一条 read 上都有什么片段（三代）；一段序列出现在哪些组装的哪里（三代）；二代位点长什么样、覆盖够不够 | `alnkit plot *`（`tracks`/`ribbon` 面向三代，`insertion`/`coverage` 面向二代） |


## 三条链路

```text
三代：参考 + 外源 → 合并参考 → 比对 → 目标 read 分类 → overlap 图 → 断点/类 → 组装 → 落点(丝带图)
二代：BAM → SR/PESR/DUP 断点证据 → 归属/同源干扰 → 位点窗口合并 → 三张表
二代出图：BAM → 表（clip 证据 / 分箱统计）→ 示意图 / 覆盖度图（`plot insertion` / `plot coverage`）
绘图：任意 FASTA/PAF → 丝带图 ；任意 BAM/PAF/段表 → 单 read 轨道图
```

## 设计原则

1. **表是接口，图是渲染**：所有分析先落成 TSV，图从 TSV 出，可重画而不必重跑比对。
2. **单一来源**：同一功能只有一处实现；同一列定义只有一处（`alnkit.core.tables`）。
3. **契约先定**：坐标系、列名、列序、排序规则写死在核心层，模块不得自己拼表头。
4. **可独立使用**：`plot` 与 `seg table` 不依赖任何上游步骤产物。

## 下一步

- [安装](install.md)
- [快速开始](quickstart.md)
- [坐标约定](concepts/coordinates.md)
- [表契约](concepts/tables.md)
- [各命令的输出表](concepts/output-tables.md)
