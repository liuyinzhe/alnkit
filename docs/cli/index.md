# CLI 参考

> 本页与 `cli/*.md` 由 `scripts/gen_cli_docs.py` 从 argparse 定义生成，**不要手改**。
> 重新生成：`python scripts/gen_cli_docs.py`。

| 命令 | 说明 |
|---|---|
| [`alnkit doctor`](alnkit-doctor.md) | 检查运行环境（Python / 依赖 / 外部命令） |
| [`alnkit selftest`](alnkit-selftest.md) | 运行内置自检（不需要输入数据） |
| [`alnkit seg table`](alnkit-seg-table.md) | BAM -> 比对段表 |
| [`alnkit plot tracks`](alnkit-plot-tracks.md) | 单条 read 的比对片段轨道图 |
| [`alnkit plot ribbon`](alnkit-plot-ribbon.md) | 共线性丝带图（参考 vs 多个查询片段） |
| [`alnkit short scan`](alnkit-short-scan.md) | BAM -> 断点证据三张表 |
| [`alnkit short insertion-plot`](alnkit-short-insertion-plot.md) | 外源插入位点示意图（单侧 / 同断点双侧 / 多断点） |
| [`alnkit short coverage`](alnkit-short-coverage.md) | 区间覆盖广度与深度（PNG + 两张表） |
| [`alnkit long graph`](alnkit-long-graph.md) | read-vs-read overlap 边与连通分量 |
| [`alnkit long classify`](alnkit-long-classify.md) | 目标 reads 分类 + 接头事件/位点/区域聚类 |
| [`alnkit long breakpoints`](alnkit-long-breakpoints.md) | 断点图方案 A：紧邻片段断点 + 共享 reads 合并成类 |
| [`alnkit long assemble`](alnkit-long-assemble.md) | 按 overlap 关系贪心延伸组装每个类 |
| [`alnkit long locate`](alnkit-long-locate.md) | 50 bp 断点锚点 + 外源查询 -> 丝带图 + 参考坐标还原 |
