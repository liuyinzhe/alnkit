# 坐标约定

!!! note "事实依据"
    本页结论来自对三条链路产物的实测审计（旧项目 + BAM/PAF 反查），不是推测。
    审计记录：机器列除二代 `Insert_Position` 外均为 0-based 半开；左右剪两个分支共 303 条比对验证。

## 三条硬规则

1. **内部计算一律 0-based 半开** `[start, end)`，长度 = `end - start`。
2. **输出列名必须带基准后缀**：`_0based` 或 `_1based`。长度列（`*_bp`、`*_len`、`*_length`、`*_count`）与
   `*_side` / `*_type` / `*_id` 不带后缀，因为它们不是坐标。
3. **区间展示串一律 1-based 闭区间** `chrom:start-end`（IGV / samtools / VCF 习惯）。

```python
from alnkit.core.coords import format_region_1based, parse_region_1based

format_region_1based("chr1", 2000000, 2000050)  # 'chr1:2000001-2000050'
parse_region_1based("chr1:2000001-2000050")  # ('chr1', 2000000, 2000050)
```

不要手写 `+1` / `-1`，只用 `alnkit.core.coords` 的入口（`to_1based` / `to_0based` /
`interval_0based_to_1based` / `interval_1based_to_0based` / `left_base_1based` / `left_base_0based`）。

## read 轴

read 坐标永远是 **read 自身 5'→3' 正向**（反向链比对段已翻转）。`strand` 只描述比对方向，不改坐标。

## 断点这类“单点”

历史产物（二代主表 `Insert_Position`）用的是**边界左侧那个碱基的 1-based 坐标**。
它与 0-based 边界数值相同，容易误读，所以规范表同时给三列：

| 列 | 含义 | 例（chr1 上某插入） |
|---|---|---|
| `insert_pos_1based` | 断点左侧碱基的 1-based 坐标（= legacy `Insert_Position`，不改值） | 1000001 |
| `insert_pos_0based` | 同一碱基的 0-based 坐标 | 1000000 |
| `insert_boundary_0based` | 断点边界（0-based，= 左侧碱基数） | 1000001 |

恒等式：`insert_pos_1based - 1 == insert_pos_0based`，且 `insert_boundary_0based == insert_pos_1based`。
最后一条是历史巧合，**不是**笔误——代码里有断言守住。

## 旧表怎么办

legacy 产物（二代三张 `.xls`、两套旧段表）**列名与数值一律不改**，用于黄金对账与既有下游。
新表用规范列名，两者可通过下表的别名对应：

| 旧列名 | 规范列名 | 基准 |
|---|---|---|
| `ref_start` / `ref_pos_start` / `target_start` / `segment_start` | `ref_start_0based` | 0-based 半开 |
| `ref_end` / `ref_pos_end` / `target_end` / `segment_end` | `ref_end_0based` | 0-based 半开 |
| `read_start` / `read_pos_start` / `M_read_pos_start` | `read_start_0based` | 0-based 半开，正向 read 轴 |
| `read_end` / `read_pos_end` / `M_read_pos_end` | `read_end_0based` | 0-based 半开，正向 read 轴 |
| `end_coordinate` / `adjacent_coordinate` | `end_boundary_0based` | 0-based 边界 |
| `boundary_min/max/median` | `boundary_min_0based` / `_max_` / `_median_` | 0-based 边界 |
| `window_50bp`（串） | `window_display_1based` + `window_start_0based` / `window_end_0based` | 串 1-based 闭；新列 0-based 半开 |
| `M_var_len` | `aln_block_bp` | 长度 |
| `Insert_Position` | `insert_pos_1based`（同值）+ `insert_pos_0based` + `insert_boundary_0based` | 见上节 |

## 两个容易踩的坑

- **`end_side` ≠ 参考侧**：三代段末端表里，边界取哪一端由 `genome_side` 决定
  （`genome_left` → 段末端，`genome_right` → 段起点）；`end_side` 只是 read 轴的左/右。
  两者混用会差一个整段。
- **`ref_aln_len` 与 CIGAR 的 M 之和可以不等**：前者是真实参考跨度（含首尾 D/N），
  后者来自“把 ≤5 bp 的 I/D 并进 M”的简化视图。两者都被保留，别互相校验。
