"""表契约：列序、列说明、读写。

规则：

* 列序只在这里定义，任何模块不得自己拼表头；
* 坐标列名**必须**带 ``_0based`` / ``_1based`` 后缀（长度列、``*_side`` 等除外）；
* 写入原子化，行序由调用方显式给出；
* 浮点固定 4 位小数（与历史产物一致，改格式会破坏黄金对账）。
"""

from __future__ import annotations

import csv
import os
import shutil
from collections.abc import Callable, Iterable, Mapping, Sequence
from types import TracebackType
from typing import Any, Literal, TextIO

__all__ = [
    "COLUMNS",
    "COLUMN_DOC",
    "columns_of",
    "is_coordinate_column",
    "read_legacy_tsv",
    "read_tsv",
    "write_legacy_tsv",
    "validate_column_names",
    "validate_columns",
    "write_tsv",
]

#: 判定“这是坐标列”的词根（按 ``_`` 切分后逐段比较）
COORDINATE_TOKENS = frozenset({"start", "end", "pos", "position", "coordinate", "boundary"})

#: 坐标列必须以此结尾
BASIS_SUFFIXES = ("_0based", "_1based")

#: 一定不是坐标的列后缀（方向、类别、编号）
NON_COORDINATE_SUFFIXES = ("_side", "_type", "_name", "_id")
#: 长度类后缀：与坐标词根同时出现时是可疑命名（如 ``ref_end_bp``），仍按坐标列要求基准后缀
LENGTH_SUFFIXES = ("_len", "_length", "_bp", "_count")

COLUMNS: dict[str, tuple[str, ...]] = {
    # 规范化的“比对段表”（三套历史表头的统一出口）
    "segments": (
        "read_name",
        "read_length",
        "segment_role",
        "chrom",
        "chrom_length",
        "ref_start_0based",
        "ref_end_0based",
        "read_start_0based",
        "read_end_0based",
        "strand",
        "mapq",
        "ref_aln_len",
        "read_aln_len",
        "cigar",
        "kept",
        "drop_reason",
        "is_secondary",
        "is_supplementary",
        "nm",
        "match_columns",
        "inserted_bases",
        "deleted_bases",
    ),
    # 二代位点规范表（机器接口）：`short scan` 与三张 legacy .xls 并行输出。
    "insert_sites": (
        "chrom",
        "insert_pos_1based",
        "insert_pos_0based",
        "insert_boundary_0based",
        "n_sr",
        "n_pesr",
        "n_dup",
        "sample",
        "query_insert",
        "nearest_other_insert",
        "other_inserts",
    ),
    # 同源位点诊断：位点是否落在“与外源同源的基因组区域”里（外源回比 / 已知 BED）
    "site_homology": (
        "sample",
        "chrom",
        "insert_pos_1based",
        "is_homology_site",
        "n_regions",
        "region_sources",
        "best_region_chrom",
        "best_region_start_1based",
        "best_region_end_1based",
        "best_region_length_bp",
        "best_region_identity",
        "best_region_query_coverage",
    ),
    # 组装支撑度（A 方案）：把类内 reads 用 minimap2 打回 contig 的可复核证据
    "assembly_support": (
        "class_id",
        "contig",
        "contig_length_bp",
        "n_reads_total",
        "n_reads_mapped",
        "covered_bp",
        "covered_fraction",
        "mean_depth",
        "depth_p10",
        "depth_cliffs",
        "ends_covered",
        "status",
    ),
    # 丝带图映射块。注意：这里登记的是**规范名**，而 `alnkit plot ribbon` 目前仍按旧表头
    # （query_start / ref_start …）写出——迁移期有意保留，等 C1 的黄金口径统一后再切。
    "ribbon_blocks": (
        "query_fragment",
        "query_start_0based",
        "query_end_0based",
        "strand",
        "reference",
        "ref_start_0based",
        "ref_end_0based",
        "block_bp",
        "identity",
    ),
    # 外源插入位点示意图（`plot insertion`）的逐 read 证据：唯一绘图输入。
    "insertion_clip_reads": (
        "sample",
        "chrom",
        "breakpoint_0based",
        "clip_side",
        "read_name",
        "read_length",
        "strand",
        "mapq",
        "ref_start_0based",
        "ref_end_0based",
        "clip_length_bp",
        "clip_sequence",
    ),
    # 同一批 clip read 逐列投票的结果：每个断点每侧一行（无 read 的侧不出现）。
    "insertion_clip_consensus": (
        "sample",
        "chrom",
        "breakpoint_0based",
        "breakpoint_1based",
        "clip_side",
        "n_reads",
        "consensus_bp",
        "consensus_fraction",
        "stop_reason",
        "consensus_sequence",
    ),
    # 覆盖度：区域一行汇总（图标题/副标题里的 reads、mapQ、全区间均值都取自这里）。
    "coverage_region": (
        "sample",
        "chrom",
        "region_start_1based",
        "region_end_1based",
        "region_length_bp",
        "n_bins",
        "n_reads_total",
        "n_reads_selected",
        "n_reads_filtered",
        "mean_mapq",
        "mean_baseq",
        "covered_bases",
        "percent_covered",
        "mean_depth",
    ),
    # 覆盖度：逐 bin 的广度与深度（`plot coverage` 的两个绘图输入之一）。
    "coverage_bins": (
        "bin_index",
        "bin_start_1based",
        "bin_end_1based",
        "bin_width_bp",
        "covered_bases",
        "depth_sum",
    ),
    # 覆盖度：逐碱基深度（`--bases-table` 才写，默认不产出）。
    "coverage_bases": ("pos_1based", "depth"),
}

COLUMN_DOC: dict[str, dict[str, str]] = {
    "segments": {
        "read_name": "read 名（与 BAM/FASTQ 一致）",
        "read_length": "read 长度（bp）；未知时留空",
        "segment_role": "段角色：genome / exog / vector / target",
        "chrom": "参考序列名",
        "chrom_length": "参考序列长度（bp）",
        "ref_start_0based": "参考起点，0-based 半开",
        "ref_end_0based": "参考终点，0-based 半开（不含）",
        "read_start_0based": "read 轴起点，已翻到正向 read，0-based 半开",
        "read_end_0based": "read 轴终点，已翻到正向 read，0-based 半开（不含）",
        "strand": "比对方向：+ / -",
        "mapq": "MAPQ",
        "ref_aln_len": "参考跨度（含首尾 D/N），与简化 CIGAR 的 M 之和可不等",
        "read_aln_len": "read 跨度",
        "cigar": "简化 CIGAR（两个 M 之间总长 <=5 bp 的 I/D 并入 M）",
        "kept": "1/0：是否参与分类；只有 `long classify` 的产物有真实含义，`seg table` 恒为 1",
        "drop_reason": "被丢弃原因；未丢弃填 -",
        "is_secondary": "是否 secondary 比对（1/0）",
        "is_supplementary": "是否 supplementary 比对（1/0）",
        "nm": "编辑距离（NM 标签；缺失时为空）",
        "match_columns": "CIGAR 里 M/=/X 的碱基总数",
        "inserted_bases": "CIGAR 里 I 的碱基数",
        "deleted_bases": "CIGAR 里 D 的碱基数",
    },
    "insert_sites": {
        "chrom": "参考序列名",
        "insert_pos_1based": "断点左侧碱基的 1-based 坐标（= legacy Insert_Position）",
        "insert_pos_0based": "同一碱基的 0-based 坐标（= insert_pos_1based - 1）",
        "insert_boundary_0based": "断点边界，0-based（数值等于 insert_pos_1based）",
        "n_sr": "SR 证据数",
        "n_pesr": "PESR 证据数",
        "n_dup": "该位点上被 markdup 标记的 read 数",
        "sample": "样本名（BAM 文件名首个 . 之前）",
        "query_insert": "本次运行针对的插入片段（= --query-chr）",
        "nearest_other_insert": (
            "除 query 外、支持该位点 read 命中最多的插入片段"
            "（并列取名字小的）；没有其它插入片段命中时填 -"
        ),
        "other_inserts": "去掉 query 与最近者之后，其余命中过的插入片段，`,` 分隔；无则 -",
    },
    "site_homology": {
        "sample": "输出前缀（--sample-name）",
        "chrom": "位点所在参考序列",
        "insert_pos_1based": "断点左侧碱基的 1-based 坐标（= Insert_Position 列）",
        "is_homology_site": "yes/no：是否落在一条同源区域里（两侧各放宽 50 bp）",
        "n_regions": "命中的同源区域个数",
        "region_sources": "命中的来源：insert:<插入名> 或 bed:<名字>，多个用分号连接",
        "best_region_chrom": "最可信的同源区域所在参考序列（无命中写 -）",
        "best_region_start_1based": "最可信同源区域的起始（1-based，闭区间）",
        "best_region_end_1based": "最可信同源区域的结束（1-based，闭区间）",
        "best_region_length_bp": "最可信同源区域的长度（bp）",
        "best_region_identity": "回比得到的一致率（0-1，按 cs 标签的真实差异算）；BED 来源写 -",
        "best_region_query_coverage": "插入序列被这段区域覆盖的比例（0-1）；BED 来源写 -",
    },
    "assembly_support": {
        "class_id": "类 ID（= 目录名）",
        "contig": "contig 名（contigs.fasta 里的记录名）",
        "contig_length_bp": "contig 长度（bp）",
        "n_reads_total": "该类 reads 总数（FASTQ 行数 / 4）",
        "n_reads_mapped": "能比回该 contig 的 read 数（MAPQ >= 10）",
        "covered_bp": "被至少一条 read 覆盖的碱基数（合并区间后）",
        "covered_fraction": "covered_bp / contig_length_bp，0-1 浮点",
        "mean_depth": "比对碱基总数 / contig 长度",
        "depth_p10": "非零窗口深度的 10% 分位（窗口 200 bp）",
        "depth_cliffs": "相邻 200 bp 窗口深度降到前一个 1/5 以下的次数（chimera 信号）",
        "ends_covered": "首尾各 100 bp 中是否至少有一半被 read 覆盖（yes/no）",
        "status": "ok / low_coverage / chimera_suspect / unsupported_ends / no_reads",
    },
    "ribbon_blocks": {
        "query_fragment": "查询片段名",
        "query_start_0based": "查询片段自身坐标起点，0-based 半开",
        "query_end_0based": "查询片段自身坐标终点，0-based 半开（不含）",
        "strand": "查询相对参考的方向：+ / -",
        "reference": "参考序列名",
        "ref_start_0based": "参考起点，0-based 半开",
        "ref_end_0based": "参考终点，0-based 半开（不含）",
        "block_bp": "块长度（= ref_end - ref_start）",
        "identity": "一致率，0-1 浮点，固定 4 位小数",
    },
    "insertion_clip_reads": {
        "sample": "样本名（BAM 文件名首个 . 之前的部分）",
        "chrom": "参考序列名",
        "breakpoint_0based": "断点边界，0-based（数值等于断点左侧碱基的 1-based 坐标）",
        "clip_side": "clip 在比对块的哪一侧：left / right",
        "read_name": "read 名（与 BAM 一致）",
        "read_length": "read 长度（含 soft clip 的碱基，bp）",
        "strand": "比对方向：+ / -",
        "mapq": "MAPQ",
        "ref_start_0based": "比对块起点，0-based 半开",
        "ref_end_0based": "比对块终点，0-based 半开（不含）",
        "clip_length_bp": "soft clip 长度（bp）",
        "clip_sequence": (
            "clip 碱基，与参考同向；left = 比对块左侧那一段（反向链 read 上是 read 的 3' 端）"
        ),
    },
    "insertion_clip_consensus": {
        "sample": "样本名（BAM 文件名首个 . 之前的部分）",
        "chrom": "参考序列名",
        "breakpoint_0based": "断点边界，0-based",
        "breakpoint_1based": "同一个断点，按“左侧碱基的 1-based 坐标”写（数值与 0based 列相同）",
        "clip_side": "left / right",
        "n_reads": "该侧参与投票的 clip read 数（候选数，不代表每个碱基都投了票）",
        "consensus_bp": "一致性序列长度（bp）",
        "consensus_fraction": "一致率 0-1（每列多数碱基数之和 / 全部条数之和）",
        "stop_reason": (
            "中断原因：all columns / coverage<2 (only N reads) / 2-read overlap mismatch / no clip"
        ),
        "consensus_sequence": (
            "一致性序列（该侧无序列写 -）；left 侧以接头端为第 0 列、向断点外侧延伸（参考方向）"
        ),
    },
    "coverage_region": {
        "sample": "样本名（BAM 文件名首个 . 之前的部分）",
        "chrom": "参考序列名",
        "region_start_1based": "统计区间起点，1-based 闭区间",
        "region_end_1based": "统计区间终点，1-based 闭区间",
        "region_length_bp": "区间长度（bp）",
        "n_bins": "分箱个数（= coverage_bins 的行数）",
        "n_reads_total": "fetch 到的比对记录数（含被过滤掉的）",
        "n_reads_selected": "通过 flag / mapQ / 长度过滤的 read 数",
        "n_reads_filtered": "被过滤掉的 read 数（= 总数 - 选中数）",
        "mean_mapq": "选中 read 的平均比对质量",
        "mean_baseq": "被覆盖位置上的平均碱基质量（只出现在终端摘要，图中不显示）",
        "covered_bases": "被覆盖的碱基数（深度 >= min_depth 的位置数）",
        "percent_covered": "覆盖广度百分比 0-100（= covered_bases / region_length_bp）",
        "mean_depth": "全区间平均深度（分子只累加被覆盖位置，分母是区间长度）",
    },
    "coverage_bins": {
        "bin_index": "第几个 bin，从 1 开始",
        "bin_start_1based": "bin 起点，1-based 闭区间",
        "bin_end_1based": "bin 终点，1-based 闭区间（相邻 bin 首尾相接）",
        "bin_width_bp": "bin 宽度（bp）；最后一个 bin 吸收除不尽的余数",
        "covered_bases": "bin 内被覆盖的碱基数",
        "depth_sum": "bin 内逐碱基深度之和（除以宽度就是该 bin 的平均深度）",
    },
    "coverage_bases": {
        "pos_1based": "碱基坐标，1-based",
        "depth": "该位置的深度（已按 --max-depth 封顶）",
    },
}


def is_coordinate_column(name: str) -> bool:
    """按词根判断是否是坐标列。

    ``read_start_0based`` 是；``end_side`` / ``junction_type`` 不是（类别后缀优先）；
    ``ref_end_bp`` 这种“坐标词根 + 长度后缀”的混搭**仍按坐标列处理**——它几乎一定是
    写错了名字，让校验要求基准后缀，作者必须显式做决定。
    """
    if name.endswith(NON_COORDINATE_SUFFIXES):
        return False
    parts = name.split("_")
    return any(part in COORDINATE_TOKENS for part in parts)


def columns_of(table: str) -> tuple[str, ...]:
    """取某张表的规范列序。"""
    try:
        return COLUMNS[table]
    except KeyError:
        raise ValueError(f"未知的表：{table}（已注册：{', '.join(sorted(COLUMNS))}）") from None


def validate_column_names(table: str, names: Sequence[str]) -> None:
    """检查一组列名是否满足基准后缀规则；违反即抛 ``AssertionError``。"""
    problems = [
        name for name in names if is_coordinate_column(name) and not name.endswith(BASIS_SUFFIXES)
    ]
    if problems:
        raise AssertionError(
            f"表 {table} 的坐标列缺少基准后缀 {BASIS_SUFFIXES}：{', '.join(problems)}"
        )


def validate_columns(table: str) -> None:
    """检查已注册表的列名是否满足基准后缀规则。"""
    validate_column_names(table, columns_of(table))


def read_tsv(path: str | os.PathLike[str]) -> tuple[list[str], list[dict[str, str]]]:
    """读 TSV：跳过 ``#`` 注释行，返回 (列名, 行字典)。"""
    lines = []
    with open(path, encoding="utf-8", newline="") as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            lines.append(line)
    if not lines:
        return [], []
    reader = csv.DictReader(lines, delimiter="\t", restkey="__extra__")
    return list(reader.fieldnames or ()), [dict(row) for row in reader]


def read_legacy_tsv(
    path: str | os.PathLike[str],
) -> tuple[list[str], list[list[str]]]:
    """读历史格式 TSV：跳过 ``#`` 注释与空行，返回 ``(表头, 行列表)``。

    与 :func:`read_tsv` 的区别是行仍是字符串列表（旧脚本按位置取列）。
    """
    header: list[str] = []
    rows: list[list[str]] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not header:
                header = fields
            else:
                rows.append(fields)
    return header, rows


def _default_file_mode() -> int:
    """按当前 umask 算出默认文件权限（mkstemp 固定 0600，写出后要改回来）。"""
    current = os.umask(0)
    os.umask(current)
    return 0o666 & ~current


class _AtomicWriter:
    """临时文件 + ``os.replace``；失败清理，成功保留目标原权限。"""

    def __init__(self, target: str) -> None:
        self.target = target
        self.temporary = ""

    def __enter__(self) -> TextIO:
        import tempfile

        directory = os.path.dirname(os.path.abspath(self.target))
        os.makedirs(directory, exist_ok=True)
        descriptor, self.temporary = tempfile.mkstemp(
            dir=directory, prefix=".alnkit-", suffix=".tmp"
        )
        self._handle = os.fdopen(descriptor, "w", encoding="utf-8", newline="")
        return self._handle

    def __exit__(
        self,
        exc_type: type[BaseException] | None,
        exc: BaseException | None,
        traceback: TracebackType | None,
    ) -> Literal[False]:
        self._handle.close()
        try:
            if exc_type is None:
                if os.path.exists(self.target):
                    shutil.copymode(self.target, self.temporary)
                else:
                    os.chmod(self.temporary, _default_file_mode())
                os.replace(self.temporary, self.target)
                return False
        finally:
            if os.path.exists(self.temporary):
                os.unlink(self.temporary)
        return False


def _format_value(value: Any) -> str:
    if value is None:
        return "-"
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


def write_tsv(
    path: str | os.PathLike[str],
    table: str,
    rows: Iterable[Mapping[str, Any]],
    *,
    comments: Sequence[str] = (),
    sort_key: Callable[[Mapping[str, Any]], Any] | None = None,
) -> int:
    """按注册列序写 TSV（原子写）；返回写入行数。

    ``sort_key`` 必须显式给出（或确认输入已排序），禁止依赖 dict/set 迭代顺序。
    未注册的列在任何字节落盘前就会报错；失败时不破坏已存在的目标文件。
    """
    columns = columns_of(table)
    column_set = set(columns)
    ordered = [dict(row) for row in rows]
    if sort_key is not None:
        ordered.sort(key=sort_key)
    unknown: set[str] = set()
    for row in ordered:
        unknown |= set(row) - column_set
    if unknown:
        raise ValueError(f"表 {table} 出现未注册的列：{', '.join(sorted(unknown))}")

    target = os.fspath(path)
    with _AtomicWriter(target) as handle:
        for comment in comments:
            handle.write(f"# {comment}\n")
        handle.write("\t".join(columns) + "\n")
        for row in ordered:
            handle.write("\t".join(_format_value(row.get(name)) for name in columns) + "\n")
    return len(ordered)


def write_legacy_tsv(
    path: str | os.PathLike[str],
    columns: Sequence[str],
    rows: Iterable[Sequence[Any]],
    comments: Sequence[str] = (),
) -> int:
    """写历史格式 TSV（原子写）：列名由调用方给、不加基准校验、值用 ``str()``。

    这是三代旧脚本（``sequence_io.write_tsv``）的语义：注释行 + 表头 + 逐行 ``str(value)``，
    不做浮点格式化、不做列注册校验。新表请用 :func:`write_tsv`。
    """
    target = os.fspath(path)
    written = 0
    with _AtomicWriter(target) as handle:
        for comment in comments:
            handle.write("# " + comment + "\n")
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(value) for value in row) + "\n")
            written += 1
    return written


for _table in COLUMNS:
    validate_columns(_table)
