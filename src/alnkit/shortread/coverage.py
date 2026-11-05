"""``short coverage`` 的数据层：BAM + 区间 → 逐碱基深度 → 分箱广度/深度 → 两张表。

数据流（全是数组运算，不含绘图与命令行逻辑）::

    BAM + Region + Filters
        → depth[]                每个碱基上有多少条 read 的比对碱基压在它上面
        → CoverageResult         区域汇总字段 + 逐碱基数组 + 分箱边界
        → coverage_region 表      一行：reads / mapQ / 覆盖广度 / 平均深度
        → coverage_bins 表        每行一个 bin：covered_bases、depth_sum

图的画法在 :mod:`alnkit.plot.coverage`，它只读这两张表。

口径（与 ``samtools coverage`` 对齐）：

* 只有 ``M`` / ``=`` / ``X`` 给深度加 1；``D``/``N`` 占参考位置但没有碱基，不计数；
  ``S``/``I``/``H`` 不占参考位置，不计数；
* ``depth >= min_depth`` 的位置才算“被覆盖”；
* 平均深度的分母是**区间长度**，不是覆盖碱基数。
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np

from ..core import coords, tables
from ..core.errors import UsageError

__all__ = [
    "MAX_REGION_LENGTH",
    "CoverageResult",
    "Filters",
    "Region",
    "accumulate_depth_from_reads",
    "compute_coverage",
    "make_bin_edges",
    "parse_flags",
    "write_coverage_tables",
]

### CIGAR 操作码：决定一个 read 的哪一段落在参考序列上 ###
# 数值来自 SAM 规范，BAM 里存的就是这几个 op code。
CIGAR_MATCH = 0  # M  比对（有碱基、占参考位置）
CIGAR_INSERTION = 1  # I  插入（有碱基、不占参考位置）
CIGAR_DELETION = 2  # D  缺失（无碱基、占参考位置）
CIGAR_REFERENCE_SKIP = 3  # N  跳跃（无碱基、占参考位置）
CIGAR_SOFT_CLIP = 4  # S  软裁剪（有碱基、不占参考位置）
CIGAR_HARD_CLIP = 5  # H  硬裁剪（序列里根本没有）
CIGAR_PADDING = 6  # P  占位（无意义）
CIGAR_SEQUENCE_MATCH = 7  # =  与参考完全一致
CIGAR_SEQUENCE_MISMATCH = 8  # X  与参考不一致

#: 一个碱基真正落在参考位置上（要计入深度）的三种 op
CIGAR_OPS_WITH_ALIGNED_BASE = (CIGAR_MATCH, CIGAR_SEQUENCE_MATCH, CIGAR_SEQUENCE_MISMATCH)
#: 占参考位置但没有碱基：缺失与跳跃都不能算覆盖
CIGAR_OPS_SKIPPING_REFERENCE = (CIGAR_DELETION, CIGAR_REFERENCE_SKIP)
#: 只消耗 read 自己的碱基：插入与软裁剪不占参考位置
CIGAR_OPS_CONSUMING_QUERY_ONLY = (CIGAR_INSERTION, CIGAR_SOFT_CLIP)

### BAM flag 位与 read 级过滤 ###
BAM_FLAG_UNMAPPED = 0x4
BAM_FLAG_SECONDARY = 0x100
BAM_FLAG_QUALITY_FAIL = 0x200
BAM_FLAG_DUPLICATE = 0x400
BAM_FLAG_SUPPLEMENTARY = 0x800

#: 默认跳过的 read：没比对上的、次要比对的、QC 不过的、PCR/测序重复的
DEFAULT_EXCLUDED_FLAGS = (
    BAM_FLAG_UNMAPPED | BAM_FLAG_SECONDARY | BAM_FLAG_QUALITY_FAIL | BAM_FLAG_DUPLICATE
)

#: ``--ff`` / ``--rf`` 里可以写名字，也可以直接写数字
FLAG_BITS_BY_NAME = {
    "PAIRED": 0x1,
    "PROPER_PAIR": 0x2,
    "UNMAP": BAM_FLAG_UNMAPPED,
    "MUNMAP": 0x8,
    "REVERSE": 0x10,
    "MREVERSE": 0x20,
    "READ1": 0x40,
    "READ2": 0x80,
    "SECONDARY": BAM_FLAG_SECONDARY,
    "QCFAIL": BAM_FLAG_QUALITY_FAIL,
    "DUP": BAM_FLAG_DUPLICATE,
    "SUPPLEMENTARY": BAM_FLAG_SUPPLEMENTARY,
}

#: 逐碱基三条数组约 12 byte/碱基；超过这个长度就先报错让用户缩小 --region
MAX_REGION_LENGTH = 20_000_000


def parse_flags(flag_text: str) -> int:
    """把 ``--ff`` / ``--rf`` 的写法转成 flag 位掩码。

    输入示例: parse_flags("UNMAP,DUP") → 输出示例: 1028（= 0x4 | 0x400）
    输入示例: parse_flags("0x504")    → 输出示例: 1284
    输入示例: parse_flags("")         → 输出示例: 0（空串表示不过滤）
    """
    if flag_text.strip() == "":
        return 0
    try:
        return int(flag_text, 0)  # 支持 '0x504' / '1284'
    except ValueError:
        pass
    flag_bits = 0
    for flag_name in flag_text.replace(",", " ").split():
        name = flag_name.upper()
        if name not in FLAG_BITS_BY_NAME:
            raise UsageError(f"未知的 flag 名：{flag_name!r}（可写数字或用名字，如 UNMAP,DUP）")
        flag_bits |= FLAG_BITS_BY_NAME[name]
    return flag_bits


### 过滤条件：哪些 read 参与统计、哪些位置算“被覆盖” ###


@dataclass(frozen=True, slots=True)
class Filters:
    """一次统计用到的全部阈值。

    示例::

        Filters()                                   # 默认：跳过 UNMAP/SECONDARY/QCFAIL/DUP
        Filters(min_mapq=20, min_depth=5)           # 只算 mapQ>=20，深度 5 以上才算覆盖
        Filters(fail_flags=0, max_depth=0)          # 不按 flag 过滤，深度不封顶
    """

    fail_flags: int = DEFAULT_EXCLUDED_FLAGS  # 带上这些 flag 的 read 跳过
    require_flags: int = 0  # 缺这些 flag 的 read 跳过
    min_mapq: int = 0  # 最小比对质量
    min_len: int = 0  # 最小 query 长度（含 soft clip）
    min_depth: int = 1  # 低于该深度的位置不算覆盖
    max_depth: int = 1_000_000  # 每个位置深度上限，0 = 不限


### 区间：内部一律 0-based 半开 ###


@dataclass(frozen=True, slots=True)
class Region:
    """要统计的区间（0-based 半开）。

    示例: Region("chr1", 0, 100) → 覆盖 chr1 的 0..99 位
          Region("chr1", 99, 200) → 命令行写 chr1:100-200
    """

    name: str
    start: int  # 0-based，含
    end: int  # 0-based，不含

    @property
    def length(self) -> int:
        """区间长度。输入示例: Region("chr1", 99, 200) → 输出示例: 101"""
        return self.end - self.start


### 统计结果：区域汇总字段 + 逐碱基数组 + 分箱边界 ###


@dataclass
class CoverageResult:
    """一次统计的全部结果，写表与画图都只读这个对象。

    示例（区间 100 bp，3 条 read 覆盖 0..74）::

        result.covered_bases        -> 75
        result.percent_covered      -> 75.0
        result.mean_coverage        -> 1.25     # Σdepth / 区间长度
        result.breadth_percent_per_bin -> array([100., 100., 100.,   0.])
        result.mean_depth_per_bin      -> array([ 2.,   2.,   1.,   0.])
    """

    region: Region
    filters: Filters
    read_count: int = 0  # fetch 到的全部记录（含被过滤掉的）
    selected_read_count: int = 0  # 通过 flag/mapQ/长度过滤的记录
    mapping_quality_sum: int = 0  # 通过过滤的记录的 mapQ 之和
    depth: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))
    base_quality_sum: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))
    base_quality_base_count: np.ndarray = field(default_factory=lambda: np.zeros(0, np.int32))
    bin_edges: np.ndarray = field(default_factory=lambda: np.zeros(1, np.int64))

    ### 区域汇总 ###
    @property
    def covered_mask(self) -> np.ndarray:
        """逐碱基布尔数组：该位置是否算“被覆盖”。

        输入示例: depth=[2, 2, 0, 1], min_depth=1  →  输出示例: [True, True, False, True]
        """
        return self.depth >= max(self.filters.min_depth, 1)

    @property
    def covered_bases(self) -> int:
        """被覆盖的碱基数。输入示例: 100 bp 区间里 75 个位置有覆盖 → 输出示例: 75"""
        return int(self.covered_mask.sum(dtype=np.int64))

    @property
    def percent_covered(self) -> float:
        """覆盖广度百分比。输入示例: 75 / 100 → 输出示例: 75.0"""
        if not self.region.length:
            return 0.0
        return 100.0 * self.covered_bases / self.region.length

    @property
    def mean_coverage(self) -> float:
        """平均深度（分母是区间长度，不是覆盖碱基数）。

        输入示例: 50 个位置深度 2、25 个位置深度 1、其余 0  →  输出示例: 1.25
        """
        if not self.region.length:
            return 0.0
        return float(self.depth[self.covered_mask].sum(dtype=np.int64)) / self.region.length

    @property
    def mean_baseq(self) -> float:
        """覆盖位置上的平均碱基质量。输入示例: 全部 Q35 → 输出示例: 35.0"""
        counted_bases = int(self.base_quality_base_count[self.covered_mask].sum(dtype=np.int64))
        if not counted_bases:
            return 0.0
        return float(self.base_quality_sum[self.covered_mask].sum(dtype=np.int64)) / counted_bases

    @property
    def mean_mapq(self) -> float:
        """平均比对质量 = mapQ 之和 / 通过过滤的 read 数。输入示例: 180 / 3 → 输出示例: 60.0"""
        if not self.selected_read_count:
            return 0.0
        return self.mapping_quality_sum / self.selected_read_count

    @property
    def filtered_read_count(self) -> int:
        """被过滤掉的 read 数。输入示例: fetch 7 条、选中 5 条 → 输出示例: 2"""
        return self.read_count - self.selected_read_count

    ### 分箱 ###
    @property
    def bin_widths(self) -> np.ndarray:
        """每个 bin 的宽度（bp）。输入示例: edges [0, 25, 50, 100] → 输出示例: [25, 25, 50]"""
        return np.diff(self.bin_edges).astype(np.int64)

    @property
    def covered_bases_per_bin(self) -> np.ndarray:
        """每个 bin 内被覆盖的碱基数。输入示例: → 输出示例: [25, 20, 0]"""
        # 显式转成 int64：numpy < 2 的 reduceat 不保证把 bool 上推成整数
        return np.add.reduceat(self.covered_mask.astype(np.int64), self.bin_edges[:-1])

    @property
    def depth_sum_per_bin(self) -> np.ndarray:
        """每个 bin 内深度之和。输入示例: → 输出示例: [50, 20, 0]"""
        return np.add.reduceat(self.depth.astype(np.int64), self.bin_edges[:-1])

    @property
    def breadth_percent_per_bin(self) -> np.ndarray:
        """上图：每个 bin 的覆盖广度 %。输入示例: bin 宽 25、20 bp 有覆盖 → 输出示例: 80.0"""
        return 100.0 * self.covered_bases_per_bin / self.bin_widths

    @property
    def mean_depth_per_bin(self) -> np.ndarray:
        """下图：每个 bin 的平均深度。输入示例: bin 宽 25、深度和 50 → 输出示例: 2.0"""
        return self.depth_sum_per_bin / self.bin_widths

    @property
    def bin_starts(self) -> np.ndarray:
        """每个 bin 的起始坐标（0-based）。输入示例: region.start=99, edges [0,5] → [99]"""
        return self.region.start + self.bin_edges[:-1]

    @property
    def bin_ends(self) -> np.ndarray:
        """每个 bin 的结束坐标（0-based 不含）。输入示例: region.start=99, edges [0,5] → [104]"""
        return self.region.start + self.bin_edges[1:]


### 分箱 ###


def make_bin_edges(region_length: int, requested_bins: int) -> np.ndarray:
    """把区间等宽切分，返回 bin 边界（相对区间起点，0-based）。

    输入示例: region_length=100, requested_bins=4  → 输出示例: array([0, 25, 50, 75, 100])
    输入示例: region_length=103, requested_bins=4  → 输出示例: array([0, 25, 50, 75, 103])
    输入示例: region_length=5,   requested_bins=50 → 输出示例: array([0, 1, 2, 3, 4, 5])

    说明: bin 宽度 = region_length // bins，前 bins-1 个宽度一致，最后一个吸收余数，
    保证整幅图覆盖整个区间（len(edges) == min(bins, region_length) + 1）。
    """
    bin_count = max(1, min(requested_bins, region_length))
    bin_width = region_length // bin_count
    edges = np.arange(bin_count, dtype=np.int64) * bin_width
    return np.append(edges, region_length)


### 深度累加：扫一遍 reads，按 CIGAR 的比对块往数组里加 ###


def accumulate_depth_from_reads(
    read_iterator: Iterable[Any],
    region: Region,
    filters: Filters,
    depth: np.ndarray,
    base_quality_sum: np.ndarray,
    base_quality_base_count: np.ndarray,
) -> tuple[int, int, int]:
    """把 reads 的比对块累加进三个逐碱基数组（原地修改），并统计 read 计数。

    输入示例::
        read_iterator            = bam.fetch("chr1", 0, 100)
        region                   = Region("chr1", 0, 100)
        depth / base_quality_*   = np.zeros(100, np.int32) 三条

    输出示例（一条 read：CIGAR 5M5D5M、起点 0-based 20、mapQ 60）::

        (1, 1, 60)                 # (read_count, selected_read_count, mapQ 之和)
        depth[20:25] == 1          # 第一段 M
        depth[25:30] == 0          # D 不带碱基 -> 不算覆盖
        depth[30:35] == 1          # 第二段 M

    边界处理: 比对块超出 region 的部分被裁掉，不会写到数组之外。
    """
    region_start, region_end = region.start, region.end
    read_count = selected_read_count = mapping_quality_sum = 0

    for read in read_iterator:
        read_count += 1

        # 1) read 级过滤：flag -> mapQ -> query 长度
        if read.flag & filters.fail_flags:
            continue
        if filters.require_flags and not read.flag & filters.require_flags:
            continue
        if read.mapping_quality < filters.min_mapq:
            continue
        cigar_blocks = read.cigartuples or ()  # 没有 CIGAR 的 read 不贡献碱基，
        if not cigar_blocks:  # 但它仍算一条被选中的 read
            selected_read_count += 1
            mapping_quality_sum += read.mapping_quality
            continue
        if filters.min_len:
            query_length = sum(
                block_length
                for cigar_operation, block_length in cigar_blocks
                if cigar_operation in CIGAR_OPS_CONSUMING_QUERY_ONLY
                or cigar_operation in CIGAR_OPS_WITH_ALIGNED_BASE
            )
            if query_length < filters.min_len:
                continue

        selected_read_count += 1
        mapping_quality_sum += read.mapping_quality

        # 2) 按 CIGAR 走 read：只有“有碱基且占参考位置”的块才计入深度
        base_qualities = None if read.query_qualities is None else np.asarray(read.query_qualities)
        reference_position = read.reference_start
        query_position = 0
        for cigar_operation, block_length in cigar_blocks:
            if cigar_operation in CIGAR_OPS_WITH_ALIGNED_BASE:
                # 与 region 求交，避免越界写入
                overlap_start = max(reference_position, region_start)
                overlap_end = min(reference_position + block_length, region_end)
                if overlap_start < overlap_end:
                    array_start = overlap_start - region_start
                    array_end = overlap_end - region_start
                    depth[array_start:array_end] += 1
                    if base_qualities is not None:
                        query_start = query_position + overlap_start - reference_position
                        query_end = query_position + overlap_end - reference_position
                        base_quality_sum[array_start:array_end] += base_qualities[
                            query_start:query_end
                        ]
                        base_quality_base_count[array_start:array_end] += 1
                reference_position += block_length
                query_position += block_length
            elif cigar_operation in CIGAR_OPS_SKIPPING_REFERENCE:
                reference_position += block_length
            elif cigar_operation in CIGAR_OPS_CONSUMING_QUERY_ONLY:
                query_position += block_length
            # CIGAR_HARD_CLIP / CIGAR_PADDING 两边都不推进

    return read_count, selected_read_count, mapping_quality_sum


### 对外入口 ###


def compute_coverage(
    bam, region: Region, filters: Filters | None = None, requested_bins: int = 200
) -> CoverageResult:
    """统计一个区间的覆盖度，返回 :class:`CoverageResult`。

    输入示例::
        bam      = pysam.AlignmentFile("reads.bam")
        region   = Region("chr1", 99, 200)          # 命令行写 chr1:100-200
        filters  = Filters(min_mapq=20)
        requested_bins = 4

    输出示例::
        result.selected_read_count     -> 12
        result.covered_bases           -> 75
        result.mean_coverage           -> 1.25
        result.breadth_percent_per_bin -> array([100., 100., 100., 0.])
        result.bin_edges               -> array([0, 25, 50, 75, 101])

    区间过长会抛 ``UsageError``（逐碱基数组要 12 byte/碱基，一次几 GB 内存不值得）。
    """
    region_length = region.length
    if filters is None:
        filters = Filters()
    if region_length > MAX_REGION_LENGTH:
        raise UsageError(
            f"区间 {region.name}:{coords.to_1based(region.start)}-{region.end} 长 "
            f"{region_length / 1e6:.0f} Mb；逐碱基数组需要 {region_length * 12 / 1e9:.1f} GB。"
            "请缩小 --region"
        )

    depth = np.zeros(region_length, np.int32)
    base_quality_sum = np.zeros(region_length, np.int32)
    base_quality_base_count = np.zeros(region_length, np.int32)

    read_iterator = bam.fetch(region.name, region.start, region.end)
    read_count, selected_read_count, mapping_quality_sum = accumulate_depth_from_reads(
        read_iterator, region, filters, depth, base_quality_sum, base_quality_base_count
    )

    if filters.max_depth > 0:  # --max-depth 对每个位置封顶
        np.minimum(depth, filters.max_depth, out=depth)

    return CoverageResult(
        region=region,
        filters=filters,
        read_count=read_count,
        selected_read_count=selected_read_count,
        mapping_quality_sum=mapping_quality_sum,
        depth=depth,
        base_quality_sum=base_quality_sum,
        base_quality_base_count=base_quality_base_count,
        bin_edges=make_bin_edges(region_length, requested_bins),
    )


### 写表 ###


def region_row(result: CoverageResult, sample: str) -> dict[str, Any]:
    """把区域汇总整理成 ``coverage_region`` 的一行。

    输入示例: 区间 100 bp、覆盖 75 bp、3 条 read、平均深度 1.25
    输出示例: {"sample": "sampleA", "chrom": "chr1", "region_start_1based": 1,
              "n_reads_selected": 3, "percent_covered": 75.0, "mean_depth": 1.25, …}
    """
    # 0-based 半开 [start, end) -> 1-based 闭区间：区间列一律由 coords 换算
    region_start_1based, region_end_1based = coords.interval_0based_to_1based(
        result.region.start, result.region.end
    )
    return {
        "sample": sample,
        "chrom": result.region.name,
        "region_start_1based": region_start_1based,
        "region_end_1based": region_end_1based,
        "region_length_bp": result.region.length,
        "n_bins": len(result.bin_edges) - 1,
        "n_reads_total": result.read_count,
        "n_reads_selected": result.selected_read_count,
        "n_reads_filtered": result.filtered_read_count,
        "mean_mapq": result.mean_mapq,
        "mean_baseq": result.mean_baseq,
        "covered_bases": result.covered_bases,
        "percent_covered": result.percent_covered,
        "mean_depth": result.mean_coverage,
    }


def bin_rows(result: CoverageResult) -> list[dict[str, Any]]:
    """把分箱结果整理成 ``coverage_bins`` 的行（图里要的广度/深度都能从这 6 列算出来）。

    输入示例: 4 个 bin，覆盖 [25, 20, 0, 30]，深度和 [50, 20, 0, 60]
    输出示例: [{"bin_index": 1, "bin_start_1based": 1, "bin_end_1based": 25,
               "bin_width_bp": 25, "covered_bases": 25, "depth_sum": 50}, …]
    """
    return [
        {
            "bin_index": index,
            "bin_start_1based": coords.to_1based(int(bin_start)),
            "bin_end_1based": int(bin_end),
            "bin_width_bp": int(bin_width),
            "covered_bases": int(covered_bases),
            "depth_sum": int(depth_sum),
        }
        for index, (bin_start, bin_end, bin_width, covered_bases, depth_sum) in enumerate(
            zip(
                result.bin_starts,
                result.bin_ends,
                result.bin_widths,
                result.covered_bases_per_bin,
                result.depth_sum_per_bin,
            ),
            start=1,
        )
    ]


def base_rows(result: CoverageResult) -> Iterator[dict[str, Any]]:
    """逐碱基深度（``--bases-table`` 用）：``pos_1based`` + ``depth``。

    输入示例: 区间 chr1:1-3，深度 [2, 0, 1]
    输出示例: [{"pos_1based": 1, "depth": 2}, {"pos_1based": 2, "depth": 0}, …]
    """
    for offset, depth_value in enumerate(result.depth):
        yield {
            "pos_1based": coords.to_1based(result.region.start + offset),
            "depth": int(depth_value),
        }


def write_coverage_tables(
    result: CoverageResult,
    sample: str,
    *,
    region_path: str | Path,
    bins_path: str | Path,
    bases_path: str | Path | None = None,
    comments: Sequence[str] = (),
) -> dict[str, int]:
    """写 ``coverage_region`` / ``coverage_bins``（以及可选的 ``coverage_bases``）。

    输入示例: region_path="out/demo.region.tsv", bins_path="out/demo.bins.tsv"
    输出示例: {"region": 1, "bins": 4}（写出的行数）
    """
    written = {
        "region": tables.write_tsv(
            region_path, "coverage_region", [region_row(result, sample)], comments=comments
        ),
        "bins": tables.write_tsv(
            bins_path,
            "coverage_bins",
            bin_rows(result),
            comments=comments,
            sort_key=lambda row: row["bin_index"],
        ),
    }
    if bases_path is not None:
        written["bases"] = tables.write_tsv(
            bases_path,
            "coverage_bases",
            base_rows(result),
            comments=comments,
            sort_key=lambda row: row["pos_1based"],
        )
    return written
