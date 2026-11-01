"""BAM / CIGAR 层：解析比对块、取断点坐标、抽取插入命中、统计比对特征。"""

from __future__ import annotations

import functools
import logging
import re
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import pysam

from .model import (
    AlignmentBlock,
    AlignmentRecord,
    CIGAR_FULL_RE,
    CIGAR_TOKEN_RE,
    ConfigurationError,
    InsertHit,
    LEFT_CLIP_RE,
    RIGHT_CLIP_RE,
)


# ------------------------------------------------------------------------------
# ### 4. CIGAR 解析与断点坐标计算 ###
# ------------------------------------------------------------------------------
def parse_alignment_blocks(alignment: "pysam.AlignedSegment") -> List[AlignmentBlock]:
    """把一条比对的 CIGAR 拆成若干"参考序列上的比对块"（只统计 M/=/X）。

    输入示例:
        alignment.cigarstring = '47S67M1I35M', alignment.reference_start = 7000000
    输出示例:
        [AlignmentBlock(reference_start=7000000, reference_end=7000100, read_start=47, read_end=114, length=67),
         AlignmentBlock(reference_start=7000101, reference_end=7000200, read_start=115, read_end=150, length=35)]

    说明（与历史实现保持一致，便于复现）:
        - reference_start/end 为 0-based 半开区间；read_start/end 为 read 内部偏移
        - soft clip（S）把 read 游标直接置成 clip 长度；insertion（I）与 hard clip（H）累加到 read 游标
        - deletion（D）/skip（N）/padding（P）只累加参考游标
        - CIGAR 拆块结果按“几何量”缓存（同一 BAM 里 CIGAR 种类很少，实测 2.6 万条比对只有 856 种）；
          **只缓存相对偏移**，绝对坐标必须再加 reference_start，否则不同位点的比对会串坐标
    """
    reference_start = alignment.reference_start
    return [
        AlignmentBlock(
            reference_start + block_start,
            reference_start + block_end,
            read_start,
            read_end,
            block_end - block_start,
        )
        for block_start, block_end, read_start, read_end in cached_block_geometry(
            tuple(alignment.cigartuples or ())
        )
    ]


@functools.lru_cache(maxsize=8192)
def cached_block_geometry(
    cigartuples: Tuple[Tuple[int, int], ...],
) -> Tuple[Tuple[int, int, int, int], ...]:
    """把 CIGAR 拆成相对几何量：(参考起点偏移, 参考终点偏移, read 起点, read 终点)。

    输入示例: ((4, 47), (0, 67), (1, 1), (0, 35))        # 47S67M1I35M
    输出示例: ((0, 67, 47, 114), (68, 103, 115, 150))

    说明: 只依赖 CIGAR，因此可以跨比对复用；绝对坐标由调用方加 reference_start。
    """
    geometry: List[Tuple[int, int, int, int]] = []
    read_cursor = 0
    reference_cursor = 0
    for operation, length in cigartuples:
        if operation == 4:  # S：soft clip
            read_cursor = length
            continue
        if operation in (0, 7, 8):  # M / = / X：比对块
            geometry.append(
                (reference_cursor, reference_cursor + length, read_cursor, read_cursor + length)
            )
            read_cursor += length
            reference_cursor += length
            continue
        if operation in (1, 5):  # I / H：消耗 read
            read_cursor += length
            continue
        if operation in (2, 3, 6):  # D / N / P：消耗参考
            reference_cursor += length
            continue
    return tuple(geometry)


def detect_clip_side(cigarstring: str) -> Tuple[str, int]:
    """判断 soft clip 出现在比对的哪一侧，并返回 clip 长度。

    输入示例: '47S67M1I35M' -> ('left', 47)
              '99M51S'      -> ('right', 51)
              '150M'        -> ('none', 0)
    输出示例: (clip_side, clip_length)，clip_side 取 'left' / 'right' / 'none'
    """
    left_match = LEFT_CLIP_RE.match(cigarstring)
    if left_match:
        return "left", int(left_match.group(1))
    right_match = RIGHT_CLIP_RE.search(cigarstring)
    if right_match:
        return "right", int(right_match.group(2))
    return "none", 0


def choose_evidence_block(
    blocks: List[AlignmentBlock], cigarstring: str, legacy_coordinates: bool
) -> AlignmentBlock:
    """选取用于计算断点的那个比对块。

    输入示例:
        blocks = [67M 块(7000000-7000100), 35M 块(7000101-7000200)]
        cigarstring = '47S67M1I35M'
    输出示例:
        legacy_coordinates=False -> 67M 块（与 clip 相邻，断点=7000000）
        legacy_coordinates=True  -> 35M 块（历史行为，断点=7000200，实测会偏 67 bp）

    说明: 历史实现固定取"最后一个 M 块"，遇到 '47S67M1I35M' 这类中间带 I/D 的 CIGAR
          会把断点算到比对内部；默认改为取与 clip 相邻的 M 块。
    """
    if legacy_coordinates:
        return blocks[-1]
    clip_side, _ = detect_clip_side(cigarstring)
    if clip_side == "left":
        return blocks[0]
    return blocks[-1]


def compute_junction_position(record: AlignmentRecord) -> int:
    """由比对记录算出断点位置（参考序列上的坐标）。

    输入示例: AlignmentRecord(cigarstring='99M51S',  block_end=4000100, clip_side='right')
    输出示例: 4000100
              AlignmentRecord(cigarstring='47S67M1I35M', block_start=7000000, clip_side='left') -> 7000000

    说明: 左侧 clip -> 取比对块起始；其余 -> 取比对块末端。
          两种写法都等于"断点左侧那个碱基的 1-based 坐标"，因此可以互相比较。
    """
    if record.clip_side == "left":
        return record.block_start
    return record.block_end


def build_alignment_record(
    alignment: "pysam.AlignedSegment", block: AlignmentBlock, mismatch_count: int
) -> AlignmentRecord:
    """把 pysam 比对对象 + 选定比对块，组装成 AlignmentRecord。

    输入示例:
        alignment: cigarstring='99M51S', is_reverse=True, query_length=150,
                   reference_name='chr1', mapping_quality=60
        block: AlignmentBlock(reference_start=4000000, reference_end=4000100,
                              read_start=0, read_end=99, length=99)
    输出示例:
        AlignmentRecord(read_start=51, read_end=150, strand='-', clip_side='right', clip_length=51, ...)
        （反向链的 read 偏移按历史实现翻转为原始 read 方向：start' = 长度 - end，end' = 长度 - start）

    说明: 输出的 12 个字段顺序与 reads_info 表的数据列一一对应。
          read 长度优先取 `query_length`（SEQ 存在时 = len(seq)）；没有 SEQ 的 BAM 里
          `query_length` 为 0，此时回退到 `infer_query_length()`（由 CIGAR 推算，pysam API），
          避免反向链的 read 偏移被算成负数。
    """
    clip_side, clip_length = detect_clip_side(alignment.cigarstring or "")
    read_start, read_end = block.read_start, block.read_end
    if alignment.is_reverse:
        query_length = alignment.query_length or alignment.infer_query_length() or 0
        read_start, read_end = query_length - block.read_end, query_length - block.read_start
    return AlignmentRecord(
        read_name=alignment.query_name,
        read_pair="read1" if alignment.is_read1 else "read2",
        cigarstring=alignment.cigarstring or "",
        is_supplementary=alignment.is_supplementary,
        chromosome=alignment.reference_name,
        block_start=block.reference_start,
        block_end=block.reference_end,
        read_start=read_start,
        read_end=read_end,
        block_length=block.length,
        strand="-" if alignment.is_reverse else "+",
        mapping_quality=alignment.mapping_quality,
        mismatch_count=mismatch_count,
        clip_side=clip_side,
        clip_length=clip_length,
        is_duplicate=alignment.is_duplicate,
    )


def split_length_pair_matches(cigarstring: str, supplementary_cigarstring: str) -> bool:
    """判断两条 CIGAR 的"M/S 长度对"是否相同（用于同源区分的历史校验）。

    输入示例: cigarstring='117M33S', supplementary_cigarstring='116S34M'
    输出示例: False   （长度对 [117,33] vs [116,34] 不同）
              cigarstring='116S34M', supplementary_cigarstring='116S34M' -> True
    说明: 历史实现在 CIGAR 不匹配 `\\d+[SM]\\d+[MS]` 时会抛异常，这里改为返回 False。
    """
    own_match = re.search(r"(\d+)[SM](\d+)[MS]", cigarstring)
    other_match = re.search(r"(\d+)[SM](\d+)[MS]", supplementary_cigarstring)
    if not own_match or not other_match:
        return False
    own_pair = sorted([own_match.group(1), own_match.group(2)])
    other_pair = sorted([other_match.group(1), other_match.group(2)])
    return own_pair == other_pair


# ------------------------------------------------------------------------------
# ### 5. 过滤规则：错配数、质量与结构 ###
# ------------------------------------------------------------------------------
def count_mismatches(
    alignment: "pysam.AlignedSegment",
    mismatch_source: str,
    on_missing_tag: str,
    statistics: Optional[Dict[str, int]] = None,
) -> Optional[int]:
    """统计一条比对的错配数。

    输入示例:
        alignment.get_tag('NM') = 0                       -> 0
        alignment 无 NM 标签, on_missing_tag='skip'        -> None（调用方丢弃该比对）
        alignment 无 NM 标签, on_missing_tag='error'       -> 抛错退出

    输出示例: 0 / 1 / None
    说明:
        nm: 直接读 NM 标签（= 错配 + 插入 + 删除碱基总数），最稳健
        md: 历史算法 —— 数 MD 标签里的 [ATGC] 字母，若含 '^'（删除）再减去 CIGAR 的删除碱基数
    """
    if mismatch_source == "nm":
        if not alignment.has_tag("NM"):
            if statistics is not None:
                statistics["drop_missing_tag"] = statistics.get("drop_missing_tag", 0) + 1
            if on_missing_tag == "error":
                raise RuntimeError(f"比对缺少 NM 标签: {alignment.query_name}")
            return None
        return int(alignment.get_tag("NM"))

    if not alignment.has_tag("MD"):
        if statistics is not None:
            statistics["drop_missing_tag"] = statistics.get("drop_missing_tag", 0) + 1
        if on_missing_tag == "error":
            raise RuntimeError(f"比对缺少 MD 标签: {alignment.query_name}")
        return None
    deletion_base_count = sum(
        length for operation, length in (alignment.cigartuples or ()) if operation == 2
    )
    md_tag = alignment.get_tag("MD")
    base_letters = re.findall(r"[ATGC]", md_tag)
    if "^" in md_tag:
        return len(base_letters) - deletion_base_count
    return len(base_letters)


def has_matching_split_alignment(alignment: "pysam.AlignedSegment") -> bool:
    """校验比对在目标序列上时，是否存在"等长 split 比对"（SA 标签）。

    输入示例:
        alignment.cigarstring = '117M33S'
        alignment.get_tag('SA') = 'target_ref,6078,-,116S34M,0,0;'
    输出示例:
        False（长度对 [117,33] 与 [116,34] 不同）
        True （存在与本 read CIGAR 等长的 SA 记录）

    说明: 只有主比对落在 --query-chr 上时才调用；没有 SA 标签时返回 True（历史行为：
          此时不做校验）。
    """
    if not alignment.has_tag("SA"):
        return True
    for supplementary in alignment.get_tag("SA").split(";"):
        if not supplementary.strip():
            continue
        fields = supplementary.split(",")
        if len(fields) < 4:
            continue
        if split_length_pair_matches(alignment.cigarstring or "", fields[3]):
            return True
    return False


# ------------------------------------------------------------------------------
# ### 6. 共有 read 的区分：插入片段归属判定 ###
# ------------------------------------------------------------------------------
def score_insert_hit(hit: InsertHit) -> Tuple[int, int, int, int]:
    """对一条“插入序列上的命中”打分，用于在两个插入片段之间裁决 read 归属。

    输入示例:
        InsertHit(chromosome='homolog_ref', mapping_quality=0,  block_length=31, mismatch_count=0, clip_length=119)
        InsertHit(chromosome='target_ref',  mapping_quality=37, block_length=150, mismatch_count=0, clip_length=0)
    输出示例: (0, 31, 0, -119) 与 (37, 150, 0, 0)

    说明（依次比较，越大越可信）: 比对质量 -> 比对块长度 -> 错配数（取负）-> clip 长度（取负）
    """
    return (
        hit.mapping_quality,
        hit.block_length,
        -hit.mismatch_count,
        -hit.clip_length,
    )


def measure_cigar_block_and_clip(cigarstring: str) -> Tuple[int, int]:
    """从 CIGAR（可来自 XA 标签）里量出比对块总长度与 soft clip 总长度。

    输入示例: '116S34M' -> (34, 116)；'81S69M' -> (69, 81)；'150M' -> (150, 0)
    输出示例: (block_length, clip_length)
    """
    block_length = 0
    clip_length = 0
    for length, operation in CIGAR_TOKEN_RE.findall(cigarstring):
        if operation in ("M", "=", "X"):
            block_length += int(length)
        elif operation == "S":
            clip_length += int(length)
    return block_length, clip_length


def parse_xa_entry(entry: str, insert_chromosomes: Sequence[str]) -> Optional[InsertHit]:
    """解析一个 XA 备选比对条目，转成 InsertHit（只保留插入序列上的命中）。

    输入示例（两种常见 XA 写法都支持）:
        'target_ref,+335,31M119S,0'                # 4 字段: rname, pos(带符号), CIGAR, NM
        'target_ref,335,+,31M119S,60,0'           # 6 字段: rname, pos, strand, CIGAR, mapQ, NM
    输出示例: InsertHit(chromosome='target_ref', mapping_quality=0, block_length=31,
                        mismatch_count=0, clip_length=119, source='xa', record=None)
    非插入序列或不含 CIGAR 时输出示例: None

    说明: 字段数不统一，因此按“哪个字段长得像 CIGAR”来定位，再取它后面的数字作 mapQ/NM。
    """
    fields = [field.strip() for field in entry.split(",")]
    if len(fields) < 4 or fields[0] not in insert_chromosomes:
        return None
    cigar_index = None
    for index, field in enumerate(fields):
        if CIGAR_FULL_RE.fullmatch(field):
            cigar_index = index
            break
    if cigar_index is None:
        return None
    block_length, clip_length = measure_cigar_block_and_clip(fields[cigar_index])
    if block_length == 0:
        return None
    trailing_numbers = [int(field) for field in fields[cigar_index + 1 :] if field.isdigit()]
    return InsertHit(
        chromosome=fields[0],
        mapping_quality=trailing_numbers[0] if len(trailing_numbers) >= 2 else 0,
        block_length=block_length,
        mismatch_count=trailing_numbers[-1] if trailing_numbers else 0,
        clip_length=clip_length,
        source="xa",
        record=None,
    )


def collect_insert_hits(
    alignment: "pysam.AlignedSegment",
    insert_chromosomes: Sequence[str],
    own_record: Optional[AlignmentRecord],
) -> List[InsertHit]:
    """收集一条比对带来的“插入序列命中”：自身命中（若落在插入序列上）+ XA 标签里的备选命中。

    输入示例:
        alignment.reference_name='chr1', alignment.mapping_quality=60
        own_record = chr1 上的 AlignmentRecord
        alignment.get_tag('XA') = 'homolog_ref,+335,31M119S,0;decoy_ref,+335,31M119S,0;'
        insert_chromosomes = ['decoy_ref', 'homolog_ref', 'target_ref']
    输出示例:
        [InsertHit('homolog_ref', 0, 31, 0, 119, 'xa', None),
         InsertHit('decoy_ref', 0, 31, 0, 119, 'xa', None)]

    说明: 自身命中带上 record 与 source（alignment / secondary），供归属打分与干扰标记使用；
          XA 命中没有对应比对记录，source='xa'、record=None。
    """
    insert_hits: List[InsertHit] = []
    if alignment.reference_name in insert_chromosomes and own_record is not None:
        insert_hits.append(
            InsertHit(
                chromosome=own_record.chromosome,
                mapping_quality=own_record.mapping_quality,
                block_length=own_record.block_length,
                mismatch_count=own_record.mismatch_count,
                clip_length=own_record.clip_length,
                source="secondary" if alignment.is_secondary else "alignment",
                record=own_record,
            )
        )

    if alignment.has_tag("XA"):
        xa_tag = str(alignment.get_tag("XA"))
        # 子串预检：XA 里没有插入序列名时，逐条解析/正则匹配都是白工（大 BAM 里绝大多数如此）
        if mentions_any_contig(xa_tag, insert_chromosomes):
            for entry in xa_tag.split(";"):
                if not entry.strip():
                    continue
                hit = parse_xa_entry(entry, insert_chromosomes)
                if hit is not None:
                    insert_hits.append(hit)
    return insert_hits


def mentions_any_contig(tag_value: object, contig_names: Sequence[str]) -> bool:
    """标签值（SA/XA 的字符串）里是否出现任一 contig 名。

    输入示例: tag_value='target_ref,+335,31M119S,0;chr1,100,+,50M100S,0;',
              contig_names=('decoy_ref','homolog_ref','target_ref')
    输出示例: True（含 'target_ref'）
              同类输入只含 'chr22' 时 -> False

    说明: 这是纯子串判定（C 级字符串搜索），比“先 split 再逐条正则”快很多；
          可能因 contig 名互为子串而多收，但只会“多看一条 read”，不会误删。
    """
    text = str(tag_value)
    return any(name in text for name in contig_names)


# ------------------------------------------------------------------------------
# ### 7. 重复 read（DUP）：只认 BAM 的 markdup 标记 ###
# ------------------------------------------------------------------------------
def read_is_duplicate(records: Sequence[AlignmentRecord]) -> bool:
    """一条 read 是否被上游 markdup 标记为重复（read 级聚合）。

    输入示例:
        records = [chr1 上的比对（is_duplicate=False）, target_ref 上的比对（is_duplicate=True）]
    输出示例: True

    说明:
        - samtools markdup 默认**只给 primary 记录打 FLAG 0x400**（补充比对要加 -S 才打），
          而一条 read 的两段比对可能分别落在基因组侧与插入侧，所以必须按 read 聚合，
          不能只看某一条记录，也不能因为“有一条没标记”就判为非重复。
        - 3.0 起不再有“按签名聚类”的判据：签名（染色体+断点+比对块）在真断点上几乎人人相同，
          会把独立片段算成重复（实测 49 vs 21 条），无法修正，已删除。
    """
    return any(record.is_duplicate for record in records)


def warn_about_alignment_features(statistics: Dict[str, int], bam_path: Path) -> None:
    """检查这个 BAM 的比对特征是否满足脚本要求（SA / soft clip / 补充比对 / XA），不满足就告警。

    输入示例（本项目 BAM 的实测值）:
        statistics = {'alignments_total': 25969, 'alignments_with_sa': 1798,
                      'alignments_with_xa': 11339, 'supplementary_alignments': 0,
                      'secondary_alignments': 906, 'hard_clipped_alignments': 0}
    输出示例（日志）:
        WARNING 本 BAM 没有 is_supplementary 记录、但有 906 条次要比对：典型是 bwa 用了 -M
                （把补充比对标成 secondary）。去掉 -M 并保留 -Y 重跑，可让脚本把补充比对
                也纳入证据与归属判定。
        INFO    比对特征检查通过：SA=1798 补充比对=0 次要比对=906 XA=11339 hard clip=0

    说明: 这些告警只影响“证据多少”，不影响脚本能否运行；参数要求见 docs/bwa_比对参数说明.md。
    """
    warnings_issued = 0
    if statistics["alignments_total"] == 0:
        logging.warning("%s 里没有可比对记录（全部未比对？）", bam_path)
        return
    if statistics["alignments_with_sa"] == 0:
        logging.warning(
            "%s 里没有任何 SA 标签：嵌合（split）证据可能全部缺失。"
            "请检查比对时是否用了 -M、或比对器不输出 SA（要求见 docs/bwa_比对参数说明.md）",
            bam_path,
        )
        warnings_issued += 1
    if statistics["supplementary_alignments"] == 0 and statistics["secondary_alignments"] > 0:
        logging.warning(
            "%s 没有 is_supplementary 记录、但有 %d 条次要比对：典型是 bwa 用了 -M"
            "（把补充比对标成 secondary）。去掉 -M 并保留 -Y 重跑，可让脚本把补充比对也纳入"
            "证据与归属判定（当前仍可用次要比对/XA 做归属判定）",
            bam_path,
            statistics["secondary_alignments"],
        )
        warnings_issued += 1
    if statistics["hard_clipped_alignments"] > 0:
        logging.warning(
            "%s 有 %d 条 hard clip 比对（补充比对未用 -Y 的典型特征）：这些记录提供不了断点/"
            "归属信息，建议比对时加 -Y",
            bam_path,
            statistics["hard_clipped_alignments"],
        )
        warnings_issued += 1
    if statistics["alignments_with_xa"] == 0:
        logging.warning(
            "%s 里没有任何 XA 标签（bwa mem 默认会输出，除非 -h 0）：同源/共有 read 的判定"
            "只能依赖次要比对记录",
            bam_path,
        )
        warnings_issued += 1
    if warnings_issued == 0:
        logging.info(
            "比对特征检查通过：SA=%d 补充比对=%d 次要比对=%d XA=%d hard clip=%d",
            statistics["alignments_with_sa"],
            statistics["supplementary_alignments"],
            statistics["secondary_alignments"],
            statistics["alignments_with_xa"],
            statistics["hard_clipped_alignments"],
        )


def report_duplicate_flags(statistics: Dict[str, int], bam_path: Path, require_flags: bool) -> None:
    """检查 BAM 是否做过 markdup（FLAG 0x400），未做就告警（或按配置报错）。

    输入示例:
        statistics = {'pass1_alignments': 25969, 'duplicate_alignments': 4373,
                      'optical_duplicate_alignments': 0}
    输出示例（日志）:
        INFO    重复标记: 4373 / 25969 条比对（16.8%），其中光学重复 0 条
        WARNING <bam> 里没有任何 FLAG 0x400（未做 markdup）：dup 两列将为 0。
                建议先 samtools markdup

    说明: 只统计第 1 遍扫到的全部比对（已关闭候选闸门时 = 全文件），用于判断“这个 BAM 有没有
          被 markdup 处理过”；它不影响 SR/PESR，只决定 dup 两列有没有意义。
    """
    duplicate_alignments = statistics["duplicate_alignments"]
    if duplicate_alignments == 0:
        message = (
            "%s 里没有任何 FLAG 0x400（未做 markdup）：dup_read_count / dup_read_name 两列"
            "将为 0。建议先 samtools markdup"
        )
        if require_flags:
            raise ConfigurationError(message % bam_path)
        logging.warning(message, bam_path)
        return
    examined = statistics["pass1_alignments"] or statistics["alignments_total"]
    logging.info(
        "重复标记: %d / %d 条比对（%.1f%%），其中光学重复 %d 条",
        duplicate_alignments,
        examined,
        100.0 * duplicate_alignments / max(examined, 1),
        statistics["optical_duplicate_alignments"],
    )
