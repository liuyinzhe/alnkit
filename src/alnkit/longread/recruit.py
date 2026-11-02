"""方案 C（借鉴 DEGAP）：围绕 contig 末端招募 read，只把“悬垂”移过来接着长。

与贪心延伸的区别：贪心每次只吃一条已经 all-vs-all 比对过的 read；这里每轮只拿 contig
末端 500 bp 当探针做一次 minimap2 招募（不是 all-vs-all），把越过端点的 read 的**悬垂**
按“同一致序列前缀出现次数最多”的簇移植过来。参考 DEGAP 的两条关键做法：

1. **围绕已知端点招募**，而不是让组装器自由发挥（实测 hifiasm 只能落回 3/8 个断点锚点）；
2. **每轮验收**：本轮必须真的延长 ≥ ``min_progress_bp``，否则丢弃（不回滚已写出的 contig，
   因为只有验收通过才修改序列）。

用 minimap2 + pysam，不引入 hifiasm/racon；失败或没有招募到东西就原样返回。
"""

from __future__ import annotations

import collections
import os
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from ..core.fasta import reverse_complement, write_fasta_records
from ..core.logging import get_logger

__all__ = ["RecruitSettings", "extend_contig_ends", "parse_alignment"]

logger = get_logger("assemble.recruit")

#: 招募到多长的悬垂才值得移植（也是“本轮有进展”的判据）
MIN_PROGRESS_BP = 100
#: 判断两个悬垂是否同源时比较的前缀长度
CLUSTER_PREFIX_BP = 50
#: 认定“比对到达端点”的容差
EDGE_TOLERANCE_BP = 5
#: 判断“悬垂里装的是重复序列”时用的 k-mer 长度
REPEAT_KMER_BP = 50
#: 悬垂前缀出现条数超过“接缝锚点出现条数”的这个倍数就算重复
REPEAT_DEPTH_FACTOR = 3.0


@dataclass(slots=True, frozen=True)
class RecruitSettings:
    """一轮招募的参数（默认值来自 DEGAP README 的量级）。"""

    edge_bp: int = 500
    min_identity: float = 0.99
    min_mapq: int = 20
    min_support: int = 2
    min_progress_bp: int = MIN_PROGRESS_BP
    max_rounds: int = 6
    threads: int = 4


@dataclass(slots=True, frozen=True)
class Alignment:
    """一条 read 打到探针上的比对（坐标都在探针坐标系里）。"""

    read_name: str
    reverse: bool
    reference_start: int
    reference_end: int
    query_start: int
    query_end: int
    identity: float
    mapq: int


_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")


def _walk_cigar(cigar: str, stop_reference: int):
    """沿 CIGAR 走到参考坐标 ``stop_reference``，返回对应的 read 偏移。

    输入示例: cigar="100S400M200M", stop_reference=400
    输出示例: 500（前 100 个 S 也占 read 偏移，但不占参考）
    """
    reference = 0
    query = 0
    for length, operation in _CIGAR.findall(cigar):
        length = int(length)
        if operation in ("M", "=", "X"):
            if reference + length >= stop_reference:
                return query + (stop_reference - reference)
            reference += length
            query += length
        elif operation in ("I", "S"):
            query += length
        elif operation in ("D", "N"):
            if reference + length >= stop_reference:
                return query
            reference += length
        # H/P 不占 read 也不占参考
    return query


def parse_alignment(line: str, edge_length: int, min_mapq: int, min_identity: float):
    """把一行 SAM 解析成 :class:`Alignment`；不合格（未比对/secondary/低质量）返回 None。

    identity 用 ``NM`` 标签算（``1 - NM / 比对长度``），这是 minimap2 自己的口径。
    """
    if line.startswith("@"):
        return None
    fields = line.rstrip("\n").split("\t")
    if len(fields) < 11:
        return None
    flag = int(fields[1])
    if flag & 0x4 or flag & 0x100 or flag & 0x800:
        return None
    mapq = int(fields[4])
    if mapq < min_mapq:
        return None
    cigar = fields[5]
    reference = 0
    for length, operation in _CIGAR.findall(cigar):
        if operation in ("M", "D", "N", "=", "X"):
            reference += int(length)
    if reference <= 0:
        return None
    matches = next((int(tag[5:]) for tag in fields[11:] if tag.startswith("NM:i:")), None)
    identity = 1.0 - (matches or 0) / reference
    if identity < min_identity:
        return None
    start = int(fields[3]) - 1
    start = max(0, min(start, edge_length - 1))
    end = min(edge_length, start + reference)
    query_end = _walk_cigar(cigar, end - start)
    return Alignment(
        read_name=fields[0],
        reverse=bool(flag & 0x10),
        reference_start=start,
        reference_end=end,
        query_start=_walk_cigar(cigar, 0),
        query_end=query_end,
        identity=identity,
        mapq=mapq,
    )


def _count_reads_containing(probe: str, read_sequences: dict[str, str]) -> int:
    """有多少条 read 含这个 k-mer（正反链都算）。"""
    complement = reverse_complement(probe)
    return sum(
        1 for sequence in read_sequences.values() if probe in sequence or complement in sequence
    )


def _looks_repetitive(hangover: str, anchor: str, read_sequences: dict[str, str]) -> bool:
    """悬垂是不是“别处拷贝”的重复序列——拿接缝锚点当本底深度比。

    实测教训：C05 左端把一段在 123/128 条 read 里都出现的重复序列当悬垂移植过来，
    结果那 11 kb 没有任何 read 能覆盖（reads 都优先比对到别处的拷贝）。
    不能简单地“出现次数多就否决”——35× 深度下真实的唯一序列本来也会出现在 ~35 条
    read 里。所以拿**接缝处 50 bp 锚点**的出现条数当本底：重复序列的拷贝数会显著更高。

    输入示例: hangover 前缀 123 条 read 命中 / anchor 35 条 → 3.5 > 3 → 判为重复
    输出示例: True
    """
    probe = hangover[:REPEAT_KMER_BP]
    anchor_probe = anchor[-REPEAT_KMER_BP:]
    if len(probe) < REPEAT_KMER_BP or len(anchor_probe) < REPEAT_KMER_BP:
        return True
    anchor_hits = max(1, _count_reads_containing(anchor_probe, read_sequences))
    return _count_reads_containing(probe, read_sequences) > anchor_hits * REPEAT_DEPTH_FACTOR


def _pick_hangover(candidates, min_support: int, min_progress_bp: int):
    """按“序列前缀相同”的簇投票，取支持数最多且够长的悬垂。

    输入示例: candidates=[(name, "ACGT…" 800 bp) …]（同源 read 的悬垂）
    输出示例: "ACGT…" 或 None（支持数不足 / 都不够长）
    """
    clusters: dict[str, list[str]] = collections.defaultdict(list)
    for _name, sequence in candidates:
        clusters[sequence[:CLUSTER_PREFIX_BP]].append(sequence)
    if not clusters:
        return None
    # 先按支持数、再按长度定序（同分时取更长的悬垂）
    ranked = sorted(
        clusters.values(),
        key=lambda group: (len(group), max(len(sequence) for sequence in group)),
        reverse=True,
    )
    best = ranked[0]
    if len(best) < min_support:
        return None
    chosen = max(best, key=len)
    if len(chosen) < min_progress_bp:
        return None
    return chosen


def _run_minimap2(command: list[str], probe_path: Path) -> str:
    logger.debug("运行：%s", " ".join(command))
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        tail = (completed.stderr or "").strip().splitlines()[-5:]
        logger.warning("探针招募失败（%s）：%s", probe_path.name, "；".join(tail))
        return ""
    return completed.stdout


def _recruit(
    edge_sequence: str,
    reads_fastq: Path,
    read_sequences: dict[str, str],
    minimap2: str,
    settings: RecruitSettings,
    side: str,
):
    """招募一轮，返回 (悬垂序列, 支持数)；没有合格悬垂时返回 None。"""
    with tempfile.TemporaryDirectory(prefix="alnkit-recruit-") as workdir:
        probe_path = Path(workdir) / "edge.fa"
        write_fasta_records({"edge": edge_sequence}, probe_path)
        command = [
            minimap2,
            "-ax",
            "map-hifi",
            "-Y",
            "-t",
            str(max(1, settings.threads)),
            str(probe_path),
            os.fspath(reads_fastq),
        ]
        sam = _run_minimap2(command, probe_path)

    edge_length = len(edge_sequence)
    candidates = []
    for line in sam.splitlines():
        alignment = parse_alignment(line, edge_length, settings.min_mapq, settings.min_identity)
        if alignment is None:
            continue
        sequence = read_sequences.get(alignment.read_name)
        if sequence is None:
            continue
        oriented = reverse_complement(sequence) if alignment.reverse else sequence
        if side == "right":
            if alignment.reference_end < edge_length - EDGE_TOLERANCE_BP:
                continue
            hangover = oriented[alignment.query_end :]
        else:
            if alignment.reference_start > EDGE_TOLERANCE_BP:
                continue
            hangover = oriented[: alignment.query_start]
        if hangover:
            candidates.append((alignment.read_name, hangover))

    chosen = _pick_hangover(candidates, settings.min_support, settings.min_progress_bp)
    if chosen is None:
        return None
    if _looks_repetitive(chosen, edge_sequence, read_sequences):
        logger.debug("悬垂前缀是重复序列（%d bp），跳过", len(chosen))
        return None
    support = sum(
        1 for _name, sequence in candidates if sequence.startswith(chosen[:CLUSTER_PREFIX_BP])
    )
    return chosen, support


def extend_contig_ends(
    sequence: str,
    reads_fastq: Path,
    read_sequences: dict[str, str],
    minimap2: str,
    settings: RecruitSettings | None = None,
) -> tuple[str, int, int]:
    """把 contig 的两端按“招募 + 悬垂移植 + 每轮验收”延长。

    返回 ``(新序列, 右端轮数, 左端轮数)``；没有任何一轮通过验收就原样返回。
    """
    settings = settings or RecruitSettings()
    right_rounds = 0
    for _ in range(settings.max_rounds):
        recruited = _recruit(
            sequence[-settings.edge_bp :], reads_fastq, read_sequences, minimap2, settings, "right"
        )
        if recruited is None:
            break
        hangover, support = recruited
        logger.debug("右端 +%d bp（%d 条 read 支持）", len(hangover), support)
        sequence = sequence + hangover
        right_rounds += 1

    left_rounds = 0
    for _ in range(settings.max_rounds):
        recruited = _recruit(
            sequence[: settings.edge_bp], reads_fastq, read_sequences, minimap2, settings, "left"
        )
        if recruited is None:
            break
        hangover, support = recruited
        logger.debug("左端 +%d bp（%d 条 read 支持）", len(hangover), support)
        sequence = hangover + sequence  # 左端悬垂是“端点之前”的序列，直接前置
        left_rounds += 1
    return sequence, right_rounds, left_rounds
