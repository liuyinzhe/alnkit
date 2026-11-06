"""同源位点诊断：位点是不是落在“与某条插入序列同源的基因组区域”里。

为什么需要：二代插入分析最常见的假阳性来自**外源与基因组内源序列同源**（例如外源
`insertA` 与参考里的内源同源片段有 ~99% 同源）。这些位点的 SR/PESR 支持同样很高，却根本不是
整合位点。现有 `is_homology_interference` 只比较"插入序列之间"，看不到这种情况。

两种证据来源，可单独或同时使用：

* ``--reference-fasta``：**对每个候选位点做局部回比**——取位点两侧各 ``--homology-window-bp``
  （默认 5 kb）的参考窗口，把插入序列比上去（``-x splice`` 兼顾 cDNA 样插入，``-x asm20``
  兼顾与基因组连续的插入）。只比小窗口，所以内存是 MB 级；
* ``--homology-bed``：用户直接给出已知的同源基因/区域（BED，0-based 半开）。

两者都只做**注释**，不改变位点计算，也不改变历史三张表；不指定这两个参数时不产生任何新输出。
"""

from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
import re
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from ..core.errors import RunError, UsageError
from ..core.logging import get_logger

__all__ = [
    "HomologyRegion",
    "annotate_sites",
    "load_bed_regions",
    "local_homology_regions",
    "write_site_homology_table",
]

logger = get_logger("shortread.homology")

#: 认定同源的最小一致率
DEFAULT_MIN_IDENTITY = 0.90
#: 认定同源的最小插入侧比对长度（bp）
DEFAULT_MIN_LENGTH = 200
#: 局部回比时位点两侧各取多长
DEFAULT_WINDOW_BP = 50000
#: 判定时向两侧放宽的碱基数（断点常落在同源区边界外几十 bp）
PADDING_BP = 50


@dataclass(slots=True, frozen=True)
class HomologyRegion:
    """基因组上一段“与外源同源”的区域（0-based 半开）。"""

    chrom: str
    start_0based: int
    end_0based: int
    source: str  # "insert:insertA" 或 "bed:homologA"
    identity: Optional[float]  # BED 来源没有一致率
    query_coverage: Optional[float] = None  # 插入序列被这段区域覆盖的比例


def load_bed_regions(paths: Sequence[Path]) -> List[HomologyRegion]:
    """读 BED（0-based 半开），每行 ``chrom start end [name]``；name 缺省用文件名。

    输入示例（known.bed）: ``chr1  1000000  1037000  homologA``
    输出示例: [HomologyRegion('chr1', 1000000, 1037000, 'bed:homologA', None)]
    """
    regions: List[HomologyRegion] = []
    for path in paths:
        if not Path(path).is_file():
            raise UsageError(f"--homology-bed 不存在：{os.fspath(path)}")
        for line_number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                raise UsageError(f"{os.fspath(path)}:{line_number}: BED 至少要有 3 列")
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise UsageError(f"{os.fspath(path)}:{line_number}: 坐标不是整数") from error
            if end <= start:
                raise UsageError(f"{os.fspath(path)}:{line_number}: end <= start")
            name = fields[3] if len(fields) > 3 and fields[3] else Path(path).stem
            regions.append(HomologyRegion(fields[0], start, end, f"bed:{name}", None))
    if regions:
        logger.info("同源区域（BED）: %d 段，来自 %d 个文件", len(regions), len(paths))
    return regions


def _read_insert_sequences(reference_fasta: Path, insert_names: Iterable[str]) -> Dict[str, str]:
    """从参考 FASTA 里取出插入序列本身（按名字）。"""
    import pysam

    sequences: Dict[str, str] = {}
    with pysam.FastaFile(str(reference_fasta)) as fasta:
        available = set(fasta.references)
        for name in insert_names:
            if name not in available:
                raise UsageError(f"--reference-fasta 里没有插入序列 {name}")
            sequences[name] = fasta.fetch(name).upper()
    return sequences


def _count_cs_differences(cs_tag: str) -> Tuple[int, int]:
    """从 minimap2 的 ``cs`` 标签数出替换数与 indel 碱基数（内含子不计）。

    输入示例: cs = "=ACGT*ag+tt=GG~gt1234ag=TT"
    输出示例: (1, 2)   # 1 个替换、2 个 indel 碱基
    """
    substitutions = len(re.findall(r"\*[a-z][a-z]", cs_tag))
    # 匹配串带着前导 +/-, 所以碱基数要减 1
    indel_bases = sum(len(match) - 1 for match in re.findall(r"[+-][a-z]+", cs_tag))
    return substitutions, indel_bases


def _best_local_hit(
    window_sequence: str,
    insert_sequences: Dict[str, str],
    minimap2: str,
    *,
    min_identity: float,
    min_length: int,
    threads: int,
) -> Optional[Tuple[str, float, int, int]]:
    """把插入序列比对到一段参考窗口，返回最佳命中 (插入名, 一致率, 窗口内起点, 终点)。

    输入示例: 窗口 10 kb、插入 insertA 1200 bp
    输出示例: ('insertA', 0.988, 3000, 3400)

    说明: 同时用 ``splice`` 与 ``asm20`` 两个预设（前者处理 cDNA 样插入，后者处理连续插入），
          取一致率最高的命中；一致率按**插入侧比对长度**算（splice 的 PAF ``alnlen`` 含内含子，
          直接除会严重低估，实测把 98.8% 算成 3.7%）。
    """
    if not insert_sequences:
        return None
    with tempfile.TemporaryDirectory(prefix="alnkit-homology-") as workdir:
        insert_path = Path(workdir) / "inserts.fa"
        with insert_path.open("w", encoding="utf-8") as handle:
            for name, sequence in insert_sequences.items():
                handle.write(f">{name}\n")
                for offset in range(0, len(sequence), 60):
                    handle.write(sequence[offset : offset + 60] + "\n")
        window_path = Path(workdir) / "window.fa"
        window_path.write_text(f">window\n{window_sequence}\n", encoding="utf-8")

        best: Optional[Tuple[str, float, int, int, float]] = None
        for preset in ("splice", "asm20"):
            command = [
                minimap2,
                "-x",
                preset,
                "-c",
                "--cs=long",
                "--secondary=no",
                "-t",
                str(max(1, threads)),
                str(window_path),
                str(insert_path),
            ]
            completed = subprocess.run(command, capture_output=True, text=True, check=False)
            if completed.returncode != 0:
                tail = (completed.stderr or "").strip().splitlines()[-5:]
                raise RunError("插入序列回比窗口失败：" + "；".join(tail))
            for line in completed.stdout.splitlines():
                fields = line.split("\t")
                if len(fields) < 12:
                    continue
                query = fields[0]
                query_start, query_end = int(fields[2]), int(fields[3])
                target_start, target_end = int(fields[7]), int(fields[8])
                query_aligned = query_end - query_start
                if query_aligned < min_length:
                    continue
                cs_tag = next((item[5:] for item in fields[11:] if item.startswith("cs:Z:")), "")
                substitutions, indel_bases = _count_cs_differences(cs_tag)
                identity = (
                    1.0 - (substitutions + indel_bases) / query_aligned if query_aligned else 0.0
                )
                if identity < min_identity:
                    continue
                insert_length = len(insert_sequences[query])
                coverage = query_aligned / insert_length if insert_length else 0.0
                if best is None or identity > best[1]:
                    best = (
                        query,
                        round(min(identity, 1.0), 4),
                        target_start,
                        target_end,
                        round(min(coverage, 1.0), 4),
                    )
    return best


def local_homology_regions(
    reference_fasta: Path,
    insert_names: Sequence[str],
    sites: Sequence[Tuple[str, int]],
    minimap2: str,
    *,
    window_bp: int = DEFAULT_WINDOW_BP,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_length: int = DEFAULT_MIN_LENGTH,
    threads: int = 4,
) -> List[HomologyRegion]:
    """对每个候选位点做“插入序列 vs 位点周围窗口”的局部回比，返回同源区域。

    输入示例: sites=[('chr1', 1000000), ('chr1', 2000229)], window_bp=5000
    输出示例: [HomologyRegion('chr1', 995000, 1000000, 'insert:insertA', 0.988)]
              （整合位点 chr1:2000229 附近比不中插入，所以不产生区域）

    说明: 只比小窗口，所以内存是 MB 级；不需要把插入比到 3 Gb 全基因组（splice 预设对
          hg38 的索引要 ~15 GB，普通机器会 OOM，实测过）。
    """
    import pysam

    if not sites:
        return []
    insert_sequences = _read_insert_sequences(reference_fasta, insert_names)
    regions: List[HomologyRegion] = []
    with pysam.FastaFile(str(reference_fasta)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        for chromosome, position in sites:
            if chromosome not in lengths:
                continue
            start = max(0, position - window_bp)
            end = min(lengths[chromosome], position + window_bp)
            if end <= start:
                continue
            window_sequence = fasta.fetch(chromosome, start, end).upper()
            hit = _best_local_hit(
                window_sequence,
                insert_sequences,
                minimap2,
                min_identity=min_identity,
                min_length=min_length,
                threads=threads,
            )
            if hit is None:
                continue
            insert_name, identity, target_start, target_end, coverage = hit
            regions.append(
                HomologyRegion(
                    chrom=chromosome,
                    start_0based=start + target_start,
                    end_0based=start + target_end,
                    source=f"insert:{insert_name}",
                    identity=identity,
                    query_coverage=coverage,
                )
            )
    logger.info(
        "同源区域（局部回比）: %d 段（identity >= %.2f, 插入侧长度 >= %d bp, 窗口 ±%d bp）",
        len(regions),
        min_identity,
        min_length,
        window_bp,
    )
    return regions


def annotate_sites(
    sites: Sequence[Tuple[str, int]], regions: Sequence[HomologyRegion]
) -> List[Dict[str, object]]:
    """给每个位点标注是否落在同源区域里，返回可直接写表的行。

    输入示例:
        sites = [('chr1', 1000000), ('chr1', 2000229)]
        regions = [HomologyRegion('chr1', 995000, 1000000, 'insert:insertA', 0.988)]
    输出示例:
        [{'chrom': 'chr1', 'insert_pos_1based': 1000000, 'is_homology_site': 'yes', ...},
         {'chrom': 'chr1', 'insert_pos_1based': 2000229, 'is_homology_site': 'no', ...}]

    说明: 位点坐标是 legacy 的 ``Insert_Position``（断点左侧碱基 1-based，数值上等于 0-based
          边界）；判定用该位置与区域的重叠，并向两侧各放宽 ``PADDING_BP``——断点常落在同源区
          边界外几十 bp（外源片段比内源短）。
    """
    by_chromosome: Dict[str, List[HomologyRegion]] = {}
    for region in regions:
        by_chromosome.setdefault(region.chrom, []).append(region)

    rows: List[Dict[str, object]] = []
    for chromosome, position in sites:
        hits = [
            region
            for region in by_chromosome.get(chromosome, [])
            if region.start_0based - PADDING_BP <= position <= region.end_0based + PADDING_BP
        ]
        best = max(
            hits,
            key=lambda item: (item.identity or 0.0, item.end_0based - item.start_0based),
            default=None,
        )
        rows.append(
            {
                "chrom": chromosome,
                "insert_pos_1based": position,
                "is_homology_site": "yes" if hits else "no",
                "n_regions": len(hits),
                "region_sources": ";".join(sorted({hit.source for hit in hits})) or "-",
                "best_region_chrom": best.chrom if best else "-",
                "best_region_start_1based": best.start_0based + 1 if best else "-",
                "best_region_end_1based": best.end_0based if best else "-",
                "best_region_length_bp": best.end_0based - best.start_0based if best else "-",
                "best_region_identity": (
                    best.identity if best and best.identity is not None else "-"
                ),
                "best_region_query_coverage": (
                    best.query_coverage if best and best.query_coverage is not None else "-"
                ),
            }
        )
    return rows


def write_site_homology_table(
    output_path: Path, rows: Sequence[Dict[str, object]], sample_name: str
) -> int:
    """写出 <sample>.site_homology.tsv（新表 `site_homology` 的产物）。"""
    from ..core import tables

    payload = [{"sample": sample_name, **row} for row in rows]
    # 表里不带 # 注释（解读见 docs/concepts/short-scan-output.md §5.1），第一行就是表头
    return tables.write_tsv(
        output_path,
        "site_homology",
        payload,
        sort_key=lambda row: (str(row["chrom"]), int(row["insert_pos_1based"])),
    )
