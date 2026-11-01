"""FASTA/FASTQ 的小工具（D1）：反向互补、按 60 bp 折行写出、FAI 长度。

只放多个步骤都要用的最小函数；参考基因组的随机访问请直接用 ``pysam.FastaFile``。
"""

from __future__ import annotations

import gzip
import os
from collections.abc import Mapping

import pysam

from .errors import UsageError

__all__ = [
    "read_fasta_records",
    "read_fastq_sequences",
    "reverse_complement",
    "write_fasta_records",
]

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")
LINE_WIDTH = 60


def read_fasta_records(path: str | os.PathLike[str]) -> dict[str, str]:
    """读小 FASTA（支持 ``.gz``），返回 ``{记录名: 大写序列}``。

    只用于外源序列 / 组装片段这类小文件；参考基因组请用 ``pysam.FastaFile`` 随机访问。
    记录名取 ``>`` 后的第一个字段；CRLF 与空行都已处理。
    """
    opener = gzip.open if str(path).endswith((".gz", ".bgz")) else open
    records: dict[str, str] = {}
    name: str | None = None
    chunks: list[str] = []
    try:
        with opener(path, "rt", encoding="utf-8") as handle:
            for line in handle:
                line = line.rstrip("\r\n")
                if line.startswith(">"):
                    if name is not None:
                        records[name] = "".join(chunks).upper()
                    name = line[1:].split()[0] if line[1:].split() else ""
                    chunks = []
                elif line:
                    chunks.append(line)
    except FileNotFoundError:
        raise UsageError(f"找不到 FASTA：{os.fspath(path)}") from None
    if name is not None:
        records[name] = "".join(chunks).upper()
    return records


def reverse_complement(sequence: str) -> str:
    """反向互补（保留大小写）。"""
    return sequence.translate(_COMPLEMENT)[::-1]


def write_fasta_records(
    records: Mapping[str, str], path: str | os.PathLike[str], width: int = LINE_WIDTH
) -> int:
    """按输入顺序写 FASTA，每行 ``width`` bp；返回写入的记录数。"""
    target = os.fspath(path)
    os.makedirs(os.path.dirname(os.path.abspath(target)), exist_ok=True)
    with open(target, "w", encoding="utf-8", newline="") as handle:
        for name, sequence in records.items():
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), width):
                handle.write(sequence[start : start + width] + "\n")
    return len(records)


def read_fastq_sequences(path: str | os.PathLike[str]) -> dict[str, str]:
    """读 FASTQ，返回 ``{read_name: 大写序列}``（用 ``pysam.FastxFile``，与旧脚本一致）。"""
    if not os.path.exists(path):
        raise UsageError(f"找不到 FASTQ：{os.fspath(path)}")
    sequences: dict[str, str] = {}
    with pysam.FastxFile(os.fspath(path)) as handle:
        for record in handle:
            if record.name:
                sequences[record.name] = (record.sequence or "").upper()
    return sequences
