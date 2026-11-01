"""核心层序列与 PAF 小工具（D1/D2）。"""

from __future__ import annotations

from pathlib import Path

import pytest

from alnkit.core import fasta, paf
from alnkit.core.errors import UsageError


def test_reverse_complement_keeps_case() -> None:
    assert fasta.reverse_complement("ACGTN") == "NACGT"
    assert fasta.reverse_complement("acgt") == "acgt"


def test_write_fasta_wraps_at_60(tmp_path: Path) -> None:
    target = tmp_path / "out.fa"
    written = fasta.write_fasta_records({"c1": "A" * 130}, target)
    assert written == 1
    lines = target.read_text(encoding="utf-8").splitlines()
    assert lines[0] == ">c1"
    assert [len(line) for line in lines[1:]] == [60, 60, 10]


def test_read_fastq_sequences_uppercases(tmp_path: Path) -> None:
    path = tmp_path / "r.fastq"
    path.write_text("@r1 desc\nacgt\n+\nIIII\n@r2\nTTTT\n+\nIIII\n", encoding="utf-8")
    assert fasta.read_fastq_sequences(path) == {"r1": "ACGT", "r2": "TTTT"}


def test_missing_fastq_raises_usage_error(tmp_path: Path) -> None:
    with pytest.raises(UsageError, match="找不到 FASTQ"):
        fasta.read_fastq_sequences(tmp_path / "nope.fastq")


def test_parse_paf_columns(tmp_path: Path) -> None:
    path = tmp_path / "x.paf"
    path.write_text(
        "# 注释\n"
        "q1\t1200\t0\t1200\t+\tt1\t32000\t500\t1700\t1180\t1200\t60\tNM:i:20\n"
        "short\t10\t0\t10\t+\tt1\t100\t0\t10\t10\t10\t60\n"
        "bad\t10\t0\t10\t+\tt1\t100\t0\t10\t1\t0\t60\n",
        encoding="utf-8",
    )
    blocks = paf.parse_paf(path)
    assert len(blocks) == 2  # 列数不足的行被跳过；block_bp=0 的第二行保留
    first = blocks[0]
    assert (first["query"], first["target"]) == ("q1", "t1")
    assert (first["query_start"], first["query_end"]) == (0, 1200)
    assert (first["target_start"], first["target_end"]) == (500, 1700)
    assert first["block_bp"] == 1200
    assert first["identity"] == pytest.approx(1180 / 1200)
    assert first["mapq"] == 60


def test_read_fasta_records_uppercases_and_handles_crlf(tmp_path: Path) -> None:
    path = tmp_path / "x.fa"
    path.write_bytes(b">seq1 desc\r\nacgt\r\n>seq2\r\nTTTT\r\n")
    assert fasta.read_fasta_records(path) == {"seq1": "ACGT", "seq2": "TTTT"}


def test_read_fasta_records_supports_gzip(tmp_path: Path) -> None:
    import gzip

    path = tmp_path / "x.fa.gz"
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(">s\nacgt\n")
    assert fasta.read_fasta_records(path) == {"s": "ACGT"}


def test_read_fasta_records_missing_file(tmp_path: Path) -> None:
    with pytest.raises(UsageError, match="找不到 FASTA"):
        fasta.read_fasta_records(tmp_path / "nope.fa")
