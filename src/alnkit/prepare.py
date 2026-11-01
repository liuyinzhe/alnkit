"""准备阶段（A1/A2）：合并参考、把 reads 比对到合并参考。

这两步原来是 shell（``01_merge_ref.sh`` / ``02_align_merged.sh``），这里改成 Python，
逐字节复刻 shell 的产物：

* ``reference.fa`` / ``merged_reference.fa``：原样字节（不重新折行），与 shell 的
  ``gzip -dc`` + ``cat`` 等价；
* ``reference_headers.tsv``：等价于
  ``awk '/^>/{name=substr($1,2); $1=""; sub(/^ /,""); print name"\\t"$0}'``
  ——注意 awk 在给字段赋值后会**用单空格重组 $0**，所以多空格会被压平；
* ``exogenous_uppercased.fa``：逐行去 ``\\r``，``>`` 行原样，序列行转大写；
* ``.fai``：用 ``pysam.faidx``（与 ``samtools faidx`` 同实现）。
"""

from __future__ import annotations

import gzip
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pysam

from .core.errors import RunError, UsageError
from .core.logging import get_logger

__all__ = ["AlignSettings", "MergeSettings", "align_reads", "merge_reference"]

logger = get_logger("prepare")

#: 长读长比对预设（与 02_align_merged.sh 一致）
LONG_READ_PRESET = "map-hifi"
#: 保留的次要比对数上限（-N），与旧脚本一致
SECONDARY_LIMIT = "50"


@dataclass(slots=True)
class MergeSettings:
    """``alnkit ref merge`` 的入参。"""

    reference: str = ""
    exog: tuple[str, ...] = ()
    outdir: str = "."


@dataclass(slots=True)
class AlignSettings:
    """``alnkit align reads`` 的入参。"""

    reference: str = ""
    reads: str = ""
    outdir: str = "."
    name: str = "reads_to_merged"
    threads: int = 12


def _open_reference(path: str):
    if path.endswith((".gz", ".bgz")):
        return gzip.open(path, "rb")
    return open(path, "rb")


def _write_headers_tsv(reference_fasta: Path, target: Path) -> int:
    """复刻 awk：``name\\tdescription``，多余空白压成单空格。"""
    written = 0
    with (
        open(reference_fasta, encoding="utf-8") as source,
        open(target, "w", encoding="utf-8", newline="") as sink,
    ):
        for line in source:
            if not line.startswith(">"):
                continue
            fields = line[1:].split()
            if not fields:
                sink.write("\t\n")
                written += 1
                continue
            sink.write(fields[0] + "\t" + " ".join(fields[1:]) + "\n")
            written += 1
    return written


def _copy_ensuring_trailing_newline(source, sink) -> None:
    """字节复制，但保证拼接边界处有换行（源文件缺尾部换行时补一个）。"""
    last = b""
    while chunk := source.read(1 << 20):
        sink.write(chunk)
        last = chunk[-1:]
    if last and last != b"\n":
        sink.write(b"\n")


def _uppercase_fasta(source: Path, sink) -> None:
    """把一份外源 FASTA 写进 sink：去 ``\\r``、``>`` 行原样、序列行大写。"""
    with open(source, encoding="utf-8", newline="") as handle:
        for line in handle:
            line = line.replace("\r", "")
            if line.startswith(">"):
                sink.write(line)
            else:
                sink.write(line.upper())


def merge_reference(settings: MergeSettings) -> dict[str, object]:
    """参考 + 外源 -> 合并参考（含 ``.fai`` 与 header 表）。"""
    reference = Path(settings.reference)
    if not reference.exists():
        raise UsageError(f"找不到参考 FASTA：{reference}")
    if not settings.exog:
        raise UsageError("至少给一个 --exog（要合并进参考的外源 FASTA）")
    for path in settings.exog:
        if not Path(path).exists():
            raise UsageError(f"找不到外源 FASTA：{path}")

    outdir = Path(settings.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    plain_reference = outdir / "reference.fa"
    exog_upper = outdir / "exogenous_uppercased.fa"
    merged = outdir / "merged_reference.fa"

    # 1) 解压（或原样复制）成未压缩 reference.fa
    with _open_reference(str(reference)) as source, open(plain_reference, "wb") as sink:
        shutil.copyfileobj(source, sink, length=1 << 20)
    pysam.faidx(str(plain_reference))

    # 2) header 表
    _write_headers_tsv(plain_reference, outdir / "reference_headers.tsv")

    # 3) 外源大写 + 拼接
    with open(exog_upper, "w", encoding="utf-8", newline="") as sink:
        for path in settings.exog:
            _uppercase_fasta(Path(path), sink)
    with (
        open(merged, "wb") as sink,
        open(plain_reference, "rb") as left,
        open(exog_upper, "rb") as right,
    ):
        shutil.copyfileobj(left, sink, length=1 << 20)
        shutil.copyfileobj(right, sink, length=1 << 20)
    pysam.faidx(str(merged))

    with open(merged, "rb") as handle:
        sequences = sum(1 for line in handle if line.startswith(b">"))
    logger.info("合并参考完成：%d 条序列", sequences)
    return {
        "outdir": str(outdir),
        "reference": str(plain_reference),
        "merged": str(merged),
        "sequences": sequences,
    }


def _check_tool(name: str) -> str:
    path = shutil.which(name)
    if path is None:
        raise UsageError(f"找不到 {name}（比对/排序需要它）：conda install -c bioconda {name}")
    return path


def _sort_and_index(sam_path: Path, bam_path: Path, threads: int) -> str:
    """SAM -> 排序 BAM + 索引，优先用 pysam 内置的 samtools（不依赖外部可执行文件）。

    pysam 自带 samtools（0.24 里是 1.23），``pysam.sort`` / ``pysam.index`` 等价于
    ``samtools sort`` / ``samtools index``；只是 ``@PG`` 里记录的调用者不同——
    对账用的是比对记录级 md5，不受影响。

    只在 pysam 封装不可用时才回退到 PATH 里的 samtools。
    """
    try:
        pysam.sort("-@", str(threads), "-o", str(bam_path), str(sam_path))
        pysam.index("-@", str(threads), str(bam_path))
        return f"pysam(samtools {pysam.__samtools_version__})"
    except (AttributeError, OSError, pysam.SamtoolsError) as error:  # pragma: no cover
        samtools = shutil.which("samtools")
        if samtools is None:
            raise RunError(
                f"samtools 排序/建索引失败（{error}），且 PATH 里没有 samtools"
            ) from None
        _run(
            [samtools, "sort", "-@", str(threads), "-o", str(bam_path), str(sam_path)],
            context="samtools sort",
        )
        _run([samtools, "index", "-@", str(threads), str(bam_path)], context="samtools index")
        return "samtools(CLI)"


def _run(command: list[str], *, context: str, stdout=None) -> None:
    logger.info("运行：%s", " ".join(command))
    completed = subprocess.run(command, stdout=stdout, stderr=subprocess.PIPE, check=False)
    if completed.returncode != 0:
        tail = completed.stderr.decode("utf-8", "replace").strip().splitlines()[-5:]
        raise RunError(f"{context} 失败（退出码 {completed.returncode}）：\n  " + "\n  ".join(tail))


def align_reads(settings: AlignSettings) -> dict[str, object]:
    """reads -> 合并参考（minimap2 比对 + pysam 排序/建索引），产物与旧 shell 等价。"""
    minimap2 = _check_tool("minimap2")
    for label, path in (("参考", settings.reference), ("reads", settings.reads)):
        if not Path(path).exists():
            raise UsageError(f"找不到{label}：{path}")

    outdir = Path(settings.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    index_path = outdir / "merged_reference.mmi"
    bam_path = outdir / f"{settings.name}.sorted.bam"

    _run(
        [
            minimap2,
            "-ax",
            LONG_READ_PRESET,
            "-t",
            str(settings.threads),
            "-d",
            str(index_path),
            settings.reference,
        ],
        context="minimap2 建索引",
        stdout=subprocess.DEVNULL,
    )

    # minimap2 写临时 SAM，再用 pysam（内置 samtools）排序/建索引（不依赖 PATH 里的 samtools）
    with tempfile.NamedTemporaryFile("wb", suffix=".sam", dir=outdir, delete=False) as sam:
        sam_path = Path(sam.name)
    try:
        with open(sam_path, "wb") as sink:
            _run(
                [
                    minimap2,
                    "-ax",
                    LONG_READ_PRESET,
                    "-Y",
                    "-N",
                    SECONDARY_LIMIT,
                    "-t",
                    str(settings.threads),
                    str(index_path),
                    settings.reads,
                ],
                context="minimap2 比对",
                stdout=sink,
            )
        sorter = _sort_and_index(sam_path, bam_path, settings.threads)
    finally:
        sam_path.unlink(missing_ok=True)

    logger.info("比对完成：%s（排序/索引由 %s 完成）", bam_path, sorter)
    return {
        "outdir": str(outdir),
        "bam": str(bam_path),
        "index": str(bam_path) + ".bai",
        "reference_index": str(index_path),
        "sorter": sorter,
    }
