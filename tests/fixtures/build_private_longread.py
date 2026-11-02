#!/usr/bin/env python3
"""准备 ``alnkit long *`` 的私有测试数据（不进公开仓库）。

    python tests/fixtures/build_private_longread.py [--source DIR] [--dest DIR] [--trust-current]

真实测序数据**不入库**，只在本机准备。布局（默认 ``~/.local/share/alnkit-testdata``）::

    longread/
    ├── input/                      # 各阶段输入（A3 产物会被拷进来供 A5 使用）
    └── golden/
        ├── classify/   PROVENANCE.txt
        ├── breakpoints/PROVENANCE.txt
        └── assemble/   PROVENANCE.txt

黄金必须由**旧脚本**生成，否则回归退化成"自己证明自己"：

* 旧脚本缺失 → 直接失败（不静默改用 alnkit）；
* ``--trust-current`` 是显式选项，仅供旧脚本已不可得时排障，PROVENANCE 记 ``generator=alnkit``，
  ``tests/test_golden_provenance.py`` 会因此判失败；
* 每次生成前先删黄金目录，避免两次来源不同的产物混在一起。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

#: 上游数据目录（用 --source 或 ALNKIT_UPSTREAM_LONGREAD 指定，仓库里不留真实路径）
DEFAULT_SOURCE = Path(os.environ.get("ALNKIT_UPSTREAM_LONGREAD", "."))
DEFAULT_DEST = Path(
    os.environ.get("ALNKIT_TEST_DATA", Path.home() / ".local/share/alnkit-testdata")
)


def write_provenance(
    directory: Path, *, generator: str, command: list[str], script: Path | None
) -> None:
    lines = [
        f"generator={generator}",
        f"generated_at={datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"command={' '.join(command)}",
    ]
    if script is not None:
        lines += [f"script={script}", f"script_sha256={hashlib.sha256(script.read_bytes()).hexdigest()}"]
    (directory / "PROVENANCE.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_stage(
    *,
    name: str,
    golden: Path,
    script: Path,
    upstream_args: list[str],
    alnkit_args: list[str],
    trust_current: bool,
    env: dict[str, str] | None = None,
) -> int:
    """用旧脚本（默认）或本包（--trust-current）生成一个阶段的黄金。"""
    if golden.exists():
        shutil.rmtree(golden)
    golden.mkdir(parents=True)

    if not trust_current:
        if not script.exists():
            return _fail(f"{name}: 找不到旧脚本 {script}；确要自证请显式加 --trust-current")
        command = [sys.executable, str(script), *upstream_args]
        generator, provenance_script = "upstream", script
    else:
        command = [sys.executable, "-m", "alnkit.cli", *alnkit_args]
        generator, provenance_script = "alnkit", None
        print(f"警告：{name} 由 alnkit 自己生成（--trust-current），不能用于对账", file=sys.stderr)

    command += ["--outdir", str(golden)]
    completed = subprocess.run(command, check=False, env={**os.environ, **(env or {})})
    if completed.returncode != 0:
        return _fail(f"{name} 生成失败（退出码 {completed.returncode}）")
    write_provenance(golden, generator=generator, command=command, script=provenance_script)
    print(f"黄金 {name}: {sorted(p.name for p in golden.iterdir())}")
    return 0


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="准备 long 链路的私有测试数据")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--dest", type=Path, default=DEFAULT_DEST)
    parser.add_argument(
        "--exog-names", default=None,
        help="合并参考里外源序列的名字，逗号分隔（真实名只写进本地 manifest，不入库）",
    )
    parser.add_argument("--query-chr", default=None, help="分类时作为 query 的插入序列名")
    parser.add_argument(
        "--trust-current", action="store_true",
        help="旧脚本不可用时用当前实现生成（失去对账意义，测试会拒绝）",
    )
    args = parser.parse_args(argv)

    source, dest = args.source, args.dest / "longread"
    scripts = source / "scripts"
    bam = source / "02_align_merged/reads_to_merged.sorted.bam"
    fastq = source / "data/rawdata/input.fastq"
    fai = source / "01_merge_ref/merged_reference.fa.fai"
    for path in (bam, fastq, fai):
        if not path.exists():
            return _fail(f"找不到 {path}")

    inputs = dest / "input"
    inputs.mkdir(parents=True, exist_ok=True)
    for path in (bam, fastq, fai):
        shutil.copy2(path, inputs / path.name)

    # 真实序列名只落在本地 manifest（不入库）：测试从它读，仓库里只出现占位名
    exog_names = [name for name in (args.exog_names or "").split(",") if name]
    if not exog_names:
        return _fail("需要 --exog-names（合并参考里外源序列名，逗号分隔）")
    manifest = {
        "exog_names": exog_names,
        "query_chr": args.query_chr or "target_ref",
    }
    (inputs / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    # --- A3 classify ---
    classify_args = [
        "--bam", str(inputs / bam.name), "--reads", str(inputs / fastq.name),
        "--merged-fai", str(inputs / fai.name), "--exog", "insertA", "--exog", "insertB",
    ]
    code = run_stage(
        name="classify", golden=dest / "golden/classify",
        script=scripts / "03_classify_reads.py",
        upstream_args=classify_args, alnkit_args=["long", "classify", *classify_args],
        trust_current=args.trust_current,
    )
    if code:
        return code

    # A5 的输入就是 A3 的产物：从黄金里拷进 input/（早先漏了这步，A5 只能靠手工准备）
    for name in ("alignment_segments.tsv", "read_classification.tsv", "related_reads.fastq"):
        shutil.copy2(dest / "golden/classify" / name, inputs / name)
    edges = source / "04_overlap_graph/overlap_edges.tsv"
    if not edges.exists():
        return _fail(f"找不到 {edges}")
    shutil.copy2(edges, inputs / edges.name)

    # --- A5 breakpoints ---
    breakpoints_args = [
        "--segments", str(inputs / "alignment_segments.tsv"),
        "--classification", str(inputs / "read_classification.tsv"),
        "--overlaps", str(inputs / edges.name),
        "--reads", str(inputs / "related_reads.fastq"),
        "--merged-fai", str(inputs / fai.name),
    ]
    code = run_stage(
        name="breakpoints", golden=dest / "golden/breakpoints",
        script=scripts / "breakpoints.py",
        upstream_args=breakpoints_args,
        alnkit_args=["long", "breakpoints", *breakpoints_args],
        trust_current=args.trust_current,
    )
    if code:
        return code
    golden_dir = dest / "golden/breakpoints"
    pools = golden_dir / "class_read_pools.tsv"
    if pools.exists():
        # 该列是 outdir 的绝对路径，随机器变；归一化后再当黄金
        pools.write_text(
            pools.read_text(encoding="utf-8").replace(str(golden_dir), "<OUTDIR>"), encoding="utf-8"
        )

    # --- A6 assemble：输入是 classes/ + overlaps.paf ---
    paf = source / "04_overlap_graph/overlaps.paf"
    classes = source / "05_breakpoints/classes"
    if not paf.exists() or not classes.exists():
        return _fail(f"缺少 {paf} 或 {classes}")
    shutil.copy2(paf, inputs / paf.name)
    if (inputs / "classes").exists():
        shutil.rmtree(inputs / "classes")
    shutil.copytree(classes, inputs / "classes")

    assemble_args = ["--classes-dir", str(inputs / "classes"), "--paf", str(inputs / paf.name)]
    code = run_stage(
        name="assemble", golden=dest / "golden/assemble",
        script=scripts / "overlap_assembly.py",
        upstream_args=assemble_args,
        alnkit_args=["long", "assemble", *assemble_args],
        trust_current=args.trust_current,
        # 上游老版本对 set 迭代顺序敏感；固定 seed 后跨机可比（上游 20b4dda 起已确定性）
        env={"PYTHONHASHSEED": "0"},
    )
    if code:
        return code

    # --- A7 locate：需要参考基因组（2.3 GB，用软链接）+ 外源 fasta + 各类组装 ---
    reference = source / "01_merge_ref/reference.fa"
    exog_files = [source / "data/exog/insertA.fa", source / "data/exog/insertB.fa"]
    if reference.exists() and all(path.exists() for path in exog_files):
        for suffix in ("", ".fai"):
            link = inputs / ("reference.fa" + suffix)
            if link.is_symlink() or link.exists():
                link.unlink()
            link.symlink_to(Path(str(reference) + suffix))
        (inputs / "exog").mkdir(exist_ok=True)
        for path in exog_files:
            shutil.copy2(path, inputs / "exog" / path.name)
        if (inputs / "assembly").exists():
            shutil.rmtree(inputs / "assembly")
        shutil.copytree(classes.parent.parent / "06_assembly", inputs / "assembly")

        locate_args = [
            "--regions", str(dest / "golden/breakpoints/breakpoint_classes.tsv"),
            "--sites", str(dest / "golden/breakpoints/breakpoint_sites.tsv"),
            "--reference-fasta", str(inputs / "reference.fa"),
            "--exog", str(inputs / "exog/insertA.fa"),
            "--exog", str(inputs / "exog/insertB.fa"),
            "--assembly-dir", str(inputs / "assembly"),
            "--threads", "8",
        ]
        code = run_stage(
            name="locate", golden=dest / "golden/locate",
            script=scripts / "ribbon_positions.py",
            upstream_args=locate_args, alnkit_args=["long", "locate", *locate_args],
            trust_current=args.trust_current,
        )
        if code:
            return code
        # 文本产物里含 outdir 绝对路径（class_position_summary / summary.md）；归一化后再当黄金
        golden_locate = dest / "golden/locate"
        for path in list(golden_locate.rglob("*.tsv")) + list(golden_locate.rglob("*.md")):
            path.write_text(
                path.read_text(encoding="utf-8").replace(str(golden_locate), "<OUTDIR>"),
                encoding="utf-8",
            )
    else:
        print("跳过 A7 黄金：缺 reference.fa 或外源 fasta", file=sys.stderr)

    print(f"\n完成：{dest}\n用法：export ALNKIT_TEST_DATA={args.dest}")
    return 0


def _fail(message: str) -> int:
    print(f"错误：{message}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
