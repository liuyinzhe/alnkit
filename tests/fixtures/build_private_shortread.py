#!/usr/bin/env python3
"""准备 ``alnkit short scan`` 的私有测试数据（不进公开仓库）。

真实测序数据**不入库**，只在本机准备（路径用 `--bam-dir` 指定）：

    python tests/fixtures/build_private_shortread.py [--source DIR] [--dest DIR]

产物布局（默认 ``~/.local/share/alnkit-testdata``）::

    shortread/
    ├── bam/
    │   ├── sampleA.sorted.bam     # 未 markdup（dup 列应为 0）
    │   ├── sampleA.markdup.bam    # 已 markdup（dup 列非零）
    │   └── ...（其余样本同理）
    └── golden/
        ├── owned/ markdup/ consensus/                # 三套黄金产物（列名/逐字节）

黄金产物由**旧脚本** ``parse_bam.py`` v3.0.0 生成，用于证明移植后输出不变。
若旧脚本已不可得，可用 ``--trust-current`` 用当前实现重新生成（会失去对账意义，
只用于本地排障）。
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

#: 上游数据目录（用 --source 或 ALNKIT_UPSTREAM_SHORTREAD 指定，仓库里不留真实路径）
DEFAULT_SOURCE = Path(os.environ.get("ALNKIT_UPSTREAM_SHORTREAD", "."))
DEFAULT_DEST = Path(
    os.environ.get("ALNKIT_TEST_DATA", Path.home() / ".local/share/alnkit-testdata")
)

#: (case 名, BAM 文件名, 额外参数)
#: (case 名, BAM 文件名, 额外参数模板)；插入序列名运行时从 manifest 取
CASES = [
    ("owned", "sampleA.sorted.bam", {"merge_window_bp": "10"}),
    ("markdup", "sampleA.markdup.bam", {}),
    ("consensus", "sampleA.sorted.bam", {"policy": "consensus", "merge_window_bp": "10"}),
]


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description="准备 short scan 的私有测试数据")
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument(
        "--bam-dir",
        type=Path,
        default=None,
        help="上游 BAM 目录（含 *.sorted.bam，子目录 markdup/ 含 *.markdup.bam）",
    )
    parser.add_argument("--dest", type=Path, default=DEFAULT_DEST)
    parser.add_argument(
        "--query-chr", default=None, help="上游 BAM 里作为 query 的插入序列名（真实名不出现在仓库）"
    )
    parser.add_argument(
        "--insert-chrs", default=None, help="逗号分隔的全部插入序列名；缺省等于 --query-chr"
    )
    parser.add_argument(
        "--trust-current",
        action="store_true",
        help="旧脚本不可用时，用当前实现生成黄金（失去对账意义）",
    )
    args = parser.parse_args(argv)

    source, dest = args.source, args.dest / "shortread"
    bam_dir = args.bam_dir or (source / "bams")
    if not bam_dir.exists():
        return _fail(
            f"找不到 {bam_dir}：用 --bam-dir 指定上游 BAM 目录"
            "（含 *.sorted.bam 与 markdup/*.markdup.bam）"
        )

    (dest / "bam").mkdir(parents=True, exist_ok=True)
    # 上游 BAM 的真实文件名不出现在本仓库：这里按排序结果映射成 sampleA/sampleB
    plain = sorted(bam_dir.glob("*.sorted.bam"))
    marked = sorted((bam_dir / "markdup").glob("*.markdup.bam"))
    if not plain or not marked:
        return _fail(f"{bam_dir} 或 {bam_dir / 'markdup'} 下没有 *.sorted.bam / *.markdup.bam")
    labels = [f"sample{chr(ord('A') + index)}" for index in range(len(plain))]
    for label, path in zip(labels, plain, strict=False):
        shutil.copy2(path, dest / "bam" / f"{label}.sorted.bam")
    for label, path in zip(labels, marked, strict=False):
        shutil.copy2(path, dest / "bam" / f"{label}.markdup.bam")

    old_script = source / "parse_bam.py"
    if not old_script.exists() and not args.trust_current:
        return _fail(f"找不到 {old_script}；确要自证请显式加 --trust-current")

    # 真实 contig 名只落在本地 manifest 里（仓库与黄金之外不出现）
    query_chr = args.query_chr or "target_ref"
    insert_chrs = args.insert_chrs or query_chr
    manifest = {
        "query_chr": query_chr,
        "insert_chrs": [name for name in insert_chrs.split(",") if name],
    }
    (dest / "manifest.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )

    for case, bam_name, extra_kwargs in CASES:
        out = dest / "golden" / case
        if out.exists():
            shutil.rmtree(out)  # 免得两次来源不同的产物混在一起
        out.mkdir(parents=True)
        extra: list[str] = ["--insert-chrs", ",".join(manifest["insert_chrs"])]
        if extra_kwargs.get("merge_window_bp"):
            extra += ["--merge-window-bp", extra_kwargs["merge_window_bp"]]
        if extra_kwargs.get("policy"):
            extra += ["--insert-origin-policy", extra_kwargs["policy"]]
        stage_args = [
            "--bam",
            str(dest / "bam" / bam_name),
            "--query-chr",
            manifest["query_chr"],
            "--log-level",
            "ERROR",
            *extra,
        ]
        if args.trust_current:
            command = [sys.executable, "-m", "alnkit.cli", "short", "scan", *stage_args]
            generator = "alnkit"
            print(
                f"警告：{case} 由 alnkit 自己生成（--trust-current），不能用于对账", file=sys.stderr
            )
        else:
            command = [sys.executable, str(old_script), *stage_args]
            generator = "upstream"
        command += ["--out-dir", str(out)]
        completed = subprocess.run(command, check=False)
        if completed.returncode != 0:
            return _fail(f"{case} 生成失败（退出码 {completed.returncode}）")
        write_provenance(
            out,
            generator=generator,
            command=command,
            script=old_script if generator == "upstream" else None,
        )
        print(f"黄金 {case}: {sorted(p.name for p in out.iterdir())}")

    print(f"\n完成。用法：export ALNKIT_TEST_DATA={args.dest}")
    return 0


def write_provenance(
    directory: Path, *, generator: str, command: list[str], script: Path | None
) -> None:
    lines = [
        f"generator={generator}",
        f"generated_at={datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        f"command={' '.join(command)}",
    ]
    if script is not None:
        lines += [
            f"script={script}",
            f"script_sha256={hashlib.sha256(script.read_bytes()).hexdigest()}",
        ]
    (directory / "PROVENANCE.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")


def _fail(message: str) -> int:
    print(f"错误：{message}", file=sys.stderr)
    return 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
