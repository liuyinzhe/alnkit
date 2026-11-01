"""PAF 解析（D2）：minimap2 原生格式，坐标一律 0-based 半开。

与绘图模块里的 PAF 解析保持同一语义，但这里保留全部 12 列（含 ``matches`` / ``mapq``），
组装与 overlap 图需要一致率与比对长度。
"""

from __future__ import annotations

import os

__all__ = ["PAF_MIN_FIELDS", "parse_paf"]

PAF_MIN_FIELDS = 12


def parse_paf(path: str | os.PathLike[str]) -> list[dict[str, object]]:
    """读 PAF（需要 ``-c`` 才有第 10/11 列的匹配数与比对长度）。

    返回每个比对块的字段：``query, query_length, query_start, query_end, strand,
    target, target_length, target_start, target_end, matches, block_bp, identity, mapq``。
    空行、``#`` 注释、列数不足、比对长度为 0 的行会被跳过。
    """
    blocks: list[dict[str, object]] = []
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < PAF_MIN_FIELDS:
                continue
            aligned_length = int(fields[10])
            if aligned_length <= 0:
                continue
            matches = int(fields[9])
            blocks.append(
                {
                    "query": fields[0],
                    "query_length": int(fields[1]),
                    "query_start": int(fields[2]),
                    "query_end": int(fields[3]),
                    "strand": fields[4],
                    "target": fields[5],
                    "target_length": int(fields[6]),
                    "target_start": int(fields[7]),
                    "target_end": int(fields[8]),
                    "matches": matches,
                    "block_bp": aligned_length,
                    "identity": matches / aligned_length,
                    "mapq": int(fields[11]),
                }
            )
    return blocks
