"""PAF 解析（D2）：minimap2 原生格式，坐标一律 0-based 半开。

与绘图模块里的 PAF 解析保持同一语义，但这里保留全部 12 列（含 ``matches`` / ``mapq``），
组装与 overlap 图需要一致率与比对长度。
"""

from __future__ import annotations

import os

from .errors import UsageError

__all__ = ["PAF_MIN_FIELDS", "parse_paf"]

PAF_MIN_FIELDS = 12


def parse_paf(path: str | os.PathLike[str]) -> list[dict[str, object]]:
    """读 PAF（需要 ``-c`` 才有第 10/11 列的匹配数与比对长度）。

    返回每个比对块的字段：``query, query_length, query_start, query_end, strand,
    target, target_length, target_start, target_end, matches, block_bp, identity, mapq``。
    空行、``#`` 注释、列数不足、比对长度为 0 的行会被跳过。
    """
    blocks: list[dict[str, object]] = []
    zero_length_rows = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < PAF_MIN_FIELDS:
                continue
            aligned_length = int(fields[10])
            if aligned_length <= 0:
                zero_length_rows += 1
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
    if not blocks and zero_length_rows:
        # minimap2 不带 -c 时 PAF 第 10/11 列是 0：照这样跑下去只会得到“0 条边/0 个 contig”
        raise UsageError(
            f"{path} 里 {zero_length_rows} 行 PAF 的比对长度为 0（第 10 列）："
            "请用 `minimap2 -c`（或 -cx <预设>）重新生成"
        )
    return blocks
