#!/usr/bin/env python3
"""共线性丝带图（C1）：「一条序列 vs 多个片段」的比对关系图。

用于看清某个外源/参考片段出现在哪些组装片段的哪里，或反过来看清某条组装片段上都有什么
（外源拷贝、参考块、类组装片段…）。两种输入方式：

* 方式 1（默认）：给 ``reference`` + 一个或多个 ``query``，内部跑 minimap2；
* 方式 2：给 ``paf``（现成比对结果），直接画图，不需要 fasta，也不依赖 minimap2。

产物：``<outdir>/<name>.png|.pdf``、``<outdir>/<name>.blocks.tsv``（坐标 0-based 半开），
可选 ``.placements.tsv`` / ``.regions.tsv``。

来源：``synteny_ribbon`` v1.10.1（git a41c840）单文件实现的忠实移植；``--paf`` 模式下
``blocks.tsv`` 与旧脚本逐字节一致（见 ``tests/test_plot_ribbon.py`` 的黄金回归）。
"""

from __future__ import annotations

import collections
import os
import shlex
import subprocess
from dataclasses import dataclass

from ..core.errors import UsageError

VERSION = "1.10.1"

# 每行一个颜色（超过则循环）
COLOR_CYCLE = [
    "#4C72B0",
    "#DD8452",
    "#55A868",
    "#C44E52",
    "#8172B3",
    "#937860",
    "#DA8BC3",
    "#8C8C8C",
    "#CCB974",
    "#64B5CD",
]

# 参考行在图里的占位名（--order / --order-file 里用这个名字指代参考行）
REFERENCE_TOKEN = "reference"

# --preset auto 的档位：从长到短，第一个满足 `长度 >= 下界` 的档位生效。
# 依据（见 README「什么时候用 --preset auto」的实测表）：asm5/asm10/asm20 带 min_dp_max=200，
# 短于约 200 bp 的比对会被直接丢弃；40–999 bp 用 -k11 -w5 的小种子参数召回率明显更高。
PRESET_TIERS = (
    (1000, ["-x", "asm20"]),
    (40, ["-x", "map-ont", "-k11", "-w5", "-n1", "-m10"]),
)
# auto 不接受更短的记录：到此长度只能用 -k8 一类的小 k，重复多的参考上假阳性会暴涨
MIN_AUTO_LENGTH = 40
# 同一文件里最长/最短记录的比值超过它就认为"跨度过大"，要求用户按长度拆分
MAX_LENGTH_SPAN_RATIO = 20
# 所有内部比对都加：
#   -c              PAF 输出匹配数/比对长度（本工具用它们算 identity），必需
#   --secondary=yes sr 预设默认 --secondary=no（同一条 query 只留一处，多拷贝会塌缩），
#                   而丝带图要画多拷贝；-x 之后写的显式选项会覆盖预设值，所以这几项必须排在预设参数后面。
ALWAYS_MINIMAP2_ARGS = ["-c", "--secondary=yes", "-N", "500"]


### 输入准备 ###


def parse_spec(spec_text, where=""):
    """解析 `标签=fasta路径[:片段名]` 形式的输入说明（片段名可省略）。

    说明：标签会作为该行在图上的显示名（所以可以把"参考基因组坐标"这类说明写进标签，便于读图），
          但标签里不能再出现 =（= 只当"标签/路径"的分隔符）。
          where 只用于报错前缀（从列表文件读时带"文件 第N行："），直接给命令行参数时留空。

    输入示例: parse_spec("insert_alpha=/data/exog.fa:insert_alpha")
    输出示例: ("insert_alpha", "/data/exog.fa", "insert_alpha")
    """
    label, _, path_part = spec_text.partition("=")
    if not label.strip():
        # 比对到旧实现：旧版允许空标签，结果是图上出现一行没有名字的轨道；这里提前拦住
        raise UsageError(
            "{0}标签不能为空（格式是 标签=fasta路径[:片段名]）：{1}".format(where, spec_text)
        )
    if not path_part:
        raise UsageError(
            "{0}输入格式应为 标签=fasta路径[:片段名]，收到：{1}".format(where, spec_text)
        )
    if "=" in path_part:
        # 常见踩坑：把标签写成 "insert_alpha = 1-1200"，于是第一个 = 之后又冒出一个 =，路径整体错位
        raise UsageError(
            "{0}说明里的标签不能包含 =（格式是 标签=fasta路径[:片段名]，= 只能当分隔符用）："
            "{1}\n  建议把标签写成 insert_alpha (1-1200) 这类形式".format(where, spec_text)
        )
    fasta_path, separator, contig = path_part.rpartition(":")
    if not separator or "/" in contig or not contig:
        fasta_path, contig = path_part, None
    return label, fasta_path, contig


def list_records(fasta_path):
    """列出 fasta 里的记录（名称与长度）。

    输入示例: list_records("reference_side.fasta")
    输出示例: [("insert_alpha", 1353), ("insert_beta", 642), ("flank_reference_region", 13762)]
    """
    records = []
    name = None
    length = 0
    with open(fasta_path) as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, length))
                name = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if name is not None:
        records.append((name, length))
    return records


def upper_case_fasta(text):
    """把 fasta 文本里的序列行转成大写（表头原样保留）。

    输入示例: upper_case_fasta(">s1\\nacgt\\n")
    输出示例: ">s1\\nACGT\\n"
    说明：本机 minimap2 2.31 不做小写软屏蔽判定（实测大小写结果一致），但较旧版本会把小写当软屏蔽
         而忽略种子；统一转大写可以跨版本/跨工具都稳定，所以脚本内一律先转大写再比对。
    """
    output_lines = []
    for line in text.splitlines(True):
        if line.startswith(">"):
            output_lines.append(line)
        else:
            output_lines.append(line.upper())
    return "".join(output_lines)


def extract_single_record(fasta_path, contig_name, output_path):
    """把一条记录单独写成 fasta（顺带转大写）。

    输入示例: extract_single_record("exog.fa", "insert_alpha", "ref.fasta")
    输出示例: 生成 ref.fasta（序列已大写），返回该序列长度
    """
    import pysam

    with pysam.FastaFile(fasta_path) as fasta:
        sequence = fasta.fetch(contig_name).upper()
    with open(output_path, "w") as output_file:
        output_file.write(">{0}\n{1}\n".format(contig_name, sequence))
    return len(sequence)


def write_upper_case_fasta(fasta_path, output_path):
    """把整个 fasta（可能多条记录）转大写写到新文件，返回记录数。

    输入示例: write_upper_case_fasta("side.fasta", "side_upper.fasta")  -> 输出示例: 5
    """
    with open(fasta_path) as input_file, open(output_path, "w") as output_file:
        text = upper_case_fasta(input_file.read())
        output_file.write(text)
    return len(list_records(output_path))


### 行顺序（长序列可放末尾 + 按列表排布）###


def read_order_file(order_file):
    """读"顺序列表文件"：每行一个片段名（`reference` 代表参考行），空行与 `#` 注释跳过。

    输入示例: read_order_file("order.txt")  其中内容为 "class1\\n# 注释\\nreference\\n"
    输出示例: ["class1", "reference"]
    """
    names = []
    with open(order_file) as order_reader:
        for line in order_reader:
            name = line.split("#", 1)[0].strip()
            if name:
                names.append(name)
    return names


def read_query_file(query_file):
    """读"查询列表文件"：每行一个 `标签=fasta路径[:片段名]`，空行与 `#` 注释跳过。

    与逐条写 --query 完全等效（顺序也一致：文件从上到下 = 命令行从左到右）。

    输入示例: read_query_file("queries.txt")  其中内容为 "a=/d/a.fa:asm_a\n# 注释\nb=/d/b.fa\n"
    输出示例: ["a=/d/a.fa:asm_a", "b=/d/b.fa"]
    """
    specs = []
    with open(query_file) as query_reader:
        for line_number, line in enumerate(query_reader, 1):
            spec = line.split("#", 1)[0].strip()
            if not spec:
                continue
            # 先按同一套规则校验，报错时能指出是哪份文件的第几行
            parse_spec(spec, "{0} 第 {1} 行：".format(query_file, line_number))
            specs.append(spec)
    if not specs:
        raise UsageError(
            "{0} 里没有查询片段：每行写一个 标签=fasta路径[:片段名]"
            "（空行与 # 注释会被跳过）".format(query_file)
        )
    return specs


def resolve_row_order(query_names, order_text=None, reference_position="top", order_file_names=()):
    """决定图中**自上而下**的行顺序（返回行名列表，参考行用 "reference" 表示）。

    输入示例:
        resolve_row_order(["assembly_a", "assembly_b"], "assembly_b,reference,assembly_a")          -> ["assembly_b", "reference", "assembly_a"]
        resolve_row_order(["assembly_a", "assembly_b"], None, "bottom")                     -> ["class1", "class2", "reference"]
        resolve_row_order(["assembly_a", "assembly_b"], None, "top", ("reference", "class2")) -> ["reference", "class2", "class1"]
    说明:
        * order_text 与 order_file_names 里都可以用行名，特殊名字 reference 代表参考行
          （同时给两者时，--order 的条目在前、列表文件的条目在后）；
        * 未列出的行按原顺序接在最后；
        * 名字写错（不在行名里、也不是 reference）会返回告警列表，不会静默丢行；
        * 给了顺序时，reference_position 只作为"参考行没被显式指定位置"时的兜底。
    """
    if order_text and order_file_names:
        raise UsageError("--order 与 --order-file 只能用一个（两种给行顺序的方式互斥）")
    tokens = []
    if order_text:
        tokens += [piece.strip() for piece in order_text.split(",") if piece.strip()]
    tokens += [str(name).strip() for name in order_file_names if str(name).strip()]

    ordered = []
    warnings = []
    for token in tokens:
        if token == REFERENCE_TOKEN:
            if REFERENCE_TOKEN not in ordered:
                ordered.append(REFERENCE_TOKEN)
        elif token in query_names:
            if token not in ordered:
                ordered.append(token)
        else:
            warnings.append(
                "顺序列表里的名字不存在，已忽略：{0}（可用行名：{1}）".format(
                    token, ", ".join(query_names) if query_names else "无"
                )
            )
    if REFERENCE_TOKEN not in ordered:
        if reference_position == "bottom":
            ordered = ordered + list(query_names) + [REFERENCE_TOKEN]
        else:
            ordered = [REFERENCE_TOKEN] + ordered + list(query_names)
    else:
        ordered = ordered + [name for name in query_names if name not in ordered]
    seen = set()
    deduped = []
    for name in ordered:
        if name not in seen:
            seen.add(name)
            deduped.append(name)
    return deduped, warnings


### 比对：内部比对 或 直接读 PAF ###


def build_minimap2_command(reference_fasta, query_fasta, threads, preset_args):
    """拼内部比对的 minimap2 命令。

    说明：-x 及其参数必须排在其它选项前面（minimap2 先应用预设，后面的显式选项再覆盖），
          所以 ALWAYS_MINIMAP2_ARGS 固定放在 preset_args 之后。

    输入示例: build_minimap2_command("ref.fa", "q.fa", 8, ["-x", "asm20"])
    输出示例: ["minimap2", "-x", "asm20", "-c", "--secondary=yes", "-N", "500", "-t", "8", "ref.fa", "q.fa"]
    """
    return (
        ["minimap2"]
        + list(preset_args)
        + ALWAYS_MINIMAP2_ARGS
        + ["-t", str(threads), reference_fasta, query_fasta]
    )


def align_file_to_reference(reference_fasta, query_fasta, threads, preset_args):
    """把一个查询 fasta（可含多条记录）用该文件对应的参数比对到参考 fasta，返回逐块结果。

    输入示例: align_file_to_reference("ref.fa", "queries.fa", 8, ["-x", "asm20"])
    输出示例: [{"record": "insert_alpha", "q_start": 0, "q_end": 1353, "r_start": 100, "r_end": 1453,
                "strand": "+", "block_bp": 1353, "identity": 1.0}]
    说明：minimap2 非 0 退出（预设名写错、文件读不了等）时直接把 stderr 带出来，
          否则会变成"没有比对块"这种难查的错。
    """
    command = build_minimap2_command(reference_fasta, query_fasta, threads, preset_args)
    finished = subprocess.run(command, capture_output=True, text=True)
    if finished.returncode != 0:
        raise UsageError(
            "minimap2 比对失败（退出码 {0}）：{1}\n  命令：{2}".format(
                finished.returncode, (finished.stderr or "").strip(), " ".join(command)
            )
        )
    return parse_paf_text(finished.stdout)


def tier_of_length(length_bp):
    """长度 → PRESET_TIERS 的档位序号；短于 MIN_AUTO_LENGTH 返回 None。

    输入示例: tier_of_length(1200) -> 0   tier_of_length(600) -> 1   tier_of_length(20) -> None
    """
    for index, (lower_bound, _args) in enumerate(PRESET_TIERS):
        if length_bp >= lower_bound:
            return index
    return None


def parse_query_preset(text):
    """解析 --query-preset 的 `标签=预设[ 额外参数]`（额外参数可省略）。

    输入示例: parse_query_preset("inserts=sr")              -> ("inserts", ["-x", "sr"])
              parse_query_preset("flank=map-ont -k11 -w5")  -> ("flank", ["-x", "map-ont", "-k11", "-w5"])
    """
    label, _, spec = text.partition("=")
    options = shlex.split(spec)
    if not options:
        raise UsageError(
            "--query-preset 应写成 标签=预设[ 额外参数]，收到：{0}\n"
            '  例：--query-preset "inserts=sr" 或 --query-preset "flank=map-ont -k11 -w5 -n1 -m10"'.format(
                text
            )
        )
    return label, ["-x"] + options


def label_text(label):
    """标签写进消息/日志时的显示形式（空标签显示成 <空>，避免出现"标签 ）"这种半句话）。"""
    return label or "<空>"


def collect_preset_overrides(preset_texts, available_labels):
    """把 --query-preset 的文本列表解成 {标签: 参数列表}；标签写错/重复就报错（不静默忽略）。

    输入示例: collect_preset_overrides(["inserts=sr"], {"region", "inserts"})
    输出示例: {"inserts": ["-x", "sr"]}
    说明：键只能是标签（--query 的 `=` 左边），不是 fasta 记录名，也不是文件名。
    """
    overrides = {}
    for text in preset_texts or []:
        label, preset_args = parse_query_preset(text)
        if label not in available_labels:
            raise UsageError(
                "--query-preset 里的标签 {0} 不存在。本次可用的查询标签：{1}\n"
                "  提示：键是 --query 左边那一段（标签），不是 fasta 里的记录名；参考侧标签也不适用。".format(
                    label_text(label),
                    ", ".join(sorted(label_text(name) for name in available_labels)) or "（无）",
                )
            )
        if label in overrides:
            raise UsageError("--query-preset 给标签 {0} 指定了两次，请只留一个".format(label))
        overrides[label] = preset_args
    return overrides


def resolve_file_preset(fasta_path, labels, lengths, preset_option, overrides, display_path=None):
    """决定一个查询 fasta 用哪套 minimap2 参数，返回 (参数列表, 来源说明)。

    优先级：--query-preset（按标签） > --preset auto（按长度分档） > 全局 --preset。
    auto 下：有记录 < MIN_AUTO_LENGTH → 报错；最长/最短 > MAX_LENGTH_SPAN_RATIO → 报错（要求拆分）。

    display_path 只用于报错信息（多记录文件比对用的是 outdir 里转大写后的副本，报错时给用户原路径更好读）。
    输入示例: resolve_file_preset("q.fa", ["inserts"], [600, 1200], "auto", {})
    输出示例: (["-x", "map-ont", "-k11", "-w5", "-n1", "-m10"], "auto：最短记录 600 bp")
    """
    picked = {label: overrides[label] for label in labels if label in overrides}
    if len({tuple(preset_args) for preset_args in picked.values()}) > 1:
        raise UsageError(
            "同一个 fasta 被多个标签引用，但 --query-preset 给了不同参数：{0}\n  {1}\n"
            "  这个文件只会比对一次，请统一参数，或把它拆成两个文件分别给 --query".format(
                display_path or fasta_path,
                "；".join(
                    "{0} → {1}".format(label, " ".join(preset_args))
                    for label, preset_args in sorted(picked.items())
                ),
            )
        )
    if picked:
        label, preset_args = sorted(picked.items())[0]
        return list(preset_args), "--query-preset（标签 {0}）".format(label_text(label))
    if preset_option != "auto":
        return ["-x", preset_option], "--preset {0}".format(preset_option)

    shown = display_path or fasta_path
    shortest, longest = min(lengths), max(lengths)
    if shortest < MIN_AUTO_LENGTH:
        raise UsageError(
            "{0} 里有记录短于 {1} bp（最短 {2} bp），--preset auto 不对这类极短片段自动选参：\n"
            "  {1} bp 以下只能用 -k8 一类的小 k 参数，重复较多的参考上假阳性会暴涨。\n"
            "  要画这类片段请显式指定参数，例如：\n"
            '    --query-preset "<标签>=map-ont -k8 -w4 -n1 -m8" --min-block-bp 20\n'
            "  也可以先用 --paf 模式自行比对，再读 PAF 画图".format(
                shown, MIN_AUTO_LENGTH, shortest
            )
        )
    if longest / shortest > MAX_LENGTH_SPAN_RATIO:
        raise UsageError(
            "{0} 里记录长度跨度过大（{1:,}–{2:,} bp，{3} 条；标签：{4}），"
            "--preset auto 只能给一个文件选一套参数：\n"
            "  请按长度拆成多个 fasta 分别给 --query（每个文件一个标签），例如：\n"
            '    --query "long=<长片段.fa>" --query "short=<短片段.fa>" --preset auto ...\n'
            "  拆文件：seqkit seq -m 1000 mixed.fa > long.fa ; seqkit seq -M 999 mixed.fa > short.fa\n"
            "  或分别比对后合并 PAF 画图：\n"
            "    minimap2 -cx asm20 ref.fa long.fa > a.paf\n"
            "    minimap2 -cx map-ont -k11 -w5 -n1 -m10 --secondary=yes ref.fa short.fa > b.paf\n"
            "    cat a.paf b.paf > all.paf ; synteny_ribbon.py --paf all.paf ...\n"
            '  如果确认想用同一套参数，用 --query-preset "<标签>=<预设>" 显式指定即可跳过这个检查'.format(
                shown, shortest, longest, len(lengths), ", ".join(sorted(set(labels)))
            )
        )
    # 文件内跨档但跨度不大时，按最短记录选档：短档（小 k）参数对长记录也安全（实测到 5 kb 仍然满召回），
    # 反之 asm 预设会丢掉/截短短记录。
    return list(PRESET_TIERS[tier_of_length(shortest)][1]), "auto：最短记录 {0:,} bp".format(
        shortest
    )


def parse_paf_text(paf_text):
    """解析 PAF 文本为比对块列表（也可用于直接读现成 PAF 文件）。

    输入示例: parse_paf_text("q\\t8000\\t0\\t8000\\t+\\tref\\t38291\\t100\\t8100\\t7900\\t8000\\t60\\n")
    输出示例: [{"record": "q", "q_length": 8000, "q_start": 0, "q_end": 8000, "r_name": "ref",
                "r_length": 38291, "r_start": 100, "r_end": 8100, "strand": "+",
                "block_bp": 8000, "identity": 0.9875, "mapq": 60}]
    说明：PAF 列依次是 qname qlen qstart qend strand tname tlen tstart tend nmatch alnlen mapq；
          identity = nmatch/alnlen（所以生成 PAF 时 **必须加 -c**，否则这两列是 0）。
    """
    blocks = []
    for line in paf_text.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            continue
        aligned_length = int(fields[10])
        blocks.append(
            {
                "record": fields[0],
                "q_length": int(fields[1]),
                "q_start": int(fields[2]),
                "q_end": int(fields[3]),
                "r_name": fields[5],
                "r_length": int(fields[6]),
                "r_start": int(fields[7]),
                "r_end": int(fields[8]),
                "strand": fields[4],
                "block_bp": aligned_length,
                "matches": int(fields[9]),
                "identity": int(fields[9]) / aligned_length if aligned_length else 0.0,
                "mapq": int(fields[11]),
            }
        )
    return blocks


def read_paf_file(paf_path):
    """读现成的 PAF 文件。输入示例: read_paf_file("aln.paf") -> 输出示例: 比对块列表"""
    with open(paf_path) as paf_reader:
        return parse_paf_text(paf_reader.read())


### 绘图 ###


def bezier_points(start, control1, control2, end, samples=48):
    """生成三次贝塞尔曲线的采样点（用于把直线梯形改成柔和的丝带）。

    输入示例: bezier_points((0, 0), (0, 1), (1, 0), (1, 1))
    输出示例: [(x, y), ...]（samples+1 个点，从 start 到 end）
    """
    points = []
    for index in range(samples + 1):
        t = index / samples
        x = (
            ((1 - t) ** 3) * start[0]
            + 3 * ((1 - t) ** 2) * t * control1[0]
            + 3 * (1 - t) * (t**2) * control2[0]
            + (t**3) * end[0]
        )
        y = (
            ((1 - t) ** 3) * start[1]
            + 3 * ((1 - t) ** 2) * t * control1[1]
            + 3 * (1 - t) * (t**2) * control2[1]
            + (t**3) * end[1]
        )
        points.append((x, y))
    return points


def build_ribbon(reference_left, reference_right, query_left, query_right, reference_y, query_y):
    """用两条贝塞尔曲线拼出一条丝带的顶点序列（上沿 + 下沿）。

    输入示例: build_ribbon(0.1, 0.3, 0.15, 0.5, 3.0, 1.0)
    输出示例: [(x, y), ...]  可直接交给 matplotlib 的 Polygon
    说明：控制点取垂直方向 45% 处（常用 sankey/丝带画法）：连接处切线竖直、中段平滑；
         参考行在下方（reference_y < query_y）时曲线方向自动反过来，仍然平滑。
    """
    vertical_span = query_y - reference_y
    control_shift = 0.45 * vertical_span
    upper_edge = bezier_points(
        (reference_left, reference_y),
        (reference_left, reference_y + control_shift),
        (query_left, query_y - control_shift),
        (query_left, query_y),
    )
    lower_edge = bezier_points(
        (query_right, query_y),
        (query_right, query_y - control_shift),
        (reference_right, reference_y + control_shift),
        (reference_right, reference_y),
    )
    return upper_edge + lower_edge


def format_block_label(block, label_coordinate="query"):
    """把一个比对块格式化成图上的坐标标签（始终带 (+/-) 方向）。

    输入示例: block={"q_start": 0, "q_end": 50, "r_start": 41755, "r_end": 41805,
                     "strand": "-"}, label_coordinate="reference"
    输出示例: "41755-41805 (-)"
    说明：
      * label_coordinate="query"     -> `查询起点-查询终点 (方向)`，默认，保持旧行为；
      * label_coordinate="reference" -> `参考起点-参考终点 (方向)`，适合“长参考 + 短查询”；
      * label_coordinate="both"      -> `查询坐标 (方向) @ 参考坐标`。
    """
    query_text = "{0}-{1}".format(block["q_start"], block["q_end"])
    reference_text = "{0}-{1}".format(block["r_start"], block["r_end"])
    if label_coordinate == "reference":
        return "{0} ({1})".format(reference_text, block["strand"])
    if label_coordinate == "both":
        return "{0} ({1}) @ {2}".format(query_text, block["strand"], reference_text)
    return "{0} ({1})".format(query_text, block["strand"])


def draw_ribbon_plot(
    row_order,
    length_of,
    blocks,
    output_path,
    title,
    dpi=100,
    max_labels_per_row=8,
    min_label_bp=300,
    label_coordinate="query",
):
    """按给定行顺序画丝带图（参考行可以不在第一行）。

    输入示例: draw_ribbon_plot(["reference", "assembly_a"], {"reference": 32000, "assembly_a": 23200}, blocks,
                              "out.png", "标题", 150, 8)
    输出示例: 生成 png（按 dpi）与同名 pdf
    说明：
      * block_bp < min_label_bp 的块不标坐标（默认 300；画 50 bp 侧翼时可设 20）；
      * 每个标签写成 `start-end (+/-)`，正负方向始终显式标出；
      * label_coordinate 决定坐标是 query（查询片段自身）/ reference（参考行，即组装 contig）/ both；
      * 一行上的标签超过 max_labels_per_row 时，只标首尾两个并加"N blocks"说明，
        避免串联多拷贝行的坐标文字互相重叠（精确坐标在 .blocks.tsv / .regions.tsv 里）。
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon, Rectangle

    query_names = [name for name in row_order if name != REFERENCE_TOKEN]
    color_of = {
        name: COLOR_CYCLE[index % len(COLOR_CYCLE)] for index, name in enumerate(query_names)
    }
    row_count = len(row_order)
    y_of = {name: (row_count - 1 - index) for index, name in enumerate(row_order)}  # 自上而下
    reference_y = y_of[REFERENCE_TOKEN]
    reference_length = length_of[REFERENCE_TOKEN]
    reference_name = length_of.get("__reference_name__", REFERENCE_TOKEN)
    figure, axes = plt.subplots(figsize=(14, 1.4 + 0.75 * row_count))
    left, right = 0.18, 0.97

    def x_of(position, row_length):
        """某一行上的坐标 → 图坐标（每行按自身长度缩放到 left..right）。"""
        if row_length <= 0:
            return left
        return left + (position / row_length) * (right - left)

    for name in row_order:
        row_y = y_of[name]
        row_length = length_of[name]
        if name == REFERENCE_TOKEN:
            axes.plot(
                [x_of(0, row_length), x_of(row_length, row_length)],
                [row_y, row_y],
                color="#333333",
                linewidth=2.2,
                zorder=4,
            )
            axes.text(
                left - 0.012,
                row_y,
                "{0}\n({1:,} bp)".format(reference_name, row_length),
                ha="right",
                va="center",
                fontsize=10,
                fontweight="bold",
            )
        else:
            axes.plot(
                [x_of(0, row_length), x_of(row_length, row_length)],
                [row_y, row_y],
                color="#BBBBBB",
                linewidth=1.4,
                zorder=2,
            )
            axes.text(
                left - 0.012,
                row_y,
                "{0}\n({1:,} bp)".format(name, row_length),
                ha="right",
                va="center",
                fontsize=9,
                color=color_of[name],
            )

    # 第一遍：画丝带与两侧矩形，同时收集"要标的标签"（同一位置相差 <=20 bp 只留一个）
    labels_of_row = {name: [] for name in query_names}
    seen_labels = set()
    for block in sorted(blocks, key=lambda item: (y_of.get(item["record"], 0), item["q_start"])):
        row_name = block["record"]
        if row_name not in y_of or row_name == REFERENCE_TOKEN:
            continue
        row_y = y_of[row_name]
        row_length = length_of[row_name]
        query_left = x_of(block["q_start"], row_length)
        query_right = x_of(block["q_end"], row_length)
        reference_left = x_of(block["r_start"], reference_length)
        reference_right = x_of(block["r_end"], reference_length)
        color = color_of[row_name]

        axes.add_patch(
            Polygon(
                build_ribbon(
                    reference_left, reference_right, query_left, query_right, reference_y, row_y
                ),
                closed=True,
                facecolor=color,
                alpha=0.22,
                edgecolor="none",
                zorder=1,
            )
        )
        axes.add_patch(
            Rectangle(
                (query_left, row_y - 0.075),
                max(query_right - query_left, 0.0004),
                0.15,
                facecolor=color,
                edgecolor="none",
                zorder=4,
            )
        )
        axes.add_patch(
            Rectangle(
                (reference_left, reference_y - 0.075),
                max(reference_right - reference_left, 0.0004),
                0.15,
                facecolor=color,
                edgecolor="none",
                zorder=5,
            )
        )
        if block["block_bp"] < min_label_bp:
            continue
        label_text = format_block_label(block, label_coordinate)
        seen_key = (
            row_name,
            round(block["q_start"] / 20),
            round(block["q_end"] / 20),
            block["strand"],
        )
        if seen_key in seen_labels:
            continue
        seen_labels.add(seen_key)
        labels_of_row[row_name].append(
            {"text": label_text, "middle": (query_left + query_right) / 2, "color": color}
        )

    # 第二遍：按策略画标签（块太多时只标首尾 + 说明，避免文字互相重叠）
    label_offsets = {}
    for row_name, labels in labels_of_row.items():
        if not labels:
            continue
        row_y = y_of[row_name]
        away = 1 if row_y > reference_y else -1  # 标签放在"远离参考行"的一侧
        chosen = labels if len(labels) <= max_labels_per_row else [labels[0], labels[-1]]
        for level, label in enumerate(chosen):
            axes.text(
                label["middle"],
                row_y + away * (0.2 + 0.17 * level),
                label["text"],
                ha="center",
                va="bottom" if away > 0 else "top",
                fontsize=7,
                color=label["color"],
            )
            label_offsets.setdefault(row_name, []).append(label["middle"])
        if len(labels) > max_labels_per_row:
            all_left = min(label["middle"] for label in labels)
            all_right = max(label["middle"] for label in labels)
            axes.text(
                (all_left + all_right) / 2,
                row_y + away * (0.2 + 0.17 * len(chosen)),
                "{0} blocks (exact coords in .regions.tsv)".format(len(labels)),
                ha="center",
                va="bottom" if away > 0 else "top",
                fontsize=7,
                color=color_of[row_name],
            )

    copy_count = collections.Counter(block["record"] for block in blocks)
    handles = [
        plt.Line2D(
            [],
            [],
            color=color_of[name],
            linewidth=6,
            label="{0}{1}".format(
                name, " ×{0}".format(copy_count[name]) if copy_count[name] > 1 else ""
            ),
        )
        for name in query_names
    ]
    if handles:
        axes.legend(
            handles=handles,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.05 - 0.02 * row_count),
            frameon=False,
            fontsize=8,
            ncol=min(5, len(query_names)),
        )
    extra_space = 0.22 * max([len(values) for values in label_offsets.values()] or [0])
    axes.set_xlim(0, 1)
    axes.set_ylim(-0.8 - extra_space, (row_count - 1) + 0.8 + extra_space)
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_title(title, fontsize=11)
    for spine in ("top", "right", "left", "bottom"):
        axes.spines[spine].set_visible(False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi)
    figure.savefig(output_path.rsplit(".", 1)[0] + ".pdf")
    plt.close(figure)


### 同源比对区域表 ###

REGION_TABLE_COLUMNS = (
    ("block_id", "区域内顺序编号：按参考坐标升序（同一参考区域再按查询坐标升序）"),
    ("ref_name", "参考序列名（坐标系的基准）"),
    ("ref_length_bp", "参考序列总长度（用于理解 ref_coverage）"),
    ("ref_start", "该同源区域在参考序列上的起点（0-based 半开）"),
    ("ref_end", "该同源区域在参考序列上的终点（0-based 半开）"),
    ("ref_region_bp", "该区域在参考序列上占的长度 = ref_end - ref_start"),
    (
        "gap_from_previous_ref_region_bp",
        "与**上一行参考区域**的间隔（bp，取距离所以 >=0）：0 = 紧邻或相互重叠"
        "（例如同一参考区域的多份拷贝）；第一行没有上一行，填占位符 -",
    ),
    ("ref_coverage", "该区域覆盖整条参考序列的比例（1.00 = 整段覆盖）"),
    ("query_order", "该查询片段在图里的行序（1 = 图最上面一行；受 --order/--order-file 影响）"),
    ("query_name", "落在该参考区域的查询片段名（图里的行名）"),
    ("query_length_bp", "该查询片段的总长度（用于理解 query_coverage）"),
    ("query_start", "同一区域在该查询片段**自身坐标**上的起点（0-based 半开）"),
    ("query_end", "同一区域在该查询片段**自身坐标**上的终点（0-based 半开）"),
    ("query_region_bp", "该区域在查询片段上占的长度 = query_end - query_start"),
    ("query_coverage", "该区域覆盖整条查询片段的比例（1.00 = 整条都落在参考上）"),
    ("strand", "查询相对参考的方向（+ 同向 / - 反向）"),
    ("identity_percent", "相似度百分比 = 匹配碱基数 / 比对长度（如 99.87%）"),
    ("ref_region", "人类可读的参考区间 `ref_name:ref_start-ref_end`（可直接喂给 samtools faidx）"),
    ("query_region", "人类可读的查询区间 `query_name:query_start-query_end`"),
)


def build_region_rows(reference_name, reference_length, row_lengths, blocks, query_order=None):
    """把比对块整理成"以参考为坐标系"的区域表行。

    输入示例: build_region_rows("ref", 12000, {"asmA": 23200}, blocks, {"asmA": 1})
    输出示例: 一行形如
        [1, "ref", 12000, 5000, 17000, 12000, "-", "1.0000", 1, "asmA", 23200, 8000, 20000, 12000,
         "0.5172", "+", "99.87%", "ref:5000-17000", "asmA:8000-20000"]
    说明：排序为"参考坐标升序 → 查询坐标升序"，所以相邻两行就是参考上相邻的两个区域；
         gap_from_previous_ref_region_bp 用上一行的 ref_end 算；第一行没有前一行，填占位符 -。
    """
    order_of_query = dict(query_order or {})
    ordered = sorted(
        blocks, key=lambda block: (block["r_start"], block["q_start"], block["record"])
    )
    rows = []
    previous_ref_end = None
    for index, block in enumerate(ordered, 1):
        query_length = row_lengths.get(block["record"], 0)
        ref_span = block["r_end"] - block["r_start"]
        query_span = block["q_end"] - block["q_start"]
        # 间隔取"距离"，最小 0：紧邻或相互重叠（例如同一参考区域的多份拷贝）都记 0
        gap = "-" if previous_ref_end is None else str(max(0, block["r_start"] - previous_ref_end))
        rows.append(
            [
                index,
                reference_name,
                reference_length,
                block["r_start"],
                block["r_end"],
                ref_span,
                gap,
                "{0:.4f}".format(ref_span / reference_length) if reference_length else "-",
                order_of_query.get(block["record"], "-"),
                block["record"],
                query_length,
                block["q_start"],
                block["q_end"],
                query_span,
                "{0:.4f}".format(query_span / query_length) if query_length else "-",
                block["strand"],
                "{0:.2f}%".format(block["identity"] * 100),
                "{0}:{1}-{2}".format(reference_name, block["r_start"], block["r_end"]),
                "{0}:{1}-{2}".format(block["record"], block["q_start"], block["q_end"]),
            ]
        )
        previous_ref_end = (
            max(previous_ref_end, block["r_end"])
            if previous_ref_end is not None
            else block["r_end"]
        )
    return rows


def merge_intervals(intervals):
    """把 (起, 止) 区间合并成不重叠区间（用来算"覆盖了多少 bp"，重叠只算一次）。

    输入示例: [(0, 10), (5, 20), (30, 40)]  -> 输出示例: [(0, 20), (30, 40)]
    """
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], end)
        else:
            merged.append([start, end])
    return [(start, end) for start, end in merged]


def build_placement_rows(
    reference_name, reference_length, row_lengths, blocks, query_order=None, full_tolerance_bp=10
):
    """**每个比对位置一行**：这段序列在参考上的某一处落在哪、方向、多相似、覆盖多少。

    输入示例: build_placement_rows("asm", 32000, {"insert_alpha": 1200}, blocks, {"insert_alpha": 2})
    输出示例: [["insert_alpha", 1200, 1, 3, "asm", 32000, 12000, 13200, 1200, 0, 1200, 1200, "+", "100.00%", "yes"],
               ["insert_alpha", 1200, 2, 3, …]]（下略）
    说明：**不做汇总** —— 一行 = 一个比对块，方向和相似度逐处给出（不再出现"混合"这种笼统值）；
          `placement_index` / `n_placements` = 这是第几处 / 这段共几处（按参考坐标升序编号）。
          `query_start/query_end` 是这一处在**该片段自身坐标**上的区间（配合 query_length_bp 就能看出是否整段命中）；
          `covers_full_query` = yes/no（该处查询侧跨度 ≥ 片段长度 − full_tolerance_bp）。
          行序：先按图中的行序（给了 query_order 就按它），同一片段内按参考坐标升序。
    """
    order_of_query = dict(query_order or {})
    grouped = collections.defaultdict(list)
    for block in blocks:
        grouped[block["record"]].append(block)
    names = sorted(grouped, key=lambda name: (order_of_query.get(name, 10**6), name))
    rows = []
    for name in names:
        query_length = row_lengths.get(name, 0)
        hits = sorted(grouped[name], key=lambda block: (block["r_start"], block["q_start"]))
        for index, block in enumerate(hits, 1):
            query_span = block["q_end"] - block["q_start"]
            rows.append(
                [
                    name,
                    query_length,
                    index,
                    len(hits),
                    reference_name,
                    reference_length,
                    block["r_start"],
                    block["r_end"],
                    block["r_end"] - block["r_start"],
                    block["q_start"],
                    block["q_end"],
                    query_span,
                    block["strand"],
                    "{0:.2f}%".format(block["identity"] * 100),
                    "yes" if query_span >= query_length - full_tolerance_bp else "no",
                ]
            )
    return rows


PLACEMENT_COLUMNS = (
    ("query_name", "查询片段名（每行一段：外源序列或参考侧片段）"),
    ("query_length_bp", "该片段自身长度（bp）"),
    ("placement_index", "这是该片段在参考上的**第几处**（按参考坐标升序编号）"),
    ("n_placements", "该片段在参考上**共几处**"),
    ("reference_name", "参考序列名（本工具的参考行；在「单条组装」用法里就是那条组装片段）"),
    ("reference_length_bp", "参考序列总长度（bp）"),
    ("ref_start", "这一处在参考上的起点（0-based 半开）"),
    ("ref_end", "这一处在参考上的终点（0-based 半开）"),
    ("ref_length_bp", "这一处在参考上占的长度 = ref_end − ref_start"),
    ("query_start", "这一处在**该片段自身坐标**上的起点（0-based 半开）"),
    ("query_end", "这一处在**该片段自身坐标**上的终点（0-based 半开）"),
    ("query_length_in_hit_bp", "这一处覆盖该片段的长度 = query_end − query_start"),
    ("strand", "这一处的方向：+ 同向 / - 反向（**逐处给出**，不汇总）"),
    ("identity_percent", "这一处的相似度（百分比）"),
    ("covers_full_query", "这一处是否整段命中：yes = 查询侧跨度 ≥ 片段长度 − 10 bp"),
)


def write_placement_table(
    path, reference_name, reference_length, row_lengths, blocks, query_order=None
):
    """写"每行一段序列在参考上的位置汇总"表（表头前带目标与逐列说明）。

    输入示例: write_placement_table("out.placements.tsv", "asm", 32000, {"insert_alpha": 1200}, blocks,
                                    {"insert_alpha": 2})
    输出示例: 表头前带 # 目标 / # 列说明 注释块；**一行一个比对位置**；返回数据行数
    """
    rows = build_placement_rows(reference_name, reference_length, row_lengths, blocks, query_order)
    with open(path, "w") as table_file:
        table_file.write(
            "# 目标：逐处列出「每一行（查询片段）在参考上的位置」——**一行一个比对位置**，\n"
        )
        table_file.write(
            "#       给出起止坐标、方向、相似度，以及「这是第几处 / 共几处、是否整段命中」\n"
        )
        table_file.write(
            "# 坐标系：全部 0-based 半开；ref_start/ref_end = 参考（本图参考行）坐标，\n"
        )
        table_file.write("#         query_start/query_end = 该片段自身坐标\n")
        table_file.write(
            "# 与 .regions.tsv 的分工：区域表以参考为主、按参考区域排列（含与上一区域的间隔）；\n"
        )
        table_file.write("#        本表以查询片段为主、按片段分组逐处列出（含第几处/共几处）\n")
        table_file.write("# 参考序列：{0}（{1:,} bp）\n".format(reference_name, reference_length))
        table_file.write("# 列说明（按表头顺序）：\n")
        for name, description in PLACEMENT_COLUMNS:
            table_file.write("#   {0}：{1}\n".format(name, description))
        table_file.write("\t".join(name for name, _ in PLACEMENT_COLUMNS) + "\n")
        for row in rows:
            table_file.write("\t".join(str(value) for value in row) + "\n")
    return len(rows)


def write_region_table(
    path, reference_name, reference_length, row_lengths, blocks, query_order=None, display_name=None
):
    """写"同源比对区域表"：以参考序列为坐标参考系，一行一个比对区域；表内自带逐列说明。

    输入示例: write_region_table("out.regions.tsv", "ref", 12000, {"asmA": 23200}, blocks, {"asmA": 1})
    输出示例: 表头前带 # 目标 / # 坐标系 / # 列说明 注释块，随后是表头与数据行；返回行数
    """
    rows = build_region_rows(reference_name, reference_length, row_lengths, blocks, query_order)
    with open(path, "w") as table_file:
        table_file.write("# 目标：输出同源比对区域坐标信息（以参考序列为坐标参考系）\n")
        table_file.write(
            "# 参考序列：{0}（{1:,} bp）\n".format(display_name or reference_name, reference_length)
        )
        table_file.write(
            "# 坐标系：全部 0-based 半开；ref_* 为参考序列坐标，query_* 为查询片段自身坐标\n"
        )
        table_file.write(
            "#        若某个查询片段本身就取自另一套坐标（例如某段基因组窗口），"
            "它在原坐标里的位置 = 该片段起点 + query_start/query_end；本表不给那套坐标\n"
        )
        table_file.write(
            "# 排序：先按参考坐标升序，同一参考区域再按查询坐标升序"
            "（同一段参考上的多个片段会相邻，便于看多拷贝/多片段命中）\n"
        )
        table_file.write("# 列说明（按表头顺序）：\n")
        for name, description in REGION_TABLE_COLUMNS:
            table_file.write("#   {0}：{1}\n".format(name, description))
        table_file.write("\t".join(name for name, _ in REGION_TABLE_COLUMNS) + "\n")
        for row in rows:
            table_file.write("\t".join(str(value) for value in row) + "\n")
    return len(rows)


### 命令行 ###


@dataclass(slots=True)
class RibbonSettings:
    """``alnkit plot ribbon`` 的全部入参（与旧脚本参数一一对应）。"""

    outdir: str = "."
    name: str = "ribbon"
    reference: str | None = None
    query: tuple[str, ...] = ()
    query_file: str | None = None
    paf: str | None = None
    reference_name: str | None = None
    query_name: tuple[str, ...] = ()
    title: str = ""
    reference_position: str = "top"
    order: str | None = None
    order_file: str | None = None
    dpi: int = 100
    min_label_bp: int = 300
    label_coordinate: str = "query"
    max_labels_per_row: int = 8
    placement_table: str | None = None
    region_table: str | None = None
    min_identity: float = 0.8
    min_block_bp: int = 200
    preset: str = "asm20"
    query_preset: tuple[str, ...] = ()
    threads: int = 8


def collect_rows_from_fasta_mode(settings, outdir_prefix):
    """方式 1（内部比对）：从 --reference/--query 准备参考序列、行长度与每个文件对应的行名。

    输入示例: collect_rows_from_fasta_mode(args, "run")
    输出示例: (参考fasta, 显示名, 记录名, 参考长度, {行名:长度}, {行名:(文件, 记录名, 标签)}, {比对文件:原路径})
    说明：显示名可以带说法（如 `insert_alpha (1-1200)`），表里一律写记录名（insert_alpha），便于与其它表对齐核对。
          标签跟着行名一起返回，因为 --query-preset 是按**标签**（不是记录名/文件名）覆盖比对参数的。
    """
    reference_label, reference_path, reference_contig = parse_spec(settings.reference)
    records = list_records(reference_path)
    if reference_contig:
        reference_record = reference_contig  # 从文件里取哪条序列
    elif len(records) == 1:
        reference_record = records[0][0]
    else:
        raise UsageError(
            "参考侧 {0} 有 {1} 条序列，请用 `标签=路径:片段名` 指定其中一条".format(
                reference_path, len(records)
            )
        )
    # 行上显示的名字：优先用标签（可以写坐标等说明，便于读图），没给标签就用记录名
    reference_name = reference_label or reference_record
    reference_fasta = os.path.join(settings.outdir, "{0}.reference.fasta".format(outdir_prefix))
    reference_length = extract_single_record(reference_path, reference_record, reference_fasta)

    row_lengths = {}
    file_of_row = {}
    source_of_file = {}  # 比对用的路径 → 用户给的路径（多记录文件会被转大写另存，日志里显示原路径）
    for spec in settings.query or []:
        label, fasta_path, contig = parse_spec(spec)
        if contig:
            row_name = label or contig
            row_lengths[row_name] = extract_single_record(
                fasta_path,
                contig,
                os.path.join(
                    settings.outdir,
                    "{0}.{1}.fasta".format(outdir_prefix, row_name.replace("/", "_")),
                ),
            )
            file_of_row[row_name] = (fasta_path, contig, label)
            continue
        # 多记录文件：先整体转大写再比对（本机 minimap2 大小写一致，转大写在旧版本/其它工具上也稳）
        upper_path = os.path.join(
            settings.outdir,
            "{0}.{1}.upper.fasta".format(outdir_prefix, (label or "query").replace("/", "_")),
        )
        write_upper_case_fasta(fasta_path, upper_path)
        query_records = list_records(upper_path)
        for record_name, record_length in query_records:
            if len(query_records) == 1 and label:
                row_name = label
            elif label:
                row_name = "{0}:{1}".format(label, record_name)
            else:
                row_name = record_name
            row_lengths[row_name] = record_length
            file_of_row[row_name] = (upper_path, record_name, label)
            source_of_file[upper_path] = fasta_path
    # reference_name = 行上显示的名字（标签优先）；reference_record = 文件里的记录名（表里用，稳定可核对）
    return (
        reference_fasta,
        reference_name,
        reference_record,
        reference_length,
        row_lengths,
        file_of_row,
        source_of_file,
    )


def collect_rows_from_paf_mode(settings):
    """方式 2（读 PAF）：从 PAF 里准备参考名/行长度，并返回全部比对块。

    输入示例: collect_rows_from_paf_mode(args)  -> 输出示例: ("ref", 38291, {"q": 8000}, [块...])
    """
    blocks = read_paf_file(settings.paf)
    if not blocks:
        raise UsageError("PAF 文件里没有可用的比对块：{0}".format(settings.paf))
    target_names = sorted({block["r_name"] for block in blocks})
    if settings.reference_name:
        if settings.reference_name not in target_names:
            raise UsageError(
                "--reference-name {0} 不在 PAF 的目标里（可用：{1}）".format(
                    settings.reference_name, ", ".join(target_names)
                )
            )
        reference_name = settings.reference_name
    elif len(target_names) == 1:
        reference_name = target_names[0]
    else:
        raise UsageError(
            "PAF 里有多条参考序列（{0}），请用 --reference-name 指定其中一条".format(
                ", ".join(target_names)
            )
        )

    row_lengths = {}
    for block in blocks:
        if block["r_name"] == reference_name:
            row_lengths[block["record"]] = block["q_length"]
    reference_length = next(
        block["r_length"] for block in blocks if block["r_name"] == reference_name
    )
    kept = [block for block in blocks if block["r_name"] == reference_name]
    return reference_name, reference_length, row_lengths, kept


def run(settings: RibbonSettings) -> dict[str, object]:
    os.makedirs(settings.outdir, exist_ok=True)

    if settings.paf:
        if settings.reference or settings.query:
            print("提示：--paf 模式下忽略 --reference/--query（片段名与长度都从 PAF 里读）")
        if settings.query_preset or settings.preset == "auto":
            print("提示：--paf 模式下 --preset/--query-preset 不生效（比对已由 PAF 决定）")
        reference_name, reference_length, row_lengths, blocks = collect_rows_from_paf_mode(settings)
        reference_record = reference_name
    else:
        if not settings.reference or not settings.query:
            raise UsageError(
                "方式 1 需要同时给 --reference 与 --query；或改用 --paf 读现成比对结果"
            )
        (
            reference_fasta,
            reference_name,
            reference_record,
            reference_length,
            row_lengths,
            file_of_row,
            source_of_file,
        ) = collect_rows_from_fasta_mode(settings, settings.name)
        rows_by_file = {}
        label_of_row = {}
        for row_name, (fasta_path, _record_name, row_label) in file_of_row.items():
            rows_by_file.setdefault(fasta_path, []).append(row_name)
            label_of_row[row_name] = row_label
        # --query-preset 先校验标签是否存在（写错了不要等到比对完才发现）
        overrides = collect_preset_overrides(settings.query_preset, set(label_of_row.values()))
        blocks = []
        for fasta_path, row_names in rows_by_file.items():
            file_lengths = [row_lengths[name] for name in row_names]
            display_path = source_of_file.get(fasta_path, fasta_path)
            preset_args, preset_source = resolve_file_preset(
                fasta_path,
                [label_of_row[name] for name in row_names],
                file_lengths,
                settings.preset,
                overrides,
                display_path,
            )
            print(
                "比对参数：{0}（{1} 条记录，{2:,}–{3:,} bp）→ {4}    [{5}]".format(
                    display_path,
                    len(row_names),
                    min(file_lengths),
                    max(file_lengths),
                    " ".join(preset_args),
                    preset_source,
                )
            )
            for block in align_file_to_reference(
                reference_fasta, fasta_path, settings.threads, preset_args
            ):
                record_name = block["record"]
                matched = next(
                    (
                        name
                        for name in row_names
                        if name == record_name
                        or name.endswith(":" + record_name)
                        or file_of_row[name][1] == record_name
                    ),
                    None,
                )
                if matched is None:
                    continue
                block["record"] = matched
                blocks.append(block)

    # 只画指定的查询片段（PAF 模式常用来挑选）
    if settings.query_name:
        wanted = set(settings.query_name)
        missing = wanted - set(row_lengths)
        if missing:
            print(
                "告警：--query-name 里这些名字不在比对结果中，已忽略：{0}".format(
                    ", ".join(sorted(missing))
                )
            )
        row_lengths = {name: length for name, length in row_lengths.items() if name in wanted}
        blocks = [block for block in blocks if block["record"] in wanted]

    # 行顺序：--order + --order-file（列表文件里每行一个片段名）
    order_file_names = read_order_file(settings.order_file) if settings.order_file else ()
    row_order, order_warnings = resolve_row_order(
        list(row_lengths), settings.order, settings.reference_position, order_file_names
    )
    for message in order_warnings:
        print("告警：{0}".format(message))

    kept = [
        block
        for block in blocks
        if block["identity"] >= settings.min_identity and block["block_bp"] >= settings.min_block_bp
    ]
    if not kept:
        # 短片段常碰到这个：块本身没问题，只是默认阈值把它们滤掉了（不给提示就是一张空图）
        print(
            "告警：没有可画的映射块（比对到 {0} 个，全被 --min-identity {1} / --min-block-bp {2} 过滤，"
            "或本来就没命中）；块很短时可把 --min-block-bp 调小，如 --min-block-bp 20".format(
                len(blocks), settings.min_identity, settings.min_block_bp
            )
        )
    length_of = dict(row_lengths)
    length_of[REFERENCE_TOKEN] = reference_length
    length_of["__reference_name__"] = reference_name
    output_path = os.path.join(settings.outdir, "{0}.png".format(settings.name))
    title = settings.title or "{0} vs {1} query fragments".format(reference_name, len(row_lengths))
    draw_ribbon_plot(
        row_order,
        length_of,
        kept,
        output_path,
        title,
        settings.dpi,
        settings.max_labels_per_row,
        settings.min_label_bp,
        settings.label_coordinate,
    )

    table_path = os.path.join(settings.outdir, "{0}.blocks.tsv".format(settings.name))
    with open(table_path, "w") as table_file:
        table_file.write(
            "query_fragment\tquery_start\tquery_end\tstrand\treference\tref_start\tref_end"
            "\tblock_bp\tidentity\n"
        )
        for block in sorted(kept, key=lambda item: (item["record"], item["q_start"])):
            table_file.write(
                "\t".join(
                    [
                        block["record"],
                        str(block["q_start"]),
                        str(block["q_end"]),
                        block["strand"],
                        reference_record,
                        str(block["r_start"]),
                        str(block["r_end"]),
                        str(block["block_bp"]),
                        "{0:.4f}".format(block["identity"]),
                    ]
                )
                + "\n"
            )
    if settings.placement_table:
        placement_rows = write_placement_table(
            settings.placement_table,
            reference_record,
            reference_length,
            row_lengths,
            kept,
            {name: index for index, name in enumerate(row_order, 1)},
        )
        print("位置汇总表 {0}（{1} 行）".format(settings.placement_table, placement_rows))
    if settings.region_table:
        region_rows = write_region_table(
            settings.region_table,
            reference_record,
            reference_length,
            row_lengths,
            kept,
            {name: index for index, name in enumerate(row_order, 1)},
            display_name=reference_name,
        )
        print("同源比对区域表 {0}（{1} 个区域）".format(settings.region_table, region_rows))
    print(
        "图 {0}（参考 {1} {2:,} bp，查询 {3} 个片段，保留 {4} 个映射块）".format(
            output_path, reference_name, reference_length, len(row_lengths), len(kept)
        )
    )
    print("行顺序（自上而下）：{0}".format(" > ".join(row_order)))

    return {
        "png": output_path,
        "blocks_tsv": table_path,
        "reference": reference_name,
        "reference_length": reference_length,
        "queries": len(row_lengths),
        "blocks": len(kept),
        "row_order": row_order,
    }
