"""CLI 骨架：命令树、退出码、帮助信息。"""

from __future__ import annotations

import pytest

from alnkit import __version__
from alnkit.cli import build_parser, main


def test_parser_exposes_commands() -> None:
    parser = build_parser()
    actions = [a for a in parser._subparsers._group_actions if hasattr(a, "choices")]  # noqa: SLF001
    choices = set(actions[0].choices)
    assert choices == {"doctor", "selftest", "seg", "plot", "short", "long"}


def test_no_arguments_prints_help(capsys) -> None:
    assert main([]) == 0
    out = capsys.readouterr().out
    assert "alnkit" in out
    assert "doctor" in out


def test_version(capsys) -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(["--version"])
    assert excinfo.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_unknown_command_exits_2() -> None:
    with pytest.raises(SystemExit) as excinfo:
        main(["nope"])
    assert excinfo.value.code == 2


def test_selftest_passes(capsys) -> None:
    assert main(["selftest"]) == 0
    out = capsys.readouterr().out
    assert "通过" in out
    assert "FAIL" not in out


def test_doctor_runs(capsys) -> None:
    assert main(["doctor"]) == 0
    out = capsys.readouterr().out
    assert "alnkit" in out
    assert "pysam" in out


def test_global_options_are_accepted_after_subcommand(capsys) -> None:
    assert main(["selftest", "--log-level", "ERROR", "--quiet"]) == 0


def test_removed_assemble_flags_are_rejected(capsys) -> None:
    """0.4.0 评估后删除的参数不能再被接受（否则用户的旧脚本会静默失效）。"""
    for flag in ("--min-junction-depth-ratio", "--recruit-edges", "--min-bridge-reads"):
        with pytest.raises(SystemExit) as excinfo:
            main(["long", "assemble", "--classes-dir", "x", "--paf", "y", "--outdir", "z", flag])
        assert excinfo.value.code == 2
