"""``alnkit doctor``：必需/可选依赖的判定与提示。"""

from __future__ import annotations

import argparse

import pytest

from alnkit.cli import doctor
from alnkit.core import errors
from alnkit.core.external import Tool


def _fake_tool(name: str, available: bool) -> Tool:
    return Tool(
        name=name,
        path="/usr/bin/" + name if available else None,
        version="9.9" if available else None,
    )


def test_reports_missing_optional_and_aligners(monkeypatch, capsys) -> None:
    monkeypatch.setattr(doctor, "_package_version", lambda name: "1.0")
    monkeypatch.setattr(doctor, "find_tool", lambda name: _fake_tool(name, False))

    assert doctor.run(argparse.Namespace()) == 0
    out = capsys.readouterr().out
    assert "pysam" in out
    assert "matplotlib" in out
    assert "警告" in out  # 比对器一个都没有


def test_all_present_has_no_warning(monkeypatch, capsys) -> None:
    monkeypatch.setattr(doctor, "_package_version", lambda name: "1.0")
    monkeypatch.setattr(doctor, "find_tool", lambda name: _fake_tool(name, True))

    assert doctor.run(argparse.Namespace()) == 0
    assert "警告" not in capsys.readouterr().out


def test_missing_required_dependency_raises(monkeypatch) -> None:
    monkeypatch.setattr(doctor, "_package_version", lambda name: None if name == "pysam" else "1.0")
    monkeypatch.setattr(doctor, "find_tool", lambda name: _fake_tool(name, True))

    with pytest.raises(errors.RunError, match="pysam"):
        doctor.run(argparse.Namespace())


def test_version_probe_falls_back_to_help_output(tmp_path, monkeypatch) -> None:
    """`bwa --version` 不被支持：要退化成"跑一次无参调用，从帮助里抠版本"。"""
    from alnkit.core import external

    fake = tmp_path / "bwa"
    fake.write_text(
        "#!/usr/bin/env bash\n"
        'if [ "${1:-}" = "--version" ]; then\n'
        "  echo \"[main] unrecognized command '--version'\" >&2\n"
        "  exit 1\n"
        "fi\n"
        "echo 'Program: bwa (alignment via Burrows-Wheeler transformation)' >&2\n"
        "echo 'Version: 0.7.19-r1273' >&2\n",
        encoding="utf-8",
    )
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{tmp_path}:{__import__('os').environ['PATH']}")

    tool = external.find_tool("bwa")
    assert tool.available is True
    assert tool.version == "0.7.19-r1273", tool.version


def test_version_probe_handles_plain_version_line() -> None:
    from alnkit.core import external

    assert external._extract_version("2.31-r1302\n") == "2.31-r1302"
    assert external._extract_version("samtools 1.24\nUsing htslib 1.24\n") == "1.24"
    assert external._extract_version("[main] unrecognized command '--version'") is None
