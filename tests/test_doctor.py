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
