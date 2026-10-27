"""Tool runs are isolated from ``~/.local`` Python packages (docs/specs/fix-pythonnousersite.md).

The Python 3.10 tool envs imported ViennaRNA/biopython from the user site because nothing set
``PYTHONNOUSERSITE``; OSTIR results differed on 6 of ~9,000 genes.
"""

from __future__ import annotations

import os
import re
import subprocess
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import pytest

from gippa.setupdb import run as sr

ROOT = Path(__file__).resolve().parents[2]
NF_CONFIG = ROOT / "workflow" / "nextflow.config"
FARMSHARE_ENV = ROOT / "scripts" / "farmshare_env.sh"
OSTIR_MAIN = ROOT / "workflow" / "modules" / "local" / "ostir" / "main.nf"


def _strip_comments(text: str) -> str:
    return re.sub(r"(^|\s)//[^\n]*", r"\1", text)  # not the // of docker:// URIs


def _block(text: str, open_at: int) -> str:
    """The text between the ``{`` at ``open_at`` and its matching ``}`` (braces balanced)."""
    assert text[open_at] == "{"
    depth = 0
    for i in range(open_at, len(text)):
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        if depth == 0:
            return text[open_at + 1 : i]
    raise AssertionError(f"unbalanced braces from offset {open_at}")


def _top_level_blocks(text: str, name: str) -> list[str]:
    """Bodies of every ``name { ... }`` scope at brace depth 0."""
    out, depth, i = [], 0, 0
    pat = re.compile(rf"\b{name}\s*\{{")
    while i < len(text):
        if depth == 0 and (m := pat.match(text, i)) and (i == 0 or not text[i - 1].isalnum()):
            body = _block(text, m.end() - 1)
            out.append(body)
            i = m.end() + len(body) + 1
            continue
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        i += 1
    return out


def _config() -> str:
    return _strip_comments(NF_CONFIG.read_text(encoding="utf-8"))


# ---------------------------------------------------------------- nextflow.config


def test_nextflow_config_sets_pythonnousersite_in_a_top_level_env_scope() -> None:
    envs = _top_level_blocks(_config(), "env")
    assert envs, "workflow/nextflow.config has no top-level env { } scope"
    assert any(re.search(r"\bPYTHONNOUSERSITE\s*=\s*['\"]1['\"]", b) for b in envs), envs


def test_env_scope_is_not_inside_a_profile() -> None:
    text = _config()
    profiles = _top_level_blocks(text, "profiles")
    assert len(profiles) == 1
    # the top-level env scope must exist outside profiles (applies under every profile)
    outside = text.replace(profiles[0], "")
    assert any("PYTHONNOUSERSITE" in b for b in _top_level_blocks(outside, "env"))


def _profile(name: str) -> str:
    profiles = _top_level_blocks(_config(), "profiles")[0]
    m = re.search(rf"(?m)^\s*{name}\s*\{{", profiles)
    assert m, f"no {name} profile"
    return _block(profiles, m.end() - 1)


def test_docker_env_whitelist_passes_pythonnousersite() -> None:
    m = re.search(r"docker\.envWhitelist\s*=\s*['\"]([^'\"]*)['\"]", _profile("docker"))
    assert m, "docker profile has no docker.envWhitelist"
    names = {n.strip() for n in m.group(1).split(",")}
    assert "PYTHONNOUSERSITE" in names, names
    assert "TMPDIR" in names, names  # the existing per-task TMPDIR stays whitelisted


# ---------------------------------------------------------------- ToolRunner.run


def _capture(monkeypatch: pytest.MonkeyPatch) -> list[dict[str, Any]]:
    calls: list[dict[str, Any]] = []

    def fake_run(cmd: Sequence[str], *args: Any, **kwargs: Any) -> subprocess.CompletedProcess:
        calls.append({"cmd": list(cmd), **kwargs})
        return subprocess.CompletedProcess(list(cmd), 7)

    monkeypatch.setattr(sr.subprocess, "run", fake_run)
    monkeypatch.setattr(sr.ToolRunner, "command", lambda self, tool, argv: ["fake", tool, *argv])
    return calls


def test_tool_runner_sets_pythonnousersite(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("PYTHONNOUSERSITE", raising=False)
    calls = _capture(monkeypatch)
    assert sr.ToolRunner().run("ostir", ["--version"]) == 7
    assert len(calls) == 1 and calls[0]["cmd"] == ["fake", "ostir", "--version"]
    env = calls[0].get("env")
    assert env is not None, "ToolRunner.run passes no env= (the child inherits the user site)"
    assert env["PYTHONNOUSERSITE"] == "1"
    assert "PYTHONNOUSERSITE" not in os.environ  # os.environ itself is untouched


def test_tool_runner_env_inherits_os_environ(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GIPPA_TEST_SENTINEL", "keep-me")
    calls = _capture(monkeypatch)
    sr.ToolRunner().run("ostir", [])
    env = calls[0].get("env")
    assert env is not None
    assert env["GIPPA_TEST_SENTINEL"] == "keep-me"
    assert env.get("PATH") == os.environ.get("PATH")
    assert {k: v for k, v in env.items() if k != "PYTHONNOUSERSITE"} == {
        k: v for k, v in os.environ.items() if k != "PYTHONNOUSERSITE"
    }


def test_tool_runner_overrides_pythonnousersite_zero(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("PYTHONNOUSERSITE", "0")
    calls = _capture(monkeypatch)
    sr.ToolRunner().run("ostir", [])
    env = calls[0].get("env")
    assert env is not None and env["PYTHONNOUSERSITE"] == "1"
    assert os.environ["PYTHONNOUSERSITE"] == "0"


# ---------------------------------------------------------------- farmshare_env.sh


def test_farmshare_env_exports_pythonnousersite_statically() -> None:
    text = FARMSHARE_ENV.read_text(encoding="utf-8")
    assert re.search(r"(?m)^\s*export\s+PYTHONNOUSERSITE=['\"]?1['\"]?\s*(#.*)?$", text), (
        "scripts/farmshare_env.sh does not export PYTHONNOUSERSITE=1"
    )


def test_sourcing_farmshare_env_exports_pythonnousersite(tmp_path: Path) -> None:
    env = {
        "HOME": str(tmp_path),
        "PATH": "/usr/bin:/bin",
        "USER": "gippa-test",
        "GIPPA_ROOT": str(tmp_path / "local"),
        "PYTHONNOUSERSITE": "0",
    }
    out = subprocess.run(
        ["bash", "--noprofile", "--norc", "-c",
         f'. "{FARMSHARE_ENV}"; env | grep "^PYTHONNOUSERSITE=" || true'],
        env=env, capture_output=True, text=True, timeout=60, check=False,
    )  # fmt: skip
    assert out.returncode == 0, out.stderr
    assert out.stdout.strip() == "PYTHONNOUSERSITE=1", out.stdout


# ---------------------------------------------------------------- OSTIR versions.yml


def _section(text: str, name: str) -> str:
    m = re.search(rf"(?m)^\s*{name}:\s*$", text)
    assert m, f"no {name}: section in {OSTIR_MAIN}"
    rest = text[m.end() :]
    nxt = re.search(r"(?m)^\s*(stub|script|output|input|when|exec):\s*$", rest)
    return rest[: nxt.start()] if nxt else rest


def test_ostir_versions_record_viennarna_binding() -> None:
    script = _section(OSTIR_MAIN.read_text(encoding="utf-8"), "script")
    versions = script[script.index("END_VERSIONS") :]
    m = re.search(r"(?m)^\s*viennarna:\s*(.+)$", versions)
    assert m, "OSTIR versions.yml has no viennarna: line"
    assert re.search(r"import RNA\b", m.group(1)) and "RNA.__version__" in m.group(1), m.group(1)
    assert "python" in m.group(1)


def test_ostir_stub_versions_record_viennarna() -> None:
    stub = _section(OSTIR_MAIN.read_text(encoding="utf-8"), "stub")
    assert re.search(r"(?m)^\s*viennarna:\s*\S+", stub), "OSTIR stub versions.yml has no viennarna:"
