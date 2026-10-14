"""build_envs scripts, Dockerfiles, NOTICE.md and import guard (M2-01; §3.8, R5, §14.6)."""

from __future__ import annotations

import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]
BUILD_SH = REPO / "scripts" / "build_envs.sh"
BUILD_SBATCH = REPO / "scripts" / "slurm" / "build_envs.sbatch"
CHECK_IMPORTS = REPO / "scripts" / "ci" / "check_imports.py"
DOCKER_TOOLS = [
    "uniop",
    "promotech",
    "promoter_calculator",
    "terminatornet",
    "bactermfinder",
    "ostir",
]


@pytest.mark.parametrize("path", [BUILD_SH, BUILD_SBATCH])
def test_scripts_exist_and_parse(path: Path) -> None:
    assert path.is_file(), path
    proc = subprocess.run(["bash", "-n", str(path)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr


def test_sbatch_header() -> None:
    text = BUILD_SBATCH.read_text()
    lines = text.splitlines()
    assert lines[0] == "#!/bin/bash"
    assert re.search(r"^#SBATCH --job-name=giap-\S+$", text, re.M)
    assert re.search(r"^#SBATCH --output=\S*/local/logs/%x\.%j\.out$", text, re.M)
    assert "#SBATCH --partition=normal" in lines
    assert "#SBATCH --cpus-per-task=4" in lines
    assert "#SBATCH --mem=16G" in lines
    assert "#SBATCH --time=08:00:00" in lines
    assert "set -euo pipefail" in lines
    assert "scripts/farmshare_env.sh" in text
    assert 'cd "$GIAP_REPO"' in text


def test_sbatch_builds_envs_and_pulls_images() -> None:
    text = BUILD_SBATCH.read_text()
    assert "build_envs.sh" in text
    assert "apptainer pull" in text
    assert "NXF_APPTAINER_CACHEDIR" in text


def test_build_envs_sh_contract() -> None:
    text = BUILD_SH.read_text()
    assert "micromamba clean -a -y" in text
    assert "--explicit" in text
    assert "workflow/envs/locks" in text
    assert "--override-channels" in text or "nodefaults" in text


def test_build_envs_sh_unknown_tool_exits_2() -> None:
    env = {k: v for k, v in os.environ.items() if k in {"PATH", "HOME", "LANG"}}
    env["GIAP_REPO"] = str(REPO)
    proc = subprocess.run(
        ["bash", str(BUILD_SH), "no_such_tool"],
        capture_output=True,
        text=True,
        cwd=REPO,
        env=env,
        timeout=60,
    )
    assert proc.returncode == 2
    out = proc.stdout + proc.stderr
    for name in ("bakta", "uniop", "terminatornet"):
        assert name in out


def test_resolve_digest_script_compiles() -> None:
    path = REPO / "scripts" / "resolve_digest.py"
    assert path.is_file()
    compile(path.read_text(), str(path), "exec")


@pytest.mark.parametrize("tool", DOCKER_TOOLS)
def test_dockerfiles_exist(tool: str) -> None:
    assert (REPO / "workflow" / "containers" / tool / "Dockerfile").is_file()


@pytest.mark.parametrize("tool", ["uniop", "promotech", "terminatornet", "bactermfinder"])
def test_dockerfile_clones_pinned_commit(tool: str) -> None:
    from giap.annotate.tools import TOOLS

    sha = TOOLS[tool].version.removeprefix("git:")
    assert sha in (REPO / "workflow" / "containers" / tool / "Dockerfile").read_text()


def test_no_tool_source_vendored() -> None:
    """R5/R11: per-tool folders are allowed; tool source or model files inside them are not."""
    source_names = {"terminatornet.py", "promotech.py", "uniop", "uniop.py", "bactermfinder.py"}
    model_suffixes = {".pkl", ".model", ".joblib", ".h5", ".pt", ".sav", ".onnx"}
    for top in ("src", "workflow", "scripts"):
        for path in (REPO / top).rglob("*"):
            if path.is_file():
                assert path.name.lower() not in source_names, path
                assert path.suffix.lower() not in model_suffixes, path
    allowed = {"Dockerfile", "main.nf", "meta.yml", ".dockerignore"}
    for tree in (REPO / "workflow" / "containers", REPO / "workflow" / "modules" / "local"):
        if not tree.is_dir():
            continue
        for path in tree.rglob("*"):
            if path.is_file():
                assert path.name in allowed or path.suffix in {".yml", ".config"}, path


# ---------------------------------------------------------------- NOTICE.md


def _notice_rows() -> dict[str, str]:
    rows = {}
    for line in (REPO / "NOTICE.md").read_text().splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) >= 2 and line.lstrip().startswith("|"):
            key = re.sub(r"[^a-z0-9]", "", cells[0].lower())
            rows[key] = cells[1]
    return rows


def test_notice_lists_every_tool_with_licence() -> None:
    from giap.annotate.tools import TOOLS

    rows = _notice_rows()
    for name, spec in TOOLS.items():
        key = re.sub(r"[^a-z0-9]", "", name.lower())
        assert key in rows, f"NOTICE.md has no row for {name}"
        assert spec.licence.lower() in rows[key].lower(), (name, spec.licence, rows[key])


def test_notice_genomad_non_commercial() -> None:
    line = next(
        ln for ln in (REPO / "NOTICE.md").read_text().splitlines() if ln.startswith("| geNOMAD")
    )
    assert "non-commercial" in line.lower()


# ---------------------------------------------------------------- check_imports.py


def _run_check(root: Path) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(CHECK_IMPORTS), str(root)], capture_output=True, text=True
    )


@pytest.mark.parametrize(
    "source", ["import terminatornet\n", "from bactermfinder.model import predict\n"]
)
def test_check_imports_forbids_separate_process_tools(tmp_path: Path, source: str) -> None:
    (tmp_path / "pkg").mkdir()
    (tmp_path / "pkg" / "mod.py").write_text(source)
    proc = _run_check(tmp_path)
    assert proc.returncode == 1
    assert "mod.py:1" in proc.stdout


def test_check_imports_forbidden_set() -> None:
    import importlib.util

    spec = importlib.util.spec_from_file_location("check_imports_m201", CHECK_IMPORTS)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    assert {"terminatornet", "bactermfinder"} <= mod.FORBIDDEN


def test_check_imports_passes_on_src() -> None:
    assert _run_check(REPO / "src").returncode == 0
