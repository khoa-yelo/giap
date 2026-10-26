"""``giap setup``'s ToolRunner outside a source checkout (#347 review B2).

In an installed wheel ``Path(__file__).parents[3]`` is ``lib/python3.x``: the pinned images must be
read from the workflow bundled in the package (M11-01), and a missing ``GIAP_ROOT`` is a clear
``ToolUnavailable`` instead of a path under site-packages.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from giap.setupdb import run as sr


def _wheel(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> Path:
    """Pretend to be an installed package whose bundled workflow lives under ``tmp_path``."""
    wf = tmp_path / "site" / "giap" / "workflow"
    (wf / "conf").mkdir(parents=True)
    (wf / "main.nf").write_text("// bundled\n")
    monkeypatch.delenv("GIAP_REPO", raising=False)
    monkeypatch.setattr(sr, "_checkout", lambda: None)
    monkeypatch.setattr("giap.runner.nextflow.workflow_main", lambda: wf / "main.nf")
    return wf


def test_checkout_is_found_from_the_source_tree(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GIAP_REPO", raising=False)
    repo = sr._checkout()
    assert repo is not None and (repo / "workflow" / "main.nf").is_file()


def test_giap_repo_overrides_the_checkout(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("GIAP_REPO", str(tmp_path))
    assert sr._checkout() == tmp_path
    assert sr._containers_config() == tmp_path / "workflow" / "conf" / "containers.config"


def test_wheel_reads_the_bundled_containers_config(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    wf = _wheel(monkeypatch, tmp_path)
    assert sr._containers_config() == wf / "conf" / "containers.config"


def test_wheel_without_giap_root_is_tool_unavailable(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    """14f2: an installed wheel no longer needs $GIAP_ROOT; with no env, no image and no docker
    the tool is unavailable for that reason, not for the unset variable."""
    _wheel(monkeypatch, tmp_path)
    monkeypatch.delenv("GIAP_ROOT", raising=False)
    nodocker = tmp_path / "nodocker"
    nodocker.mkdir()
    monkeypatch.setenv("PATH", str(nodocker))
    with pytest.raises(sr.ToolUnavailable, match="build_envs") as exc:
        sr.ToolRunner().command("bakta", ["bakta_db", "download"])
    assert "GIAP_ROOT is not set" not in str(exc.value)


def test_wheel_with_giap_root_uses_its_env(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    _wheel(monkeypatch, tmp_path)
    root = tmp_path / "root"
    (root / "micromamba" / "envs" / "giap-bakta").mkdir(parents=True)
    monkeypatch.setenv("GIAP_ROOT", str(root))
    monkeypatch.delenv("MAMBA_ROOT_PREFIX", raising=False)
    cmd = sr.ToolRunner().command("bakta", ["bakta_db", "download"])
    assert cmd[:4] == ["micromamba", "run", "-n", "giap-bakta"]


def test_wheel_image_lookup_reads_the_bundled_pin(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    wf = _wheel(monkeypatch, tmp_path)
    ref = "quay.io/biocontainers/bakta@sha256:" + "a" * 64
    (wf / "conf" / "containers.config").write_text(
        f"process {{\n  withName: 'BAKTA' {{ container = 'docker://{ref}' }}\n}}\n"
    )
    root = tmp_path / "root"
    img = root / "containers" / ("quay.io-biocontainers-bakta@sha256-" + "a" * 64 + ".img")
    img.parent.mkdir(parents=True)
    img.write_text("")
    monkeypatch.setenv("GIAP_ROOT", str(root))
    monkeypatch.delenv("MAMBA_ROOT_PREFIX", raising=False)
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR", raising=False)
    cmd = sr.ToolRunner().command("bakta", ["bakta_db", "download"])
    assert cmd[0] == "apptainer" and cmd[-3] == str(img)
