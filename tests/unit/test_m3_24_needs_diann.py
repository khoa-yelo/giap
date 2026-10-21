"""The ``needs_diann`` capability probe (M3-24 lead decision 12; GIAP_SPEC.md §3.6).

Khoa supplied DIA-NN as an Apptainer image, so the probe of ``tests/conftest.py`` must accept
``${DIANN_SIF:-<decision 1 path>}`` as well as the unpacked ``$GIAP_TOOLS_DIR/diann-2.7.0*``
build the spec expects. Nothing here runs DIA-NN: only ``Path.is_file`` is consulted.

The repository ``conftest.py`` is loaded as a fresh module under a private name so its
environment-dependent constants (``GIAP_TOOLS_DIR``) can be set per test; the plugin instance
pytest itself uses is never touched.
"""

from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest
from helpers.pilot_kit import DIANN_SIF, REPO

ROOT_CONFTEST = REPO / "tests" / "conftest.py"


def _conftest(monkeypatch: pytest.MonkeyPatch, tools: Path, sif: str | None) -> ModuleType:
    monkeypatch.setenv("GIAP_TOOLS_DIR", str(tools))
    if sif is None:
        monkeypatch.delenv("DIANN_SIF", raising=False)
    else:
        monkeypatch.setenv("DIANN_SIF", sif)
    spec = importlib.util.spec_from_file_location("giap_conftest_probe", ROOT_CONFTEST)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _probe(module: ModuleType) -> bool:
    return bool(module.CHECKS["needs_diann"]())


def test_the_image_satisfies_the_capability(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    image = tmp_path / "diann_2.3.1.sif"
    image.write_text("stub sif\n")
    module = _conftest(monkeypatch, tmp_path / "tools", str(image))
    assert _probe(module) is True


def test_the_capability_is_missing_without_image_or_binary(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _conftest(monkeypatch, tmp_path / "tools", str(tmp_path / "absent.sif"))
    assert _probe(module) is False


def test_an_unpacked_2_7_0_build_still_satisfies_it(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Decision 12: the ``$GIAP_TOOLS_DIR`` probe stays for a future unpacked build."""
    tools = tmp_path / "tools"
    (tools / "diann-2.7.0").mkdir(parents=True)
    (tools / "diann-2.7.0" / "diann").write_text("stub binary\n")
    module = _conftest(monkeypatch, tools, str(tmp_path / "absent.sif"))
    assert _probe(module) is True


def test_the_default_is_the_image_of_decision_1(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    module = _conftest(monkeypatch, tmp_path / "tools", None)
    assert module.DIANN_SIF_DEFAULT == DIANN_SIF
    assert _probe(module) is Path(DIANN_SIF).is_file()


def test_the_marker_is_declared_and_describes_the_image() -> None:
    markers = [
        line
        for line in (REPO / "pyproject.toml").read_text().splitlines()
        if line.strip().startswith('"needs_diann:')
    ]
    assert len(markers) == 1, markers
    assert "DIANN_SIF" in markers[0], markers[0]
