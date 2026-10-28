"""``workflow_main()`` in a checkout and in an installed package (M11-01 edge cases).

Both layouts are faked in ``tmp_path`` by pointing the module's ``__file__`` at a fake
``gippa/runner/nextflow.py``. Also: nothing in the shipped workflow writes into ``projectDir``
(site-packages may be read-only).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

import gippa.runner.nextflow as nf

REPO = Path(__file__).resolve().parents[2]


def _touch(p: Path) -> Path:
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text("x\n")
    return p


def _module_at(monkeypatch: pytest.MonkeyPatch, pkg: Path) -> None:
    monkeypatch.setattr(nf, "__file__", str(_touch(pkg / "runner" / "nextflow.py")))


def test_checkout_layout_prefers_repo_workflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    repo = tmp_path.resolve() / "repo"
    main = _touch(repo / "workflow" / "main.nf")
    _module_at(monkeypatch, repo / "src" / "gippa")
    assert nf.workflow_main() == main


def test_installed_layout_uses_package_workflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    site = tmp_path.resolve() / "venv" / "lib" / "python3.12" / "site-packages"
    main = _touch(site / "gippa" / "workflow" / "main.nf")
    _module_at(monkeypatch, site / "gippa")
    got = nf.workflow_main()
    assert got == main and got.is_file()
    assert got.parent == site / "gippa" / "workflow"


def test_checkout_wins_when_both_exist(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """An editable install has the package copy too (if built in place): the checkout still wins."""
    repo = tmp_path.resolve() / "repo"
    main = _touch(repo / "workflow" / "main.nf")
    _touch(repo / "src" / "gippa" / "workflow" / "main.nf")
    _module_at(monkeypatch, repo / "src" / "gippa")
    assert nf.workflow_main() == main


def test_this_checkout_resolves_its_own_workflow() -> None:
    assert nf.workflow_main() == (REPO / "workflow" / "main.nf").resolve()


# --- read-only projectDir ---------------------------------------------------------------------

_SHIPPED = [
    p
    for p in (REPO / "workflow").rglob("*")
    if p.is_file()
    and "tests" not in p.relative_to(REPO / "workflow").parts
    and p.suffix in {".nf", ".config", ".groovy", ".sh", ".R", ".py"}
]
_WRITES = [
    re.compile(r"publishDir\b[^\n]*projectDir"),
    re.compile(r"\bworkDir\s*=[^\n]*projectDir"),
    re.compile(r"(>>?|\btee\b|\bmkdir\b|\btouch\b|\brm\b)\s+(-\S+\s+)*[\"']?\$\{?projectDir"),
    re.compile(r"\b(cp|mv|ln)\b[^\n]*\s[\"']?\$\{?projectDir[^\s\"']*[\"']?\s*$"),
    re.compile(
        r"projectDir[^\n]*\.(text\s*=|append\b|write\b|mkdirs?\b|copyTo\b|moveTo\b|delete\b)"
    ),
    re.compile(r"projectDir[^\n]*<<"),
]
_BIND = re.compile(r"-B\s+\$\{projectDir\}(\S+?)(\s|\"|$)")


def test_shipped_workflow_files_found() -> None:
    names = {p.name for p in _SHIPPED}
    assert {"main.nf", "nextflow.config"} <= names


@pytest.mark.parametrize("path", _SHIPPED, ids=lambda p: p.relative_to(REPO).as_posix())
def test_no_writes_into_projectdir(path: Path) -> None:
    bad = []
    for n, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if "projectDir" not in line:
            continue
        if any(r.search(line) for r in _WRITES):
            bad.append((n, line.strip()))
        for m in _BIND.finditer(line):
            if m.group(1) != "/assets/ps/ps:/usr/bin/ps:ro":
                bad.append((n, line.strip()))
    assert bad == []
