"""M11-15g: the step images' hashed Python lock agrees with ``uv.lock`` (item 8 AC 5).

Spec docs/specs/M11-15-portable-clusters.md item 8a: ``pip install --require-hashes -r
requirements.lock`` "(from ``uv export --frozen --no-dev --no-emit-project``, hashes kept,
committed, checked against ``uv.lock`` by a test)"; AC 5 "re-export with ``--no-emit-project``,
hashes kept".

Readings chosen where the spec leaves room:

* The lock is ``workflow/containers/step/requirements.lock``, uv's export format, with uv's header
  naming the command (``--frozen``, ``--no-dev``, ``--no-emit-project``; never ``--no-hashes``, no
  ``--extra``/``--all-extras``: the image carries the runtime dependencies only).
* Static check (always): every requirement is ``name==version`` (one line per uv.lock fork,
  e.g. numpy for Python < 3.12 and >= 3.12) with at least one ``--hash=sha256:``; the set of
  (name, version) is exactly the closure of the ``gippa`` package's runtime ``dependencies`` in
  ``uv.lock`` (extras followed only where a dependency asks for one; markers ignored, as uv export
  keeps every marker); each hash set equals the ``uv.lock`` hashes of that package (sdist +
  wheels). No editable, path, URL or ``gippa`` line.
* Re-export (when a ``uv`` binary is on PATH): ``uv export --frozen --no-dev --no-emit-project``
  offline in the repository gives the same (name, version, marker, hashes) entries.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import tomllib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
LOCK = REPO / "workflow" / "containers" / "step" / "requirements.lock"
UV_LOCK = REPO / "uv.lock"
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
EXPORT = ["uv", "export", "--frozen", "--no-dev", "--no-emit-project"]

Entry = tuple[str, str, str, frozenset[str]]  # name, version, marker, hashes


def _norm(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def _parse(text: str) -> dict[tuple[str, str], Entry]:
    """uv/pip requirements with hashes -> {(normalised name, version): entry}.

    A name may appear once per version: uv.lock forks by marker (e.g. ``numpy`` for
    ``python_full_version < '3.12'`` and ``>= '3.12'``).
    """
    joined = re.sub(r"\\[ \t]*\n", " ", text)
    out: dict[tuple[str, str], Entry] = {}
    for raw in joined.splitlines():
        line = raw.split(" #", 1)[0].strip()
        if not line or line.startswith("#"):
            continue
        hashes = frozenset(re.findall(r"--hash=sha256:([0-9a-f]{64})", line))
        req = re.sub(r"\s--hash=\S+", "", line).strip()
        m = re.fullmatch(r"([A-Za-z0-9][A-Za-z0-9._\-]*)==(\S+?)\s*(?:;\s*(.*))?", req)
        assert m, f"not a pinned requirement: {raw!r}"
        key = (_norm(m.group(1)), m.group(2))
        assert key not in out, f"duplicate requirement {key}"
        out[key] = (key[0], key[1], (m.group(3) or "").strip(), hashes)
    return out


def _lock_text() -> str:
    assert LOCK.is_file(), f"{LOCK.relative_to(REPO)} is missing (M11-15g)"
    return LOCK.read_text(encoding="utf-8")


def _uv_packages() -> dict[tuple[str, str], dict]:
    lock = tomllib.loads(UV_LOCK.read_text(encoding="utf-8"))
    return {(_norm(p["name"]), p.get("version", "")): p for p in lock["package"]}


def _resolve(dep: dict, pkgs: dict[tuple[str, str], dict]) -> tuple[str, str]:
    """A ``uv.lock`` dependency -> its package key (``version`` is given when a name forks)."""
    name = _norm(dep["name"])
    if "version" in dep:
        return name, dep["version"]
    (key,) = [k for k in pkgs if k[0] == name]
    return key


def _runtime_closure() -> set[tuple[str, str]]:
    pkgs = _uv_packages()
    (root,) = [p for k, p in pkgs.items() if k[0] == _norm(PKG)]
    seen: set[tuple[tuple[str, str], tuple[str, ...]]] = set()
    stack = [(_resolve(d, pkgs), tuple(d.get("extra", ()))) for d in root.get("dependencies", [])]
    while stack:
        key, extras = stack.pop()
        if (key, extras) in seen:
            continue
        seen.add((key, extras))
        pkg = pkgs[key]
        deps = list(pkg.get("dependencies", []))
        for e in extras:
            deps += pkg.get("optional-dependencies", {}).get(e, [])
        stack += [(_resolve(d, pkgs), tuple(d.get("extra", ()))) for d in deps]
    return {k for k, _ in seen}


def _uv_hashes(pkg: dict) -> frozenset[str]:
    files = ([pkg["sdist"]] if "sdist" in pkg else []) + list(pkg.get("wheels", []))
    return frozenset(f["hash"].removeprefix("sha256:") for f in files if "hash" in f)


# ---------------------------------------------------------------- the committed file


def test_lock_exists() -> None:
    assert LOCK.is_file(), f"{LOCK.relative_to(REPO)} is missing (M11-15g)"


def test_header_names_the_export_command() -> None:
    header = [ln for ln in _lock_text().splitlines()[:5] if ln.startswith("#")]
    cmd = next((ln.lstrip("# ").strip() for ln in header if "uv export" in ln), None)
    assert cmd, f"no uv export command in the header: {header}"
    words = cmd.split()
    for flag in ("--frozen", "--no-dev", "--no-emit-project"):
        assert flag in words, (flag, cmd)
    for flag in ("--no-hashes", "--all-extras", "--extra", "--all-groups", "--group", "-e"):
        assert flag not in words, (flag, cmd)


def test_every_requirement_is_pinned_and_hashed() -> None:
    entries = _parse(_lock_text())
    assert entries, "empty lock"
    for name, _, _, hashes in entries.values():
        assert hashes, f"{name} has no --hash (pip --require-hashes would refuse it)"


def test_no_project_editable_path_or_url_line() -> None:
    text = _lock_text()
    assert not re.search(r"^\s*-e\s", text, re.M), "editable line"
    assert "file://" not in text and not re.search(r"^\s*\.{0,2}/", text, re.M)
    assert not re.search(r"^\S+\s+@\s+\S+", text, re.M), "direct URL requirement"
    names = {n for n, _ in _parse(text)}
    assert _norm(PKG) not in names, f"{PKG} itself is in the lock (--no-emit-project)"


def test_names_are_the_runtime_closure_of_the_package() -> None:
    got = set(_parse(_lock_text()))
    want = _runtime_closure()
    assert got == want, {"extra": sorted(got - want), "missing": sorted(want - got)}


def test_no_dev_or_docs_only_package() -> None:
    got = {n for n, _ in _parse(_lock_text())}
    for dev in ("pytest", "hypothesis", "ruff", "mypy", "mkdocs-material", "towncrier"):
        assert dev not in got, dev


def test_versions_and_hashes_equal_uv_lock() -> None:
    pkgs = _uv_packages()
    for name, version, _, hashes in _parse(_lock_text()).values():
        assert (name, version) in pkgs, f"{name}=={version} not in uv.lock"
        assert hashes == _uv_hashes(pkgs[name, version]), f"{name}: hashes differ from uv.lock"


# ---------------------------------------------------------------- re-export (AC 5)


def test_reexport_from_uv_lock_is_identical() -> None:
    uv = shutil.which("uv")
    if uv is None:
        pytest.skip("no uv binary on PATH; the static uv.lock comparison above still holds")
    env = {**os.environ, "UV_OFFLINE": "1", "UV_NO_SYNC": "1"}
    done = subprocess.run(
        [uv, *EXPORT[1:]], cwd=REPO, env=env, capture_output=True, text=True, timeout=120
    )
    assert done.returncode == 0, done.stderr
    assert _parse(_lock_text()) == _parse(done.stdout), "requirements.lock is stale: re-export it"
