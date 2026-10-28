"""M11-15g review round 1, B1 and a nit: hashed build backend for the sdist-only dependencies.

Review of PR #436, B1: ``pandarallel`` and ``upsetplot`` have only an sdist in
``workflow/containers/step/requirements.lock``. Under pip's default build isolation their build
backend (setuptools) is fetched by a separate pip run that does not inherit ``--require-hashes``,
so an unhashed setuptools enters the build (ADR-0060 decision 1, ADR-0065). Required: the
hashed ``setuptools`` pinned in the lock is installed first, and the lock is installed with
``--no-build-isolation``. Static: the Dockerfiles are only read; nothing is built.

Readings chosen where the review leaves room:

* pip accepts ``--hash`` only inside a requirements file, so the setuptools bootstrap is a
  ``pip install ... -r <file>`` whose ``<file>`` is a committed build-context file (bind-mounted
  or COPYed, like the lock). Its only requirement is ``setuptools==<the lock's version>`` carrying
  exactly the lock's hashes for setuptools; the command (or the file) enables hash checking.
* The bootstrap runs before the lock install (an earlier segment, same or earlier ``RUN``), and
  the lock install carries ``--no-build-isolation``.
* Nit: ``STEP_IMAGES`` in ``scripts/release/image_key.py`` and ``scripts/release/readiness.py``
  name the same images (the two copies are kept in sync by this test).
"""

from __future__ import annotations

import importlib.util
import re
import sys
from pathlib import Path
from types import ModuleType

import pytest

REPO = Path(__file__).resolve().parents[2]
CONTAINERS = REPO / "workflow" / "containers"
IMAGES = {
    "step": CONTAINERS / "step" / "Dockerfile",
    "step-r": CONTAINERS / "step-r" / "Dockerfile",
}
LOCK = "workflow/containers/step/requirements.lock"
PIP_INSTALL = re.compile(r"\bpip3?\s+install\b")
REQ_FLAG = re.compile(r"(?:^|\s)(?:-r|--requirement)(?:\s+|=)(\S+)")


# ---------------------------------------------------------------- helpers


def _text(path: Path) -> str:
    assert path.is_file(), f"{path.relative_to(REPO)} is missing (M11-15g)"
    return path.read_text(encoding="utf-8")


def _instructions(text: str) -> list[tuple[str, str]]:
    lines = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    joined = re.sub(r"\\[ \t]*\n", " ", "\n".join(lines))
    out = []
    for line in joined.splitlines():
        words = line.strip().split(None, 1)
        if len(words) == 2:
            out.append((words[0].upper(), words[1].strip()))
    return out


def _context_files(path: Path) -> dict[str, str]:
    """In-image path -> build-context source, from ``COPY`` (no ``--from``) and bind mounts."""
    out: dict[str, str] = {}
    for op, rest in _instructions(_text(path)):
        if op == "COPY":
            words = rest.split()
            flags = [w for w in words if w.startswith("--")]
            paths = [w for w in words if not w.startswith("--")]
            if any(f.startswith("--from") for f in flags) or len(paths) < 2:
                continue
            for src in paths[:-1]:
                src = src.removeprefix("./")
                dest = paths[-1]
                out[dest + Path(src).name if dest.endswith("/") else dest] = src
        elif op == "RUN":
            for m in re.finditer(r"--mount=(\S+)", rest):
                kv = dict(p.split("=", 1) for p in m.group(1).split(",") if "=" in p)
                if kv.get("type", "bind") == "bind" and "source" in kv and "target" in kv:
                    out[kv["target"]] = kv["source"].removeprefix("./")
    return out


def _pip(path: Path) -> list[str]:
    out = []
    for op, rest in _instructions(_text(path)):
        if op == "RUN":
            body = re.sub(r"--mount=\S+", " ", rest)
            out += [s.strip() for s in re.split(r"&&|\|\||;|\|", body) if s.strip()]
    return [s for s in out if PIP_INSTALL.search(s)]


def _requirements(text: str) -> dict[str, tuple[str, set[str]]]:
    """``name -> (version, hashes)`` of a hashed requirements file (comments dropped)."""
    joined = re.sub(r"\\[ \t]*\n", " ", text)
    out: dict[str, tuple[str, set[str]]] = {}
    for line in joined.splitlines():
        line = line.split("#", 1)[0].strip()
        m = re.match(r"([A-Za-z0-9][A-Za-z0-9._\-]*)==(\S+)", line)
        if m:
            name = re.sub(r"[-_.]+", "-", m.group(1)).lower()
            out[name] = (m.group(2), set(re.findall(r"--hash=sha256:([0-9a-f]{64})", line)))
    return out


def _locked_setuptools() -> tuple[str, set[str]]:
    reqs = _requirements((REPO / LOCK).read_text(encoding="utf-8"))
    assert "setuptools" in reqs, "setuptools is not in requirements.lock"
    version, hashes = reqs["setuptools"]
    assert hashes, "setuptools has no hash in requirements.lock"
    return version, hashes


def _lock_cmd(path: Path) -> str:
    ctx = _context_files(path)
    cmds = [s for s in _pip(path) if any(ctx.get(t) == LOCK for t in REQ_FLAG.findall(s))]
    assert len(cmds) == 1, f"{path.relative_to(REPO)}: expected one pip install of {LOCK}: {cmds}"
    return cmds[0]


def _bootstraps(path: Path) -> list[tuple[str, str, str]]:
    """``(pip command, context source, its text)`` for each pip install of a committed file
    whose only requirement is setuptools."""
    ctx = _context_files(path)
    out = []
    for cmd in _pip(path):
        for target in REQ_FLAG.findall(cmd):
            src = ctx.get(target)
            if src is None or src == LOCK or not (REPO / src).is_file():
                continue
            text = (REPO / src).read_text(encoding="utf-8")
            if set(_requirements(text)) == {"setuptools"}:
                out.append((cmd, src, text))
    return out


# ---------------------------------------------------------------- B1


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_lock_is_installed_without_build_isolation(image: str) -> None:
    cmd = _lock_cmd(IMAGES[image])
    assert "--no-build-isolation" in cmd.split(), (
        f"{image}: the sdist-only locked packages (pandarallel, upsetplot) are built in an "
        f"isolated env that fetches an unhashed setuptools; add --no-build-isolation: {cmd}"
    )
    assert "--require-hashes" in cmd.split(), cmd


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_hashed_locked_setuptools_is_installed_before_the_lock(image: str) -> None:
    path = IMAGES[image]
    boots = _bootstraps(path)
    assert len(boots) == 1, (
        f"{image}: expected one pip install -r <committed file pinning only setuptools> before "
        f"the lock, got {[b[0] for b in boots]}"
    )
    cmd, src, text = boots[0]
    version, hashes = _locked_setuptools()
    got_version, got_hashes = _requirements(text)["setuptools"]
    assert got_version == version, f"{src}: setuptools=={got_version}, the lock has {version}"
    assert got_hashes == hashes, f"{src}: setuptools hashes differ from requirements.lock"
    assert "--require-hashes" in cmd.split() or "--hash=" in text, cmd
    pip = _pip(path)
    assert pip.index(cmd) < pip.index(_lock_cmd(path)), pip


def test_both_images_bootstrap_from_the_same_file() -> None:
    srcs = {image: [b[1] for b in _bootstraps(path)] for image, path in IMAGES.items()}
    assert srcs["step"] and srcs["step"] == srcs["step-r"], srcs


# ---------------------------------------------------------------- nit: one STEP_IMAGES


def _load(name: str) -> ModuleType:
    path = REPO / "scripts" / "release" / f"{name}.py"
    spec = importlib.util.spec_from_file_location(f"_m11_15g_review_{name}", path)
    assert spec is not None and spec.loader is not None, path
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # dataclasses look their module up while it executes
    try:
        spec.loader.exec_module(mod)
    finally:
        sys.modules.pop(spec.name, None)
    return mod


def test_step_images_constant_is_the_same_in_image_key_and_readiness() -> None:
    ik, rd = _load("image_key"), _load("readiness")
    assert set(ik.STEP_IMAGES) == set(rd.STEP_IMAGES), (ik.STEP_IMAGES, rd.STEP_IMAGES)
    assert set(ik.STEP_IMAGES) == set(IMAGES), ik.STEP_IMAGES
