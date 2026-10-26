"""Locked, checksummed Dockerfiles and the pin-equality test (M11-14b).

Spec docs/specs/M11-14-container-install.md item 2 ("Locked, checked builds"), the edge case "key
unchanged, upstream source changed", AC 1 (a ``GIT_SHA`` or env pin differing from
``annotate/tools.py`` fails) and AC 2 (every ``ADD`` of a URL has ``--checksum``); Q-L12
(BacTermFinder locks); ADR-0060 §1. Static: the Dockerfiles, env files and locks are only read.

Readings chosen where the spec leaves room:

* "Installs from the explicit lock": every ``micromamba install|create`` in the Dockerfile reads
  ``-f <the COPY destination of workflow/envs/locks/<tool>.explicit.txt>``; none reads the env
  ``.yml``. When a pip lock exists it is copied and installed by a ``RUN`` that contains
  ``pip install`` and ``--no-deps`` (nothing is resolved at build time; the ``RUN`` may filter the
  ``name @ file://`` lines that are conda's own packages first).
* ``--checksum=sha256:<64 hex>`` may be literal or come from an ``ARG`` default.
* uniop is out of scope for the ``ADD``/lock rules here (14u1 replaces its ``ADD`` with ``COPY``);
  its ``GIT_SHA`` still has to equal the registry when present.
* A lock is coherent with its env file: every ``name=version`` conda pin of the ``.yml`` is in the
  explicit lock, every ``name==version`` pip pin in the pip lock; the pip lock never contradicts
  the explicit lock (same package, other version); a ``.yml`` without a ``pip:`` section has a pip
  lock with no pip-only package. (The committed python-3.10 pip locks hold user-site packages,
  e.g. ``ViennaRNA==2.7.1`` next to conda ``viennarna 2.6.4`` for ostir: a ``pip freeze`` without
  ``PYTHONNOUSERSITE=1``. Installing them in an image would replace the pinned tool.)
"""

from __future__ import annotations

import re
import shutil
from pathlib import Path

import pytest
import yaml

from giap.annotate.tools import TOOLS

REPO = Path(__file__).resolve().parents[2]
GHCR = ["bactermfinder", "ostir", "promoter_calculator", "promotech", "terminatornet"]
GIT = sorted(t for t in [*GHCR, "uniop"] if TOOLS[t].version.startswith("git:"))
PYPI_PKG = {"ostir": "ostir", "promoter_calculator": "promotercalculator"}
HEX40 = re.compile(r"\b[0-9a-f]{40}\b")


def _norm(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def _instructions(text: str) -> list[tuple[str, str]]:
    """``(INSTRUCTION, rest)`` per instruction; comments dropped, continuations joined."""
    lines = [ln for ln in text.splitlines() if not ln.lstrip().startswith("#")]
    joined = re.sub(r"\\[ \t]*\n", " ", "\n".join(lines))
    out = []
    for line in joined.splitlines():
        words = line.strip().split(None, 1)
        if len(words) == 2:
            out.append((words[0].upper(), words[1].strip()))
    return out


def _args(instr: list[tuple[str, str]]) -> dict[str, str]:
    out = {}
    for op, rest in instr:
        if op == "ARG" and "=" in rest:
            k, v = rest.split("=", 1)
            out[k.strip()] = v.strip().strip('"')
    return out


def _expand(s: str, args: dict[str, str]) -> str:
    return re.sub(r"\$\{(\w+)\}|\$(\w+)", lambda m: args.get(m[1] or m[2], m[0]), s)


def _split_flags(rest: str) -> tuple[list[str], list[str]]:
    words = rest.split()
    flags = [w for w in words if w.startswith("--")]
    return flags, [w for w in words if not w.startswith("--")]


def _dockerfile(tool: str, root: Path | None = None) -> str:
    return ((root or REPO) / "workflow" / "containers" / tool / "Dockerfile").read_text()


def _copies(tool: str) -> dict[str, str]:
    """COPY source -> destination."""
    out = {}
    for op, rest in _instructions(_dockerfile(tool)):
        if op == "COPY":
            _, paths = _split_flags(rest)
            for src in paths[:-1]:
                out[src.removeprefix("./")] = paths[-1]
    return out


def _runs(tool: str) -> list[str]:
    return [rest for op, rest in _instructions(_dockerfile(tool)) if op == "RUN"]


def _url_adds(tool: str) -> list[tuple[list[str], list[str]]]:
    instr = _instructions(_dockerfile(tool))
    args = _args(instr)
    out = []
    for op, rest in instr:
        if op == "ADD":
            flags, paths = _split_flags(_expand(rest, args))
            if any(re.match(r"https?://", p) for p in paths[:-1]):
                out.append((flags, paths))
    return out


def _lock(tool: str, kind: str) -> Path:
    return REPO / "workflow" / "envs" / "locks" / f"{tool}.{kind}.txt"


def _env(tool: str) -> dict[str, object]:
    data: dict[str, object] = yaml.safe_load(
        (REPO / "workflow" / "envs" / f"{tool}.yml").read_text()
    )
    return data


def _yml_pins(tool: str) -> tuple[dict[str, str], dict[str, str]]:
    conda, pip = {}, {}
    for dep in _env(tool)["dependencies"]:  # type: ignore[union-attr]
        if isinstance(dep, str):
            name, _, ver = dep.partition("=")
            conda[_norm(name)] = ver.lstrip("=")
        else:
            for p in dep["pip"]:
                name, _, ver = p.partition("==")
                pip[_norm(name)] = ver
    return conda, pip


def _explicit_pkgs(tool: str) -> dict[str, str]:
    text = _lock(tool, "explicit").read_text()
    pat = r"/([^/\s]+)-([^-/\s]+)-([^-/\s]+)\.(?:conda|tar\.bz2)"
    return {_norm(m[1]): m[2] for m in re.finditer(pat, text)}


def _pip_pkgs(tool: str) -> dict[str, str]:
    out = {}
    for line in _lock(tool, "pip").read_text().splitlines():
        if "==" in line and not line.lstrip().startswith("#"):
            name, ver = line.strip().split("==", 1)
            out[_norm(name)] = ver.strip()
    return out


def _pin_mismatches(tool: str, root: Path) -> list[str]:
    """Pins in the Dockerfile and env file of ``tool`` that differ from ``annotate/tools.py``."""
    spec = TOOLS[tool]
    bad = []
    docker = _dockerfile(tool, root)
    env = (root / "workflow" / "envs" / f"{tool}.yml").read_text()
    if spec.version.startswith("git:"):
        sha = spec.version.removeprefix("git:")
        got = _args(_instructions(docker)).get("GIT_SHA")
        if got is not None and got != sha:
            bad.append(f"Dockerfile GIT_SHA {got}")
        for where, text in (("Dockerfile", docker), ("env", env)):
            bad += [f"{where} commit {h}" for h in HEX40.findall(text) if h != sha]
    if tool in PYPI_PKG:
        pins = re.findall(rf"(?im)^\s*-\s*{PYPI_PKG[tool]}==(\S+)", env)
        if pins != [spec.version]:
            bad.append(f"env pin {pins}")
    return bad


# --- install from the locks ---------------------------------------------------------------


@pytest.mark.parametrize("tool", GHCR)
def test_dockerfile_copies_and_installs_the_explicit_lock(tool: str) -> None:
    copies = _copies(tool)
    src = f"workflow/envs/locks/{tool}.explicit.txt"
    assert src in copies, f"{tool}: COPY {src} missing"
    mamba = [r for r in _runs(tool) if re.search(r"micromamba\s+(install|create)", r)]
    assert mamba, f"{tool}: no micromamba install/create"
    for r in mamba:
        assert re.search(rf"-f\s+{re.escape(copies[src])}(\s|$)", r), r


@pytest.mark.parametrize("tool", GHCR)
def test_dockerfile_never_solves_the_env_file(tool: str) -> None:
    copies = _copies(tool)
    yml = f"workflow/envs/{tool}.yml"
    for r in _runs(tool):
        if re.search(r"micromamba\s+(install|create)", r):
            assert "env.yml" not in r and ".yml" not in r, r
            if yml in copies:
                assert copies[yml] not in r, r


@pytest.mark.parametrize("tool", GHCR)
def test_dockerfile_installs_the_pip_lock_without_resolving(tool: str) -> None:
    if not _lock(tool, "pip").is_file():
        pytest.skip(f"{tool} has no pip lock")
    copies = _copies(tool)
    src = f"workflow/envs/locks/{tool}.pip.txt"
    assert src in copies, f"{tool}: COPY {src} missing"
    hits = [r for r in _runs(tool) if copies[src] in r and "pip install" in r]
    assert hits, f"{tool}: no RUN pip install of {copies[src]}"
    assert all("--no-deps" in r for r in hits), hits


# --- ADD --checksum ---------------------------------------------------------------------------


@pytest.mark.parametrize("tool", GHCR)
def test_every_url_add_has_a_sha256_checksum(tool: str) -> None:
    for flags, paths in _url_adds(tool):
        sums = [f for f in flags if f.startswith("--checksum=")]
        assert len(sums) == 1, f"{tool}: ADD {paths[0]} without --checksum"
        assert re.fullmatch(r"--checksum=sha256:[0-9a-f]{64}", sums[0]), sums[0]


@pytest.mark.parametrize("tool", [t for t in GIT if t != "uniop"])
def test_git_tool_adds_exactly_the_pinned_archive(tool: str) -> None:
    sha = TOOLS[tool].version.removeprefix("git:")
    adds = _url_adds(tool)
    assert len(adds) == 1, adds
    assert adds[0][1][0] == f"{TOOLS[tool].source_repo}/archive/{sha}.tar.gz"


@pytest.mark.parametrize("tool", sorted(set(GHCR) - set(GIT)))
def test_pypi_tools_fetch_no_url(tool: str) -> None:
    assert _url_adds(tool) == []


def test_archive_checksums_are_distinct() -> None:
    sums = [f for t in GHCR for flags, _ in _url_adds(t) for f in flags if "--checksum" in f]
    assert len(sums) == len(set(sums)) == len([t for t in GIT if t != "uniop"]), sums


# --- pin equality (AC 1) ----------------------------------------------------------------------


@pytest.mark.parametrize("tool", sorted({*GHCR, "uniop"}))
def test_dockerfile_and_env_pins_equal_the_registry(tool: str) -> None:
    assert _pin_mismatches(tool, REPO) == []


@pytest.mark.parametrize("tool", [t for t in GIT if t != "uniop"])
def test_git_dockerfile_declares_git_sha(tool: str) -> None:
    sha = TOOLS[tool].version.removeprefix("git:")
    assert _args(_instructions(_dockerfile(tool))).get("GIT_SHA") == sha


def _copy_tool(tool: str, tmp_path: Path) -> Path:
    for rel in (f"workflow/containers/{tool}/Dockerfile", f"workflow/envs/{tool}.yml"):
        (tmp_path / rel).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(REPO / rel, tmp_path / rel)
    return tmp_path


def test_a_changed_git_sha_is_caught(tmp_path: Path) -> None:
    root = _copy_tool("promotech", tmp_path)
    p = root / "workflow/containers/promotech/Dockerfile"
    sha = TOOLS["promotech"].version.removeprefix("git:")
    p.write_text(p.read_text().replace(f"GIT_SHA={sha}", "GIT_SHA=" + "a" * 40))
    assert _pin_mismatches("promotech", root)


def test_a_changed_env_pin_is_caught(tmp_path: Path) -> None:
    root = _copy_tool("ostir", tmp_path)
    p = root / "workflow/envs/ostir.yml"
    p.write_text(p.read_text().replace("ostir==1.1.3", "ostir==1.1.4"))
    assert _pin_mismatches("ostir", root) == ["env pin ['1.1.4']"]


# --- locks: BacTermFinder (Q-L12) and coherence -------------------------------------------------


def test_bactermfinder_has_both_locks() -> None:
    for kind in ("explicit", "pip"):
        assert _lock("bactermfinder", kind).is_file(), f"locks/bactermfinder.{kind}.txt missing"


@pytest.mark.parametrize("tool", GHCR)
def test_explicit_lock_is_a_hash_pinned_url_list(tool: str) -> None:
    lines = [ln for ln in _lock(tool, "explicit").read_text().splitlines() if ln.strip()]
    body = [ln for ln in lines if not ln.startswith("#")]
    assert body[0] == "@EXPLICIT", body[:1]
    assert len(body) > 1
    for ln in body[1:]:
        assert re.fullmatch(r"https://\S+\.(conda|tar\.bz2)#[0-9a-f]{32,64}", ln), ln


@pytest.mark.parametrize("tool", GHCR)
def test_env_conda_pins_are_in_the_explicit_lock(tool: str) -> None:
    conda, _ = _yml_pins(tool)
    got = _explicit_pkgs(tool)
    assert {n: got.get(n) for n in conda} == conda


@pytest.mark.parametrize("tool", GHCR)
def test_env_pip_pins_are_in_the_pip_lock(tool: str) -> None:
    _, pip = _yml_pins(tool)
    if not pip:
        return
    got = _pip_pkgs(tool)
    assert {n: got.get(n) for n in pip} == pip


@pytest.mark.parametrize("tool", GHCR)
def test_pip_lock_never_contradicts_the_explicit_lock(tool: str) -> None:
    if not _lock(tool, "pip").is_file():
        pytest.skip(f"{tool} has no pip lock")
    conda = _explicit_pkgs(tool)
    clash = {n: (v, conda[n]) for n, v in _pip_pkgs(tool).items() if n in conda and conda[n] != v}
    assert clash == {}, f"{tool}: pip lock vs explicit lock (user-site packages?): {clash}"


@pytest.mark.parametrize("tool", GHCR)
def test_no_pip_only_packages_without_a_pip_section(tool: str) -> None:
    _, pip = _yml_pins(tool)
    if pip or not _lock(tool, "pip").is_file():
        return
    conda = _explicit_pkgs(tool)
    extra = sorted(n for n in _pip_pkgs(tool) if n not in conda)
    assert extra == [], f"{tool}.yml has no pip section, but its pip lock adds {extra}"


@pytest.mark.parametrize("tool", GHCR)
def test_pip_lock_install_drops_conda_file_lines(tool: str) -> None:
    """Review item 9 (M13): `name @ file://…` lines are conda's own packages; pip cannot fetch them.

    The RUN that pip-installs the lock filters them out first (grep -v / sed on `@ file://`)."""
    if not _lock(tool, "pip").is_file():
        pytest.skip(f"{tool} has no pip lock")
    dest = _copies(tool)[f"workflow/envs/locks/{tool}.pip.txt"]
    hits = [r for r in _runs(tool) if dest in r and "pip install" in r]
    assert hits, f"{tool}: no RUN pip install of {dest}"
    for r in hits:
        assert re.search(r"(grep\s+-v|sed\b)[^;&|]*@ file://", r), r
        filt = r.index("@ file://")
        assert filt < r.index("pip install"), "filter the lock before pip reads it"
