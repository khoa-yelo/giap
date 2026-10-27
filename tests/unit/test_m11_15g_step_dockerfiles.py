"""M11-15g: the step images ``<pkg>-step`` and ``<pkg>-step-r`` (static Dockerfile tests).

Spec docs/specs/M11-15-portable-clusters.md item 8a (with the B9/B12 fixes: procps, the hashed
lock, the marker file) and item 8 AC 4: "base by digest; ``procps`` installed; ``pip install
--require-hashes -r requirements.lock``, then the wheel with ``--no-deps``; the marker file; labels;
no conda, R or tool installs (step-r: exactly the fgsea lock added)". Static: the Dockerfiles are
only read; nothing is built.

Readings chosen where the spec leaves room:

* Paths: ``workflow/containers/step/Dockerfile`` (named by the spec) and
  ``workflow/containers/step-r/Dockerfile`` (the directory is the image suffix, as for the tool
  images). Both are built from the repository root, like every other image. The shared lock is
  ``workflow/containers/step/requirements.lock``.
* A build-context file reaches a ``RUN`` either by ``COPY src dest`` or by
  ``RUN --mount=type=bind,source=src,target=dest``; both are accepted.
* "As a separate invocation": the lock and the wheel are two different ``pip install`` commands
  (they may share one ``RUN`` joined by ``&&``); the lock command comes first, carries
  ``--require-hashes`` and ``-r``; the wheel command carries ``--no-deps`` and names a ``.whl`` (or
  a ``$…WHEEL…`` variable); no third ``pip install`` exists.
* The base: some ``FROM python:3.12-slim[-<debian>]@sha256:<64 hex>``; every ``FROM`` is pinned by
  digest or names an earlier stage. Both images use the same python digest (same interpreter as
  the image the lock was checked for).
* apt: every ``apt-get install`` has ``--no-install-recommends``; the union of apt packages is
  ``{procps, ca-certificates}`` with ``procps`` present (``ca-certificates`` may already be in the
  base).
* The marker ``/opt/gippa/STEP_IMAGE``: written by a shell command or COPYed from a committed
  file. In step its content is empty (``touch``, ``: >``, ``true >``, ``printf '' >``,
  ``echo -n >`` or a bare ``> file``); in step-r it is ``/opt/conda`` (trailing newline allowed).
* Labels: ``org.opencontainers.image.{version,source,revision,licenses}`` and
  ``org.gippa.wheel.sha256`` present in both. step: ``licenses="MIT"`` exactly. step-r aggregates R
  (Q-K37): its ``licenses`` names ``MIT`` and ``GPL-2`` (as the fgsea image does).
* step-r adds exactly the fgsea explicit lock (``micromamba install|create … -f <its dest>``), no
  other conda lock, no ``.yml`` solve, no R package install, no pip lock but requirements.lock.
* The step images are not ADR-0060 keyed tool images: ``scripts/release/image_key.ghcr_tools``
  and the readiness H4 images check do not list them (their tag is ``<version>``, the Q-K39 ADR).
"""

from __future__ import annotations

import re
import shlex
import sys
import tomllib
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
CONTAINERS = REPO / "workflow" / "containers"
STEP = CONTAINERS / "step" / "Dockerfile"
STEP_R = CONTAINERS / "step-r" / "Dockerfile"
IMAGES = {"step": STEP, "step-r": STEP_R}
LOCK = "workflow/containers/step/requirements.lock"
FGSEA_LOCK = "workflow/envs/locks/fgsea.explicit.txt"
MARKER = "/opt/gippa/STEP_IMAGE"
DIGEST = re.compile(r"@sha256:[0-9a-f]{64}$")
PYBASE = re.compile(
    r"^(?:docker\.io/)?(?:library/)?python:3\.12-slim(?:-[a-z]+)?@sha256:([0-9a-f]{64})$"
)
PIP_INSTALL = re.compile(r"\bpip3?\s+install\b")
CONDA = re.compile(r"\b(micromamba|mamba|conda)\b|/opt/conda")
R_WORDS = re.compile(r"\bRscript\b|\br-base\b|\bR\s+-e\b|install\.packages|BiocManager")


# ---------------------------------------------------------------- helpers


def _text(path: Path) -> str:
    assert path.is_file(), f"{path.relative_to(REPO)} is missing (M11-15g)"
    return path.read_text(encoding="utf-8")


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


def _of(path: Path, op: str) -> list[str]:
    return [rest for o, rest in _instructions(_text(path)) if o == op]


def _split_flags(rest: str) -> tuple[list[str], list[str]]:
    words = rest.split()
    return [w for w in words if w.startswith("--")], [w for w in words if not w.startswith("--")]


def _context_files(path: Path) -> dict[str, str]:
    """Build-context source -> in-image path, from ``COPY`` (no ``--from``) and bind mounts."""
    out: dict[str, str] = {}
    for op, rest in _instructions(_text(path)):
        if op == "COPY":
            flags, paths = _split_flags(rest)
            if any(f.startswith("--from") for f in flags) or len(paths) < 2:
                continue
            for src in paths[:-1]:
                out[src.removeprefix("./")] = paths[-1]
        elif op == "RUN":
            for m in re.finditer(r"--mount=(\S+)", rest):
                kv = dict(p.split("=", 1) for p in m.group(1).split(",") if "=" in p)
                if kv.get("type", "bind") == "bind" and "source" in kv and "target" in kv:
                    out[kv["source"].removeprefix("./")] = kv["target"]
    return out


def _segments(path: Path) -> list[str]:
    """Every shell command of every ``RUN``, in order (split on ``&&``, ``||``, ``;``, ``|``)."""
    out = []
    for rest in _of(path, "RUN"):
        body = re.sub(r"--mount=\S+", " ", rest)
        out += [s.strip() for s in re.split(r"&&|\|\||;|\|", body) if s.strip()]
    return out


def _pip(path: Path) -> list[str]:
    return [s for s in _segments(path) if PIP_INSTALL.search(s)]


def _lock_cmds(path: Path) -> list[str]:
    return [s for s in _pip(path) if re.search(r"(^|\s)(-r|--requirement)(\s|=)", s)]


def _wheel_cmds(path: Path) -> list[str]:
    wheel = re.compile(r"\.whl\b|\$\{?\w*WHEEL\w*\}?", re.IGNORECASE)
    return [s for s in _pip(path) if s not in _lock_cmds(path) and wheel.search(s)]


def _labels(path: Path) -> dict[str, str]:
    out = {}
    for rest in _of(path, "LABEL"):
        for k, v in re.findall(r'([\w.\-]+)=("[^"]*"|\S+)', rest):
            out[k] = v.strip('"')
    return out


def _froms(path: Path) -> list[tuple[str, str | None]]:
    """``(image, stage alias)`` per ``FROM``."""
    out = []
    for rest in _of(path, "FROM"):
        words = [w for w in rest.split() if not w.startswith("--")]
        alias = words[2] if len(words) >= 3 and words[1].upper() == "AS" else None
        out.append((words[0], alias))
    return out


def _python_digest(path: Path) -> str:
    found = [m.group(1) for img, _ in _froms(path) if (m := PYBASE.match(img))]
    assert found, f"{path.relative_to(REPO)}: no FROM python:3.12-slim@sha256:<digest>"
    return found[0]


def _marker_content(path: Path) -> str | None:
    """The marker's content as the Dockerfile writes it; ``None`` if it is never written."""
    ctx = _context_files(path)
    for src, dest in ctx.items():
        if dest == MARKER or (dest.endswith("/") and dest.rstrip("/") + "/STEP_IMAGE" == MARKER):
            return (REPO / src).read_text(encoding="utf-8").rstrip("\n")
    for seg in _segments(path):
        m = re.search(rf"^(.*?)>\s*{re.escape(MARKER)}\s*$", seg)
        if m is None:
            if re.fullmatch(rf"touch\s+{re.escape(MARKER)}", seg):
                return ""
            continue
        head = m.group(1).strip()
        if head in {"", ":", "true", "printf ''", 'printf ""', "echo -n", "echo -n ''"}:
            return ""
        words = shlex.split(head)
        if words and words[0] in {"printf", "echo"}:
            args = [w for w in words[1:] if w not in {"-n", "%s", "%s\\n"}]
            return "".join(args).replace("\\n", "").strip()
        return head
    return None


# ---------------------------------------------------------------- both images


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_dockerfile_exists(image: str) -> None:
    assert IMAGES[image].is_file(), f"{IMAGES[image].relative_to(REPO)} is missing (M11-15g)"


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_every_base_is_pinned_by_digest_or_an_earlier_stage(image: str) -> None:
    froms = _froms(IMAGES[image])
    assert froms, "no FROM"
    stages: set[str] = set()
    for img, alias in froms:
        assert img in stages or DIGEST.search(img), f"{image}: base not pinned by digest: {img}"
        if alias:
            stages.add(alias)


def test_both_images_use_the_same_python_slim_digest() -> None:
    assert _python_digest(STEP) == _python_digest(STEP_R)


def test_step_image_is_a_single_python_slim_stage() -> None:
    """step: "FROM python:3.12-slim@sha256:<pinned>" and nothing else to copy from."""
    froms = _froms(STEP)
    assert PYBASE.match(froms[-1][0]), froms
    assert all(PYBASE.match(img) for img, _ in froms), froms


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_apt_installs_only_procps_and_ca_certificates(image: str) -> None:
    path = IMAGES[image]
    apt = [s for s in _segments(path) if re.search(r"\bapt-get\s+(-\S+\s+)*install\b", s)]
    assert apt, f"{image}: no apt-get install (procps is needed for Nextflow's task metrics)"
    pkgs: set[str] = set()
    for seg in apt:
        assert "--no-install-recommends" in seg, seg
        after = re.split(r"\binstall\b", seg, maxsplit=1)[1].split()
        pkgs |= {w for w in after if not w.startswith("-")}
    assert "procps" in pkgs, pkgs
    assert pkgs <= {"procps", "ca-certificates"}, pkgs
    for seg in _segments(path):
        assert not re.search(r"\b(apt|apk|dnf|yum)\s+(add|install)\b", seg), seg


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_installs_the_hashed_lock_with_require_hashes(image: str) -> None:
    path = IMAGES[image]
    ctx = _context_files(path)
    assert LOCK in ctx, f"{image}: {LOCK} is neither COPYed nor bind-mounted"
    cmds = _lock_cmds(path)
    assert len(cmds) == 1, f"{image}: expected one pip install -r, got {cmds}"
    cmd = cmds[0]
    assert "--require-hashes" in cmd, cmd
    dest = ctx[LOCK]
    target = dest if not dest.endswith("/") else dest + "requirements.lock"
    assert re.search(rf"(-r|--requirement)(\s+|=){re.escape(target)}(\s|$)", cmd), (cmd, target)
    assert ".whl" not in cmd, f"the wheel must be a separate pip invocation: {cmd}"


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_installs_the_wheel_separately_with_no_deps_after_the_lock(image: str) -> None:
    path = IMAGES[image]
    wheels = _wheel_cmds(path)
    assert len(wheels) == 1, f"{image}: expected one pip install of the wheel, got {wheels}"
    assert "--no-deps" in wheels[0].split(), wheels[0]
    pip = _pip(path)
    assert pip.index(_lock_cmds(path)[0]) < pip.index(wheels[0]), pip
    assert len(pip) == 2, f"{image}: a pip install other than the lock and the wheel: {pip}"
    srcs = [s for s in _context_files(path) if s.startswith("dist/")]
    assert srcs, f"{image}: the wheel does not come from dist/ in the build context"
    assert all(s.endswith(".whl") or "WHEEL" in s.upper() for s in srcs), srcs


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_home_is_tmp(image: str) -> None:
    envs = " ".join(_of(IMAGES[image], "ENV"))
    assert re.search(r'(^|\s)HOME(=|\s+)"?/tmp"?(\s|$)', envs), envs


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_labels(image: str) -> None:
    labels = _labels(IMAGES[image])
    for key in (
        "org.opencontainers.image.version",
        "org.opencontainers.image.source",
        "org.opencontainers.image.revision",
        "org.opencontainers.image.licenses",
        "org.gippa.wheel.sha256",
    ):
        assert key in labels, f"{image}: label {key} missing ({sorted(labels)})"


def test_step_licence_label_is_mit() -> None:
    assert _labels(STEP)["org.opencontainers.image.licenses"] == "MIT"


def test_step_r_licence_label_names_mit_and_r_gpl2() -> None:
    lic = _labels(STEP_R)["org.opencontainers.image.licenses"]
    assert "MIT" in lic and "GPL-2" in lic, lic


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_fetches_no_url(image: str) -> None:
    for op, rest in _instructions(_text(IMAGES[image])):
        if op == "ADD":
            _, paths = _split_flags(rest)
            assert not [p for p in paths[:-1] if re.match(r"https?://", p)], rest
    for seg in _segments(IMAGES[image]):
        assert not re.search(r"\b(curl|wget|git\s+clone)\b", seg), seg


@pytest.mark.parametrize("image", sorted(IMAGES))
def test_installs_no_bioinformatics_tool(image: str) -> None:
    from gippa.annotate.tools import TOOLS

    allowed = {"codon_bias"} | ({"fgsea"} if image == "step-r" else set())
    text = " ".join(_of(IMAGES[image], "RUN") + _of(IMAGES[image], "COPY"))
    for tool in sorted(set(TOOLS) - allowed):
        assert not re.search(rf"\b{re.escape(tool)}\b", text, re.IGNORECASE), (image, tool)
    assert not [s for s in _context_files(IMAGES[image]) if s.startswith("workflow/envs/")] or (
        image == "step-r"
    )


# ---------------------------------------------------------------- step: no conda, no R


def test_step_has_no_conda_and_no_r() -> None:
    text = " ".join(
        rest
        for op, rest in _instructions(_text(STEP))
        if op in {"FROM", "RUN", "COPY", "ADD", "ENV"}
    )
    assert not CONDA.search(text), CONDA.search(text)
    assert not R_WORDS.search(text), R_WORDS.search(text)


def test_step_marker_is_empty() -> None:
    content = _marker_content(STEP)
    assert content is not None, f"step does not write {MARKER}"
    assert content == "", f"step marker must be empty, got {content!r}"


# ---------------------------------------------------------------- step-r: exactly the fgsea lock


def test_step_r_marker_is_the_conda_prefix() -> None:
    content = _marker_content(STEP_R)
    assert content is not None, f"step-r does not write {MARKER}"
    assert content == "/opt/conda", content


def test_step_r_adds_exactly_the_fgsea_explicit_lock() -> None:
    ctx = _context_files(STEP_R)
    conda_locks = sorted(s for s in ctx if s.startswith("workflow/envs/"))
    assert conda_locks == [FGSEA_LOCK], conda_locks
    mamba = [s for s in _segments(STEP_R) if re.search(r"micromamba\s+(install|create)", s)]
    assert len(mamba) == 1, mamba
    assert re.search(rf"-f\s+{re.escape(ctx[FGSEA_LOCK])}(\s|$)", mamba[0]), mamba[0]
    assert ".yml" not in mamba[0], mamba[0]
    for seg in _segments(STEP_R):
        assert not re.search(r"install\.packages|BiocManager|\bR\s+-e\b", seg), seg
        assert not re.search(r"\b(conda|mamba)\s+(install|create)\b", seg), seg


def test_step_r_conda_env_lives_at_the_marker_prefix() -> None:
    """The fgsea env must be at /opt/conda (micromamba's root prefix, or copied there)."""
    mamba = next(s for s in _segments(STEP_R) if re.search(r"micromamba\s+(install|create)", s))
    prefix = re.search(r"(?:-p|--prefix)\s+(\S+)", mamba)
    named = re.search(r"(?:-n|--name)\s+(\S+)", mamba)
    assert (prefix and prefix.group(1).rstrip("/") == "/opt/conda") or (
        named and named.group(1) == "base"
    ), mamba
    multi = len(_froms(STEP_R)) > 1
    if multi:
        copies = [r for r in _of(STEP_R, "COPY") if r.startswith("--from")]
        assert any(re.search(r"/opt/conda/?\s+/opt/conda/?$", c) for c in copies), copies


# ---------------------------------------------------------------- not ADR-0060 keyed tool images


def test_step_images_are_not_keyed_ghcr_tools() -> None:
    assert STEP.is_file() and STEP_R.is_file(), "step Dockerfiles missing (M11-15g)"
    sys.path.insert(0, str(REPO / "scripts" / "release"))
    try:
        import image_key
    finally:
        sys.path.pop(0)
    tools = image_key.ghcr_tools(REPO)
    assert "step" not in tools and "step-r" not in tools, tools


def test_readiness_h4_does_not_ask_images_json_for_step_images() -> None:
    assert STEP.is_file() and STEP_R.is_file(), "step Dockerfiles missing (M11-15g)"
    sys.path.insert(0, str(REPO / "scripts" / "release"))
    try:
        import readiness
    finally:
        sys.path.pop(0)
    reasons = readiness._h4_images(REPO)
    assert not [r for r in reasons if re.search(r"\bstep(-r)?\b", r)], reasons
