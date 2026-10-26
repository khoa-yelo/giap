"""Content-addressed image keys and tags, and the image cache (M11-14 items 2-3; ADR-0060 §1).

``root`` is the directory holding ``workflow/``. It defaults, at call time, to the checkout or to
the installed package that bundles ``workflow/`` (``giap.runner.nextflow.workflow_main()``), so a
wheel computes the same keys as the repository.
"""

from __future__ import annotations

import hashlib
import json
import re
import shlex
from collections.abc import Mapping
from pathlib import Path

__all__ = ["image_cache", "key", "key_inputs", "key_text", "tag", "uncovered_copy_sources"]

_GIT_VERSION = re.compile(r"git:([0-9a-f]{40})")
_TAG = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_.-]{0,127}")
_EXTRA_INPUTS: dict[str, tuple[str, ...]] = {"fgsea": ("workflow/bin/giap_fgsea.R",)}


def _root(root: Path | None) -> Path:
    if root is not None:
        return Path(root)
    from giap.runner import nextflow

    return nextflow.workflow_main().parent.parent


def key_inputs(tool: str, root: Path | None = None) -> list[str]:
    """Repo-relative POSIX paths of every file the image of ``tool`` is built from, sorted."""
    base = _root(root)
    found = [p for p in (base / "workflow" / "containers" / tool).rglob("*") if p.is_file()]
    optional = [
        f"workflow/envs/{tool}.yml",
        f"workflow/envs/locks/{tool}.explicit.txt",
        f"workflow/envs/locks/{tool}.pip.txt",
        *_EXTRA_INPUTS.get(tool, ()),
    ]
    found += [base / rel for rel in optional if (base / rel).is_file()]
    return sorted(p.relative_to(base).as_posix() for p in found)


def key_text(tool: str, root: Path | None = None) -> str:
    """One ``<path>\\t<mode>\\t<sha256>\\n`` line per input (git's two file modes)."""
    base = _root(root)
    inputs = key_inputs(tool, base)
    if not inputs:
        raise ValueError(f"unknown tool {tool!r}: no image inputs under {base}")
    lines = []
    for rel in inputs:
        p = base / rel
        mode = "100755" if p.stat().st_mode & 0o111 else "100644"
        lines.append(f"{rel}\t{mode}\t{hashlib.sha256(p.read_bytes()).hexdigest()}\n")
    return "".join(lines)


def key(tool: str, root: Path | None = None) -> str:
    """The first 12 hex of sha256 over ``key_text``."""
    return hashlib.sha256(key_text(tool, root).encode("utf-8")).hexdigest()[:12]


def tag(version: str, image_key: str) -> str:
    """``<version>-<key>`` with the version made tag-safe; ValueError if no valid tag results."""
    m = _GIT_VERSION.fullmatch(version)
    base = f"git-{m.group(1)[:8]}" if m else re.sub(r"[^A-Za-z0-9_.-]", "_", version)
    out = f"{base}-{image_key}"
    if not _TAG.fullmatch(out):
        raise ValueError(f"not a valid image tag: {out!r}")
    return out


def _copy_sources(dockerfile: str) -> list[str]:
    text = re.sub(r"\\\r?\n", " ", dockerfile)
    sources: list[str] = []
    for line in text.splitlines():
        words = line.strip().split(None, 1)
        if len(words) < 2 or words[0].upper() != "COPY":
            continue
        rest = words[1].strip()
        flags: list[str] = []
        while rest.startswith("--"):
            flag, _, rest = rest.partition(" ")
            flags.append(flag)
            rest = rest.strip()
        if any(f.startswith("--from") for f in flags):
            continue
        args = json.loads(rest) if rest.startswith("[") else shlex.split(rest)
        sources += [str(a) for a in args[:-1]]
    return sources


def uncovered_copy_sources(tool: str, root: Path | None = None) -> list[str]:
    """Build-context ``COPY`` sources of the tool's Dockerfile that are not key inputs.

    A directory source is covered only when every file under it is a key input."""
    base = _root(root)
    inputs = set(key_inputs(tool, base))
    dockerfile = (base / "workflow" / "containers" / tool / "Dockerfile").read_text()
    out = []
    for src in _copy_sources(dockerfile):
        rel = src.removeprefix("./").rstrip("/")
        path = base / rel
        files = [path] if not path.is_dir() else [p for p in path.rglob("*") if p.is_file()]
        rels = {p.relative_to(base).as_posix() for p in files}
        if not rels or not rels <= inputs:
            out.append(src)
    return out


def image_cache(environ: Mapping[str, str], db_dir: Path) -> Path:
    """``NXF_APPTAINER_CACHEDIR``, else ``$GIAP_ROOT/containers``, else ``<db_dir>/containers``."""
    if nxf := environ.get("NXF_APPTAINER_CACHEDIR"):
        return Path(nxf)
    if giap_root := environ.get("GIAP_ROOT"):
        return Path(giap_root) / "containers"
    return Path(db_dir) / "containers"
