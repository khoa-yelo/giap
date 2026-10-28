"""Fake container runtimes that refuse network pulls in tests (``tests/conftest.py``).

Unit tests of ``gippa setup`` once ran the host's real ``apptainer pull`` (about 6 GB per image)
and filled ``/tmp`` (fixed one by one on #437 and #441). ``gippa setup`` and ``gippa run`` find
the runtimes on ``PATH`` only: ``subprocess.run(["apptainer", ...])``/``["docker", ...]`` in
``gippa.setupdb.images``, ``shutil.which("apptainer") or shutil.which("singularity")`` and
``shutil.which("docker")`` in ``gippa.setupdb.uniop_build`` and ``gippa.preflight.environment``;
no environment variable names a runtime binary. So the autouse fixture ``_no_container_pulls``
puts the scripts written here first on ``PATH``, one per runtime the host has.

Each script refuses (exit ``REFUSED_EXIT``, ``REFUSED`` on stderr, the real binary never runs):

- ``apptainer``/``singularity``: ``pull``; any argument that is a remote image URI
  (``docker://``, ``library://``, ``oras://``, ``shub://``, ``http(s)://``), so ``exec``/``run``
  of a remote image too; ``build`` of a definition file whose ``Bootstrap:`` is not local
  (``localimage``, ``scratch``, ``oci-archive``, ``docker-archive``).
- ``docker``: ``pull``, ``push``, ``build``, ``buildx``, ``image pull|push|build``,
  ``builder build``, ``compose pull|push|build|up|create|run``; ``run``/``create`` (also as
  ``container run|create``) unless ``--pull=never`` / ``--pull never``.

Anything else runs the real binary found on ``PATH`` at session start (capability probes such as
``apptainer --version`` keep working). A runtime the host lacks gets no script, so
``shutil.which`` answers as without the guard. A test's own fake, prepended to ``PATH`` after
the autouse fixture ran (or a monkeypatched ``subprocess.run``), wins over these.
"""

from __future__ import annotations

import os
import shlex
import shutil
from collections.abc import Mapping
from pathlib import Path

GUARD_ENV = "GIPPA_TEST_CONTAINER_GUARD"
NAMES = ("apptainer", "singularity", "docker")
REFUSED = "network pull refused in tests"
REFUSED_EXIT = 97
OPT_OUT = "real_container_pull"

_HEAD = """#!/bin/sh
# Written by tests/helpers/container_guard.py (tests/conftest.py _no_container_pulls).
name={name}
real={real}
refuse() {{
    echo "gippa test guard: {refused}: $name $* (tests/conftest.py _no_container_pulls; put \
your own fake first on PATH or mark the test {opt_out})" >&2
    exit {code}
}}
remote() {{
    case "$1" in
    docker://*|library://*|oras://*|shub://*|http://*|https://*) return 0;;
    esac
    return 1
}}
"""

_APPTAINER = """
sub=""
for a in "$@"; do
    case "$a" in -*) ;; *) sub=$a; break;; esac
done
[ "$sub" = pull ] && refuse "$@"
for a in "$@"; do
    remote "$a" && refuse "$@"
done
if [ "$sub" = build ]; then
    for a in "$@"; do
        [ -f "$a" ] || continue
        boot=$(head -c 65536 "$a" 2>/dev/null | LC_ALL=C tr -d '\\000' \\
            | LC_ALL=C sed -n 's/^[[:space:]]*[Bb][Oo][Oo][Tt][Ss][Tt][Rr][Aa][Pp][[:space:]]*:\
[[:space:]]*\\([A-Za-z_-]*\\).*/\\1/p' | head -n 1 | tr 'A-Z' 'a-z')
        case "$boot" in
        ""|localimage|scratch|oci-archive|docker-archive) ;;
        *) refuse "$@";;
        esac
    done
fi
"""

_DOCKER = """
sub=""; sub2=""; skip=0
for a in "$@"; do
    if [ "$skip" = 1 ]; then skip=0; continue; fi
    case "$a" in
    -H|--host|-c|--context|--config|-l|--log-level|--tlscacert|--tlscert|--tlskey) skip=1;;
    -*) ;;
    *) if [ -z "$sub" ]; then sub=$a; else sub2=$a; break; fi;;
    esac
done
never=0; prev=""
for a in "$@"; do
    [ "$a" = --pull=never ] && never=1
    [ "$prev" = --pull ] && [ "$a" = never ] && never=1
    prev=$a
done
case "$sub" in
pull|push|build|buildx) refuse "$@";;
image) case "$sub2" in pull|push|build) refuse "$@";; esac;;
builder) [ "$sub2" = build ] && refuse "$@";;
compose) case "$sub2" in pull|push|build|up|create|run) refuse "$@";; esac;;
run|create) [ "$never" = 1 ] || refuse "$@";;
container) case "$sub2" in run|create) [ "$never" = 1 ] || refuse "$@";; esac;;
esac
"""

_TAIL = """
if [ -z "$real" ]; then
    echo "gippa test guard: no real $name on PATH" >&2
    exit 127
fi
exec "$real" "$@"
"""


def guard_script(name: str, real: str | None) -> str:
    """The ``sh`` script shadowing ``name``; ``real`` is the binary it hands safe calls to."""
    head = _HEAD.format(name=shlex.quote(name), real=shlex.quote(real or ""), refused=REFUSED,
                        opt_out=OPT_OUT, code=REFUSED_EXIT)  # fmt: skip
    body = _DOCKER if name == "docker" else _APPTAINER
    return head + body + _TAIL


def real_binaries(path: str, skip: tuple[str, ...] = ()) -> dict[str, str]:
    """``{name: absolute path}`` of the runtimes on ``path``, ignoring the ``skip`` directories
    (an outer guard, so a child pytest never wraps a guard in a guard)."""
    dirs = [d for d in path.split(os.pathsep) if d]
    if skip:
        gone = {Path(s).resolve() for s in skip}
        dirs = [d for d in dirs if Path(d).resolve() not in gone]
    search = os.pathsep.join(dirs)
    return {n: os.path.abspath(p) for n in NAMES if (p := shutil.which(n, path=search))}


def write_guard(bindir: Path, reals: Mapping[str, str | None]) -> Path:
    """Write one guard script per name of ``reals`` into ``bindir`` (mode 0755)."""
    bindir.mkdir(parents=True, exist_ok=True)
    for name, real in reals.items():
        exe = bindir / name
        exe.write_text(guard_script(name, real), encoding="utf-8")
        exe.chmod(0o755)
    return bindir
