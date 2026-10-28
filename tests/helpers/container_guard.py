"""Fake container runtimes that refuse network pulls in tests (``tests/conftest.py``).

Unit tests of ``gippa setup`` once ran the host's real ``apptainer pull`` (about 6 GB per image)
and filled ``/tmp`` (fixed one by one on #437 and #441). ``gippa setup`` and ``gippa run`` find
the runtimes on ``PATH`` only: ``subprocess.run(["apptainer", ...])``/``["docker", ...]`` in
``gippa.setupdb.images``, ``shutil.which("apptainer") or shutil.which("singularity")`` and
``shutil.which("docker")`` in ``gippa.setupdb.uniop_build`` and ``gippa.preflight.environment``;
no environment variable names a runtime binary. So the autouse fixture ``_no_container_pulls``
puts the scripts written here first on ``PATH``, one per runtime the host has.

Each script refuses (exit ``REFUSED_EXIT``, ``REFUSED`` on stderr, the real binary never runs):

- ``apptainer``/``singularity`` (the subcommand is the first word that is neither a flag nor
  the value of ``-c``/``--config``):
  - ``pull``;
  - ``build`` with any remote image URI argument, a definition on stdin (``build out.sif -``:
    stdin cannot be inspected without consuming it), or a definition file whose ``Bootstrap:``
    is not local (``localimage``, ``scratch``, ``oci-archive``, ``docker-archive``);
  - any other subcommand (``exec``/``run``/``shell``/``instance start``/...) whose image is a
    remote URI. Remote schemes (``docker://``, ``library://``, ``oras://``, ``shub://``,
    ``http(s)://``) match in any case (``DOCKER://``). Only the words up to the image are
    scanned: the image is the first positional that is a local path or a non-remote URI
    (``instance://``), and the words after it are the container's own command, so
    ``exec local.sif curl https://...`` runs. Values of the known value flags (``--env``,
    ``--bind``, ``--home``, ...; ``_APPTAINER_VALUE_FLAGS``) are never scanned or taken for the
    image: ``--env URL=https://...`` is an argument, not a pull. When no word is a local path
    (the image does not exist) every word is scanned, so ``exec missing.sif curl https://...``
    is refused: a false refusal, never a pull.
- ``docker``: ``pull``, ``push``, ``build``, ``buildx``, ``image pull|push|build``,
  ``builder build``, ``import``/``image import`` of a URL, ``plugin install|upgrade|push``,
  ``service create|update``, ``stack deploy|up``, ``compose`` with any later word
  ``pull|push|build|up|create|run``; ``run``/``create`` (also as ``container run|create``)
  unless ``--pull=never``/``--pull never`` comes among the options before the image (the
  container's own arguments do not count; ``_DOCKER_RUN_VALUE_FLAGS`` names the flags whose
  value is skipped, and an unknown flag's value is taken for the image: a refusal again).

Anything else runs the real binary found on ``PATH`` at session start (capability probes such as
``apptainer --version`` keep working). A runtime the host lacks gets no script, so
``shutil.which`` answers as without the guard. A test's own fake, prepended to ``PATH`` after
the autouse fixture ran (or a monkeypatched ``subprocess.run``), wins over these.

Known limits:

- A subprocess whose environment names its own ``PATH`` skips the guard. Helpers that build
  such an environment call ``guarded_path`` (``tests/helpers/pilot_runs_m10_04.py``,
  ``tests/helpers/a3_run_m10_10.py``): with the guard active (``GUARD_ENV`` set) it prepends
  guard scripts for the runtimes *that* ``PATH`` holds, forwarding to them, so the subprocess
  finds the same runtimes as without the guard. Some unit tests still set ``PATH=/usr/bin:/bin``
  inline for scripts that never start a runtime.
- An unknown apptainer flag whose value is an existing path ends the scan early, so a remote
  image after it passes; a runtime path resolved at import or collection time skips the guard.
- Other registry traffic is not refused (``docker search``, ``manifest``, ``login``,
  ``apptainer remote``, ``apptainer push``): none of it downloads an image.
"""

from __future__ import annotations

import hashlib
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
remote() {{  # a remote image URI; schemes are case-insensitive (DOCKER://, Https://)
    case "$1" in *://*) ;; *) return 1;; esac
    low=$(printf '%s' "$1" | LC_ALL=C tr 'A-Z' 'a-z')
    case "$low" in
    docker://*|library://*|oras://*|shub://*|http://*|https://*) return 0;;
    esac
    return 1
}}
"""

# apptainer exec/run/shell/... flags that take a value: their value is never the image and
# never scanned (``--env URL=https://...`` is an argument, not a pull)
_APPTAINER_VALUE_FLAGS = "|".join([
    "-B", "--bind", "--mount", "-H", "--home", "--env", "--env-file", "-W", "--workdir", "-S",
    "--scratch", "--pwd", "--security", "--network", "--network-args", "--hostname", "-o",
    "--overlay", "--app", "--add-caps", "--drop-caps", "--dns", "--fusemount", "--apply-cgroups",
    "--cpus", "--cpu-shares", "--cpuset-cpus", "--cpuset-mems", "--memory",
    "--memory-reservation", "--memory-swap", "--pids-limit", "--blkio-weight", "--device",
    "--pid-file", "--passphrase", "--pem-path", "--authfile",
])  # fmt: skip

_APPTAINER = """
sub=""; skip=0
for a in "$@"; do
    if [ "$skip" = 1 ]; then skip=0; continue; fi
    case "$a" in -c|--config) skip=1;; -*) ;; *) sub=$a; break;; esac
done
[ "$sub" = pull ] && refuse "$@"
if [ "$sub" = build ]; then
    for a in "$@"; do
        [ "$a" = - ] && refuse "$@"  # definition on stdin: not inspectable without consuming it
        remote "$a" && refuse "$@"
        [ -f "$a" ] || continue
        boot=$(head -c 65536 "$a" 2>/dev/null | LC_ALL=C tr -d '\\000' \\
            | LC_ALL=C sed -n 's/^[[:space:]]*[Bb][Oo][Oo][Tt][Ss][Tt][Rr][Aa][Pp][[:space:]]*:\
[[:space:]]*\\([A-Za-z_-]*\\).*/\\1/p' | head -n 1 | tr 'A-Z' 'a-z')
        case "$boot" in
        ""|localimage|scratch|oci-archive|docker-archive) ;;
        *) refuse "$@";;
        esac
    done
else
    # exec/run/shell/instance start/...: scan the words up to the image, the first positional
    # that is a local path or a non-remote URI (instance://, docker-daemon://); the words after
    # it are the container's own command (curl https://... is not a pull)
    seen=0; skip=0
    for a in "$@"; do
        if [ "$skip" = 1 ]; then skip=0; continue; fi
        if [ "$seen" = 0 ]; then
            case "$a" in -c|--config) skip=1;; -*) ;; *) seen=1;; esac
            continue
        fi
        case "$a" in
        @VALUE_FLAGS@) skip=1;;
        -*) ;;
        *) remote "$a" && refuse "$@"
           case "$a" in *://*) break;; esac
           [ -e "$a" ] && break;;
        esac
    done
fi
""".replace("@VALUE_FLAGS@", _APPTAINER_VALUE_FLAGS)

_DOCKER_GLOBAL_VALUE_FLAGS = (
    "-H|--host|-c|--context|--config|-l|--log-level|--tlscacert|--tlscert|--tlskey"
)
# docker run/create flags that take a value (a value never counts as the image); an unknown
# flag's value is taken for the image, so a later --pull=never is missed: a refusal, never a pull
_DOCKER_RUN_VALUE_FLAGS = "|".join([
    "-a", "--attach", "-c", "--cpu-shares", "-e", "--env", "-h", "--hostname", "-l", "--label",
    "-m", "--memory", "-p", "--publish", "-u", "--user", "-v", "--volume", "-w", "--workdir",
    "--add-host", "--annotation", "--blkio-weight", "--cap-add", "--cap-drop", "--cgroup-parent",
    "--cgroupns", "--cidfile", "--cpu-period", "--cpu-quota", "--cpus", "--cpuset-cpus",
    "--cpuset-mems", "--detach-keys", "--device", "--dns", "--dns-option", "--dns-search",
    "--domainname", "--entrypoint", "--env-file", "--expose", "--gpus", "--group-add",
    "--health-cmd", "--health-interval", "--health-retries", "--health-timeout", "--ip", "--ip6",
    "--ipc", "--isolation", "--label-file", "--link", "--log-driver", "--log-opt",
    "--mac-address", "--memory-reservation", "--memory-swap", "--mount", "--name", "--network",
    "--net", "--network-alias", "--pid", "--platform", "--restart", "--runtime",
    "--security-opt", "--shm-size", "--stop-signal", "--stop-timeout", "--storage-opt",
    "--sysctl", "--tmpfs", "--ulimit", "--userns", "--uts", "--volumes-from",
])  # fmt: skip

_DOCKER = """
sub=""; sub2=""; skip=0
for a in "$@"; do
    if [ "$skip" = 1 ]; then skip=0; continue; fi
    case "$a" in
    @GLOBAL@) skip=1;;
    -*) ;;
    *) if [ -z "$sub" ]; then sub=$a; else sub2=$a; break; fi;;
    esac
done
# --pull=never / --pull never among the run/create options only (before the image), never in
# the container's own command
never=0; phase=0; skip=""
for a in "$@"; do
    if [ -n "$skip" ]; then
        [ "$skip" = pull ] && [ "$a" = never ] && never=1
        skip=""; continue
    fi
    if [ "$phase" = 0 ]; then
        case "$a" in @GLOBAL@) skip=1;; -*) ;; run|create) phase=1;; esac
        continue
    fi
    case "$a" in
    --pull=never) never=1;;
    --pull) skip=pull;;
    @RUN@) skip=1;;
    -*) ;;
    *) break;;
    esac
done
any_remote() {
    for a in "$@"; do remote "$a" && return 0; done
    return 1
}
case "$sub" in
pull|push|build|buildx) refuse "$@";;
image) case "$sub2" in
    pull|push|build) refuse "$@";;
    import) any_remote "$@" && refuse "$@";;
    esac;;
import) any_remote "$@" && refuse "$@";;
builder) [ "$sub2" = build ] && refuse "$@";;
plugin) case "$sub2" in install|upgrade|push) refuse "$@";; esac;;
service) case "$sub2" in create|update) refuse "$@";; esac;;
stack) case "$sub2" in deploy|up) refuse "$@";; esac;;
compose)  # any later word: compose's value flags (-f c.yml, -p name, ...) hide the subcommand
    for a in "$@"; do
        case "$a" in pull|push|build|up|create|run) refuse "$@";; esac
    done;;
run|create) [ "$never" = 1 ] || refuse "$@";;
container) case "$sub2" in run|create) [ "$never" = 1 ] || refuse "$@";; esac;;
esac
""".replace("@GLOBAL@", _DOCKER_GLOBAL_VALUE_FLAGS).replace("@RUN@", _DOCKER_RUN_VALUE_FLAGS)

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


def guarded_path(path: str, environ: Mapping[str, str] | None = None) -> str:
    """``path`` for a subprocess environment, with guard scripts first while the guard is active.

    Without ``GUARD_ENV`` in ``environ`` (default ``os.environ``; an opted-out test or no
    runtime on the host) or when ``path`` holds no runtime, ``path`` is returned unchanged. Else
    scripts for exactly the runtimes ``path`` holds are written once under the guard directory
    (they forward to ``path``'s binaries) and their directory is put first."""
    env = os.environ if environ is None else environ
    guard = env.get(GUARD_ENV)
    if not guard or not Path(guard).is_dir():
        return path
    reals = real_binaries(path, (guard,))
    if not reals:
        return path
    key = hashlib.sha256(repr(sorted(reals.items())).encode()).hexdigest()[:12]
    bindir = Path(guard) / f"path-{key}"
    if not all((bindir / n).is_file() for n in reals):
        write_guard(bindir, reals)
    return f"{bindir}{os.pathsep}{path}"
