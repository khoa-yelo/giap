"""Fake Slurm / GitHub executables and a sandbox GIAP checkout for the M2-21 shell-script tests.

Nothing here reaches Slurm, GitHub, the network, the real ``$GIAP_ROOT`` or the real repository:

* ``Sandbox.repo`` is a fresh git repository holding a copy of this checkout's ``scripts/`` and
  ``ci/`` (``scripts/farmshare_env.sh`` replaced by a stub that only exports sandbox paths), with
  a bare ``origin`` next to it, so ``git fetch origin`` / ``worktree add origin/main`` stay local.
* ``Sandbox.bin`` is first on a minimal ``PATH`` (``bin:/usr/bin:/bin``) and holds fakes for
  ``sbatch squeue scontrol scancel srun sacct gh hostname nf-test giap nextflow`` and a ``python``
  wrapper that fakes ``python -m pytest`` and runs everything else with the real interpreter.
  Every fake logs its argv (JSON lines under ``Sandbox.state``).
* The scripts are run only after ``assert_honours_giap_repo`` found no hard-coded path into the
  real repository outside ``${GIAP_REPO:-…}`` defaults and ``--output``/``--error`` arguments.
"""

from __future__ import annotations

import json
import os
import re
import shutil
import stat
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
REAL_REPO = "/home/classes/bios/270/khoa/giap"
SLURM_TOOLS = ("sbatch", "squeue", "scontrol", "scancel", "srun", "sacct", "scrontab")

_SBATCH = r"""#!{python}
import json, os, sys
state = os.environ["FAKE_SLURM_DIR"]
log = os.path.join(state, "sbatch.jsonl")
counter = os.path.join(state, "next_id")
n = int(open(counter).read()) if os.path.exists(counter) else 4242
calls = sum(1 for _ in open(log)) if os.path.exists(log) else 0
env = {k: v for k, v in os.environ.items()
       if k.startswith(("NIGHTLY_", "WEEKLY_", "GIAP_")) or k == "MODE"}
script = None
for a in sys.argv[1:]:
    if not a.startswith("-") and os.path.isfile(a):
        script = open(a).read(); break
rec = {"argv": sys.argv[1:], "env": env, "cwd": os.getcwd(), "id": n, "script": script}
fail = os.environ.get("FAKE_SBATCH_FAIL", "")
failed = fail == "all" or str(calls) in fail.split(",")
rec["failed"] = failed
with open(log, "a") as fh:
    fh.write(json.dumps(rec) + "\n")
if failed:
    print("sbatch: error: Batch job submission failed: fake failure", file=sys.stderr)
    sys.exit(1)
open(counter, "w").write(str(n + 1))
print(n if any(a.startswith("--parsable") for a in sys.argv[1:]) else f"Submitted batch job {n}")
"""

_SACCT = r"""#!{python}
import json, os, sys
state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
with open(os.path.join(state, "sacct.jsonl"), "a") as fh:
    fh.write(json.dumps(argv) + "\n")
rc = int(os.environ.get("FAKE_SACCT_RC", "0"))
if rc:
    print("sacct: error: fake failure", file=sys.stderr); sys.exit(rc)
path = os.path.join(state, "sacct_rows.json")
rows = json.load(open(path)) if os.path.exists(path) else []
fields = ["JobID", "JobName", "State", "ExitCode"]
jobs, parsable, noheader, trailing = [], False, False, False
i = 0
while i < len(argv):
    a = argv[i]
    if a in ("--format", "-o"):
        fields = argv[i + 1].split(","); i += 1
    elif a.startswith("--format="):
        fields = a.split("=", 1)[1].split(",")
    elif a in ("-j", "--jobs"):
        jobs += argv[i + 1].split(","); i += 1
    elif a.startswith("--jobs="):
        jobs += a.split("=", 1)[1].split(",")
    elif a.startswith("-j") and len(a) > 2:
        jobs += a[2:].split(",")
    elif a in ("-P", "--parsable2"):
        parsable = True
    elif a in ("-p", "--parsable"):
        parsable = trailing = True
    elif a in ("-n", "--noheader"):
        noheader = True
    elif a.startswith("-") and not a.startswith("--") and set(a[1:]) <= set("nPpXa"):
        parsable = parsable or "P" in a or "p" in a
        trailing = trailing or "p" in a
        noheader = noheader or "n" in a
    i += 1
fields = [f.split("%")[0] for f in fields]
def match(jid):
    return not jobs or any(jid == j or jid.startswith(j + "_") or jid.startswith(j + ".")
                           for j in jobs)
sep = "|" if parsable else "  "
if not noheader:
    print(sep.join(fields))
for r in rows:
    if match(r.get("JobID", "")):
        print(sep.join(str(r.get(f, "")) for f in fields) + ("|" if trailing else ""))
"""

_GH = r"""#!{python}
import json, os, sys
state = os.environ["FAKE_SLURM_DIR"]
argv = sys.argv[1:]
spath = os.path.join(state, "gh_state.json")
st = json.load(open(spath)) if os.path.exists(spath) else {"issues": [], "next": 77}
stdin = None
if "--body-file" in argv and argv[argv.index("--body-file") + 1] == "-":
    stdin = sys.stdin.read()
elif "--body-file" in argv:
    stdin = open(argv[argv.index("--body-file") + 1]).read()
elif "--body" in argv:
    stdin = argv[argv.index("--body") + 1]
with open(os.path.join(state, "gh.jsonl"), "a") as fh:
    fh.write(json.dumps({"argv": argv, "body": stdin}) + "\n")
fail = os.environ.get("FAKE_GH_FAIL", "")
if argv[:2] == ["issue", "list"]:
    print(json.dumps([{"number": i["number"], "title": i["title"]} for i in st["issues"]]))
elif argv[:2] == ["issue", "create"]:
    if fail == "create":
        print("gh: fake failure", file=sys.stderr); sys.exit(1)
    title = argv[argv.index("--title") + 1] if "--title" in argv else ""
    st["issues"].append({"number": st["next"], "title": title})
    print(f"https://github.invalid/khoa-yelo/giap/issues/{st['next']}")
    st["next"] += 1
elif argv[:2] == ["issue", "comment"]:
    if fail == "comment":
        print("gh: fake failure", file=sys.stderr); sys.exit(1)
    print(f"https://github.invalid/khoa-yelo/giap/issues/{argv[2]}#issuecomment-1")
json.dump(st, open(spath, "w"))
"""

_LOGGER = r"""#!/bin/bash
printf '%s\n' "$(basename "$0") $*" >> "$FAKE_SLURM_DIR/calls.log"
exit 0
"""

_HOSTNAME = r"""#!/bin/bash
h="${FAKE_HOSTNAME:-barley-01}"
case "${1:-}" in -s|--short) echo "${h%%.*}" ;; *) echo "$h" ;; esac
"""

_PYTHON = r"""#!/bin/bash
if [ "${1:-}" = "-m" ] && [ "${2:-}" = "pytest" ]; then
  printf '%s\n' "$*" >> "$FAKE_SLURM_DIR/pytest.argv"
  case "${FAKE_PYTEST:-ok}" in
    term) kill -TERM 0; sleep 5; echo "1 failed in 5.00s"; exit 1 ;;
    fail) echo "1 failed in 0.01s"; exit 1 ;;
    *) echo "1 passed in 0.01s"; exit 0 ;;
  esac
fi
exec {python} "$@"
"""

_ENV_STUB = """# Sandbox stand-in for scripts/farmshare_env.sh (tests/helpers/slurm_sandbox.py).
: "${GIAP_ROOT:?sandbox GIAP_ROOT}"
export GIAP_REPO GIAP_ROOT
export UV_PROJECT_ENVIRONMENT="$GIAP_ROOT/venv"
export GIAP_DB_DIR="$GIAP_ROOT/db" GIAP_DATA_DIR="$GIAP_ROOT/data" GIAP_LOG_DIR="$GIAP_ROOT/logs"
"""


FAKE_GH = _GH


def write_exe(path: Path, text: str) -> Path:
    path.write_text(text)
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


def assert_honours_giap_repo(text: str, what: str) -> None:
    """Refuse to run a script that would touch the real repository from the sandbox."""
    t = re.sub(r"\$\{GIAP_REPO:-[^}]*\}", "", text)
    t = re.sub(r"--(output|error)[= ]\S+", "", t)
    bad = [ln for ln in t.splitlines()
           if REAL_REPO in ln and not ln.lstrip().startswith("#")]  # fmt: skip
    assert not bad, (
        f"{what}: hard-coded path into the real repository (use ${{GIAP_REPO:-…}}): {bad}"
    )
    assert not re.search(r"/usr/bin/(sbatch|sacct|squeue|scontrol|scancel|scrontab)\b", t), what


@dataclass
class Sandbox:
    tmp: Path
    repo: Path
    root: Path
    bin: Path
    state: Path
    env: dict[str, str] = field(default_factory=dict)

    @property
    def head(self) -> str:
        cmd = ["git", "-C", str(self.repo), "rev-parse", "HEAD"]
        return subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.strip()

    def jsonl(self, name: str) -> list[Any]:
        p = self.state / name
        return [json.loads(ln) for ln in p.read_text().splitlines()] if p.is_file() else []

    def sbatch_calls(self) -> list[dict[str, Any]]:
        return self.jsonl("sbatch.jsonl")

    def gh_calls(self) -> list[dict[str, Any]]:
        return self.jsonl("gh.jsonl")

    def set_sacct(self, rows: list[dict[str, str]]) -> None:
        (self.state / "sacct_rows.json").write_text(json.dumps(rows))

    def tsv_rows(self, name: str = "nightly.tsv") -> list[list[str]]:
        p = self.root / "logs" / name
        if not p.is_file():
            return []
        return [ln.split("\t") for ln in p.read_text().splitlines() if ln.strip()]

    def run(self, cmd: list[str], *, env: dict[str, str] | None = None, cwd: Path | None = None,
            timeout: float = 180) -> subprocess.CompletedProcess[str]:  # fmt: skip
        return subprocess.run(
            cmd,
            env={**self.env, **(env or {})},
            cwd=str(cwd or self.tmp),
            capture_output=True,
            text=True,
            timeout=timeout,
            start_new_session=True,  # `kill -TERM 0` in the fake pytest hits only the script
            check=False,
        )


def _git(*args: str, cwd: Path, env: dict[str, str]) -> None:
    subprocess.run(["git", *args], cwd=str(cwd), env=env, check=True, capture_output=True)


def make_sandbox(tmp: Path) -> Sandbox:
    repo, root, state = tmp / "repo", tmp / "root", tmp / "state"
    home = tmp / "home"
    for d in (repo, root / "logs", state, home):
        d.mkdir(parents=True)
    ignore = shutil.ignore_patterns("__pycache__", "*.pyc")
    shutil.copytree(REPO / "scripts", repo / "scripts", ignore=ignore)
    shutil.copytree(REPO / "ci", repo / "ci", ignore=ignore)
    (repo / "scripts" / "farmshare_env.sh").write_text(_ENV_STUB)
    bindir = root / "venv" / "bin"
    bindir.mkdir(parents=True)
    py = sys.executable
    write_exe(bindir / "sbatch", _SBATCH.replace("{python}", py))
    write_exe(bindir / "sacct", _SACCT.replace("{python}", py))
    write_exe(bindir / "gh", _GH.replace("{python}", py))
    write_exe(bindir / "hostname", _HOSTNAME)
    write_exe(bindir / "python", _PYTHON.replace("{python}", py))
    write_exe(bindir / "python3", _PYTHON.replace("{python}", py))
    for name in ("squeue", "scontrol", "scancel", "srun", "scrontab", "nf-test", "giap",
                 "nextflow", "module"):  # fmt: skip
        write_exe(bindir / name, _LOGGER)
    env = {
        "PATH": f"{bindir}:/usr/bin:/bin",
        "HOME": str(home),
        "LANG": "C.UTF-8",
        "USER": os.environ.get("USER", "giap"),
        "GIAP_REPO": str(repo),
        "GIAP_ROOT": str(root),
        "FAKE_SLURM_DIR": str(state),
        "GIT_AUTHOR_NAME": "sandbox",
        "GIT_AUTHOR_EMAIL": "sandbox@example.invalid",
        "GIT_COMMITTER_NAME": "sandbox",
        "GIT_COMMITTER_EMAIL": "sandbox@example.invalid",
        "GIT_CONFIG_GLOBAL": "/dev/null",
        "GIT_CONFIG_NOSYSTEM": "1",
        "GH_TOKEN": "fake-token-never-used",
        "GH_HOST": "github.invalid",
    }
    _git("init", "-q", "-b", "main", cwd=repo, env=env)
    _git("add", "-A", cwd=repo, env=env)
    _git("commit", "-q", "-m", "sandbox", cwd=repo, env=env)
    _git("clone", "-q", "--bare", str(repo), str(tmp / "origin.git"), cwd=tmp, env=env)
    _git("remote", "add", "origin", str(tmp / "origin.git"), cwd=repo, env=env)
    _git("fetch", "-q", "origin", cwd=repo, env=env)
    return Sandbox(tmp=tmp, repo=repo, root=root, bin=bindir, state=state, env=env)


def export_pairs(argv: list[str]) -> dict[str, str]:
    """``--export=ALL,K=V,…`` (or ``--export K=V``) of one sbatch call → ``{K: V}``."""
    out: dict[str, str] = {}
    for i, a in enumerate(argv):
        val = None
        if a.startswith("--export="):
            val = a.split("=", 1)[1]
        elif a == "--export" and i + 1 < len(argv):
            val = argv[i + 1]
        if val is None:
            continue
        for part in val.split(","):
            k, sep, v = part.partition("=")
            if sep:
                out[k] = v
    return out


def option(argv: list[str], long: str, short: str | None = None) -> str | None:
    """Value of ``--long=V`` / ``--long V`` / ``-s V`` / ``-sV`` in one argv (last wins)."""
    val = None
    for i, a in enumerate(argv):
        if a.startswith(f"--{long}="):
            val = a.split("=", 1)[1]
        elif i + 1 < len(argv) and (a == f"--{long}" or (short and a == f"-{short}")):
            val = argv[i + 1]
        elif short and a.startswith(f"-{short}") and len(a) > 2 and not a.startswith("--"):
            val = a[2:]
    return val
