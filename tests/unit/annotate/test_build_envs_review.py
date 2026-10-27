"""build_envs.sh against a fake micromamba (M2-01 lead decisions after code review, 1-4 and 7).

No conda is needed: ``micromamba`` and ``git`` are small shell stubs in ``$GIPPA_ROOT/bin`` (which
``farmshare_env.sh`` puts first on PATH). The stub logs its argv, creates env prefixes, prints
per-env version output and can be told to fail a step. The script runs against a throw-away
``GIPPA_REPO`` whose ``src`` and ``scripts`` link to this checkout, so locks and markers land in
``tmp_path``.

Decision 5 (Apptainer cache file name) is not covered here: the naming is inline in
``scripts/slurm/build_envs.sbatch``, not a callable function.
"""

from __future__ import annotations

import hashlib
import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[3]

FAKE_MICROMAMBA = r"""#!/bin/bash
# Fake micromamba for build_envs.sh tests. Controlled by files under $FAKE_DIR.
echo "micromamba $*" >> "$FAKE_DIR/calls.log"
envs="$MAMBA_ROOT_PREFIX/envs"
name_of() { while [ "$#" -gt 0 ]; do if [ "$1" = "-n" ]; then echo "$2"; return; fi; shift; done; }
case "$1" in
  env)
    case "$2" in
      list)
        echo "# conda environments:"
        for d in "$envs"/*/; do [ -d "$d" ] && echo "$(basename "$d")  $d"; done
        exit 0 ;;
      export)
        n="$(name_of "$@")"
        [ -e "$FAKE_DIR/fail_export_$n" ] && exit 1
        echo "@EXPLICIT"
        echo "https://conda.anaconda.org/conda-forge/linux-64/fake-1.0-0.conda#$n"
        exit 0 ;;
      remove)
        n="$(name_of "$@")"; rm -rf "${envs:?}/$n"; exit 0 ;;
    esac
    exit 0 ;;
  create)
    n="$(name_of "$@")"
    [ -e "$FAKE_DIR/fail_create_$n" ] && exit 1
    mkdir -p "$envs/$n/conda-meta" "$envs/$n/bin"
    exit 0 ;;
  remove)
    n="$(name_of "$@")"; rm -rf "${envs:?}/$n"; exit 0 ;;
  run)
    shift
    n=""
    while [ "$#" -gt 0 ]; do
      case "$1" in
        -n) n="$2"; shift 2 ;;
        -*) shift ;;
        *) break ;;
      esac
    done
    [ -d "$envs/$n" ] || { echo "no env $n" >&2; exit 1; }
    case " $* " in
      *" pip "*freeze*|*"pip freeze"*)
        [ -e "$FAKE_DIR/fail_pip_$n" ] && { echo "pip exploded" >&2; exit 1; }
        echo "fakepkg==1.0"
        echo "otherpkg==2.3.4"
        exit 0 ;;
    esac
    [ -e "$FAKE_DIR/fail_version_$n" ] && { echo "version exploded" >&2; exit 1; }
    if [ -e "$FAKE_DIR/version_$n" ]; then cat "$FAKE_DIR/version_$n"; exit 0; fi
    CONDA_PREFIX="$envs/$n" PATH="$envs/$n/bin:$PATH" exec "$@" ;;
  clean) exit 0 ;;
esac
exit 0
"""

FAKE_GIT = r"""#!/bin/bash
echo "git $*" >> "$FAKE_DIR/calls.log"
if [ "$1" = "clone" ]; then
  for a in "$@"; do last="$a"; done
  mkdir -p "$last"
fi
exit 0
"""

# Tools `all` must build (no pinned container, not opt-in) and must skip (decision 4).
ALL_BUILDS = {"uniop", "promotech", "promoter_calculator", "terminatornet", "ostir", "fgsea"}
ALL_SKIPS = {"bakta", "transtermhp", "genomad", "antismash", "miniprot", "blast", "bactermfinder"}


class Sandbox:
    def __init__(self, tmp: Path) -> None:
        self.repo = tmp / "repo"
        self.root = tmp / "root"
        self.fake = tmp / "fake"
        for d in (self.repo / "workflow" / "envs", self.repo / "workflow" / "conf", self.fake):
            d.mkdir(parents=True)
        (self.repo / "src").symlink_to(REPO / "src")
        (self.repo / "scripts").symlink_to(REPO / "scripts")
        for yml in (REPO / "workflow" / "envs").glob("*.yml"):
            shutil.copy(yml, self.repo / "workflow" / "envs" / yml.name)
        # frozen-lock tools (fgsea, M8-17a) build only from a committed lock: here the one the
        # fake micromamba exports, so the post-build lock check passes
        locks = self.repo / "workflow" / "envs" / "locks"
        locks.mkdir()
        (locks / "fgsea.explicit.txt").write_text(
            "@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/fake-1.0-0.conda#gippa-fgsea\n"
        )
        shutil.copy(
            REPO / "workflow" / "conf" / "containers.config",
            self.repo / "workflow" / "conf" / "containers.config",
        )
        bin_dir = self.root / "bin"
        bin_dir.mkdir(parents=True)
        fake_r = "#!/bin/bash\necho '[1] \u20181.36.2\u2019'\n"  # fgsea's version command (M8-17a)
        for name, body in (("micromamba", FAKE_MICROMAMBA), ("git", FAKE_GIT), ("R", fake_r)):
            path = bin_dir / name
            path.write_text(body)
            path.chmod(0o755)
        (self.fake / "calls.log").write_text("")
        self.envs = self.root / "micromamba" / "envs"
        self.locks = self.repo / "workflow" / "envs" / "locks"

    def version(self, tool: str, text: str) -> None:
        (self.fake / f"version_gippa-{tool}").write_text(text + "\n")

    def fail(self, step: str, tool: str) -> None:
        (self.fake / f"fail_{step}_gippa-{tool}").write_text("")

    def run(self, *tools: str) -> subprocess.CompletedProcess[str]:
        env = {k: v for k, v in os.environ.items() if k in {"PATH", "LANG", "TMPDIR"}}
        env.update(
            GIPPA_REPO=str(self.repo),
            GIPPA_ROOT=str(self.root),
            GIPPA_PYTHON=sys.executable,
            FAKE_DIR=str(self.fake),
            HOME=str(self.root / "home"),
        )
        return subprocess.run(
            ["bash", str(self.repo / "scripts" / "build_envs.sh"), *tools],
            capture_output=True,
            text=True,
            cwd=self.repo,
            env=env,
            timeout=300,
        )

    def calls(self) -> list[str]:
        return (self.fake / "calls.log").read_text().splitlines()

    def created(self) -> set[str]:
        out = set()
        for line in self.calls():
            words = line.split()
            if words[:2] == ["micromamba", "create"] and "-n" in words:
                out.add(words[words.index("-n") + 1].removeprefix("gippa-"))
        return out

    def reset_log(self) -> None:
        (self.fake / "calls.log").write_text("")

    def marker(self, tool: str) -> Path:
        return self.envs / f"gippa-{tool}" / ".gippa-built"

    def yml_sha(self, tool: str) -> str:
        data = (self.repo / "workflow" / "envs" / f"{tool}.yml").read_bytes()
        return hashlib.sha256(data).hexdigest()


@pytest.fixture
def sb(tmp_path: Path) -> Sandbox:
    box = Sandbox(tmp_path)
    box.version("bakta", "bakta 1.12.2")
    box.version("miniprot", "0.18")
    box.version("promoter_calculator", "1.2.4")
    box.version("ostir", "1.1.3")
    return box


def _ok(proc: subprocess.CompletedProcess[str]) -> None:
    assert proc.returncode == 0, proc.stdout + proc.stderr


# ---------------------------------------------------------------- decision 1: pip freeze lock


def test_pip_freeze_lock_written(sb: Sandbox) -> None:
    _ok(sb.run("bakta"))
    pip_lock = sb.locks / "bakta.pip.txt"
    assert pip_lock.is_file()
    assert "fakepkg==1.0" in pip_lock.read_text()
    assert any(
        "run" in c.split() and "-n" in c and "gippa-bakta" in c and "pip" in c and "freeze" in c
        for c in sb.calls()
    ), sb.calls()
    assert "@EXPLICIT" in (sb.locks / "bakta.explicit.txt").read_text()


# ---------------------------------------------------------------- decision 2: build marker


def test_marker_holds_yml_sha256(sb: Sandbox) -> None:
    _ok(sb.run("bakta"))
    marker = sb.marker("bakta")
    assert marker.is_file()
    assert sb.yml_sha("bakta") in marker.read_text()


def test_matching_marker_skips_rebuild(sb: Sandbox) -> None:
    _ok(sb.run("bakta"))
    assert sb.created() == {"bakta"}
    sb.reset_log()
    proc = sb.run("bakta")
    _ok(proc)
    assert sb.created() == set()
    assert "gippa-bakta exists" in proc.stdout


def test_missing_marker_triggers_rebuild(sb: Sandbox) -> None:
    _ok(sb.run("bakta"))
    sb.marker("bakta").unlink()  # half-built env: prefix present, no marker
    sb.reset_log()
    _ok(sb.run("bakta"))
    assert sb.created() == {"bakta"}
    assert sb.yml_sha("bakta") in sb.marker("bakta").read_text()


def test_stale_marker_triggers_rebuild(sb: Sandbox) -> None:
    _ok(sb.run("bakta"))
    sb.marker("bakta").write_text("0" * 64 + "\n")
    sb.reset_log()
    _ok(sb.run("bakta"))
    assert sb.created() == {"bakta"}
    assert sb.yml_sha("bakta") in sb.marker("bakta").read_text()


def test_changed_yml_triggers_rebuild(sb: Sandbox) -> None:
    _ok(sb.run("bakta"))
    yml = sb.repo / "workflow" / "envs" / "bakta.yml"
    yml.write_text(yml.read_text() + "# changed\n")
    sb.reset_log()
    _ok(sb.run("bakta"))
    assert sb.created() == {"bakta"}
    assert sb.yml_sha("bakta") in sb.marker("bakta").read_text()


def test_pip_failure_leaves_no_marker(sb: Sandbox) -> None:
    sb.fail("pip", "bakta")
    proc = sb.run("bakta")
    assert proc.returncode != 0
    assert not sb.marker("bakta").exists()
    assert (
        not (sb.locks / "bakta.pip.txt").exists()
        or not (sb.locks / "bakta.pip.txt").read_text().strip()
    )


def test_create_failure_leaves_no_marker(sb: Sandbox) -> None:
    sb.fail("create", "bakta")
    assert sb.run("bakta").returncode != 0
    assert not sb.marker("bakta").exists()


# ---------------------------------------------------------------- decision 3: version step


def _no_locks(sb: Sandbox, tool: str) -> None:
    assert not (sb.locks / f"{tool}.explicit.txt").exists()
    assert not (sb.locks / f"{tool}.pip.txt").exists()


def test_failing_version_command_fails_tool(sb: Sandbox) -> None:
    sb.fail("version", "bakta")
    assert sb.run("bakta").returncode != 0
    _no_locks(sb, "bakta")


def test_version_mismatch_fails_tool(sb: Sandbox) -> None:
    sb.version("bakta", "bakta 1.11.0")
    proc = sb.run("bakta")
    assert proc.returncode != 0
    _no_locks(sb, "bakta")


def test_version_mismatch_fails_only_that_tool(sb: Sandbox) -> None:
    sb.version("miniprot", "0.17")
    proc = sb.run("bakta", "miniprot")
    assert proc.returncode != 0
    _no_locks(sb, "miniprot")
    assert (sb.locks / "bakta.explicit.txt").is_file()
    assert (sb.locks / "bakta.pip.txt").is_file()


def test_git_commit_mismatch_fails_tool(sb: Sandbox) -> None:
    sb.version("uniop", "0" * 40)
    assert sb.run("uniop").returncode != 0
    _no_locks(sb, "uniop")


def test_git_commit_match_passes(sb: Sandbox) -> None:
    from gippa.annotate.tools import TOOLS

    _ok(sb.run("uniop"))  # stub runs the real version_cmd, which reads GIPPA_COMMIT
    sha = TOOLS["uniop"].version.removeprefix("git:")
    assert (sb.envs / "gippa-uniop" / "share" / "uniop" / "GIPPA_COMMIT").read_text().strip() == sha
    assert (sb.locks / "uniop.explicit.txt").is_file()
    assert (sb.locks / "uniop.pip.txt").is_file()


# ---------------------------------------------------------------- decision 4: `all` scope


def test_all_partition_matches_registry() -> None:
    from gippa.annotate.tools import TOOLS, read_containers_config

    pinned = {
        k.lower() for k in read_containers_config(REPO / "workflow" / "conf" / "containers.config")
    }
    opt_in = {n for n, t in TOOLS.items() if t.module == "opt_in"}
    conda = {n for n, t in TOOLS.items() if not t.in_process}
    assert conda - pinned - opt_in == ALL_BUILDS
    assert (pinned | opt_in) & conda == ALL_SKIPS


def test_all_skips_container_and_opt_in_tools(sb: Sandbox) -> None:
    _ok(sb.run("all"))
    assert sb.created() == ALL_BUILDS
    for tool in ALL_SKIPS:
        assert not (sb.envs / f"gippa-{tool}").exists(), tool
        assert not (sb.locks / f"{tool}.explicit.txt").exists(), tool


@pytest.mark.parametrize("tool", ["bakta", "bactermfinder"])
def test_explicit_name_builds_skipped_tool(sb: Sandbox, tool: str) -> None:
    _ok(sb.run(tool))
    assert sb.created() == {tool}
    assert sb.marker(tool).is_file()


# ---------------------------------------------------------------- decision 7: comment wording


def test_farmshare_env_comment_says_gippa_root_must_be_exported() -> None:
    lines = (REPO / "scripts" / "farmshare_env.sh").read_text().splitlines()
    comments = [ln.split("#", 1)[1] for ln in lines if "#" in ln]
    hits = [
        c
        for c in comments
        if "worktree" in c.lower()
        and "GIPPA_ROOT" in c
        and "export" in c.lower()
        and "only" in c.lower()
    ]
    assert hits, "farmshare_env.sh must say worktrees share local/ only when GIPPA_ROOT is exported"
    assert not any(
        "worktrees share the main checkout's local/" in c and "only" not in c for c in comments
    )
