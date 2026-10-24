"""Env files, ``PILOT_TOOLS`` and the ``build_envs.sh`` lookup for the A3 pilot (M10-08a).

Spec: docs/specs/M10-08-a3-feasibility.md, deliverables 1-3, acceptance criterion 10, lead
decisions Q4 (earlier: not in ``TOOLS``), Q10 (pyteomics pinned to 4.7.5) and the Cascadia 0.0.7
facts table (PyPI only, Apache-2.0, wheel sha256 ``50b52043…f092a3``, Python 3.10,
PyTorch 2.0.x with CUDA 11.8). Route B2: TIMSCONVERT at commit ``3cb96aef…`` (tag v2.0.0).

``build_envs.sh`` runs against a throw-away ``GIAP_REPO`` whose ``src``/``scripts`` link to this
checkout and a fake ``micromamba`` (bash 3.2: no mapfile, no associative arrays), so nothing is
solved, cloned or downloaded. No statistical claim, so no calibration test (spec §8.4).
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
import yaml
from helpers.a3_kit import (
    CASCADIA_VERSION,
    CASCADIA_WHEEL_SHA_HEAD,
    CASCADIA_WHEEL_SHA_TAIL,
    ENVS,
    TIMSCONVERT_COMMIT,
    require,
)
from helpers.slurm_sandbox import REPO

PILOT = ("cascadia", "timsconvert")


def _yml(stem: str) -> dict[str, object]:
    return yaml.safe_load(require(ENVS / f"{stem}.yml", "deliverables 1-2").read_text())


def _header(stem: str) -> str:
    text = require(ENVS / f"{stem}.yml", "deliverables 1-2").read_text()
    return "\n".join(ln for ln in text.splitlines() if ln.lstrip().startswith("#"))


def _deps(stem: str) -> tuple[list[str], list[str]]:
    conda: list[str] = []
    pip: list[str] = []
    for d in _yml(stem).get("dependencies") or []:  # type: ignore[union-attr]
        if isinstance(d, dict):
            pip += [str(p) for p in d.get("pip") or []]
        else:
            conda.append(str(d).replace(" ", ""))
    return conda, pip


def _conda(stem: str, pkg: str) -> list[str]:
    pat = re.compile(rf"^{re.escape(pkg)}(?=[=<>!~]|$)", re.I)
    return [d for d in _deps(stem)[0] if pat.match(d)]


# ------------------------------------------------------------------ deliverable 1: cascadia.yml


def test_cascadia_env_name_and_channels() -> None:
    data = _yml("cascadia")
    assert data["name"] == "giap-cascadia"
    assert data["channels"] == ["conda-forge"], "conda-forge only (deliverable 1)"


@pytest.mark.parametrize(
    ("pkg", "pin"),
    [
        ("python", r"=3\.10(\.|=|\*|$)"),
        ("pytorch", r"=2\.0(\.|=|\*|$)"),
        ("pyteomics", r"=4\.7\.5$"),  # Q10: pinned, not a range
        ("psims", r"=1\.4\.0$"),
    ],
)
def test_cascadia_conda_pins(pkg: str, pin: str) -> None:
    hits = _conda("cascadia", pkg)
    assert len(hits) == 1, (pkg, _deps("cascadia")[0])
    assert re.search(pin, hits[0][len(pkg) :]), hits[0]


def test_cascadia_numpy_below_2() -> None:
    (dep,) = _conda("cascadia", "numpy")
    rest = dep[len("numpy") :]
    assert rest.startswith("<2") or re.match(r"^=1\.", rest), dep


def test_cascadia_sqlalchemy_2() -> None:
    (dep,) = _conda("cascadia", "sqlalchemy")
    rest = dep[len("sqlalchemy") :]
    assert rest.startswith(">=2") or re.match(r"^=2\.", rest), dep


@pytest.mark.parametrize("pkg", ["lxml", "rdkit", "numba", "h5py"])
def test_cascadia_other_dependencies(pkg: str) -> None:
    assert _conda("cascadia", pkg), pkg


def test_cascadia_cuda_11_8_build() -> None:
    conda, _ = _deps("cascadia")
    cuda = [d for d in conda if "cuda" in d.lower()]
    assert any("11.8" in d or "cuda118" in d.lower() for d in cuda), conda


def test_cascadia_pip_pin() -> None:
    _, pip = _deps("cascadia")
    assert "cascadia==0.0.7" in pip, pip


def test_cascadia_header_records_source_licence_sha_and_date() -> None:
    h = _header("cascadia")
    assert "PyPI" in h
    assert "github.com/Noble-Lab/cascadia" in h
    assert "Apache-2.0" in h
    shas = re.findall(r"\b[0-9a-f]{64}\b", h)
    assert any(s.startswith(CASCADIA_WHEEL_SHA_HEAD) and s.endswith(CASCADIA_WHEEL_SHA_TAIL)
               for s in shas), "the 0.0.7 wheel sha256"  # fmt: skip
    assert re.search(r"\b20\d\d-\d\d-\d\d\b", h), "a date"
    assert "4.7.5" in h, "Q10: the pyteomics pin is recorded in the header"


# ------------------------------------------------------------------ deliverable 2: timsconvert.yml


def test_timsconvert_env() -> None:
    data = _yml("timsconvert")
    assert data["name"] == "giap-timsconvert"
    assert data["channels"] == ["conda-forge"]
    (py,) = _conda("timsconvert", "python")
    assert re.match(r"^python=3\.11(\.|=|\*|$)", py), py
    assert _conda("timsconvert", "pip")
    _, pip = _deps("timsconvert")
    want = f"timsconvert @ git+https://github.com/gtluu/timsconvert@{TIMSCONVERT_COMMIT}"
    assert [p.replace(" ", "") for p in pip] == [want.replace(" ", "")], pip


def test_timsconvert_header() -> None:
    h = _header("timsconvert")
    assert "B2" in h
    assert "TDF-SDK" in h
    assert re.search(r"ADR-(M10-B|\d{4})", h)
    assert re.search(r"never redistribut", h, re.I)


@pytest.mark.parametrize("stem", PILOT)
def test_no_defaults_channel(stem: str) -> None:
    text = require(ENVS / f"{stem}.yml", "deliverables 1-2").read_text()
    assert "defaults" not in text.split("dependencies:")[0]
    assert "nodefaults" not in text  # channels are exactly conda-forge


@pytest.mark.parametrize("stem", PILOT)
def test_locks_when_present(stem: str) -> None:
    pip_lock = ENVS / "locks" / f"{stem}.pip.txt"
    if not pip_lock.is_file():
        pytest.skip(f"{stem} locks come from the FarmShare env build (acceptance 14)")
    text = pip_lock.read_text()
    if stem == "cascadia":
        assert re.search(r"^cascadia==0\.0\.7$", text, re.M)
        assert re.search(r"^pyteomics==4\.7\.5$", text, re.M)
    else:
        assert TIMSCONVERT_COMMIT in text
    assert (ENVS / "locks" / f"{stem}.explicit.txt").is_file()


# ------------------------------------------------------------------ deliverable 3: PILOT_TOOLS


def test_pilot_tools_entries() -> None:
    from giap.annotate.tools import PILOT_TOOLS, TOOLS, ToolSpec

    assert set(PILOT_TOOLS) == set(PILOT)
    assert not set(PILOT_TOOLS) & set(TOOLS), "Q4: never in TOOLS"
    for name, spec in PILOT_TOOLS.items():
        assert isinstance(spec, ToolSpec)
        assert spec.name == name
        assert spec.env_name == f"giap-{name}"
        assert spec.env_file == f"workflow/envs/{name}.yml"
        assert (REPO / spec.env_file).is_file()
        assert spec.licence == "Apache-2.0"
        assert not spec.in_process
        assert spec.version_cmd
    assert PILOT_TOOLS["cascadia"].version == CASCADIA_VERSION
    assert PILOT_TOOLS["timsconvert"].version in {"2.0.0", f"git:{TIMSCONVERT_COMMIT}"}


def test_pilot_tools_is_read_only() -> None:
    from giap.annotate.tools import PILOT_TOOLS

    with pytest.raises(TypeError):
        PILOT_TOOLS["x"] = PILOT_TOOLS["cascadia"]  # type: ignore[index]


@pytest.mark.parametrize(
    "overrides",
    [
        [],
        ["annotation.terminators.bactermfinder=true",
         "annotation.terminators.cross_check_transtermhp=true"],
    ],
)  # fmt: skip
def test_required_tools_never_returns_pilot_tools(overrides: list[str]) -> None:
    from giap.annotate.tools import required_tools
    from giap.config import resolve_config

    cfg = resolve_config(overrides=overrides)
    for abundance in (False, True):
        names = {t.name for t in required_tools(cfg, abundance=abundance)}
        assert not names & set(PILOT), names


# ------------------------------------------------------------------ build_envs.sh lookup

FAKE_MICROMAMBA = r"""#!/bin/bash
# Fake micromamba (bash 3.2): logs argv, creates env prefixes, prints $FAKE_DIR/version_<env>.
echo "micromamba $*" >> "$FAKE_DIR/calls.log"
envs="$MAMBA_ROOT_PREFIX/envs"
name_of() { while [ "$#" -gt 0 ]; do if [ "$1" = "-n" ]; then echo "$2"; return; fi; shift; done; }
case "$1" in
  env)
    case "$2" in
      export)
        echo "@EXPLICIT"
        echo "https://conda.anaconda.org/conda-forge/linux-64/fake-1.0-0.conda"
        exit 0 ;;
      remove) n="$(name_of "$@")"; rm -rf "${envs:?}/$n"; exit 0 ;;
    esac
    exit 0 ;;
  create)
    n="$(name_of "$@")"; mkdir -p "$envs/$n/conda-meta" "$envs/$n/bin"; exit 0 ;;
  run)
    shift; n=""
    while [ "$#" -gt 0 ]; do
      case "$1" in -n) n="$2"; shift 2 ;; -*) shift ;; *) break ;; esac
    done
    [ -d "$envs/$n" ] || { echo "no env $n" >&2; exit 1; }
    case " $* " in *" pip "*freeze*) echo "fakepkg==1.0"; exit 0 ;; esac
    if [ -e "$FAKE_DIR/version_$n" ]; then cat "$FAKE_DIR/version_$n"; exit 0; fi
    echo "no version for $n" >&2; exit 1 ;;
esac
exit 0
"""

FAKE_GIT = r"""#!/bin/bash
echo "git $*" >> "$FAKE_DIR/calls.log"
if [ "$1" = "clone" ]; then for a in "$@"; do last="$a"; done; mkdir -p "$last"; fi
exit 0
"""


class BuildBox:
    def __init__(self, tmp: Path) -> None:
        self.repo, self.root, self.fake = tmp / "repo", tmp / "root", tmp / "fake"
        for d in (self.repo / "workflow" / "envs" / "locks", self.repo / "workflow" / "conf",
                  self.fake, self.root / "bin"):  # fmt: skip
            d.mkdir(parents=True)
        (self.repo / "src").symlink_to(REPO / "src")
        (self.repo / "scripts").symlink_to(REPO / "scripts")
        for yml in ENVS.glob("*.yml"):
            shutil.copy(yml, self.repo / "workflow" / "envs" / yml.name)
        shutil.copy(REPO / "workflow" / "conf" / "containers.config",
                    self.repo / "workflow" / "conf" / "containers.config")  # fmt: skip
        for name, body in (("micromamba", FAKE_MICROMAMBA), ("git", FAKE_GIT)):
            p = self.root / "bin" / name
            p.write_text(body)
            p.chmod(0o755)
        (self.fake / "calls.log").write_text("")

    def version(self, tool: str, text: str) -> None:
        (self.fake / f"version_giap-{tool}").write_text(text + "\n")

    def run(self, *tools: str) -> subprocess.CompletedProcess[str]:
        env = {k: v for k, v in os.environ.items() if k in {"PATH", "LANG", "TMPDIR"}}
        env.update(GIAP_REPO=str(self.repo), GIAP_ROOT=str(self.root), GIAP_PYTHON=sys.executable,
                   FAKE_DIR=str(self.fake), HOME=str(self.root / "home"))  # fmt: skip
        return subprocess.run(["bash", str(self.repo / "scripts" / "build_envs.sh"), *tools],
                              capture_output=True, text=True, cwd=self.repo, env=env,
                              timeout=300, check=False)  # fmt: skip

    def calls(self) -> list[str]:
        return (self.fake / "calls.log").read_text().splitlines()


@pytest.mark.parametrize("tool", PILOT)
def test_build_envs_resolves_pilot_tool(tmp_path: Path, tool: str) -> None:
    require(ENVS / f"{tool}.yml", "deliverables 1-2")
    from giap.annotate.tools import PILOT_TOOLS

    box = BuildBox(tmp_path)
    box.version(tool, PILOT_TOOLS[tool].version.removeprefix("git:"))
    res = box.run(tool)
    assert res.returncode == 0, res.stdout + res.stderr
    creates = [c for c in box.calls() if c.startswith("micromamba create")]
    assert len(creates) == 1, box.calls()
    assert f"-n giap-{tool}" in creates[0]
    assert f"workflow/envs/{tool}.yml" in creates[0]
    marker = box.root / "micromamba" / "envs" / f"giap-{tool}" / ".giap-built"
    assert marker.is_file()


def test_build_envs_all_lists_no_pilot_tool(tmp_path: Path) -> None:
    require(ENVS / "cascadia.yml", "deliverable 1")
    require(ENVS / "timsconvert.yml", "deliverable 2")
    from giap.annotate.tools import TOOLS

    box = BuildBox(tmp_path)
    for name, spec in TOOLS.items():
        box.version(name, spec.version.removeprefix("git:"))
    box.version("fgsea", "[1] '1.36.2'")
    res = box.run("all")
    line = next((ln for ln in res.stdout.splitlines() if ln.startswith("all:")), None)
    assert line is not None, res.stdout + res.stderr
    assert not set(line.split()[1:]) & set(PILOT), line
    assert not [c for c in box.calls() if "giap-cascadia" in c or "giap-timsconvert" in c]


def test_build_envs_usage_lists_pilot_tools() -> None:
    require(ENVS / "cascadia.yml", "deliverable 1")
    env = {k: v for k, v in os.environ.items() if k in {"PATH", "HOME", "LANG"}}
    env["GIAP_REPO"] = str(REPO)
    cmd = ["bash", str(REPO / "scripts" / "build_envs.sh"), "no_such_tool"]
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=60, check=False)
    assert res.returncode == 2
    assert "cascadia" in res.stderr and "timsconvert" in res.stderr
