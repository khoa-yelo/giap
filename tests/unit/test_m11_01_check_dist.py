"""``scripts/release/check_dist.py`` on fake artefacts (issue M11-01 item 4, examples 4-7).

Every wheel and sdist here is built in ``tmp_path`` with ``zipfile``/``tarfile``: no network,
no hatchling, no install. ``_wheel``/``_sdist`` write a *valid* artefact (every required file of
item 4, including each ``workflow/envs/*.yml`` and ``workflow/conf/*.config`` of the checkout);
a test then adds, drops or edits entries to provoke one rule.
"""

from __future__ import annotations

import io
import re
import subprocess
import sys
import tarfile
import zipfile
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "check_dist.py"
_V = re.search(r'__version__\s*=\s*"([^"]+)"', (REPO / "src/giap/__init__.py").read_text())
assert _V is not None
V = _V.group(1)  # giap.__version__ of the checkout (the default --version)
MB = 1024 * 1024

ENVS = sorted(p.name for p in (REPO / "workflow/envs").glob("*.yml"))
CONFS = sorted(p.name for p in (REPO / "workflow/conf").glob("*.config"))

WHEEL_REQUIRED = [
    "giap/__init__.py",
    "giap/cli.py",
    "giap/py.typed",
    "giap/config/defaults.yaml",
    "giap/schemas/datapackage.json",
    "giap/data/calibration.yaml",
    "giap/data/db_sources.yaml",
    "giap/workflow/main.nf",
    "giap/workflow/nextflow.config",
    *(f"giap/workflow/envs/{n}" for n in ENVS),
    *(f"giap/workflow/conf/{n}" for n in CONFS),
]
# allowed look-alikes: modules/local/ is not top-level local/; the two uniop/ dirs are exempt
WHEEL_EXTRA = [
    "giap/workflow/modules/local/bakta/main.nf",
    "giap/workflow/modules/local/uniop/main.nf",
    "giap/workflow/containers/uniop/Dockerfile",
    "giap/workflow/bin/giap_fgsea.R",
    "giap/workflow/assets/ps/ps",
]
SDIST_REQUIRED = [
    "pyproject.toml",
    "README.md",
    "LICENSE",
    "NOTICE.md",
    "CHANGELOG.md",
    "CITATION.cff",
    "src/giap/__init__.py",
    "workflow/main.nf",
]
SDIST_EXTRA = [
    "src/giap/cli.py",
    "workflow/nextflow.config",
    "workflow/modules/local/uniop/main.nf",
    "workflow/containers/uniop/Dockerfile",
    "workflow/bin/giap_fgsea.R",
]


def _whl_name(v: str = V) -> str:
    return f"giap-{v}-py3-none-any.whl"


def _sdist_name(v: str = V) -> str:
    return f"giap-{v}.tar.gz"


def _metadata(v: str, license_line: str) -> bytes:
    return f"Metadata-Version: 2.4\nName: giap\nVersion: {v}\n{license_line}\n".encode()


def _wheel(
    dist: Path,
    *,
    name_v: str = V,
    meta_v: str = V,
    add: dict[str, bytes] | None = None,
    drop: frozenset[str] = frozenset(),
    symlinks: dict[str, str] | None = None,
    license_line: str = "License-Expression: MIT",
    entry_points: str = "[console_scripts]\ngiap = giap.cli:main\n",
    license_file: bool = True,
) -> Path:
    path = dist / _whl_name(name_v)
    di = f"giap-{meta_v}.dist-info"
    files: dict[str, bytes] = {p: b"x\n" for p in WHEEL_REQUIRED + WHEEL_EXTRA}
    files[f"{di}/METADATA"] = _metadata(meta_v, license_line)
    files[f"{di}/WHEEL"] = (
        b"Wheel-Version: 1.0\nGenerator: fake\nRoot-Is-Purelib: true\nTag: py3-none-any\n"
    )
    files[f"{di}/entry_points.txt"] = entry_points.encode()
    if license_file:
        files[f"{di}/licenses/LICENSE"] = b"MIT License\n"
    files.update(add or {})
    with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as zf:
        for arc, data in files.items():
            if arc in drop:
                continue
            info = zipfile.ZipInfo(arc, date_time=(2026, 10, 24, 0, 0, 0))
            mode = 0o755 if "/bin/" in arc or arc.endswith("/ps/ps") else 0o644
            info.external_attr = (0o100000 | mode) << 16
            zf.writestr(info, data)
        for arc, target in (symlinks or {}).items():
            info = zipfile.ZipInfo(arc, date_time=(2026, 10, 24, 0, 0, 0))
            info.external_attr = (0o120000 | 0o777) << 16
            zf.writestr(info, target)
        zf.writestr(f"{di}/RECORD", "")
    return path


def _sdist(
    dist: Path,
    *,
    name_v: str = V,
    meta_v: str = V,
    add: dict[str, bytes] | None = None,
    drop: frozenset[str] = frozenset(),
    symlinks: dict[str, str] | None = None,
) -> Path:
    """``add`` keys are full archive paths (``giap-V/...``); ``drop`` keys are root-relative."""
    path = dist / _sdist_name(name_v)
    root = f"giap-{name_v}"
    files: dict[str, bytes] = {
        f"{root}/{p}": b"x\n" for p in SDIST_REQUIRED + SDIST_EXTRA if p not in drop
    }
    files[f"{root}/PKG-INFO"] = _metadata(meta_v, "License-Expression: MIT")
    files.update(add or {})
    with tarfile.open(path, "w:gz") as tf:
        for arc, data in files.items():
            info = tarfile.TarInfo(arc)
            info.size = len(data)
            info.mode = 0o755 if "/bin/" in arc else 0o644
            tf.addfile(info, io.BytesIO(data))
        for arc, target in (symlinks or {}).items():
            info = tarfile.TarInfo(arc)
            info.type = tarfile.SYMTYPE
            info.linkname = target
            tf.addfile(info)
    return path


def _run(dist: Path, *args: str) -> tuple[int, list[str]]:
    assert SCRIPT.is_file(), f"{SCRIPT.relative_to(REPO)} does not exist (M11-01 item 4)"
    p = subprocess.run(
        [sys.executable, str(SCRIPT), str(dist), *args],
        cwd=REPO,
        capture_output=True,
        text=True,
        timeout=120,
    )
    return p.returncode, (p.stdout + p.stderr).splitlines()


def _fails(lines: list[str]) -> list[str]:
    return [ln for ln in lines if ln.startswith("FAIL ")]


@pytest.fixture
def dist(tmp_path: Path) -> Path:
    d = tmp_path / "dist"
    d.mkdir()
    return d


# --- valid artefacts --------------------------------------------------------------------------


def test_valid_pair_exits_0_with_one_ok_line_per_artefact(dist: Path) -> None:
    _wheel(dist)
    _sdist(dist)
    (dist / ".gitignore").write_text("*\n")  # uv build writes this into dist/
    code, lines = _run(dist)
    assert code == 0, lines
    ok = [ln for ln in lines if ln.startswith("ok")]
    assert len(ok) == 2, lines
    assert any(_whl_name() in ln for ln in ok)
    assert any(_sdist_name() in ln for ln in ok)
    assert not _fails(lines)


def test_explicit_matching_version_exits_0(dist: Path) -> None:
    _wheel(dist)
    _sdist(dist)
    code, lines = _run(dist, "--version", V)
    assert code == 0, lines


def test_license_header_alternative_accepted(dist: Path) -> None:
    _wheel(dist, license_line="License: MIT")
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 0, lines


# --- examples 4-7 -----------------------------------------------------------------------------


def test_example4_file_over_5mb(dist: Path) -> None:
    _wheel(dist, add={"giap/data/huge.bin": b"\0" * (6 * MB)})
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    assert f"FAIL {_whl_name()}: file > 5 MB: giap/data/huge.bin" in lines
    assert len(_fails(lines)) == 1, lines


def test_file_under_5mb_is_fine(dist: Path) -> None:
    _wheel(dist, add={"giap/data/big.bin": b"\0" * (4 * 1000 * 1000)})
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 0, lines


def test_example5_sdist_forbidden_tests_path(dist: Path) -> None:
    _wheel(dist)
    _sdist(dist, add={f"giap-{V}/tests/data/x.tsv.gz": b"x"})
    code, lines = _run(dist)
    assert code == 1, lines
    assert f"FAIL {_sdist_name()}: forbidden path: giap-{V}/tests/data/x.tsv.gz" in lines
    assert len(_fails(lines)) == 1, lines


def test_example6_version_mismatch_one_fail_per_artefact(dist: Path) -> None:
    other = "1.0.0" if V != "1.0.0" else "9.9.9"
    _wheel(dist)
    _sdist(dist)
    code, lines = _run(dist, "--version", other)
    assert code == 1, lines
    fails = _fails(lines)
    assert len([f for f in fails if f.startswith(f"FAIL {_whl_name()}: version")]) == 1, lines
    assert len([f for f in fails if f.startswith(f"FAIL {_sdist_name()}: version")]) == 1, lines
    assert len(fails) == 2, lines


def test_example7_empty_dir_exits_2(dist: Path) -> None:
    code, lines = _run(dist)
    assert code == 2, lines
    assert any("no wheel and no sdist in" in ln and dist.name in ln for ln in lines), lines


@pytest.mark.parametrize("which", ["wheel", "sdist"])
def test_missing_one_kind_exits_2(dist: Path, which: str) -> None:
    (_wheel if which == "wheel" else _sdist)(dist)
    code, lines = _run(dist)
    assert code == 2, lines


# --- version rule -----------------------------------------------------------------------------


def test_default_version_is_checkout_version(dist: Path) -> None:
    _wheel(dist, name_v="7.7.7", meta_v="7.7.7")
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    assert [f for f in _fails(lines) if f.startswith("FAIL giap-7.7.7-py3-none-any.whl: version")]
    assert not [f for f in _fails(lines) if f.startswith(f"FAIL {_sdist_name()}")], lines


def test_metadata_version_disagreeing_with_name_fails(dist: Path) -> None:
    _wheel(dist, meta_v="7.7.7")
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    assert [f for f in _fails(lines) if f.startswith(f"FAIL {_whl_name()}: version")], lines


# --- required content -------------------------------------------------------------------------


@pytest.mark.parametrize(
    "missing",
    [
        "giap/py.typed",
        "giap/cli.py",
        "giap/config/defaults.yaml",
        "giap/data/db_sources.yaml",
        "giap/workflow/main.nf",
        "giap/workflow/nextflow.config",
        "giap/workflow/envs/bakta.yml",
        "giap/workflow/conf/test.config",
    ],
)
def test_wheel_missing_required_file(dist: Path, missing: str) -> None:
    _wheel(dist, drop=frozenset({missing}))
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    fails = _fails(lines)
    assert len(fails) == 1, lines
    assert fails[0].startswith(f"FAIL {_whl_name()}: ") and fails[0].endswith(missing), lines


def test_every_checkout_env_and_conf_is_required() -> None:
    assert "bakta.yml" in ENVS and "test.config" in CONFS  # the parametrisation above stays real


def test_wheel_without_console_script_fails(dist: Path) -> None:
    _wheel(dist, entry_points="[console_scripts]\nother = giap.cli:main\n")
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    assert [f for f in _fails(lines) if f.startswith(f"FAIL {_whl_name()}: ")], lines


def test_wheel_non_mit_license_fails(dist: Path) -> None:
    _wheel(dist, license_line="License-Expression: GPL-3.0-only")
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    assert [f for f in _fails(lines) if f.startswith(f"FAIL {_whl_name()}: ")], lines


def test_wheel_without_license_file_fails(dist: Path) -> None:
    _wheel(dist, license_file=False)
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    fails = _fails(lines)
    assert [f for f in fails if f.startswith(f"FAIL {_whl_name()}: ") and "LICENSE" in f], lines


@pytest.mark.parametrize(
    "missing",
    ["pyproject.toml", "README.md", "LICENSE", "NOTICE.md", "CHANGELOG.md", "CITATION.cff"],
)
def test_sdist_missing_required_file(dist: Path, missing: str) -> None:
    _wheel(dist)
    _sdist(dist, drop=frozenset({missing}))
    code, lines = _run(dist)
    assert code == 1, lines
    fails = _fails(lines)
    assert len(fails) == 1, lines
    assert fails[0].startswith(f"FAIL {_sdist_name()}: ") and fails[0].endswith(missing), lines


@pytest.mark.parametrize("tree", ["src", "workflow"])
def test_sdist_without_src_or_workflow_fails(dist: Path, tree: str) -> None:
    _wheel(dist)
    _sdist(dist, drop=frozenset(p for p in SDIST_REQUIRED + SDIST_EXTRA if p.startswith(tree)))
    code, lines = _run(dist)
    assert code == 1, lines
    assert [f for f in _fails(lines) if f.startswith(f"FAIL {_sdist_name()}: ")], lines


# --- forbidden content ------------------------------------------------------------------------

R = f"giap-{V}"
FORBIDDEN_WHEEL = [
    "local/notes.txt",
    "tests/test_x.py",
    "docs/index.md",
    "giap/workflow/tests/main.nf.test",
    "giap/data/model.sif",
    "giap/data/run.raw",
    "giap/data/run.mzML",
    "giap/data/table.parquet",
    "giap/cli.pyc",
    "giap/__pycache__/cli.cpython-312.pyc",
    "giap/data/RegulonDB_network.tsv",
    "giap/data/x_regulondb.txt",
    "giap/data/uniop/op.py",
    "giap/workflow/uniop/run.py",
]
FORBIDDEN_SDIST = [
    f"{R}/local/notes.txt",
    f"{R}/tests/data/x.tsv.gz",
    f"{R}/validation/a.tsv",
    f"{R}/.github/workflows/ci.yml",
    f"{R}/.claude/agents/x.md",
    f"{R}/docs/index.md",
    f"{R}/benchmarks/baseline.json",
    f"{R}/workflow/tests/main.nf.test",
    f"{R}/workflow/assets/run.d",
    f"{R}/src/giap/data/regulondb.tsv",
    f"{R}/src/giap/uniop/op.py",
]


@pytest.mark.parametrize("path", FORBIDDEN_WHEEL)
def test_wheel_forbidden(dist: Path, path: str) -> None:
    _wheel(dist, add={path: b"x"})
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    fails = _fails(lines)
    assert fails and all(f.startswith(f"FAIL {_whl_name()}: ") for f in fails), lines
    assert {f.rsplit(": ", 1)[-1] for f in fails} == {path}, lines


@pytest.mark.parametrize("path", FORBIDDEN_SDIST)
def test_sdist_forbidden(dist: Path, path: str) -> None:
    _wheel(dist)
    _sdist(dist, add={path: b"x"})
    code, lines = _run(dist)
    assert code == 1, lines
    fails = _fails(lines)
    assert fails and all(f.startswith(f"FAIL {_sdist_name()}: ") for f in fails), lines
    assert {f.rsplit(": ", 1)[-1] for f in fails} == {path}, lines


def test_symlink_in_wheel_rejected(dist: Path) -> None:
    _wheel(dist, symlinks={"giap/workflow/link.nf": "main.nf"})
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    assert f"FAIL {_whl_name()}: symlink: giap/workflow/link.nf" in lines


def test_symlink_in_sdist_rejected(dist: Path) -> None:
    _wheel(dist)
    _sdist(dist, symlinks={f"{R}/workflow/link.nf": "main.nf"})
    code, lines = _run(dist)
    assert code == 1, lines
    assert f"FAIL {_sdist_name()}: symlink: {R}/workflow/link.nf" in lines


def test_all_violations_reported_not_the_first(dist: Path) -> None:
    _wheel(dist, add={"giap/data/huge.bin": b"\0" * (6 * MB), "giap/data/model.sif": b"x"})
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    fails = _fails(lines)
    assert len(fails) == 2, lines
    assert {f.rsplit(": ", 1)[-1] for f in fails} == {"giap/data/huge.bin", "giap/data/model.sif"}


# --- property: the flagged set is exactly the injected set ------------------------------------

_INJECT = [("wheel", p) for p in FORBIDDEN_WHEEL[:7]] + [("sdist", p) for p in FORBIDDEN_SDIST[:7]]


@settings(
    max_examples=12, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(chosen=st.sets(st.sampled_from(_INJECT), max_size=6))
def test_property_every_violation_reported_once_per_path(
    tmp_path_factory: pytest.TempPathFactory, chosen: set[tuple[str, str]]
) -> None:
    d = tmp_path_factory.mktemp("dist")
    _wheel(d, add={p: b"x" for k, p in chosen if k == "wheel"})
    _sdist(d, add={p: b"x" for k, p in chosen if k == "sdist"})
    code, lines = _run(d)
    assert code == (1 if chosen else 0), lines
    names = {"wheel": _whl_name(), "sdist": _sdist_name()}
    flagged = set()
    for f in _fails(lines):
        art, _, rest = f[len("FAIL ") :].partition(": ")
        flagged.add((art, rest.rsplit(": ", 1)[-1]))
    assert flagged == {(names[k], p) for k, p in chosen}, lines


# --- review #330 S2 (lead): boundaries the hand mutants survived ------------------------------


@pytest.mark.parametrize(("size", "ok"), [(5 * MB, True), (5 * MB + 1, False), (5_000_001, True)])
def test_5mb_limit_is_strictly_greater_than_5_mib(dist: Path, size: int, ok: bool) -> None:
    """The R11 limit is 5 MiB (``MB`` here is 2**20): exactly 5 MiB passes (``>``, not ``>=``), and
    5,000,001 bytes (over a decimal 5 MB) passes too."""
    _wheel(dist, add={"giap/data/edge.bin": b"\0" * size})
    _sdist(dist)
    code, lines = _run(dist)
    assert (code == 0) is ok, lines


def test_hardlink_in_sdist_rejected(dist: Path) -> None:
    sdist = _sdist(dist)
    with tarfile.open(sdist, "r:gz") as tf:
        members = [(m, tf.extractfile(m).read() if m.isfile() else None)  # type: ignore[union-attr]
                   for m in tf.getmembers()]  # fmt: skip
    link = tarfile.TarInfo(f"giap-{V}/src/giap/hard.py")
    link.type, link.linkname = tarfile.LNKTYPE, f"giap-{V}/README.md"
    with tarfile.open(sdist, "w:gz") as tf:
        for m, data in members:
            tf.addfile(m, io.BytesIO(data) if data is not None else None)
        tf.addfile(link)
    _wheel(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    assert f"FAIL {_sdist_name()}: symlink: giap-{V}/src/giap/hard.py" in lines, lines


def test_metadata_with_two_version_headers_fails(dist: Path) -> None:
    di = f"giap-{V}.dist-info"
    meta = _metadata(V, "License-Expression: MIT") + f"Version: {V}\n".encode()
    _wheel(dist, add={f"{di}/METADATA": meta})
    _sdist(dist)
    code, lines = _run(dist)
    assert code == 1, lines
    assert [ln for ln in _fails(lines) if ": version: " in ln], lines
