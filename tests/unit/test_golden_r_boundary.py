"""The R licence boundary and the golden-provenance guard (M8-12 "Licences and the boundary",
"Regeneration"; M8-index §Licences; GIAP_SPEC.md R5, §11.5, §14.6).

Interface assumed (the spec leaves the names open):

- ``scripts/ci/check_imports.py`` has ``rpy2`` in ``FORBIDDEN`` (spec example 5).
- ``scripts/ci/check_no_rscript.py [ROOT]`` (default ``src``): exit 1 when any file under ROOT
  (every file, not only ``*.py``) contains the text ``Rscript``; prints one
  ``<path>:<line>: ...`` line per hit. Importable: ``violations(root: Path) -> list[str]``.
- ``scripts/ci/check_golden_provenance.py [BASE]`` (default ``origin/main``): exit 1 when the
  diff ``BASE...HEAD`` touches ``tests/data/golden/r/`` without touching
  ``tests/data/golden/r/PROVENANCE.md``. Importable: ``violations(changed: Iterable[str]) ->
  list[str]`` over repo-relative POSIX paths.
- ``.github/workflows/ci.yml`` runs both new checks.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from collections.abc import Iterable
from pathlib import Path
from types import ModuleType

import pytest
from hypothesis import assume, given
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
CI = REPO / "scripts" / "ci"
CHECK_IMPORTS = CI / "check_imports.py"
CHECK_RSCRIPT = CI / "check_no_rscript.py"
CHECK_PROV = CI / "check_golden_provenance.py"
CI_YML = REPO / ".github" / "workflows" / "ci.yml"
GOLDEN_R = "tests/data/golden/r"
PROVENANCE = f"{GOLDEN_R}/PROVENANCE.md"


def _load(path: Path) -> ModuleType:
    assert path.is_file(), f"missing {path.relative_to(REPO)}"
    spec = importlib.util.spec_from_file_location(path.stem, path)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run(script: Path, *args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    assert script.is_file(), f"missing {script.relative_to(REPO)}"
    return subprocess.run(
        [sys.executable, str(script), *args], capture_output=True, text=True, cwd=cwd
    )


# --- rpy2 in check_imports.py (spec example 5) ---------------------------------------------


def test_rpy2_is_forbidden() -> None:
    assert "rpy2" in _load(CHECK_IMPORTS).FORBIDDEN


@pytest.mark.parametrize(
    "line",
    [
        "import rpy2\n",
        "import rpy2.robjects as ro\n",
        "from rpy2.robjects import r\n",
        "from rpy2 import robjects\n",
    ],
)
def test_check_imports_exits_1_on_rpy2(tmp_path: Path, line: str) -> None:
    src = tmp_path / "src" / "giap"
    src.mkdir(parents=True)
    (src / "bridge.py").write_text(line)
    proc = _run(CHECK_IMPORTS, str(tmp_path / "src"))
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert "rpy2" in proc.stdout


def test_check_imports_still_passes_on_the_real_src() -> None:
    proc = _run(CHECK_IMPORTS, str(REPO / "src"))
    assert proc.returncode == 0, proc.stdout


# --- the Rscript guard ---------------------------------------------------------------------


@pytest.mark.parametrize(
    ("name", "text"),
    [
        ("run.py", 'subprocess.run(["Rscript", "x.R"])\n'),
        ("run.py", "cmd = 'Rscript --vanilla golden.R'\n"),
        ("main.nf", "script:\n  Rscript make.R\n"),
        ("helper.sh", "#!/bin/bash\nRscript -e 'library(fgsea)'\n"),
        ("data.yaml", "tool: Rscript\n"),
    ],
)
def test_rscript_guard_exits_1_on_any_file_mentioning_rscript(
    tmp_path: Path, name: str, text: str
) -> None:
    root = tmp_path / "src"
    (root / "giap" / "sub").mkdir(parents=True)
    (root / "giap" / "sub" / name).write_text(text)
    (root / "giap" / "clean.py").write_text("x = 1\n")
    proc = _run(CHECK_RSCRIPT, str(root))
    assert proc.returncode == 1, proc.stdout + proc.stderr
    assert name in proc.stdout
    assert "clean.py" not in proc.stdout


def test_rscript_guard_reports_the_line_number(tmp_path: Path) -> None:
    root = tmp_path / "src"
    root.mkdir()
    (root / "a.py").write_text("a = 1\nb = 2\nc = 'Rscript'\n")
    found = _load(CHECK_RSCRIPT).violations(root)
    assert len(found) == 1
    assert found[0].startswith(f"{root / 'a.py'}:3:")


def test_rscript_guard_passes_a_clean_tree(tmp_path: Path) -> None:
    root = tmp_path / "src"
    root.mkdir()
    (root / "a.py").write_text("# R is not used here; fgsea is ported to Python\nx = 1\n")
    (root / "blob.bin").write_bytes(bytes(range(256)))  # binary files must not crash the guard
    proc = _run(CHECK_RSCRIPT, str(root))
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert proc.stdout.strip() == ""


def test_rscript_guard_passes_the_real_src() -> None:
    proc = _run(CHECK_RSCRIPT, str(REPO / "src"))
    assert proc.returncode == 0, proc.stdout


def test_rscript_guard_defaults_to_src() -> None:
    proc = _run(CHECK_RSCRIPT, cwd=REPO)
    assert proc.returncode == 0, proc.stdout + proc.stderr


_SAFE_TEXT = st.text(
    alphabet=st.characters(blacklist_categories=("Cs",), blacklist_characters="\x00"),
    max_size=200,
).filter(lambda s: "Rscript" not in s)


@given(prefix=_SAFE_TEXT, suffix=_SAFE_TEXT, planted=st.booleans())
def test_rscript_guard_flags_exactly_the_files_that_mention_it(
    tmp_path_factory: pytest.TempPathFactory, prefix: str, suffix: str, planted: bool
) -> None:
    assume(planted or "Rscript" not in prefix + suffix)  # PR #192 B3: the join can spell it
    root = tmp_path_factory.mktemp("src")
    body = prefix + ("Rscript" if planted else "") + suffix
    (root / "f.txt").write_text(body, encoding="utf-8")
    found = _load(CHECK_RSCRIPT).violations(root)
    assert bool(found) == planted


# --- the provenance guard ------------------------------------------------------------------


def _prov_violations(changed: Iterable[str]) -> list[str]:
    result: list[str] = _load(CHECK_PROV).violations(list(changed))
    return result


@pytest.mark.parametrize(
    ("changed", "fails"),
    [
        ([f"{GOLDEN_R}/outputs/fgsea.tsv.gz"], True),
        ([f"{GOLDEN_R}/VERSIONS.txt"], True),
        ([f"{GOLDEN_R}/inputs/gsea_sets.tsv.gz", "src/giap/x.py"], True),
        ([f"{GOLDEN_R}/outputs/fgsea.tsv.gz", PROVENANCE], False),
        ([PROVENANCE], False),
        (["src/giap/stats/fgsea.py", "tests/data/golden/sim50k/mapping/a.tsv"], False),
        ([], False),
        (["tests/data/golden/rx/a.tsv"], False),  # a sibling directory, not golden/r/
    ],
)
def test_provenance_guard_examples(changed: list[str], fails: bool) -> None:
    assert bool(_prov_violations(changed)) == fails


_PATHS = st.sampled_from(
    [
        f"{GOLDEN_R}/outputs/fgsea.tsv.gz",
        f"{GOLDEN_R}/outputs/goseq.tsv.gz",
        f"{GOLDEN_R}/inputs/ora_genes.tsv.gz",
        f"{GOLDEN_R}/VERSIONS.txt",
        PROVENANCE,
        "src/giap/stats/wallenius.py",
        "docs/specs/M8-12-golden-r-env.md",
        "tests/data/golden/sim50k/m4/x.tsv",
    ]
)


@given(changed=st.lists(_PATHS, unique=True))
def test_provenance_guard_property(changed: list[str]) -> None:
    touches = any(p.startswith(f"{GOLDEN_R}/") and p != PROVENANCE for p in changed)
    expected = touches and PROVENANCE not in changed
    assert bool(_prov_violations(changed)) == expected


def _git(repo: Path, *args: str) -> None:
    subprocess.run(
        ["git", "-c", "user.name=t", "-c", "user.email=t@example.org", *args],
        cwd=repo,
        check=True,
        capture_output=True,
    )


@pytest.mark.parametrize("update_provenance", [False, True])
def test_provenance_guard_cli_on_a_git_diff(tmp_path: Path, update_provenance: bool) -> None:
    repo = tmp_path / "repo"
    out = repo / GOLDEN_R / "outputs"
    out.mkdir(parents=True)
    (out / "fgsea.tsv.gz").write_bytes(b"old")
    (repo / PROVENANCE).write_text("job 1\n")
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "base")
    _git(repo, "checkout", "-q", "-b", "feature")
    (out / "fgsea.tsv.gz").write_bytes(b"new")
    if update_provenance:
        (repo / PROVENANCE).write_text("job 2\n")
    _git(repo, "commit", "-q", "-am", "regenerate")
    proc = _run(CHECK_PROV, "main", cwd=repo)
    assert proc.returncode == (0 if update_provenance else 1), proc.stdout + proc.stderr


# --- wired into CI -------------------------------------------------------------------------


@pytest.mark.parametrize("script", ["check_no_rscript.py", "check_golden_provenance.py"])
def test_ci_yml_runs_the_new_checks(script: str) -> None:
    assert f"scripts/ci/{script}" in CI_YML.read_text()
