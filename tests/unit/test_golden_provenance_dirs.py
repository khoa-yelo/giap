"""The provenance guard covers both golden directories (docs/specs/M8-12b-golden-median-exact.md
§4.3, B9, example 10).

``scripts/ci/check_golden_provenance.py``: ``GOLDEN_R`` becomes the tuple
``("tests/data/golden/r/", "tests/data/golden/wallenius/")``; a changed file under either
directory needs **that** directory's ``PROVENANCE.md`` in the same diff. ``violations(changed)``
keeps its signature; the M8-12 cases in tests/unit/test_golden_r_boundary.py still hold.
"""

from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path
from types import ModuleType

import pytest
from hypothesis import given
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
CHECK = REPO / "scripts" / "ci" / "check_golden_provenance.py"
CI_YML = REPO / ".github" / "workflows" / "ci.yml"
R = "tests/data/golden/r/"
W = "tests/data/golden/wallenius/"
R_PROV, W_PROV = R + "PROVENANCE.md", W + "PROVENANCE.md"
W_FILE = W + "exact_tails.tsv.gz"


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("check_golden_provenance_dirs", CHECK)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _violations(changed: list[str]) -> list[str]:
    result: list[str] = _load().violations(changed)
    return result


def test_golden_dirs_tuple() -> None:
    assert tuple(_load().GOLDEN_R) == (R, W)


@pytest.mark.parametrize(
    ("changed", "fails"),
    [
        ([W_FILE], True),
        ([W_FILE, W_PROV], False),
        ([W_FILE, R_PROV], True),  # the other directory's PROVENANCE does not count
        ([R + "outputs/fgsea.tsv.gz", W_PROV], True),
        ([R + "outputs/fgsea.tsv.gz", R_PROV, W_FILE], True),
        ([R + "outputs/fgsea.tsv.gz", R_PROV, W_FILE, W_PROV], False),
        ([W_PROV], False),
        (["tests/data/golden/wallenius_old/a.tsv"], False),  # a sibling, not golden/wallenius/
        (["tests/data/golden/rx/a.tsv"], False),
    ],
)
def test_example_10_cases(changed: list[str], fails: bool) -> None:
    assert bool(_violations(changed)) == fails


def test_violation_names_the_file() -> None:
    found = _violations([W_FILE, R_PROV])
    assert len(found) == 1
    assert W_FILE in found[0] and W_PROV in found[0]


_PATHS = st.sampled_from(
    [
        R + "outputs/goseq.tsv.gz",
        R + "inputs/gsea_stats.tsv.gz",
        R_PROV,
        W_FILE,
        W + "notes.txt",
        W_PROV,
        "src/giap/function/ora.py",
        "tests/data/golden/sim50k/m4/x.tsv",
    ]
)


@given(changed=st.lists(_PATHS, unique=True))
def test_each_directory_needs_its_own_provenance(changed: list[str]) -> None:
    expected = set()
    for root, prov in ((R, R_PROV), (W, W_PROV)):
        if prov not in changed:
            expected |= {p for p in changed if p.startswith(root) and p != prov}
    found = _violations(changed)
    assert len(found) == len(expected)
    for p in expected:
        assert any(p in line for line in found), p


def _git(repo: Path, *args: str) -> None:
    subprocess.run(
        ["git", "-c", "user.name=t", "-c", "user.email=t@example.org", *args],
        cwd=repo,
        check=True,
        capture_output=True,
    )


@pytest.mark.parametrize("update_provenance", [False, True])
def test_cli_on_a_wallenius_diff(tmp_path: Path, update_provenance: bool) -> None:
    repo = tmp_path / "repo"
    (repo / W).mkdir(parents=True)
    (repo / W_FILE).write_bytes(b"old")
    (repo / W_PROV).write_text("job 1\n")
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-q", "-m", "base")
    _git(repo, "checkout", "-q", "-b", "feature")
    (repo / W_FILE).write_bytes(b"new")
    if update_provenance:
        (repo / W_PROV).write_text("job 2\n")
    _git(repo, "commit", "-q", "-am", "regenerate")
    proc = subprocess.run(
        [sys.executable, str(CHECK), "main"], capture_output=True, text=True, cwd=repo
    )
    assert proc.returncode == (0 if update_provenance else 1), proc.stdout + proc.stderr
    if not update_provenance:
        assert W_FILE in proc.stdout


def test_ci_step_name_mentions_both_issues() -> None:
    assert "golden data carry a PROVENANCE.md update (M8-12, M8-12b)" in CI_YML.read_text()
