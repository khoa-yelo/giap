"""Evidence reuse predicate ``reusable(C, R)`` and its two path lists (M11-14 14g; ADR-0060 §4).

Spec docs/specs/M11-14-container-install.md item 4 "Evidence reuse", examples 13b, 13c and 13d,
and AC 13 (version-line exception; ``NON_STATS_PATHS`` and the never list disjoint, every entry
exists). Every git test builds its own throw-away repository; nothing touches the network.

Interface fixed by these tests (in ``scripts/release/readiness.py``, which holds the lists):

``NON_STATS_PATHS: tuple[str, ...]``
    exactly the 14 entries of the spec; a directory entry ends with ``/``, a file entry does not.
``NEVER_PATHS: tuple[str, ...]``
    the never list, exactly the 15 entries of the spec (same notation).
``non_stats_path(path: str) -> bool``
    ``path`` (repo-relative POSIX) lies under an entry of ``NON_STATS_PATHS``, matched on whole
    path components: a directory entry ``a/b/`` matches ``a/b/x`` but not ``a/bx``; a file entry
    matches only itself. Pure, no git.
``reusable(repo: Path, c: str, r: str) -> bool``
    ``D = git diff --name-only --no-renames C R``; true iff every path of ``D`` is a
    ``non_stats_path``, except ``src/giap/__init__.py`` when ``git diff -U0 C R -- <it>`` is one
    hunk with one ``-`` and one ``+`` line, both matching ``^__version__ = "[^"]+"$``; always
    false when ``D`` contains ``scripts/release/readiness.py``. Tests pass full 40-hex SHAs.

Not a statistical method, so there is no calibration test (spec §8.4); the predicate decides
when calibration evidence is reused, hence the exact-list pins below (Q-L15: ``stats`` review).
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest
from helpers.release_readiness import REPO, load_module
from hypothesis import given, settings
from hypothesis import strategies as st

SPEC_NON_STATS = {
    "docs/", "changelog.d/", "README.md", "CHANGELOG.md", "mkdocs.yml", "PROGRESS.md",
    "AGENTS.md", "CLAUDE.md", ".claude/", ".github/", "scripts/dev/", "scripts/docs/",
    "src/giap/report/", "src/giap/viz/",
}  # fmt: skip
SPEC_NEVER = {
    "src/giap/preflight/", "src/giap/runner/", "src/giap/setupdb/", "src/giap/cli.py",
    "src/giap/data/", "src/giap/config/", "src/giap/rng.py", "src/giap/io/", "workflow/",
    "pyproject.toml", "uv.lock", "tests/sim/", "tests/calibration/", "ci/",
    "scripts/release/readiness.py",
}  # fmt: skip

_GIT_ENV = {
    "GIT_AUTHOR_NAME": "Fixture", "GIT_AUTHOR_EMAIL": "fixture@example.invalid",
    "GIT_COMMITTER_NAME": "Fixture", "GIT_COMMITTER_EMAIL": "fixture@example.invalid",
    "GIT_CONFIG_GLOBAL": os.devnull, "GIT_CONFIG_SYSTEM": os.devnull,
}  # fmt: skip

INIT = '"""GIAP."""\n\nfrom __future__ import annotations\n\n__version__ = "1.0.0"\n'


@pytest.fixture(scope="module")
def rd():  # type: ignore[no-untyped-def]
    return load_module()


# ------------------------------------------------------------------ a throw-away repository


class Repo:
    def __init__(self, path: Path) -> None:
        self.path = path
        path.mkdir()
        self.git("init", "-q", "-b", "main")

    def git(self, *args: str) -> str:
        out = subprocess.run(
            ["git", "-C", str(self.path), *args],
            env={**os.environ, **_GIT_ENV}, check=True, capture_output=True, text=True,
        )  # fmt: skip
        return out.stdout.strip()

    def write(self, rel: str, text: str) -> None:
        p = self.path / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")

    def commit(self, msg: str) -> str:
        self.git("add", "-A")
        self.git("commit", "-q", "--allow-empty", "-m", msg)
        return self.git("rev-parse", "HEAD")


@pytest.fixture
def repo(tmp_path: Path) -> tuple[Repo, str]:
    """A repository at C: version 1.0.0, a stats module, a viz module, docs and a readiness.py."""
    r = Repo(tmp_path / "repo")
    r.write("src/giap/__init__.py", INIT)
    r.write("src/giap/cli.py", '"""CLI."""\n')
    r.write("src/giap/stats/x.py", '"""stats."""\n')
    r.write("src/giap/report/card.py", '"""card."""\n')
    r.write("src/giap/viz/locus.py", '"""locus."""\n')
    r.write("docs/index.md", "# GIAP\n")
    r.write("README.md", "# GIAP\n")
    r.write("CITATION.cff", "version: 1.0.0\n")
    r.write("scripts/release/readiness.py", "# readiness\n")
    return r, r.commit("C")


def _bump(r: Repo, new: str = "1.0.1") -> None:
    r.write("src/giap/__init__.py", INIT.replace('"1.0.0"', f'"{new}"'))


# ------------------------------------------------------------------ the predicate on git diffs


def test_same_commit_is_reusable(rd, repo: tuple[Repo, str]) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    assert rd.reusable(r.path, c, c) is True


def test_example12_viz_docs_and_version_line_reusable(rd, repo) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    r.write("src/giap/viz/locus.py", '"""locus, restyled."""\n')
    r.write("docs/x.md", "# x\n")
    _bump(r)
    assert rd.reusable(r.path, c, r.commit("R")) is True


@pytest.mark.parametrize(
    "rel",
    ["docs/methods/new.md", "changelog.d/+fix.bugfix.md", "README.md", "CHANGELOG.md",
     "mkdocs.yml", "PROGRESS.md", "AGENTS.md", "CLAUDE.md", ".claude/agents/x.md",
     ".github/workflows/ci.yml", "scripts/dev/waitmerge.sh", "scripts/docs/gen.py",
     "src/giap/report/card.py", "src/giap/viz/locus.py"],
)  # fmt: skip
def test_each_non_stats_entry_alone_is_reusable(rd, repo, rel: str) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    r.write(rel, "changed\n")
    assert rd.reusable(r.path, c, r.commit(f"touch {rel}")) is True


@pytest.mark.parametrize(
    "rel",
    ["src/giap/cli.py", "src/giap/stats/fdr.py", "src/giap/rng.py", "src/giap/io/x.py",
     "workflow/main.nf", "pyproject.toml", "uv.lock", "tests/sim/giapsim/nulls.py",
     "tests/calibration/test_x.py", "ci/gates.yaml", "benchmarks/baseline.json",
     "tests/helpers/x.py", "CITATION.cff", "src/giap/preflight/x.py"],
)  # fmt: skip
def test_example13_any_other_path_is_not_reusable(rd, repo, rel: str) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    r.write(rel, "changed\n")
    r.write("docs/x.md", "# x\n")
    assert rd.reusable(r.path, c, r.commit(f"touch {rel}")) is False


def test_example13b_readiness_change_never_reusable(rd, repo) -> None:  # type: ignore[no-untyped-def]
    """Self-protection: a commit widening NON_STATS_PATHS cannot certify itself."""
    r, c = repo
    r.write("scripts/release/readiness.py", "# readiness\nNON_STATS_PATHS = ('src/',)\n")
    r.write("docs/x.md", "# x\n")
    assert rd.reusable(r.path, c, r.commit("widen the list")) is False


def test_example13b_holds_even_if_readiness_lists_itself(rd, repo, monkeypatch) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    r.write("scripts/release/readiness.py", "# readiness, edited\n")
    rr = r.commit("R")
    monkeypatch.setattr(rd, "NON_STATS_PATHS", (*rd.NON_STATS_PATHS, "scripts/release/"))
    assert rd.reusable(r.path, c, rr) is False


def test_example13c_move_out_of_stats_lists_both_paths(rd, repo) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    r.git("mv", "src/giap/stats/x.py", "docs/x.py")
    rr = r.commit("move stats into docs")
    names = r.git("diff", "--name-only", "--no-renames", c, rr).splitlines()
    assert sorted(names) == ["docs/x.py", "src/giap/stats/x.py"]
    assert rd.reusable(r.path, c, rr) is False


def test_a_deleted_stats_file_counts_as_changed(rd, repo) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    r.git("rm", "-q", "src/giap/stats/x.py")
    assert rd.reusable(r.path, c, r.commit("delete stats module")) is False


def test_a_move_inside_non_stats_paths_is_reusable(rd, repo) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    r.git("mv", "src/giap/viz/locus.py", "docs/locus.py")
    assert rd.reusable(r.path, c, r.commit("move viz into docs")) is True


def test_example13d_component_match(rd, repo) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    r.write("src/giap/report_stats.py", '"""not under report/."""\n')
    assert rd.reusable(r.path, c, r.commit("report_stats")) is False


# ------------------------------------------------------------------ version-line exception


def test_version_line_alone_is_reusable(rd, repo) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    _bump(r, "1.0.1rc1")
    assert rd.reusable(r.path, c, r.commit("version")) is True


@pytest.mark.parametrize(
    "new",
    [
        # a second hunk elsewhere in the file
        INIT.replace('"1.0.0"', '"1.0.1"').replace('"""GIAP."""', '"""GIAP, edited."""'),
        # one hunk, but a second added line
        INIT.replace('__version__ = "1.0.0"\n', '__version__ = "1.0.1"\nSEED_SALT = 1\n'),
        # one hunk, but a second removed line
        INIT.replace('from __future__ import annotations\n\n__version__ = "1.0.0"',
                     '__version__ = "1.0.1"'),
        # one -/+ pair that is not a version line
        INIT.replace("from __future__ import annotations", "import os"),
        # the version written with single quotes does not match ^__version__ = "[^"]+"$
        INIT.replace('__version__ = "1.0.0"', "__version__ = '1.0.1'"),
        # a trailing comment
        INIT.replace('__version__ = "1.0.0"', '__version__ = "1.0.1"  # bump'),
        # only an added line
        INIT + "SEED_SALT = 1\n",
    ],
    ids=["second-hunk", "added-line", "removed-line", "not-version", "single-quotes",
         "comment", "append"],
)  # fmt: skip
def test_any_other_change_to_init_is_not_reusable(rd, repo, new: str) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    r.write("src/giap/__init__.py", new)
    assert rd.reusable(r.path, c, r.commit("init")) is False


def test_version_line_with_a_stats_change_is_not_reusable(rd, repo) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    _bump(r)
    r.write("src/giap/stats/x.py", '"""stats, changed."""\n')
    assert rd.reusable(r.path, c, r.commit("R")) is False


def test_two_version_hunks_across_commits_still_one_hunk(rd, repo) -> None:  # type: ignore[no-untyped-def]
    """The diff is C..R, not per commit: 1.0.0 -> 1.0.1rc1 -> 1.0.1 is one -/+ pair."""
    r, c = repo
    _bump(r, "1.0.1rc1")
    r.commit("rc1")
    _bump(r, "1.0.1")
    assert rd.reusable(r.path, c, r.commit("final")) is True


def test_unknown_commit_is_not_reusable(rd, repo) -> None:  # type: ignore[no-untyped-def]
    r, c = repo
    assert rd.reusable(r.path, "0" * 40, c) is False


# ------------------------------------------------------------------ the two lists


def test_non_stats_paths_are_the_spec_list(rd) -> None:  # type: ignore[no-untyped-def]
    assert isinstance(rd.NON_STATS_PATHS, tuple)
    assert len(rd.NON_STATS_PATHS) == len(set(rd.NON_STATS_PATHS))
    assert set(rd.NON_STATS_PATHS) == SPEC_NON_STATS


def test_never_list_is_the_spec_list(rd) -> None:  # type: ignore[no-untyped-def]
    assert isinstance(rd.NEVER_PATHS, tuple)
    assert len(rd.NEVER_PATHS) == len(set(rd.NEVER_PATHS))
    assert set(rd.NEVER_PATHS) == SPEC_NEVER
    assert "scripts/release/readiness.py" in rd.NEVER_PATHS


def _under(path: str, entry: str) -> bool:
    return path.startswith(entry) if entry.endswith("/") else path == entry


def test_lists_are_disjoint_on_components(rd) -> None:  # type: ignore[no-untyped-def]
    for n in rd.NEVER_PATHS:
        assert not rd.non_stats_path(n), n
        assert not rd.non_stats_path(n + "x" if n.endswith("/") else n), n
        for e in rd.NON_STATS_PATHS:
            assert not _under(n, e) and not _under(e, n), (n, e)


@pytest.mark.parametrize("entry", sorted(SPEC_NON_STATS | SPEC_NEVER))
def test_every_entry_exists_in_the_checkout(entry: str) -> None:
    p = REPO / entry.rstrip("/")
    if entry.endswith("/"):
        assert p.is_dir(), f"{entry}: no such directory"
    else:
        assert p.is_file(), f"{entry}: no such file"


@pytest.mark.parametrize("name", ["NON_STATS_PATHS", "NEVER_PATHS"])
def test_entries_in_the_module_exist(rd, name: str) -> None:  # type: ignore[no-untyped-def]
    for entry in getattr(rd, name):
        p = REPO / entry.rstrip("/")
        assert p.is_dir() if entry.endswith("/") else p.is_file(), (name, entry)


# ------------------------------------------------------------------ component matching (pure)


def test_component_examples(rd) -> None:  # type: ignore[no-untyped-def]
    yes = ["docs/a.md", "docs/methods/x/y.md", "src/giap/report/card.py", "README.md",
           "src/giap/viz/locus.py", ".github/workflows/release.yml"]  # fmt: skip
    no = ["src/giap/report_stats.py", "src/giap/reports/x.py", "docsx/a.md", "README.md.bak",
          "README.md/x", "src/giap/vizx.py", "scripts/devtools/x.sh", "src/giap/cli.py",
          "scripts/release/readiness.py", "doc/a.md", ""]  # fmt: skip
    assert [p for p in yes if not rd.non_stats_path(p)] == []
    assert [p for p in no if rd.non_stats_path(p)] == []


SEGMENT = st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789_.-", min_size=1, max_size=8)
SUFFIX = st.lists(SEGMENT, min_size=1, max_size=4).map("/".join)
DIR_ENTRIES = sorted(e for e in SPEC_NON_STATS if e.endswith("/"))
FILE_ENTRIES = sorted(e for e in SPEC_NON_STATS if not e.endswith("/"))
NEVER_DIRS = sorted(e for e in SPEC_NEVER if e.endswith("/"))


@settings(max_examples=200, deadline=None)
@given(entry=st.sampled_from(DIR_ENTRIES), rest=SUFFIX)
def test_property_everything_under_a_dir_entry_matches(rd, entry: str, rest: str) -> None:  # type: ignore[no-untyped-def]
    assert rd.non_stats_path(entry + rest)


@settings(max_examples=200, deadline=None)
@given(entry=st.sampled_from(DIR_ENTRIES), glue=SEGMENT, rest=SUFFIX)
def test_property_a_sibling_name_never_matches(rd, entry: str, glue: str, rest: str) -> None:  # type: ignore[no-untyped-def]
    sibling = entry.rstrip("/") + glue  # "docs" + "x" -> "docsx/...", one component, not docs/
    assert not rd.non_stats_path(f"{sibling}/{rest}")
    assert not rd.non_stats_path(sibling)


@settings(max_examples=100, deadline=None)
@given(entry=st.sampled_from(FILE_ENTRIES), glue=SEGMENT, rest=SUFFIX)
def test_property_a_file_entry_matches_only_itself(rd, entry: str, glue: str, rest: str) -> None:  # type: ignore[no-untyped-def]
    assert rd.non_stats_path(entry)
    assert not rd.non_stats_path(entry + glue)
    assert not rd.non_stats_path(f"{entry}/{rest}")
    assert not rd.non_stats_path(f"{rest}/{entry}")


@settings(max_examples=200, deadline=None)
@given(entry=st.sampled_from(NEVER_DIRS), rest=SUFFIX)
def test_property_nothing_under_the_never_list_matches(rd, entry: str, rest: str) -> None:  # type: ignore[no-untyped-def]
    assert not rd.non_stats_path(entry + rest)
