"""``scripts/ci/check_changelog.py BASE`` (M11-12 definition item 1, examples 1-8, edge cases).

Every case builds a throwaway git repository in ``tmp_path``: ``main`` holds the base commit, a
``feature`` branch holds the PR commits, and the script runs from the repository with ``main`` as
BASE (the ``check_roles.py`` pattern). Labels come from a fake ``GITHUB_EVENT_PATH`` JSON.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts/ci/check_changelog.py"
MISSING = "changelog fragment missing"
EXEMPT_LABEL = "exempt: label no-changelog"
EXEMPT_TYPES = "exempt: refactor/test/ci/chore only"
TYPES = ("feature", "bugfix", "doc", "removal", "misc")


def _git_env() -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not k.startswith(("GIT_", "GITHUB_"))}
    env.update(
        GIT_AUTHOR_NAME="t",
        GIT_AUTHOR_EMAIL="t@example.invalid",
        GIT_COMMITTER_NAME="t",
        GIT_COMMITTER_EMAIL="t@example.invalid",
        GIT_CONFIG_GLOBAL=os.devnull,
        GIT_CONFIG_NOSYSTEM="1",
    )
    return env


@dataclass
class Repo:
    path: Path

    def git(self, *args: str) -> str:
        return subprocess.run(
            ["git", "-c", "commit.gpgsign=false", *args],
            cwd=self.path,
            env=_git_env(),
            capture_output=True,
            text=True,
            check=True,
        ).stdout

    def write(self, files: Mapping[str, str]) -> None:
        for rel, text in files.items():
            p = self.path / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(text)

    def commit(self, subject: str, files: Mapping[str, str] | None = None) -> None:
        if files:
            self.write(files)
        self.git("add", "-A")
        self.git("commit", "--allow-empty", "-q", "-m", subject, "-m", "Agent-Role: implementer")


@dataclass
class Result:
    code: int
    out: str


def _repo(tmp_path: Path, base_files: Mapping[str, str] | None = None) -> Repo:
    root = tmp_path / "repo"
    root.mkdir()
    r = Repo(root)
    r.git("init", "-q", "-b", "main")
    r.commit(
        "chore: base",
        {
            "src/gippa/cli.py": "x = 1\n",
            "src/gippa/x.py": "y = 1\n",
            "docs/user-guide/run.md": "# run\n",
            "workflow/main.nf": "// main\n",
            "workflow/tests/old.nf.test": "// t\n",
            "changelog.d/+m00-scaffold.feature.md": "Scaffold.\n",
            **(base_files or {}),
        },
    )
    r.git("checkout", "-q", "-b", "feature")
    return r


def _run(
    repo: Repo,
    labels: list[str] | None = None,
    event: dict[str, object] | None = None,
    base: str = "main",
) -> Result:
    env = _git_env()
    if labels is not None or event is not None:
        payload = (
            event
            if event is not None
            else {"pull_request": {"number": 341, "labels": [{"name": n} for n in labels or []]}}
        )
        ev = repo.path.parent / "event.json"
        ev.write_text(json.dumps(payload))
        env["GITHUB_EVENT_PATH"] = str(ev)
    p = subprocess.run(
        [sys.executable, str(SCRIPT), base],
        cwd=repo.path,
        env=env,
        capture_output=True,
        text=True,
    )
    return Result(p.returncode, p.stdout + p.stderr)


def _pr(tmp_path: Path, files: Mapping[str, str], subject: str = "feat: change") -> Repo:
    r = _repo(tmp_path)
    r.commit(subject, files)
    return r


def test_script_exists() -> None:
    assert SCRIPT.is_file(), "scripts/ci/check_changelog.py is missing (M11-12 item 1)"


# ---- the spec's examples 1-8 -------------------------------------------------------------------


def test_example1_code_with_fragment_passes(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/cli.py": "x = 2\n", "changelog.d/341.feature.md": "New flag.\n"})
    res = _run(r)
    assert res.code == 0, res.out


def test_example2_code_without_fragment_fails_with_message(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/cli.py": "x = 2\n"})
    res = _run(r)
    assert res.code == 1, res.out
    assert MISSING in res.out
    assert "src/gippa/cli.py" in res.out
    assert "changelog.d/" in res.out
    assert "no-changelog" in res.out


def test_example3_docs_only_passes(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"docs/user-guide/run.md": "# run\nmore\n"}, subject="docs: run page")
    res = _run(r)
    assert res.code == 0, res.out
    assert MISSING not in res.out


def test_example4_workflow_tests_only_passes(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"workflow/tests/x.nf.test": "// x\n"}, subject="feat: nf-test")
    res = _run(r)
    assert res.code == 0, res.out


def test_example5_label_no_changelog_exempts(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/cli.py": "x = 2\n"})
    res = _run(r, labels=["enhancement", "no-changelog"])
    assert res.code == 0, res.out
    assert EXEMPT_LABEL in res.out


def test_example6_refactor_only_commits_exempt(tmp_path: Path) -> None:
    r = _repo(tmp_path)
    r.commit("refactor: split x", {"src/gippa/x.py": "y = 2\n"})
    r.commit("refactor(cli): rename", {"src/gippa/x.py": "y = 3\n"})
    res = _run(r)
    assert res.code == 0, res.out
    assert EXEMPT_TYPES in res.out


def test_example7_empty_fragment_fails(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n", "changelog.d/342.bugfix.md": ""})
    res = _run(r)
    assert res.code == 1, res.out
    assert "empty fragment" in res.out


def test_example8_misnamed_fragment_fails(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n", "changelog.d/notes.txt": "Notes.\n"})
    res = _run(r)
    assert res.code == 1, res.out
    assert MISSING in res.out


# ---- needs a fragment ---------------------------------------------------------------------------


def test_workflow_code_needs_a_fragment(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"workflow/main.nf": "// changed\n"})
    assert _run(r).code == 1


def test_workflow_code_with_fragment_passes(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"workflow/main.nf": "// changed\n", "changelog.d/7.misc.md": "Tweak.\n"})
    assert _run(r).code == 0


def test_tests_scripts_and_docs_need_no_fragment(tmp_path: Path) -> None:
    r = _pr(
        tmp_path,
        {
            "tests/unit/test_a.py": "def test_a():\n    pass\n",
            "scripts/ci/x.py": "print(1)\n",
            "README.md": "# r\n",
            "workflow/tests/old.nf.test": "// changed\n",
        },
        subject="feat: tooling",
    )
    assert _run(r).code == 0


def test_no_changes_passes(tmp_path: Path) -> None:
    r = _repo(tmp_path)
    assert _run(r).code == 0


def test_deleted_src_file_needs_a_fragment(tmp_path: Path) -> None:
    r = _repo(tmp_path)
    r.git("rm", "-q", "src/gippa/x.py")
    r.commit("feat: drop x")
    res = _run(r)
    assert res.code == 1, res.out
    assert "src/gippa/x.py" in res.out


def test_rename_into_src_counts_at_new_path(tmp_path: Path) -> None:
    r = _repo(tmp_path, base_files={"tools/helper.py": "def h():\n    return 1\n" * 5})
    r.git("mv", "tools/helper.py", "src/gippa/helper.py")
    r.commit("feat: move helper")
    res = _run(r)
    assert res.code == 1, res.out
    assert "src/gippa/helper.py" in res.out


# ---- has a fragment -----------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name", [f"341.{t}.md" for t in TYPES] + ["+m11-release.feature.md", "+fix-x.bugfix.md"]
)
def test_valid_fragment_names_pass(tmp_path: Path, name: str) -> None:
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n", f"changelog.d/{name}": "Text.\n"})
    res = _run(r)
    assert res.code == 0, res.out


@pytest.mark.parametrize(
    "name",
    [
        "341.feat.md",  # unknown type
        "341.feature.txt",  # wrong suffix
        "abc.feature.md",  # neither a PR number nor +slug
        "341.md",  # no type
        "README.md",
        "sub/341.feature.md",  # not directly under changelog.d/
    ],
)
def test_invalid_fragment_names_fail(tmp_path: Path, name: str) -> None:
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n", f"changelog.d/{name}": "Text.\n"})
    res = _run(r)
    assert res.code == 1, res.out


def test_fragment_outside_changelog_dir_does_not_count(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n", "docs/341.feature.md": "Text.\n"})
    assert _run(r).code == 1


def test_whitespace_only_fragment_is_empty(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n", "changelog.d/342.bugfix.md": " \n\n"})
    res = _run(r)
    assert res.code == 1, res.out
    assert "empty fragment" in res.out


def test_one_valid_fragment_beside_a_misnamed_file_passes(tmp_path: Path) -> None:
    r = _pr(
        tmp_path,
        {
            "src/gippa/x.py": "y = 2\n",
            "changelog.d/notes.txt": "n\n",
            "changelog.d/343.feature.md": "Text.\n",
        },
    )
    assert _run(r).code == 0


def test_fragment_already_on_base_does_not_count(tmp_path: Path) -> None:
    # the base already holds +m00-scaffold.feature.md; the PR must *add* one
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n"})
    assert _run(r).code == 1


def test_fragment_in_a_later_commit_counts(tmp_path: Path) -> None:
    r = _repo(tmp_path)
    r.commit("feat: x", {"src/gippa/x.py": "y = 2\n"})
    r.commit("docs: fragment", {"changelog.d/344.feature.md": "Text.\n"})
    assert _run(r).code == 0


def test_fragment_added_then_deleted_does_not_count(tmp_path: Path) -> None:
    r = _repo(tmp_path)
    r.commit("feat: x", {"src/gippa/x.py": "y = 2\n", "changelog.d/344.feature.md": "Text.\n"})
    r.git("rm", "-q", "changelog.d/344.feature.md")
    r.commit("feat: oops")
    assert _run(r).code == 1


# ---- exemptions ---------------------------------------------------------------------------------


def test_other_labels_do_not_exempt(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n"})
    res = _run(r, labels=["documentation", "changelog"])
    assert res.code == 1, res.out
    assert EXEMPT_LABEL not in res.out


def test_event_without_pull_request_is_not_exempt(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n"})
    res = _run(r, event={"ref": "refs/heads/main", "labels": [{"name": "no-changelog"}]})
    assert res.code == 1, res.out


def test_missing_event_file_is_not_exempt(tmp_path: Path) -> None:
    r = _pr(tmp_path, {"src/gippa/x.py": "y = 2\n"})
    env_unset = _run(r)  # GITHUB_EVENT_PATH removed by _git_env
    assert env_unset.code == 1, env_unset.out


@pytest.mark.parametrize(
    "subjects",
    [
        ["test: add cases"],
        ["ci: guard"],
        ["chore: tidy"],
        ["refactor: split"],
        ["test(m11): a", "ci(guards): b", "chore(deps): c", "refactor(cli): d"],
    ],
)
def test_maintenance_only_subjects_exempt(tmp_path: Path, subjects: list[str]) -> None:
    r = _repo(tmp_path)
    for i, s in enumerate(subjects):
        r.commit(s, {"src/gippa/x.py": f"y = {i + 10}\n"})
    res = _run(r)
    assert res.code == 0, res.out
    assert EXEMPT_TYPES in res.out


@pytest.mark.parametrize(
    "subjects",
    [
        ["refactor: split", "feat: new flag"],
        ["fix: bug"],
        ["docs: page"],
        ["perf: faster"],
        ["Refactor the cli"],  # not a conventional subject
        ["refactoring: x"],  # a different type
        ["chore release 1.0.0"],  # no colon
    ],
)
def test_other_subjects_do_not_exempt(tmp_path: Path, subjects: list[str]) -> None:
    r = _repo(tmp_path)
    for i, s in enumerate(subjects):
        r.commit(s, {"src/gippa/x.py": f"y = {i + 10}\n"})
    res = _run(r)
    assert res.code == 1, res.out
    assert EXEMPT_TYPES not in res.out


def test_release_prep_deleting_fragments_is_exempt(tmp_path: Path) -> None:
    """Edge case: release prep drops all fragments, bumps the version (``chore(release):``)."""
    r = _repo(tmp_path, base_files={"src/gippa/__init__.py": '__version__ = "0.0.1"\n'})
    r.git("rm", "-q", "changelog.d/+m00-scaffold.feature.md")
    r.commit(
        "chore(release): 1.0.0",
        {"src/gippa/__init__.py": '__version__ = "1.0.0"\n', "CHANGELOG.md": "# 1.0.0\n"},
    )
    res = _run(r)
    assert res.code == 0, res.out
    assert EXEMPT_TYPES in res.out


# ---- merge commits from main --------------------------------------------------------------------


def test_main_changes_merged_into_the_branch_do_not_count(tmp_path: Path) -> None:
    """Edge case: ``BASE...HEAD`` uses the merge base; main's own src/ changes do not count."""
    r = _repo(tmp_path)
    r.commit("docs: page", {"docs/user-guide/run.md": "# run\nmore\n"})
    r.git("checkout", "-q", "main")
    r.commit("feat: main moves on", {"src/gippa/cli.py": "x = 99\n"})
    r.git("checkout", "-q", "feature")
    r.git("merge", "-q", "--no-ff", "--no-edit", "main")
    res = _run(r)
    assert res.code == 0, res.out


def test_merge_commit_subject_does_not_break_the_type_exemption(tmp_path: Path) -> None:
    r = _repo(tmp_path)
    r.commit("refactor: split", {"src/gippa/x.py": "y = 2\n"})
    r.git("checkout", "-q", "main")
    r.commit("docs: main moves on", {"docs/other.md": "o\n"})
    r.git("checkout", "-q", "feature")
    r.git("merge", "-q", "--no-ff", "--no-edit", "main")  # subject "Merge branch 'main' ..."
    res = _run(r)
    assert res.code == 0, res.out
    assert EXEMPT_TYPES in res.out


# ---- properties ---------------------------------------------------------------------------------

_SLUG = st.from_regex(r"[a-z0-9][a-z0-9-]{0,15}", fullmatch=True)
_PROP = settings(
    max_examples=15, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)


@_PROP
@given(
    stem=st.one_of(st.integers(1, 99999).map(str), _SLUG.map(lambda s: "+" + s)),
    kind=st.sampled_from(TYPES),
)
def test_property_any_rule_conforming_fragment_passes(
    tmp_path_factory: pytest.TempPathFactory, stem: str, kind: str
) -> None:
    r = _pr(
        tmp_path_factory.mktemp("p"),
        {"src/gippa/x.py": "y = 2\n", f"changelog.d/{stem}.{kind}.md": "Text.\n"},
    )
    res = _run(r)
    assert res.code == 0, res.out


@_PROP
@given(
    number=st.integers(1, 99999),
    kind=st.text(alphabet="abcdefghijklmnopqrstuvwxyz", min_size=1, max_size=8).filter(
        lambda k: k not in TYPES
    ),
)
def test_property_unknown_type_fails(
    tmp_path_factory: pytest.TempPathFactory, number: int, kind: str
) -> None:
    r = _pr(
        tmp_path_factory.mktemp("p"),
        {"src/gippa/x.py": "y = 2\n", f"changelog.d/{number}.{kind}.md": "Text.\n"},
    )
    res = _run(r)
    assert res.code == 1, res.out


@_PROP
@given(
    sub=st.from_regex(r"[a-z][a-z0-9_]{0,8}", fullmatch=True),
    top=st.sampled_from(
        ["docs", "tests/unit", "scripts", "workflow/tests", "validation", "changelog.d"]
    ),
)
def test_property_non_code_paths_never_need_a_fragment(
    tmp_path_factory: pytest.TempPathFactory, sub: str, top: str
) -> None:
    r = _pr(tmp_path_factory.mktemp("p"), {f"{top}/{sub}.txt": "z\n"}, subject="feat: other")
    res = _run(r)
    assert res.code == 0, res.out
