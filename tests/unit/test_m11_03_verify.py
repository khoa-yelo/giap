"""``readiness.py verify``: the commit rule and the re-check of the ``repo`` rows (issue M11-03,
examples 11-12; Q-L4; edge case "Detached HEAD").

Spec: docs/specs/M11-03-release-readiness.md, "Commit rule (verify)". ``verify`` passes when
HEAD = ``commit`` or every path changed in ``commit..HEAD`` is release paperwork
(``validation/release/**``, ``PROGRESS.md``, ``docs/**`` but not ``docs/methods/**``,
``README.md``, ``CITATION.cff``, ``CHANGELOG.md``, ``changelog.d/**``); any other path gives
``fail (stale: <first path>)`` and exit 1. ``paperwork_path`` is the module's pure predicate for
one repository-relative path (property tests below).
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from helpers.release_readiness import Fixture, git, load_module, metrics_doc, progress_md
from hypothesis import given
from hypothesis import strategies as st


@pytest.fixture
def computed(tmp_path: Path) -> Fixture:
    fx = Fixture.build(tmp_path)
    proc, doc = fx.compute()
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert doc is not None and doc["commit"] == fx.r
    return fx


def _out(p: object) -> str:
    return getattr(p, "stdout", "") + getattr(p, "stderr", "")


def test_verify_head_is_commit(computed: Fixture) -> None:
    proc = computed.verify()
    assert proc.returncode == 0, _out(proc)


def test_example11_paperwork_commit_passes(computed: Fixture) -> None:
    md = computed.out / "READINESS.md"
    md.write_text(md.read_text(encoding="utf-8") + "\nSigned off.\n", encoding="utf-8")
    computed.commit("readiness paperwork")
    proc = computed.verify()
    assert proc.returncode == 0, _out(proc)


def test_verify_paperwork_mix_passes(computed: Fixture) -> None:
    computed.write_text("docs/user-guide/install.md", "# Install\n\nEdited.\n")
    computed.write_text("CHANGELOG.md", (computed.repo / "CHANGELOG.md").read_text() + "\nMore.\n")
    computed.write_text("README.md", "# GIPPA 1.0\n")
    computed.commit("docs")
    assert computed.verify().returncode == 0


def test_example12_code_commit_is_stale(computed: Fixture) -> None:
    computed.write_text("src/gippa/cli.py", '"""CLI (fixture), changed."""\n')
    computed.commit("touch cli")
    proc = computed.verify()
    assert proc.returncode == 1
    assert "fail (stale: src/gippa/cli.py)" in _out(proc)


def test_methods_note_is_not_paperwork(computed: Fixture) -> None:
    computed.write_text("docs/methods/bh.md", "# BH\n")
    computed.commit("methods note")
    proc = computed.verify()
    assert proc.returncode == 1 and "stale: docs/methods/bh.md" in _out(proc)


def test_stale_path_over_two_commits(computed: Fixture) -> None:
    computed.write_text("PROGRESS.md", (computed.repo / "PROGRESS.md").read_text() + "\nx\n")
    computed.commit("progress")
    computed.write_text("workflow/main.nf", "// changed\n")
    computed.commit("workflow")
    proc = computed.verify()
    assert proc.returncode == 1 and "stale: workflow/main.nf" in _out(proc)


def test_verify_rechecks_repo_rows(computed: Fixture) -> None:
    # PROGRESS.md is paperwork, so the commit rule passes, but H1 now fails on the new HEAD
    computed.write_text("PROGRESS.md", progress_md(range(1, 9)))
    computed.commit("drop gates")
    proc = computed.verify()
    assert proc.returncode == 1, _out(proc)
    assert "H1" in _out(proc)


def test_verify_does_not_trust_tampered_json(tmp_path: Path) -> None:
    fx = Fixture.build(tmp_path)
    fx.write_json("validation/pilot/metrics/metrics.json", metrics_doc(share=0.9949))
    proc, out = fx.compute()
    assert proc.returncode == 1 and out is not None
    for r in out["rows"]:
        if r["id"] == "R1":
            r["status"] = "pass"
    out["overall"] = "pass"
    p = fx.out / "readiness.json"
    p.write_text(json.dumps(out), encoding="utf-8")
    assert fx.verify(p).returncode == 1


def test_verify_tag_created_after_compute_fails(computed: Fixture) -> None:
    git(computed.repo, "tag", "v1.0.0")
    proc = computed.verify()
    assert proc.returncode == 1 and "tag v1.0.0 exists" in _out(proc)


def test_verify_unknown_commit_fails(computed: Fixture) -> None:
    p = computed.out / "readiness.json"
    doc = json.loads(p.read_text(encoding="utf-8"))
    doc["commit"] = "0" * 40
    p.write_text(json.dumps(doc), encoding="utf-8")
    assert computed.verify(p).returncode == 1


# ------------------------------------------------------------------ paperwork_path (pure)

ALLOWED = [
    "validation/release/1.0.0/READINESS.md",
    "validation/release/1.0.0/readiness.json",
    "PROGRESS.md",
    "docs/release.md",
    "docs/user-guide/install.md",
    "docs/adr/ADR-0060-x.md",
    "README.md",
    "CITATION.cff",
    "CHANGELOG.md",
    "changelog.d/+m11.feature.md",
]
DENIED = [
    "src/gippa/cli.py",
    "docs/methods/bh.md",
    "docs/methods/sub/x.md",
    "validation/pilot/metrics/metrics.json",
    "validation/releases/x",
    "tests/unit/test_x.py",
    "workflow/main.nf",
    "pyproject.toml",
    "PROGRESS.md.bak",
    "sub/README.md",
    "docsx/a.md",
    "ci/gates.yaml",
]


@pytest.mark.parametrize("path", ALLOWED)
def test_paperwork_allowed(path: str) -> None:
    assert load_module().paperwork_path(path) is True


@pytest.mark.parametrize("path", DENIED)
def test_paperwork_denied(path: str) -> None:
    assert load_module().paperwork_path(path) is False


_seg = st.text(alphabet="abcdefghijklmnopqrstuvwxyz0123456789_-.", min_size=1, max_size=12).filter(
    lambda s: s not in (".", "..")
)
_rel = st.lists(_seg, min_size=1, max_size=4).map("/".join)


@given(_rel)
def test_property_validation_release_and_changelog_d_trees(rest: str) -> None:
    mod = load_module()
    assert mod.paperwork_path(f"validation/release/{rest}")
    assert mod.paperwork_path(f"changelog.d/{rest}")


@given(_rel)
def test_property_docs_tree_except_methods(rest: str) -> None:
    mod = load_module()
    assert mod.paperwork_path(f"docs/methods/{rest}") is False
    if rest.split("/")[0] != "methods":
        assert mod.paperwork_path(f"docs/{rest}") is True


@given(_rel)
def test_property_code_trees_never_paperwork(rest: str) -> None:
    mod = load_module()
    for top in ("src", "workflow", "scripts", "tests", "ci", "benchmarks", ".github"):
        assert mod.paperwork_path(f"{top}/{rest}") is False
