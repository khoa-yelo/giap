"""M11-14h: the schema diff test and ``fetch-depth: 0`` in CI.

docs/specs/M11-14-container-install.md §4 "Schema diff test (14h)" and acceptance criterion 14:
"a removed column with an unchanged major fails, and an added column with an unchanged schema
minor fails (fixture trees). No ``v*`` tag skips with ``no previous release``; a shallow fixture
repo fails with ``fetch tags``"; "14h sets ``fetch-depth: 0`` and ``fetch-tags: true`` on the CI
job that runs it".

Interface pinned here (test-author reading, beyond the spec text):

- ``scripts/release/schema_diff.py`` (loaded by path, like ``readiness.py``) exposes
  ``check(root: Path) -> list[str]``: the problems, ``[]`` when HEAD is compatible with the
  newest ``v*`` tag reachable from HEAD. It reads the old tree with ``git show <tag>:<path>``
  only (no network, no checkout). It raises ``NoPreviousRelease`` (message ``no previous
  release``) when no ``v*`` tag is reachable, and returns a problem containing ``fetch tags:
  actions/checkout with fetch-depth: 0`` when ``git rev-parse --is-shallow-repository`` prints
  ``true`` (checked first: a shallow clone has no tags either).
- Compared: the ``fields[].name`` of every ``src/giap/schemas/tables/*.schema.json`` (a table
  file that disappears removes its columns, a new one adds them) and the top-level
  ``properties`` keys of ``src/giap/schemas/run_manifest.schema.json``. The package version is
  ``__version__`` of ``src/giap/__init__.py``; ``output_schema_version`` is read from
  ``src/giap/schemas/datapackage.json``.
- A removal (a rename is a removal plus an addition) needs a larger ``__version__`` major; an
  addition needs a larger ``output_schema_version`` (minor, or major). Each problem names the
  table (file stem, or ``run_manifest``) and the column.

The last tests run the check on this repository (skipped until the first ``v*`` tag) and read
``.github/workflows/ci.yml``.
"""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "release" / "schema_diff.py"
CI_YML = REPO / ".github" / "workflows" / "ci.yml"
FETCH_TAGS = "fetch tags: actions/checkout with fetch-depth: 0"

_GIT_ENV = {
    "GIT_AUTHOR_NAME": "Fixture",
    "GIT_AUTHOR_EMAIL": "fixture@example.invalid",
    "GIT_COMMITTER_NAME": "Fixture",
    "GIT_COMMITTER_EMAIL": "fixture@example.invalid",
    "GIT_CONFIG_GLOBAL": os.devnull,
    "GIT_CONFIG_SYSTEM": os.devnull,
}

GENES = ["locus_tag", "contig", "start", "end"]
MANIFEST_KEYS = ["schema_version", "output_schema_version", "giap_version", "operon_source"]


def _load() -> ModuleType:
    spec = importlib.util.spec_from_file_location("giap_release_schema_diff", SCRIPT)
    assert spec is not None and spec.loader is not None, f"{SCRIPT} missing"
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod
    spec.loader.exec_module(mod)
    return mod


def git(repo: Path, *args: str) -> str:
    out = subprocess.run(
        ["git", "-C", str(repo), *args],
        env={**os.environ, **_GIT_ENV},
        check=True,
        capture_output=True,
        text=True,
    )
    return out.stdout.strip()


def _write(root: Path, rel: str, text: str) -> None:
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_tree(
    root: Path,
    *,
    version: str = "1.0.0",
    output: str = "1.0",
    tables: dict[str, list[str]] | None = None,
    manifest: list[str] | None = None,
) -> None:
    """The four files the check reads; ``tables`` maps a stem to its column names."""
    tables = {"genes": GENES, "islands": ["island_id", "p_value"]} if tables is None else tables
    _write(root, "src/giap/__init__.py", f'"""GIAP."""\n\n__version__ = "{version}"\n')
    tdir = root / "src" / "giap" / "schemas" / "tables"
    if tdir.is_dir():
        for old in tdir.glob("*.schema.json"):
            old.unlink()
    for stem, cols in tables.items():
        doc = {"name": stem, "fields": [{"name": c, "type": "string"} for c in cols]}
        _write(root, f"src/giap/schemas/tables/{stem}.schema.json", json.dumps(doc, indent=2))
    props = {k: {"type": "string"} for k in (MANIFEST_KEYS if manifest is None else manifest)}
    _write(
        root,
        "src/giap/schemas/run_manifest.schema.json",
        json.dumps({"type": "object", "required": sorted(props), "properties": props}),
    )
    package = {"name": "giap-run", "output_schema_version": output, "resources": []}
    _write(root, "src/giap/schemas/datapackage.json", json.dumps(package, indent=2))


def commit(root: Path, msg: str, tag: str | None = None) -> None:
    git(root, "add", "-A")
    git(root, "commit", "-q", "--allow-empty", "-m", msg)
    if tag:
        git(root, "tag", tag)


def released(tmp: Path, **kw: Any) -> Path:
    """A repo whose ``v1.0.0`` tag holds the baseline tree; HEAD = the tag."""
    root = tmp / "repo"
    root.mkdir()
    git(root, "init", "-q", "-b", "main")
    write_tree(root, **kw)
    commit(root, "release 1.0.0", tag="v1.0.0")
    return root


def head(root: Path, **kw: Any) -> list[str]:
    write_tree(root, **kw)
    commit(root, "change")
    problems = _load().check(root)
    assert isinstance(problems, list)
    return [str(p) for p in problems]


# ---------------------------------------------------------------- fixture trees


def test_unchanged_tree_passes(tmp_path: Path) -> None:
    root = released(tmp_path)
    assert head(root, version="1.0.1") == []


def test_removed_column_with_unchanged_major_fails(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(
        root,
        version="1.1.0",
        output="1.1",
        tables={"genes": GENES[:-1], "islands": ["island_id", "p_value"]},
    )
    assert out, "a removed column must fail without a major bump"
    assert any("genes" in p and "end" in p for p in out), out


def test_removed_column_with_a_major_bump_passes(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(
        root,
        version="2.0.0",
        output="2.0",
        tables={"genes": GENES[:-1], "islands": ["island_id", "p_value"]},
    )
    assert out == []


def test_renamed_column_with_unchanged_major_fails(tmp_path: Path) -> None:
    root = released(tmp_path)
    cols = [*GENES[:-1], "stop"]
    out = head(
        root,
        version="1.1.0",
        output="1.1",
        tables={"genes": cols, "islands": ["island_id", "p_value"]},
    )
    assert any("genes" in p and "end" in p for p in out), out


def test_removed_table_with_unchanged_major_fails(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(root, version="1.1.0", output="1.1", tables={"genes": GENES})
    assert any("islands" in p for p in out), out


def test_added_column_with_unchanged_schema_minor_fails(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(
        root,
        version="1.1.0",
        output="1.0",
        tables={"genes": [*GENES, "operon_source"], "islands": ["island_id", "p_value"]},
    )
    assert out, "an added column must fail without an output_schema_version minor bump"
    assert any("genes" in p and "operon_source" in p for p in out), out


def test_added_column_with_a_schema_minor_bump_passes(tmp_path: Path) -> None:
    root = released(tmp_path)
    out = head(
        root,
        version="1.1.0",
        output="1.1",
        tables={"genes": [*GENES, "operon_source"], "islands": ["island_id", "p_value"]},
    )
    assert out == []


def test_added_table_needs_a_schema_minor_bump(tmp_path: Path) -> None:
    root = released(tmp_path)
    tables = {"genes": GENES, "islands": ["island_id", "p_value"], "tus": ["tu_id"]}
    out = head(root, version="1.1.0", output="1.0", tables=tables)
    assert any("tus" in p for p in out), out


def test_manifest_key_removed_or_added(tmp_path: Path) -> None:
    """The manifest schema counts like a table (14d's ``operon_source`` is such an addition)."""
    root = released(tmp_path, manifest=MANIFEST_KEYS[:-1])
    added = head(root, version="1.1.0", output="1.0", manifest=MANIFEST_KEYS)
    assert any("operon_source" in p for p in added), added
    removed = head(root, version="1.1.0", output="1.1", manifest=MANIFEST_KEYS[:-2])
    assert any("output_schema_version" in p for p in removed), removed


def test_compares_with_the_newest_reachable_tag(tmp_path: Path) -> None:
    """v1.1.0 already has the new column: HEAD against v1.1.0 passes; a tag on a side branch
    that HEAD cannot reach is ignored."""
    root = released(tmp_path)
    write_tree(
        root,
        version="1.1.0",
        output="1.1",
        tables={"genes": [*GENES, "x"], "islands": ["island_id", "p_value"]},
    )
    commit(root, "release 1.1.0", tag="v1.1.0")
    git(root, "checkout", "-q", "-b", "side", "v1.0.0")
    write_tree(root, version="9.0.0", output="9.0", tables={"genes": ["only"]})
    commit(root, "side", tag="v9.0.0")
    git(root, "checkout", "-q", "main")
    out = head(
        root,
        version="1.1.1",
        output="1.1",
        tables={"genes": [*GENES, "x"], "islands": ["island_id", "p_value"]},
    )
    assert out == []


def test_no_v_tag_skips_with_no_previous_release(tmp_path: Path) -> None:
    mod = _load()
    root = tmp_path / "repo"
    root.mkdir()
    git(root, "init", "-q", "-b", "main")
    write_tree(root)
    commit(root, "first", tag="not-a-release")
    with pytest.raises(mod.NoPreviousRelease, match="no previous release"):
        mod.check(root)


def test_shallow_clone_fails_with_fetch_tags(tmp_path: Path) -> None:
    """Shallow and no tags: a failure naming the CI fix, never the skip."""
    src = released(tmp_path)
    write_tree(src, version="1.0.1")
    commit(src, "after the release")
    clone = tmp_path / "shallow"
    subprocess.run(
        ["git", "clone", "-q", "--depth", "1", "--no-tags", src.as_uri(), str(clone)],
        env={**os.environ, **_GIT_ENV},
        check=True,
        capture_output=True,
    )
    assert git(clone, "rev-parse", "--is-shallow-repository") == "true"
    out = [str(p) for p in _load().check(clone)]
    assert any(FETCH_TAGS in p for p in out), out


# ---------------------------------------------------------------- this repository and CI


def test_this_repository_is_schema_compatible_with_its_last_release() -> None:
    """The schema diff test itself: skipped before the first ``v*`` tag."""
    mod = _load()
    try:
        problems = mod.check(REPO)
    except mod.NoPreviousRelease as exc:
        pytest.skip(str(exc))
    assert problems == []


def _jobs() -> dict[str, Any]:
    doc = yaml.safe_load(CI_YML.read_text(encoding="utf-8"))
    jobs = doc["jobs"]
    assert isinstance(jobs, dict)
    return jobs


@pytest.mark.parametrize("job", ["test", "test-compat"])
def test_ci_jobs_running_unit_tests_fetch_full_history_and_tags(job: str) -> None:
    """The jobs that run ``tests/unit`` (and so the schema diff test) check out with
    ``fetch-depth: 0`` and ``fetch-tags: true``."""
    steps = _jobs()[job]["steps"]
    checkout = [s for s in steps if str(s.get("uses", "")).startswith("actions/checkout@")]
    assert len(checkout) == 1, steps
    opts = checkout[0].get("with") or {}
    assert opts.get("fetch-depth") == 0, opts
    assert opts.get("fetch-tags") is True, opts
