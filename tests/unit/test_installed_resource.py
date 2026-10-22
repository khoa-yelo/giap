"""M8-06a item 7: ``giap.setupdb.installed.resource(db_dir, name)`` (the ``cog_table`` rule of
M4-15 decision 5, generalised; docs/specs/M8-06-function-step-and-run.md).

``resource`` returns ``db_dir / SOURCES[name].path`` when the setup manifest lists ``name`` and the
file exists, else None; ``cog_table(db) == resource(db, "cog2020_def")``. The test DB is
``helpers.function_db`` (never the real ``$GIAP_DB_DIR``).
"""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest
from helpers.function_db import FUNCTION_RESOURCES, ITEMS, PATHS, manifest_doc, write_function_db
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st


def _xrefs() -> pd.DataFrame:
    return pd.DataFrame({"locus_tag": ["g1", "g1", "g1"], "db": ["COG", "GO", "KEGG"],
                         "id": ["COG0001", "GO:0001201", "K00101"]})  # fmt: skip


def test_function_resources_constant() -> None:
    from giap.setupdb.installed import FUNCTION_RESOURCES as LIVE

    assert LIVE == FUNCTION_RESOURCES


def test_paths_come_from_db_sources_yaml() -> None:
    """The helper's layout is ``db_sources.yaml``'s (what ``resource`` joins onto ``db_dir``)."""
    from giap.setupdb.items import SOURCES

    for name in ITEMS:
        assert SOURCES[name].path == PATHS[name], name


@pytest.mark.parametrize("name", [*FUNCTION_RESOURCES, "cog2020_def"])
def test_listed_and_present_gives_the_path(name: str, tmp_path: Path) -> None:
    from giap.setupdb.installed import resource
    from giap.setupdb.items import SOURCES

    db = write_function_db(tmp_path, _xrefs())
    got = resource(tmp_path, name)
    assert got == tmp_path / SOURCES[name].path
    assert got == db.paths[name]


@pytest.mark.parametrize("name", FUNCTION_RESOURCES)
def test_listed_but_absent_gives_none(name: str, tmp_path: Path) -> None:
    """Degenerate case: listed in the manifest, file deleted -> not installed."""
    from giap.setupdb.installed import resource

    db = write_function_db(tmp_path, _xrefs())
    path = db.paths[name]
    assert path is not None
    path.unlink()
    assert resource(tmp_path, name) is None


@pytest.mark.parametrize("name", FUNCTION_RESOURCES)
def test_present_but_not_listed_gives_none(name: str, tmp_path: Path) -> None:
    from giap.setupdb.installed import resource

    others = [n for n in ITEMS if n != name]
    write_function_db(tmp_path, _xrefs(), listed=others)
    assert (tmp_path / PATHS[name]).is_file()
    assert resource(tmp_path, name) is None


def test_no_db_dir_or_no_manifest_gives_none(tmp_path: Path) -> None:
    from giap.setupdb.installed import resource

    assert resource(None, "go_basic") is None
    (tmp_path / "go").mkdir()
    (tmp_path / PATHS["go_basic"]).write_text("format-version: 1.2\n", encoding="utf-8")
    assert resource(tmp_path, "go_basic") is None  # no db_manifest.json


def test_provisional_manifest_form_is_read(tmp_path: Path) -> None:
    """``read_manifest``'s provisional ``{"<name>": {"version": ...}}`` form counts too."""
    from giap.setupdb.installed import resource

    write_function_db(tmp_path, _xrefs(), items=["cog2020_fun"])
    doc = {"cog2020_fun": {"version": "2020", "path": PATHS["cog2020_fun"]}}
    (tmp_path / "db_manifest.json").write_text(json.dumps(doc) + "\n", encoding="utf-8")
    assert resource(tmp_path, "cog2020_fun") == tmp_path / PATHS["cog2020_fun"]
    assert resource(tmp_path, "go_basic") is None


def test_status_not_ok_is_not_installed(tmp_path: Path) -> None:
    from giap.setupdb.installed import resource

    write_function_db(tmp_path, _xrefs())
    doc = manifest_doc(ITEMS)
    for item in doc["items"]:  # type: ignore[attr-defined]
        if item["name"] == "go_basic":
            item["status"] = "failed"
    (tmp_path / "db_manifest.json").write_text(json.dumps(doc) + "\n", encoding="utf-8")
    assert resource(tmp_path, "go_basic") is None
    assert resource(tmp_path, "goslim_generic") is not None


@pytest.mark.parametrize("case", ["listed", "absent", "unlisted", "none"])
def test_cog_table_is_resource_cog2020_def(case: str, tmp_path: Path) -> None:
    from giap.setupdb.installed import cog_table, resource

    db_dir: Path | None = tmp_path
    if case == "listed":
        write_function_db(tmp_path, _xrefs())
    elif case == "absent":
        db = write_function_db(tmp_path, _xrefs())
        assert db.paths["cog2020_def"] is not None
        db.paths["cog2020_def"].unlink()
    elif case == "unlisted":
        write_function_db(tmp_path, _xrefs(), listed=FUNCTION_RESOURCES)
    else:
        db_dir = None
    assert cog_table(db_dir) == resource(db_dir, "cog2020_def")
    assert (cog_table(db_dir) is not None) == (case == "listed")


@settings(
    max_examples=30, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(written=st.sets(st.sampled_from(ITEMS)), listed=st.sets(st.sampled_from(ITEMS)))
def test_resource_is_listed_and_present(
    written: set[str], listed: set[str], tmp_path_factory: pytest.TempPathFactory
) -> None:
    """Property: ``resource(db, n)`` is a path iff ``n`` is listed and its file exists, and then
    it is that file."""
    from giap.setupdb.installed import resource

    root = tmp_path_factory.mktemp("prop")
    write_function_db(root, _xrefs(), items=written, listed=listed)
    for name in ITEMS:
        got = resource(root, name)
        if name in written and name in listed:
            assert got == root / PATHS[name]
            assert got.is_file()
        else:
            assert got is None, name
