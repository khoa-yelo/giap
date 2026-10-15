"""`db_manifest.json` (issue M2-08; GIAP_SPEC.md §7.9.1, §6.7): schema, determinism, merging,
atomic writes, and the readers (`environment.tools` in preflight, `scripts/probe_env.py`).

Ambiguities resolved here:
- The manifest is ``{"schema_version": "1.0", "items": [DbRecord as a JSON object, ...]}`` sorted
  by name; a record's JSON object has exactly the ``DbRecord`` fields.
- ``write_manifest`` replaces a record with the same name and keeps the others.
- Atomic means write-to-temp then ``os.replace``/``os.rename``: if that rename fails the old
  manifest is intact.
- Preflight maps manifest names to its requirement names (``bakta_light``/``bakta_full`` ->
  ``bakta_db``, ``genomad`` -> ``genomad_db``, ``antismash`` -> ``antismash_db``,
  ``promotech_models`` unchanged) and counts only ``status: "ok"`` records as found.
"""

from __future__ import annotations

import dataclasses
import importlib.util
import json
import os
from pathlib import Path
from typing import Any

import jsonschema
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st
from setup_guard import setup_guard  # noqa: F401  (autouse: no network, no tools; decision 3)

ROOT = Path(__file__).resolve().parents[2]
SCHEMA = ROOT / "src" / "giap" / "schemas" / "db_manifest.schema.json"
FIELDS = {
    "name",
    "version",
    "url",
    "checksum",
    "size_bytes",
    "retrieved",
    "path",
    "status",
    "message",
}


def rec(name: str, status: str = "ok", version: str = "6.0", **kw: Any) -> Any:
    from giap.setupdb import DbRecord

    base: dict[str, Any] = {
        "name": name,
        "version": version,
        "url": f"https://zenodo.org/records/1/{name}.tar.gz",
        "checksum": "md5:" + "0" * 32,
        "size_bytes": 123,
        "retrieved": "2026-10-15",
        "path": name,
        "status": status,
    }
    base.update(kw)
    return DbRecord(**base)


def read(db: Path) -> dict[str, Any]:
    data = json.loads((db / "db_manifest.json").read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def schema() -> dict[str, Any]:
    data = json.loads(SCHEMA.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


# ---------------------------------------------------------------- schema


def test_schema_is_a_valid_draft_2020_12_schema() -> None:
    s = schema()
    assert s["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    jsonschema.Draft202012Validator.check_schema(s)


def test_written_manifest_validates_and_has_the_documented_shape(tmp_path: Path) -> None:
    from giap.setupdb import write_manifest

    write_manifest(tmp_path, [rec("genomad", version="1.9"), rec("bakta_light")])
    data = read(tmp_path)
    jsonschema.Draft202012Validator(schema()).validate(data)
    assert set(data) == {"schema_version", "items"} and data["schema_version"] == "1.0"
    assert [i["name"] for i in data["items"]] == ["bakta_light", "genomad"]
    for item in data["items"]:
        assert set(item) == FIELDS


@pytest.mark.parametrize(
    "bad",
    [
        {"schema_version": "1.0"},
        {"schema_version": "1.0", "items": [{"name": "x"}]},
        {"schema_version": "1.0", "items": [dict.fromkeys(FIELDS, "x") | {"size_bytes": 1}]},
    ],
)
def test_schema_rejects_malformed_manifests(bad: dict[str, Any]) -> None:
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema()).validate(bad)


def test_schema_status_is_an_enum() -> None:
    good = dataclasses.asdict(rec("bakta_light"))
    v = jsonschema.Draft202012Validator(schema())
    v.validate({"schema_version": "1.0", "items": [good]})
    with pytest.raises(jsonschema.ValidationError):
        v.validate({"schema_version": "1.0", "items": [good | {"status": "maybe"}]})


# ---------------------------------------------------------------- merge, determinism, atomicity


def test_write_manifest_merges_by_name(tmp_path: Path) -> None:
    from giap.setupdb import write_manifest

    write_manifest(tmp_path, [rec("bakta_light"), rec("genomad", version="1.9")])
    write_manifest(tmp_path, [rec("genomad", status="mismatch", version="1.9")])
    got = {i["name"]: i for i in read(tmp_path)["items"]}
    assert set(got) == {"bakta_light", "genomad"}
    assert got["genomad"]["status"] == "mismatch" and got["bakta_light"]["status"] == "ok"


def test_write_manifest_leaves_no_temporary_files(tmp_path: Path) -> None:
    from giap.setupdb import write_manifest

    write_manifest(tmp_path, [rec("bakta_light")])
    write_manifest(tmp_path, [rec("genomad")])
    assert [p.name for p in tmp_path.iterdir() if p.name != ".giap_setup.lock"] == [
        "db_manifest.json"
    ]


def test_write_manifest_is_atomic(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from giap.setupdb import write_manifest

    write_manifest(tmp_path, [rec("bakta_light")])
    before = (tmp_path / "db_manifest.json").read_bytes()

    def boom(*_a: Any, **_k: Any) -> None:
        raise OSError("disk full")

    monkeypatch.setattr(os, "replace", boom)
    monkeypatch.setattr(os, "rename", boom)
    with pytest.raises(OSError):
        write_manifest(tmp_path, [rec("genomad")])
    assert (tmp_path / "db_manifest.json").read_bytes() == before


_names = st.lists(
    st.sampled_from(["bakta_light", "bakta_full", "genomad", "antismash", "promotech_models"]),
    unique=True,
    min_size=1,
)


@settings(max_examples=40, deadline=None)
@given(names=_names, data=st.data())
def test_manifest_bytes_do_not_depend_on_record_order(
    tmp_path_factory: pytest.TempPathFactory, names: list[str], data: st.DataObject
) -> None:
    from giap.setupdb import write_manifest

    records = [
        rec(n, status=data.draw(st.sampled_from(["ok", "absent", "mismatch"]))) for n in names
    ]
    shuffled = data.draw(st.permutations(records))
    a, b = tmp_path_factory.mktemp("a"), tmp_path_factory.mktemp("b")
    write_manifest(a, records)
    write_manifest(b, shuffled)
    assert (a / "db_manifest.json").read_bytes() == (b / "db_manifest.json").read_bytes()
    items = read(a)["items"]
    assert [i["name"] for i in items] == sorted(names)
    assert items == [dataclasses.asdict(r) for r in sorted(records, key=lambda r: r.name)]


@settings(max_examples=40, deadline=None)
@given(first=_names, second=_names)
def test_merge_is_union_with_later_records_winning(
    tmp_path_factory: pytest.TempPathFactory, first: list[str], second: list[str]
) -> None:
    from giap.setupdb import write_manifest

    d = tmp_path_factory.mktemp("m")
    write_manifest(d, [rec(n, message="first") for n in first])
    write_manifest(d, [rec(n, message="second") for n in second])
    got = {i["name"]: i["message"] for i in read(d)["items"]}
    assert set(got) == set(first) | set(second)
    assert all(got[n] == "second" for n in second)


# ---------------------------------------------------------------- readers


def _setup_manifest(db: Path) -> None:
    from giap.setupdb import write_manifest

    write_manifest(
        db,
        [
            rec("bakta_full", path="bakta/db"),
            rec("genomad", version="1.9"),
            rec("promotech_models", version="RF-HOT"),
            rec("antismash", status="mismatch", version="8.0"),
        ],
    )


def test_preflight_reads_the_setup_manifest(tmp_path: Path) -> None:
    from giap.preflight.environment import read_manifest

    _setup_manifest(tmp_path)
    found = read_manifest(str(tmp_path))
    assert found is not None
    assert found.get("bakta_db") == "6.0"
    assert found.get("genomad_db") == "1.9"
    assert found.get("promotech_models") is not None
    assert found.get("antismash_db") is None, "a mismatched item is not found"


def test_environment_tools_sees_setup_databases(tmp_path: Path) -> None:
    from giap.preflight.environment import (
        EnvProbe,
        environment_checks,
        read_manifest,
        requirements,
    )

    _setup_manifest(tmp_path)
    found = read_manifest(str(tmp_path))
    assert found is not None
    reqs = [r for r in requirements(["mobile", "bgc"], "targeted") if r.kind in ("db", "model")]
    tools = environment_checks(reqs, EnvProbe(found, 0, 0), 0, 0)[0]
    assert tools.level == "error"
    value = tools.value
    assert isinstance(value, dict)
    assert value.get("missing") == ["antismash_db"], value
    assert "wrong_version" not in value


def test_probe_env_reads_the_manifest(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    spec = importlib.util.spec_from_file_location(
        "probe_env_m208", ROOT / "scripts" / "probe_env.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path))
    without = json.dumps(mod.databases())
    _setup_manifest(tmp_path)
    with_manifest = json.dumps(mod.databases())
    assert "mismatch" not in without and "1.9" not in without
    assert "mismatch" in with_manifest and "1.9" in with_manifest, with_manifest
