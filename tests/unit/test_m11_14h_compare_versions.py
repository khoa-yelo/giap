"""M11-14h: ``giap compare`` version, digest, database and operon-source warnings.

docs/specs/M11-14-container-install.md item 1 option (c) "Compare (14h) warns
``compare.operon_source: run A uniop, run B distance_model`` when the sources differ", item 1b
(a local UniOP build: warn on ``source_sha256``, the key or the runtime, never on its
``container_digest``), §4 "Compatibility check (14h)" (the table of effects), examples 10 and 11,
acceptance criterion 12 ("two fake run directories whose manifests differ only in the listed
fields, including the ``operon_source`` warning"); ADR-0060 §6.

Interface pinned here (test-author reading, beyond the spec text):

- ``giap.compare.versions.compare_versions(manifests)`` is pure. ``manifests`` maps each run
  label (in run order) to its parsed ``run_manifest.json`` (``{}`` when the run has none). It
  returns the warning lines, each starting with its event name (``compare.cross_version: ...``),
  one line per differing tool. It reads every key with ``.get``: a manifest written before 14d
  or 14h (no ``operon_source``, no ``output_schema_version``, no ``tools``) never raises, and a
  key absent from both runs never warns.
- The CLI prints every line on stderr (``CliRunner`` output), writes each as a ``WARNING``
  record of ``CMP/logs/events.jsonl`` (``step: compare``, ``event`` = the line's prefix, the full
  line somewhere in the record) and shows it in ``CMP/report/compare.html``; the token
  ``compare_versions`` appears in both the page and the event log.
- A differing ``output_schema_version`` major is checked before any table is written: exit 1
  with ``compare: output schema 1.x and 2.x differ; rerun one with the same major``; with
  ``--allow-cross-major`` the run goes on (exit 0) and ``output schema 1.x and 2.x differ`` stays
  as a warning.

Runs are small ``OUT/`` trees from ``tests/helpers/compare_data.py`` with a synthetic
``run_manifest.json``. Never touches the real ``$GIAP_DB_DIR`` or the network.
"""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from click.testing import CliRunner, Result
from helpers.compare_data import SHA_A, gene, write_run
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.exitcodes import ExitCode

FAST = ["--set", "stats.bootstrap=50"]
N = 8
TABLES = ("compare_detection", "compare_abundance", "compare_stability")

OSTIR_A = "sha256:3f9c0a1b2d4e" + "0" * 52
OSTIR_B = "sha256:77a0c1d2e3f4" + "1" * 52


@pytest.fixture(autouse=True)
def _no_db(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    (tmp_path / "db_empty").mkdir()
    monkeypatch.setenv("GIAP_DB_DIR", str(tmp_path / "db_empty"))
    monkeypatch.setenv("GIAP_ROOT", str(tmp_path / "db_empty"))


# ---------------------------------------------------------------- fixtures


def _genes(shift: float) -> list[dict[str, object]]:
    vals = np.sin(np.arange(N) * 1.7) * 3 + np.arange(N) * 0.4
    return [
        gene(f"g{i:02d}", start=1000 * i, log2=float(v + shift * (i % 3)), tu=f"TU_{i // 2}")
        for i, v in enumerate(vals)
    ]


def tool(name: str, version: str, digest: str | None = None, **extra: Any) -> dict[str, Any]:
    return {
        "name": name,
        "version": version,
        "container_digest": digest,
        "conda_lock_hash": None,
        "processes": [name.upper()],
        **extra,
    }


def db(name: str, version: str | None) -> dict[str, Any]:
    return {
        "name": name,
        "version": version,
        "checksum": None,
        "retrieved": None,
        "size_bytes": None,
    }


def manifest(**over: Any) -> dict[str, Any]:
    """The fields compare reads, equal between runs unless overridden; ``None`` drops a key."""
    doc: dict[str, Any] = {
        "schema_version": "1.0",
        "output_schema_version": "1.0",
        "giap_version": "1.0.0",
        "inputs": [
            {
                "role": "genome",
                "file": "genome.fna",
                "path": "/x/genome.fna",
                "sha256": SHA_A,
                "size": 1000,
            }
        ],
        "tools": [tool("ostir", "1.1.3", OSTIR_A), tool("bakta", "1.9.4")],
        "databases": [db("bakta_db", "6.0")],
        "operon_source": "uniop",
        "operon_source_reason": None,
    }
    for key, value in over.items():
        if value is None:
            doc.pop(key, None)
        else:
            doc[key] = value
    return doc


def _run(path: Path, doc: Mapping[str, Any] | None, shift: float) -> Path:
    out = write_run(path, _genes(shift), manifest_sha=SHA_A if doc is not None else None)
    if doc is not None:
        (out / "run_manifest.json").write_text(json.dumps(doc, indent=2) + "\n", encoding="utf-8")
    return out


def _pair(tmp: Path, a: Mapping[str, Any] | None, b: Mapping[str, Any] | None) -> tuple[Path, Path]:
    return _run(tmp / "runA", a, 0.0), _run(tmp / "runB", b, 0.5)


def _compare(tmp: Path, a: Path, b: Path, *extra: str) -> Result:
    from giap.cli import cli

    args = [str(a), str(b), "--outdir", str(tmp / "CMP"), "--labels", "A,B", *FAST, *extra]
    res = CliRunner().invoke(cli, ["compare", *args])
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    return res


def _events(tmp: Path) -> list[dict[str, Any]]:
    path = tmp / "CMP" / "logs" / "events.jsonl"
    assert path.is_file(), f"no {path}"
    return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]


def _html(tmp: Path) -> str:
    return (tmp / "CMP" / "report" / "compare.html").read_text(encoding="utf-8")


def _warning_lines(text: str) -> list[str]:
    return [m.group(0) for m in re.finditer(r"compare\.[a-z_]+: [^\n]*", text)]


# ---------------------------------------------------------------- the operon-source warning


def test_operon_source_differs_warns_in_output_events_and_report(tmp_path: Path) -> None:
    """Item 1 (c): the exact warning line of the spec; compare still writes its tables."""
    a, b = _pair(
        tmp_path,
        manifest(),
        manifest(operon_source="distance_model", operon_source_reason="UniOP not installed"),
    )
    res = _compare(tmp_path, a, b)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    line = "compare.operon_source: run A uniop, run B distance_model"
    assert line in res.output
    assert all((tmp_path / "CMP" / "tables" / f"{t}.tsv").is_file() for t in TABLES)
    hits = [e for e in _events(tmp_path) if e.get("event") == "compare.operon_source"]
    assert len(hits) == 1, _events(tmp_path)
    assert hits[0]["level"] == "WARNING" and hits[0]["step"] == "compare"
    assert line in json.dumps(hits[0])
    assert line in _html(tmp_path)


def test_same_operon_source_does_not_warn(tmp_path: Path) -> None:
    a, b = _pair(tmp_path, manifest(), manifest())
    res = _compare(tmp_path, a, b)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert "compare.operon_source" not in res.output
    assert not any(e.get("event") == "compare.operon_source" for e in _events(tmp_path))


@pytest.mark.parametrize(
    ("src_a", "src_b"),
    [
        ("uniop", "distance_model"),
        ("distance_model", "uniop"),
        ("uniop", "precomputed"),
        ("precomputed", "distance_model"),
    ],
)
def test_operon_source_line_names_each_run_in_order(src_a: str, src_b: str) -> None:
    from giap.compare.versions import compare_versions

    out = compare_versions({"A": manifest(operon_source=src_a), "B": manifest(operon_source=src_b)})
    assert f"compare.operon_source: run A {src_a}, run B {src_b}" in out


def test_operon_source_reason_alone_does_not_warn() -> None:
    """The table lists ``operon_source``; the reason is free text."""
    from giap.compare.versions import compare_versions

    out = compare_versions(
        {
            "A": manifest(
                operon_source="distance_model", operon_source_reason="UniOP not installed"
            ),
            "B": manifest(
                operon_source="distance_model",
                operon_source_reason="UniOP disabled (annotation.uniop.enabled=false)",
            ),
        }
    )
    assert out == []


# ---------------------------------------------------------------- example 10 and the table rows


def test_example_10_cross_version_and_tool_digest(tmp_path: Path) -> None:
    """Example 10: a 1.0.0 run and a 1.1.0 run whose OSTIR digests differ: exit 0, two warnings."""
    a, b = _pair(
        tmp_path,
        manifest(),
        manifest(
            giap_version="1.1.0", tools=[tool("ostir", "1.1.4", OSTIR_B), tool("bakta", "1.9.4")]
        ),
    )
    res = _compare(tmp_path, a, b)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert "compare.cross_version: runs from 1.0.0 and 1.1.0" in res.output
    assert (
        "compare.tool_versions: ostir 1.1.3 sha256:3f9c0a1b2d4e vs 1.1.4 sha256:77a0c1d2e3f4"
        in res.output
    )
    events = {e.get("event") for e in _events(tmp_path) if e.get("level") == "WARNING"}
    assert {"compare.cross_version", "compare.tool_versions"} <= events
    html = _html(tmp_path)
    assert "compare.cross_version: runs from 1.0.0 and 1.1.0" in html
    assert "compare_versions" in html
    assert any("compare_versions" in json.dumps(e) for e in _events(tmp_path))


def test_package_patch_only_differs_gives_nothing() -> None:
    from giap.compare.versions import compare_versions

    assert (
        compare_versions({"A": manifest(giap_version="1.0.2"), "B": manifest(giap_version="1.0.7")})
        == []
    )


def test_package_minor_differs_gives_the_spec_line() -> None:
    from giap.compare.versions import compare_versions

    out = compare_versions(
        {"A": manifest(giap_version="1.0.2"), "B": manifest(giap_version="1.1.0")}
    )
    assert out == ["compare.cross_version: runs from 1.0.2 and 1.1.0"]


def test_same_tool_version_different_digest_warns() -> None:
    """Any tool ``version`` or ``container_digest`` differs: the digest alone is enough."""
    from giap.compare.versions import compare_versions

    out = compare_versions(
        {
            "A": manifest(),
            "B": manifest(tools=[tool("ostir", "1.1.3", OSTIR_B), tool("bakta", "1.9.4")]),
        }
    )
    assert out == [
        "compare.tool_versions: ostir 1.1.3 sha256:3f9c0a1b2d4e vs 1.1.3 sha256:77a0c1d2e3f4"
    ]


def test_one_line_per_differing_tool() -> None:
    from giap.compare.versions import compare_versions

    out = compare_versions(
        {
            "A": manifest(),
            "B": manifest(tools=[tool("ostir", "1.1.4", OSTIR_B), tool("bakta", "1.10.0")]),
        }
    )
    lines = [x for x in out if x.startswith("compare.tool_versions: ")]
    assert len(lines) == 2, out
    assert sum(" ostir 1.1.3 " in f" {x.split(': ', 1)[1]} " for x in lines) == 1
    assert sum(x.split(": ", 1)[1].startswith("bakta 1.9.4") for x in lines) == 1
    assert all(" vs " in x for x in lines)


def test_database_version_differs() -> None:
    from giap.compare.versions import compare_versions

    out = compare_versions({"A": manifest(), "B": manifest(databases=[db("bakta_db", "6.1")])})
    assert out == ["compare.db_versions: bakta_db 6.0 vs 6.1"]


def test_manifests_equal_in_every_listed_field_give_no_warning(tmp_path: Path) -> None:
    a, b = _pair(tmp_path, manifest(), manifest())
    res = _compare(tmp_path, a, b)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert _warning_lines(res.output) == []


# ---------------------------------------------------------------- item 1b: local UniOP builds


def _uniop(**over: Any) -> dict[str, Any]:
    base = {"local_build": True, "runtime": "docker", "source_sha256": "c" * 64}
    return tool("uniop", "git:3ed14dfc", "sha256:" + "d" * 64, **{**base, **over})


def _with_uniop(**over: Any) -> dict[str, Any]:
    return manifest(tools=[tool("ostir", "1.1.3", OSTIR_A), _uniop(**over)])


def test_local_uniop_rebuild_with_another_digest_does_not_warn() -> None:
    """Item 1b: a rebuild of the same key on another host is not a tool change."""
    from giap.compare.versions import compare_versions

    out = compare_versions(
        {"A": _with_uniop(), "B": _with_uniop(container_digest="sha256:" + "e" * 64)}
    )
    assert out == []


@pytest.mark.parametrize("over", [{"source_sha256": "f" * 64}, {"runtime": "apptainer"}])
def test_local_uniop_source_or_runtime_change_warns(over: dict[str, Any]) -> None:
    from giap.compare.versions import compare_versions

    out = compare_versions({"A": _with_uniop(), "B": _with_uniop(**over)})
    lines = [x for x in out if x.startswith("compare.tool_versions: uniop")]
    assert len(lines) == 1, out


# ---------------------------------------------------------------- example 11: schema majors


def test_example_11_cross_major_exits_1_and_writes_no_table(tmp_path: Path) -> None:
    a, b = _pair(
        tmp_path, manifest(output_schema_version="1.3"), manifest(output_schema_version="2.0")
    )
    res = _compare(tmp_path, a, b)
    assert res.exit_code == ExitCode.INVALID_INPUT, res.output
    assert "compare: output schema 1.x and 2.x differ; rerun one with the same major" in res.output
    assert not any((tmp_path / "CMP" / "tables" / f"{t}.tsv").exists() for t in TABLES)


def test_example_11_allow_cross_major_keeps_the_warning(tmp_path: Path) -> None:
    a, b = _pair(
        tmp_path, manifest(output_schema_version="1.3"), manifest(output_schema_version="2.0")
    )
    res = _compare(tmp_path, a, b, "--allow-cross-major")
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert "output schema 1.x and 2.x differ" in res.output
    assert all((tmp_path / "CMP" / "tables" / f"{t}.tsv").is_file() for t in TABLES)


def test_schema_minor_only_differs_is_not_an_error(tmp_path: Path) -> None:
    a, b = _pair(
        tmp_path, manifest(output_schema_version="1.0"), manifest(output_schema_version="1.3")
    )
    res = _compare(tmp_path, a, b)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert "output schema" not in res.output


# ---------------------------------------------------------------- manifests from before 14d/14h


OLD_KEYS = ("operon_source", "operon_source_reason", "output_schema_version")


def test_manifests_without_the_new_fields_still_compare(tmp_path: Path) -> None:
    """A pre-14d manifest (no ``operon_source*``, no ``output_schema_version``) is read with
    ``.get``: exit 0, tables written, no operon-source warning."""
    old = manifest(**dict.fromkeys(OLD_KEYS))
    assert not set(OLD_KEYS) & set(old)
    a, b = _pair(tmp_path, old, dict(old))
    res = _compare(tmp_path, a, b)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert _warning_lines(res.output) == []
    assert all((tmp_path / "CMP" / "tables" / f"{t}.tsv").is_file() for t in TABLES)


def test_minimal_manifests_and_preflight_only_runs_still_compare(tmp_path: Path) -> None:
    """The M9-22 fixtures: a manifest of ``schema_version`` + ``inputs`` only, and a run with
    no manifest at all (genome from ``preflight.json``)."""
    a = write_run(tmp_path / "runA", _genes(0.0), manifest_sha=SHA_A)
    b = write_run(tmp_path / "runB", _genes(0.5))
    assert not (b / "run_manifest.json").exists()
    res = _compare(tmp_path, a, b)
    assert res.exit_code == ExitCode.SUCCESS, res.output
    assert "Traceback" not in res.output


def test_old_and_new_manifest_mix_does_not_fail(tmp_path: Path) -> None:
    """One run from before 14h (no ``output_schema_version``), one with ``2.0``: no major can be
    compared, so no exit 1."""
    a, b = _pair(
        tmp_path, manifest(**dict.fromkeys(OLD_KEYS)), manifest(output_schema_version="2.0")
    )
    res = _compare(tmp_path, a, b)
    assert res.exit_code == ExitCode.SUCCESS, res.output


# ---------------------------------------------------------------- properties


_VER = st.tuples(st.integers(0, 3), st.integers(0, 3), st.integers(0, 3))
_KEYS = (
    "giap_version",
    "tools",
    "databases",
    "operon_source",
    "operon_source_reason",
    "output_schema_version",
    "inputs",
    "schema_version",
)


@settings(max_examples=60, deadline=None)
@given(va=_VER, vb=_VER)
def test_cross_version_iff_major_or_minor_differs(
    va: tuple[int, int, int], vb: tuple[int, int, int]
) -> None:
    from giap.compare.versions import compare_versions

    a, b = ".".join(map(str, va)), ".".join(map(str, vb))
    out = compare_versions({"A": manifest(giap_version=a), "B": manifest(giap_version=b)})
    cross = [x for x in out if x.startswith("compare.cross_version: ")]
    if va[:2] == vb[:2]:
        assert cross == []
    else:
        assert cross == [f"compare.cross_version: runs from {a} and {b}"]


@settings(max_examples=60, deadline=None)
@given(drop_a=st.sets(st.sampled_from(_KEYS)), drop_b=st.sets(st.sampled_from(_KEYS)))
def test_any_missing_key_never_raises_and_warns_only_with_both_values(
    drop_a: set[str], drop_b: set[str]
) -> None:
    """Old manifests: any subset of keys may be absent. Every value present in both runs is
    equal here, so a key missing from one run never produces a warning."""
    from giap.compare.versions import compare_versions

    a = manifest(**dict.fromkeys(drop_a))
    b = manifest(**dict.fromkeys(drop_b))
    out = compare_versions({"A": a, "B": b})
    assert out == []
    assert compare_versions({"A": a, "B": dict(a)}) == []
    assert compare_versions({"A": {}, "B": {}}) == []
