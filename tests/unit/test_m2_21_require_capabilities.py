"""``--require-capabilities`` with the nightly gate (M2-21 example 2; GIAP_SPEC.md §3.6).

"same with ``local/db/genomad_db`` renamed away → the job fails (not skips): ``required capability
needs_genomad_db is missing``". The capability plugin (``tests/conftest.py``) probes the database
directory live (``$GIAP_DB_DIR``), not ``.giap/capabilities.json``; the test builds a tmp DB
directory holding every nightly database layout (and a tmp ``capabilities.json`` claiming
everything present, which must not mask the live probe) and runs a pytester session with the full
``farmshare_nightly`` list of ``ci/gates.yaml``.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
CONFTEST = (REPO / "tests" / "conftest.py").read_text()
GATES = yaml.safe_load((REPO / "ci" / "gates.yaml").read_text())["required_capabilities"]
DB_CAPS = ("needs_bakta_db", "needs_bakta_full", "needs_genomad_db", "needs_promotech_models")

pytest_plugins = ["pytester"]

TESTS = """
import pytest

@pytest.mark.needs_genomad_db
def test_genomad():
    pass

@pytest.mark.needs_bakta_full
def test_bakta_full():
    pass

@pytest.mark.needs_promotech_models
def test_promotech():
    pass
"""


def _db(root: Path, genomad: bool) -> Path:
    db = root / "db"
    (db / "bakta" / "db").mkdir(parents=True)
    (db / "bakta" / "db-light").mkdir(parents=True)
    (db / "promotech" / "models").mkdir(parents=True)
    (db / "promotech" / "models" / "RF-HOT.model").write_text("x")
    if genomad:
        (db / "genomad_db").mkdir()
    return db


def _session(pytester: pytest.Pytester, monkeypatch: pytest.MonkeyPatch, genomad: bool,
             *args: str) -> pytest.RunResult:  # fmt: skip
    db = _db(pytester.path, genomad)
    caps = pytester.path / ".giap" / "capabilities.json"
    caps.parent.mkdir()
    caps.write_text(json.dumps({c: True for c in GATES["farmshare_nightly"]}))
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    monkeypatch.setenv("GIAP_ROOT", str(pytester.path))
    pytester.makeconftest(CONFTEST)
    pytester.makeini("[pytest]\nmarkers =\n" + "".join(f"    {c}: x\n" for c in DB_CAPS))
    pytester.makepyfile(TESTS)
    return pytester.runpytest_subprocess("-p", "no:cacheprovider", "-rs", *args)


def test_nightly_gate_list_is_known_to_the_plugin() -> None:
    import importlib.util

    spec = importlib.util.spec_from_file_location(
        "m2_21_root_conftest", REPO / "tests" / "conftest.py"
    )
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    CHECKS = mod.CHECKS
    assert set(GATES["farmshare_nightly"]) <= set(CHECKS)
    assert set(GATES["farmshare_weekly"]) <= set(CHECKS)
    assert set(DB_CAPS) <= set(GATES["farmshare_nightly"])


def test_all_databases_present_passes(
    pytester: pytest.Pytester, monkeypatch: pytest.MonkeyPatch
) -> None:
    res = _session(pytester, monkeypatch, True, "--require-capabilities=" + ",".join(DB_CAPS))
    res.assert_outcomes(passed=3)


def test_example_2_genomad_db_renamed_away_fails_not_skips(
    pytester: pytest.Pytester, monkeypatch: pytest.MonkeyPatch
) -> None:
    res = _session(pytester, monkeypatch, False, "--require-capabilities=" + ",".join(DB_CAPS))
    out = res.stdout.str()
    assert "required capability needs_genomad_db is missing" in out, out
    outcomes = res.parseoutcomes()
    assert outcomes.get("skipped", 0) == 0, out
    assert outcomes.get("passed", 0) == 2, out
    assert outcomes.get("failed", 0) + outcomes.get("errors", outcomes.get("error", 0)) == 1, out
    assert res.ret != 0


def test_without_the_option_a_missing_database_is_a_skip(
    pytester: pytest.Pytester, monkeypatch: pytest.MonkeyPatch
) -> None:
    res = _session(pytester, monkeypatch, False)
    res.assert_outcomes(passed=2, skipped=1)
    res.stdout.fnmatch_lines(["*missing capability: needs_genomad_db*"])
