"""MVP-04 gate 2 matches the gene-level detection warning tokens exactly
(docs/specs/MVP-04-ecoli-evidence-run.md gate 2; statistician diagnosis of run 20261023_a1,
local/stat/mvp04-gates-2-12.md).

Gate 2 means the gene-level (detection step) warnings ``decoy_unscored``, ``decoy_pairing_low``,
``q_floor_above_alpha`` and ``no_eligible_genes``. The novel_orfs step's per-class warning
``novel_orfs.class_q_floor_above_alpha`` (src/gippa/steps/novel.py, manifest entry
``novel_orfs:novel_orfs.class_q_floor_above_alpha``) is a different warning and must not fail the
gate; the scorer found it by substring. The same holds for a ``class_`` form of each other token.

An exact token still fails the gate wherever it appears: a manifest ``warnings`` entry
(``<step>:<event>``), a bare manifest flag, or an events.jsonl line (as the event name's last
dotted part or as a ``data.flag`` value).

Uses the all-pass fixture of tests/helpers/mvp04_summary_data.py unchanged.
"""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest
from helpers import mvp04_summary_data as S

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "dev" / "mvp04_summary.py"
TOKENS = ("decoy_unscored", "decoy_pairing_low", "q_floor_above_alpha", "no_eligible_genes")


def mod() -> ModuleType:
    spec = importlib.util.spec_from_file_location("mvp04_summary_gate2", SCRIPT)
    assert spec is not None and spec.loader is not None, SCRIPT
    m = importlib.util.module_from_spec(spec)
    sys.modules["mvp04_summary_gate2"] = m
    spec.loader.exec_module(m)
    return m


def gate2(d: S.Dirs) -> dict[str, Any]:
    res = mod().gate2(d.out, d.inputs)
    assert isinstance(res, dict) and "warnings_found" in res, res
    return res


@pytest.fixture
def d(tmp_path: Path) -> S.Dirs:
    return S.build(tmp_path)


def _manifest_warning(d: S.Dirs, text: str) -> None:
    S.edit_json(d.out / "run_manifest.json", lambda m: m["warnings"].append(text))


def _event(d: S.Dirs, record: dict[str, Any]) -> None:
    with (d.out / "logs" / "events.jsonl").open("a") as fh:
        fh.write(json.dumps(record) + "\n")


def _novel_class_event(token: str) -> dict[str, Any]:
    """The novel_orfs per-class warning as src/gippa/steps/novel.py logs it (``class_<token>``)."""
    return {
        "event": f"novel_orfs.class_{token}",
        "level": "WARNING",
        "step": "novel_orfs",
        "message": "novel_orfs: class intergenic (picked_decoy) has 6 targets; its q floor "
        "0.1667 > alpha 0.01, so none can be called",
        "data": {"regime": "picked_decoy", "class": "intergenic", "n_targets": 6,
                 "q_floor": 1 / 6, "alpha": 0.01},
    }  # fmt: skip


def test_fixture_passes_gate2(d: S.Dirs) -> None:
    res = gate2(d)
    assert res["warnings_found"] == [] and res["gate_ok"] is True


# ---------------------------------------------------------------- class_ forms do not count


@pytest.mark.parametrize("token", TOKENS)
@pytest.mark.parametrize(
    "form",
    ["novel_orfs:novel_orfs.class_{t}", "novel_orfs.class_{t}", "class_{t}"],
    ids=["step_colon_event", "event_only", "bare"],
)
def test_a_class_token_in_the_manifest_does_not_fail_gate2(
    d: S.Dirs, token: str, form: str
) -> None:
    _manifest_warning(d, form.format(t=token))
    res = gate2(d)
    assert res["warnings_found"] == [], res
    assert res["gate_ok"] is True and res["status"] == "pass", res


@pytest.mark.parametrize("token", TOKENS)
def test_a_class_token_event_does_not_fail_gate2(d: S.Dirs, token: str) -> None:
    _event(d, _novel_class_event(token))
    res = gate2(d)
    assert res["warnings_found"] == [], res
    assert res["gate_ok"] is True and res["status"] == "pass", res


def test_run_20261023_a1_novel_warnings_pass_gate2(d: S.Dirs) -> None:
    """The run's own manifest warnings and one class warning per picked class (diagnosis)."""
    _manifest_warning(d, "novel_orfs:novel_orfs.class_q_floor_above_alpha")
    _manifest_warning(d, "region_layer:region_layer.layer_missing")
    for cls, t in (("antisense", 14), ("cds_outframe", 11), ("intergenic", 6),
                   ("nterm_extension", 7)):  # fmt: skip
        ev = _novel_class_event("q_floor_above_alpha")
        ev["data"] |= {"class": cls, "n_targets": t, "q_floor": 1 / t}
        _event(d, ev)
    res = gate2(d)
    assert res["value"] == "decoy_picked"
    assert res["warnings_found"] == [] and res["status"] == "pass", res


# ---------------------------------------------------------------- exact tokens still fail


@pytest.mark.parametrize("token", TOKENS)
def test_an_exact_token_manifest_warning_fails_gate2(d: S.Dirs, token: str) -> None:
    _manifest_warning(d, f"gene_table:detection.{token}")
    res = gate2(d)
    assert res["warnings_found"] == [token] and res["status"] == "fail", res


@pytest.mark.parametrize("token", TOKENS)
def test_an_exact_token_manifest_flag_fails_gate2(d: S.Dirs, token: str) -> None:
    S.edit_json(d.out / "run_manifest.json", lambda m: m["flags"].append(token))
    res = gate2(d)
    assert res["warnings_found"] == [token] and res["status"] == "fail", res


@pytest.mark.parametrize("token", TOKENS)
def test_an_exact_token_event_name_fails_gate2(d: S.Dirs, token: str) -> None:
    _event(d, {"event": f"detection.{token}", "level": "WARNING", "step": "gene_table",
               "data": {}})  # fmt: skip
    res = gate2(d)
    assert res["warnings_found"] == [token] and res["status"] == "fail", res


@pytest.mark.parametrize("token", TOKENS)
def test_an_exact_token_event_flag_fails_gate2(d: S.Dirs, token: str) -> None:
    _event(d, {"event": "gene_table.detection", "level": "WARNING", "step": "gene_table",
               "data": {"flag": token}})  # fmt: skip
    res = gate2(d)
    assert res["warnings_found"] == [token] and res["status"] == "fail", res


@pytest.mark.parametrize("token", TOKENS)
def test_exact_token_found_next_to_its_class_form(d: S.Dirs, token: str) -> None:
    """Both present: the exact one is found (once), the class form adds nothing."""
    _manifest_warning(d, f"novel_orfs:novel_orfs.class_{token}")
    _event(d, _novel_class_event(token))
    _manifest_warning(d, f"gene_table:detection.{token}")
    res = gate2(d)
    assert res["warnings_found"] == [token] and res["status"] == "fail", res


def test_only_the_exact_tokens_are_reported(d: S.Dirs) -> None:
    """Class forms of every token plus the exact ``decoy_pairing_low``: only that one is found."""
    for t in TOKENS:
        _manifest_warning(d, f"novel_orfs:novel_orfs.class_{t}")
        _event(d, _novel_class_event(t))
    _event(d, {"event": "gene_table.detection", "level": "WARNING", "step": "gene_table",
               "data": {"flag": "decoy_pairing_low"}})  # fmt: skip
    res = gate2(d)
    assert res["warnings_found"] == ["decoy_pairing_low"] and res["status"] == "fail", res


# --------------------------------------------------------------- #326 review S1 (surviving mutants)


@pytest.mark.parametrize("token", TOKENS)
@pytest.mark.parametrize("form", ["gene_table:{t}", "/{t}", '"{t}"'])
def test_a_colon_slash_or_quote_bounded_token_fails_gate2(d: S.Dirs, token: str, form: str) -> None:
    _manifest_warning(d, form.format(t=token))
    res = gate2(d)
    assert res["warnings_found"] == [token] and res["status"] == "fail", res


@pytest.mark.parametrize("token", TOKENS)
@pytest.mark.parametrize("form", ["{t}_x", "novel_{t}", "gene_{t}", "{t}2", "x{t}"])
def test_a_token_inside_a_longer_word_does_not_fail_gate2(d: S.Dirs, token: str, form: str) -> None:
    _manifest_warning(d, f"novel_orfs:{form.format(t=token)}")
    res = gate2(d)
    assert res["warnings_found"] == [] and res["status"] == "pass", res
