"""``READINESS.md`` rendering, ``overall`` and the threshold constants (issue M11-03, acceptance
criteria 3-4; edge case "Row-level data"; R3, R11).

Module interface (``scripts/release/readiness.py`` loaded by path): ``render(doc) -> str`` is the
only producer of ``READINESS.md`` (a pure function of ``readiness.json``); ``overall(statuses)``
aggregates row statuses (``fail`` > ``pending`` > ``pass``; ``deferred`` and ``waived`` count as
pass); ``TE1_MIN``, ``MAPPING_FDR_MAX``, ``RUNTIME_MIN_S``, ``RUNTIME_MAX_S`` are the thresholds.
"""

from __future__ import annotations

import ast
import json
import re
from pathlib import Path

import pytest
from helpers.release_readiness import REPO, SCRIPT, Fixture, load_module, metrics_doc
from hypothesis import given
from hypothesis import strategies as st

PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")  # the M10-05 guard rule


@pytest.fixture
def computed(tmp_path: Path) -> Fixture:
    fx = Fixture.build(tmp_path)
    fx.compute()
    return fx


def test_readiness_md_is_render_of_json(computed: Fixture) -> None:
    mod = load_module()
    doc = json.loads((computed.out / "readiness.json").read_text(encoding="utf-8"))
    md = (computed.out / "READINESS.md").read_text(encoding="utf-8")
    first, second = mod.render(doc), mod.render(json.loads(json.dumps(doc)))
    assert first == second == md


def test_render_shows_every_row_and_overall(computed: Fixture) -> None:
    md = (computed.out / "READINESS.md").read_text(encoding="utf-8")
    doc = json.loads((computed.out / "readiness.json").read_text(encoding="utf-8"))
    for r in doc["rows"]:
        assert re.search(rf"\|\s*{r['id']}\s*\|", md), r["id"]
    assert doc["commit"] in md and "1.0.0" in md and doc["overall"] in md


def test_render_depends_on_json_only(computed: Fixture) -> None:
    mod = load_module()
    doc = json.loads((computed.out / "readiness.json").read_text(encoding="utf-8"))
    for r in doc["rows"]:
        if r["id"] == "R1":
            r["value"] = "12.34 %"
    assert "12.34 %" in mod.render(doc)


def test_no_peptide_like_strings(tmp_path: Path) -> None:
    fx = Fixture.build(tmp_path)
    doc = metrics_doc(share=0.9949)
    doc["te1"]["misses"] = [{"peptide": "PEPTIDEKLMNQR", "n": 1}]  # row-level, must not leak
    fx.write_json("validation/pilot/metrics/metrics.json", doc)
    fx.compute()
    md = (fx.out / "READINESS.md").read_text(encoding="utf-8")
    assert not PEPTIDE_LIKE.search(md), PEPTIDE_LIKE.findall(md)


@given(st.lists(st.sampled_from(["pass", "fail", "pending", "deferred", "waived"]), max_size=20))
def test_property_overall(statuses: list[str]) -> None:
    got = load_module().overall(statuses)
    if "fail" in statuses:
        assert got == "fail"
    elif "pending" in statuses:
        assert got == "pending"
    else:
        assert got == "pass"


# ------------------------------------------------------------------ thresholds (criterion 4, R3)


def test_threshold_constants() -> None:
    mod = load_module()
    assert mod.TE1_MIN == 0.995
    assert mod.MAPPING_FDR_MAX == 0.01
    assert mod.RUNTIME_MIN_S == 960
    assert mod.RUNTIME_MAX_S == 3300


def test_thresholds_are_literal_module_constants() -> None:
    tree = ast.parse(SCRIPT.read_text(encoding="utf-8"))
    lit: dict[str, object] = {}
    for node in tree.body:
        if isinstance(node, ast.Assign | ast.AnnAssign):
            targets = node.targets if isinstance(node, ast.Assign) else [node.target]
            for t in targets:
                if isinstance(t, ast.Name) and isinstance(node.value, ast.Constant):
                    lit[t.id] = node.value.value
    assert lit.get("TE1_MIN") == 0.995 and lit.get("MAPPING_FDR_MAX") == 0.01
    assert lit.get("RUNTIME_MIN_S") == 960 and lit.get("RUNTIME_MAX_S") == 3300


def test_thresholds_not_read_from_config() -> None:
    src = SCRIPT.read_text(encoding="utf-8")
    for cfg in ("gates.yaml", "tolerances.yaml", "calibration.yaml", "defaults.yaml"):
        assert cfg not in src, f"readiness.py names {cfg}: thresholds must be constants (R3)"


def test_script_exists_in_scripts_release() -> None:
    assert SCRIPT == REPO / "scripts" / "release" / "readiness.py"
    assert SCRIPT.is_file()
