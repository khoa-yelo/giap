"""The ``novel_orfs`` manifest block and the class-FDR warnings
(docs/specs/M7-11b-class-fdr-wiring.md rules 9-10, Q4, Q7; degenerate rows "resume duplicates
a class event" and "--strict with a floor or decoys_not_makedb warning").

* ``novel_orfs = {"classes": [...]}``: the de-duplicated ``novel_orfs.class_fdr`` event data
  (``_event_data``), sorted by ``(regime, class)``, NaN -> null; ``None`` for a run without
  the step (read here as: no such event). ``run_manifest.schema.json`` declares the block.
* The two WARNING events reach ``warnings`` as ``novel_orfs:<event>``; they are no run-flag
  tokens (Q4), so ``flags`` does not gain them.
"""

from __future__ import annotations

import json
import math
from typing import Any

from helpers.manifest_data import build, event, schema, validate

CLASS_EVENT = "novel_orfs.class_fdr"
FLOOR_EVENT = "novel_orfs.class_q_floor_above_alpha"
MAKEDB_EVENT = "novel_orfs.decoys_not_makedb"


def _picked(cls: str = "intergenic", n: int = 150, **over: Any) -> dict[str, Any]:
    data = {"regime": "picked_decoy", "class": cls, "estimated": True, "reason": None,
            "n_targets": n, "n_target_wins": n, "n_decoys": 30, "paired_share": 1.0,
            "c": math.nan, "q_floor": 1 / n, "floor_above_alpha": 1 / n > 0.01,
            "alpha": 0.01, "n_called": n}  # fmt: skip
    data.update(over)
    return event("novel_orfs", CLASS_EVENT, data)


def _denovo(cls: str = "antisense", n: int = 40, c: float = 1.0) -> dict[str, Any]:
    return event("novel_orfs", CLASS_EVENT, {
        "regime": "denovo_subgroup", "class": cls, "estimated": True, "reason": None,
        "n_targets": n, "n_target_wins": None, "n_decoys": 0, "paired_share": None, "c": c,
        "q_floor": c / n, "floor_above_alpha": c / n > 0.01, "alpha": 0.01, "n_called": 0,
    })  # fmt: skip


def test_schema_declares_the_block() -> None:
    prop = schema()["properties"].get("novel_orfs")
    assert prop is not None
    assert "required" not in prop or "classes" in prop.get("required", ["classes"])


def test_no_class_event_gives_null_block() -> None:
    m = build()
    assert "novel_orfs" in m and m["novel_orfs"] is None
    validate(m)


def test_block_sorted_by_regime_then_class_and_nan_is_null() -> None:
    evs = [_denovo("antisense"), _picked("intergenic"), _picked("cds_outframe", 120),
           _denovo("intergenic", 120, 0.5)]  # fmt: skip
    m = build(events=evs)
    got = m["novel_orfs"]["classes"]
    assert [(d["regime"], d["class"]) for d in got] == [
        ("denovo_subgroup", "antisense"), ("denovo_subgroup", "intergenic"),
        ("picked_decoy", "cds_outframe"), ("picked_decoy", "intergenic")]  # fmt: skip
    pk = got[-1]
    assert pk["c"] is None  # NaN -> null
    assert pk["q_floor"] == 1 / 150 and pk["n_called"] == 150
    assert set(pk) == {"regime", "class", "estimated", "reason", "n_targets", "n_target_wins",
                       "n_decoys", "paired_share", "c", "q_floor", "floor_above_alpha", "alpha",
                       "n_called"}  # fmt: skip
    validate(m)
    json.dumps(m, allow_nan=False)


def test_resume_duplicate_collapses() -> None:
    one = build(events=[_picked(), _denovo()])["novel_orfs"]
    two = build(events=[_picked(), _denovo(), _picked(), _denovo()])["novel_orfs"]
    assert one == two and len(one["classes"]) == 2


def test_warnings_reach_the_manifest_but_not_the_flags() -> None:
    evs = [
        _denovo(),
        event("novel_orfs", FLOOR_EVENT, {"regime": "denovo_subgroup", "class": "antisense",
                                          "n_targets": 40, "q_floor": 0.025, "alpha": 0.01},
              level="WARNING"),
        event("novel_orfs", MAKEDB_EVENT, {"n_decoy_rows": 7}, level="WARNING"),
    ]  # fmt: skip
    m = build(events=evs)
    assert f"novel_orfs:{FLOOR_EVENT}" in m["warnings"]
    assert f"novel_orfs:{MAKEDB_EVENT}" in m["warnings"]
    assert not any("q_floor" in f or "decoys_not_makedb" in f for f in m["flags"]), m["flags"]
    validate(m)


def test_block_with_unestimated_picked_class_validates() -> None:
    ev = _picked("antisense", 3, estimated=False, reason="decoys_unscored", n_target_wins=0,
                 n_decoys=0, q_floor=None, floor_above_alpha=False, n_called=0)  # fmt: skip
    m = build(events=[ev])
    assert m["novel_orfs"]["classes"][0]["reason"] == "decoys_unscored"
    assert m["novel_orfs"]["classes"][0]["q_floor"] is None
    validate(m)
