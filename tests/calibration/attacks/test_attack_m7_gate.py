"""Teeth of the M7-13a variants gate: each check has a mutant that must fail
(docs/specs/M7-13-calibration-gate.md § Attacks: "Each attack has a mutant that must fail";
attack 6 row "A6c doubles as the teeth of A6b"; GIAP_SPEC.md §8.5 "A failing attack blocks the
merge"). The mutants reuse the gate's own replicates (``helpers.m7_gate_records``: one GIAP run
feeds the check and its mutant, rule 1) and are evaluated through the same criteria, logged to a
private directory (``teeth_log``), never the nightly suite log.

* G1 with the confidence filter removed (every substitution row predicted): a confirmed ``fail``
  of the pooled F1 criterion.
* G2's check on the isobaric attack's de novo reads (error-laden, no strain): violations in at
  least one replicate.
* A6a with ``mass_equivalent`` stripped from every row: unflagged rows in at least one replicate.
* A6b with every swap row taken as ``high``: a confirmed FWER ``fail``.

The novel-ORF attacks (a, b, d, GC, 5i) are M7-13b's.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import _criteria as C
import pytest
import yaml
from helpers import m7_gate_records as M

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("variants")]
NIGHTLY = pytest.mark.calibration_nightly
SUITE = "variants"


def test_smoke_one_attack_replicate_record_is_complete() -> None:
    rec = M.a6_record(0)
    assert set(rec["a6a"]) == set(rec["a6a_mutant"]) == {"n_rows", "n_unflagged"}
    assert set(rec["a6b"]) == set(rec["a6b_mutant"]) == {"folded", "unfolded"}
    assert rec["a6a_mutant"]["n_rows"] == rec["a6a"]["n_rows"]
    assert rec["a6b_mutant"]["folded"] >= rec["a6b"]["folded"]
    assert isinstance(rec["g2_mutant_violations"], int)


@NIGHTLY
def test_g1_mutant_without_confidence_filter_fails(teeth_log: Path) -> None:
    R = C.replicates(SUITE, False)
    bound = float(yaml.safe_load(C.TOLERANCES_PATH.read_text())["f1_min"]["variants_high"])
    rec = C.calibrate_pooled(SUITE, M.G1_BASE, "mutant_any_confidence", bound,
                             lambda i: M.g1_record(i, "100k"), R, log_dir=teeth_log,
                             guard={"n_truth": 50, "n_pred": 50})  # fmt: skip
    assert rec["status"] == "fail", rec


def _a6(R: int) -> list[dict[str, Any]]:
    return [M.a6_record(i) for i in range(R)]


@NIGHTLY
def test_g2_check_finds_internal_high_rows_on_erroneous_reads() -> None:
    got = _a6(C.replicates(SUITE, False))
    assert any(r["g2_mutant_violations"] > 0 for r in got), "G2's check has no teeth"


@NIGHTLY
def test_a6a_mutant_without_mass_equivalent_fails() -> None:
    got = _a6(C.replicates(SUITE, False))
    assert any(r["a6a_mutant"]["n_unflagged"] > 0 for r in got), "A6a's check has no teeth"


@NIGHTLY
def test_a6b_mutant_every_swap_high_fails(teeth_log: Path) -> None:
    R, stream = C.replicates(SUITE, False), "variants.a6b.mutant"
    index = {C.suite_seed(stream, i): i for i in range(2 * R)}

    def fn(rng: Any) -> float:
        i = index[int(rng.bit_generator.seed_seq.entropy)]
        return float(M.a6_record(i)["a6b_mutant"]["folded"] > 0)

    rec = C.calibrate(SUITE, "fwer", 0.05, fn, R, stream=stream, log_dir=teeth_log)
    assert rec["status"] == "fail", rec
