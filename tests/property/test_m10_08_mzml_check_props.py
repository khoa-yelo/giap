"""Hypothesis properties of ``a3_feasibility.py mzml-check`` (M10-08a; spec "The mzML check").

Invariants (from the rule table, which counts and bounds but never depends on spectrum order):

* a synthetic file built to satisfy M1-M8 for a random acquisition stub passes, every rule true;
* J's verdict fields do not depend on the order of the spectra in the file (the rules are
  counts, sets and maxima);
* widening every spectrum of one window above half the acquisition span fails M4 and only M4.

``main(argv) -> int`` is called in-process (test-author convention, ``helpers/a3_kit.py``).
"""

from __future__ import annotations

import json
import random
import tempfile
from functools import cache
from pathlib import Path
from typing import Any

from helpers.a3_kit import A3_FEASIBILITY, RULES, Spec, require, write_acquisition, write_mzml
from helpers.pxd_fixtures import load
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st


@cache
def _main() -> Any:
    return load(require(A3_FEASIBILITY, "deliverable 9, mzml-check")).main


_KEYS = ("rules", "pass", "failed", "empty_ms2_share")


@st.composite
def valid_case(draw: st.DrawFn) -> tuple[list[Spec], dict[str, Any]]:
    cycles = draw(st.integers(1, 5))
    n_win = draw(st.integers(1, 6))
    lo_mz = draw(st.floats(300, 600))
    span = draw(st.floats(100, 900))
    width = draw(st.floats(0.5, span / 2 / n_win - 0.1))
    step = span / n_win
    wins = [(round(lo_mz + step * (k + 0.5), 2), round(width / 2, 2), round(width / 2, 2))
            for k in range(n_win)]  # fmt: skip
    gradient = draw(st.floats(5, 2000))
    sp: list[Spec] = []
    for c in range(cycles):
        t = gradient * (c + 1) / cycles * 0.97
        sp.append(Spec(1, t - gradient * 0.01))
        sp += [Spec(2, t - gradient * 0.001 * (j + 1), (w,), n_peaks=draw(st.integers(1, 9)))
               for j, w in enumerate(wins)]  # fmt: skip
    groups = draw(st.integers(1, n_win))
    acq = {
        "msms_type_counts": {"0": cycles, "9": cycles},
        "window_groups": groups,
        "dia_window_rows": n_win,
        "mz_range": [lo_mz, lo_mz + span],
        "gradient_length_s": gradient,
    }
    return sp, acq


def _check(sp: list[Spec], acq: dict[str, Any], d: Path, name: str) -> tuple[int, dict[str, Any]]:
    mzml = write_mzml(d / f"{name}.mzML", sp)
    acq_path = write_acquisition(d / f"{name}.acq.json", {"WT_1": acq})
    out = d / f"{name}.json"
    rc = _main()(["mzml-check", str(mzml), "--acquisition", str(acq_path), "--run", "WT_1",
                  "--out", str(out)])  # fmt: skip
    return rc, json.loads(out.read_text())


_SETTINGS = settings(max_examples=25, deadline=None,
                     suppress_health_check=[HealthCheck.too_slow])  # fmt: skip


@_SETTINGS
@given(valid_case(), st.randoms(use_true_random=False))
def test_valid_files_pass_and_order_does_not_matter(
    case: tuple[list[Spec], dict[str, Any]], rnd: random.Random
) -> None:
    sp, acq = case
    with tempfile.TemporaryDirectory() as tmp:
        rc, j = _check(sp, acq, Path(tmp), "a")
        assert rc == 0, j
        assert j["rules"] == dict.fromkeys(RULES, True)
        shuffled = list(sp)
        rnd.shuffle(shuffled)
        rc2, j2 = _check(shuffled, acq, Path(tmp), "b")
        assert rc2 == 0
        assert {k: j2[k] for k in _KEYS} == {k: j[k] for k in _KEYS}


@_SETTINGS
@given(valid_case(), st.data())
def test_widening_one_window_fails_only_m4(
    case: tuple[list[Spec], dict[str, Any]], data: st.DataObject
) -> None:
    sp, acq = case
    span = acq["mz_range"][1] - acq["mz_range"][0]
    wins = sorted({s.windows[0] for s in sp if s.ms_level == 2})  # type: ignore[index]
    target = data.draw(st.sampled_from(wins))
    half = round(span / 4 + 1.0, 2)  # width 2 x half > span / 2
    wide = (target[0], half, half)
    sp = [Spec(s.ms_level, s.rt, (wide,), s.mode, s.n_peaks)
          if s.ms_level == 2 and s.windows[0] == target else s for s in sp]  # fmt: skip
    with tempfile.TemporaryDirectory() as tmp:
        rc, j = _check(sp, acq, Path(tmp), "w")
        assert rc == 1
        want = dict.fromkeys(RULES, True)
        want["M4"] = False
        assert j["rules"] == want
