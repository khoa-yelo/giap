"""``a3_feasibility.py mzml-check``: the 8 structural rules tied to ``acquisition.json`` (M10-08a).

Spec: docs/specs/M10-08-a3-feasibility.md, "The mzML check" (rules M1-M8), example 11, edge case
"an mzML with 0 MS2 spectra fails M2", exit codes (0 pass, 1 failed check, 2 precondition) and
acceptance criterion 4 ("passes example 11 and fails each of M1-M8 on a fixture that breaks only
that rule; every failure exits 1 and is named in J").

Every mzML here is synthetic and written in ``tmp_path`` (``helpers.a3_kit.write_mzml``); the
binary arrays are not valid base64, so a checker that decodes them fails. The ``J`` layout is the
test author's convention documented in ``helpers/a3_kit.py``. No statistical claim is made here
(the check is deterministic), so there is no calibration test (spec §8.4).
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from helpers.a3_kit import (
    A3_FEASIBILITY,
    ACQUISITION,
    EX11_ACQ,
    EX11_WINDOWS,
    RULES,
    Spec,
    example11_spectra,
    is_canonical_json,
    ms2_indices,
    require,
    with_window,
    write_acquisition,
    write_mzml,
    wt1_spectra,
)
from helpers.slurm_sandbox import REPO


def _check(
    mzml: Path, acq: Path, out: Path, run: str = "WT_1", *extra: str
) -> subprocess.CompletedProcess[str]:
    require(A3_FEASIBILITY, "deliverable 9, mzml-check")
    cmd = [sys.executable, str(A3_FEASIBILITY), "mzml-check", str(mzml), "--acquisition",
           str(acq), "--run", run, "--out", str(out), *extra]  # fmt: skip
    env = {"PATH": "/usr/bin:/bin", "PYTHONPATH": str(REPO / "src"), "LANG": "C.UTF-8"}
    return subprocess.run(cmd, capture_output=True, text=True, timeout=300, env=env, check=False)


def _run_case(
    tmp: Path, spectra: list[Spec], acq: dict[str, Any] | None = None, **kw: Any
) -> tuple[subprocess.CompletedProcess[str], dict[str, Any] | None]:
    mzml = write_mzml(tmp / "WT_1.mzML", spectra, **kw)
    acq_path = write_acquisition(tmp / "acquisition.json", {"WT_1": acq or EX11_ACQ})
    out = tmp / "mzml_check.json"
    res = _check(mzml, acq_path, out)
    return res, (json.loads(out.read_text()) if out.is_file() else None)


def _rules(j: dict[str, Any] | None) -> dict[str, bool]:
    assert j is not None, "mzml-check wrote no J"
    assert set(j["rules"]) == set(RULES), j["rules"]
    return j["rules"]


# ------------------------------------------------------------------ example 11: the pass


@pytest.mark.parametrize("indexed", [False, True], ids=["mzML", "indexedmzML"])
def test_example11_passes_with_every_rule_true(tmp_path: Path, indexed: bool) -> None:
    res, j = _run_case(tmp_path, example11_spectra(), indexed=indexed)
    assert res.returncode == 0, res.stdout + res.stderr
    assert _rules(j) == dict.fromkeys(RULES, True)
    assert j is not None and j["pass"] is True and j["failed"] == []
    assert j["empty_ms2_share"] == 0


def test_example11_j_is_canonical_and_deterministic(tmp_path: Path) -> None:
    mzml = write_mzml(tmp_path / "WT_1.mzML", example11_spectra())
    acq = write_acquisition(tmp_path / "acq.json", {"WT_1": EX11_ACQ})
    a, b = tmp_path / "a.json", tmp_path / "b.json"
    assert _check(mzml, acq, a).returncode == 0
    assert _check(mzml, acq, b).returncode == 0
    assert is_canonical_json(a.read_text()), "J: sorted keys, indent 2, trailing newline"
    assert a.read_bytes() == b.read_bytes()


def test_scan_start_time_in_minutes_is_converted(tmp_path: Path) -> None:
    # msconvert writes minutes (UO:0000031); M7 compares seconds with gradient_length_s.
    res, j = _run_case(tmp_path, example11_spectra(), rt_unit="minute")
    assert res.returncode == 0, res.stdout + res.stderr
    assert _rules(j)["M7"] is True


def test_window_triples_are_rounded_to_0_01(tmp_path: Path) -> None:
    sp = example11_spectra()
    t, lo, up = EX11_WINDOWS[0]
    with_window(sp, ms2_indices(sp)[4], (t + 0.0004, lo - 0.0003, up + 0.0002))
    res, j = _run_case(tmp_path, sp)
    assert res.returncode == 0, res.stdout + res.stderr
    assert _rules(j)["M5"] is True


# ------------------------------------------------------------------ example 11: the M4 failure


def test_example11_one_ms2_widened_to_100_fails_m4(tmp_path: Path) -> None:
    sp = example11_spectra()
    t = EX11_WINDOWS[3][0]
    with_window(sp, ms2_indices(sp)[3], (t, 50.0, 50.0))  # width 100 > 50 = half the span
    res, j = _run_case(tmp_path, sp)
    assert res.returncode == 1, res.stdout + res.stderr
    assert _rules(j)["M4"] is False
    assert j is not None and "M4" in j["failed"] and j["pass"] is False


# ------------------------------------------------------------------ one fixture per rule


def _break(rule: str) -> tuple[list[Spec], dict[str, Any]]:
    """Example 11 changed so that exactly ``rule`` fails."""
    sp = example11_spectra()
    acq = json.loads(json.dumps(EX11_ACQ))
    ms2 = ms2_indices(sp)
    if rule == "M1":  # a 4th MS1 spectrum; acquisition says 3
        sp.insert(5, Spec(1, 3.0))
    elif rule == "M2":  # 13 MS2 > 3 frames x 4 window rows
        sp.append(Spec(2, 9.2, (EX11_WINDOWS[0],)))
        sp.sort(key=lambda s: s.rt)
    elif rule == "M3":  # one MS2 spectrum without an isolation window
        sp[ms2[5]] = replace(sp[ms2[5]], windows=(None,))
    elif rule == "M4":  # every spectrum of one window widened to 100 (still 4 distinct triples)
        t = EX11_WINDOWS[3][0]
        for i in ms2:
            if sp[i].windows[0] == EX11_WINDOWS[3]:
                with_window(sp, i, (t, 50.0, 50.0))
    elif rule == "M5":  # one distinct window < window_groups = 2
        for i in ms2:
            with_window(sp, i, EX11_WINDOWS[0])
    elif rule == "M6":  # one spectrum flagged profile
        sp[ms2[2]] = replace(sp[ms2[2]], mode="profile")
    elif rule == "M7":  # last RT 8.5 s < 0.9 x 10 s
        sp = [replace(s, rt=s.rt - 1.0) if s.rt >= 8.0 else s for s in sp]
    elif rule == "M8":  # one empty MS2 spectrum
        sp[ms2[7]] = replace(sp[ms2[7]], n_peaks=0)
    else:  # pragma: no cover
        raise AssertionError(rule)
    return sp, acq


@pytest.mark.parametrize("rule", RULES)
def test_fixture_breaking_only_one_rule(tmp_path: Path, rule: str) -> None:
    sp, acq = _break(rule)
    res, j = _run_case(tmp_path, sp, acq)
    assert res.returncode == 1, f"{rule}: {res.stdout}{res.stderr}"
    want = dict.fromkeys(RULES, True)
    want[rule] = False
    assert _rules(j) == want
    assert j is not None and j["failed"] == [rule] and j["pass"] is False


def test_m2_lower_bound(tmp_path: Path) -> None:
    acq = {**EX11_ACQ, "msms_type_counts": {"0": 3, "9": 13}}  # 12 MS2 < 13 frames
    res, j = _run_case(tmp_path, example11_spectra(), acq)
    assert res.returncode == 1
    assert _rules(j)["M2"] is False


def test_zero_ms2_spectra_fails_m2(tmp_path: Path) -> None:
    res, j = _run_case(tmp_path, [s for s in example11_spectra() if s.ms_level == 1])
    assert res.returncode == 1
    assert _rules(j)["M2"] is False


@pytest.mark.parametrize(
    "variant", ["two_precursors", "no_precursor_list", "zero_lower_offset", "zero_upper_offset"]
)
def test_m3_variants(tmp_path: Path, variant: str) -> None:
    sp = example11_spectra()
    i = ms2_indices(sp)[6]
    t, lo, up = sp[i].windows[0]  # type: ignore[misc]
    if variant == "two_precursors":
        with_window(sp, i, (t, lo, up), (t, lo, up))
    elif variant == "no_precursor_list":
        with_window(sp, i)
    elif variant == "zero_lower_offset":
        with_window(sp, i, (t, 0.0, up))
    else:
        with_window(sp, i, (t, lo, 0.0))
    res, j = _run_case(tmp_path, sp)
    assert res.returncode == 1
    assert _rules(j)["M3"] is False
    assert j is not None and "M3" in j["failed"]


def test_m5_too_many_windows(tmp_path: Path) -> None:
    sp = example11_spectra()
    with_window(sp, ms2_indices(sp)[0], (450.0, 12.5, 12.5))  # 5 distinct > dia_window_rows 4
    res, j = _run_case(tmp_path, sp)
    assert res.returncode == 1
    assert _rules(j)["M5"] is False


def test_m6_unflagged_spectrum_fails(tmp_path: Path) -> None:
    sp = example11_spectra()
    sp[0] = replace(sp[0], mode=None)  # an MS1 spectrum without centroid flag
    res, j = _run_case(tmp_path, sp)
    assert res.returncode == 1
    assert _rules(j)["M6"] is False


def test_m8_reports_the_empty_share(tmp_path: Path) -> None:
    sp = example11_spectra()
    ms2 = ms2_indices(sp)
    for i in ms2[:3]:
        sp[i] = replace(sp[i], n_peaks=0)
    res, j = _run_case(tmp_path, sp)
    assert res.returncode == 1
    assert _rules(j)["M8"] is False
    assert j is not None and j["empty_ms2_share"] == pytest.approx(3 / 12)


def test_several_failures_are_all_named(tmp_path: Path) -> None:
    sp, acq = _break("M6")
    sp = [replace(s, rt=s.rt - 1.0) if s.rt >= 8.0 else s for s in sp]  # and M7
    res, j = _run_case(tmp_path, sp, acq)
    assert res.returncode == 1
    assert j is not None and j["failed"] == ["M6", "M7"]


# ------------------------------------------------------------------ the committed WT_1 values


@pytest.fixture(scope="module")
def wt1_mzml(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return write_mzml(tmp_path_factory.mktemp("wt1") / "WT_1.mzML", wt1_spectra(), indexed=True)


def test_wt1_shaped_file_passes_against_committed_acquisition(
    tmp_path: Path, wt1_mzml: Path
) -> None:
    out = tmp_path / "J.json"
    res = _check(wt1_mzml, ACQUISITION, out, "WT_1")
    assert res.returncode == 0, res.stdout + res.stderr
    assert json.loads(out.read_text())["rules"] == dict.fromkeys(RULES, True)


def test_rules_use_the_named_run(tmp_path: Path, wt1_mzml: Path) -> None:
    # WT_2 has 8,805 diaPASEF frames: 8,804 MS2 spectra are too few for it (M2).
    out = tmp_path / "J.json"
    res = _check(wt1_mzml, ACQUISITION, out, "WT_2")
    assert res.returncode == 1, res.stdout + res.stderr
    j = json.loads(out.read_text())
    assert j["failed"] == ["M2"]


# ------------------------------------------------------------------ preconditions (exit 2)


def test_missing_mzml_exits_2(tmp_path: Path) -> None:
    acq = write_acquisition(tmp_path / "acq.json", {"WT_1": EX11_ACQ})
    res = _check(tmp_path / "absent.mzML", acq, tmp_path / "J.json")
    assert res.returncode == 2, res.stdout + res.stderr


def test_missing_acquisition_exits_2(tmp_path: Path) -> None:
    mzml = write_mzml(tmp_path / "WT_1.mzML", example11_spectra())
    res = _check(mzml, tmp_path / "absent.json", tmp_path / "J.json")
    assert res.returncode == 2, res.stdout + res.stderr


def test_unknown_run_exits_2(tmp_path: Path) -> None:
    mzml = write_mzml(tmp_path / "WT_1.mzML", example11_spectra())
    acq = write_acquisition(tmp_path / "acq.json", {"WT_1": EX11_ACQ})
    res = _check(mzml, acq, tmp_path / "J.json", "WT_9")
    assert res.returncode == 2, res.stdout + res.stderr


def test_unknown_option_exits_2(tmp_path: Path) -> None:
    mzml = write_mzml(tmp_path / "WT_1.mzML", example11_spectra())
    acq = write_acquisition(tmp_path / "acq.json", {"WT_1": EX11_ACQ})
    res = _check(mzml, acq, tmp_path / "J.json", "WT_1", "--no-such-option")
    assert res.returncode == 2, res.stdout + res.stderr


def test_malformed_xml_is_an_error_not_a_pass(tmp_path: Path) -> None:
    mzml = write_mzml(tmp_path / "WT_1.mzML", example11_spectra())
    mzml.write_text(mzml.read_text()[: len(mzml.read_text()) // 2])  # truncated file
    acq = write_acquisition(tmp_path / "acq.json", {"WT_1": EX11_ACQ})
    res = _check(mzml, acq, tmp_path / "J.json")
    assert res.returncode in (1, 2), res.stdout + res.stderr
    assert "Traceback" not in res.stderr
