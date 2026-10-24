"""``validation/pilot/acquisition_mode.py``: the acquisition mode of every run (M10-01).

Spec: docs/specs/M10-01-pxd078459-download.md, "Acquisition mode", examples 9-10, edge cases
(no ``analysis.tdf``, missing ``analysis.tdf_bin``, diaPASEF variants, empty ``Frames``),
acceptance criterion 5; question Q4 (all three runs are checked).

``MsMsType`` values (0 = MS1, 8 = ddaPASEF MS/MS, 9 = diaPASEF MS/MS) are taken from the AlphaTims
documentation of the TDF format; the fixtures are SQLite files built at test time by
``helpers.pxd_fixtures.make_tdf``, never real data, and no TDF reader is imported.

Interface settled by the test author:

* ``python -I acquisition_mode.py RUN.d [RUN.d ...] --out FILE.json``. The JSON is
  ``{"dataset_mode": ..., "runs": [...]}`` with one object per argument, in argument order, holding
  at least ``run`` (directory name without ``.d``), ``mode``, ``msms_type_counts`` (``{"0": n,
  "9": n}``: only the types present), ``window_groups`` (distinct ``WindowGroup`` values),
  ``mz_range`` (``[min(IsolationMz - IsolationWidth/2), max(IsolationMz + IsolationWidth/2)]``),
  ``gradient_length_s`` (max ``Frames.Time``), ``metadata`` (the ``GlobalMetadata`` keys of the
  spec that are present) and ``incomplete_run`` (bool).
* Exit 0 only for ``dataset_mode = diaPASEF``; 3 for ``ddaPASEF`` (message names §17 and M10-02);
  1 otherwise. The JSON is written in every case.

No statistical claim is made in this module, so it has no calibration test (spec §8.4).
"""

from __future__ import annotations

import hashlib
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import pytest
from helpers.pxd_fixtures import ACQUISITION_MODE, make_run
from hypothesis import given, settings
from hypothesis import strategies as st

DIA = {"msms": {0: 1200, 9: 10800}, "dia_groups": 12, "windows_per_group": 2}  # example 9
DDA = {"msms": {0: 600, 8: 2400}, "dia_groups": 0, "pasef_rows": 40}


def _run(*args: str | Path) -> subprocess.CompletedProcess[str]:
    assert ACQUISITION_MODE.is_file(), f"{ACQUISITION_MODE} does not exist yet (deliverable 4)"
    cmd = [sys.executable, "-I", str(ACQUISITION_MODE), *map(str, args)]
    return subprocess.run(cmd, capture_output=True, text=True, timeout=120, check=False)


def _check(
    tmp: Path, runs: dict[str, dict[str, Any]]
) -> tuple[subprocess.CompletedProcess[str], dict[str, Any]]:
    dirs = []
    for name, kw in runs.items():
        d = tmp / f"{name}.d"
        make_run(d, **kw)
        dirs.append(d)
    out = tmp / "acquisition.json"
    res = _run(*dirs, "--out", out)
    assert out.is_file(), f"no JSON written\n{res.stdout}\n{res.stderr}"
    return res, json.loads(out.read_text())


def _modes(doc: dict[str, Any]) -> list[str]:
    return [r["mode"] for r in doc["runs"]]


def test_example9_three_diapasef_runs_exit_0(tmp_path: Path) -> None:
    res, doc = _check(tmp_path, {f"WT_{i}": DIA for i in (1, 2, 3)})
    assert res.returncode == 0, res.stdout + res.stderr
    assert doc["dataset_mode"] == "diaPASEF"
    assert [r["run"] for r in doc["runs"]] == ["WT_1", "WT_2", "WT_3"]
    assert _modes(doc) == ["diaPASEF"] * 3


def test_example9_window_scheme_and_metadata(tmp_path: Path) -> None:
    d = tmp_path / "WT_1.d"
    info = make_run(d, **DIA)
    out = tmp_path / "a.json"
    res = _run(d, "--out", out)
    assert res.returncode == 0, res.stdout + res.stderr
    run = json.loads(out.read_text())["runs"][0]
    assert run["msms_type_counts"] == {"0": 1200, "9": 10800}
    assert run["window_groups"] == 12
    assert run["mz_range"] == pytest.approx(info["mz_range"])
    assert run["gradient_length_s"] == pytest.approx(info["gradient_length_s"])
    for key in ("InstrumentName", "InstrumentFamily", "AcquisitionSoftwareVersion",
                "AcquisitionDateTime"):  # fmt: skip
        assert run["metadata"][key] == info["metadata"][key], key
    assert run["incomplete_run"] is False


def test_example10_one_dda_run_makes_the_dataset_inconsistent(tmp_path: Path) -> None:
    res, doc = _check(tmp_path, {"WT_1": DIA, "WT_2": DIA, "WT_3": DDA})
    assert _modes(doc) == ["diaPASEF", "diaPASEF", "ddaPASEF"]
    assert doc["dataset_mode"] == "inconsistent"
    assert res.returncode == 1, res.stdout + res.stderr


def test_all_ddapasef_exits_3_and_names_section_17(tmp_path: Path) -> None:
    res, doc = _check(tmp_path, {f"WT_{i}": DDA for i in (1, 2, 3)})
    assert doc["dataset_mode"] == "ddaPASEF"
    assert res.returncode == 3, res.stdout + res.stderr
    msg = res.stdout + res.stderr
    assert "§17" in msg and "M10-02" in msg, msg


def test_ms1_only(tmp_path: Path) -> None:
    res, doc = _check(tmp_path, {f"WT_{i}": {"msms": {0: 50}, "dia_groups": 0} for i in (1, 2, 3)})
    assert _modes(doc) == ["ms1_only"] * 3
    assert doc["dataset_mode"] == "ms1_only"
    assert res.returncode == 1


def test_mixed_run(tmp_path: Path) -> None:
    mixed = {"msms": {0: 10, 8: 10, 9: 10}, "dia_groups": 4, "pasef_rows": 5}
    res, doc = _check(tmp_path, {f"WT_{i}": mixed for i in (1, 2, 3)})
    assert _modes(doc) == ["mixed"] * 3
    assert doc["dataset_mode"] == "mixed"
    assert res.returncode == 1


def test_type9_without_windows_is_unknown(tmp_path: Path) -> None:
    runs = {f"WT_{i}": {"msms": {0: 10, 9: 30}, "dia_groups": 0} for i in (1, 2, 3)}
    res, doc = _check(tmp_path, runs)
    assert _modes(doc) == ["unknown"] * 3
    assert res.returncode == 1


def test_empty_frames_is_unknown(tmp_path: Path) -> None:
    res, doc = _check(tmp_path, {f"WT_{i}": {"msms": {}, "dia_groups": 0} for i in (1, 2, 3)})
    assert _modes(doc) == ["unknown"] * 3
    assert res.returncode == 1


def test_run_without_analysis_tdf_is_not_a_bruker_run(tmp_path: Path) -> None:
    dirs = []
    for i in (1, 2):
        make_run(tmp_path / f"WT_{i}.d", **DIA)
        dirs.append(tmp_path / f"WT_{i}.d")
    bad = tmp_path / "WT_3.d"
    bad.mkdir()
    (bad / "WT_3.mzML").write_text("<mzML/>\n")
    dirs.append(bad)
    out = tmp_path / "a.json"
    res = _run(*dirs, "--out", out)
    assert res.returncode == 1, res.stdout + res.stderr
    assert "not a Bruker TDF run" in res.stdout + res.stderr
    doc = json.loads(out.read_text())
    assert _modes(doc)[2] == "unknown"
    assert doc["dataset_mode"] == "inconsistent"


def test_missing_tdf_bin_is_incomplete_and_exits_1(tmp_path: Path) -> None:
    dirs = []
    for i in (1, 2, 3):
        make_run(tmp_path / f"WT_{i}.d", tdf_bin=i != 2, **DIA)
        dirs.append(tmp_path / f"WT_{i}.d")
    out = tmp_path / "a.json"
    res = _run(*dirs, "--out", out)
    assert res.returncode == 1, res.stdout + res.stderr
    doc = json.loads(out.read_text())
    assert _modes(doc) == ["diaPASEF"] * 3, "the metadata still read without the binary"
    assert [r["incomplete_run"] for r in doc["runs"]] == [False, True, False]


def test_never_writes_to_the_database(tmp_path: Path) -> None:
    """Acceptance 5: mtime and sha256 unchanged; no journal or WAL file appears."""
    dirs = []
    for i in (1, 2, 3):
        make_run(tmp_path / f"WT_{i}.d", **DIA)
        dirs.append(tmp_path / f"WT_{i}.d")
    tdf = dirs[0] / "analysis.tdf"
    before = (tdf.stat().st_mtime_ns, hashlib.sha256(tdf.read_bytes()).hexdigest())
    listing = sorted(p.name for p in dirs[0].iterdir())
    tdf.chmod(0o444)  # a writer would fail loudly on a read-only file
    try:
        res = _run(*dirs, "--out", tmp_path / "a.json")
    finally:
        tdf.chmod(0o644)
    assert res.returncode == 0, res.stdout + res.stderr
    after = (tdf.stat().st_mtime_ns, hashlib.sha256(tdf.read_bytes()).hexdigest())
    assert after == before
    assert sorted(p.name for p in dirs[0].iterdir()) == listing


def test_opens_the_database_read_only() -> None:
    assert ACQUISITION_MODE.is_file(), ACQUISITION_MODE
    code = ACQUISITION_MODE.read_text()
    assert re.search(r"mode=ro", code), "sqlite3.connect('file:…?mode=ro', uri=True) (spec)"
    imports = re.findall(r"^\s*(?:from|import)\s+([\w.]+)", code, re.M)
    roots = {i.split(".")[0] for i in imports}
    non_std = sorted(r for r in roots if r not in sys.stdlib_module_names and r != "__future__")
    assert not non_std, f"deliverable 4 is stdlib only: {non_std}"


# ------------------------------------------------------------------ property: the mode rule


def _expected(counts: dict[int, int], dia_rows: int, pasef_rows: int) -> str:
    present = {t for t, n in counts.items() if n}
    if 8 in present and 9 in present:
        return "mixed"
    if 9 in present and dia_rows:
        return "diaPASEF"
    if 8 in present and pasef_rows:
        return "ddaPASEF"
    if present == {0}:
        return "ms1_only"
    return "unknown"


@settings(max_examples=25, deadline=None)
@given(
    st.integers(0, 4),
    st.integers(0, 4),
    st.integers(0, 4),
    st.integers(0, 3),
    st.integers(0, 3),
)
def test_mode_rule(n0: int, n8: int, n9: int, groups: int, pasef: int) -> None:
    counts = {0: n0, 8: n8, 9: n9}
    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp) / "WT_1.d"
        make_run(d, msms={t: n for t, n in counts.items() if n}, dia_groups=groups,
                 windows_per_group=1, pasef_rows=pasef)  # fmt: skip
        out = Path(tmp) / "a.json"
        res = _run(d, "--out", out)
        doc = json.loads(out.read_text())
    want = _expected(counts, groups, pasef)
    assert doc["runs"][0]["mode"] == want, (counts, groups, pasef)
    assert doc["dataset_mode"] == want
    assert res.returncode == {"diaPASEF": 0, "ddaPASEF": 3}.get(want, 1)
