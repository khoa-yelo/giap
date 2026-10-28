"""M9-25 example 2: changing ``style.DETECTION["detected"]`` to another hex makes snapshots fail
(RMS > 2) (docs/specs/M9-25-snapshot-baselines.md; GIPPA_SPEC.md §10.1).

The spec lists "V01, V02, V03, ..."; V03 colours genes by the cividis abundance scale and never
uses ``DETECTION["detected"]`` (gippa.viz.locus), so the plots checked are V01 (circos status
track), V02 (detected markers) and V16 (UpSet bars) (test-author's reading, reported to the lead).

Every snapshot instance of these plots on ``sim_plots`` is drawn in a fresh interpreter in which
``gippa.viz.style`` is loaded from its source with only the ``"detected"`` entry of ``DETECTION``
changed (``STATUS_COLOUR`` and every module that copies the colour at import pick it up), then
compared with the committed baselines by ``matplotlib.testing.compare.compare_images``; each
plot's most affected instance must exceed RMS 2. With ``--mpl`` (the ubuntu 3.12 leg) the
unchanged colour is also drawn and every instance must match (the control).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest
from helpers import snapshot_fixture as sf

PIDS = ("V01", "V02", "V16")
OLD, NEW = "#0072B2", "#D55E00"

_SCRIPT = r"""
import importlib.util, json, sys
from pathlib import Path
old, new, dest = sys.argv[1], sys.argv[2], Path(sys.argv[3])
pids = sys.argv[4:]
import gippa.viz
src = Path(gippa.viz.__file__).with_name("style.py").read_text(encoding="utf-8")
needle = f'"detected": "{old}"'
assert src.count(needle) == 1, "DETECTION['detected'] literal not found once in style.py"
if new != old:
    src = src.replace(needle, f'"detected": "{new}"')
spec = importlib.util.spec_from_loader("gippa.viz.style", loader=None)
mod = importlib.util.module_from_spec(spec)
mod.__file__ = gippa.viz.__file__.replace("__init__.py", "style.py")
sys.modules["gippa.viz.style"] = mod
exec(compile(src, mod.__file__, "exec"), mod.__dict__)
gippa.viz.style = mod
from helpers import snapshot_fixture as sf
cases = {(c.pid, c.stem): c for c in sf.cases() if c.fixture == sf.RICH and c.error is None}
out = {}
import matplotlib.style
for (pid, _), case in cases.items():
    if pid not in pids:
        continue
    with matplotlib.style.context("default", after_reset=True):
        mod.apply()
        sf.figure(case).savefig(dest / case.filename, dpi=sf.SNAPSHOT_DPI)
    out.setdefault(pid, []).append(case.filename)
(dest / "drawn.json").write_text(json.dumps(out), encoding="utf-8")
"""


def _draw(colour: str, dest: Path) -> dict[str, list[str]]:
    dest.mkdir(parents=True)
    env = {**os.environ, "MPLBACKEND": "Agg"}
    paths = [str(sf.REPO / "src"), str(sf.REPO / "tests" / "sim"), str(sf.REPO / "tests")]
    env["PYTHONPATH"] = os.pathsep.join([*paths, env.get("PYTHONPATH", "")])
    cmd = [sys.executable, "-c", _SCRIPT, OLD, colour, str(dest), *PIDS]
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=600)
    assert res.returncode == 0, res.stderr[-4000:]
    drawn: dict[str, list[str]] = json.loads((dest / "drawn.json").read_text(encoding="utf-8"))
    return drawn


def _rms(actual: Path, expected: Path) -> float:
    from matplotlib.testing.compare import compare_images

    res = compare_images(str(expected), str(actual), tol=0, in_decorator=True)
    return 0.0 if res is None else float(res["rms"])


def _require() -> None:
    why = sf.missing(sf.RICH)
    if why is not None:
        pytest.fail(why)
    if not sf.BASELINE_DIR.is_dir():
        pytest.fail("tests/snapshot/baseline is not committed")


def test_detected_colour_change_breaks_v01_v02_v16(tmp_path: Path) -> None:
    _require()
    drawn = _draw(NEW, tmp_path / "new")
    assert sorted(drawn) == list(PIDS)
    # per plot, its most affected instance
    rms = {p: max(_rms(tmp_path / "new" / f, sf.BASELINE_DIR / f) for f in fs)
           for p, fs in drawn.items()}  # fmt: skip
    assert {p: r for p, r in rms.items() if r <= sf.TOLERANCE} == {}, rms


def test_unchanged_colour_matches_the_baselines(
    tmp_path: Path, pytestconfig: pytest.Config
) -> None:
    if not pytestconfig.getoption("--mpl", default=False):
        pytest.skip("the control compares pixels: --mpl leg only (index decision 22)")
    _require()
    drawn = _draw(OLD, tmp_path / "old")
    rms = {f: _rms(tmp_path / "old" / f, sf.BASELINE_DIR / f) for fs in drawn.values() for f in fs}
    assert {f: r for f, r in rms.items() if r > sf.TOLERANCE} == {}, rms
