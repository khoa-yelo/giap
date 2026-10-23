"""M9-26 example 2 and work item 4: two renders give byte-identical ``.svg``, ``.pdf`` and ``.png``
files, and two ``giap report`` calls give a byte-identical ``report.html`` and ``plots/``
(docs/specs/M9-26-scale-and-determinism.md; M9-01 item 7; GIAP_SPEC.md R7/TE7, §7.9.4).

Identity is asserted within one environment only (degenerate row "PNG bytes differ across
matplotlib patch versions": expected).

- Every registered plot is rendered on ``helpers.byte_identity.golden_out`` (the committed
  giapsim sim50k golden tables) in two fresh interpreters with different ``PYTHONHASHSEED``s;
  the trees of images and data TSVs must be equal byte for byte. The M9-25 snapshot fixture runs
  the same comparison in tests/snapshot/test_snapshot_byte_identity.py once it is committed.
- V07 (the spec's example 2) is also rendered twice in one process.
- A Hypothesis property: for small giapsim layouts (``helpers.scale_fixture.stress_out``), V1 and
  V2 rendered twice in one process are byte-identical in all three formats (the layouts carry
  two origin-crossing genes; no shrink phase, M9-26 lead decision 4).
- ``giap report OUT`` is run in two fresh interpreters (different hash seeds) on two copies of
  one ``OUT/``: ``report/report.html`` and every file under ``plots/`` are equal. Marked
  ``slow``.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest
from helpers.byte_identity import (
    FORMATS,
    REPO,
    differing,
    golden_out,
    render_all,
    render_in_subprocess,
    tree_bytes,
)
from hypothesis import HealthCheck, Phase, given, settings
from hypothesis import strategies as st

#: plots the golden sim50k tables are rich enough to draw (each must render, not skip)
GOLDEN_RENDERED = ("V01", "V02", "V03", "V05", "V06", "V07", "V08", "V10", "V17", "V18")


@pytest.fixture(scope="module")
def golden(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return golden_out(tmp_path_factory.mktemp("golden50k"))


def test_every_plot_byte_identical_across_two_interpreters(golden: Path, tmp_path: Path) -> None:
    a = render_in_subprocess(golden, tmp_path / "a", hashseed=1)
    b = render_in_subprocess(golden, tmp_path / "b", hashseed=2)
    assert a == b, "render status differs between the two runs"
    missing = [p for p in GOLDEN_RENDERED if a.get(p) != "rendered"]
    assert not missing, f"golden fixture no longer renders {missing}: {a}"
    ta, tb = tree_bytes(tmp_path / "a"), tree_bytes(tmp_path / "b")
    for pid in GOLDEN_RENDERED:
        for fmt in FORMATS:
            assert any(k.startswith(f"{pid}_") and k.endswith(f".{fmt}") for k in ta), (pid, fmt)
    assert differing(ta, tb) == []


def test_example_2_v07_two_renders_identical_in_all_formats(golden: Path, tmp_path: Path) -> None:
    assert render_all(golden, tmp_path / "a", ["V07"]) == {"V07": "rendered"}
    assert render_all(golden, tmp_path / "b", ["V07"]) == {"V07": "rendered"}
    ta, tb = tree_bytes(tmp_path / "a"), tree_bytes(tmp_path / "b")
    assert {k.rsplit(".", 1)[-1] for k in ta} == {*FORMATS, "tsv"}
    assert differing(ta, tb) == []


@settings(
    max_examples=6,
    deadline=None,
    # no shrink phase: a failing example re-renders V1/V2 per shrink step (≈ 333 s, review N6)
    phases=[Phase.explicit, Phase.generate],
    suppress_health_check=[HealthCheck.function_scoped_fixture, HealthCheck.too_slow],
)
@given(n_genes=st.integers(20, 300), seed=st.integers(0, 2**31 - 1))
def test_property_v01_v02_two_renders_identical(
    tmp_path_factory: pytest.TempPathFactory, n_genes: int, seed: int
) -> None:
    from helpers.scale_fixture import stress_out

    tmp = tmp_path_factory.mktemp("prop")
    run = stress_out(tmp, seed=seed, n_genes=n_genes)
    first = render_all(run.out, tmp / "a", ["V01", "V02"])
    second = render_all(run.out, tmp / "b", ["V01", "V02"])
    assert first == second == {"V01": "rendered", "V02": "rendered"}
    assert differing(tree_bytes(tmp / "a"), tree_bytes(tmp / "b")) == []


def _giap_report(out: Path, hashseed: int, db: Path) -> None:
    env = {**os.environ, "PYTHONHASHSEED": str(hashseed), "MPLBACKEND": "Agg"}
    env["GIAP_DB_DIR"] = str(db)
    paths = [str(REPO / "src"), str(REPO / "tests" / "sim"), str(REPO / "tests")]
    env["PYTHONPATH"] = os.pathsep.join([*paths, env.get("PYTHONPATH", "")])
    cmd = [sys.executable, "-c", "from giap.cli import cli; cli()", "report", str(out)]
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=1800)
    assert res.returncode == 0, f"{res.stdout[-4000:]}\n{res.stderr[-4000:]}"


@pytest.mark.slow
def test_giap_report_twice_byte_identical_report_and_images(golden: Path, tmp_path: Path) -> None:
    db = tmp_path / "no_db"  # never the real $GIAP_DB_DIR; giap report reads no database
    db.mkdir()
    outs = []
    for i, seed in enumerate((11, 12)):
        out = tmp_path / f"run{i}" / "OUT"
        shutil.copytree(golden, out)
        _giap_report(out, seed, db)
        outs.append(out)
    pages = [(o / "report" / "report.html").read_bytes() for o in outs]
    assert pages[0], "empty report.html"
    assert pages[0] == pages[1], "report.html differs between two giap report calls"
    pa, pb = tree_bytes(outs[0] / "plots"), tree_bytes(outs[1] / "plots")
    assert any(k.endswith(".pdf") for k in pa) and any(k.endswith(".png") for k in pa)
    assert differing(pa, pb) == []
