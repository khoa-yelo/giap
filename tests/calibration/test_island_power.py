"""M5-07: the island power recorder ``scripts/dev/island_power.py`` (GIPPA_SPEC.md §7.3.6 "Power
(recorded, not gated): islands of 5, 10 and 20 TUs with effects of 1, 1.5 and 2 SD"; spec
docs/specs/M5-07-island-calibration-power.md "Public interface", Examples 4-5, lead decision 30).

Power is recorded, never gated: these tests check the recorder's grid, seeds, metric definitions
and JSON lines, not power values. The one claim asserted on values is the structural zero of the
low direction (lead decision 28 / ADR-0021 Amendment 1, M5-06 decision 50): a low-direction cell
never has more same-direction hits than the null allows.

Interface assumed of the script (implementer; importable as a module from its path):

``N_GENES = 4500``, ``RHO = 0.5``, ``R_DEFAULT = 100``
``CELLS``: the 18 ``(n_tus, effect_sd, direction)`` of {5, 10, 20} x {1, 1.5, 2} x {high, low};
``cell_seed(index) -> int``: ``calibration_base + index`` (tolerances.yaml ``seeds``);
``scenario(n_tus, effect_sd, direction, *, detection="covariates_only") -> (LayoutConfig,
ExpressionConfig)``: ``scenarios.island_power`` on ``N_GENES`` at ``RHO`` with the effect signed by
the direction, and the planted island kept off every region (spec "Edge cases");
``replicate_metrics(islands, planted, direction) -> (hit, jaccard, false_island)``: ``islands`` a
sequence of ``(direction, frozenset of TU IDs)``, ``planted`` the planted TU IDs; ``hit`` = some
same-direction island shares >= 1 TU with the planted one, ``jaccard`` = TU Jaccard of the best
such island (NaN without a hit), ``false_island`` = some island shares no TU with it;
``run_cell(n_tus, effect_sd, direction, R, seed, *, detection="covariates_only") -> dict``;
``main(argv) -> int`` with ``--replicates``, ``--cells`` (comma-separated indices into ``CELLS``),
``--detection``, ``--out`` (default ``local/logs/calibration/islands_power.jsonl``), ``--run-id``;
one JSON line per cell with ``n_tus``, ``effect_sd``, ``direction``, ``detection``, ``R``,
``power``, ``mean_jaccard``, ``false_island_rate``, ``commit``, ``run_id``.
"""

from __future__ import annotations

import importlib.util
import json
import math
import sys
from functools import cache
from pathlib import Path
from types import ModuleType

import _criteria as C
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "dev" / "island_power.py"
FIELDS = {
    "n_tus", "effect_sd", "direction", "detection", "R", "power", "mean_jaccard",
    "false_island_rate", "commit", "run_id",
}  # fmt: skip

pytestmark = pytest.mark.calsuite("islands")


@cache
def power_module() -> ModuleType:
    assert SCRIPT.is_file(), f"{SCRIPT} is missing (M5-07 power recorder)"
    spec = importlib.util.spec_from_file_location("island_power", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules["island_power"] = mod
    spec.loader.exec_module(mod)
    return mod


# --------------------------------------------------------------------------- grid and seeds


def test_grid_is_the_spec_grid() -> None:
    m = power_module()
    assert (m.N_GENES, m.RHO, m.R_DEFAULT) == (4500, 0.5, 100)
    cells = [(int(n), float(e), str(d)) for n, e, d in m.CELLS]
    assert len(cells) == 18
    assert set(cells) == {
        (n, e, d) for n in (5, 10, 20) for e in (1.0, 1.5, 2.0) for d in ("high", "low")
    }


def test_cell_seeds_are_calibration_base_plus_index() -> None:
    m = power_module()
    base = C.load_tolerances().calibration_base
    assert [m.cell_seed(i) for i in range(len(m.CELLS))] == [base + i for i in range(18)]


@pytest.mark.parametrize("direction", ["high", "low"])
@pytest.mark.parametrize("detection", ["covariates_only", "abundance"])
def test_scenario_is_island_power_with_signed_effect(direction: str, detection: str) -> None:
    from gippasim import scenarios as sc

    m = power_module()
    lcfg, ecfg = m.scenario(10, 1.5, direction, detection=detection)
    (isl,) = ecfg.islands
    assert isl.replicon == "chr1"
    assert isl.n_tus == 10
    assert isl.effect == (1.5 if direction == "high" else -1.5)
    assert isl.first_tu == sc.island_first_tu(4500, lcfg.mean_tu_size)
    assert ecfg.ar1_rho == 0.5
    assert lcfg.n_genes_chromosome == 4500
    assert ecfg.detection == detection


@pytest.mark.parametrize("seed", [0, 1, 2])
def test_planted_island_avoids_regions(seed: int) -> None:
    """Spec "Edge cases": the planted island shares no gene with any region (region effects must
    not add power). Checked on drawn layouts of the largest island (20 TUs)."""
    import numpy as np

    import gippasim.layout as la

    m = power_module()
    lcfg, ecfg = m.scenario(20, 2.0, "high")
    (isl,) = ecfg.islands
    lay = la.simulate_layout(lcfg, np.random.default_rng(seed))
    chr_tus = lay.tus[lay.tus["contig"] == "chr1"].sort_values("start", kind="mergesort")
    planted = chr_tus.iloc[isl.first_tu : isl.first_tu + isl.n_tus]
    genes = {g for gs in planted["genes"] for g in gs}
    assert len(planted) == 20
    assert genes.isdisjoint(set(lay.region_genes["locus_tag"]))


# --------------------------------------------------------------------------- metric definitions

P = frozenset({"T1", "T2", "T3", "T4"})


def test_metrics_no_island() -> None:
    hit, jac, false = power_module().replicate_metrics([], P, "high")
    assert (hit, false) == (False, False)
    assert math.isnan(jac)


def test_metrics_same_direction_overlap_is_a_hit() -> None:
    hit, jac, false = power_module().replicate_metrics(
        [("high", frozenset({"T3", "T4", "T5", "T6"}))], P, "high"
    )
    assert hit and not false
    assert jac == pytest.approx(2 / 6)


def test_metrics_opposite_direction_is_not_a_hit() -> None:
    hit, jac, false = power_module().replicate_metrics([("low", P)], P, "high")
    assert not hit and not false  # it touches the planted island, so it is not a false island
    assert math.isnan(jac)


def test_metrics_best_of_several_and_a_false_island() -> None:
    islands = [
        ("high", frozenset({"T4", "T5"})),
        ("high", frozenset({"T1", "T2", "T3"})),
        ("high", frozenset({"T90", "T91"})),
    ]
    hit, jac, false = power_module().replicate_metrics(islands, P, "high")
    assert hit and false
    assert jac == pytest.approx(3 / 4)


tu_sets = st.frozensets(st.sampled_from([f"T{i}" for i in range(12)]), min_size=1, max_size=8)


@settings(max_examples=200, deadline=None)
@given(
    islands=st.lists(st.tuples(st.sampled_from(["high", "low"]), tu_sets), max_size=5),
    planted=tu_sets,
    direction=st.sampled_from(["high", "low"]),
)
def test_metrics_properties(
    islands: list[tuple[str, frozenset[str]]], planted: frozenset[str], direction: str
) -> None:
    """hit iff a same-direction island meets the planted TUs; Jaccard in (0, 1] on hits and the
    maximum over them; false iff some island misses the planted TUs."""
    hit, jac, false = power_module().replicate_metrics(islands, planted, direction)
    same = [s for d, s in islands if d == direction and s & planted]
    assert hit == bool(same)
    assert false == any(not (s & planted) for _, s in islands)
    if hit:
        assert jac == pytest.approx(max(len(s & planted) / len(s | planted) for s in same))
        assert 0.0 < jac <= 1.0
    else:
        assert math.isnan(jac)


# --------------------------------------------------------------------------- running cells


def _lines(path: Path) -> list[dict[str, object]]:
    return [json.loads(x) for x in path.read_text().splitlines() if x.strip()]


@pytest.mark.calibration
def test_main_writes_one_line_per_cell(tmp_path: Path) -> None:
    m = power_module()
    out = tmp_path / "islands_power.jsonl"
    rc = m.main(["--replicates", "2", "--cells", "0,17", "--out", str(out), "--run-id", "t1"])
    assert rc == 0
    lines = _lines(out)
    assert len(lines) == 2
    for line, idx in zip(lines, (0, 17), strict=True):
        assert set(line) >= FIELDS
        n, e, d = m.CELLS[idx]
        assert (line["n_tus"], line["effect_sd"], line["direction"]) == (n, e, d)
        assert line["R"] == 2
        assert line["detection"] == "covariates_only"
        assert line["run_id"] == "t1"
        assert isinstance(line["commit"], str) and line["commit"]
        for key in ("power", "false_island_rate"):
            assert 0.0 <= float(line[key]) <= 1.0  # type: ignore[arg-type]


@pytest.mark.calibration
def test_main_detection_abundance_grid(tmp_path: Path) -> None:
    """Lead decision 30: the grid also runs with ``detection = "abundance"``."""
    m = power_module()
    out = tmp_path / "p.jsonl"
    assert m.main(["--replicates", "1", "--cells", "5", "--detection", "abundance",
                   "--out", str(out)]) == 0  # fmt: skip
    (line,) = _lines(out)
    assert line["detection"] == "abundance"


@pytest.mark.calibration
def test_run_cell_same_seed_same_record() -> None:
    m = power_module()
    keys = ("power", "mean_jaccard", "false_island_rate", "R")
    a = m.run_cell(10, 2.0, "high", 2, m.cell_seed(0))
    b = m.run_cell(10, 2.0, "high", 2, m.cell_seed(0))
    assert [a[k] for k in keys] == pytest.approx([b[k] for k in keys], nan_ok=True)


@pytest.mark.calibration
def test_teeth_recorder_counts_a_huge_planted_island() -> None:
    """Teeth for the recorder: a 20-TU island at +6 SD (far beyond the grid) is found in every
    replicate, so a recorder that never counts a hit (or counts the wrong direction) fails."""
    m = power_module()
    rec = m.run_cell(20, 6.0, "high", 3, m.cell_seed(0))
    print(f"teeth power: {rec}")
    assert rec["power"] == 1.0
    assert 0.0 < float(rec["mean_jaccard"]) <= 1.0


@pytest.mark.calibration
def test_low_direction_is_a_structural_zero() -> None:
    """Lead decision 28 / decision 50: the low row is a structural zero, not low power. The
    largest low cell (20 TUs, -2 SD) has no more same-direction hits than the null FWER bound."""
    m = power_module()
    R = 20  # fast tier (< 60 s); the full grid at R = 100 is the Slurm recorder run
    rec = m.run_cell(20, 2.0, "low", R, m.cell_seed(17))
    print(f"low 20 TUs 2 SD: {rec}")
    assert float(rec["power"]) <= C.level_bound(0.05, R)
