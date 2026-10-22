"""M5-07 island power recorder (GIAP_SPEC.md §7.3.6 "Power (recorded, not gated)"; spec
docs/specs/M5-07-island-calibration-power.md "Public interface", lead decision 30).

Grid: n_tus in {5, 10, 20} x effect in {1, 1.5, 2} SD x direction in {high, low} on the 4,500-gene
``scenarios.island_power`` layout at rho 0.5, R = 100 per cell, cell seed ``calibration_base`` +
cell index; replicate r draws from ``giap.rng.for_step(f"islands_power.r{r}", cell seed)``. The
truth (layout + expression) is the simulator's, the scan is the real M5-06 step
(``giap.context.islands``, ``stats.bootstrap = stats.permutations = 1000``, default
``context.islands``). One JSON line per cell goes to ``local/logs/calibration/islands_power.jsonl``.
Runs inside a Slurm job (scripts/slurm/), never in the session allocation.

    PYTHONPATH=src:tests/sim python scripts/dev/island_power.py [--replicates R] [--cells 0,5]
        [--detection covariates_only|abundance] [--out PATH] [--run-id ID]
"""

from __future__ import annotations

import argparse
import copy
import dataclasses
import json
import math
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
for _p in ("src", "tests/sim", "tests/helpers", "tests/calibration"):
    if str(REPO / _p) not in sys.path:
        sys.path.insert(0, str(REPO / _p))

import _criteria as C  # noqa: E402
import islands_data as idat  # noqa: E402

from giap.rng import for_step  # noqa: E402

N_GENES = 4500
RHO = 0.5
R_DEFAULT = 100
B = PERM = 1000  # as the FWER gate: stats.bootstrap = stats.permutations = 1000
ALPHA = 0.05
STREAM = "islands_power"
CELLS: tuple[tuple[int, float, str], ...] = tuple(
    (n, e, d) for d in ("high", "low") for n in (5, 10, 20) for e in (1.0, 1.5, 2.0)
)
DEFAULT_OUT = REPO / "local" / "logs" / "calibration" / "islands_power.jsonl"

Island = tuple[str, frozenset[str]]


def cell_seed(index: int) -> int:
    """``calibration_base + index`` (tolerances.yaml ``seeds``)."""
    return C.load_tolerances().calibration_base + index


def scenario(
    n_tus: int, effect_sd: float, direction: str, *, detection: str = "covariates_only"
) -> Any:
    """``scenarios.island_power`` at ``N_GENES`` / ``RHO`` with the effect signed by ``direction``.

    The layout carries no region (spec "Edge cases": region effects must not add power, and the
    default regions are placed at random, so only an empty region set keeps every layout drawn
    off the island).
    """
    from giapsim import scenarios as sc

    if direction not in ("high", "low"):
        raise ValueError(f"unknown direction {direction!r}")
    signed = effect_sd if direction == "high" else -effect_sd
    lcfg, ecfg = sc.island_power(N_GENES, RHO, n_tus, signed, detection=detection)  # type: ignore[arg-type]
    return dataclasses.replace(lcfg, regions=()), ecfg


def replicate_metrics(
    islands: Sequence[Island], planted: frozenset[str], direction: str
) -> tuple[bool, float, bool]:
    """(hit, Jaccard of the best same-direction overlapping island or NaN, false island)."""
    jac = [
        len(s & planted) / len(s | planted) for d, s in islands if d == direction and s & planted
    ]
    false = any(not (s & planted) for _, s in islands)
    return bool(jac), (max(jac) if jac else math.nan), false


def _replicate(rng: Any, pair: Any) -> tuple[list[Island], frozenset[str]]:
    """Reported islands (direction, TU IDs) and the planted TU IDs of one replicate."""
    import giapsim.expression as ex
    import giapsim.layout as la

    truth_rng = copy.deepcopy(rng)  # the same draws scenario_contigs makes, to read the truth
    lay = la.simulate_layout(pair[0], truth_rng)
    eg = ex.simulate_expression(lay, pair[1], truth_rng).genes
    planted_genes = set(eg.loc[eg["island_id"].notna(), "locus_tag"].astype(str))
    planted = frozenset(
        str(t) for t, gs in zip(lay.tus["tu_id"], lay.tus["genes"], strict=True)
        if planted_genes & set(gs)
    )  # fmt: skip
    contigs = idat.scenario_contigs(rng, pair)
    cfg = idat.config(bootstrap=B, permutations=PERM, alpha=ALPHA)
    genes, tus, replicons = idat.frames(contigs)
    out, _ = idat.step(genes, idat.build_series(genes, tus, replicons, cfg), cfg, rng)
    tu_of = dict(zip(genes["locus_tag"].astype(str), genes["tu_id"].astype(str), strict=True))
    found = [
        (str(d), frozenset(tu_of[str(g)] for g in gs))
        for d, gs in zip(out["direction"], out["genes"], strict=True)
    ]
    return found, planted


def run_cell(
    n_tus: int,
    effect_sd: float,
    direction: str,
    R: int,
    seed: int,
    *,
    detection: str = "covariates_only",
) -> dict[str, Any]:
    """Power, mean Jaccard over hits and false-island rate of one cell over ``R`` replicates."""
    pair = scenario(n_tus, effect_sd, direction, detection=detection)
    hits, jacs, falses = 0, [], 0
    for r in range(R):
        found, planted = _replicate(for_step(f"{STREAM}.r{r}", seed), pair)
        hit, jac, false = replicate_metrics(found, planted, direction)
        hits += hit
        falses += false
        if hit:
            jacs.append(jac)
    return {
        "n_tus": n_tus,
        "effect_sd": effect_sd,
        "direction": direction,
        "detection": detection,
        "R": R,
        "seed": seed,
        "power": hits / R,
        "mean_jaccard": sum(jacs) / len(jacs) if jacs else math.nan,
        "false_island_rate": falses / R,
    }


def _commit() -> str:
    try:
        res = subprocess.run(["git", "-C", str(REPO), "rev-parse", "HEAD"],
                             capture_output=True, text=True, check=True)  # fmt: skip
        return res.stdout.strip() or "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0] if __doc__ else None)
    ap.add_argument("--replicates", type=int, default=R_DEFAULT)
    ap.add_argument("--cells", default=None, help="comma-separated indices into CELLS (all)")
    ap.add_argument("--detection", default="covariates_only",
                    choices=("covariates_only", "abundance"))  # fmt: skip
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--run-id", default="")
    args = ap.parse_args(argv)
    idx = range(len(CELLS)) if args.cells is None else [int(x) for x in args.cells.split(",")]
    commit = _commit()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    for i in idx:
        n, e, d = CELLS[i]
        rec = run_cell(n, e, d, args.replicates, cell_seed(i), detection=args.detection)
        rec.update(cell=i, commit=commit, run_id=args.run_id)
        line = json.dumps({k: (None if isinstance(v, float) and math.isnan(v) else v)
                           for k, v in rec.items()})  # fmt: skip
        with args.out.open("a") as fh:
            fh.write(line + "\n")
        print(line, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
