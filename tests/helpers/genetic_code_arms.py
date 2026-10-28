"""Genetic-code arms of the existing calibration suites (M11-16f; docs/specs/
M11-16-translation-tables-4-25.md item 9, decision D14).

Arm 11 stays in each suite's own file, unchanged and first. Arms 4 and 25 reuse that suite's
tolerance keys, replicate counts and seeds (R3: ``tolerances.yaml`` and the seed lists are not
touched):

- ``sim_at(t)`` runs every ``gippasim.simulate`` call of the block at genetic code ``t`` and
  writes ``annotation.translation_table: t`` into the simulator's ``config.yaml`` (gippasim does
  not, see ``stamp_config``), so helpers such as ``placement_sims.run_chain`` need no change;
- ``arm_seeds`` maps the replicate seeds the harness draws for an arm's own stream name (its
  records stay apart from arm 11's) to the seeds arm 11 draws for the same replicate index.
"""

from __future__ import annotations

import dataclasses
import sys
from collections.abc import Iterator
from contextlib import contextmanager
from pathlib import Path
from typing import Any

import pytest
import yaml

#: arms in the fast tier and in ``calibration_full`` (item 9 table, D14)
RECALL_FAST = (4,)
RECALL_FULL = (4, 25)
SENSE_TABLES = (4, 25)
#: TGA's residue at each table (NCBI)
TGA_AA = {4: "W", 25: "G"}


def stamp_config(config_yaml: Path | str | None, table: int) -> None:
    """Set ``annotation.translation_table`` in the simulator's run config (no-op at 11)."""
    if config_yaml is None or table == 11:
        return
    path = Path(config_yaml)
    data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    data.setdefault("annotation", {})["translation_table"] = table
    path.write_text(yaml.safe_dump(data, sort_keys=False), encoding="utf-8")


@contextmanager
def sim_at(table: int) -> Iterator[None]:
    """Every ``gippasim.simulate`` call inside the block simulates at genetic code ``table``."""
    import gippasim.simulate  # noqa: F401  (the module; the package attribute is the function)

    mod = sys.modules["gippasim.simulate"]
    orig = mod.simulate

    def at_table(config: Any, seed: int, outdir: Any) -> Any:
        res = orig(dataclasses.replace(config, genetic_code=table), seed, outdir)
        stamp_config(res.config_yaml, table)
        return res

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(mod, "simulate", at_table)
        yield


def arm_seeds(base: str, arm_base: str, n: int) -> dict[int, int]:
    """Replicate seed of stream ``arm_base`` -> replicate seed of stream ``base``, index < ``n``
    (the draw of ``_placement_suite.run_fdr_suite`` and ``run_balance``)."""
    import _criteria as C

    def draw(stream: str, i: int) -> int:
        return int(C.replicate_rng(stream, i).integers(0, 2**31 - 1))

    out = {draw(arm_base, i): draw(base, i) for i in range(n)}
    assert len(out) == n, "replicate seed collision between indices"
    return out


def in_frame_tga(nt: str) -> int:
    """Number of in-frame TGA codons of a coding stretch."""
    return sum(nt[i : i + 3] == "TGA" for i in range(0, len(nt) - 2, 3))
