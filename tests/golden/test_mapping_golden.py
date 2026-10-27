"""Golden mapping tables for the 50-kb fixture (issue M3-17 acceptance; GIPPA_SPEC.md §7.2.9
"Golden: the hit tables for the 50-kb fixture", §11.4; M3-21 "golden" row).

The fixture: ``simulate(SimConfig(proteomics=ProteomicsConfig(), write_precomputed=True),
GOLDEN_SEED)`` (50 kb, M3-19 defaults), ``gippa preflight`` with genome, peptides and Bakta-like
JSON, then the mapping steps through ``run_step`` (helpers.mapping_steps.drive) with
``seed = GOLDEN_SEED``. The four mapping TSVs must equal ``tests/data/golden/sim50k/mapping/``
byte for byte.

Regenerating (only on purpose, reviewed in the PR diff): ``GIPPA_UPDATE_GOLDEN=1 pytest
tests/golden/test_mapping_golden.py`` writes the files and skips.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path

import pytest
from helpers.mapping_steps import (
    MAPPING_TSVS,
    drive,
    fresh_steps,
    module2_files,
    preflight,
    roots_of,
)

from gippa.config import resolve_config

GOLDEN = Path(__file__).resolve().parents[1] / "data" / "golden" / "sim50k" / "mapping"


def test_mapping_tsvs_equal_the_committed_golden_files(tmp_path: Path, golden_seed: int) -> None:
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SimConfig, simulate

    res = simulate(
        SimConfig(proteomics=ProteomicsConfig(), write_precomputed=True),
        golden_seed,
        tmp_path / "sim",
    )
    run = preflight(res, tmp_path / "pf")
    cfg = resolve_config(run.config_yaml, [f"seed={golden_seed}"])
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, tmp_path / "work", roots_of(run))
    files = module2_files(tmp_path / "work", avail)
    if os.environ.get("GIPPA_UPDATE_GOLDEN") == "1":
        GOLDEN.mkdir(parents=True, exist_ok=True)
        for f in MAPPING_TSVS:
            shutil.copy(files[f], GOLDEN / f)
        pytest.skip(f"golden files written to {GOLDEN}")
    missing = [f for f in MAPPING_TSVS if not (GOLDEN / f).is_file()]
    assert not missing, f"golden files missing in {GOLDEN}: {missing}"
    for f in MAPPING_TSVS:
        got, want = files[f].read_text(), (GOLDEN / f).read_text()
        if got != want:
            diff = [
                (i + 1, a, b)
                for i, (a, b) in enumerate(zip(want.splitlines(), got.splitlines(), strict=False))
                if a != b
            ][:5]
            pytest.fail(f"{f} differs from the golden file; first lines (n, want, got): {diff}")
        assert files[f].read_bytes() == (GOLDEN / f).read_bytes(), f
