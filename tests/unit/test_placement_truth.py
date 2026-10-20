"""The M3-21 truth oracle (tests/helpers/placement_truth.py, placement_sims.py): PR #113 review B1
and the overlap rule of lead decision 2."""

from __future__ import annotations

from pathlib import Path

import pytest
from helpers.placement_sims import dbstrain_input
from helpers.placement_truth import Truth


@pytest.fixture(scope="module")
def dbstrain(tmp_path_factory: pytest.TempPathFactory):  # type: ignore[no-untyped-def]
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SimConfig, simulate
    from giapsim.strain import StrainConfig

    tmp = tmp_path_factory.mktemp("dbstrain")
    cfg = SimConfig(
        genome_length=100_000,
        strain=StrainConfig(snv_rate_per_kb=5.0),
        proteomics=ProteomicsConfig(pi_wrong=0.01, search_db="other", n_runs=1),
    )
    res = simulate(cfg, 20261021, tmp / "sim")
    return dbstrain_input(res, tmp)


def test_every_database_strain_peptide_has_an_origin(dbstrain) -> None:  # type: ignore[no-untyped-def]
    """B1: each database peptide that replaced a sample peptide has the sample peptide's origin,
    including pairs equal under I = L (an L/I SNV)."""
    _path, truth, back = dbstrain
    assert back
    missing = [p for p in back if not truth.origins.get(p)]
    assert not missing, f"database peptides without an origin: {missing[:10]}"


def test_absent_only_when_different_under_i_equals_l(dbstrain) -> None:  # type: ignore[no-untyped-def]
    _path, truth, back = dbstrain
    il = {
        p for p, ss in back.items() if all(p.replace("I", "L") == s.replace("I", "L") for s in ss)
    }
    assert not (truth.absent & il)
    differ = {p for p, ss in back.items() if any(p.replace("I", "L") != s.replace("I", "L")
                                                 for s in ss)}  # fmt: skip
    assert truth.absent == differ
    assert differ, "a 5-SNV/kb strain must change some detected peptides"


def _truth() -> Truth:
    t = Truth({"c": 1000}, {"c": True})
    t.add("P", [("c", 100, 130, "+")], ["G1"])
    t.add("M", [("c", 400, 430, "-")], ["G2"])
    t.add("O", [("c", 990, 1020, "+")], ["G3"])  # across the origin
    return t


@pytest.mark.parametrize(
    ("hit", "pep", "ok"),
    [
        (("c", 100, 130, "+"), "P", True),
        (("c", 115, 145, "+"), "P", True),  # exactly 50 %, same phase
        (("c", 118, 148, "+"), "P", False),  # < 50 %
        (("c", 101, 131, "+"), "P", False),  # other phase
        (("c", 100, 130, "-"), "P", False),  # other strand
        (("c", 397, 430, "-"), "M", True),  # minus: phase from the end
        (("c", 400, 431, "-"), "M", False),
        (("c", 990, 1020, "+"), "O", True),
        (("c", -10, 20, "+"), "O", True),  # wrapped coordinates
        (("c", 100, 130, "+"), "X", False),  # no origin
    ],
)
def test_overlap_rule(hit: tuple[str, int, int, str], pep: str, ok: bool) -> None:
    assert _truth().overlaps(hit, pep) is ok


def test_dbstrain_path_written(dbstrain) -> None:  # type: ignore[no-untyped-def]
    path, _t, _b = dbstrain
    assert Path(path).is_file()
