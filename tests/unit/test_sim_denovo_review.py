"""gippasim de novo: the lead decisions after code review of PR #66 (issue M3-20).

docs/specs/M3-20-sim-denovo-errors.md, "Lead decisions after code review", items 1-2:

1. Redraw rule: a drawn error type without a site is redrawn among the remaining positive-weight
   types; only when none has a site does the read fall back to ``substitution`` (even at weight 0).
2. Mutant killers: wrong-ID spectra get a correct ``true_peptide`` that differs from the database
   peptide and varies across spectra; weights {block: 0.5, swap: 0.5} never give ``substitution``;
   every reversed-stretch length occurs; both residue orders of a block are read and written;
   de novo ``protein_ids`` are NA; the de novo RNG stream is its own (one pinned hash).

Convention pinned here: the decision text says stretch lengths "2-6", but the spec (§11.3 list in
the issue, "3-6 residues reversed") and ``reverse_stretch`` use 3-6; a 2-residue reversal is an
adjacent ``swap``. The test pins 3-6.
"""

from __future__ import annotations

import csv
import hashlib
from collections import Counter
from pathlib import Path
from typing import Any

import numpy as np
import pytest

import gippasim

SEED = 11
BASE = {"genome_length": 60_000}
#: sha256 of truth/denovo.tsv for SimConfig(genome_length=60_000,
#: proteomics=ProteomicsConfig(pi_wrong=0.05), denovo=DenovoConfig()), seed 11 (stream key 4).
DENOVO_TRUTH_SHA256 = "319e261f77cec3487bebedba788636d60aba9af55ea3f2922a77dec4f6762ad5"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def run_sim(out: Path, **dn_kw: Any) -> Path:
    from gippasim.denovo import DenovoConfig
    from gippasim.proteomics import ProteomicsConfig

    cfg = gippasim.SimConfig(
        proteomics=ProteomicsConfig(pi_wrong=0.05), denovo=DenovoConfig(**dn_kw), **BASE
    )
    gippasim.simulate(cfg, SEED, out)
    return out


@pytest.fixture(scope="module")
def default(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return run_sim(tmp_path_factory.mktemp("default"))


def wrong_map(out: Path) -> dict[str, bool]:
    w: dict[str, bool] = {}
    for t in rows(out / "truth" / "peptides.tsv"):
        if t["origin"] != "decoy":
            w[t["peptide"]] = w.get(t["peptide"], False) or t["is_wrong"].lower() == "true"
    return w


# ---------------------------------------------------------------- item 2: wrong-ID spectra


def test_wrong_id_true_peptide_is_another_correct_peptide(default: Path) -> None:
    """Kills M01 (wrong IDs keep the DB peptide as truth)."""
    wrong = wrong_map(default)
    db = {(r["run"], r["spectrum_id"]): r["peptide"] for r in rows(default / "peptides.tsv")}
    n = 0
    for t in rows(default / "truth" / "denovo.tsv"):
        pep = db[(t["run"], t["spectrum_id"])]
        if wrong[pep]:
            n += 1
            assert t["true_peptide"] != pep, t
            assert wrong.get(t["true_peptide"]) is False, t
    assert n >= 20, n


def test_wrong_id_true_peptides_vary(default: Path) -> None:
    """Kills M02 (every wrong-ID spectrum gets the same correct peptide)."""
    wrong = wrong_map(default)
    db = {(r["run"], r["spectrum_id"]): r["peptide"] for r in rows(default / "peptides.tsv")}
    trues = [
        t["true_peptide"]
        for t in rows(default / "truth" / "denovo.tsv")
        if wrong[db[(t["run"], t["spectrum_id"])]]
    ]
    assert len(trues) >= 20
    assert len(set(trues)) >= len(trues) // 2, (len(set(trues)), len(trues))


# ---------------------------------------------------------------- item 1: redraw rule


def test_block_swap_weights_never_give_substitution(tmp_path: Path) -> None:
    """Kills M05 (no redraw: a peptide without a block site falls back to substitution)."""
    out = run_sim(tmp_path, error_weights={"block": 0.5, "swap": 0.5})
    types = Counter(t["error_type"] for t in rows(out / "truth" / "denovo.tsv"))
    assert types["substitution"] == 0, types
    assert set(types) <= {"none", "block", "swap"}, types
    assert types["block"] > 0 and types["swap"] > 0, types


def test_no_site_among_positive_weights_falls_back_to_substitution(tmp_path: Path) -> None:
    out = run_sim(tmp_path, error_weights={"block": 1.0})
    types = Counter(t["error_type"] for t in rows(out / "truth" / "denovo.tsv"))
    assert set(types) <= {"none", "block", "substitution"}, types
    assert types["block"] > 0 and types["substitution"] > 0, types


# ---------------------------------------------------------------- item 2: operator coverage


def test_every_stretch_length_occurs() -> None:
    """Kills M07 (stretches capped at 5). Distinct residues: the changed span is the length."""
    from gippasim.denovo import apply_error

    pep = "ACDEFGHKMNPQRSTVWY"
    rng = np.random.default_rng(0)
    lengths = set()
    for _ in range(400):
        out = apply_error(pep, "reversed_stretch", rng)
        d = [i for i, (a, b) in enumerate(zip(pep, out, strict=True)) if a != b]
        lengths.add(d[-1] - d[0] + 1)
    assert lengths == {3, 4, 5, 6}, lengths


def test_block_read_in_reverse_order() -> None:
    """Kills M10 (only the listed order of a block is searched): GA (reverse of AG) reads as Q."""
    from gippasim.denovo import apply_error

    pep = "PEPKGAEKR"  # only block site: GA
    rng = np.random.default_rng(0)
    assert {apply_error(pep, "block", rng) for _ in range(20)} == {"PEPKQEKR"}


def test_block_written_in_both_orders() -> None:
    """Kills M09 (the new block is always written in one order): Q reads as AG and as GA."""
    from gippasim.denovo import apply_error

    pep = "PEPKQEKR"  # only block site: Q
    rng = np.random.default_rng(0)
    assert {apply_error(pep, "block", rng) for _ in range(200)} == {"PEPKAGEKR", "PEPKGAEKR"}


# ---------------------------------------------------------------- item 2: output columns, RNG


def test_denovo_protein_ids_are_na(default: Path) -> None:
    """Kills M21 (de novo rows keep the database protein_ids)."""
    dn = rows(default / "peptides_denovo.tsv")
    assert dn and "protein_ids" in dn[0]
    assert {r["protein_ids"] for r in dn} == {"NA"}


def test_denovo_stream_pinned(default: Path) -> None:
    """Kills M16 (de novo draws on the proteomics stream key 3)."""
    got = hashlib.sha256((default / "truth" / "denovo.tsv").read_bytes()).hexdigest()
    assert got == DENOVO_TRUTH_SHA256
