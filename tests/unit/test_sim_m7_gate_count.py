"""giapsim M7-12 lead decision 26: the gate count of ``novel_orf_power`` (class q floor; PR #225,
statistician round 3; GIAP_SPEC.md §7.6.4, §7.6.6, §8.4; docs/methods/target-decoy.md).

The picked class q ``(D + 1) / T`` cannot fall below ``1 / T``, so a class with fewer than
``1 / alpha`` = 100 target ORFs (``stats.alpha.novel_orfs`` = 0.01) calls nothing. M7-13 G3/G4
therefore plant ``NOVEL_PER_CLASS_GATE`` ORFs per class at 1 Mb; the default sizes keep
``NOVEL_PER_CLASS`` and every existing byte pin.

Interface choices pinned here:

* ``giapsim.scenarios.novel_orf_power(seed, size="100k", n_per_class: int | None = None)``;
  ``None`` means ``NOVEL_PER_CLASS[size]`` and gives a config and output byte-identical to the
  call without the argument; an int replaces the per-class count of all four classes
  (``n_novel_orfs`` and ``n_nterm_extensions``) and nothing else.
* ``giapsim.scenarios.NOVEL_PER_CLASS_GATE == 300`` (exported in ``__all__``).
* **Shortfall report.** With ``m7_truth`` and novel ORFs requested, the simulator writes
  ``truth/novel_orfs.json``: ``{"requested": {class: int}, "realised": {class: int},
  "shortfall": {class: int}}`` over the four classes ``intergenic``, ``cds_outframe``,
  ``antisense``, ``nterm_extension``; ``realised`` equals the class's row count in
  ``truth/novel_orfs.tsv`` and ``shortfall = requested - realised`` (>= 0). Intergenic planting
  stops at the number of usable gaps (decision 15), about 100-145 at 1 Mb, so at the gate count
  its shortfall is positive and reported (decision 26: accepted and documented, not forced).
* At 1 Mb and the gate count, ORFs with ``n_true_peptides >= 1`` number >= 100 in
  ``cds_outframe``, ``antisense`` and ``nterm_extension``, and >= 90 in ``intergenic`` (seed
  ``SEED``; the intergenic count is seed-dependent: seed 1 gives 88 of 104 planted).

The 1-Mb cases take about 8 s each and are ``slow``.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Any

import pytest
from helpers import m7_sim_oracle as O

import giapsim

SEED = 20261022
ALPHA_NOVEL = 0.01  # stats.alpha.novel_orfs (src/giap/config/defaults.yaml)
CLASSES4 = ("intergenic", "cds_outframe", "antisense", "nterm_extension")


def sc() -> Any:
    import giapsim.scenarios

    return giapsim.scenarios


def _digest(out: Path) -> dict[str, str]:
    return {
        str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(out.rglob("*"))
        if p.is_file()
    }


def _requested(cfg: Any) -> dict[str, int]:
    return {**dict(cfg.n_novel_orfs), "nterm_extension": cfg.n_nterm_extensions}


# ---------------------------------------------------------------- the constant and the argument


def test_gate_constant_is_300() -> None:
    assert sc().NOVEL_PER_CLASS_GATE == 300
    assert "NOVEL_PER_CLASS_GATE" in sc().__all__


def test_gate_count_reaches_one_over_alpha() -> None:
    assert sc().NOVEL_PER_CLASS_GATE * ALPHA_NOVEL >= 1


@pytest.mark.parametrize("size", ["100k", "300k", "1M"])
def test_none_config_equals_default(size: str) -> None:
    default = sc().novel_orf_power(SEED, size)
    assert sc().novel_orf_power(SEED, size, n_per_class=None) == default
    explicit = sc().novel_orf_power(SEED, size, n_per_class=sc().NOVEL_PER_CLASS[size])
    assert explicit == default


@pytest.mark.parametrize("n", [1, 5, 300])
def test_n_per_class_sets_all_four_classes_only(n: int) -> None:
    default = sc().novel_orf_power(SEED, "100k")
    cfg = sc().novel_orf_power(SEED, "100k", n_per_class=n)
    assert _requested(cfg) == dict.fromkeys(CLASSES4, n)
    same = dataclasses.replace(
        cfg, n_novel_orfs=default.n_novel_orfs, n_nterm_extensions=default.n_nterm_extensions
    )
    assert same == default


def test_none_output_byte_identical_at_100k(tmp_path: Path) -> None:
    a, b = tmp_path / "default", tmp_path / "none"
    giapsim.simulate(sc().novel_orf_power(SEED, "100k"), SEED, a)
    giapsim.simulate(sc().novel_orf_power(SEED, "100k", n_per_class=None), SEED, b)
    da, db = _digest(a), _digest(b)
    assert da == db
    assert "truth/novel_orfs.tsv" in da


# ---------------------------------------------------------------- shortfall report (100 kb)


def _report(out: Path) -> dict[str, dict[str, int]]:
    rep = json.loads((out / "truth" / "novel_orfs.json").read_text())
    rows = Counter(r["class"] for r in O.rows(out / "truth" / "novel_orfs.tsv"))
    for key in ("requested", "realised", "shortfall"):
        assert set(rep[key]) == set(CLASSES4), (key, rep[key])
    for c in CLASSES4:
        assert rep["realised"][c] == rows[c], (c, rep, rows)
        assert rep["shortfall"][c] == rep["requested"][c] - rep["realised"][c] >= 0, (c, rep)
    return rep


def test_report_written_at_default_100k(tmp_path: Path) -> None:
    cfg = sc().novel_orf_power(SEED, "100k")
    giapsim.simulate(cfg, SEED, tmp_path)
    rep = _report(tmp_path)
    assert rep["requested"] == dict.fromkeys(CLASSES4, sc().NOVEL_PER_CLASS["100k"])


def test_report_records_intergenic_shortfall_at_100k(tmp_path: Path) -> None:
    """100 kb has about 10 usable gaps (seed ``SEED``), so 20 per class shows the shortfall
    without a 1-Mb run (300 per class does not fit in 100 kb)."""
    giapsim.simulate(sc().novel_orf_power(SEED, "100k", n_per_class=20), SEED, tmp_path)
    rep = _report(tmp_path)
    assert rep["requested"] == dict.fromkeys(CLASSES4, 20)
    assert rep["shortfall"]["intergenic"] > 0, rep
    for c in ("cds_outframe", "antisense", "nterm_extension"):
        assert rep["shortfall"][c] == 0, rep


# ---------------------------------------------------------------- the gate count at 1 Mb


@pytest.fixture(scope="module")
def gate_1m(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("m7gate1M")
    cfg = sc().novel_orf_power(SEED, "1M", n_per_class=sc().NOVEL_PER_CLASS_GATE)
    giapsim.simulate(cfg, SEED, out)
    return out


def _with_peptide(out: Path) -> Counter[str]:
    n: Counter[str] = Counter()
    for r in O.rows(out / "truth" / "novel_orfs.tsv"):
        n[r["class"]] += int(r["n_true_peptides"]) >= 1
    return n


@pytest.mark.slow
@pytest.mark.parametrize("cls", ["cds_outframe", "antisense", "nterm_extension"])
def test_gate_count_1m_class_reaches_one_over_alpha(gate_1m: Path, cls: str) -> None:
    n = _with_peptide(gate_1m)
    assert n[cls] >= round(1 / ALPHA_NOVEL), (cls, dict(n))


@pytest.mark.slow
def test_gate_count_1m_intergenic_at_least_90_and_reported(gate_1m: Path) -> None:
    n = _with_peptide(gate_1m)
    assert n["intergenic"] >= 90, dict(n)
    rep = _report(gate_1m)
    assert rep["requested"] == dict.fromkeys(CLASSES4, sc().NOVEL_PER_CLASS_GATE)
    assert rep["shortfall"]["intergenic"] > 0, rep
    for c in ("cds_outframe", "antisense", "nterm_extension"):
        assert rep["shortfall"][c] == 0, rep
