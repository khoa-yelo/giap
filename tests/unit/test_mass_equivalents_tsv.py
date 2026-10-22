"""M7-05 the committed ``mass_equivalents.tsv`` and its generator (docs/specs/M7-05-variant-flags.md
rule 1, lead decision 2, acceptance criterion 1; GIAP_SPEC.md §7.5.3, Appendix B.5).

Interface pinned here (test-author reading of the spec):

* ``python -m giap.data.build_mass_equivalents`` writes ``src/giap/data/mass_equivalents.tsv``
  (the package data file); ``--out PATH`` writes the same bytes to ``PATH`` instead, which is how
  this test regenerates it without touching the committed file.
* Format: ``#`` comment lines first (at least one, naming the generator
  ``build_mass_equivalents`` and Appendix ``B.5``), then a header row whose first four columns are
  ``left``, ``right``, ``delta_da``, ``kind`` (the ``MassEquivalence`` fields), then one row per
  entry of ``giap.mapping.masses.mass_equivalents(0.05)`` in that function's order. ``delta_da``
  is a decimal that parses to the entry's value within 5e-7 Da.
* The run time never reads the file (lead decision 2); ``test_flags.py`` checks the flag follows
  the configured tolerance.
"""

from __future__ import annotations

import os
import subprocess
import sys
from importlib.resources import files
from pathlib import Path

import pytest

from giap.mapping.masses import mass_equivalents

MODULE = "giap.data.build_mass_equivalents"


def _committed() -> str:
    res = files("giap.data").joinpath("mass_equivalents.tsv")
    assert res.is_file(), "src/giap/data/mass_equivalents.tsv is not committed"
    return res.read_text(encoding="utf-8")


def _rows(text: str) -> list[dict[str, str]]:
    lines = [ln for ln in text.splitlines() if ln and not ln.startswith("#")]
    head = lines[0].split("\t")
    assert head[:4] == ["left", "right", "delta_da", "kind"], head
    return [dict(zip(head, ln.split("\t"), strict=True)) for ln in lines[1:]]


@pytest.fixture(scope="module")
def rows() -> list[dict[str, str]]:
    return _rows(_committed())


def _multiset(side: str) -> str:
    return side if "+" in side else "".join(sorted(side))


def _has(rows: list[dict[str, str]], a: str, b: str, delta: float | None = None) -> bool:
    want = {_multiset(a), _multiset(b)}
    for r in rows:
        if {_multiset(r["left"]), _multiset(r["right"])} == want:
            # B.5 prints rounded deltas (M + Oxidation vs F is 0.033005, printed 0.0329)
            return delta is None or abs(float(r["delta_da"]) - delta) < 2e-4
    return False


def test_generator_regenerates_the_committed_file(tmp_path: Path) -> None:
    out = tmp_path / "mass_equivalents.tsv"
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join(p for p in sys.path if p)
    proc = subprocess.run(
        [sys.executable, "-m", MODULE, "--out", str(out)],
        capture_output=True, text=True, env=env, check=False,
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    assert out.read_text(encoding="utf-8") == _committed()


def test_header_comment_names_the_generator_and_b5() -> None:
    comments = "\n".join(ln for ln in _committed().splitlines() if ln.startswith("#"))
    assert "build_mass_equivalents" in comments
    assert "B.5" in comments


def test_tsv_equals_the_in_process_list_at_the_default_tolerance(
    rows: list[dict[str, str]],
) -> None:
    """Lead decision 2: the committed list and ``mass_equivalents(0.05)`` agree row for row."""
    want = mass_equivalents(0.05)
    assert len(rows) == len(want)
    for r, e in zip(rows, want, strict=True):
        assert (r["left"], r["right"], r["kind"]) == (e.left, e.right, e.kind)
        assert abs(float(r["delta_da"]) - e.delta_da) <= 5e-7, (r, e)


# The twelve checks of Appendix B.5 "the generated list MUST pass" at 0.05 Da.
B5_CHECKS: list[tuple[str, str, float]] = [
    ("N", "GG", 0.0),
    ("Q", "AG", 0.0),
    ("GT", "AS", 0.0),
    ("AV", "GL", 0.0),
    ("AN", "GQ", 0.0),
    ("W", "AD", 0.0153),
    ("W", "EG", 0.0153),
    ("R", "GV", 0.0112),
    ("K", "Q", 0.0364),
    ("A+Oxidation", "S", 0.0),
    ("N+Deamidated", "D", 0.0),
    ("G+Methyl", "A", 0.0),
    ("M+Oxidation", "F", 0.0329),
]


@pytest.mark.parametrize(("a", "b", "delta"), B5_CHECKS, ids=[f"{a}~{b}" for a, b, _ in B5_CHECKS])
def test_b5_checks_pass_on_the_committed_file(
    rows: list[dict[str, str]], a: str, b: str, delta: float
) -> None:
    assert _has(rows, a, b, delta), f"{a} vs {b} (Δ {delta}) missing from mass_equivalents.tsv"


def test_every_row_is_within_the_default_tolerance(rows: list[dict[str, str]]) -> None:
    assert rows
    for r in rows:
        assert 0.0 <= float(r["delta_da"]) <= 0.05, r
        assert r["kind"] in {"block", "near_block", "substitution", "modification"}, r


def test_no_entry_for_g_phospho_vs_h(rows: list[dict[str, str]]) -> None:
    """Δ 0.0711: outside the default list (the 0.1 Da flag in test_flags.py is recomputed)."""
    assert not _has(rows, "G+Phospho", "H")
