"""Properties of gippasim's de novo error model (issue M3-20; GIPPA_SPEC.md §11.3, §7.2.4).

Spec acceptance criterion: every error type, applied to a peptide, changes the sequence (beyond
I/L) and keeps it within 7-40 aa; ``block``/``near_block`` errors change the monoisotopic mass by
at most 0.05 Da. Also: the result is explained by the independent oracle
(tests/helpers/denovo_oracle.py) for the claimed type, ``redraw_il`` only touches I/L, the
operators are deterministic in the generator, and the simulated truth labels agree with the
oracle for random seeds and error weights.

Conventions as in tests/unit/test_sim_denovo.py.
"""

from __future__ import annotations

import csv
import tempfile
from pathlib import Path

import numpy as np
from helpers.denovo_oracle import AA20, ERROR_TYPES, explains, il, mass
from hypothesis import given, settings
from hypothesis import strategies as st

import gippasim

#: Contains a site for every error type in both length directions (N/GG, W/AD, R, A/N/G).
MOTIF = "NGGWRAD"
peptides = st.text(alphabet=st.sampled_from(AA20), min_size=7, max_size=40)


@st.composite
def with_sites(draw: st.DrawFn) -> str:
    n = draw(st.integers(7, 40))
    body = draw(st.text(alphabet=st.sampled_from(AA20), min_size=n, max_size=n))
    at = draw(st.integers(0, n - len(MOTIF)))
    return body[:at] + MOTIF + body[at + len(MOTIF) :]


def check(pep: str, out: str, et: str) -> None:
    assert isinstance(out, str)
    assert il(out) != il(pep), (pep, out, et)
    assert 7 <= len(out) <= 40, (pep, out, et)
    assert set(out) <= set(AA20), out
    if et in ("block", "near_block"):
        assert abs(mass(out) - mass(pep)) <= 0.05, (pep, out, et)
    if et == "random":
        assert 7 <= len(out) <= 20, out
    assert explains(pep, out, et), (pep, out, et)


@settings(max_examples=300, deadline=None)
@given(pep=with_sites(), et=st.sampled_from(ERROR_TYPES), seed=st.integers(0, 2**32 - 1))
def test_every_error_type_changes_peptide_within_bounds(pep: str, et: str, seed: int) -> None:
    from gippasim.denovo import apply_error

    check(pep, apply_error(pep, et, np.random.default_rng(seed)), et)


@settings(max_examples=300, deadline=None)
@given(pep=peptides, et=st.sampled_from(ERROR_TYPES), seed=st.integers(0, 2**32 - 1))
def test_any_peptide_error_valid_or_rejected(pep: str, et: str, seed: int) -> None:
    from gippasim.denovo import apply_error

    try:
        out = apply_error(pep, et, np.random.default_rng(seed))
    except ValueError:
        # only allowed when the type has no applicable site (never for these two)
        assert et not in ("substitution", "random"), (pep, et)
        return
    check(pep, out, et)


@settings(max_examples=200, deadline=None)
@given(pep=with_sites(), et=st.sampled_from(ERROR_TYPES), seed=st.integers(0, 2**32 - 1))
def test_apply_error_deterministic(pep: str, et: str, seed: int) -> None:
    from gippasim.denovo import apply_error

    a = apply_error(pep, et, np.random.default_rng(seed))
    b = apply_error(pep, et, np.random.default_rng(seed))
    assert a == b


@settings(max_examples=300, deadline=None)
@given(pep=peptides, seed=st.integers(0, 2**32 - 1))
def test_redraw_il_only_touches_il(pep: str, seed: int) -> None:
    from gippasim.denovo import redraw_il

    out = redraw_il(pep, np.random.default_rng(seed))
    assert len(out) == len(pep) and il(out) == il(pep)
    for a, b in zip(pep, out, strict=True):
        assert a == b or {a, b} == {"I", "L"}


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


@settings(max_examples=6, deadline=None)
@given(
    seed=st.integers(0, 10_000),
    weights=st.lists(st.floats(0.0, 1.0), min_size=6, max_size=6),
    w_sub=st.floats(0.05, 1.0),
)
def test_simulated_truth_agrees_with_oracle(seed: int, weights: list[float], w_sub: float) -> None:
    # substitution (always applicable) keeps a positive weight, so fallbacks stay in-weight
    weights = [w_sub, *weights]
    from gippasim.denovo import DenovoConfig
    from gippasim.proteomics import ProteomicsConfig

    cfg = gippasim.SimConfig(
        genome_length=30_000,
        proteomics=ProteomicsConfig(pi_wrong=0.05, n_runs=1),
        denovo=DenovoConfig(error_weights=dict(zip(ERROR_TYPES, weights, strict=True))),
    )
    with tempfile.TemporaryDirectory() as d:
        out = Path(d)
        gippasim.simulate(cfg, seed, out)
        truth = _rows(out / "truth" / "denovo.tsv")
        denovo = _rows(out / "peptides_denovo.tsv")
    assert len(truth) == len(denovo) > 0
    for t in truth:
        et = t["error_type"]
        assert et == "none" or weights[ERROR_TYPES.index(et)] > 0, t
        assert explains(t["true_peptide"], t["peptide"], et), t
        assert (t["correct"].lower() == "true") == (il(t["peptide"]) == il(t["true_peptide"]))
        assert 7 <= len(t["peptide"]) <= 40
