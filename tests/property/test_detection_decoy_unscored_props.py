"""Hypothesis properties of M4-24 (docs/specs/M4-24-decoy-unscored.md, acceptance criterion 4;
GIAP_SPEC.md §7.3.2 picked steps 1-3, §8 preamble / R4).

(a) ``decoy_unscored`` in ``warnings`` iff the pairing guard passes, the target side (``accepted``)
    is non-empty and S != C_T or S is empty;
(b) setting every value of a key column outside S to NA on both sides leaves ``tier``, ``table``
    and the ``q_threshold`` rows unchanged. Scope: S non-empty (with S empty, NA-ing the target q
    can move the no-decoy chain from ``mayu`` to ``two_peptide``, which rule 3 intends);
(c) when every paired decoy is complete in every column of C_T (C_T non-empty), there is no
    ``decoy_unscored`` and the output equals the output with the decoy-only columns set to NA;
(d) row-permutation invariance of tier, warnings and table with partly NA key columns.

Every column of each side is drawn ``full`` / ``none`` / ``partial`` so S and C_T take all their
shapes. Target and decoy q stay <= ``q_value_max`` (0.01) or NA, so every target is accepted and
the decoy q filter (rule 6) removes nothing: pairing = every decoy row.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import detection_tiers_fixtures as fx
from detection_tiers_fixtures import Builder, cfg

FAST = settings(max_examples=60, deadline=None, suppress_health_check=[HealthCheck.too_slow])
KEYS = ("q_value", "pep", "score")
VALUES = {
    "q_value": st.sampled_from([0.0, 0.0005, 0.001, 0.005, 0.01]),
    "pep": st.sampled_from([0.001, 0.01, 0.05, 0.2]),
    "score": st.sampled_from([10.0, 20.0, 30.0, 40.0]),
}
MODE = st.sampled_from(["full", "none", "partial"])
UNSCORED = "decoy_unscored"


@st.composite
def _rows(draw: st.DrawFn, n_genes: int, max_size: int, modes: dict[str, str]) -> list[Any]:
    rows = []
    for _ in range(draw(st.integers(0, max_size))):
        g = draw(st.integers(0, n_genes - 1))
        vals = []
        for k in KEYS:
            m = modes[k]
            keep = m == "full" or (m == "partial" and draw(st.booleans()))
            vals.append(draw(VALUES[k]) if keep else None)
        rows.append((g, *vals))
    return rows


@st.composite
def scenario(draw: st.DrawFn) -> dict[str, Any]:
    n = draw(st.integers(1, 6))
    tmodes = {k: draw(MODE) for k in KEYS}
    dmodes = {k: draw(MODE) for k in KEYS}
    return {
        "n": n,
        "partial": draw(st.sets(st.integers(0, n - 1), max_size=2)),
        "targets": draw(_rows(n, 10, tmodes)),
        "decoys": draw(_rows(n, 6, dmodes)),
        "unpaired": draw(st.integers(0, 2)),
        "alpha": draw(st.sampled_from([0.01, 0.1, 0.5, 1.0])),
    }


def _na(rows: list[Any], cols: set[str]) -> list[Any]:
    return [(r[0], *(None if k in cols else v for k, v in zip(KEYS, r[1:], strict=True)))
            for r in rows]  # fmt: skip


def _run(sc: dict[str, Any], shuffle: np.random.Generator | None = None) -> fx.Run:
    b = Builder(sc["n"], partial={fx.tag(i) for i in sc["partial"]})
    for g, q, pep, score in sc["targets"]:
        b.target(fx.tag(g), q, pep, score)
    for g, q, pep, score in sc["decoys"]:
        b.decoy([f"rev_{fx.tag(g)}"], q, pep, score)
    for i in range(sc["unpaired"]):
        b.decoy([f"rev_UNKNOWN{i}"], 0.001, 0.01, 10.0)
    c = cfg(stats={"alpha": {"identification": sc["alpha"]}})
    return b.run(c, shuffle=shuffle)


def _complete(rows: list[Any]) -> tuple[str, ...]:
    return tuple(k for j, k in enumerate(KEYS) if all(r[1 + j] is not None for r in rows))


def _expected(sc: dict[str, Any], r: fx.Run) -> tuple[bool, tuple[str, ...], tuple[str, ...]]:
    """(guard passes and the target side is non-empty, C_T, S)."""
    n_paired, n_all = len(sc["decoys"]), len(sc["decoys"]) + sc["unpaired"]
    guard = n_all > 0 and n_paired / n_all >= 0.8
    acc = r.accepted
    ct = tuple(k for k in KEYS if len(acc) and bool(acc[k].astype(float).notna().any()))
    comp = set(_complete(sc["decoys"]))
    s = tuple(k for k in ct if k in comp)
    return guard and len(acc) > 0, ct, s


def _same_output(a: fx.Run, b: fx.Run) -> None:
    assert a.det.tier == b.det.tier
    pd.testing.assert_frame_equal(a.det.table, b.det.table)
    pd.testing.assert_frame_equal(fx.curve_rows(a.det).reset_index(drop=True),
                                  fx.curve_rows(b.det).reset_index(drop=True))  # fmt: skip


@FAST
@given(sc=scenario())
def test_a_warning_iff_shared_key_differs_from_target_key(sc: dict[str, Any]) -> None:
    r = _run(sc)
    active, ct, s = _expected(sc, r)
    expect = active and (s != ct or not s)
    assert (UNSCORED in r.det.warnings) == expect, (ct, s, r.det.tier, r.det.warnings)
    if active and not s:
        assert r.det.tier != "decoy_picked", "S empty: the no-decoy chain (rule 3)"
    elif active:
        assert r.det.tier == "decoy_picked"


@FAST
@given(sc=scenario())
def test_b_columns_outside_s_do_not_matter(sc: dict[str, Any]) -> None:
    r = _run(sc)
    active, _ct, s = _expected(sc, r)
    if not active or not s:
        return
    out = set(KEYS) - set(s)
    na = {**sc, "targets": _na(sc["targets"], out), "decoys": _na(sc["decoys"], out)}
    _same_output(r, _run(na))


@FAST
@given(sc=scenario())
def test_c_complete_decoys_never_warn_and_decoy_only_columns_are_ignored(
    sc: dict[str, Any],
) -> None:
    r = _run(sc)
    active, ct, _s = _expected(sc, r)
    if not active or not ct or not set(ct) <= set(_complete(sc["decoys"])):
        return
    assert UNSCORED not in r.det.warnings
    assert r.det.tier == "decoy_picked"
    decoy_only = set(KEYS) - set(ct)
    _same_output(r, _run({**sc, "decoys": _na(sc["decoys"], decoy_only)}))


@FAST
@given(sc=scenario(), seed=st.integers(0, 2**32 - 1))
def test_d_row_permutation_invariant(sc: dict[str, Any], seed: int) -> None:
    a = _run(sc)
    b = _run(sc, shuffle=np.random.default_rng(seed))
    assert a.det.tier == b.det.tier
    assert sorted(a.det.warnings) == sorted(b.det.warnings)
    pd.testing.assert_frame_equal(a.det.table, b.det.table)
