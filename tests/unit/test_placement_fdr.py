"""Placement FDR (G3): strata, pi_wrong, q-values, ``mapping_fdr`` rows (issue M3-10; GIAP_SPEC.md
§7.2.5 steps 1-6, §8.1, §9.3, Appendix B.6).

Conventions pinned here (docs/specs/M3-10-placement-fdr.md, lead decisions in M3-index.md):
- ``stratum_qvalues(e_target, e_decoy, pi_decoy, pi_plus)`` (ADR-0010, supersedes M3-10 Q2):
  T(e) = #targets with E <= e, W(e) = sum of the decoy weights ``pi_decoy`` over decoys with
  E <= e (a decoy tied with a target counts at that threshold), c = ``pi_plus``,
  FDR(e) = (W(e) + c) / T(e) evaluated at every target E, q_i = min over e >= E_i, clipped to
  <= 1; output in input order. Constant pi: weights pi, c = pi -> pi (D(e) + 1) / T(e) (§7.2.5).
- ``placement``: each decoy hit carries its own peptide's pi; c = max pi over the searched
  peptides of the decoy type (all its rows except ``same_db`` exact peptides in stratum 0; 1 if
  none; m >= 1 -> 1).
- pi_wrong: database m = 0 -> ``pi_wrong_database`` (mean PEP over peptides with PEP, else max
  q <= ``input.q_value_max``, else 1); de novo m = 0 -> the peptide's ``pi_denovo`` (1 if missing);
  m >= 1 -> 1.
- Lists are per decoy type: database peptides vs reverse decoys, de novo vs shuffle decoys.
- ``search_db = same`` with ``m_target = 0``: placement_q = 0, flag ``same_db``, left out of the
  stratum-0 target and decoy lists.
- No target hit: no placement_q (absent from the per-peptide table, or NaN and not accepted).
- ``mapping_fdr`` rows: decoy_type x n_diff x length class (labels ``7``, ``8``, ``9``, ``10-11``,
  ``12-14``, ``15-19``, ``20-40``) plus ``all``; ``decoys`` = D, ``fdr_stratum`` = (W + c)/T
  (W over the class's decoy hits, c of the stratum; NA when T = 0), ``pi_wrong`` = (W + c)/(D + 1);
  ``section = placement``.
- Accept when q <= ``stats.alpha.placement`` (inclusive).
- ``low_complexity`` placement flag (lead decision 11, M3-03 statistician review): set on a
  dipeptide repeat / homopolymer peptide, never on an ordinary peptide.
"""

from __future__ import annotations

import dataclasses
import math
from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
import pandas as pd
import pytest

from giap.config import GiapConfig
from giap.mapping.chance import ResidueFreqs, expected_hits
from giap.mapping.hits import hits_frame
from giap.mapping.placement_fdr import (
    PeptideCosts,
    best_costs,
    pi_wrong_database,
    placement,
    stratum_qvalues,
)
from giap.stats.tda import tda_qvalues

ALPHABET19 = "ACDEFGHKLMNPQRSTVWY"
PREFIX = "MDEFGHKLNPQRSTVWYACDEFGHKLMNPQRSTVWY"
CLASS_LABELS = ["7", "8", "9", "10-11", "12-14", "15-19", "20-40"]
FDR_COLUMNS = {
    "length_class",
    "n_diff",
    "targets",
    "decoys",
    "pi_wrong",
    "fdr_stratum",
    "n_accepted",
    "decoy_type",
}


# ------------------------------------------------------------------------------------ helpers


def pep_of(i: int, length: int = 9) -> str:
    """Distinct, ordinary (>= 4 distinct residues) peptide number ``i`` of ``length`` aa."""
    assert length >= 7 and i < 19**3
    digits = "".join(ALPHABET19[(i // 19**k) % 19] for k in (2, 1, 0))
    return PREFIX[: length - 3] + digits


def costs_frame(rows: Iterable[dict[str, Any]]) -> pd.DataFrame:
    """Costs table built through ``PeptideCosts`` (defaults: database, search_db other)."""
    out = []
    for r in rows:
        r = dict(r)
        r.setdefault("length", len(r["peptide"]))
        r.setdefault("source", "database")
        r.setdefault("search_db", "other")
        for k in ("m_target", "m_decoy", "e_target", "e_decoy", "pep", "q_value", "pi_denovo"):
            r.setdefault(k, None)
        out.append(dataclasses.asdict(PeptideCosts(**r)))
    return pd.DataFrame(out, columns=[f.name for f in dataclasses.fields(PeptideCosts)])


def targets(
    es: Sequence[float], m: int, *, start: int = 0, length: int = 9, **kw: Any
) -> list[dict[str, Any]]:
    return [
        {"peptide": pep_of(start + i, length), "m_target": m, "e_target": e, **kw}
        for i, e in enumerate(es)
    ]


def decoys(
    es: Sequence[float], m: int, *, start: int = 5000, length: int = 9, **kw: Any
) -> list[dict[str, Any]]:
    """Peptides with a decoy hit only (no target hit)."""
    return [
        {"peptide": pep_of(start + i, length), "m_decoy": m, "e_decoy": e, **kw}
        for i, e in enumerate(es)
    ]


def cfg(alpha: float = 0.01, q_value_max: float = 0.01) -> GiapConfig:
    return GiapConfig.model_validate(
        {"stats": {"alpha": {"placement": alpha}}, "input": {"q_value_max": q_value_max}}
    )


def flags_of(x: object) -> set[str]:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return set()
    if isinstance(x, str):
        return {t for t in x.replace(",", "|").replace(";", "|").split("|") if t}
    return {str(t) for t in x}  # type: ignore[attr-defined]


def run(rows: Iterable[dict[str, Any]], **kw: Any) -> tuple[pd.DataFrame, pd.DataFrame]:
    per, fdr = placement(costs_frame(rows), cfg(**kw))
    return per, fdr


def by_peptide(per: pd.DataFrame) -> dict[str, dict[str, Any]]:
    assert per["peptide"].is_unique
    return {r["peptide"]: r for r in per.to_dict("records")}


def q_of(per: pd.DataFrame, peptides: Sequence[str]) -> np.ndarray:
    d = by_peptide(per)
    return np.array([float(d[p]["placement_q"]) for p in peptides])


def fdr_row(fdr: pd.DataFrame, m: int, cls: str, decoy_type: str = "reverse") -> dict[str, Any]:
    sel = fdr[
        (fdr["decoy_type"] == decoy_type)
        & (fdr["n_diff"].astype(int) == m)
        & (fdr["length_class"].astype(str) == cls)
    ]
    assert len(sel) == 1, (m, cls, decoy_type, fdr)
    row: dict[str, Any] = sel.iloc[0].to_dict()
    return row


def unplaced(per: pd.DataFrame, peptide: str) -> bool:
    d = by_peptide(per)
    if peptide not in d:
        return True
    r = d[peptide]
    return bool(pd.isna(r["placement_q"])) and not bool(r["accepted"])


# ------------------------------------------------------------------ stratum_qvalues (examples)


def test_example1_b6_150_targets_before_first_decoy() -> None:
    et = np.linspace(1e-6, 1e-3, 150)
    q = stratum_qvalues(et, np.array([0.5]), np.ones(1), 1.0)
    np.testing.assert_allclose(q, 1 / 150, rtol=1e-12)
    assert round(float(q[0]), 5) == 0.00667


def test_example2_b6_50_targets_no_decoys() -> None:
    q = stratum_qvalues(np.linspace(1e-6, 1e-3, 50), np.array([]), np.array([]), 1.0)
    np.testing.assert_allclose(q, 0.02, rtol=1e-12)


def test_example3_stratum1_pi_one() -> None:
    q = stratum_qvalues(np.array([0.001, 0.002, 0.003]), np.array([0.0025]), np.ones(1), 1.0)
    np.testing.assert_allclose(q, [0.5, 0.5, 2 / 3], rtol=1e-12)


def test_example3_input_order_preserved() -> None:
    q = stratum_qvalues(np.array([0.003, 0.001, 0.002]), np.array([0.0025]), np.ones(1), 1.0)
    np.testing.assert_allclose(q, [2 / 3, 0.5, 0.5], rtol=1e-12)


def test_example4_stratum0_pi_002() -> None:
    et, ed = np.array([0.001, 0.002, 0.003]), np.array([0.0025])
    q = stratum_qvalues(et, ed, np.full(1, 0.02), 0.02)
    np.testing.assert_allclose(q, [0.01, 0.01, 0.02 * 2 / 3], rtol=1e-12)


def test_tie_between_target_and_decoy_counts_the_decoy() -> None:
    # at e = 2: T = 4, D = 1 (tied decoy counted) -> 2/4; without the tie it would be 1/4
    q = stratum_qvalues(np.array([1.0, 2.0, 2.0, 2.0]), np.array([2.0]), np.ones(1), 1.0)
    np.testing.assert_allclose(q, [0.5, 0.5, 0.5, 0.5], rtol=1e-12)


def test_tied_targets_share_one_threshold() -> None:
    q = stratum_qvalues(np.array([1.0, 1.0, 3.0]), np.array([2.0]), np.ones(1), 1.0)
    # FDR(1) = 1/2, FDR(3) = 2/3
    np.testing.assert_allclose(q, [0.5, 0.5, 2 / 3], rtol=1e-12)


def test_fdr_above_one_is_clipped() -> None:
    q = stratum_qvalues(np.array([1.0]), np.array([0.1, 0.2, 0.3]), np.ones(3), 1.0)
    np.testing.assert_allclose(q, [1.0])


def test_empty_target_list() -> None:
    q = stratum_qvalues(np.array([]), np.array([0.1, 0.2]), np.ones(2), 1.0)
    assert q.shape == (0,)


def test_per_target_pi_uses_mean_over_targets_up_to_e() -> None:
    """ADR-0010 (supersedes M3-10 Q2): decoys carry their own pi, FDR(e) = (W(e) + c) / T(e)."""
    q = stratum_qvalues(np.array([1.0, 2.0, 3.0]), np.array([2.5]), np.array([0.5]), 0.6)
    # FDR(1) = 0.6/1, FDR(2) = 0.6/2, FDR(3) = (0.5 + 0.6)/3
    np.testing.assert_allclose(q, [0.3, 0.3, 1.1 / 3], rtol=1e-12)
    assert np.round(q, 6).tolist() == [0.3, 0.3, 0.366667]


def test_decoy_weights_sum_up_to_e_and_c_is_added_once() -> None:
    # decoys E = 0.5 (w 0.2), 1.5 (w 0.4), 9 (w 1.0); c = 1
    # FDR(1) = (0.2 + 1)/1, FDR(2) = (0.6 + 1)/2 = 0.8, FDR(3) = 1.6/3
    q = stratum_qvalues(
        np.array([1.0, 2.0, 3.0]), np.array([9.0, 0.5, 1.5]), np.array([1.0, 0.2, 0.4]), 1.0
    )
    np.testing.assert_allclose(q, [1.6 / 3, 1.6 / 3, 1.6 / 3], rtol=1e-12)
    q = stratum_qvalues(np.array([1.0, 2.0]), np.array([0.5, 1.5]), np.array([0.2, 0.4]), 0.5)
    # FDR(1) = 0.7, FDR(2) = 1.1/2 = 0.55
    np.testing.assert_allclose(q, [0.55, 0.55], rtol=1e-12)


def test_q_depends_on_e_only_through_the_ranking() -> None:
    """E is a ranking score only (M3-09 lead decision): monotone transforms of E change no q."""
    et, ed, w = np.array([0.001, 0.002, 0.003]), np.array([0.0025]), np.full(1, 0.02)
    q = stratum_qvalues(et, ed, w, 0.02)
    np.testing.assert_allclose(stratum_qvalues(np.log(et), np.log(ed), w, 0.02), q, rtol=1e-12)
    np.testing.assert_allclose(stratum_qvalues(1e6 * et, 1e6 * ed, w, 0.02), q, rtol=1e-12)


@pytest.mark.parametrize("pi", [1.0, 0.3, 0.02])
def test_constant_pi_equals_m4_tda_qvalues(pi: float) -> None:
    rng = np.random.default_rng(7)
    et = np.round(rng.exponential(1.0, 300), 2)  # rounding makes ties
    ed = np.round(rng.exponential(2.0, 200), 2)
    q = stratum_qvalues(et, ed, np.full(ed.size, pi), pi)
    np.testing.assert_allclose(q, tda_qvalues(-et, -ed, pi=pi, plus_one=True), rtol=1e-12)


# ------------------------------------------------------------------------- pi_wrong_database


def _pep_table(pep: Sequence[float | None], q: Sequence[float | None]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide": [pep_of(i) for i in range(len(pep))],
            "pep": pd.array([np.nan if v is None else v for v in pep], dtype="float64"),
            "q_value": pd.array([np.nan if v is None else v for v in q], dtype="float64"),
        }
    )


def test_example5_mean_pep() -> None:
    t = _pep_table([0.001, 0.003, 0.008], [0.001, 0.001, 0.001])
    assert pi_wrong_database(t, 0.01) == pytest.approx(0.004, rel=1e-12)


def test_example6_max_q_within_cap() -> None:
    t = _pep_table([None, None, None], [0.001, 0.0098, 0.03])
    assert pi_wrong_database(t, 0.01) == pytest.approx(0.0098, rel=1e-12)


def test_example7_neither_pep_nor_q() -> None:
    t = _pep_table([None, None], [None, None])
    assert pi_wrong_database(t, 0.01) == 1.0


def test_mean_pep_over_peptides_with_pep_only() -> None:
    t = _pep_table([0.002, None, 0.004], [None, 0.005, None])
    assert pi_wrong_database(t, 0.01) == pytest.approx(0.003, rel=1e-12)


def test_no_q_within_cap_falls_back_to_one() -> None:
    t = _pep_table([None, None], [0.02, 0.03])
    assert pi_wrong_database(t, 0.01) == 1.0


def test_empty_table_gives_one() -> None:
    assert pi_wrong_database(_pep_table([], []), 0.01) == 1.0


# --------------------------------------------------------------------------------- best_costs


def _uniform_freqs() -> ResidueFreqs:
    return ResidueFreqs(f=np.full(19, 1 / 19), n_positions=1_000_000)


def _hit(peptide: str, search: str, n_diff: int, start: int) -> dict[str, object]:
    return {
        "peptide": peptide,
        "query": peptide,
        "search": search,
        "decoy_type": None if search == "target" else "reverse",
        "contig": "chr1",
        "frame_id": 0,
        "aa_start": start,
        "aa_end": start + len(peptide),
        "start": 3 * start,
        "end": 3 * (start + len(peptide)),
        "strand": "+",
        "frame": 1,
        "route": "exact" if n_diff == 0 else "tolerant",
        "n_diff": n_diff,
        "diff_types": "" if n_diff == 0 else "substitution",
        "blocks": None,
    }


def _peptides(names: Sequence[str]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "peptide": list(names),
            "q_value": [0.001] * len(names),
            "pep": [np.nan] * len(names),
            "is_decoy": [False] * len(names),
            "source": ["database"] * len(names),
            "search_db": ["other"] * len(names),
            "pi_denovo": [np.nan] * len(names),
        }
    )


def test_best_costs_counts_each_peptide_once_at_best_cost() -> None:
    a, b = "PEPTIDEKLR", "MDEFGHKAC"
    hits = hits_frame(
        [
            _hit(a, "target", 1, 10),
            _hit(a, "target", 0, 40),
            _hit(a, "target", 2, 70),
            _hit(a, "decoy", 2, 5),
            _hit(a, "decoy", 1, 9),
            _hit(b, "decoy", 1, 3),
        ]
    )
    freqs = _uniform_freqs()
    out = best_costs(hits, _peptides([a, b]), freqs)
    assert out["peptide"].is_unique
    rows = {r["peptide"]: r for r in out.to_dict("records")}
    ra, rb = rows[a], rows[b]
    assert int(ra["m_target"]) == 0 and int(ra["m_decoy"]) == 1
    assert int(ra["length"]) == len(a)
    assert ra["e_target"] == pytest.approx(expected_hits(a, 0, freqs), rel=1e-12)
    assert ra["e_decoy"] == pytest.approx(expected_hits(a, 1, freqs), rel=1e-12)
    # example 9 input: decoy hit only
    assert pd.isna(rb["m_target"]) and pd.isna(rb["e_target"])
    assert int(rb["m_decoy"]) == 1
    assert rb["e_decoy"] == pytest.approx(expected_hits(b, 1, freqs), rel=1e-12)
    assert ra["source"] == "database" and ra["search_db"] == "other"
    assert {f.name for f in dataclasses.fields(PeptideCosts)} <= set(out.columns)


def test_best_costs_order_independent() -> None:
    a, b = "PEPTIDEKLR", "MDEFGHKAC"
    rows = [_hit(a, "target", 1, 10), _hit(a, "decoy", 0, 5), _hit(b, "target", 2, 3)]
    freqs = _uniform_freqs()
    one = best_costs(hits_frame(rows), _peptides([a, b]), freqs)
    two = best_costs(hits_frame(rows[::-1]), _peptides([b, a]), freqs)
    key = ["peptide"]
    pd.testing.assert_frame_equal(
        one.sort_values(key).reset_index(drop=True),
        two.sort_values(key).reset_index(drop=True)[one.columns],
    )


def test_best_costs_drops_input_decoy_rows() -> None:
    """Input decoy rows (``is_decoy``) are never mapped (M3-17) and never enter the lists."""
    a, rev = "PEPTIDEKLR", "RLKEDITPEP"
    hits = hits_frame([_hit(a, "target", 0, 10), _hit(rev, "target", 0, 30)])
    peps = _peptides([a, rev])
    peps.loc[peps["peptide"] == rev, "is_decoy"] = True
    out = best_costs(hits, peps, _uniform_freqs())
    assert out["peptide"].tolist() == [a]


def test_peptide_costs_is_frozen() -> None:
    pc = PeptideCosts(
        peptide="PEPTIDEK", length=8, source="database", search_db="other", m_target=0,
        m_decoy=None, e_target=0.1, e_decoy=None, pep=None, q_value=None, pi_denovo=None,
    )  # fmt: skip
    with pytest.raises(dataclasses.FrozenInstanceError):
        pc.m_target = 1  # type: ignore[misc]


# ------------------------------------------------------------------------- placement examples


def test_placement_example1_b6() -> None:
    rows = targets(np.linspace(1e-6, 1e-3, 150), 1) + decoys([0.5], 1)
    per, fdr = run(rows)
    names = [r["peptide"] for r in rows[:150]]
    np.testing.assert_allclose(q_of(per, names), 1 / 150, rtol=1e-12)
    d = by_peptide(per)
    assert all(bool(d[p]["accepted"]) for p in names)
    tot = fdr_row(fdr, 1, "all")
    assert int(tot["targets"]) == 150 and int(tot["decoys"]) == 1
    assert int(tot["n_accepted"]) == 150
    assert float(tot["pi_wrong"]) == 1.0
    assert float(tot["fdr_stratum"]) == pytest.approx(2 / 150, rel=1e-12)


def test_placement_example2_b6_none_accepted() -> None:
    rows = targets(np.linspace(1e-6, 1e-3, 50), 1)
    per, fdr = run(rows)
    names = [r["peptide"] for r in rows]
    np.testing.assert_allclose(q_of(per, names), 0.02, rtol=1e-12)
    assert not per.set_index("peptide").loc[names, "accepted"].astype(bool).any()
    tot = fdr_row(fdr, 1, "all")
    assert int(tot["targets"]) == 50 and int(tot["decoys"]) == 0
    assert float(tot["fdr_stratum"]) == pytest.approx(0.02, rel=1e-12)
    assert int(tot["n_accepted"]) == 0


def test_placement_example3_stratum1_ignores_pep() -> None:
    # PEP present: pi_wrong = 1 for m >= 1 all the same
    rows = targets([0.001, 0.002, 0.003], 1, pep=0.02) + decoys([0.0025], 1, pep=0.02)
    per, fdr = run(rows)
    names = [r["peptide"] for r in rows[:3]]
    np.testing.assert_allclose(q_of(per, names), [0.5, 0.5, 2 / 3], rtol=1e-12)
    assert float(fdr_row(fdr, 1, "9")["pi_wrong"]) == 1.0


def test_placement_example4_stratum0_pi_from_pep_inclusive_alpha() -> None:
    rows = targets([0.001, 0.002, 0.003], 0, pep=0.02) + decoys([0.0025], 0, pep=0.02)
    per, fdr = run(rows)
    names = [r["peptide"] for r in rows[:3]]
    np.testing.assert_allclose(q_of(per, names), [0.01, 0.01, 0.02 * 2 / 3], rtol=1e-12)
    acc = per.set_index("peptide").loc[names, "accepted"].astype(bool).tolist()
    assert acc == [True, True, False]
    row = fdr_row(fdr, 0, "9")
    assert float(row["pi_wrong"]) == pytest.approx(0.02, rel=1e-12)
    assert int(row["targets"]) == 3 and int(row["decoys"]) == 1
    assert float(row["fdr_stratum"]) == pytest.approx(0.02 * 2 / 3, rel=1e-12)
    assert int(row["n_accepted"]) == 2


def test_placement_example5_mean_pep_enters_stratum0() -> None:
    peps = [0.001, 0.003, 0.008]
    rows = [
        {"peptide": pep_of(i), "m_target": 0, "e_target": 0.001 * (i + 1), "pep": p}
        for i, p in enumerate(peps)
    ]
    per, fdr = run(rows)
    # no decoys: FDR(e_k) = 0.004 / k -> q = 0.004/3 for all
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in rows]), 0.004 / 3, rtol=1e-12)
    assert float(fdr_row(fdr, 0, "all")["pi_wrong"]) == pytest.approx(0.004, rel=1e-12)


def test_placement_example6_max_q_enters_stratum0() -> None:
    rows = [
        {"peptide": pep_of(i), "m_target": 0, "e_target": 0.001 * (i + 1), "q_value": q}
        for i, q in enumerate([0.001, 0.0098, 0.03])
    ]
    per, fdr = run(rows, q_value_max=0.01)
    assert float(fdr_row(fdr, 0, "all")["pi_wrong"]) == pytest.approx(0.0098, rel=1e-12)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in rows]), 0.0098 / 3, rtol=1e-12)


def test_placement_example7_no_pep_no_q_pi_one() -> None:
    rows = targets([0.001, 0.002], 0)
    per, fdr = run(rows)
    assert float(fdr_row(fdr, 0, "all")["pi_wrong"]) == 1.0
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in rows]), 0.5, rtol=1e-12)


def test_placement_example8_same_db() -> None:
    same_exact = {
        "peptide": pep_of(900), "search_db": "same", "m_target": 0, "e_target": 0.001,
        "m_decoy": 0, "e_decoy": 0.0005,
    }  # fmt: skip
    twin = {"peptide": pep_of(901), "search_db": "same", "m_target": 1, "e_target": 0.001}
    others = targets([0.002, 0.004], 0)
    per, fdr = run([same_exact, twin, *others])
    d = by_peptide(per)
    assert float(d[same_exact["peptide"]]["placement_q"]) == 0.0
    assert "same_db" in flags_of(d[same_exact["peptide"]]["flags"])
    assert bool(d[same_exact["peptide"]]["accepted"])
    # neither its target nor its decoy hit enters stratum 0: q = 1/1, 1/2 -> 0.5, 0.5
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in others]), [0.5, 0.5], rtol=1e-12)
    tot0 = fdr_row(fdr, 0, "all")
    assert int(tot0["targets"]) == 2 and int(tot0["decoys"]) == 0
    # the 1-difference twin goes through steps 1-4 (alone in stratum 1: q = 1/1)
    assert "same_db" not in flags_of(d[twin["peptide"]]["flags"])
    assert float(d[twin["peptide"]]["placement_q"]) == pytest.approx(1.0)
    assert not bool(d[twin["peptide"]]["accepted"])
    assert int(fdr_row(fdr, 1, "all")["targets"]) == 1


def test_placement_example9_decoy_only_peptide() -> None:
    t = targets([0.001, 0.003], 1)
    dec = decoys([0.002], 1)
    per, fdr = run(t + dec)
    # with the decoy: FDR = 1/1, 2/2 -> q = 1, 1 (without it would be 0.5, 0.5)
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in t]), [1.0, 1.0], rtol=1e-12)
    assert unplaced(per, dec[0]["peptide"])
    assert int(fdr_row(fdr, 1, "all")["decoys"]) == 1
    assert int(fdr_row(fdr, 0, "all")["decoys"]) == 0
    assert int(fdr_row(fdr, 2, "all")["decoys"]) == 0


def test_placement_example10_length_class_rows() -> None:
    rows = [
        {"peptide": pep_of(1, 8), "m_target": 1, "e_target": 0.001},
        {"peptide": pep_of(2, 12), "m_target": 1, "e_target": 0.002},
    ]
    _, fdr = run(rows)
    assert set(fdr.columns) >= FDR_COLUMNS
    assert int(fdr_row(fdr, 1, "8")["targets"]) == 1
    assert int(fdr_row(fdr, 1, "12-14")["targets"]) == 1
    assert int(fdr_row(fdr, 1, "all")["targets"]) == 2
    for cls in ("7", "9", "10-11", "15-19", "20-40"):
        r = fdr_row(fdr, 1, cls)
        assert int(r["targets"]) == 0 and pd.isna(r["fdr_stratum"])


# --------------------------------------------------------------------------- mapping_fdr shape


def test_mapping_fdr_rows_complete_and_sorted() -> None:
    rows = targets([0.001], 0) + targets([0.01], 2, start=10, length=15) + decoys([0.1], 1)
    _, fdr = run(rows)
    if "section" in fdr.columns:
        assert set(fdr["section"]) == {"placement"}
    rev = fdr[fdr["decoy_type"] == "reverse"]
    assert len(rev) == 3 * (len(CLASS_LABELS) + 1)
    for m in (0, 1, 2):
        labels = rev[rev["n_diff"].astype(int) == m]["length_class"].astype(str).tolist()
        assert sorted(labels) == sorted([*CLASS_LABELS, "all"])
        assert [x for x in labels if x != "all"] == CLASS_LABELS
    assert set(fdr["decoy_type"]) <= {"reverse", "shuffle"}
    key = list(zip(fdr["decoy_type"].astype(str), fdr["n_diff"].astype(int), strict=True))
    assert key == sorted(key)
    # whole-class counts, no threshold
    assert int(fdr_row(fdr, 2, "15-19")["targets"]) == 1
    assert int(fdr_row(fdr, 1, "9")["decoys"]) == 1


def test_stratum_with_only_decoys_has_na_fdr() -> None:
    _, fdr = run(decoys([0.001, 0.002], 0))
    r = fdr_row(fdr, 0, "all")
    assert int(r["targets"]) == 0 and int(r["decoys"]) == 2
    assert pd.isna(r["fdr_stratum"]) and int(r["n_accepted"]) == 0


# ------------------------------------------------------------------- decoy types and pi_denovo


def test_lists_are_per_decoy_type() -> None:
    db = targets([0.001, 0.002], 0)
    dn_decoy = decoys([0.0005], 0, source="denovo")
    per, fdr = run(db + dn_decoy)
    # the shuffle decoy does not enter the database (reverse) list: q = 1/1, 1/2 -> 0.5, 0.5
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in db]), [0.5, 0.5], rtol=1e-12)
    assert int(fdr_row(fdr, 0, "all", "reverse")["decoys"]) == 0
    sh = fdr_row(fdr, 0, "all", "shuffle")
    assert int(sh["decoys"]) == 1 and int(sh["targets"]) == 0


def test_denovo_stratum0_uses_per_peptide_pi() -> None:
    """ADR-0010: the decoy hit weighs its own pi (0.5), c = max pi of the searched de novo
    peptides (0.6); the targets' pi enters only through c."""
    t = [
        {"peptide": pep_of(i, 7), "source": "denovo", "m_target": 0, "e_target": e, "pi_denovo": p}
        for i, (e, p) in enumerate([(1.0, 0.5), (2.0, 0.6), (3.0, 0.3)])
    ]
    dec = decoys([2.5], 0, length=7, source="denovo", pi_denovo=0.5)
    per, fdr = run(t + dec)
    # FDR = 0.6/1, 0.6/2, (0.5 + 0.6)/3
    np.testing.assert_allclose(
        q_of(per, [r["peptide"] for r in t]), [0.3, 0.3, 1.1 / 3], rtol=1e-12
    )
    for cls in ("7", "all"):
        r = fdr_row(fdr, 0, cls, "shuffle")
        assert int(r["targets"]) == 3 and int(r["decoys"]) == 1
        assert float(r["fdr_stratum"]) == pytest.approx(1.1 / 3, rel=1e-12)
        assert float(r["pi_wrong"]) == pytest.approx(0.55, rel=1e-12)


def test_denovo_c_is_max_pi_over_all_searched_peptides() -> None:
    """c counts every searched de novo peptide, also one with no stratum-0 hit (here a target
    hit at m = 1 with pi 0.9); stratum 1 has pi = 1 and c = 1 regardless."""
    t = [
        {"peptide": pep_of(i, 7), "source": "denovo", "m_target": 0, "e_target": e, "pi_denovo": p}
        for i, (e, p) in enumerate([(1.0, 0.1), (2.0, 0.2), (3.0, 0.1), (4.0, 0.1)])
    ]
    other = {
        "peptide": pep_of(50, 7), "source": "denovo", "m_target": 1, "e_target": 1.0,
        "pi_denovo": 0.9,
    }  # fmt: skip
    per, fdr = run([*t, other])
    # no stratum-0 decoys: FDR(e_k) = 0.9 / k -> q = 0.9/4
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in t]), 0.225, rtol=1e-12)
    r0 = fdr_row(fdr, 0, "all", "shuffle")
    assert float(r0["pi_wrong"]) == pytest.approx(0.9, rel=1e-12)
    assert float(r0["fdr_stratum"]) == pytest.approx(0.225, rel=1e-12)
    assert float(by_peptide(per)[other["peptide"]]["placement_q"]) == pytest.approx(1.0)


def test_denovo_c_is_one_when_a_searched_peptide_lacks_pi() -> None:
    t = targets([1.0, 2.0, 3.0, 4.0], 0, length=7, source="denovo", pi_denovo=0.1)
    no_pi = {"peptide": pep_of(60, 7), "source": "denovo", "m_decoy": 1, "e_decoy": 1.0}
    per, fdr = run([*t, no_pi])
    np.testing.assert_allclose(q_of(per, [r["peptide"] for r in t]), 0.25, rtol=1e-12)
    assert float(fdr_row(fdr, 0, "all", "shuffle")["pi_wrong"]) == pytest.approx(1.0)


def test_mapping_fdr_class_rows_use_class_w_and_stratum_c() -> None:
    """Rows: fdr_stratum = (W + c)/T with W over the class's decoy hits and the stratum's c;
    pi_wrong = (W + c)/(D + 1)."""
    t7 = targets([1.0, 3.0], 0, length=7, source="denovo", pi_denovo=0.2)
    t8 = targets([2.0], 0, start=10, length=8, source="denovo", pi_denovo=0.8)
    d7 = decoys([1.5], 0, length=7, source="denovo", pi_denovo=0.4)
    per, fdr = run(t7 + t8 + d7)
    c = 0.8
    r7, r8, ra = (fdr_row(fdr, 0, k, "shuffle") for k in ("7", "8", "all"))
    assert float(r7["fdr_stratum"]) == pytest.approx((0.4 + c) / 2, rel=1e-12)
    assert float(r7["pi_wrong"]) == pytest.approx((0.4 + c) / 2, rel=1e-12)
    assert float(r8["fdr_stratum"]) == pytest.approx(c / 1, rel=1e-12)
    assert float(r8["pi_wrong"]) == pytest.approx(c, rel=1e-12)
    assert float(ra["fdr_stratum"]) == pytest.approx((0.4 + c) / 3, rel=1e-12)
    assert float(ra["pi_wrong"]) == pytest.approx((0.4 + c) / 2, rel=1e-12)
    for r in (r7, r8, ra):
        assert float(r["fdr_stratum"]) == pytest.approx(
            float(r["pi_wrong"]) * (int(r["decoys"]) + 1) / int(r["targets"]), rel=1e-12
        )
    # one list per stratum: FDR = 0.8/1, (0.4 + 0.8)/2, 1.2/3 -> q = 0.4, 0.4, 0.4
    names = [r["peptide"] for r in t7 + t8]
    np.testing.assert_allclose(q_of(per, names), 0.4, rtol=1e-12)


def test_denovo_missing_pi_is_one_and_m1_ignores_pi() -> None:
    t0 = [{"peptide": pep_of(1), "source": "denovo", "m_target": 0, "e_target": 1.0}]
    t1 = [
        {"peptide": pep_of(2), "source": "denovo", "m_target": 1, "e_target": 1.0,
         "pi_denovo": 0.1},
    ]  # fmt: skip
    per, fdr = run(t0 + t1)
    assert float(by_peptide(per)[pep_of(1)]["placement_q"]) == pytest.approx(1.0)
    assert float(by_peptide(per)[pep_of(2)]["placement_q"]) == pytest.approx(1.0)
    assert float(fdr_row(fdr, 1, "all", "shuffle")["pi_wrong"]) == 1.0


# ----------------------------------------------------------------------------- flags and misc


def test_low_complexity_flag() -> None:
    rows = [
        {"peptide": "AKAKAKAKAK", "m_target": 1, "e_target": 0.001},
        {"peptide": "PPPPPPPPP", "m_target": 0, "e_target": 0.001},
        {"peptide": "PEPTIDEKLR", "m_target": 1, "e_target": 0.002},
    ]
    per, _ = run(rows)
    d = by_peptide(per)
    assert "low_complexity" in flags_of(d["AKAKAKAKAK"]["flags"])
    assert "low_complexity" in flags_of(d["PPPPPPPPP"]["flags"])
    assert "low_complexity" not in flags_of(d["PEPTIDEKLR"]["flags"])


def test_per_peptide_columns_and_m_target() -> None:
    rows = targets([0.001], 2)
    per, _ = run(rows)
    assert {"peptide", "m_target", "placement_q", "accepted", "flags"} <= set(per.columns)
    assert int(by_peptide(per)[rows[0]["peptide"]]["m_target"]) == 2


def test_q_values_within_unit_interval_when_d_exceeds_t() -> None:
    per, _ = run(targets([1.0], 0) + decoys([0.1, 0.2, 0.3], 0))
    assert float(per["placement_q"].dropna().iloc[0]) == 1.0


def test_alpha_from_config() -> None:
    rows = targets([0.001, 0.002, 0.003], 1) + decoys([0.0025], 1)
    per, fdr = run(rows, alpha=0.5)
    acc = per.set_index("peptide").loc[[r["peptide"] for r in rows[:3]], "accepted"]
    assert acc.astype(bool).tolist() == [True, True, False]
    assert int(fdr_row(fdr, 1, "all")["n_accepted"]) == 2


def test_empty_input() -> None:
    per, _ = placement(costs_frame([]), cfg())
    assert len(per) == 0


def _mixed_rows() -> list[dict[str, Any]]:
    rng = np.random.default_rng(3)
    rows: list[dict[str, Any]] = []
    for i in range(120):
        src = "denovo" if i % 3 == 0 else "database"
        L = int(rng.integers(7, 25))
        r: dict[str, Any] = {"peptide": pep_of(i, L), "source": src, "pep": 0.01}
        if i % 11 == 0 and src == "database":
            r["search_db"] = "same"
        if rng.random() < 0.8:
            r["m_target"] = int(rng.integers(0, 3))
            r["e_target"] = float(np.round(rng.exponential(1.0), 1))
        if rng.random() < 0.4:
            r["m_decoy"] = int(rng.integers(0, 3))
            r["e_decoy"] = float(np.round(rng.exponential(1.0), 1))
        if src == "denovo":
            r["pi_denovo"] = float(rng.uniform(0.05, 1.0))
        rows.append(r)
    return rows


def test_deterministic() -> None:
    a = run(_mixed_rows())
    b = run(_mixed_rows())
    pd.testing.assert_frame_equal(a[0], b[0])
    pd.testing.assert_frame_equal(a[1], b[1])


def test_order_independence() -> None:
    rows = _mixed_rows()
    per1, fdr1 = run(rows)
    perm = np.random.default_rng(11).permutation(len(rows))
    per2, fdr2 = run([rows[i] for i in perm])
    s1 = per1.sort_values("peptide").reset_index(drop=True)
    s2 = per2.sort_values("peptide").reset_index(drop=True)
    pd.testing.assert_frame_equal(s1, s2)
    pd.testing.assert_frame_equal(fdr1.reset_index(drop=True), fdr2.reset_index(drop=True))


def test_mixed_input_invariants() -> None:
    rows = _mixed_rows()
    per, fdr = run(rows)
    d = by_peptide(per)
    for r in rows:
        if r.get("m_target") is None:
            assert unplaced(per, r["peptide"])
            continue
        q = float(d[r["peptide"]]["placement_q"])
        assert 0.0 <= q <= 1.0
        assert bool(d[r["peptide"]]["accepted"]) == (q <= 0.01)
        if r.get("search_db") == "same" and r["m_target"] == 0:
            assert q == 0.0
    # 'all' rows sum the class rows
    for (dt, m), g in fdr.groupby(["decoy_type", "n_diff"]):
        cls = g[g["length_class"].astype(str) != "all"]
        tot = g[g["length_class"].astype(str) == "all"].iloc[0]
        assert int(tot["targets"]) == int(cls["targets"].sum()), (dt, m)
        assert int(tot["decoys"]) == int(cls["decoys"].sum()), (dt, m)
        assert int(tot["n_accepted"]) == int(cls["n_accepted"].sum()), (dt, m)
