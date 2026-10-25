"""V3 locus order within each §8.1 family, not across scales
(docs/specs/V03-order-within-kind.md, with its "Statistician review" and "Lead decisions" of
2026-10-25; GIAP_SPEC.md §10.2 V3, §8.1).

The rule under test: each locus gets the key ``(tier, family_rank, family_index, locus_id)``.

- ``tier`` is 0 when ``error <= alpha(family)`` (plain float ``<=``), 1 when above, 2 when NA.
  ``alpha`` is ``stats.alpha.context`` for prophage, bgc and island, and
  ``stats.alpha.operon_flags`` for both flag families, read from the ``cfg`` passed in.
- The families, in order: prophage, bgc, island, operon_discordance (``uniop_class = OP``
  flags) and missed_operon (``uniop_class = NO`` flags). ``Locus.kind`` stays ``operon_flag``
  for both flag families; the tests know each flag's family from the rows they write.
- ``family_rank`` is the 1-based position among the loci of the same family and tier sorted by
  ``(error, locus_id)``; 0 in tier 2.
- Operon flags merge by ``(family, contig, span)``.

Errors are written to the TSVs as ``repr`` strings so that values at ``alpha`` and their
``math.nextafter`` neighbours survive the round trip (``format_value`` would print ``%.6g``).
"""

from __future__ import annotations

import functools
import itertools
import math
import random
import re
import tempfile
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

import pandas as pd
import pytest
from helpers import snapshot_fixture as sf
from helpers import viz_locus as vl
from hypothesis import HealthCheck, assume, given, settings
from hypothesis import strategies as st
from matplotlib.text import Text

from giap.config import GiapConfig, resolve_config
from giap.io.run_view import RunView
from giap.viz import base, style
from giap.viz.base import Instance, Skip
from giap.viz.locus import Locus, select_loci

FAMILIES = ("prophage", "bgc", "island", "operon_discordance", "missed_operon")
FLAG_FAMILIES = ("operon_discordance", "missed_operon")
KIND = {
    "prophage": "prophage",
    "bgc": "bgc",
    "island": "island",
    "operon_discordance": "operon_flag",
    "missed_operon": "operon_flag",
}
ESTIMATE = {
    "prophage": "q",
    "bgc": "q",
    "island": "FWER p",
    "operon_discordance": "q",
    "missed_operon": "decoy-FDR q",
}
HEADERS = ("significant at each family's own level", "not significant", "no estimate")
CAPTION = (
    "Order alternates between test families; values in different kinds are on different "
    "scales and are not compared."
)

#: one pool entry: (family, locus_id, error); error None = NA
Entry = tuple[str, str, float | None]
Key = tuple[str, str]  # (kind, locus_id): unique in the pool, the same ID may be in two kinds

PROPS = settings(
    max_examples=40,
    deadline=None,
    suppress_health_check=[HealthCheck.too_slow, HealthCheck.filter_too_much],
)


@pytest.fixture(autouse=True)
def _style() -> None:
    style.apply()


@functools.lru_cache(maxsize=64)
def _cfg(*overrides: str) -> GiapConfig:
    return resolve_config(overrides=list(overrides))


def _v03() -> base.PlotSpec:
    import giap.viz.locus  # noqa: F401  (registers V03)

    return base.registered()["V03"]


def _run(out: Path, cfg: GiapConfig) -> tuple[Instance, list[Instance]]:
    result = _v03().fn(RunView.from_outdir(out), cfg)
    assert not isinstance(result, Skip), result
    index = [i for i in result if i.suffix is None]
    loci = [i for i in result if i.suffix is not None]
    assert len(index) == 1, [i.suffix for i in result]
    return index[0], loci


def _ids(loci: Sequence[Locus]) -> list[str]:
    return [lo.locus_id for lo in loci]


def _keys(loci: Sequence[Locus]) -> list[Key]:
    return [(lo.kind, lo.locus_id) for lo in loci]


def _alphas(cfg: GiapConfig) -> dict[str, float]:
    a = cfg.stats.alpha
    return {f: (a.operon_flags if f in FLAG_FAMILIES else a.context) for f in FAMILIES}


# --- building an OUT/ tree from entries ------------------------------------------------------


def _cell(e: float | None) -> str | None:
    return None if e is None else repr(float(e))


def _flag_row(fam: str, pid: str, a: str, b: str, e: float | None) -> dict[str, object]:
    if fam == "operon_discordance":
        return vl.flag(pid, a, b, _cell(e))  # type: ignore[arg-type]
    row = vl.flag(pid, a, b, _cell(e), flag="novel_operon_candidate")  # type: ignore[arg-type]
    row.update(uniop_class="NO", abundance_pattern="coabundant")
    return row


def _build(tmp: Path, entries: Sequence[Entry], shuffle_seed: int | None = None) -> Path:
    """One OUT/ with every entry on chr1 in its own 1 kb slot; each flag has its own TU of two
    genes, so no two flags merge. ``shuffle_seed`` permutes the rows of the three tables."""
    regions: list[dict[str, object]] = []
    islands: list[dict[str, object]] = []
    flags: list[dict[str, object]] = []
    genes: list[dict[str, object]] = []
    tus: list[tuple] = []
    for i, (fam, lid, e) in enumerate(entries):
        s = i * 1_500 + 1
        if fam in ("prophage", "bgc"):
            regions.append(vl.region_stat(lid, fam, "chr1", s, s + 999, _cell(e)))  # type: ignore[arg-type]
        elif fam == "island":
            islands.append(vl.island(lid, "chr1", s, s + 999, _cell(e)))  # type: ignore[arg-type]
        else:
            ga, gb, tid = f"g{i}a", f"g{i}b", f"TU{i + 1:04d}"
            genes += [
                vl.gene(ga, "chr1", s, s + 399, tu_id=tid, log2=20.0),
                vl.gene(gb, "chr1", s + 500, s + 999, tu_id=tid, log2=21.0),
            ]
            tus.append(vl.tu(tid, "chr1", s, s + 999, "+", [ga, gb]))
            flags.append(_flag_row(fam, lid, ga, gb, e))
    if shuffle_seed is not None:
        rng = random.Random(shuffle_seed)
        for rows in (regions, islands, flags):
            rng.shuffle(rows)
    return vl.build_out(tmp, genes=genes, tus=tus, regions_stats=regions, islands=islands,
                        flags=flags)  # fmt: skip


def _select(
    entries: Sequence[Entry], cfg: GiapConfig, shuffle_seed: int | None = None
) -> list[Key]:
    with tempfile.TemporaryDirectory() as tmp:
        out = _build(Path(tmp), entries, shuffle_seed)
        return _keys(select_loci(RunView.from_outdir(out), cfg))


# --- the oracle pieces (computed from the entries, never from the implementation) -----------


def _family_of(entries: Sequence[Entry]) -> dict[Key, str]:
    return {(KIND[f], lid): f for f, lid, _ in entries}


def _error_of(entries: Sequence[Entry]) -> dict[Key, float | None]:
    return {(KIND[f], lid): e for f, lid, e in entries}


def _in_pool(entries: Sequence[Entry], cfg: GiapConfig) -> list[Entry]:
    a = cfg.stats.alpha.operon_flags
    return [t for t in entries if t[0] not in FLAG_FAMILIES or (t[2] is not None and t[2] <= a)]


def _tier(fam: str, e: float | None, alphas: Mapping[str, float]) -> int:
    if e is None:
        return 2
    return 0 if e <= alphas[fam] else 1


def _family_ranks(entries: Sequence[Entry], alphas: Mapping[str, float]) -> dict[Key, int]:
    """1-based rank within (family, tier) by (error, locus_id); 0 in tier 2."""
    ranks: dict[Key, int] = {}
    for fam in FAMILIES:
        for tier in (0, 1):
            mine = sorted(
                ((e, lid) for f, lid, e in entries if f == fam and _tier(f, e, alphas) == tier),
                key=lambda t: (t[0], t[1]),
            )
            for r, (_, lid) in enumerate(mine, start=1):
                ranks[(KIND[fam], lid)] = r
        for f, lid, e in entries:
            if f == fam and e is None:
                ranks[(KIND[f], lid)] = 0
    return ranks


def _expected(entries: Sequence[Entry], cfg: GiapConfig) -> list[Key]:
    """The full key ``(tier, family_rank, family_index, locus_id)`` on the pool."""
    pool, alphas = _in_pool(entries, cfg), _alphas(cfg)
    ranks = _family_ranks(pool, alphas)
    return [
        (KIND[f], lid)
        for f, lid, e in sorted(
            pool,
            key=lambda t: (_tier(t[0], t[2], alphas), ranks[(KIND[t[0]], t[1])],
                           FAMILIES.index(t[0]), t[1]),
        )
    ]  # fmt: skip


# --- generators ------------------------------------------------------------------------------

A_CTX, A_FLAG = 0.05, 0.10  # the resolve_config() defaults, asserted below
CTX_SPECIAL = (None, 0.0, 0.001, 0.01, A_CTX, math.nextafter(A_CTX, 0.0),
               math.nextafter(A_CTX, 1.0), 0.06, 0.5, 1.0)  # fmt: skip
FLAG_SPECIAL = (0.0, 0.01, 0.05, A_FLAG, math.nextafter(A_FLAG, 0.0),
                math.nextafter(A_FLAG, 1.0), 0.2, None)  # fmt: skip
#: shared by regions and islands, with flag-like IDs, so the same ID turns up in two kinds
REGION_IDS = ("PAIR_1", "PAIR_2", "PAIR_10", "X_1", "X_2", "X_10", "prophage_1")
FLAG_IDS = tuple(f"PAIR_{i}" for i in (1, 2, 3, 10, 11, 12, 20))

ctx_errors = st.one_of(
    st.sampled_from(CTX_SPECIAL), st.floats(0.0, 1.0, allow_nan=False, allow_infinity=False)
)
flag_errors = st.one_of(
    st.sampled_from(FLAG_SPECIAL), st.floats(0.0, A_FLAG, allow_nan=False, allow_infinity=False)
)


@st.composite
def pools(draw: st.DrawFn, min_size: int = 2) -> list[Entry]:
    entries: list[Entry] = []
    for fam in ("prophage", "bgc", "island"):
        rows = draw(st.lists(st.tuples(st.sampled_from(REGION_IDS), ctx_errors), max_size=4,
                             unique_by=lambda t: t[0]))  # fmt: skip
        entries += [(fam, lid, e) for lid, e in rows]
    rows3 = draw(
        st.lists(st.tuples(st.sampled_from(FLAG_IDS), st.sampled_from(FLAG_FAMILIES), flag_errors),
                 max_size=6, unique_by=lambda t: t[0])
    )  # fmt: skip
    entries += [(fam, lid, e) for lid, fam, e in rows3]
    assume(len(_in_pool(entries, _cfg())) >= min_size)
    return draw(st.permutations(entries))


def test_generator_defaults_match_resolve_config() -> None:
    assert (_cfg().stats.alpha.context, _cfg().stats.alpha.operon_flags) == (A_CTX, A_FLAG)


# --- order-preserving maps (the metamorphic property) ----------------------------------------

#: (shape, parameter): "scale" multiplies the distance to the tier's anchor by c in (0, 1];
#: "power" maps the tier's interval onto itself by x -> x**k
maps = st.tuples(st.sampled_from(("scale", "power")), st.floats(0.01, 1.0), st.floats(0.25, 4.0))


def _g(e: float | None, alpha: float, m: tuple[str, float, float]) -> float | None:
    if e is None:
        return None
    shape, c, k = m
    if e <= alpha:  # [0, alpha] onto [0, alpha]
        return e * c if shape == "scale" else alpha * (e / alpha) ** k if alpha else e
    if shape == "scale":  # (alpha, 1] towards alpha
        return alpha + (e - alpha) * c
    return alpha + (1 - alpha) * ((e - alpha) / (1 - alpha)) ** k


def _valid_map(values: Sequence[float | None], g: Callable[[float | None], float | None],
               alpha: float) -> bool:  # fmt: skip
    """``g`` is strictly increasing on ``values``, keeps ``e <= alpha`` and NA, stays in [0, 1]."""
    vs = sorted({v for v in values if v is not None})
    if any(g(v) is None for v in vs) or any(v is None and g(v) is not None for v in values):
        return False
    mapped = [float(g(v)) for v in vs]  # type: ignore[arg-type]
    if any(not 0.0 <= x <= 1.0 for x in mapped):
        return False
    if any(b <= a for a, b in itertools.pairwise(mapped)):
        return False
    return all((v <= alpha) == (x <= alpha) for v, x in zip(vs, mapped, strict=True))


def _apply(entries: Sequence[Entry], fam: str,
           g: Callable[[float | None], float | None]) -> list[Entry]:  # fmt: skip
    return [(f, lid, g(e) if f == fam else e) for f, lid, e in entries]


# --- examples (spec table) -------------------------------------------------------------------

EXAMPLES: dict[str, tuple[list[Entry], list[str]]] = {
    "2": (
        [("prophage", "prophage_1", 0.04), ("bgc", "bgc_1", 0.001), ("island", "ISL_1", 0.0005),
         ("island", "ISL_2", 0.002), ("operon_discordance", "PAIR_1", 0.08)],
        ["prophage_1", "bgc_1", "ISL_1", "PAIR_1", "ISL_2"],
    ),
    "3": (
        [("prophage", "prophage_1", 0.06), ("operon_discordance", "PAIR_1", 0.09)],
        ["PAIR_1", "prophage_1"],
    ),
    "5": (
        [("island", "ISL_2", 0.005), ("island", "ISL_1", 0.005), ("prophage", "prophage_1", 0.01)],
        ["prophage_1", "ISL_1", "ISL_2"],
    ),
    "7": (
        [("bgc", "bgc_1", 0.05), ("prophage", "prophage_2", 0.05),
         ("prophage", "prophage_1", 0.05), ("island", "ISL_1", 0.05),
         ("operon_discordance", "PAIR_1", 0.05)],
        ["prophage_1", "bgc_1", "ISL_1", "PAIR_1", "prophage_2"],
    ),
    "8": (
        [("prophage", "prophage_1", 0.01), ("prophage", "prophage_2", 0.02),
         ("prophage", "prophage_3", 0.5), ("prophage", "prophage_4", 0.7),
         ("island", "ISL_1", 0.001), ("bgc", "bgc_1", 0.3)],
        ["prophage_1", "ISL_1", "prophage_2", "prophage_3", "bgc_1", "prophage_4"],
    ),
    "11": (
        [("prophage", "prophage_1", 0.01), ("operon_discordance", "PAIR_1", 0.08),
         ("missed_operon", "PAIR_2", 0.01), ("operon_discordance", "PAIR_3", 0.02)],
        ["prophage_1", "PAIR_3", "PAIR_2", "PAIR_1"],
    ),
    "only flags, both families": (
        [("operon_discordance", "PAIR_1", 0.01), ("operon_discordance", "PAIR_2", 0.02),
         ("missed_operon", "PAIR_3", 0.03)],
        ["PAIR_1", "PAIR_3", "PAIR_2"],
    ),
}  # fmt: skip


@pytest.mark.parametrize("name", list(EXAMPLES))
def test_spec_examples(name: str) -> None:
    entries, expected = EXAMPLES[name]
    got = _select(entries, _cfg())
    assert [lid for _, lid in got] == expected
    assert got == _expected(entries, _cfg())  # the test oracle agrees with the spec table


def test_example4_flags_do_not_starve_regions_under_the_cap(tmp_path: Path) -> None:
    entries: list[Entry] = [("operon_discordance", f"PAIR_{i}", i / 10_000) for i in range(1, 61)]
    entries += [("prophage", "prophage_1", 0.03), ("bgc", "bgc_1", 0.04)]
    cfg = _cfg()
    assert cfg.report.max_locus_plots == 50
    flags = [f"PAIR_{i}" for i in range(1, 61)]
    expected = ["prophage_1", "bgc_1", *flags]
    out = _build(tmp_path, entries)
    assert _ids(select_loci(RunView.from_outdir(out), cfg)) == expected
    index, rendered = _run(out, cfg)
    assert [i.suffix for i in rendered] == ["prophage_1", "bgc_1", *flags[:48]]
    assert len(index.data) == 62
    assert list(index.data["rank"]) == list(range(1, 63))
    assert list(index.data["locus_id"]) == expected
    assert [bool(v) for v in index.data["rendered"]] == [True] * 50 + [False] * 12


def test_example9_on_the_sim_plots_fixture() -> None:
    assert sf.missing(sf.RICH) is None
    view = RunView.from_outdir(sf.fixture_out(sf.RICH))
    loci = select_loci(view, sf.fixture_cfg(sf.RICH))
    assert _ids(loci) == ["prophage_1", "ISL_1", "bgc_1"]  # PAIR_176 (q 0.881) is dropped


# --- tier boundaries at alpha ----------------------------------------------------------------

_ABOVE_CTX, _BELOW_CTX = math.nextafter(A_CTX, 1.0), math.nextafter(A_CTX, 0.0)


@pytest.mark.parametrize(
    ("q", "expected"),
    [
        (A_CTX, ["prophage_1", "ISL_1", "ISL_2"]),  # = alpha: tier 0
        (_BELOW_CTX, ["prophage_1", "ISL_1", "ISL_2"]),
        (_ABOVE_CTX, ["ISL_1", "ISL_2", "prophage_1"]),  # one ulp above: tier 1
    ],
    ids=["at_alpha", "ulp_below", "ulp_above"],
)
def test_region_tier_is_a_plain_float_comparison_with_alpha(q: float, expected: list[str]) -> None:
    entries: list[Entry] = [("prophage", "prophage_1", q), ("island", "ISL_1", 0.01),
                            ("island", "ISL_2", 0.02)]  # fmt: skip
    assert [lid for _, lid in _select(entries, _cfg())] == expected


@pytest.mark.parametrize(
    ("q", "expected"),
    [
        (A_FLAG, ["prophage_1", "PAIR_1", "prophage_2"]),
        (math.nextafter(A_FLAG, 0.0), ["prophage_1", "PAIR_1", "prophage_2"]),
        (math.nextafter(A_FLAG, 1.0), ["prophage_1", "prophage_2"]),  # dropped from the pool
    ],
    ids=["at_alpha", "ulp_below", "ulp_above"],
)
def test_flag_at_alpha_is_tier0_and_above_is_dropped(q: float, expected: list[str]) -> None:
    entries: list[Entry] = [("prophage", "prophage_1", 0.001), ("prophage", "prophage_2", 0.002),
                            ("operon_discordance", "PAIR_1", q)]  # fmt: skip
    assert [lid for _, lid in _select(entries, _cfg())] == expected


def test_same_id_in_two_kinds_breaks_by_family_order() -> None:
    entries: list[Entry] = [("island", "PAIR_1", 0.01), ("missed_operon", "PAIR_1", 0.01),
                            ("prophage", "PAIR_1", 0.01), ("bgc", "PAIR_1", 0.01)]  # fmt: skip
    want = [("prophage", "PAIR_1"), ("bgc", "PAIR_1"), ("island", "PAIR_1"),
            ("operon_flag", "PAIR_1")]  # fmt: skip
    assert _select(entries, _cfg()) == want


# --- the levels come from cfg ----------------------------------------------------------------


def test_example10_context_level_is_read_from_cfg() -> None:
    entries: list[Entry] = [("prophage", "prophage_1", 0.03), ("island", "ISL_1", 0.005),
                            ("island", "ISL_2", 0.02)]  # fmt: skip
    low = _select(entries, _cfg("stats.alpha.context=0.01"))
    assert [lid for _, lid in low] == ["ISL_1", "prophage_1", "ISL_2"]
    default = _select(entries, _cfg())
    assert [lid for _, lid in default] == ["prophage_1", "ISL_1", "ISL_2"]


def test_a_higher_context_level_moves_a_region_into_tier0() -> None:
    entries: list[Entry] = [("prophage", "prophage_1", 0.15), ("island", "ISL_1", 0.01),
                            ("island", "ISL_2", 0.02)]  # fmt: skip
    assert [lid for _, lid in _select(entries, _cfg())] == ["ISL_1", "ISL_2", "prophage_1"]
    high = _select(entries, _cfg("stats.alpha.context=0.2"))
    assert [lid for _, lid in high] == ["prophage_1", "ISL_1", "ISL_2"]


def test_operon_flags_level_is_read_from_cfg() -> None:
    entries: list[Entry] = [("prophage", "prophage_1", 0.02), ("prophage", "prophage_2", 0.03),
                            ("operon_discordance", "PAIR_1", 0.04),
                            ("missed_operon", "PAIR_2", 0.06)]  # fmt: skip
    got = _select(entries, _cfg("stats.alpha.operon_flags=0.05"))
    # PAIR_2 (q 0.06 > 0.05) leaves the pool; PAIR_1 (q 0.04) stays tier 0 and takes turns
    assert [lid for _, lid in got] == ["prophage_1", "PAIR_1", "prophage_2"]


# --- operon flags merge by (family, contig, span) --------------------------------------------


def test_flags_of_two_families_in_one_tu_stay_two_loci(tmp_path: Path) -> None:
    genes = [
        vl.gene("g3", "chr1", 40_001, 41_000, tu_id="TU0001", log2=20.0),
        vl.gene("g4", "chr1", 41_101, 42_000, tu_id="TU0001", log2=22.0),
    ]
    tus = [vl.tu("TU0001", "chr1", 40_001, 42_000, "+", ["g3", "g4"])]
    flags = [
        _flag_row("operon_discordance", "PAIR_7", "g3", "g4", 0.02),
        _flag_row("operon_discordance", "PAIR_8", "g4", "g3", 0.05),
        _flag_row("missed_operon", "PAIR_9", "g3", "g4", 0.01),
    ]
    out = vl.build_out(tmp_path, genes=genes, tus=tus, flags=flags)
    loci = select_loci(RunView.from_outdir(out), _cfg())
    # one OP locus (PAIR_7 + PAIR_8 merge in their family) then the NO locus: family order
    assert _ids(loci) == ["PAIR_7", "PAIR_9"]
    op, no = loci
    assert [lo.kind for lo in loci] == ["operon_flag", "operon_flag"]
    assert op.error == pytest.approx(0.02)
    assert no.error == pytest.approx(0.01)
    assert (op.start0, op.end0) == (no.start0, no.end0) == (40_000, 42_000)
    assert "PAIR_9" not in op.label and "novel_operon_candidate" not in op.label
    assert "internal_promoter_supported between g3 and g4 (q = 0.02)" in op.label
    assert "internal_promoter_supported between g4 and g3 (q = 0.05)" in op.label
    assert no.label == "novel_operon_candidate between g3 and g4 (q = 0.01)"


# --- properties ------------------------------------------------------------------------------


def _check_pool(entries: Sequence[Entry], got: Sequence[Key], cfg: GiapConfig) -> None:
    """Pool membership is unchanged: exactly the in-pool entries, once each."""
    want = {(KIND[f], lid) for f, lid, _ in _in_pool(entries, cfg)}
    assert len(got) == len(set(got)) == len(want)
    assert set(got) == want


@PROPS
@given(pools())
def test_property_tiers_come_in_order(entries: list[Entry]) -> None:
    cfg, alphas = _cfg(), _alphas(_cfg())
    got = _select(entries, cfg)
    _check_pool(entries, got, cfg)
    fam, err = _family_of(entries), _error_of(entries)
    tiers = [_tier(fam[k], err[k], alphas) for k in got]
    assert tiers == sorted(tiers), list(zip(got, tiers, strict=True))
    assert all(t == 0 for k, t in zip(got, tiers, strict=True) if fam[k] in FLAG_FAMILIES)


@PROPS
@given(pools())
def test_property_order_within_each_family(entries: list[Entry]) -> None:
    cfg = _cfg()
    got = _select(entries, cfg)
    _check_pool(entries, got, cfg)
    fam, err = _family_of(entries), _error_of(entries)
    for f in FAMILIES:
        sub = [k for k in got if fam[k] == f]
        e_key = [(err[k] is None, err[k] or 0.0, k[1]) for k in sub]
        assert e_key == sorted(e_key), (f, sub)


@PROPS
@given(pools(min_size=3))
def test_property_exact_rotation_rounds(entries: list[Entry]) -> None:
    """Within tiers 0 and 1, the maximal runs of constant family rank are the rounds 1, 2, ...;
    round r holds one locus from each family with at least r loci in the tier, in family order.
    Tier 2 sorts by family, then ID. The ranks are computed here, from the entries."""
    cfg, alphas = _cfg(), _alphas(_cfg())
    got = _select(entries, cfg)
    _check_pool(entries, got, cfg)
    pool = _in_pool(entries, cfg)
    fam, err = _family_of(pool), _error_of(pool)
    ranks = _family_ranks(pool, alphas)
    for tier in (0, 1):
        block = [k for k in got if _tier(fam[k], err[k], alphas) == tier]
        counts = {f: sum(1 for k in block if fam[k] == f) for f in FAMILIES}
        runs = [list(g) for _, g in itertools.groupby(block, key=lambda k: ranks[k])]
        assert [ranks[r[0]] for r in runs] == list(range(1, len(runs) + 1)), (tier, block)
        for r, run in enumerate(runs, start=1):
            assert [fam[k] for k in run] == [f for f in FAMILIES if counts[f] >= r], (tier, r, run)
        # the readable balance check: families with loci left differ by at most 1 per prefix
        for n in range(len(block) + 1):
            seen = {f: sum(1 for k in block[:n] if fam[k] == f) for f in FAMILIES}
            left = [f for f in FAMILIES if seen[f] < counts[f]]
            assert max((seen[f] for f in left), default=0) - min(
                (seen[f] for f in left), default=0) <= 1  # fmt: skip
    tail = [k for k in got if err[k] is None]
    assert tail == sorted(tail, key=lambda k: (FAMILIES.index(fam[k]), k[1]))
    assert got == _expected(entries, cfg)


@PROPS
@given(pools(), st.integers(0, 2**32 - 1))
def test_property_input_row_order_does_not_matter(entries: list[Entry], seed: int) -> None:
    cfg = _cfg()
    got = _select(entries, cfg)
    _check_pool(entries, got, cfg)
    assert _select(entries, cfg, shuffle_seed=seed) == got


@PROPS
@given(pools(), st.sampled_from(FAMILIES), maps)
def test_property_metamorphic_one_family(entries: list[Entry], f: str,
                                         m: tuple[str, float, float]) -> None:  # fmt: skip
    """A strictly increasing map of one family's errors that keeps ``e <= alpha`` and NA leaves
    the order of the whole pool unchanged."""
    cfg, alpha = _cfg(), _alphas(_cfg())[f]
    assume(any(fam == f and e is not None for fam, _, e in entries))
    g = functools.partial(_g, alpha=alpha, m=m)
    assume(_valid_map([e for fam, _, e in entries if fam == f], g, alpha))
    got = _select(entries, cfg)
    _check_pool(entries, got, cfg)
    assert _select(_apply(entries, f, g), cfg) == got


@PROPS
@given(pools(), st.tuples(*([maps] * len(FAMILIES))))
def test_property_metamorphic_every_family_at_once(
    entries: list[Entry], ms: tuple[tuple[str, float, float], ...]
) -> None:
    cfg, alphas = _cfg(), _alphas(_cfg())
    mapped = list(entries)
    for f, m in zip(FAMILIES, ms, strict=True):
        g = functools.partial(_g, alpha=alphas[f], m=m)
        assume(_valid_map([e for fam, _, e in entries if fam == f], g, alphas[f]))
        mapped = _apply(mapped, f, g)
    got = _select(entries, cfg)
    _check_pool(entries, got, cfg)
    assert _select(mapped, cfg) == got


@PROPS
@given(pools(), st.sampled_from(FAMILIES), st.sampled_from((0, 1)),
       st.sampled_from(("below", "above")))  # fmt: skip
def test_property_cross_family_stress(entries: list[Entry], f: str, tier: int, side: str) -> None:
    """Within one tier, push every error of one family below (then above) every error of the
    other families; the order cannot change if no step compares errors across families."""
    cfg, alphas = _cfg(), _alphas(_cfg())
    alpha, eps = alphas[f], 1e-6
    in_tier = [e for fam, _, e in entries if fam == f and _tier(fam, e, alphas) == tier]
    assume(in_tier)

    def g(e: float | None) -> float | None:
        if e is None or _tier(f, e, alphas) != tier:
            return e
        if tier == 0:
            return e * eps if side == "below" else alpha - (alpha - e) * eps
        return alpha + (e - alpha) * eps if side == "below" else 1.0 - (1.0 - e) * eps

    assume(_valid_map([e for fam, _, e in entries if fam == f], g, alpha))
    got = _select(entries, cfg)
    _check_pool(entries, got, cfg)
    assert _select(_apply(entries, f, g), cfg) == got


@PROPS
@given(pools(), st.sampled_from(("context", "operon_flags")), st.floats(0.0, 1.0))
def test_property_changing_an_alpha_that_no_locus_crosses(entries: list[Entry], key: str,
                                                          t: float) -> None:  # fmt: skip
    fams = FLAG_FAMILIES if key == "operon_flags" else ("prophage", "bgc", "island")
    alpha = A_FLAG if key == "operon_flags" else A_CTX
    vals = [e for f, _, e in entries if f in fams and e is not None]
    lo = max((v for v in vals if v <= alpha), default=0.0)
    hi = min((v for v in vals if v > alpha), default=1.0)
    new = lo + t * (hi - lo)
    assume(0.0 < new < 1.0 and lo <= new < hi)
    cfg_new = _cfg(f"stats.alpha.{key}={new!r}")
    assert getattr(cfg_new.stats.alpha, key) == new
    with tempfile.TemporaryDirectory() as tmp:
        view = RunView.from_outdir(_build(Path(tmp), entries))
        got = _keys(select_loci(view, _cfg()))
        _check_pool(entries, got, _cfg())
        assert _keys(select_loci(view, cfg_new)) == got


# --- the index figure ------------------------------------------------------------------------

_INDEX_ENTRIES: list[Entry] = [
    ("prophage", "prophage_1", 0.01), ("bgc", "bgc_1", 0.3), ("island", "ISL_1", 0.005),
    ("operon_discordance", "PAIR_1", 0.02), ("missed_operon", "PAIR_2", 0.03),
    ("prophage", "prophage_2", None), ("island", "ISL_2", None),
]  # fmt: skip
_INDEX_ORDER = ["prophage_1", "ISL_1", "PAIR_1", "PAIR_2", "bgc_1", "prophage_2", "ISL_2"]


def _texts_top_down(inst: Instance) -> list[str]:
    """Every non-empty visible text of the figure, top to bottom (display y)."""
    found = []
    for t in inst.figure.findobj(Text):
        s = t.get_text()
        if s.strip() and t.get_visible():
            y = t.get_transform().transform(t.get_position())[1]
            found.append((-y, s))
    return [s for _, s in sorted(found, key=lambda p: p[0])]


def _index_figure(tmp: Path, entries: Sequence[Entry]) -> list[str]:
    index, _ = _run(_build(tmp, entries), _cfg("report.max_locus_plots=0"))
    return _texts_top_down(index)


def test_index_lines_name_each_family_estimate(tmp_path: Path) -> None:
    texts = _index_figure(tmp_path, _INDEX_ENTRIES)
    lines = [s for s in texts if re.match(r"\d+\. ", s)]
    assert [re.match(r"\d+\. (\S+)", s).group(1) for s in lines] == _INDEX_ORDER  # type: ignore[union-attr]
    fam = {lid: f for f, lid, _ in _INDEX_ENTRIES}
    err = {lid: e for _, lid, e in _INDEX_ENTRIES}
    for rank, (lid, line) in enumerate(zip(_INDEX_ORDER, lines, strict=True), start=1):
        head = f"{rank}. {lid} [{KIND[fam[lid]]}] "
        assert line.startswith(head), line
        rest = line[len(head) :]
        if err[lid] is None:
            assert re.match(r"(?:(?:q|FWER p|decoy-FDR q) )?NA\b", rest), line
        else:
            assert rest.startswith(f"{ESTIMATE[fam[lid]]} {base.sig(err[lid], 3)}"), line
    assert "decoy-FDR" not in lines[_INDEX_ORDER.index("PAIR_1")]
    assert "FWER" not in lines[_INDEX_ORDER.index("prophage_1")]


def test_index_prints_tier_headers_and_the_caption(tmp_path: Path) -> None:
    texts = _index_figure(tmp_path, _INDEX_ENTRIES)
    assert sum(CAPTION in s for s in texts) == 1, texts
    markers: list[str] = []
    for s in texts:
        hit = [h for h in HEADERS if s.strip() == h or s.strip().startswith(h)]
        if hit:
            markers.append(f"H:{hit[0]}")
        elif m := re.match(r"\d+\. (\S+)", s):
            markers.append(m.group(1))
    assert markers == [
        f"H:{HEADERS[0]}", "prophage_1", "ISL_1", "PAIR_1", "PAIR_2",
        f"H:{HEADERS[1]}", "bgc_1",
        f"H:{HEADERS[2]}", "prophage_2", "ISL_2",
    ]  # fmt: skip


def test_index_prints_headers_only_for_non_empty_tiers(tmp_path: Path) -> None:
    texts = _index_figure(tmp_path, [("prophage", "prophage_1", 0.01), ("island", "ISL_1", 0.02)])
    joined = "\n".join(texts)
    assert HEADERS[0] in joined
    assert sum(CAPTION in s for s in texts) == 1
    for h in HEADERS[1:]:
        assert not any(s.strip().startswith(h) for s in texts), texts


def test_index_tsv_follows_the_new_order(tmp_path: Path) -> None:
    index, _ = _run(_build(tmp_path, _INDEX_ENTRIES), _cfg("report.max_locus_plots=2"))
    frame: pd.DataFrame = index.data
    assert list(frame["locus_id"]) == _INDEX_ORDER
    assert list(frame["rank"]) == list(range(1, 8))
    assert [bool(v) for v in frame["rendered"]] == [True, True] + [False] * 5
    assert list(frame["kind"]) == [KIND[f] for f in ("prophage", "island", "operon_discordance",
                                                     "missed_operon", "bgc", "prophage",
                                                     "island")]  # fmt: skip
