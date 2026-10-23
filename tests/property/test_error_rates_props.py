"""M7-10 properties (docs/specs/M7-10-class-error-rates.md rules 1-4): the class error rates of
``giap.variants.error_rates`` against a brute-force oracle, and their invariants.

Interface as pinned in tests/unit/variants/test_error_rates.py. Random inputs: de novo rows on a
few runs and spectrum IDs (some NA), reference rows (q around the 0.01 threshold, some decoys,
some NA spectrum IDs, at most one per (run, spectrum_id) as preflight guarantees, ADR-0005), and
ALT support rows over a small peptide alphabet that includes I/L spellings.

* The per-class ``(n, errors)`` equal the oracle; ``0 <= errors <= n`` and ``rate = errors / n``.
* I/L spelling of the reference peptides never changes a result (rule 3).
* Row order of any input never changes a result.
* Raising ``q_max`` never lowers a class's ``n`` (rule 1 only adds pairs).
* ``annotate_sites`` gives each row exactly its class's rate (or ``unestimated``) and changes no
  other column.
"""

from __future__ import annotations

import importlib
from types import ModuleType
from typing import Any

import pandas as pd
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

UNEST = "unestimated"
PEPTIDES = ("AAISAAK", "AALSAAK", "CCCWCCK", "GGNPEPK", "NNPEPK", "DDDSDDK", "EELEEK")
RUNS = ("r1", "r2")
SIDS = ("s1", "s2", "s3", "s4", "s5", None)
CLASSES = (("substitution", False), ("substitution", True), ("block", True), ("swap", False),
           ("block", False), ("frameshift", False))  # fmt: skip
SETTINGS = settings(max_examples=150, deadline=None,
                    suppress_health_check=[HealthCheck.too_slow])  # fmt: skip


def er() -> ModuleType:
    return importlib.import_module("giap.variants.error_rates")


def fold(p: str) -> str:
    return p.replace("I", "L")


dn_rows = st.lists(
    st.tuples(st.sampled_from(RUNS), st.sampled_from(SIDS), st.sampled_from(PEPTIDES)),
    max_size=25,
)
ref_rows = st.dictionaries(
    st.tuples(st.sampled_from(RUNS), st.sampled_from(SIDS)),
    st.tuples(
        st.sampled_from(PEPTIDES),
        st.sampled_from([0.0001, 0.005, 0.01, 0.0100001, 0.02, 0.2]),
        st.booleans(),
    ),
    max_size=12,
)
sup_rows = st.lists(
    st.tuples(
        st.sampled_from([f"VAR_{i}" for i in range(1, 7)]),
        st.sampled_from(PEPTIDES),
        st.sampled_from(CLASSES),
    ),
    max_size=12,
)


def frames(dn: Any, ref: Any, sup: Any) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    d = pd.DataFrame(
        {"run": pd.Series([r[0] for r in dn], dtype=object),
         "spectrum_id": pd.Series([r[1] for r in dn], dtype=object),
         "peptide": pd.Series([r[2] for r in dn], dtype=object),
         "source": pd.Series(["denovo"] * len(dn), dtype=object)}
    )  # fmt: skip
    items = list(ref.items())
    r = pd.DataFrame(
        {"run": pd.Series([k[0] for k, _ in items], dtype=object),
         "spectrum_id": pd.Series([k[1] for k, _ in items], dtype=object),
         "peptide": pd.Series([v[0] for _, v in items], dtype=object),
         "source": pd.Series(["database"] * len(items), dtype=object),
         "q_value": pd.Series([v[1] for _, v in items], dtype="float64"),
         "is_decoy": pd.Series([v[2] for _, v in items], dtype=bool)}
    )  # fmt: skip
    s = pd.DataFrame(
        {"variant_id": pd.Series([x[0] for x in sup], dtype=object),
         "peptide": pd.Series([x[1] for x in sup], dtype=object),
         "alt_type": pd.Series([x[2][0] for x in sup], dtype=object),
         "mass_equivalent": pd.Series([x[2][1] for x in sup], dtype=bool)}
    )  # fmt: skip
    return d, r, s


def oracle(dn: Any, ref: Any, sup: Any, q_max: float = 0.01) -> dict[tuple[str, bool], tuple]:
    classes: dict[str, set[tuple[str, bool]]] = {}
    for _vid, pep, cls in sup:
        classes.setdefault(pep, set()).add(cls)
    out: dict[tuple[str, bool], list[int]] = {}
    for run, sid, pep in dn:
        if sid is None or pep not in classes:
            continue
        hit = ref.get((run, sid))
        if hit is None or hit[2] or not hit[1] <= q_max:
            continue
        err = fold(pep) != fold(hit[0])
        for cls in classes[pep]:
            c = out.setdefault(cls, [0, 0])
            c[0] += 1
            c[1] += int(err)
    return {k: (v[0], v[1]) for k, v in out.items()}


def measured(d: pd.DataFrame, r: pd.DataFrame, s: pd.DataFrame, **kw: Any) -> dict:
    m = er()
    out: dict[tuple[str, bool], tuple[int, int]] = {}
    for x in m.class_rates(m.paired_observations(d, r, s, **kw)):
        assert 0 <= x.errors <= x.n
        if x.n == 0:
            assert x.rate is None
            continue
        assert abs(x.rate - x.errors / x.n) < 1e-12
        key = (str(x.alt_type), bool(x.mass_equivalent))
        assert key not in out
        out[key] = (int(x.n), int(x.errors))
    return out


@SETTINGS
@given(dn_rows, ref_rows, sup_rows)
def test_rates_equal_the_brute_force_oracle(dn: Any, ref: Any, sup: Any) -> None:
    assert measured(*frames(dn, ref, sup)) == oracle(dn, ref, sup)


@SETTINGS
@given(dn_rows, ref_rows, sup_rows, st.randoms(use_true_random=False))
def test_il_spelling_of_the_reference_is_irrelevant(dn: Any, ref: Any, sup: Any, rnd: Any) -> None:
    def respell(p: str) -> str:
        return "".join(rnd.choice("IL") if c in "IL" else c for c in p)

    ref2 = {k: (respell(v[0]), v[1], v[2]) for k, v in ref.items()}
    assert measured(*frames(dn, ref2, sup)) == measured(*frames(dn, ref, sup))


@SETTINGS
@given(dn_rows, ref_rows, sup_rows, st.randoms(use_true_random=False))
def test_row_order_is_irrelevant(dn: Any, ref: Any, sup: Any, rnd: Any) -> None:
    d, r, s = frames(dn, ref, sup)
    shuffled = [f.sample(frac=1.0, random_state=rnd.randrange(2**31)).reset_index(drop=True)
                for f in (d, r, s)]  # fmt: skip
    assert measured(*shuffled) == measured(d, r, s)


@SETTINGS
@given(dn_rows, ref_rows, sup_rows)
def test_raising_q_max_only_adds_pairs(dn: Any, ref: Any, sup: Any) -> None:
    f = frames(dn, ref, sup)
    low, high = measured(*f, q_max=0.01), measured(*f, q_max=0.05)
    for cls, (n, _e) in low.items():
        assert cls in high and high[cls][0] >= n


@SETTINGS
@given(
    st.lists(st.sampled_from(CLASSES), max_size=10),
    st.dictionaries(
        st.sampled_from(CLASSES), st.tuples(st.integers(0, 30), st.integers(0, 30)), max_size=6
    ),
)
def test_annotate_sites_gives_each_row_its_class_rate(rows: Any, table: Any) -> None:
    m = er()
    rates = [m.ClassRate(t, me, n, min(e, n), (min(e, n) / n) if n else None)
             for (t, me), (n, e) in table.items()]  # fmt: skip
    sites = pd.DataFrame(
        {"variant_id": [f"VAR_{i}" for i in range(1, len(rows) + 1)],
         "alt_type": [c[0] for c in rows],
         "flags": [("mass_equivalent", "snv") if c[1] else ("multi_nt",) for c in rows],
         "n_alt": list(range(len(rows)))}
    )  # fmt: skip
    out = m.annotate_sites(sites, rates)
    assert out["variant_id"].tolist() == sites["variant_id"].tolist()
    assert out["n_alt"].tolist() == sites["n_alt"].tolist()
    for i, (t, me) in enumerate(rows):
        n, e = table.get((t, me), (0, 0))
        e = min(e, n)
        got_r, got_n = out["class_error_rate"].iloc[i], out["class_error_n"].iloc[i]
        if n == 0:
            assert got_r == UNEST and pd.isna(got_n)
        else:
            assert got_r != UNEST and abs(float(got_r) - e / n) < 1e-12 and int(got_n) == n
