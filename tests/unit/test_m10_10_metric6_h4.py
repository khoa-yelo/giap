"""Metric 6 and hypothesis 4 as pure functions of ``validation/pilot/metrics.py`` (M10-10b;
docs/specs/M10-10-a3-run.md, acceptance criteria 8 and 9, examples 9-12, "Edge cases", the
statistician review S1-S6 and the lead decisions).

Frames are written under their schemas and read back with the reader ``metrics.py`` uses
(``helpers.a3_metrics_m10_10b.frame``); Wilson intervals are checked against an independent
textbook implementation (``helpers.a3_metrics_m10_10b.wilson``), never ``gippa``'s ``wilson_ci``.

Interface fixed here (test-author choices where the spec names a function but not its layout;
each is reversible):

* ``novel_summary(df) -> {"n_rows", "by_class": {class: {called, unestimated, not_called}},
  "totals": {called, unestimated, not_called}}``; every class of the schema enum is present.
* ``variant_summary(df) -> {"alt": {alt_type: {high, candidate}}, "mixed",
  "class_error": {"numeric", "unestimated"}}``; ``alt`` counts ``call = ALT`` rows only, every
  type x confidence present; ``class_error`` counts every row of the table (ALT and mixed).
* ``overlap(a2, a3, key) -> {n_a2, n_a3, n_shared, union, jaccard, a3_confirmed_by_a2,
  a2_found_by_a3}``: ``a2``/``a3`` are the frames of one set (already filtered), ``key`` the
  key columns; keys are de-duplicated; ``jaccard`` None with ``reason = "both_empty"`` when the
  union is 0; the shares are M10-05's ``share`` dicts.
* ``metric6(tables, compare_stability) -> dict``: ``tables`` maps an arm label to
  ``{"novel_orfs": df | None, "variants": df | None}``; without an ``"A3"`` entry the overlap
  is ``{"reason": "a3_not_run"}`` (lead decision after the M10-15a review). Output keys: one
  block per arm (``{"novel_orfs": …, "variants": …}``, an absent table ``{"reason":
  "missing:<table>.tsv"}``), ``overlap_a2_a3`` = ``{"novel_orfs": {"called": …, "any_row": …,
  "orf_id_jaccard": …}, "alt_sites": {"alt": …, "high": …}}`` and ``warnings`` (dicts with
  ``code``). ``called`` also carries ``class_agree`` and ``by_class_a3`` ({class: {n, shared}}).
  The ALT cross-check against the (A2, A3) ``alt_sites`` row raises ``Inconsistent``; the ORF
  cross-check warns (``novel_orf_key_differs``, ``duplicate_orf_key``).
* ``denovo_queries(P3, cfg3) -> DataFrame`` (``peptide``, ``score``; one row per de novo
  query, definition 1).
* ``h4_bin(n_b, k_b, genome, n_reportable) -> {"L_b", "L_b_adj", "U_b", "relation_95",
  "relation_adj"}`` (definition 7 for one bin; ``genome`` = {precision, ci_low, ci_high, …}).
  Definitions 7 and 8 follow "Lead decision: the H4 fallback (2026-10-25)": ``relation_adj``
  compares ``genome.ci_high`` with ``L_b_adj.ci_low``; ``U_b`` is descriptive only;
  ``not_supported`` needs two consecutive reportable bins below.
* ``hypothesis4(queries, d2, mapping_fdr, *, searchable=None, n_denovo_logged=None) -> dict``
  with ``verdict``, ``B``, ``bins`` (a list in ``score_bin`` order; per bin ``score_bin``,
  ``n_b``, ``k_b``, ``L_b``, ``L_b_adj``, ``U_b``, ``genome``, ``reportable``, ``relation_95``,
  ``relation_adj``, ``nested_share``, ``L_b_searchable``), ``mean_genome_precision``,
  ``pooled_L``, ``caveats`` (strings) and ``warnings``. ``hypothesis4`` folds I→L on both
  sides itself. ``n_denovo_logged`` is the "(M de novo)" count of the placement log line
  (S6): a different count raises ``Inconsistent``; None gives ``denovo_count_unchecked``.
"""

from __future__ import annotations

import itertools
import re
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.a3_metrics_m10_10b import (
    ALT_TYPES,
    BUCKETS,
    C1,
    C2,
    CONFIDENCE,
    NOVEL_CLASSES,
    PLACEMENT_ROW,
    VERDICTS,
    dn_row,
    equal_bins,
    example9_novel,
    example10_variants,
    example_stability,
    fdr_row,
    frame,
    orf,
    queries_frame,
    synthetic_peptides,
    upper,
    variant,
    wilson,
    write_denovo_table,
    z_ratio,
)
from helpers.pilot_metrics_m10 import load_metrics
from hypothesis import HealthCheck, Phase, assume, given, settings
from hypothesis import strategies as st

ORF_KEY = ("contig", "strand", "end")
# The whole-rule properties call hypothesis4 on up to 1,500 peptides per example; shrinking a
# failure takes minutes, so they report the first failing example (inputs in the message).
NO_SHRINK = (Phase.explicit, Phase.reuse, Phase.generate)
ALT_KEY = ("contig", "codon_start", "strand", "alt_aa")


@pytest.fixture(scope="module")
def M() -> Any:
    return load_metrics()


def fn(M: Any, name: str) -> Any:
    """The function ``name`` of metrics.py; a clear failure while it does not exist."""
    f = getattr(M, name, None)
    if f is None:
        pytest.fail(f"validation/pilot/metrics.py has no {name} (M10-10b)")
    return f


def codes(warnings: Any) -> list[str]:
    return [w if isinstance(w, str) else str(w.get("code")) for w in warnings or []]


def val(x: Any) -> Any:
    """A share dict's ``value``, or the number itself."""
    return x["value"] if isinstance(x, dict) else x


def approx4(x: float) -> Any:
    return pytest.approx(x, abs=5e-5)


def _metric6(
    M: Any,
    tmp: Path,
    novel: dict[str, list[Any]],
    variants: dict[str, list[Any]],
    stability: list[dict[str, object]],
) -> dict[str, Any]:
    tables = {
        arm: {
            "novel_orfs": frame("novel_orfs", novel[arm], tmp),
            "variants": frame("variants", variants[arm], tmp),
        }
        for arm in novel
    }
    out: dict[str, Any] = fn(M, "metric6")(tables, frame("compare_stability", stability, tmp))
    return out


def _ex9_10(M: Any, tmp: Path, **kw: Any) -> dict[str, Any]:
    n2, n3 = example9_novel()
    v2, v3 = example10_variants()
    novel = {"A1": [], "A2": n2, "A3": n3}
    variants = {"A1": [], "A2": v2, "A3": v3}
    return _metric6(M, tmp, novel, variants, kw.pop("stability", example_stability()))


# ------------------------------------------------------------------ per-arm summaries


def test_novel_summary_buckets_by_class(M: Any, tmp_path: Path) -> None:
    rows = [
        orf(C1, "+", 100, "intergenic", called=True, class_fdr="0.01"),
        orf(C1, "+", 400, "intergenic", called=True, class_fdr="unestimated"),  # called wins
        orf(C1, "+", 700, "intergenic", called=False, class_fdr="unestimated"),
        orf(C1, "-", 900, "sorf", called=False, class_fdr="0.2"),
        orf(C1, "-", 1200, "antisense", called=False, class_fdr="unestimated"),
        orf(C2, "+", 50, "nterm_extension", called=True, class_fdr="0.03"),
    ]
    out = fn(M, "novel_summary")(frame("novel_orfs", rows, tmp_path))
    assert out["n_rows"] == 6
    assert set(out["by_class"]) == set(NOVEL_CLASSES)
    for cls in NOVEL_CLASSES:
        assert set(out["by_class"][cls]) == set(BUCKETS), cls
    assert out["by_class"]["intergenic"] == {"called": 2, "unestimated": 1, "not_called": 0}
    assert out["by_class"]["sorf"] == {"called": 0, "unestimated": 0, "not_called": 1}
    assert out["by_class"]["antisense"] == {"called": 0, "unestimated": 1, "not_called": 0}
    assert out["by_class"]["nterm_extension"] == {"called": 1, "unestimated": 0, "not_called": 0}
    assert out["by_class"]["pseudogene"] == {"called": 0, "unestimated": 0, "not_called": 0}
    assert out["totals"] == {"called": 3, "unestimated": 2, "not_called": 1}


def test_novel_summary_empty_table(M: Any, tmp_path: Path) -> None:
    out = fn(M, "novel_summary")(frame("novel_orfs", [], tmp_path))
    assert out["n_rows"] == 0
    assert out["totals"] == dict.fromkeys(BUCKETS, 0)
    assert all(v == dict.fromkeys(BUCKETS, 0) for v in out["by_class"].values())


def test_variant_summary_alt_by_type_and_confidence(M: Any, tmp_path: Path) -> None:
    rows = [
        variant(C1, 10, "+", "N", i=1, alt_type="substitution", confidence="high"),
        variant(
            C1,
            13,
            "+",
            "E",
            i=2,
            alt_type="substitution",
            confidence="candidate",
            class_error_rate="0.02",
        ),
        variant(C1, 40, "-", "W", i=3, alt_type="swap", confidence="high"),
        variant(C1, 70, "+", "Q", i=4, alt_type="frameshift", confidence="candidate"),
        variant(
            C1, 90, "+", "S", i=5, call="mixed", confidence="candidate", class_error_rate="0.1"
        ),
    ]
    out = fn(M, "variant_summary")(frame("variants", rows, tmp_path))
    assert set(out["alt"]) == set(ALT_TYPES)
    for t in ALT_TYPES:
        assert set(out["alt"][t]) == set(CONFIDENCE), t
    assert out["alt"]["substitution"] == {"high": 1, "candidate": 1}
    assert out["alt"]["swap"] == {"high": 1, "candidate": 0}
    assert out["alt"]["block"] == {"high": 0, "candidate": 0}
    assert out["alt"]["frameshift"] == {"high": 0, "candidate": 1}
    assert out["mixed"] == 1  # the mixed row is not in "alt"
    assert out["class_error"] == {"numeric": 2, "unestimated": 3}


def test_metric6_per_arm_blocks_and_missing_table(M: Any, tmp_path: Path) -> None:
    """A table absent in an arm (module off): ``reason = missing:<table>`` for that arm only."""
    n2, n3 = example9_novel()
    tables = {
        "A1": {"novel_orfs": None, "variants": frame("variants", [], tmp_path)},
        "A2": {
            "novel_orfs": frame("novel_orfs", n2, tmp_path),
            "variants": frame("variants", [], tmp_path),
        },
        "A3": {
            "novel_orfs": frame("novel_orfs", n3, tmp_path),
            "variants": frame("variants", [], tmp_path),
        },
    }
    stab_rows = example_stability(alt=[0, 0, 0, None, "both_empty"])
    out = fn(M, "metric6")(tables, frame("compare_stability", stab_rows, tmp_path))
    assert out["A1"]["novel_orfs"]["reason"] == "missing:novel_orfs.tsv"
    assert out["A1"]["variants"]["alt"]["substitution"] == {"high": 0, "candidate": 0}
    assert out["A2"]["novel_orfs"]["n_rows"] == 2 and out["A3"]["novel_orfs"]["n_rows"] == 3


def test_metric6_without_a3_is_a3_not_run(M: Any, tmp_path: Path) -> None:
    """Lead decision (M10-15a review): no ``--a3`` → the overlap block carries ``a3_not_run``."""
    n2, _ = example9_novel()
    v2, _ = example10_variants()
    tables = {
        arm: {
            "novel_orfs": frame("novel_orfs", n2, tmp_path),
            "variants": frame("variants", v2, tmp_path),
        }
        for arm in ("A1", "A2")
    }
    out = fn(M, "metric6")(tables, None)
    assert out["overlap_a2_a3"]["reason"] == "a3_not_run"
    assert out["A1"]["novel_orfs"]["n_rows"] == 2 and out["A2"]["variants"]["mixed"] == 0
    if "A3" in out:
        assert out["A3"].get("reason") == "a3_not_run"


def test_a3_numeric_class_error_rate_warns(M: Any, tmp_path: Path) -> None:
    """§7.5.4: A3 has no spectrum pairing, so every ALT row should be ``unestimated``."""
    n2, n3 = example9_novel()
    v2, v3 = example10_variants()
    v3[1]["class_error_rate"] = "0.05"
    out = _metric6(
        M,
        tmp_path,
        {"A1": [], "A2": n2, "A3": n3},
        {"A1": [], "A2": v2, "A3": v3},
        example_stability(),
    )
    assert "a3_class_error_estimated" in codes(out["warnings"])
    assert out["A3"]["variants"]["class_error"] == {"numeric": 1, "unestimated": 1}


def test_a3_all_unestimated_no_class_error_warning(M: Any, tmp_path: Path) -> None:
    out = _ex9_10(M, tmp_path)
    assert "a3_class_error_estimated" not in codes(out["warnings"])


# ------------------------------------------------------------------ overlap: examples 9 and 10


def test_example_9_novel_orf_overlap(M: Any, tmp_path: Path) -> None:
    out = _ex9_10(M, tmp_path)
    c = out["overlap_a2_a3"]["novel_orfs"]["called"]
    assert (c["n_a2"], c["n_a3"], c["n_shared"], c["union"]) == (2, 3, 1, 4)
    assert c["jaccard"] == pytest.approx(0.25)
    a3c, a2f = c["a3_confirmed_by_a2"], c["a2_found_by_a3"]
    assert a3c["value"] == pytest.approx(1 / 3) and (a3c["k"], a3c["n"]) == (1, 3)
    assert (a3c["ci_low"], a3c["ci_high"]) == (approx4(0.0615), approx4(0.7923))
    assert (a3c["ci_low"], a3c["ci_high"]) == pytest.approx(wilson(1, 3), abs=1e-9)
    assert a2f["value"] == pytest.approx(0.5)
    assert (a2f["ci_low"], a2f["ci_high"]) == (approx4(0.0945), approx4(0.9055))
    assert (a2f["ci_low"], a2f["ci_high"]) == pytest.approx(wilson(1, 2), abs=1e-9)
    assert c["class_agree"] == 0  # intergenic in A2, nterm_extension in A3
    by = c["by_class_a3"]
    assert by["nterm_extension"] == {"n": 1, "shared": 1}
    assert by["intergenic"] == {"n": 1, "shared": 0}
    assert by["antisense"] == {"n": 1, "shared": 0}
    assert all(
        v["n"] == 0
        for k, v in by.items()
        if k not in ("nterm_extension", "intergenic", "antisense")
    )


def test_example_9_orf_id_key_differs_warns(M: Any, tmp_path: Path) -> None:
    """compare keys on ``orf_id`` (0 of 5 shared) while the stop-codon key shares one ORF: a
    warning with both values, not an exit."""
    out = _ex9_10(M, tmp_path)
    nov = out["overlap_a2_a3"]["novel_orfs"]
    assert nov["orf_id_jaccard"] == pytest.approx(0.0)
    w = [
        x
        for x in out["warnings"]
        if isinstance(x, dict) and x.get("code") == "novel_orf_key_differs"
    ]
    assert len(w) == 1
    nums = [float(v) for k, v in w[0].items() if k != "code" and isinstance(v, int | float)]
    assert any(x == pytest.approx(0.25) for x in nums), w[0]  # the stop-codon key's Jaccard
    assert any(x == pytest.approx(0.0) for x in nums), w[0]  # compare's orf_id Jaccard


def test_orf_cross_check_count_difference_is_duplicate_orf_key(M: Any, tmp_path: Path) -> None:
    """Per-arm counts must agree with compare's (``n_a2 = n_a``, ``n_a3 = n_b``); a difference
    is the warning ``duplicate_orf_key``, never an exit."""
    rows = example_stability()
    nov = next(r for r in rows if r["family"] == "novel_orfs" and r["run_a"] == "A2")
    nov["n_b"] = 4  # compare saw 4 called ORF ids in A3
    out = _ex9_10(M, tmp_path, stability=rows)
    assert "duplicate_orf_key" in codes(out["warnings"])


def test_duplicate_orf_key_in_one_arm_counted_once(M: Any, tmp_path: Path) -> None:
    n2, n3 = example9_novel()
    n3.append(orf(C2, "+", 3011, "intergenic", lowest=2650))  # same stop codon, other orf_id
    v2, v3 = example10_variants()
    out = _metric6(
        M,
        tmp_path,
        {"A1": [], "A2": n2, "A3": n3},
        {"A1": [], "A2": v2, "A3": v3},
        example_stability(),
    )
    c = out["overlap_a2_a3"]["novel_orfs"]["called"]
    assert (c["n_a3"], c["n_shared"]) == (3, 1)
    assert "duplicate_orf_key" in codes(out["warnings"])


def test_any_row_set_includes_uncalled_rows(M: Any, tmp_path: Path) -> None:
    n2, n3 = example9_novel()
    n2.append(orf(C2, "+", 3011, "intergenic", called=False, class_fdr="unestimated"))
    v2, v3 = example10_variants()
    out = _metric6(
        M,
        tmp_path,
        {"A1": [], "A2": n2, "A3": n3},
        {"A1": [], "A2": v2, "A3": v3},
        example_stability(),
    )
    nov = out["overlap_a2_a3"]["novel_orfs"]
    assert (nov["called"]["n_a2"], nov["called"]["n_shared"]) == (2, 1)
    a = nov["any_row"]
    assert (a["n_a2"], a["n_a3"], a["n_shared"], a["union"]) == (3, 3, 2, 4)
    assert a["jaccard"] == pytest.approx(0.5)


def test_example_10_alt_sites_match_compare(M: Any, tmp_path: Path) -> None:
    out = _ex9_10(M, tmp_path)
    alt = out["overlap_a2_a3"]["alt_sites"]["alt"]
    assert (alt["n_a2"], alt["n_a3"], alt["n_shared"], alt["union"]) == (1, 2, 1, 2)
    assert alt["jaccard"] == pytest.approx(0.5)
    hi = out["overlap_a2_a3"]["alt_sites"]["high"]
    assert (hi["n_a2"], hi["n_a3"], hi["n_shared"]) == (1, 2, 1)


def test_example_10_doctored_compare_jaccard_is_inconsistent(M: Any, tmp_path: Path) -> None:
    with pytest.raises(M.Inconsistent):
        _ex9_10(M, tmp_path, stability=example_stability(alt=[1, 2, 1, 1.0, None]))


def test_doctored_compare_counts_are_inconsistent(M: Any, tmp_path: Path) -> None:
    for alt in ([1, 2, 0, 0.5, None], [2, 2, 1, 0.5, None], [1, 3, 1, 0.5, None]):
        with pytest.raises(M.Inconsistent):
            _ex9_10(M, tmp_path, stability=example_stability(alt=alt))


def test_both_empty_must_agree_with_compare_flag(M: Any, tmp_path: Path) -> None:
    n2, n3 = example9_novel()
    novel = {"A1": [], "A2": n2, "A3": n3}
    variants: dict[str, list[Any]] = {"A1": [], "A2": [], "A3": []}
    ok = _metric6(
        M, tmp_path, novel, variants, example_stability(alt=[0, 0, 0, None, "both_empty"])
    )
    alt = ok["overlap_a2_a3"]["alt_sites"]["alt"]
    assert alt["jaccard"] is None and alt["reason"] == "both_empty"
    assert alt["a3_confirmed_by_a2"]["value"] is None
    assert alt["a3_confirmed_by_a2"]["reason"] == "empty_denominator"
    with pytest.raises(M.Inconsistent):  # compare says both empty, metric 6 sees one site
        _metric6(
            M,
            tmp_path,
            novel,
            {"A1": [], "A2": [], "A3": example10_variants()[1][:1]},
            example_stability(alt=[0, 0, 0, None, "both_empty"]),
        )


def test_alt_jaccard_compared_after_6g_formatting(M: Any, tmp_path: Path) -> None:
    """S5: compare writes ``%.6g`` (1/3 → 0.333333); no raw float ``==``."""
    a2 = [variant(C1, 10, "+", "N", i=1), variant(C1, 40, "+", "E", i=2)]
    a3 = [variant(C1, 10, "+", "N", i=1), variant(C1, 70, "-", "Q", i=2)]
    n2, n3 = example9_novel()
    out = _metric6(
        M,
        tmp_path,
        {"A1": [], "A2": n2, "A3": n3},
        {"A1": [], "A2": a2, "A3": a3},
        example_stability(alt=[2, 2, 1, 1 / 3, None]),
    )
    assert out["overlap_a2_a3"]["alt_sites"]["alt"]["jaccard"] == pytest.approx(1 / 3)


def test_high_set_is_each_arms_own_high_rows(M: Any, tmp_path: Path) -> None:
    """Statistician nit: ``high`` = A2's ``high`` rows vs A3's ``high`` rows; a site high in A2
    and candidate in A3 is in A2's set only."""
    a2 = [variant(C1, 10, "+", "N", i=1), variant(C1, 40, "+", "E", i=2)]
    a3 = [variant(C1, 10, "+", "N", i=1), variant(C1, 40, "+", "E", i=2, confidence="candidate")]
    a3.append(variant(C1, 90, "+", "S", i=3, call="mixed"))  # mixed: in neither set
    n2, n3 = example9_novel()
    out = _metric6(
        M,
        tmp_path,
        {"A1": [], "A2": n2, "A3": n3},
        {"A1": [], "A2": a2, "A3": a3},
        example_stability(alt=[2, 2, 2, 1.0, None]),
    )
    hi = out["overlap_a2_a3"]["alt_sites"]["high"]
    assert (hi["n_a2"], hi["n_a3"], hi["n_shared"]) == (2, 1, 1)
    assert hi["jaccard"] == pytest.approx(0.5)


def test_edge_case_fixtures_origin_linear_minus_strand(M: Any, tmp_path: Path) -> None:
    """The "Edge cases": an origin-spanning ORF in both arms (its ``end`` as written, small, with
    ``orf_id``s near the contig end) is shared; an ORF on a linear contig counts; two ORFs on
    opposite strands with one ``end`` are different keys; a minus-strand ALT keys on
    ``codon_start`` and ``strand``."""
    a2 = [
        orf(C1, "+", 2, "intergenic", lowest=99701),
        orf("pLIN", "+", 800, "sorf"),
        orf(C1, "+", 5000, "intergenic"),
    ]
    a3 = [
        orf(C1, "+", 2, "intergenic", lowest=99650),
        orf(C1, "-", 5000, "antisense"),
        orf("pLIN", "+", 800, "sorf", lowest=650),
    ]
    v2 = [variant(C1, 7702, "-", "D", i=1), variant(C1, 7702, "+", "D", i=2)]
    v3 = [variant(C1, 7702, "-", "D", i=1), variant(C1, 7702, "-", "E", i=2)]
    out = _metric6(
        M,
        tmp_path,
        {"A1": [], "A2": a2, "A3": a3},
        {"A1": [], "A2": v2, "A3": v3},
        example_stability(novel_j=0.0, alt=[2, 2, 1, 1 / 3, None]),
    )
    c = out["overlap_a2_a3"]["novel_orfs"]["called"]
    assert (c["n_a2"], c["n_a3"], c["n_shared"], c["union"]) == (3, 3, 2, 4)
    alt = out["overlap_a2_a3"]["alt_sites"]["alt"]
    assert (alt["n_shared"], alt["union"]) == (1, 3)  # two ALT residues at one codon: two keys


# ------------------------------------------------------------------ overlap properties

_key = st.tuples(
    st.sampled_from(["c1", "c2"]),
    st.sampled_from(["+", "-"]),
    st.integers(min_value=1, max_value=12),
)


def _orf_frame(M: Any, keys: set[tuple[str, str, int]], tmp: Path) -> pd.DataFrame:
    return frame("novel_orfs", [orf(c, s, e) for c, s, e in sorted(keys)], tmp)


@settings(
    max_examples=40, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(a=st.sets(_key, max_size=8), b=st.sets(_key, max_size=8))
def test_overlap_set_identities(M: Any, tmp_path: Path, a: set[Any], b: set[Any]) -> None:
    ov = fn(M, "overlap")
    x = ov(_orf_frame(M, a, tmp_path), _orf_frame(M, b, tmp_path), ORF_KEY)
    assert (x["n_a2"], x["n_a3"], x["n_shared"]) == (len(a), len(b), len(a & b))
    assert x["union"] == len(a | b) == x["n_a2"] + x["n_a3"] - x["n_shared"]
    if a | b:
        assert x["jaccard"] == pytest.approx(len(a & b) / len(a | b))
    else:
        assert x["jaccard"] is None and x["reason"] == "both_empty"
    for name, n in (("a3_confirmed_by_a2", len(b)), ("a2_found_by_a3", len(a))):
        s = x[name]
        if n == 0:
            assert s["value"] is None and s["reason"] == "empty_denominator"
        else:
            lo, hi = wilson(len(a & b), n)
            assert s["value"] == pytest.approx(len(a & b) / n)
            assert (s["ci_low"], s["ci_high"]) == pytest.approx((lo, hi), abs=1e-9)
    y = ov(_orf_frame(M, b, tmp_path), _orf_frame(M, a, tmp_path), ORF_KEY)  # swap the arms
    assert (y["n_a2"], y["n_a3"], y["jaccard"]) == (x["n_a3"], x["n_a2"], x["jaccard"])
    assert val(y["a3_confirmed_by_a2"]) == val(x["a2_found_by_a3"])


@settings(
    max_examples=30, deadline=None, suppress_health_check=[HealthCheck.function_scoped_fixture]
)
@given(
    rows=st.lists(
        st.tuples(
            st.sampled_from(NOVEL_CLASSES),
            st.booleans(),
            st.sampled_from(["unestimated", "0.01", "0.5"]),
        ),
        max_size=15,
    )
)
def test_novel_summary_buckets_partition_rows(M: Any, tmp_path: Path, rows: list[Any]) -> None:
    df = frame(
        "novel_orfs",
        [orf(C1, "+", 10 + 3 * i, c, called=k, class_fdr=f) for i, (c, k, f) in enumerate(rows)],
        tmp_path,
    )
    out = fn(M, "novel_summary")(df)
    assert out["n_rows"] == len(rows) == sum(out["totals"].values())
    for b in BUCKETS:
        assert out["totals"][b] == sum(v[b] for v in out["by_class"].values())
    assert out["totals"]["called"] == sum(k for _, k, _ in rows)


# ------------------------------------------------------------------ hypothesis 4: definition 7
#
# Definition 7 as amended by "Lead decision: the H4 fallback (2026-10-25)": ``relation_adj`` is
# ``below_lower_bound`` when ``genome.ci_high < L_b^adj.ci_low`` (the unrescaled 95 % genome
# upper limit against the two-sided Wilson at 1 - 0.05/B); ``U_b`` is still written, as a
# descriptive column, and no longer enters any relation.


def test_example_11_bin_7(M: Any) -> None:
    b = fn(M, "h4_bin")(200, 120, {"precision": 0.71, "ci_low": 0.64, "ci_high": 0.78}, 10)
    assert b["L_b"]["value"] == pytest.approx(0.6)
    assert (b["L_b"]["ci_low"], b["L_b"]["ci_high"]) == (approx4(0.5308), approx4(0.6654))
    assert (b["L_b_adj"]["ci_low"], b["L_b_adj"]["ci_high"]) == (approx4(0.5008), approx4(0.6917))
    assert (b["L_b_adj"]["ci_low"], b["L_b_adj"]["ci_high"]) == pytest.approx(
        wilson(120, 200, 0.995), abs=1e-9
    )
    assert b["U_b"] == approx4(0.8103) and b["U_b"] == pytest.approx(upper(0.71, 0.78, 10))
    assert b["relation_95"] == "consistent" and b["relation_adj"] == "consistent"  # 0.78 >= 0.5008


def test_example_11_bin_3_below_at_95_but_not_adjusted(M: Any) -> None:
    b = fn(M, "h4_bin")(150, 90, {"precision": 0.41, "ci_low": 0.33, "ci_high": 0.49}, 10)
    assert (b["L_b"]["ci_low"], b["L_b"]["ci_high"]) == (approx4(0.5200), approx4(0.6750))
    assert (b["L_b_adj"]["ci_low"], b["L_b_adj"]["ci_high"]) == (approx4(0.4855), approx4(0.7046))
    assert b["U_b"] == approx4(0.5246)  # still written (descriptive)
    assert b["relation_95"] == "below_lower_bound"  # 0.49 < 0.52
    assert b["relation_adj"] == "consistent"  # adopted fallback: ci_high 0.49 >= 0.4855


def test_example_11_bin_5_below_adjusted(M: Any) -> None:
    b = fn(M, "h4_bin")(150, 90, {"precision": 0.31, "ci_low": 0.24, "ci_high": 0.37}, 10)
    assert b["U_b"] == approx4(0.3959)
    assert b["relation_95"] == "below_lower_bound"
    assert b["relation_adj"] == "below_lower_bound"  # ci_high 0.37 < 0.4855


def test_relation_adj_compares_ci_high_not_u_b(M: Any) -> None:
    """The fallback's discriminating case: ci_high 0.48 < L_b^adj.ci_low 0.4855 <= U_b 0.5146.
    The one-bin rule (U_b against L_b^adj) read ``consistent``; the adopted rule reads below."""
    g = {"precision": 0.40, "ci_low": 0.32, "ci_high": 0.48}
    b = fn(M, "h4_bin")(150, 90, g, 10)
    lo_adj = wilson(90, 150, 1 - 0.05 / 10)[0]
    assert 0.48 < lo_adj <= upper(0.40, 0.48, 10)  # the oracle agrees the case discriminates
    assert b["U_b"] == pytest.approx(upper(0.40, 0.48, 10))  # U_b still computed and written
    assert b["relation_adj"] == "below_lower_bound"


def test_z_ratio_of_example_11() -> None:
    assert z_ratio(10) == approx4(1.4322)  # the oracle itself agrees with the spec's digits


def test_b_equal_1_relations_coincide(M: Any) -> None:
    g = {"precision": 0.41, "ci_low": 0.33, "ci_high": 0.49}
    b = fn(M, "h4_bin")(150, 90, g, 1)
    assert b["U_b"] == pytest.approx(0.49)
    assert (b["L_b_adj"]["ci_low"], b["L_b_adj"]["ci_high"]) == pytest.approx(
        (b["L_b"]["ci_low"], b["L_b"]["ci_high"])
    )
    assert b["relation_95"] == b["relation_adj"] == "below_lower_bound"


def test_u_b_is_clipped_at_1(M: Any) -> None:
    b = fn(M, "h4_bin")(50, 10, {"precision": 0.95, "ci_low": 0.85, "ci_high": 1.0}, 10)
    assert b["U_b"] == 1.0


_genome = st.tuples(
    st.floats(min_value=0.0, max_value=1.0), st.floats(min_value=0.0, max_value=1.0)
).map(sorted)


@settings(max_examples=200, deadline=None)
@given(
    n=st.integers(min_value=20, max_value=3000),
    frac=st.floats(min_value=0, max_value=1),
    g=_genome,
    B=st.integers(min_value=1, max_value=30),
)
def test_relation_adj_implies_relation_95(
    M: Any, n: int, frac: float, g: list[float], B: int
) -> None:
    """Attack plan: ``relation_95 ⊇ relation_adj`` (a bin below at the adjusted level is below
    at 95 %: L_b^adj.ci_low <= L_b.ci_low); U_b >= ci_high and U_b = the descriptive rescaling;
    ``relation_adj`` from ``ci_high`` against ``L_b^adj.ci_low`` (the adopted fallback)."""
    k = round(frac * n)
    p, hi = g
    b = fn(M, "h4_bin")(n, k, {"precision": p, "ci_low": max(0.0, p - 0.05), "ci_high": hi}, B)
    assert b["U_b"] >= hi - 1e-12 and b["U_b"] <= 1.0
    assert b["U_b"] == pytest.approx(upper(p, hi, B), abs=1e-9)
    lo95, _ = wilson(k, n)
    lo_adj, _ = wilson(k, n, 1 - 0.05 / B)
    assert b["L_b"]["ci_low"] == pytest.approx(lo95, abs=1e-9)
    assert b["L_b_adj"]["ci_low"] == pytest.approx(lo_adj, abs=1e-9)
    assert b["relation_95"] == ("below_lower_bound" if hi < lo95 else "consistent")
    assert b["relation_adj"] == ("below_lower_bound" if hi < lo_adj else "consistent")
    if b["relation_adj"] == "below_lower_bound":
        assert b["relation_95"] == "below_lower_bound"


# ------------------------------------------------------------------ hypothesis 4: the whole rule
#
# Definition 8 as amended: ``not_supported`` needs two consecutive reportable bins (adjacent in
# ``score_bin`` order after the non-reportable bins are removed; labels skipped by merged ties do
# not break adjacency) both ``below_lower_bound``; with B = 1 it is impossible.

LOW = (0.31, 0.24, 0.37)  # example 11's bin 5: below L_b^adj (0.4855 at n = 150, k = 90)


def _ex11_frame(
    below: tuple[int, ...] = (), unreported: tuple[int, ...] = ()
) -> tuple[pd.DataFrame, set[str], list[dict[str, object]]]:
    """Ten bins of 150 peptides, k = 90 in each (example 11's bin 3 and 5 counts); genome rows
    as example 11 (bin 7 0.71 [0.64, 0.78], bin 3 0.41 [0.33, 0.49]; the bins of ``below`` 0.31
    [0.24, 0.37] as the variant bin 5; the others 0.9 [0.85, 0.95]); the bins of ``unreported``
    carry no precision (flag ``unreported``)."""
    peps, scores, by_bin = equal_bins(10, 150)
    d2 = {p for b in by_bin for p in b[:90]}
    g: dict[int, tuple[float | None, ...]] = {b: (0.9, 0.85, 0.95) for b in range(10)}
    g[7], g[3] = (0.71, 0.64, 0.78), (0.41, 0.33, 0.49)
    for b in below:
        g[b] = LOW
    rows = [
        fdr_row(b, None, None, None, flags="unreported") if b in unreported else fdr_row(b, *g[b])
        for b in range(10)
    ]
    return queries_frame(peps, scores), d2, [PLACEMENT_ROW, *rows]


def _h4(
    M: Any, tmp: Path, q: pd.DataFrame, d2: Any, fdr: list[dict[str, object]], **kw: Any
) -> dict[str, Any]:
    kw.setdefault("n_denovo_logged", len(q))
    out: dict[str, Any] = fn(M, "hypothesis4")(q, d2, frame("mapping_fdr", fdr, tmp), **kw)
    return out


def _below(out: dict[str, Any]) -> list[int]:
    return [b["score_bin"] for b in out["bins"] if b["relation_adj"] == "below_lower_bound"]


def test_example_11_verdict_consistent_with_lower_bound(M: Any, tmp_path: Path) -> None:
    q, d2, fdr = _ex11_frame()
    out = _h4(M, tmp_path, q, d2, fdr)
    assert out["B"] == 10
    assert [b["score_bin"] for b in out["bins"]] == list(range(10))
    b3, b7 = out["bins"][3], out["bins"][7]
    assert (b3["n_b"], b3["k_b"]) == (150, 90)
    assert b3["relation_95"] == "below_lower_bound" and b3["relation_adj"] == "consistent"
    assert b7["relation_95"] == b7["relation_adj"] == "consistent"
    assert b3["U_b"] == pytest.approx(upper(0.41, 0.49, 10))  # written, descriptive
    assert b3["genome"]["precision"] == pytest.approx(0.41)
    assert b3["genome"]["ci_high"] == pytest.approx(0.49)
    assert out["verdict"] == "consistent_with_lower_bound"


def test_example_11_variant_bin_5_alone_is_consistent(M: Any, tmp_path: Path) -> None:
    """Under the adopted fallback, bin 5 below on its own reads ``consistent_with_lower_bound``
    (the one-bin rule read ``not_supported``)."""
    q, d2, fdr = _ex11_frame(below=(5,))
    out = _h4(M, tmp_path, q, d2, fdr)
    assert _below(out) == [5]
    assert out["verdict"] == "consistent_with_lower_bound"


def test_example_11_bins_5_and_6_not_supported(M: Any, tmp_path: Path) -> None:
    """Example 11: bin 6 also as bin 5 (two consecutive reportable bins) -> ``not_supported``."""
    q, d2, fdr = _ex11_frame(below=(5, 6))
    out = _h4(M, tmp_path, q, d2, fdr)
    assert _below(out) == [5, 6]
    assert out["verdict"] == "not_supported"


def test_example_11_bins_5_and_7_not_adjacent(M: Any, tmp_path: Path) -> None:
    """Example 11: bins 5 and 7 below with bin 6 consistent -> ``consistent_with_lower_bound``."""
    q, d2, fdr = _ex11_frame(below=(5, 7))
    out = _h4(M, tmp_path, q, d2, fdr)
    assert _below(out) == [5, 7] and out["bins"][6]["relation_adj"] == "consistent"
    assert out["verdict"] == "consistent_with_lower_bound"


@pytest.mark.parametrize("below", [(0, 9), (0, 2, 4, 6, 8), (1, 3, 5, 7, 9)])
def test_below_bins_never_adjacent_are_consistent(
    M: Any, tmp_path: Path, below: tuple[int, ...]
) -> None:
    """Two (or five) reportable bins below, no two consecutive: no refutation."""
    q, d2, fdr = _ex11_frame(below=below)
    out = _h4(M, tmp_path, q, d2, fdr)
    assert _below(out) == list(below)
    assert out["verdict"] == "consistent_with_lower_bound"


def test_adjacent_across_a_non_reportable_bin(M: Any, tmp_path: Path) -> None:
    """Bins 5 and 7 below and bin 6 unreported: after removing bin 6 the reportable bins 5 and 7
    are consecutive -> ``not_supported``."""
    q, d2, fdr = _ex11_frame(below=(5, 7), unreported=(6,))
    out = _h4(M, tmp_path, q, d2, fdr)
    assert out["B"] == 9 and out["bins"][6]["reportable"] is False
    assert _below(out) == [5, 7]
    assert out["verdict"] == "not_supported"


def test_adjacent_across_two_non_reportable_bins(M: Any, tmp_path: Path) -> None:
    """Bins 5 and 8 below, bins 6 and 7 unreported: reportable bins 5 and 8 are consecutive."""
    q, d2, fdr = _ex11_frame(below=(5, 8), unreported=(6, 7))
    out = _h4(M, tmp_path, q, d2, fdr)
    assert out["B"] == 8
    assert out["verdict"] == "not_supported"


def _tied_bins() -> tuple[pd.DataFrame, set[str], list[int], dict[int, int]]:
    """60 distinct low scores, 40 tied at 0.9 and 100 distinct high scores: ``score_bins``
    merges edges at the ties and leaves label 4 empty (labels 0-3, 5-9; bin 3 holds the 40
    ties, every other bin 20). Every peptide is confirmed (k_b = n_b)."""
    import numpy as np

    from gippa.mapping.denovo_precision import score_bins

    scores = [
        *(0.80 + 0.05 * i / 61 for i in range(1, 61)),
        *([0.9] * 40),
        *(0.95 + 0.05 * i / 100 for i in range(1, 101)),
    ]
    peps = synthetic_peptides(len(scores))
    bins = score_bins(np.asarray(scores))
    labels = sorted({int(b) for b in bins})
    sizes = {lab: int((bins == lab).sum()) for lab in labels}
    assert labels == [0, 1, 2, 3, 5, 6, 7, 8, 9], labels  # the fixture really skips label 4
    assert min(sizes.values()) >= 20, sizes
    return queries_frame(peps, scores), set(peps), labels, sizes


@pytest.mark.parametrize(
    ("below", "want"), [((3, 5), "not_supported"), ((3, 6), "consistent_with_lower_bound")]
)
def test_adjacent_across_a_skipped_tie_label(
    M: Any, tmp_path: Path, below: tuple[int, int], want: str
) -> None:
    """Merged ties skip label 4: reportable bins 3 and 5 are consecutive (skipped labels do not
    break adjacency); bins 3 and 6 are not (bin 5 between them is reportable and consistent)."""
    q, d2, labels, _ = _tied_bins()
    fdr = [PLACEMENT_ROW] + [fdr_row(b, *(LOW if b in below else (0.95, 0.9, 1.0))) for b in labels]
    out = _h4(M, tmp_path, q, d2, fdr)
    assert [b["score_bin"] for b in out["bins"]] == labels and out["B"] == 9
    assert _below(out) == list(below)
    assert out["verdict"] == want


def test_b_equal_1_single_bin_below_is_consistent(M: Any, tmp_path: Path) -> None:
    """B = 1 (one tied bin): the bin is below, but ``not_supported`` needs two bins."""
    peps = synthetic_peptides(150)
    q = queries_frame(peps, [0.9] * 150)
    out = _h4(M, tmp_path, q, set(peps[:90]), [fdr_row(0, *LOW)])
    assert out["B"] == 1 and _below(out) == [0]
    assert out["verdict"] == "consistent_with_lower_bound"


def test_b_equal_1_among_unreported_bins_is_consistent(M: Any, tmp_path: Path) -> None:
    """B = 1 with nine unreported bins: the one reportable bin below reads consistent."""
    q, d2, fdr = _ex11_frame(below=(5,), unreported=(0, 1, 2, 3, 4, 6, 7, 8, 9))
    out = _h4(M, tmp_path, q, d2, fdr)
    assert out["B"] == 1 and _below(out) == [5]
    assert out["verdict"] == "consistent_with_lower_bound"


def test_descriptive_extras(M: Any, tmp_path: Path) -> None:
    """Definition 9: the n-weighted mean genome precision and the pooled share Σk/Σn."""
    q, d2, fdr = _ex11_frame()
    out = _h4(M, tmp_path, q, d2, fdr)
    want = (8 * 0.9 + 0.71 + 0.41) / 10  # equal n per bin
    assert val(out["mean_genome_precision"]) == pytest.approx(want)
    pooled = out["pooled_L"]
    assert (pooled["k"], pooled["n"]) == (900, 1500)
    assert (pooled["ci_low"], pooled["ci_high"]) == pytest.approx(wilson(900, 1500), abs=1e-9)


def test_example_12_small_bins_not_evaluable(M: Any, tmp_path: Path) -> None:
    peps, scores, by_bin = equal_bins(10, 19)  # n_b = 19 < 20 in every bin
    d2 = {p for b in by_bin for p in b[:10]}
    fdr = [fdr_row(b, 0.6, 0.4, 0.8) for b in range(10)]
    out = _h4(M, tmp_path, queries_frame(peps, scores), d2, fdr)
    assert out["verdict"] == "not_evaluable" and out["B"] == 0
    assert len(out["bins"]) == 10  # every bin still listed with its numbers
    for b in out["bins"]:
        assert (b["n_b"], b["k_b"]) == (19, 10) and b["reportable"] is False
        assert b["L_b"]["value"] == pytest.approx(10 / 19)
        assert b["relation_adj"] != "below_lower_bound"


def test_example_12_c_assumed_not_evaluable(M: Any, tmp_path: Path) -> None:
    q, d2, fdr = _ex11_frame(below=(5, 6))  # would be not_supported
    for r in fdr[1:]:
        r["flags"] = "c_assumed"
    out = _h4(M, tmp_path, q, d2, fdr)
    assert out["verdict"] == "not_evaluable"
    assert all("n_b" in b and "L_b" in b for b in out["bins"])


def test_unreported_bin_is_not_reportable(M: Any, tmp_path: Path) -> None:
    q, d2, fdr = _ex11_frame(unreported=(3,))
    out = _h4(M, tmp_path, q, d2, fdr)
    assert out["B"] == 9 and out["bins"][3]["reportable"] is False
    assert out["bins"][7]["U_b"] == pytest.approx(upper(0.71, 0.78, 9))


def test_d2_empty_warns_and_never_refutes(M: Any, tmp_path: Path) -> None:
    q, _, fdr = _ex11_frame(below=(5, 6))
    out = _h4(M, tmp_path, q, set(), fdr)
    assert "a2_empty" in codes(out["warnings"])
    for b in out["bins"]:
        assert b["k_b"] == 0 and b["L_b"]["value"] == 0.0 and b["L_b"]["ci_low"] == 0.0
        assert b["L_b"]["ci_high"] == pytest.approx(wilson(0, b["n_b"])[1], abs=1e-9)
    assert out["verdict"] == "consistent_with_lower_bound"


def test_i_l_folded_on_both_sides(M: Any, tmp_path: Path) -> None:
    """Cascadia writes L only; DIA-NN may write I: a confirmed peptide matches after folding."""
    peps = ["GADLEFGHK", *synthetic_peptides(20)]
    q = queries_frame(peps, [0.9] * len(peps))  # one tied bin
    out = _h4(M, tmp_path, q, {"GADIEFGHK"}, [fdr_row(0, 0.6, 0.4, 0.8)])
    assert len(out["bins"]) == 1 and out["bins"][0]["k_b"] == 1


def test_bin_labels_cross_check_is_inconsistent(M: Any, tmp_path: Path) -> None:
    """Definition 2: the bin labels must equal ``mapping_fdr.tsv``'s ``score_bin`` values."""
    q, d2, fdr = _ex11_frame()
    with pytest.raises(M.Inconsistent):
        _h4(M, tmp_path, q, d2, [*fdr, fdr_row(10, 0.5, 0.4, 0.6)])
    with pytest.raises(M.Inconsistent):
        _h4(M, tmp_path, q, d2, fdr[:-1])  # bin 9 missing


def test_denovo_count_cross_check(M: Any, tmp_path: Path) -> None:
    """S6 / lead decision: |S| against the placement log's "(M de novo)"."""
    q, d2, fdr = _ex11_frame()
    with pytest.raises(M.Inconsistent):
        _h4(M, tmp_path, q, d2, fdr, n_denovo_logged=len(q) - 1)
    out = _h4(M, tmp_path, q, d2, fdr, n_denovo_logged=None)
    assert "denovo_count_unchecked" in codes(out["warnings"])
    ok = _h4(M, tmp_path, q, d2, fdr)
    assert "denovo_count_unchecked" not in codes(ok["warnings"])


def test_nested_share(M: Any, tmp_path: Path) -> None:
    """S3: the share of a bin's peptides that are a substring of another query (I→L folded)."""
    peps = ["PEPTLDEKAGG", "PEPTLDEK", "GHWWWGGNK"]
    q = queries_frame(["PEPTLDEKAGG", "PEPTIDEK", "GHWWWGGNK"], [0.9, 0.9, 0.9])
    out = _h4(M, tmp_path, q, set(peps), [fdr_row(0, 0.6, 0.4, 0.8)])
    assert val(out["bins"][0]["nested_share"]) == pytest.approx(1 / 3)


def test_l_b_searchable(M: Any, tmp_path: Path) -> None:
    """Definition 9 / S4: the A2 confirmation rate among the searchable peptides; without the
    searchable set (no ``--strain-db``) the field is ``not_evaluable``."""
    peps = synthetic_peptides(30)
    q = queries_frame(peps, [0.9] * 30)
    d2, searchable = set(peps[:6]), set(peps[:10])
    out = _h4(M, tmp_path, q, d2, [fdr_row(0, 0.6, 0.4, 0.8)], searchable=searchable)
    s = out["bins"][0]["L_b_searchable"]
    assert (s["k"], s["n"]) == (6, 10) and s["value"] == pytest.approx(0.6)
    none = _h4(M, tmp_path, q, d2, [fdr_row(0, 0.6, 0.4, 0.8)])["bins"][0]["L_b_searchable"]
    assert none == "not_evaluable" or none.get("reason") == "not_evaluable"


def test_caveats_rendered_with_the_result(M: Any, tmp_path: Path) -> None:
    """Definition 10 and S3: (0.8, 1] only; c_value from A2 (shared dependence); ALT error rates
    unestimated; dependent rows; contaminants in c_l; A2's own FDR."""
    q, d2, fdr = _ex11_frame()
    cav = _h4(M, tmp_path, q, d2, fdr)["caveats"]
    assert isinstance(cav, list) and all(isinstance(c, str) and c for c in cav)
    text = " ".join(cav).lower()
    for needle in ("0.8", "c_value", "unestimated", "dependen", "contaminant", "fdr"):
        assert needle in text, needle


def test_caveat_one_bin_failure_reads_consistent(M: Any, tmp_path: Path) -> None:
    """Lead decision (H4 fallback): the report's H4 caveats state the cost, that a failure
    confined to one score bin reads ``consistent_with_lower_bound`` (two consecutive reportable
    bins are needed). One caveat string carries both the one-bin statement and the verdict."""
    q, d2, fdr = _ex11_frame()
    cav = _h4(M, tmp_path, q, d2, fdr)["caveats"]
    hits = [
        c
        for c in cav
        if "consistent_with_lower_bound" in c
        and re.search(r"\b(one|single)\b[^.]*\bbin\b", c.lower())
        and re.search(r"consecutive|adjacent", c.lower())
    ]
    assert hits, cav


def _rule(below_in_order: list[bool]) -> bool:
    """The oracle of definition 8: some two consecutive entries (reportable bins in score_bin
    order) are both below."""
    return any(a and b for a, b in itertools.pairwise(below_in_order))


@settings(
    max_examples=40,
    deadline=None,
    phases=NO_SHRINK,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    below=st.sets(st.integers(min_value=0, max_value=9)),
    unreported=st.sets(st.integers(min_value=0, max_value=9)),
)
def test_verdict_follows_the_two_consecutive_bins_rule(
    M: Any, tmp_path: Path, below: set[int], unreported: set[int]
) -> None:
    """Any pattern of below / unreported bins on the example 11 frame: the verdict is
    ``not_supported`` iff two reportable bins, consecutive once the unreported ones are removed,
    are both below (derived here, not read from the output's relations)."""
    q, d2, fdr = _ex11_frame(below=tuple(sorted(below)), unreported=tuple(sorted(unreported)))
    rep = [b for b in range(10) if b not in unreported]
    out = _h4(M, tmp_path, q, d2, fdr)
    assert out["B"] == len(rep)
    if not rep:
        assert out["verdict"] == "not_evaluable"
        return
    lo_adj = wilson(90, 150, 1 - 0.05 / len(rep))[0]
    ci_high = {b: LOW[2] if b in below else {3: 0.49, 7: 0.78}.get(b, 0.95) for b in rep}
    flags = [ci_high[b] < lo_adj for b in rep]
    want = "not_supported" if _rule(flags) else "consistent_with_lower_bound"
    assert out["verdict"] == want, (sorted(below), sorted(unreported))


@settings(
    max_examples=40,
    deadline=None,
    phases=NO_SHRINK,
    suppress_health_check=[HealthCheck.function_scoped_fixture],
)
@given(
    per_bin=st.integers(min_value=15, max_value=40),
    ks=st.lists(st.floats(min_value=0, max_value=1), min_size=10, max_size=10),
    gs=st.lists(_genome, min_size=10, max_size=10),
    unreported=st.lists(st.booleans(), min_size=10, max_size=10),
    c_assumed=st.booleans(),
)
def test_verdict_never_supported_and_follows_rule(
    M: Any,
    tmp_path: Path,
    per_bin: int,
    ks: list[float],
    gs: list[list[float]],
    unreported: list[bool],
    c_assumed: bool,
) -> None:
    """The verdict re-derived from definitions 6-8 with an independent oracle: reportable =
    precision reported and n_b >= 20; below = genome ci_high < Wilson(k_b, n_b, 1 - 0.05/B)
    lower limit (textbook formula); ``not_supported`` iff two consecutive reportable bins are
    below; ``not_evaluable`` with no reportable bin or any ``c_assumed``."""
    peps, scores, by_bin = equal_bins(10, per_bin)
    k = [round(f * per_bin) for f in ks]
    d2 = {p for j, b in enumerate(by_bin) for p in b[: k[j]]}
    fdr = [
        fdr_row(
            j,
            None if unreported[j] else g[0],
            None if unreported[j] else max(0.0, g[0] - 0.05),
            None if unreported[j] else g[1],
            flags="unreported"
            if unreported[j]
            else ("c_assumed" if c_assumed and j == 0 else None),
        )
        for j, g in enumerate(gs)
    ]
    rep = [j for j in range(10) if not unreported[j] and per_bin >= 20]
    B = len(rep)
    lo_adj = {j: wilson(k[j], per_bin, 1 - 0.05 / B)[0] for j in rep}
    assume(all(abs(gs[j][1] - lo_adj[j]) > 1e-6 for j in rep))  # no float near-ties
    below = [gs[j][1] < lo_adj[j] for j in rep]
    out = _h4(M, tmp_path, queries_frame(peps, scores), d2, fdr)
    assert out["verdict"] in VERDICTS and out["verdict"] != "supported"
    assert out["B"] == B == sum(b["reportable"] for b in out["bins"])
    by_label = {b["score_bin"]: b for b in out["bins"]}
    for j, x in zip(rep, below, strict=True):
        assert by_label[j]["relation_adj"] == ("below_lower_bound" if x else "consistent")
    if (c_assumed and not unreported[0]) or not rep:
        assert out["verdict"] == "not_evaluable"
    elif _rule(below):
        assert out["verdict"] == "not_supported"
    else:
        assert out["verdict"] == "consistent_with_lower_bound"
    if B == 1:
        assert out["verdict"] != "not_supported"


# ------------------------------------------------------------------ hypothesis 4: definition 1


def _cfg() -> Any:
    from gippa.config import GippaConfig

    return GippaConfig(seed=7)


def test_denovo_queries_use_the_maximum_score(M: Any, tmp_path: Path) -> None:
    """S1: one peptide on three rows with scores 0.85, 0.99, 0.90 (first row lowest) is one
    query with score 0.99 (``aggregate_peptides`` before ``select_queries``)."""
    from gippa.io.peptides import read_peptide_table

    others = synthetic_peptides(12)
    rows = [
        dn_row("GHWNQSTAK", 0.85, scan=1),
        dn_row("GHWNQSTAK", 0.99, run="WT_2", scan=2),
        dn_row("GHWNQSTAK", 0.90, run="WT_3", charge=3, scan=3),
    ]
    rows += [dn_row(p, 0.81 + 0.01 * i, scan=10 + i) for i, p in enumerate(others)]
    rows += [dn_row("GHWNQK", 0.95, scan=50)]  # 6 aa: below the 7-aa minimum, not a query
    P3, _ = read_peptide_table(write_denovo_table(tmp_path / "a3.tsv.gz", rows))
    q = fn(M, "denovo_queries")(P3, _cfg())
    assert sorted(q["peptide"].astype(str)) == sorted(["GHWNQSTAK", *others])
    assert len(q) == 13
    score = dict(zip(q["peptide"].astype(str), q["score"].astype(float), strict=True))
    assert score["GHWNQSTAK"] == pytest.approx(0.99)


def test_multi_row_peptide_lands_in_the_bin_of_its_maximum(M: Any, tmp_path: Path) -> None:
    """S1 end to end: with the first-row score (0.80) the peptide would sit in bin 0; it must
    be counted in the top bin, whose score labels GIPPA's mapping step wrote."""
    from gippa.io.peptides import read_peptide_table

    others = synthetic_peptides(39)
    rows = [dn_row("GHWNQSTAK", 0.80, scan=1), dn_row("GHWNQSTAK", 0.999, run="WT_2", scan=2)]
    rows += [dn_row(p, 0.81 + 0.004 * i, scan=10 + i) for i, p in enumerate(others)]
    P3, _ = read_peptide_table(write_denovo_table(tmp_path / "a3.tsv.gz", rows))
    q = fn(M, "denovo_queries")(P3, _cfg())
    fdr = [fdr_row(b, 0.6, 0.4, 0.8) for b in range(10)]
    out = _h4(M, tmp_path, q, {"GHWNQSTAK"}, fdr)
    assert [b["k_b"] for b in out["bins"]] == [0] * 9 + [1]
