"""M7-14 rules 1-2: the pure TE2 comparison ``te2_result`` and the ``te2`` CLI of
validation/proxies/te_proxies.py (docs/specs/M7-14-te2-proxy-and-id-mixups.md; GIAP_SPEC.md §11.6
TE2 "no ``high`` ALT calls from A2 except at sites flagged ``start_codon`` or ``stop_codon``").

The spec's acceptance names tests/unit/test_te_proxies.py; test authors do not edit existing test
files, so the ``te2_result`` cases live here.

Interface pinned here (test-author reading)::

    validation.proxies.te_proxies.te2_result(sites: pd.DataFrame) -> ProxyResult
    validation.proxies.te_proxies.te2_proxy(genome: Path, bakta_json: Path) -> ProxyResult
    python validation/proxies/te_proxies.py te2 GENOME BAKTA_JSON [--out FILE]

* ``sites`` is the ``variants.tsv`` frame of ``run_variants`` (``flags`` a tuple of §9.4 names).
* ``ProxyResult.name == "TE2"``; ``counts`` holds exactly the rule 1 keys ``alt_sites`` (rows),
  ``high`` (``confidence = high``), ``high_start_codon`` / ``high_stop_codon`` (high rows with
  that flag; a row with both counts in both) and ``violations`` (high rows with neither).
* ``passed`` is ``violations == 0`` (a zero-count assertion: no ALT site at all passes).
* ``exceptions``: the violations in ``sites`` row order, at most 50, each with ``locus_tag``
  (``gene``), ``aa_pos``, ``ref`` (``ref_aa``) and ``alt`` (``alt_aa``).
* Rule 2: with ``high_stop_codon == 0`` the result JSON says why (the Sec sites stay unmapped
  until M7-16): ``to_json()`` mentions ``M7-16``.
* The CLI writes ``to_json()`` (sorted keys) to ``--out`` (default
  ``validation/proxies/<YYYY-MM-DD>_te2.json``) and exits 0 iff ``passed``, as ``te1``/``te8``.
"""

from __future__ import annotations

import importlib
import json
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
te: Any = importlib.import_module("validation.proxies.te_proxies")

COUNT_KEYS = {"alt_sites", "high", "high_start_codon", "high_stop_codon", "violations"}
FLAGS = ("mass_equivalent", "snv", "multi_nt", "strain_difference", "database_derived",
         "start_codon", "stop_codon")  # fmt: skip


def _row(i: int, conf: str, flags: tuple[str, ...], *, ref: str = "N", alt: str = "D",
         gene: str | None = None, aa_pos: int | None = None) -> dict[str, Any]:  # fmt: skip
    from giap.variants.vcf_records import VARIANTS_COLUMNS

    row: dict[str, Any] = dict.fromkeys(VARIANTS_COLUMNS, pd.NA)
    row.update(variant_id=f"VAR_{i}", contig="chr1", codon_start=3 * i, strand="+",
               gene=gene or f"NISSLE_{i:05d}", aa_pos=aa_pos if aa_pos is not None else 10 + i,
               ref_codon="AAC", ref_aa=ref, alt_aa=alt, alt_type="substitution", call="ALT",
               confidence=conf, n_ref=0, n_alt=2 if conf == "high" else 1, source="database",
               flags=tuple(flags))  # fmt: skip
    return row


def _sites(rows: list[dict[str, Any]]) -> pd.DataFrame:
    from giap.variants.vcf_records import VARIANTS_COLUMNS

    return pd.DataFrame(rows, columns=list(VARIANTS_COLUMNS))


def _check_shape(r: Any) -> None:
    c = r.counts
    assert r.name == "TE2"
    assert set(c) == COUNT_KEYS
    assert all(type(v) is int for v in c.values()), c
    assert c["high"] <= c["alt_sites"]
    assert c["violations"] <= c["high"]
    assert c["high_start_codon"] <= c["high"] and c["high_stop_codon"] <= c["high"]
    assert r.passed is (c["violations"] == 0)
    assert len(r.exceptions) == min(50, c["violations"])
    for ex in r.exceptions:
        assert {"locus_tag", "aa_pos", "ref", "alt"} <= set(ex), ex


# ------------------------------------------------------------------ examples


def test_only_start_codon_high_sites_pass() -> None:
    """Rule 2's expected exceptions: ALT M at GTG/TTG starts (``start_codon``) are allowed."""
    rows = [_row(1, "high", ("snv", "start_codon"), ref="V", alt="M", aa_pos=1),
            _row(2, "high", ("snv", "start_codon"), ref="L", alt="M", aa_pos=1),
            _row(3, "candidate", ("snv",))]  # fmt: skip
    r = te.te2_result(_sites(rows))
    _check_shape(r)
    assert r.counts == {"alt_sites": 3, "high": 2, "high_start_codon": 2, "high_stop_codon": 0,
                        "violations": 0}  # fmt: skip
    assert r.passed is True
    assert r.exceptions == []


def test_a_high_site_without_start_or_stop_flag_fails_and_names_the_gene() -> None:
    """Degenerate row "TE2 violation | the proxy fails; ... the example rows name the genes"."""
    rows = [
        _row(1, "high", ("snv", "start_codon"), ref="V", alt="M", aa_pos=1),
        _row(2, "high", ("snv",), ref="N", alt="D", gene="NISSLE_00042", aa_pos=117),
    ]
    r = te.te2_result(_sites(rows))
    _check_shape(r)
    assert r.counts["violations"] == 1
    assert r.passed is False
    ex = r.exceptions[0]
    assert (ex["locus_tag"], int(ex["aa_pos"]), ex["ref"], ex["alt"]) == (
        "NISSLE_00042",
        117,
        "N",
        "D",
    )


def test_candidate_sites_are_not_violations() -> None:
    """TE2 is the ``high``-only test (lead decision 1); candidates without a flag pass."""
    rows = [_row(i, "candidate", ("snv",)) for i in range(1, 6)]
    r = te.te2_result(_sites(rows))
    _check_shape(r)
    assert r.counts["alt_sites"] == 5 and r.counts["high"] == 0
    assert r.counts["violations"] == 0 and r.passed is True


def test_stop_codon_high_site_is_allowed_and_counted() -> None:
    rows = [_row(1, "high", ("multi_nt", "stop_codon"), ref="*", alt="U")]
    r = te.te2_result(_sites(rows))
    _check_shape(r)
    assert r.counts["high_stop_codon"] == 1 and r.counts["violations"] == 0
    assert r.passed is True


def test_a_site_with_both_flags_counts_in_both() -> None:
    rows = [_row(1, "high", ("start_codon", "stop_codon"), ref="*", alt="M")]
    r = te.te2_result(_sites(rows))
    _check_shape(r)
    assert r.counts["high_start_codon"] == 1 and r.counts["high_stop_codon"] == 1
    assert r.counts["violations"] == 0


@pytest.mark.parametrize("flag", ["strain_difference", "database_derived", "mass_equivalent"])
def test_other_flags_do_not_exempt_a_high_site(flag: str) -> None:
    r = te.te2_result(_sites([_row(1, "high", ("snv", flag))]))
    _check_shape(r)
    assert r.counts["violations"] == 1 and r.passed is False


def test_no_alt_site_passes_with_zero_counts() -> None:
    r = te.te2_result(_sites([]))
    _check_shape(r)
    assert r.counts == dict.fromkeys(COUNT_KEYS, 0)
    assert r.passed is True


def test_violation_examples_capped_at_50_in_row_order() -> None:
    rows = [_row(i, "high", ("snv",)) for i in range(1, 81)]
    r = te.te2_result(_sites(rows))
    _check_shape(r)
    assert r.counts["violations"] == 80
    assert len(r.exceptions) == 50
    assert [e["locus_tag"] for e in r.exceptions] == [f"NISSLE_{i:05d}" for i in range(1, 51)]


def test_rule_2_zero_stop_codon_sites_are_explained() -> None:
    """Spec: "no ``stop_codon`` site is expected yet; the result records 0 and says why"."""
    rows = [_row(1, "high", ("snv", "start_codon"), ref="V", alt="M", aa_pos=1)]
    r = te.te2_result(_sites(rows))
    assert r.counts["high_stop_codon"] == 0
    text = json.dumps(r.to_json())
    assert "M7-16" in text, text


def test_to_json_round_trip() -> None:
    rows = [_row(1, "high", ("snv",), gene="NISSLE_00007", aa_pos=3)]
    r = te.te2_result(_sites(rows))
    doc = json.loads(json.dumps(r.to_json(), sort_keys=True))
    assert doc["proxy"] == "TE2"
    assert doc["passed"] is False
    assert doc["counts"] == r.counts
    assert doc["exceptions"][0]["locus_tag"] == "NISSLE_00007"


# ------------------------------------------------------------------ properties


@st.composite
def _table(draw: st.DrawFn) -> pd.DataFrame:
    n = draw(st.integers(0, 70))
    rows = []
    for i in range(1, n + 1):
        conf = draw(st.sampled_from(["high", "candidate"]))
        flags = tuple(draw(st.lists(st.sampled_from(FLAGS), unique=True, max_size=3)))
        rows.append(_row(i, conf, flags))
    return _sites(rows)


@settings(max_examples=100, deadline=None)
@given(_table())
def test_property_counts_follow_rule_1(sites: pd.DataFrame) -> None:
    r = te.te2_result(sites)
    _check_shape(r)
    high = sites[sites["confidence"] == "high"]
    fl = [set(f) for f in high["flags"]]
    assert r.counts["alt_sites"] == len(sites)
    assert r.counts["high"] == len(high)
    assert r.counts["high_start_codon"] == sum("start_codon" in f for f in fl)
    assert r.counts["high_stop_codon"] == sum("stop_codon" in f for f in fl)
    bad = [g for g, f in zip(high["gene"], fl, strict=True)
           if not {"start_codon", "stop_codon"} & f]  # fmt: skip
    assert r.counts["violations"] == len(bad)
    assert [e["locus_tag"] for e in r.exceptions] == bad[:50]


@settings(max_examples=60, deadline=None)
@given(_table())
def test_property_flagging_every_site_start_codon_passes(sites: pd.DataFrame) -> None:
    flagged = sites.copy()
    flagged["flags"] = [(*f, "start_codon") for f in sites["flags"]]
    r = te.te2_result(flagged)
    assert r.passed is True and r.counts["violations"] == 0
    assert r.counts["high_start_codon"] == r.counts["high"]


# ------------------------------------------------------------------ CLI


def _fake(passed: bool) -> Any:
    return te.ProxyResult("TE2", {"alt_sites": 3, "high": 2, "high_start_codon": 2,
                                  "high_stop_codon": 0, "violations": 0 if passed else 1},
                          {}, passed, [])  # fmt: skip


@pytest.mark.parametrize("passed", [True, False])
def test_cli_te2_writes_the_result_json_and_exits_on_pass(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, passed: bool
) -> None:
    seen: list[tuple[Path, Path]] = []

    def fake(genome: Path, bakta_json: Path) -> Any:
        seen.append((genome, bakta_json))
        return _fake(passed)

    monkeypatch.setattr(te, "te2_proxy", fake)
    out = tmp_path / "te2.json"
    rc = te.main(["te2", str(tmp_path / "g.fna.gz"), str(tmp_path / "b.json.gz"), "--out",
                  str(out)])  # fmt: skip
    assert rc == (0 if passed else 1)
    assert seen == [(tmp_path / "g.fna.gz", tmp_path / "b.json.gz")]
    doc = json.loads(out.read_text(encoding="utf-8"))
    assert doc["proxy"] == "TE2" and doc["passed"] is passed
    assert doc["counts"]["violations"] == (0 if passed else 1)


def test_cli_te2_default_output_is_dated_under_validation_proxies(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Acceptance: "the result JSON under ``validation/proxies/`` (dates only there)"."""
    import datetime as dt

    monkeypatch.setattr(te, "te2_proxy", lambda g, b: _fake(True))
    monkeypatch.setattr(te, "HERE", tmp_path)
    assert te.main(["te2", "g.fna", "b.json"]) == 0
    assert (tmp_path / f"{dt.date.today().isoformat()}_te2.json").is_file()
