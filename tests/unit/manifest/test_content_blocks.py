"""MVP-03c content blocks of the pure builder: ``sensitivity``, ``context`` (``replicons``,
``oric``), ``operons``, ``mapping``, ``calibration`` and the ``tools``/``databases`` rows
(docs/specs/MVP-03-manifest-finalize.md D-C rows marked c, examples 7 and 8, edge cases "Ties",
"db_manifest.json in the provisional form", "No --db-dir", "Circular/linear replicons", acceptance
03c; lead decisions 9, 10; backlog M2-15, M4-19 d7, MVP-01 pass 2, M6-06/M6-08; ADR-0021
Amendment 1; ADR-0017 KEGG date read by M9-16).

Interface assumed (test-author choice; the spec names the fields and their sources only). New
keyword fields of ``giap.manifest.ManifestSources``, all already parsed (the builder reads no
file):

``sensitivity``         rows of ``tables/sensitivity.tsv`` as dicts (``family``, ``variant``,
                        ``n_primary``, ``n_variant``, ``n_shared``, ``jaccard`` (float or None
                        for NA), ``agreement``), ``[]`` for a header-only table, None when the
                        file is absent. Default None.
``island_diagnostics``  rows of ``context/island_diagnostics.parquet`` (``pyarrow``
                        ``to_pylist()``: ``contig``, ``scanned``, ``block_rule``, ``l_hat``,
                        ``block_len``, ``morans_i_lag1`` (NaN when not scanned), ``flags``
                        (list)), C1 order, or None when absent. Default None.
``calibration``         the mapping ``giap.data.calibration.load()`` returns, or None. Default
                        None: then no ``uncalibrated:`` flag is added (the run tail and the
                        registry step always pass ``load()``; the 03a builder tests stay exact).

Blocks pinned (every key always present in the manifest):

* ``sensitivity`` = ``{"level": cfg.output.sensitivity, "n_rows": len(rows) (0 when None),
  "families": {family: {variant: {n_primary, n_variant, n_shared, jaccard}}}}``.
* ``context`` = ``{"replicons": [...], "oric": [...]}``: ``replicons`` = the diagnostics rows
  in their order with the seven ``DIAGNOSTICS_COLUMNS`` keys (NaN → None); ``oric`` = the
  ``data`` of each unique ``context_stats.oriter.multiple_oric`` event, sorted by ``contig``.
  Both ``[]`` when there is no source.
* ``operons`` = ``{"pi0", "bin_counts"}`` from the ``operons.missed_operons`` event with every
  ±∞ (and NaN) → None; None when no such event.
* ``mapping`` = ``{"search_db", "n_queries"}`` from the ``mapping.search_db_resolved`` event;
  None when no such event.
* ``calibration`` = ``{key: {"calibrated": bool, "provenance": str | None}}`` for every
  top-level key of the mapping, and ``flags`` gains ``uncalibrated:<key>`` for each
  ``calibrated: false`` key (still sorted unique).
* ``databases`` = one ``{name, version, checksum, retrieved, size_bytes}`` per ``status: ok``
  item, sorted by name, ``retrieved`` copied verbatim (M9-16 prints it for ``kegg*``).
"""

from __future__ import annotations

import copy
import json
import math
import random
from typing import Any

import pytest
from helpers.manifest_data import build, event, layer_missing, schema, validate
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig

KEYS_03C = {"sensitivity", "context", "operons", "mapping", "calibration"}
ORIC_EVENT = "context_stats.oriter.multiple_oric"
MISSED_EVENT = "operons.missed_operons"
SEARCH_EVENT = "mapping.search_db_resolved"
DIAG_COLS = ("contig", "scanned", "block_rule", "l_hat", "block_len", "morans_i_lag1", "flags")


def _cfg(**over: Any) -> GiapConfig:
    return GiapConfig.model_validate({"seed": 7, **over})


def _bins(lower0: float = -math.inf, upper_last: float = math.inf) -> list[dict[str, Any]]:
    return [
        {"bin": 0, "lower_bp": lower0, "upper_bp": 0.0, "N_NO": 4, "N_cv": 3, "N_sw": 1},
        {"bin": 1, "lower_bp": 0.0, "upper_bp": 50.0, "N_NO": 7, "N_cv": 5, "N_sw": 2},
        {"bin": 2, "lower_bp": 50.0, "upper_bp": upper_last, "N_NO": 2, "N_cv": 0, "N_sw": 6},
    ]


def _missed(lower0: float = -math.inf, upper_last: float = math.inf) -> dict[str, Any]:
    return event("operons", MISSED_EVENT, {"pi0": 0.8, "bin_counts": _bins(lower0, upper_last)})


def _oric(contig: str, mids: list[float], accepted: bool = True) -> dict[str, Any]:
    oric = sum(mids) / len(mids) if accepted else None
    data = {"contig": contig, "midpoints": mids, "spread": max(mids) - min(mids),
            "accepted": accepted, "oric": oric}  # fmt: skip
    return event("context_stats", ORIC_EVENT, data)


def _search(value: str = "same", n: int = 12) -> dict[str, Any]:
    return event("placement", SEARCH_EVENT, {"key": "search_db", "search_db": value,
                                             "n_queries": n})  # fmt: skip


def _sens_row(family: str, variant: str, n_p: int, n_v: int, n_s: int,
              jac: float | None) -> dict[str, Any]:  # fmt: skip
    return {"family": family, "variant": variant, "n_primary": n_p, "n_variant": n_v,
            "n_shared": n_s, "jaccard": jac, "agreement": None}  # fmt: skip


def _diag_rows() -> list[dict[str, Any]]:
    return [
        {"contig": "chr1", "scanned": True, "block_rule": "auto", "l_hat": 3, "block_len": 7,
         "morans_i_lag1": 0.21, "flags": []},
        {"contig": "pA", "scanned": True, "block_rule": "auto", "l_hat": 1, "block_len": 3,
         "morans_i_lag1": -0.05, "flags": []},
        {"contig": "lin1", "scanned": False, "block_rule": "not_scanned", "l_hat": None,
         "block_len": None, "morans_i_lag1": math.nan, "flags": ["few_values"]},
    ]  # fmt: skip


CALIBRATION = {
    "operon_distance": {"d0": 50, "s": 20, "provenance": "GIAP_SPEC §7.1.3 default",
                        "calibrated": False},
    "terminatornet_yes": {"value": 0.5, "calibrated": False},
}  # fmt: skip


# ----------------------------------------------------------------------------- shape


def test_every_03c_block_is_present_and_validates() -> None:
    m = build()
    missing = KEYS_03C - set(m)
    assert not missing, f"manifest lacks {sorted(missing)}"
    validate(m)


@pytest.mark.parametrize("field", sorted(KEYS_03C))
def test_every_03c_block_is_declared_in_the_schema(field: str) -> None:
    prop = schema()["properties"].get(field)
    assert prop is not None, field
    assert prop != {"type": "object", "description": "MVP-03c."}, f"{field} still a stub"


def test_full_03c_manifest_validates_and_is_strict_json() -> None:
    m = build(
        events=[_missed(), _oric("chr1", [100.0, 300.0]), _search()],
        sensitivity=[_sens_row("regions", "missing_min", 3, 4, 2, 0.4)],
        island_diagnostics=_diag_rows(),
        calibration=copy.deepcopy(CALIBRATION),
    )
    validate(m)
    text = json.dumps(m, allow_nan=False, sort_keys=True)
    assert "Infinity" not in text and "NaN" not in text
    assert json.loads(text) == m


# ----------------------------------------------------------------------------- example 7


def test_example_7_minus_infinity_bin_edge_becomes_null() -> None:
    m = build(events=[_missed()])
    bins = m["operons"]["bin_counts"]
    assert bins[0]["lower_bp"] is None
    assert bins[-1]["upper_bp"] is None
    json.dumps(m, allow_nan=False)  # raises on a ±∞ that slipped through
    validate(m)


def test_operons_block_copies_pi0_and_finite_counts() -> None:
    m = build(events=[_missed()])
    want = _bins()
    want[0]["lower_bp"] = None
    want[-1]["upper_bp"] = None
    assert m["operons"] == {"pi0": 0.8, "bin_counts": want}


def test_operons_block_with_finite_edges_is_unchanged() -> None:
    m = build(events=[_missed(lower0=-100.0, upper_last=500.0)])
    assert m["operons"]["bin_counts"] == _bins(-100.0, 500.0)


def test_operons_nan_count_field_also_becomes_null() -> None:
    ev = _missed()
    ev["data"]["bin_counts"][1]["upper_bp"] = math.nan
    m = build(events=[ev])
    assert m["operons"]["bin_counts"][1]["upper_bp"] is None
    json.dumps(m, allow_nan=False)


def test_no_operons_event_gives_null_block() -> None:
    m = build()
    assert m["operons"] is None
    validate(m)


def test_operons_resume_duplicate_collapses() -> None:
    assert build(events=[_missed(), _missed()])["operons"] == build(events=[_missed()])["operons"]


def test_builder_does_not_mutate_the_event() -> None:
    ev = _missed()
    build(events=[ev])
    assert ev["data"]["bin_counts"][0]["lower_bp"] == -math.inf


# ----------------------------------------------------------------------------- example 8


@pytest.mark.parametrize("rows", [None, []], ids=["absent", "header_only"])
def test_example_8_sensitivity_off(rows: list[dict[str, Any]] | None) -> None:
    m = build(cfg=_cfg(output={"sensitivity": "off"}), sensitivity=rows)
    assert m["sensitivity"] == {"level": "off", "n_rows": 0, "families": {}}
    validate(m)


@pytest.mark.parametrize("level", ["basic", "full"])
def test_sensitivity_level_follows_config(level: str) -> None:
    m = build(cfg=_cfg(output={"sensitivity": level}), sensitivity=[])
    assert m["sensitivity"]["level"] == level


def test_sensitivity_default_level_is_basic() -> None:
    assert build()["sensitivity"] == {"level": "basic", "n_rows": 0, "families": {}}


def test_sensitivity_families_from_rows() -> None:
    rows = [
        _sens_row("regions", "missing_min", 3, 4, 2, 0.4),
        _sens_row("regions", "top3", 3, 3, 3, 1.0),
        _sens_row("islands", "missing_min", 0, 0, 0, None),
    ]
    m = build(sensitivity=rows)
    assert m["sensitivity"]["n_rows"] == 3
    assert m["sensitivity"]["families"] == {
        "regions": {
            "missing_min": {"n_primary": 3, "n_variant": 4, "n_shared": 2, "jaccard": 0.4},
            "top3": {"n_primary": 3, "n_variant": 3, "n_shared": 3, "jaccard": 1.0},
        },
        "islands": {
            "missing_min": {"n_primary": 0, "n_variant": 0, "n_shared": 0, "jaccard": None},
        },
    }
    validate(m)


def test_sensitivity_nan_jaccard_becomes_null() -> None:
    m = build(sensitivity=[_sens_row("regions", "top3", 0, 0, 0, math.nan)])
    assert m["sensitivity"]["families"]["regions"]["top3"]["jaccard"] is None
    json.dumps(m, allow_nan=False)


def test_sensitivity_off_flag_not_invented() -> None:
    m = build(cfg=_cfg(output={"sensitivity": "off"}))
    assert not [f for f in m["flags"] if "sensitivity" in f]


# ----------------------------------------------------------------------------- context


def test_context_replicons_equal_the_diagnostics_rows_in_order() -> None:
    rows = _diag_rows()
    m = build(island_diagnostics=copy.deepcopy(rows))
    got = m["context"]["replicons"]
    assert [r["contig"] for r in got] == ["chr1", "pA", "lin1"]  # C1 order, never sorted
    want = copy.deepcopy(rows)
    want[2]["morans_i_lag1"] = None  # NaN → null (strict JSON)
    assert got == want
    for r in got:
        assert tuple(sorted(r)) == tuple(sorted(DIAG_COLS))
    validate(m)


def test_context_not_scanned_l_hat_stays_null_never_zero() -> None:
    """ADR-0021 Amendment 1: "not estimated" is NA, never 0."""
    lin = build(island_diagnostics=_diag_rows())["context"]["replicons"][2]
    assert lin["l_hat"] is None and lin["block_len"] is None
    assert lin["block_rule"] == "not_scanned"
    assert lin["flags"] == ["few_values"]


def test_context_without_sources_is_empty() -> None:
    m = build()
    assert m["context"] == {"replicons": [], "oric": []}
    validate(m)


def test_context_oric_from_multiple_oric_events() -> None:
    a = _oric("chr1", [3884032.0, 3925819.0])
    b = _oric("chr2", [1000.0, 300000.0], accepted=False)
    m = build(events=[b, a, a])  # a resume duplicate and out of order
    assert m["context"]["oric"] == [a["data"], b["data"]]
    assert m["context"]["oric"][1]["oric"] is None
    validate(m)


def test_context_oric_events_are_not_warnings_or_flags() -> None:
    m = build(events=[_oric("chr1", [10.0, 20.0])])
    assert m["flags"] == [] and m["warnings"] == []


# ----------------------------------------------------------------------------- mapping


def test_mapping_block_from_the_search_db_event() -> None:
    m = build(events=[_search("other", 40)])
    assert m["mapping"] == {"search_db": "other", "n_queries": 40}
    validate(m)


def test_mapping_resume_duplicate_collapses() -> None:
    assert build(events=[_search(), _search()])["mapping"] == {"search_db": "same",
                                                               "n_queries": 12}  # fmt: skip


def test_mapping_stub_event_without_values_still_validates() -> None:
    """A ``GIAP_STEP`` stub writes ``mapping.search_db_resolved`` with ``data = {key, stub}``
    only (workflow/modules/local/giap_step/main.nf ``stubEvent``): no crash, a valid manifest,
    and no invented value."""
    stub = event("placement", SEARCH_EVENT, {"key": "search_db", "stub": True})
    m = build(events=[stub])
    validate(m)
    assert m["mapping"] in (None, {"search_db": None, "n_queries": None})


def test_no_mapping_event_gives_null_block() -> None:
    m = build()
    assert m["mapping"] is None
    validate(m)


# ----------------------------------------------------------------------------- calibration


def test_calibration_lists_both_bundled_keys_and_flags_them() -> None:
    """Acceptance 03c / backlog M2-15: every calibration.yaml key, ``uncalibrated:<key>`` only
    while ``calibrated: false`` (all three keys are calibrated since ADR-0022, M6-12 part C)."""
    from giap.data.calibration import load

    cal = load()
    m = build(calibration=cal)
    assert set(m["calibration"]) == {"operon_distance", "uniop_pair_threshold", "terminatornet_yes"}
    for key, entry in cal.items():
        assert m["calibration"][key] == {"calibrated": entry["calibrated"],
                                         "provenance": entry.get("provenance")}  # fmt: skip
        assert (f"uncalibrated:{key}" in m["flags"]) is (entry["calibrated"] is False)
    validate(m)


def test_calibration_block_exact() -> None:
    m = build(calibration=copy.deepcopy(CALIBRATION))
    assert m["calibration"] == {
        "operon_distance": {"calibrated": False, "provenance": "GIAP_SPEC §7.1.3 default"},
        "terminatornet_yes": {"calibrated": False, "provenance": None},
    }
    assert m["flags"] == ["uncalibrated:operon_distance", "uncalibrated:terminatornet_yes"]


def test_calibrated_key_gets_no_flag() -> None:
    cal = copy.deepcopy(CALIBRATION)
    cal["operon_distance"]["calibrated"] = True
    cal["operon_distance"]["provenance"] = "M6 calibration run"
    m = build(calibration=cal)
    assert m["calibration"]["operon_distance"] == {"calibrated": True,
                                                   "provenance": "M6 calibration run"}  # fmt: skip
    assert m["flags"] == ["uncalibrated:terminatornet_yes"]


def test_uncalibrated_flags_merge_sorted_with_event_flags() -> None:
    m = build(calibration=copy.deepcopy(CALIBRATION),
              events=[layer_missing("gene_table", "cog"),
                      event("x", "x.run_flag", {"flag": "c_assumed"})])  # fmt: skip
    assert m["flags"] == ["c_assumed", "layer_missing:cog", "uncalibrated:operon_distance",
                          "uncalibrated:terminatornet_yes"]  # fmt: skip
    assert m["layers_missing"] == ["cog"]


def test_uncalibrated_flag_is_accepted() -> None:
    """Needs Khoa item 1, answered by ADR-0041: ``uncalibrated:<key>`` and the three other
    run-scope tokens are registered run flags."""
    from giap.manifest import RUN_FLAGS, RUN_FLAGS_PENDING

    assert "uncalibrated:" in RUN_FLAGS
    assert "uncalibrated:" not in RUN_FLAGS_PENDING
    for tok in ("input_prefiltered", "cai_reference_missing", "no_detection_limit_model"):
        assert tok in RUN_FLAGS and tok not in RUN_FLAGS_PENDING


def test_uncalibrated_flags_match_the_schema_pattern() -> None:
    validate(build(calibration=copy.deepcopy(CALIBRATION)))


def test_calibration_is_not_in_the_provenance_digest() -> None:
    """D-D: the digest covers no flags; the calibration block follows from git_commit."""
    a = build(calibration=copy.deepcopy(CALIBRATION))
    b = build()
    assert a["provenance_digest"] == b["provenance_digest"]


def test_content_blocks_do_not_change_the_digest() -> None:
    a = build()
    b = build(events=[_missed(), _search(), _oric("chr1", [1.0, 2.0])],
              sensitivity=[_sens_row("regions", "top3", 1, 1, 1, 1.0)],
              island_diagnostics=_diag_rows())  # fmt: skip
    assert a["provenance_digest"] == b["provenance_digest"]


# ----------------------------------------------------------------------------- databases


def _item(name: str, version: str, retrieved: str, n: int, status: str = "ok") -> dict[str, Any]:
    return {"name": name, "version": version, "url": f"https://example.invalid/{name}",
            "checksum": "sha256:" + format(n, "x") * 64, "size_bytes": 100 * n,
            "retrieved": retrieved, "path": name, "status": status, "message": ""}  # fmt: skip


SETUP_ITEMS = [
    _item("kegg_pathway_names", "2026-09-30", "2026-10-02", 3),
    _item("kegg_ko_pathway", "2026-09-30", "2026-10-01", 2),
    _item("bakta_light", "6.1", "2026-09-01", 1),
    _item("cog2020_def", "2020", "2026-08-15", 4),
    _item("genomad", "1.9", "2026-09-03", 5, status="failed"),
]


def test_databases_pin_retrieved_per_item_sorted() -> None:
    m = build(db_manifest={"schema_version": "1.0", "items": copy.deepcopy(SETUP_ITEMS)})
    got = m["databases"]
    assert [d["name"] for d in got] == ["bakta_light", "cog2020_def", "kegg_ko_pathway",
                                        "kegg_pathway_names"]  # fmt: skip
    want = {i["name"]: i for i in SETUP_ITEMS if i["status"] == "ok"}
    for d in got:
        src = want[d["name"]]
        assert d == {k: src[k] for k in ("name", "version", "checksum", "retrieved",
                                          "size_bytes")}  # fmt: skip
    validate(m)


def test_kegg_retrieval_dates_reach_the_manifest() -> None:
    """M9-16 prints ``retrieved`` of every ``databases[]`` entry whose name starts with ``kegg``
    (ADR-0017 date)."""
    m = build(db_manifest={"schema_version": "1.0", "items": copy.deepcopy(SETUP_ITEMS)})
    kegg = {d["name"]: d["retrieved"] for d in m["databases"] if d["name"].startswith("kegg")}
    assert kegg == {"kegg_ko_pathway": "2026-10-01", "kegg_pathway_names": "2026-10-02"}


def test_failed_items_are_left_out() -> None:
    m = build(db_manifest={"schema_version": "1.0", "items": copy.deepcopy(SETUP_ITEMS)})
    assert "genomad" not in {d["name"] for d in m["databases"]}


def test_databases_are_in_the_provenance_core_and_digest() -> None:
    from giap.manifest import provenance_core

    doc = {"schema_version": "1.0", "items": copy.deepcopy(SETUP_ITEMS)}
    a = build(db_manifest=doc)
    assert provenance_core(a)["databases"] == a["databases"]
    later = copy.deepcopy(doc)
    later["items"][1]["retrieved"] = "2026-11-01"  # a KEGG refresh
    b = build(db_manifest=later)
    assert a["provenance_digest"] != b["provenance_digest"]


def test_databases_provisional_form_has_null_checksum_and_retrieved() -> None:
    m = build(db_manifest={"kegg_ko_pathway": {"version": "2026-09"},
                           "bakta_db": {"version": "6.0"}})  # fmt: skip
    assert [d["name"] for d in m["databases"]] == ["bakta_db", "kegg_ko_pathway"]
    for d in m["databases"]:
        assert d["checksum"] is None and d["retrieved"] is None
    assert m["databases"][1]["version"] == "2026-09"
    validate(m)


def test_empty_setup_manifest_gives_no_databases() -> None:
    assert build(db_manifest={"schema_version": "1.0", "items": []})["databases"] == []


# ----------------------------------------------------------------------------- tools


def _report(process: str, name: str, version: str, digest: str | None = None,
            lock: str | None = None) -> dict[str, Any]:  # fmt: skip
    return {"process": process, "name": name, "version": version,
            "container_digest": digest, "conda_lock_hash": lock}  # fmt: skip


def test_tools_rows_exact_and_sorted() -> None:
    reports = [
        _report("NOVEL:TBLASTN", "blast", "2.16.0", "sha256:" + "a" * 64),
        _report("ANNOTATE:BAKTA", "bakta", "1.11.4", "sha256:" + "b" * 64),
        _report("ANNOTATE:GIAP_STEP", "giap", "0.1.0"),
    ]
    m = build(executor="nextflow", trace=[], local_steps=(), tool_reports=reports)
    assert m["tools"] == [
        {"name": "bakta", "version": "1.11.4", "container_digest": "sha256:" + "b" * 64,
         "conda_lock_hash": None, "processes": ["ANNOTATE:BAKTA"]},
        {"name": "blast", "version": "2.16.0", "container_digest": "sha256:" + "a" * 64,
         "conda_lock_hash": None, "processes": ["NOVEL:TBLASTN"]},
    ]  # fmt: skip
    validate(m)


def test_tools_same_process_twice_listed_once() -> None:
    """Several tasks of one process (chunks) report the same tool line."""
    rep = _report("ANNOTATE:PROMOTECH", "promotech", "1.0", None, "e" * 64)
    m = build(executor="nextflow", trace=[], local_steps=(), tool_reports=[rep, dict(rep)])
    assert m["tools"] == [{"name": "promotech", "version": "1.0", "container_digest": None,
                           "conda_lock_hash": "e" * 64,
                           "processes": ["ANNOTATE:PROMOTECH"]}]  # fmt: skip


def test_tools_are_not_in_the_digest() -> None:
    a = build(executor="nextflow", trace=[], local_steps=(),
              tool_reports=[_report("A:BAKTA", "bakta", "1.11.4")])  # fmt: skip
    b = build(executor="nextflow", trace=[], local_steps=())
    assert a["provenance_digest"] == b["provenance_digest"]


# ----------------------------------------------------------------------------- properties


_EVENTS_03C = [
    _missed(),
    _oric("chr1", [100.0, 300.0]),
    _oric("chr2", [5.0, 9000.0], accepted=False),
    _search(),
    layer_missing("gene_table", "cog"),
    event("x", "x.run_flag", {"flag": "c_assumed"}, level="WARNING"),
    event("function", "function.run_flags", {"flags": ["null_degenerate", "p_detect_clipped"]}),
]


@settings(max_examples=60, deadline=None)
@given(st.permutations(list(range(len(_EVENTS_03C)))), st.lists(
    st.integers(0, len(_EVENTS_03C) - 1), max_size=4))  # fmt: skip
def test_builder_invariant_under_event_order_and_duplicates(
    order: list[int], dups: list[int]
) -> None:
    """Test idea: ``build_manifest`` is invariant under any permutation of the event lines (and
    a resume duplicate changes nothing)."""
    base = build(events=copy.deepcopy(_EVENTS_03C), calibration=copy.deepcopy(CALIBRATION))
    evs = [copy.deepcopy(_EVENTS_03C[i]) for i in order + dups]
    assert build(events=evs, calibration=copy.deepcopy(CALIBRATION)) == base


_EDGE = st.one_of(
    st.floats(allow_nan=True, allow_infinity=True),
    st.integers(-(10**6), 10**6).map(float),
)


@settings(max_examples=80, deadline=None)
@given(st.lists(st.tuples(_EDGE, _EDGE, st.integers(0, 50)), min_size=1, max_size=6),
       st.floats(0.0, 1.0))  # fmt: skip
def test_operons_block_is_always_strict_json(
    edges: list[tuple[float, float, int]], pi0: float
) -> None:
    """Whatever the bin edges (±∞, NaN), the manifest serialises with ``allow_nan=False`` and the
    finite values are kept exactly."""
    bins = [{"bin": i, "lower_bp": lo, "upper_bp": hi, "N_NO": n, "N_cv": n, "N_sw": n}
            for i, (lo, hi, n) in enumerate(edges)]  # fmt: skip
    m = build(events=[event("operons", MISSED_EVENT, {"pi0": pi0, "bin_counts": bins})])
    json.dumps(m, allow_nan=False)
    for got, src in zip(m["operons"]["bin_counts"], bins, strict=True):
        for key in ("lower_bp", "upper_bp"):
            v = src[key]
            assert got[key] == (v if math.isfinite(v) else None)


@settings(max_examples=40, deadline=None)
@given(st.lists(
    st.tuples(st.sampled_from(["regions", "islands", "determinants"]),
              st.sampled_from(["missing_min", "missing_exclude", "single_peptide", "top3"]),
              st.integers(0, 200), st.integers(0, 200),
              st.one_of(st.none(), st.floats(0.0, 1.0))),
    max_size=8, unique_by=lambda t: (t[0], t[1])))  # fmt: skip
def test_sensitivity_n_rows_and_families_partition_the_rows(
    rows: list[tuple[str, str, int, int, float | None]],
) -> None:
    table = [_sens_row(f, v, p, q, min(p, q), j) for f, v, p, q, j in rows]
    shuffled = list(table)
    random.Random(1).shuffle(shuffled)
    m = build(sensitivity=table)
    s = m["sensitivity"]
    assert s["n_rows"] == len(table)
    assert sum(len(v) for v in s["families"].values()) == len(table)
    for r in table:
        cell = s["families"][r["family"]][r["variant"]]
        assert cell["n_primary"] == r["n_primary"] and cell["jaccard"] == r["jaccard"]
    assert build(sensitivity=shuffled)["sensitivity"] == s
    validate(m)
