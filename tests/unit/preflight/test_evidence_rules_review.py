"""Review follow-ups for the evidence preflight rules (issue M1-10; GIAP_SPEC.md §7.0, §7.3.2).

Pins the "Lead decisions after code review" in docs/specs/M1-10-preflight-evidence-rules.md and
kills mutations that survived the first test file:

- peptides given but every row dropped and no proteins file → inputs error "no usable evidence";
- ``decoy_picked`` needs decoy rows and a ranking column (q_value, pep or score); otherwise the
  tier falls back and ``peptides.content`` warns with the reason;
- "no row passes ``input.q_value_max``" counts target rows only;
- no ``peptides.decoy_share`` rule (the share is a value inside ``peptides.content``);
- ``peptides.schema`` holds errors only; ``unknown_column`` and warn findings go to content;
- protein records whose groups share a member are Info under ``proteins.ids`` (capped);
- ``abundance_checks`` resolves the scale once and reports that one value.

API assumptions made here (the decisions do not name them): ``inputs_checks`` takes a keyword
``n_usable_peptides`` (rows left after validation); the inputs check is found by its input, not
its rule id; ``abundance_checks`` may keep or drop its ``scale`` parameter (both are accepted).
"""

from __future__ import annotations

import inspect
import json
from collections import Counter
from collections.abc import Iterator, Sequence

import numpy as np
import pandas as pd
import pytest

from giap.config import GiapConfig, resolve_config
from giap.io.abundance import AbundanceTable, resolve_scale
from giap.io.issues import Issue
from giap.io.peptides import validate_peptides
from giap.io.protein_fasta import (
    COUNTED_AA,
    STANDARD_AA,
    LetterCounts,
    ProteinRecord,
    resolve_protein_id,
)
from giap.preflight.checks import LIST_CAP, Check
from giap.preflight.evidence import (
    abundance_checks,
    candidate_tier,
    inputs_checks,
    peptide_checks,
    protein_checks,
)

PEPLEN = (7, 40)

# ---------------------------------------------------------------- helpers


def cfg(*overrides: str) -> GiapConfig:
    return resolve_config(overrides=list(overrides))


def by_id(checks: Sequence[Check]) -> dict[str, Check]:
    ids = [c.rule_id for c in checks]
    assert len(ids) == len(set(ids)), f"rule reported twice: {ids}"
    for c in checks:
        json.dumps(c.value)
    return {c.rule_id: c for c in checks}


def prot(header: str, seq: str) -> ProteinRecord:
    """An in-memory record shaped like ``read_protein_fasta`` output."""
    group_id, members = resolve_protein_id(header.split(maxsplit=1)[0])
    seq = seq.upper()
    if seq.endswith("*"):
        seq = seq[:-1]
    counts = Counter(seq)
    stops = counts.pop("*", 0)
    return ProteinRecord(
        id=group_id,
        members=members,
        header=header,
        sequence=seq,
        internal_stops=stops,
        nonstandard=LetterCounts(sorted((c, n) for c, n in counts.items() if c in COUNTED_AA)),
        other=LetterCounts(
            sorted(
                (c, n) for c, n in counts.items() if c not in STANDARD_AA and c not in COUNTED_AA
            )
        ),
    )


def pchecks(records: Sequence[ProteinRecord]) -> dict[str, Check]:
    return by_id(protein_checks(records, None, PEPLEN))


_FLOAT = ("q_value", "pep", "score", "intensity")


def raw_peptides(**cols: Sequence[object]) -> pd.DataFrame:
    data: dict[str, pd.Series] = {}
    for name, values in cols.items():
        if name in _FLOAT:
            data[name] = pd.Series([np.nan if v is None else v for v in values], dtype="float64")
        elif name == "is_decoy":
            data[name] = pd.Series(list(values), dtype="boolean")
        else:
            data[name] = pd.Series(list(values), dtype=pd.StringDtype())
    return pd.DataFrame(data)


def validated(config: GiapConfig, **cols: Sequence[object]) -> tuple[pd.DataFrame, list[Issue]]:
    return validate_peptides(raw_peptides(**cols), config.input)


def pepchecks(config: GiapConfig, **cols: Sequence[object]) -> dict[str, Check]:
    df, issues = validated(config, **cols)
    return by_id(peptide_checks(issues, df, config.input))


def content_value(c: Check) -> dict[str, object]:
    assert isinstance(c.value, dict)
    return c.value


def the_inputs_check(checks: Sequence[Check]) -> Check:
    found = [c for c in checks if c.input == "inputs"]
    assert len(found) == 1, checks
    return found[0]


def abund(
    table: AbundanceTable,
    proteins: Sequence[ProteinRecord] | None = None,
    *,
    requested: str = "auto",
    stale_scale: str = "linear",
) -> dict[str, Check]:
    """Call ``abundance_checks``; ``stale_scale`` is passed only if the parameter still exists."""
    kwargs: dict[str, object] = {
        "table": table,
        "issues": [],
        "proteins": proteins,
        "requested_scale": requested,
        "peptides_have_intensity": False,
        "cfg": cfg(),
    }
    if "scale" in inspect.signature(abundance_checks).parameters:
        kwargs["scale"] = stale_scale
    return by_id(abundance_checks(**kwargs))  # type: ignore[arg-type]


def table(ids: Sequence[str], values: Sequence[float | None]) -> AbundanceTable:
    return AbundanceTable(
        id_column="protein_id",
        value_column="ibaq",
        raw_ids=tuple(ids),
        ids=tuple(ids),
        values=tuple(values),
        bad_rows=tuple(False for _ in ids),
        ignored_columns=(),
    )


def prots(*ids: str) -> list[ProteinRecord]:
    return [prot(i, "M" + "K" * 50) for i in ids]


def lists_in(v: object) -> Iterator[list[object]]:
    if isinstance(v, list):
        yield v
        for x in v:
            yield from lists_in(x)
    elif isinstance(v, dict):
        for x in v.values():
            yield from lists_in(x)


# ---------------------------------------------------------------- decision 1: no usable evidence


def test_all_peptide_rows_dropped_without_proteins_is_an_error() -> None:
    c = the_inputs_check(inputs_checks(True, False, n_usable_peptides=0))
    assert c.level == "error"
    assert "no usable evidence" in c.message


def test_all_peptide_rows_dropped_with_proteins_is_ok() -> None:
    assert the_inputs_check(inputs_checks(True, True, n_usable_peptides=0)).level == "ok"


def test_usable_peptide_rows_without_proteins_is_ok() -> None:
    assert the_inputs_check(inputs_checks(True, False, n_usable_peptides=3)).level == "ok"


def test_protein_list_tier_means_proteins_and_no_usable_peptides() -> None:
    df, _ = validated(cfg(), peptide=["PEPK", "SAMK"])  # both short → every row dropped
    assert df.empty
    assert candidate_tier(df, has_proteins=True) == "protein_list"


# ---------------------------------------------------------------- decision 2: ranking column


def test_decoys_without_ranking_column_fall_back_with_a_warning() -> None:
    cols = {"peptide": ["PEPTIDEK", "SAMPLERK"], "protein_ids": ["P1", "rev_P1"]}
    df, issues = validated(cfg(), **cols)
    assert bool(df["is_decoy"].any())
    assert candidate_tier(df, has_proteins=False) == "two_peptide"
    c = by_id(peptide_checks(issues, df, cfg().input))["peptides.content"]
    assert c.level == "warn"
    msg = c.message.lower()
    assert "decoy" in msg
    assert "rank" in msg or ("q_value" in msg and "score" in msg)


@pytest.mark.parametrize("col", ["q_value", "pep", "score"])
def test_decoys_with_a_ranking_column_are_decoy_picked(col: str) -> None:
    cols: dict[str, Sequence[object]] = {
        "peptide": ["PEPTIDEK", "SAMPLERK"],
        "protein_ids": ["P1", "rev_P1"],
        col: [0.001, 0.002],
    }
    df, issues = validated(cfg(), **cols)
    assert candidate_tier(df, has_proteins=False) == "decoy_picked"
    assert by_id(peptide_checks(issues, df, cfg().input))["peptides.content"].level == "ok"


# ---------------------------------------------------------------- decision 3: target rows only


def test_q_value_max_counts_target_rows_only() -> None:
    # only the decoy passes the cut-off → no target passes → warn
    c = pepchecks(
        cfg(),
        peptide=["PEPTIDEK", "SAMPLERK", "ELVISLIVESK"],
        protein_ids=["P1", "P2", "rev_P1"],
        q_value=[0.02, 0.03, 0.001],
    )["peptides.content"]
    assert c.level == "warn"
    assert "q_value" in c.message


def test_failing_decoys_do_not_trigger_the_q_warning() -> None:
    c = pepchecks(
        cfg(),
        peptide=["PEPTIDEK", "SAMPLERK"],
        protein_ids=["P1", "rev_P1"],
        q_value=[0.001, 0.5],
    )["peptides.content"]
    assert c.level == "ok"


def test_q_value_equal_to_max_passes() -> None:
    qmax = cfg().input.q_value_max
    c = pepchecks(cfg(), peptide=["PEPTIDEK", "SAMPLERK"], q_value=[qmax, 0.5])["peptides.content"]
    assert c.level == "ok"


def test_q_value_just_above_max_fails() -> None:
    qmax = cfg().input.q_value_max
    above = float(np.nextafter(qmax, 1.0))
    c = pepchecks(cfg(), peptide=["PEPTIDEK", "SAMPLERK"], q_value=[above, 0.5])["peptides.content"]
    assert c.level == "warn"


# ---------------------------------------------------------------- decision 4: no decoy_share rule


def test_no_separate_decoy_share_rule() -> None:
    checks = pepchecks(
        cfg(),
        peptide=["PEPTIDEK", "SAMPLERK", "ELVISLIVESK", "QWERTYK"],
        protein_ids=["P1", "P2", "rev_P1", "P3"],
        q_value=[0.001] * 4,
    )
    assert "peptides.decoy_share" not in checks
    assert set(checks) == {"peptides.schema", "peptides.content"}
    c = checks["peptides.content"]
    assert c.level == "ok"
    assert content_value(c)["decoy_share"] == pytest.approx(0.25)


# ---------------------------------------------------------------- decision 5: schema = errors only


def test_unknown_column_goes_to_content_not_schema() -> None:
    df, issues = validated(cfg(), peptide=["PEPTIDEK", "SAMPLERK"], q_value=[0.001, 0.002])
    unknown = Issue("unknown_column", "info", "unknown columns kept: ['foo']")
    checks = by_id(peptide_checks([unknown, *issues], df, cfg().input))
    schema = checks["peptides.schema"]
    assert schema.level == "ok"
    assert "unknown_column" not in json.dumps(schema.value)
    assert "unknown_column" not in schema.message
    assert "unknown_column" in json.dumps(checks["peptides.content"].value)


def test_warn_findings_go_to_content_not_schema() -> None:
    checks = pepchecks(
        cfg(), peptide=["PEPTIDEK", "PEPK", "SAMPLERK"], q_value=[0.001, 0.001, 0.001]
    )
    assert checks["peptides.schema"].level == "ok"
    assert "short_peptide" not in json.dumps(checks["peptides.schema"].value)
    assert checks["peptides.content"].level == "warn"
    assert "short_peptide" in json.dumps(checks["peptides.content"].value)


def test_schema_level_is_only_ok_or_error() -> None:
    df, issues = validated(cfg(), peptide=["PEPTIDEK", "PEPK"], is_decoy=[True, False])
    extra = [Issue("unknown_column", "info", "unknown columns kept: ['x']")]
    checks = by_id(peptide_checks([*extra, *issues], df, cfg().input))
    assert checks["peptides.schema"].level in ("ok", "error")


# ---------------------------------------------------------------- decision 6: shared members


@pytest.mark.parametrize(
    "headers", [("P1;P2", "P2;P3"), ("P1", "P1;P2")], ids=["two_groups", "single_and_group"]
)
def test_groups_sharing_a_member_are_info(headers: tuple[str, str]) -> None:
    records = [prot(headers[0], "M" + "A" * 50), prot(headers[1], "M" + "C" * 50)]
    c = pchecks(records)["proteins.ids"]
    assert c.level == "info"
    shared = "P2" if headers[0] == "P1;P2" else "P1"
    assert shared in json.dumps(c.value) or shared in c.message


def test_disjoint_groups_are_ok() -> None:
    records = [prot("P1;P2", "M" + "A" * 50), prot("P3;P4", "M" + "C" * 50)]
    assert pchecks(records)["proteins.ids"].level == "ok"


def test_duplicate_ids_stay_an_error_with_shared_members() -> None:
    records = [
        prot("P1", "M" + "A" * 50),
        prot("P1", "M" + "C" * 50),
        prot("P2;P3", "M" + "D" * 50),
        prot("P3;P4", "M" + "E" * 50),
    ]
    assert pchecks(records)["proteins.ids"].level == "error"


def test_shared_member_listing_is_capped() -> None:
    n = LIST_CAP + 5
    records = [prot(f"Q{i};S", "M" + "A" * 50) for i in range(n)]
    c = pchecks(records)["proteins.ids"]
    assert c.level == "info"
    lists = list(lists_in(c.value))
    assert lists, c.value
    assert all(len(x) <= LIST_CAP for x in lists)
    assert str(LIST_CAP) in c.message


# ---------------------------------------------------------------- decision 7: one resolved scale


@pytest.mark.parametrize(
    ("values", "requested", "stale"),
    [
        ([21.3, 22.0], "auto", "linear"),
        ([-1.0, 100.0], "auto", "linear"),
        ([1e6, 2e5], "auto", "log2"),
        ([21.3, 22.0], "linear", "log2"),
        ([1e6, 2e5], "log2", "linear"),
    ],
)
def test_scale_is_resolved_once_and_reported(
    values: list[float], requested: str, stale: str
) -> None:
    expected = resolve_scale(values, requested)  # type: ignore[arg-type]
    checks = abund(table(["P1", "P2"], values), requested=requested, stale_scale=stale)
    scale = content_value(checks["abundance.scale"])
    assert scale["scale"] == expected
    negative_error = "negative_on_linear" in json.dumps(checks["abundance.values"].value)
    assert negative_error == (expected == "linear" and min(values) < 0)


def test_undetermined_scale_reports_linear() -> None:
    checks = abund(table(["P1", "P2"], [None, None]), stale_scale="log2")
    assert content_value(checks["abundance.scale"])["scale"] == "linear"


# ---------------------------------------------------------------- surviving mutations


def test_empty_peptides_after_validation_still_report_content() -> None:
    # No schema error and every row dropped: content is reported over zero rows.
    checks = pepchecks(cfg(), peptide=["PEPK", "SAMK"], q_value=[0.001, 0.001])
    assert checks["peptides.schema"].level == "ok"
    c = checks["peptides.content"]
    v = content_value(c)
    assert v["n_rows"] == 0
    assert v["n_dropped"] == 2
    assert v["decoy_share"] is None
    assert c.level == "warn"


def test_schema_error_with_no_rows_short_circuits() -> None:
    df, issues = validated(cfg(), sequence=["PEPTIDEK"])
    assert df.empty
    out = peptide_checks(issues, df, cfg().input)
    assert [c.rule_id for c in out] == ["peptides.schema"]
    assert out[0].level == "error"


def test_schema_error_with_rows_left_still_reports_content() -> None:
    # The q_value=1.2 row is dropped with an error; the other row is still checked.
    df, issues = validated(cfg(), peptide=["PEPTIDEK", "SAMPLERK"], q_value=[1.2, 0.001])
    assert len(df) == 1
    checks = by_id(peptide_checks(issues, df, cfg().input))
    assert checks["peptides.schema"].level == "error"
    assert content_value(checks["peptides.content"])["n_rows"] == 1


def test_n_dropped_counts_only_dropping_issues() -> None:
    # decoy_without_protein_ids is a warn with n_rows but drops nothing.
    df, issues = validated(
        cfg(), peptide=["PEPTIDEK", "SAMPLERK", "PEPK"], is_decoy=[True, False, False]
    )
    assert any(i.code == "decoy_without_protein_ids" and i.n_rows for i in issues)
    extra = Issue("mixed_source", "warn", "peptide has both database and denovo rows", 7)
    c = by_id(peptide_checks([*issues, extra], df, cfg().input))["peptides.content"]
    assert content_value(c)["n_dropped"] == 1


@pytest.mark.parametrize(
    ("seq", "other"),
    [("MKOVL", {"O": 1}), ("MK1VL", {"1": 1}), ("MK-VL", {"-": 1}), ("MKO-9VL", None)],
)
def test_other_letters_warn_in_alphabet(seq: str, other: dict[str, int] | None) -> None:
    c = pchecks([prot("P1", seq)])["proteins.alphabet"]
    assert c.level == "warn"
    v = content_value(c)
    assert v["nonstandard"] == {}
    assert v["internal_stops"] == 0
    assert v["other"] == (other if other is not None else {"-": 1, "9": 1, "O": 1})


def test_join_with_empty_ids_and_no_match_is_deferred() -> None:
    c = abund(table(["", "locus_1"], [1e6, 2e5]), proteins=prots("P1", "P2"))["abundance.join"]
    assert c.level == "info"
    assert isinstance(c.value, str) and c.value.startswith("deferred")


def test_join_with_empty_ids_and_a_match_runs() -> None:
    c = abund(table(["", "P1"], [1e6, 2e5]), proteins=prots("P1", "P2"))["abundance.join"]
    assert c.level == "warn"
    assert c.value == {"abundance_not_in_proteins": [], "proteins_not_in_abundance": ["P2"]}
