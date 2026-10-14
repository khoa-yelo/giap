"""Preflight rules part B: inputs, proteins, peptides, abundance, evidence tier (issue M1-10;
GIAP_SPEC.md §1.1, §7.0, §7.3.2, §2.1 decision 3; ADR-0005 items 2-3).

T-PRE-11 … T-PRE-22 each have a passing and a failing example. All inputs are in-memory literals
(records, DataFrames, Issues); the rules are pure functions.

Lead decisions used here: ADR-0005 item 4's decoy-share threshold was declined, so the decoy share
is reported as its own Info check (``peptides.decoy_share``, value = the share) and never warns;
"no row passes ``input.q_value_max``" warns (existing key). The abundance join runs only when
≥ 1 abundance ID matches a protein (else info "deferred", ADR-0005 item 3).
"""

from __future__ import annotations

import json
from collections import Counter
from collections.abc import Sequence

import numpy as np
import pandas as pd
import pytest

from giap.config import GiapConfig, resolve_config
from giap.io.abundance import AbundanceTable
from giap.io.issues import Issue
from giap.io.peptides import validate_peptides
from giap.io.protein_fasta import (
    COUNTED_AA,
    STANDARD_AA,
    LetterCounts,
    ProteinRecord,
    resolve_protein_id,
)
from giap.preflight.checks import Check
from giap.preflight.evidence import (
    abundance_checks,
    candidate_tier,
    inputs_checks,
    peptide_checks,
    protein_checks,
    tier_check,
)

LEVELS = ("ok", "info", "warn", "error")
PEPLEN = (7, 40)
DEFERRED = "deferred: no proteins.faa; IDs are matched to locus tags after annotation"

# ---------------------------------------------------------------- helpers


def cfg(*overrides: str) -> GiapConfig:
    return resolve_config(overrides=list(overrides))


def _sane(checks: Sequence[Check], input_name: str) -> dict[str, Check]:
    ids = [c.rule_id for c in checks]
    assert len(ids) == len(set(ids)), f"rule reported twice: {ids}"
    for c in checks:
        assert c.input == input_name, c
        assert c.level in LEVELS
        assert c.rule, f"{c.rule_id}: empty human rule"
        json.dumps(c.value)
    return {c.rule_id: c for c in checks}


def prot(header: str, seq: str) -> ProteinRecord:
    """An in-memory record shaped like ``read_protein_fasta`` output."""
    token = header.split(maxsplit=1)[0] if header.split() else ""
    group_id, members = resolve_protein_id(token)
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


def pchecks(
    records: Sequence[ProteinRecord],
    parse_error: str | None = None,
    peptide_len: tuple[int, int] = PEPLEN,
) -> dict[str, Check]:
    return _sane(protein_checks(records, parse_error, peptide_len), "proteins")


_FLOAT = ("q_value", "pep", "score", "intensity")


def raw_peptides(**cols: Sequence[object]) -> pd.DataFrame:
    """A raw frame typed as ``read_peptide_table`` would type it."""
    data: dict[str, pd.Series] = {}
    for name, values in cols.items():
        if name in _FLOAT:
            data[name] = pd.Series([np.nan if v is None else v for v in values], dtype="float64")
        elif name == "is_decoy":
            data[name] = pd.Series(list(values), dtype="boolean")
        elif name == "charge":
            data[name] = pd.Series(list(values), dtype="Int64")
        else:
            data[name] = pd.Series(list(values), dtype=pd.StringDtype())
    return pd.DataFrame(data)


def validated(config: GiapConfig, **cols: Sequence[object]) -> tuple[pd.DataFrame, list[Issue]]:
    return validate_peptides(raw_peptides(**cols), config.input)


def pepchecks(config: GiapConfig, **cols: Sequence[object]) -> dict[str, Check]:
    df, issues = validated(config, **cols)
    return _sane(peptide_checks(issues, df, config.input), "peptides")


def table(
    ids: Sequence[str],
    values: Sequence[float | None],
    bad: Sequence[bool] | None = None,
    ignored: tuple[str, ...] = (),
) -> AbundanceTable:
    return AbundanceTable(
        id_column="protein_id",
        value_column="ibaq",
        raw_ids=tuple(ids),
        ids=tuple(resolve_protein_id(i)[0] for i in ids),
        values=tuple(values),
        bad_rows=tuple(bad) if bad is not None else tuple(False for _ in ids),
        ignored_columns=ignored,
    )


def achecks(
    tbl: AbundanceTable | None,
    issues: Sequence[Issue] = (),
    proteins: Sequence[ProteinRecord] | None = None,
    *,
    scale: str = "linear",
    requested: str = "auto",
    pep_intensity: bool = False,
    config: GiapConfig | None = None,
) -> dict[str, Check]:
    config = config if config is not None else cfg()
    out = abundance_checks(
        tbl,
        list(issues),
        proteins,
        scale,
        requested,
        pep_intensity,
        config,  # type: ignore[arg-type]
    )
    return _sane(out, "abundance")


def prots(*ids: str) -> list[ProteinRecord]:
    return [prot(i, "M" + "K" * 50) for i in ids]


# ---------------------------------------------------------------- T-PRE-11 inputs.evidence_present


def test_T_PRE_11_pass() -> None:
    for has_pep, has_prot in ((True, False), (False, True), (True, True)):
        c = _sane(inputs_checks(has_pep, has_prot), "inputs")["inputs.evidence_present"]
        assert c.level == "ok"


def test_T_PRE_11_fail() -> None:
    c = _sane(inputs_checks(False, False), "inputs")["inputs.evidence_present"]
    assert c.level == "error"


# ---------------------------------------------------------------- T-PRE-12 proteins.ids


def test_T_PRE_12_pass() -> None:
    records = [prot("sp|P1|X", "M" + "A" * 50), prot("P2;P3", "M" + "C" * 50)]
    assert [r.id for r in records] == ["P1", "P2;P3"]
    checks = pchecks(records)
    assert set(checks) >= {
        "proteins.ids",
        "proteins.alphabet",
        "proteins.length_classes",
        "proteins.duplicate_sequences",
    }
    assert checks["proteins.ids"].level == "ok"


def test_T_PRE_12_fail_duplicate_id() -> None:
    records = [prot("P1", "M" + "A" * 50), prot("sp|P1|Y", "M" + "C" * 50)]
    c = pchecks(records)["proteins.ids"]
    assert c.level == "error"
    assert "P1" in json.dumps(c.value)


def test_T_PRE_12_fail_empty_header() -> None:
    # read_protein_fasta raises FastaFormatError("header has no protein ID") on ">".
    c = pchecks([], parse_error="line 1: header has no protein ID")["proteins.ids"]
    assert c.level == "error"


# ---------------------------------------------------------------- T-PRE-13 proteins.alphabet


def test_T_PRE_13_pass() -> None:
    assert pchecks([prot("P1", "MKUVL*")])["proteins.alphabet"].level == "ok"


def test_T_PRE_13_fail() -> None:
    c = pchecks([prot("P1", "MKX*VL")])["proteins.alphabet"]
    assert c.level == "warn"
    assert c.value == {"nonstandard": {"X": 1}, "internal_stops": 1}


# ---------------------------------------------------------------- T-PRE-14 proteins.length_classes


def test_T_PRE_14_pass() -> None:
    c = pchecks([prot("P1", "M" + "A" * 40), prot("P2", "M" + "C" * 19)])["proteins.length_classes"]
    assert c.level == "info"
    assert c.value == {"protein": 1, "peptide": 1, "dropped": 0}


def test_T_PRE_14_fail() -> None:
    c = pchecks([prot("P1", "M" + "A" * 40), prot("SHORT1", "MKVLAC")])["proteins.length_classes"]
    assert c.level == "warn"
    assert isinstance(c.value, dict)
    assert c.value["dropped"] == 1
    assert "SHORT1" in json.dumps(c.value)


def test_T_PRE_14_uses_peptide_len() -> None:
    c = pchecks([prot("P1", "MKVLAC")], peptide_len=(5, 40))["proteins.length_classes"]
    assert c.level == "info"
    assert isinstance(c.value, dict) and c.value["dropped"] == 0 and c.value["peptide"] == 1


# ---------------------------------------------------------------- T-PRE-15 duplicate_sequences


def test_T_PRE_15_pass() -> None:
    c = pchecks([prot("P1", "MKV" + "A" * 40), prot("P2", "MKV" + "C" * 40)])[
        "proteins.duplicate_sequences"
    ]
    assert c.level == "ok"


def test_T_PRE_15_fail() -> None:
    seq = "MKV" + "A" * 40
    c = pchecks([prot("P1", seq), prot("P2", seq), prot("P3", "MKV" + "C" * 40)])[
        "proteins.duplicate_sequences"
    ]
    assert c.level == "info"
    assert c.value == [["P1", "P2"]]


# ---------------------------------------------------------------- T-PRE-16 peptides.schema


def test_T_PRE_16_pass() -> None:
    c = pepchecks(cfg(), peptide=["PEPTIDEK"], q_value=[0.01])["peptides.schema"]
    assert c.level == "ok"


@pytest.mark.parametrize(
    ("cols", "code"),
    [
        ({"sequence": ["PEPTIDEK"]}, "missing_peptide_column"),
        ({"peptide": ["PEPTIDEK", "PEPTIDER"], "q_value": [1.2, 0.01]}, "out_of_range:q_value"),
        ({"peptide": ["PEPTIDEK", "PEPTIDER"], "source": ["dna", "database"]}, "bad_enum:source"),
    ],
    ids=["no_peptide_column", "q_value_1.2", "source_dna"],
)
def test_T_PRE_16_fail(cols: dict[str, list[object]], code: str) -> None:
    c = pepchecks(cfg(), **cols)["peptides.schema"]
    assert c.level == "error"
    assert code in json.dumps(c.value) or code in c.message


def test_T_PRE_16_reader_bad_type_issue_is_schema_error() -> None:
    df, issues = validated(cfg(), peptide=["PEPTIDEK"])
    bad = Issue("bad_type:q_value", "error", "q_value: value is not a valid float64 (1 rows)", 1)
    c = _sane(peptide_checks([bad, *issues], df, cfg().input), "peptides")["peptides.schema"]
    assert c.level == "error"


def test_peptide_checks_without_peptides_returns_nothing() -> None:
    assert peptide_checks([], None, cfg().input) == []


# ---------------------------------------------------------------- T-PRE-17 peptides.content


def test_T_PRE_17_pass() -> None:
    checks = pepchecks(
        cfg(), peptide=["PEPTIDEK", "ELVISLIVESK", "SAMPLERK"], q_value=[0.001, 0.005, 0.01]
    )
    assert checks["peptides.content"].level == "ok"
    share = checks["peptides.decoy_share"]
    assert share.level == "info"
    assert share.value == 0.0


def test_T_PRE_17_fail_dropped_rows() -> None:
    c = pepchecks(
        cfg(),
        peptide=["PEPTIDEK", "PEPK", "PEPT(UniMod:21)IDEK", "SAMPLERK"],
        q_value=[0.001, 0.001, 0.001, 0.02],
    )["peptides.content"]
    assert c.level == "warn"
    assert isinstance(c.value, dict)
    assert c.value["n_dropped"] == 2
    assert c.value["q_min"] == pytest.approx(0.001)
    assert c.value["q_max"] == pytest.approx(0.02)


def test_T_PRE_17_fail_no_row_passes_q_value_max() -> None:
    cols = {"peptide": ["PEPTIDEK", "SAMPLERK"], "q_value": [0.02, 0.03]}
    c = pepchecks(cfg(), **cols)["peptides.content"]
    assert c.level == "warn"
    assert isinstance(c.value, dict)
    assert c.value["n_dropped"] == 0
    assert c.value["q_min"] == pytest.approx(0.02)
    # input.q_value_max is read from config: raising it makes the same rows pass.
    assert pepchecks(cfg("input.q_value_max=0.05"), **cols)["peptides.content"].level == "ok"


def test_T_PRE_17_no_q_column_is_not_a_q_warning() -> None:
    c = pepchecks(cfg(), peptide=["PEPTIDEK", "SAMPLERK"])["peptides.content"]
    assert c.level == "ok"


@pytest.mark.parametrize(
    ("decoys", "share"),
    [([False, False, False, True], 0.25), ([True, True, False, False], 0.5), ([True] * 4, 1.0)],
)
def test_T_PRE_17_decoy_share_is_info_and_never_warns(decoys: list[bool], share: float) -> None:
    # Khoa declined ADR-0005 item 4: no decoy-share threshold, always reported as info.
    checks = pepchecks(
        cfg(),
        peptide=["PEPTIDEK", "ELVISLIVESK", "SAMPLERK", "QWERTYK"],
        q_value=[0.001] * 4,
        is_decoy=decoys,
        protein_ids=["P1"] * 4,
    )
    c = checks["peptides.decoy_share"]
    assert c.level == "info"
    assert c.value == pytest.approx(share)
    assert checks["peptides.content"].level == "ok"


def test_decoy_share_is_over_validated_rows() -> None:
    # The short decoy row is dropped by validation; share = 1 decoy / 2 kept rows.
    c = pepchecks(
        cfg(),
        peptide=["PEPTIDEK", "SAMPLERK", "PEPK"],
        is_decoy=[True, False, True],
        protein_ids=["P1", "P1", "P1"],
    )["peptides.decoy_share"]
    assert c.value == pytest.approx(0.5)


def test_no_decoy_share_key_in_config() -> None:
    assert not hasattr(cfg().preflight, "decoy_share_warn")


# ---------------------------------------------------------------- T-PRE-18 evidence tier


def _tier(**cols: Sequence[object]) -> str:
    df, _ = validated(cfg(), **cols)
    return candidate_tier(df, has_proteins=False)


def test_T_PRE_18_decoy_picked() -> None:
    tier = _tier(
        peptide=["PEPTIDEK", "SAMPLERK"],
        protein_ids=["P1", "rev_P1"],  # decoy by prefix
        q_value=[0.001, 0.2],
    )
    assert tier == "decoy_picked"
    c = tier_check("decoy_picked")
    assert (c.rule_id, c.input, c.level) == ("peptides.evidence_tier", "peptides", "info")
    assert isinstance(c.value, dict)
    assert c.value["tier"] == "decoy_picked"
    assert c.value["provisional"] is True


def test_T_PRE_18_decoy_without_protein_ids_is_not_picked() -> None:
    tier = _tier(peptide=["PEPTIDEK", "SAMPLERK"], is_decoy=[False, True], q_value=[0.001, 0.2])
    assert tier == "mayu"


@pytest.mark.parametrize("col", ["q_value", "pep"])
def test_T_PRE_18_mayu(col: str) -> None:
    assert _tier(peptide=["PEPTIDEK", "SAMPLERK"], **{col: [0.001, None]}) == "mayu"
    c = tier_check("mayu")
    assert c.level == "info" and isinstance(c.value, dict) and c.value["tier"] == "mayu"
    assert c.value.get("provisional", False) is False


def test_T_PRE_18_two_peptide() -> None:
    assert _tier(peptide=["PEPTIDEK", "SAMPLERK"]) == "two_peptide"
    # an all-NA q_value column is treated as absent
    assert _tier(peptide=["PEPTIDEK", "SAMPLERK"], q_value=[None, None]) == "two_peptide"
    c = tier_check("two_peptide")
    assert c.level == "info" and isinstance(c.value, dict) and c.value["tier"] == "two_peptide"


def test_T_PRE_18_protein_list() -> None:
    assert candidate_tier(None, has_proteins=True) == "protein_list"
    c = tier_check("protein_list")
    assert c.level == "info" and isinstance(c.value, dict) and c.value["tier"] == "protein_list"


# ---------------------------------------------------------------- T-PRE-19 abundance.shape


def test_T_PRE_19_pass() -> None:
    c = achecks(table(["P1", "P2"], [1e6, 2e5]))["abundance.shape"]
    assert c.level == "ok"


def test_T_PRE_19_fail_one_column() -> None:
    issue = Issue(
        "no_value_column", "error", "abundance table needs an ID column and a value column"
    )
    c = achecks(None, [issue])["abundance.shape"]
    assert c.level == "error"


def test_T_PRE_19_fail_zero_rows() -> None:
    issue = Issue("empty_table", "error", "abundance table has a header but no rows")
    assert achecks(None, [issue])["abundance.shape"].level == "error"


def test_T_PRE_19_warn_extra_columns() -> None:
    tbl = table(["P1", "P2"], [1e6, 2e5], ignored=("w",))
    issue = Issue("extra_columns", "warn", "only one value column is used; ignored columns: w")
    c = achecks(tbl, [issue])["abundance.shape"]
    assert c.level == "warn"
    assert "w" in json.dumps(c.value)


# ---------------------------------------------------------------- T-PRE-20 abundance.join


def test_T_PRE_20_pass() -> None:
    c = achecks(table(["P1", "P2"], [1e6, 2e5]), proteins=prots("P1", "P2"))["abundance.join"]
    assert c.level == "ok"


def test_T_PRE_20_fail() -> None:
    c = achecks(table(["P1", "P9"], [1e6, 2e5]), proteins=prots("P1", "P2"))["abundance.join"]
    assert c.level == "warn"
    assert c.value == {"abundance_not_in_proteins": ["P9"], "proteins_not_in_abundance": ["P2"]}


def test_join_deferred_without_proteins() -> None:
    c = achecks(table(["P1", "P2"], [1e6, 2e5]), proteins=None)["abundance.join"]
    assert c.level == "info"
    assert c.value == DEFERRED


def test_join_deferred_when_no_id_matches() -> None:
    # ADR-0005 item 3: no ID matches → table assumed keyed by locus tags; join deferred.
    c = achecks(table(["locus_1", "locus_2"], [1e6, 2e5]), proteins=prots("P1", "P2"))[
        "abundance.join"
    ]
    assert c.level == "info"
    assert isinstance(c.value, str) and c.value.startswith("deferred")


def test_join_group_ids_match_on_any_member() -> None:
    c = achecks(table(["P3", "P1;P7"], [1e6, 2e5]), proteins=prots("P2;P3", "P1"))["abundance.join"]
    assert c.level == "ok"


def test_join_runs_when_peptides_have_intensity() -> None:
    c = achecks(table(["P1", "P9"], [1e6, 2e5]), proteins=prots("P1", "P2"), pep_intensity=True)[
        "abundance.join"
    ]
    assert c.level == "warn"


def test_two_kinds_of_missing_abundance_kept_apart() -> None:
    # P2 identified but not quantified (NA row); P3 never in the table (not identified).
    tbl = table(["P1", "P2"], [1e6, None])
    checks = achecks(tbl, proteins=prots("P1", "P2", "P3"))
    join = checks["abundance.join"]
    assert isinstance(join.value, dict)
    assert join.value["proteins_not_in_abundance"] == ["P3"]
    values = checks["abundance.values"]
    assert isinstance(values.value, dict)
    assert values.value["n_na"] == 1


# ---------------------------------------------------------------- T-PRE-21 abundance.values


def test_T_PRE_21_pass() -> None:
    c = achecks(table(["P1", "P2"], [1e6, 2e5]))["abundance.values"]
    assert c.level == "ok"
    assert isinstance(c.value, dict)
    assert c.value["n_na"] == 0 and c.value["n_zero"] == 0


def test_T_PRE_21_fail_non_numeric_from_bad_rows() -> None:
    tbl = table(["P1", "P2", "P3"], [1e6, None, None], bad=[False, True, False])
    issue = Issue("non_numeric_values", "error", "non-numeric value in 'ibaq' (1 rows)", 1, ("2",))
    c = achecks(tbl, [issue])["abundance.values"]
    assert c.level == "error"
    assert isinstance(c.value, dict)
    assert c.value["n_bad"] == 1
    assert c.value["n_na"] == 1  # the bad cell is not counted as missing


def test_T_PRE_21_fail_negative_on_forced_linear() -> None:
    tbl = table(["P1", "P2"], [-1.0, 100.0])
    c = achecks(tbl, scale="linear", requested="linear")["abundance.values"]
    assert c.level == "error"
    # the same values auto-detect as log2, where negatives are allowed
    assert achecks(tbl, scale="log2", requested="auto")["abundance.values"].level == "ok"


def test_T_PRE_21_fail_too_many_missing() -> None:
    ids = [f"P{i}" for i in range(10)]
    tbl = table(ids, [None, None, None, *([1e6] * 7)])
    c = achecks(tbl)["abundance.values"]
    assert c.level == "warn"
    assert isinstance(c.value, dict) and c.value["n_na"] == 3
    loose = cfg("preflight.abundance_missing_warn_frac=0.5")
    assert achecks(tbl, config=loose)["abundance.values"].level == "ok"


def test_T_PRE_21_zeros_counted_not_missing() -> None:
    c = achecks(table(["P1", "P2"], [0.0, 1e6]))["abundance.values"]
    assert c.level == "ok"
    assert isinstance(c.value, dict) and c.value["n_zero"] == 1 and c.value["n_na"] == 0


# ---------------------------------------------------------------- T-PRE-22 abundance.scale


def test_T_PRE_22_pass_auto() -> None:
    c = achecks(table(["P1", "P2"], [21.3, 22.0]), scale="log2", requested="auto")[
        "abundance.scale"
    ]
    assert c.level == "info"
    assert isinstance(c.value, dict)
    assert c.value["scale"] == "log2" and c.value["source"] == "auto"
    assert c.value.get("used", True) is True


def test_T_PRE_22_cli_override() -> None:
    c = achecks(table(["P1", "P2"], [21.3, 22.0]), scale="linear", requested="linear")[
        "abundance.scale"
    ]
    assert c.level == "info"
    assert isinstance(c.value, dict)
    assert c.value["scale"] == "linear" and c.value["source"] == "cli"


def test_T_PRE_22_unused_when_peptides_have_intensity() -> None:
    # §2.1 decision 3: protein-level values are used only without peptide intensities.
    c = achecks(table(["P1", "P2"], [1e6, 2e5]), pep_intensity=True)["abundance.scale"]
    assert c.level == "info"
    assert isinstance(c.value, dict) and c.value["used"] is False
