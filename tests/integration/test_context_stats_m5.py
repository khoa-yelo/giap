"""M5-09 examples 1 and 4-8 and the schema rules, on the extended ``context_stats`` step
(docs/specs/M5-09-context-tables-step.md, acceptance item
``tests/integration/test_context_stats_m5.py``).

GIAP_SPEC.md §6.2 (the step: gene table → ``regions_stats``, ``autocorrelation``, ``islands``,
``detection_summary``, ori-ter), §9.3 (the two column dictionaries and the genes.tsv context
group), §9.5 (Frictionless), §6.6 / R7 (determinism), R4 (every claim carries FWER/q or an
``unestimated`` label).

Interface assumed
-----------------
``giap.steps.context_stats`` gains the optional input ``features`` (C15, flat layout
``work_dir/features.parquet``) and the outputs ``autocorrelation`` → ``tables/
autocorrelation.tsv``, ``islands`` → ``tables/islands.tsv``, ``gene_context`` →
``work/context/gene_context.parquet`` (spec "Public interface"); ``run_context_stats`` keeps its
signature and finds the new input and outputs in the local layout. The per-replicon island
diagnostics of ADR-0021 Amendment 1 land in ``work/context/island_diagnostics.parquet``
(``tests/unit/context/test_island_diagnostics.py`` fixes that frame; the lead has not named the
artifact yet — see the test-author report).

The rows themselves come from the merged M5-05 / M5-06 producers and from M5-08 (PR #132, not
merged): the step only assembles, sorts and writes them, so every check here compares the written
table with ``context_tables_data``'s oracle, which calls those producers on the step's own random
streams (``giap.rng.for_step("context_stats.<part>", cfg.seed)``).

The fixture is the M4-18 integration fixture plus an oriC and a planted island; M4-20's
``build_m4_fixture`` (which the spec names) is not merged — see ``context_tables_data``.
"""

from __future__ import annotations

import json
import logging
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import context_tables_data as cd

TABLE_DIR = Path(__file__).resolve().parents[2] / "src" / "giap" / "schemas" / "tables"
#: an island is only reachable on a small fixture with a fixed block and a wide detrend span
#: (M5-06 decision 50: left censoring caps the high direction too); the ADR-0021 calibration does
#: not apply to a run with a fixed ``block_len``, which is fine for a plumbing fixture.
ISLAND_KW: dict[str, Any] = {"span_frac": 0.30, "block_len": 1, "bootstrap": 199}


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    """Four work directories, each already run through ``run_context_stats``."""
    base = tmp_path_factory.mktemp("context_stats_m5")
    out: dict[str, Any] = {}
    out["plain"] = _run(base / "plain", cd.config())
    out["again"] = _run(base / "again", cd.config())
    out["island"] = _run(base / "island", cd.config(**ISLAND_KW),
                         islands=(cd.island_spec(20, 6, 6.0),))  # fmt: skip
    out["no_features"] = _run(base / "no_features", cd.config(), features=False)
    # same gene table and the same shared ``stats`` knobs as "plain" (``stats.bootstrap`` is
    # also M4-16's): only the M5-exclusive knobs differ — example 8
    out["m5knobs"] = _run(
        base / "m5knobs", cd.config(moran_lags=(1, 3), span_frac=0.30, block_len=1)
    )
    return out


class Run:
    """One finished step run: the fixture, the config and the written outputs."""

    def __init__(self, fixture: Any, cfg: Any, records: list[logging.LogRecord]) -> None:
        self.fixture = fixture
        self.cfg = cfg
        self.records = records
        self.tables = fixture.work / "tables"
        self.autocorrelation = self.tables / "autocorrelation.tsv"
        self.islands = self.tables / "islands.tsv"
        self.gene_context = fixture.work / "context" / "gene_context.parquet"
        self.diagnostics = fixture.work / "context" / "island_diagnostics.parquet"


def _run(work: Path, cfg: Any, *, islands: Any = (), features: bool = True) -> Run:
    from giap.steps.context_stats import run_context_stats

    fixture = cd.build(work, islands=islands, features=features)
    handler = _Recorder()
    log = logging.getLogger("giap.steps.context_stats")
    log.addHandler(handler)
    try:
        run_context_stats(work, cfg, None)
    finally:
        log.removeHandler(handler)
    return Run(fixture, cfg, handler.records)


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _schema(name: str) -> dict[str, Any]:
    found = [p for p in (TABLE_DIR / f"{name}.json", TABLE_DIR / f"{name}.schema.json")
             if p.exists()]  # fmt: skip
    assert len(found) == 1, found
    data = json.loads(found[0].read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


# --------------------------------------------------------------------------- example 1


def test_example_1_both_tables_are_written(runs: dict[str, Run]) -> None:
    run = runs["island"]
    assert run.autocorrelation.exists() and run.islands.exists()
    assert run.gene_context.exists()


@pytest.mark.parametrize("name", ["plain", "island", "no_features"])
@pytest.mark.parametrize("table", ["autocorrelation", "islands"])
def test_example_1_headers_equal_the_module_constants(runs: dict[str, Run], name: str,
                                                      table: str) -> None:  # fmt: skip
    from giap.context.autocorrelation import AUTOCORR_COLUMNS
    from giap.context.islands import ISLAND_COLUMNS

    want = AUTOCORR_COLUMNS if table == "autocorrelation" else ISLAND_COLUMNS
    df = cd.read_table(getattr(runs[name], table))
    assert tuple(df.columns) == tuple(want)


@pytest.mark.parametrize("name", ["plain", "island", "no_features"])
@pytest.mark.parametrize("table", ["autocorrelation", "islands"])
def test_example_1_frictionless_validate(runs: dict[str, Run], name: str, table: str) -> None:
    frictionless = pytest.importorskip("frictionless")

    path = getattr(runs[name], table)
    schema = frictionless.Schema.from_descriptor(
        {k: v for k, v in _schema(table).items() if k != "coordinates"}
    )
    resource = frictionless.Resource(path=path.name, basepath=str(path.parent), schema=schema)
    report = resource.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def test_example_1_the_island_fixture_really_has_an_island(runs: dict[str, Run]) -> None:
    """Fixture guard: the planted island is found, so the island rows below are exercised."""
    assert len(cd.read_table(runs["island"].islands)) >= 1


# --------------------------------------------------------------------------- row assembly


def test_autocorrelation_rows_are_the_two_producers(runs: dict[str, Run]) -> None:
    """Spec: autocorrelation = M5-05 rows plus M5-08 rows, in the spec's sort order."""
    run = runs["plain"]
    want = cd.oracle_autocorrelation(run.fixture, run.cfg)
    got = cd.read_table(run.autocorrelation)
    keys = ("contig", "statistic", "null_model", "lag", "stat")
    assert [cd.row_key(r, keys) for _, r in got.iterrows()] == [
        cd.row_key(r, keys) for _, r in want.iterrows()
    ]


def test_islands_rows_are_the_m5_06_frame(runs: dict[str, Run]) -> None:
    run = runs["island"]
    want, _ = cd.oracle_islands(run.fixture, run.cfg)
    got = cd.read_table(run.islands)
    assert list(got["island_id"]) == list(want["island_id"])
    assert [int(v) for v in got["start"]] == [int(v) + 1 for v in want["start"]]
    assert [int(v) for v in got["end"]] == [int(v) for v in want["end"]]


def test_every_replicon_has_autocorrelation_rows(runs: dict[str, Run]) -> None:
    """Edge case: "the table is never missing a replicon" — even one with no eligible gene."""
    run = runs["plain"]
    got = cd.read_table(run.autocorrelation)
    assert set(got["contig"]) >= set(run.fixture.replicons["contig"].astype(str))


def test_each_replicon_has_the_full_permutation_grid(runs: dict[str, Run]) -> None:
    """Edge case: a replicon with no eligible gene still contributes its (NA) rows, so every
    replicon has both statistics under both nulls at every lag."""
    run = runs["plain"]
    lo, hi = (int(v) for v in run.cfg.context.moran_lags)
    df = cd.read_table(run.autocorrelation)
    permuted = df[df["statistic"] != "oriter_slope"]
    for contig in run.fixture.replicons["contig"].astype(str):
        rows = permuted[permuted["contig"] == contig]
        assert len(rows) == 2 * 2 * (hi - lo + 1), contig
        assert set(rows["null_model"]) == {"gene_label", "tu_block"}
        assert set(rows["statistic"]) == {"morans_i", "join_count_bb"}


def test_rows_are_sorted_as_the_spec_says(runs: dict[str, Run]) -> None:
    """Contig in C1 order, ``statistic`` and ``null_model`` in §9.3 enum order, ``lag`` (NA
    last), ``stat``."""
    run = runs["plain"]
    got = cd.read_table(run.autocorrelation)
    order = {c: i for i, c in enumerate(run.fixture.replicons["contig"].astype(str))}
    keys = [
        (
            order[str(r["contig"])],
            cd.STATISTIC_ORDER.index(str(r["statistic"])),
            cd.NULL_ORDER.index(str(r["null_model"])),
            np.inf if cd.is_na(r["lag"]) else float(r["lag"]),
            str(r["stat"]),
        )
        for _, r in got.iterrows()
    ]
    assert keys == sorted(keys)


def test_islands_are_sorted_by_contig_start_end_id(runs: dict[str, Run]) -> None:
    run = runs["island"]
    got = cd.read_table(run.islands)
    order = {c: i for i, c in enumerate(run.fixture.replicons["contig"].astype(str))}
    keys = [(order[str(r["contig"])], int(r["start"]), int(r["end"]),
             int(str(r["island_id"]).rsplit("_", 1)[-1])) for _, r in got.iterrows()]  # fmt: skip
    assert keys == sorted(keys)


# --------------------------------------------------------------------------- no row loses its label


@pytest.mark.parametrize("name", ["plain", "island"])
def test_no_autocorrelation_row_loses_its_flags(runs: dict[str, Run], name: str) -> None:
    """The step is plumbing: every row reaches the TSV with exactly the flags its producer gave
    it (the statistician's check for this issue)."""
    run = runs[name]
    want = cd.oracle_autocorrelation(run.fixture, run.cfg)
    got = cd.read_table(run.autocorrelation)
    assert len(got) == len(want)
    for (_, g), (_, w) in zip(got.iterrows(), want.iterrows(), strict=True):
        assert cd.flag_set(g["flags"]) == cd.flag_set(w["flags"]), dict(w)


def test_no_island_row_loses_its_flags(runs: dict[str, Run]) -> None:
    run = runs["island"]
    want, _ = cd.oracle_islands(run.fixture, run.cfg)
    got = cd.read_table(run.islands)
    assert len(got) == len(want)
    for (_, g), (_, w) in zip(got.iterrows(), want.iterrows(), strict=True):
        assert cd.flag_set(g["flags"]) == cd.flag_set(w["flags"]), dict(w)


@pytest.mark.parametrize("name", ["plain", "island"])
def test_no_row_loses_its_missingness(runs: dict[str, Run], name: str) -> None:
    """A statistic that does not exist stays NA on the way through the step; an estimable one
    keeps its value (``%.6g``)."""
    run = runs[name]
    want = cd.oracle_autocorrelation(run.fixture, run.cfg)
    got = cd.read_table(run.autocorrelation)
    for column in cd.AUTOCORR_FLOATS:
        for (_, g), (_, w) in zip(got.iterrows(), want.iterrows(), strict=True):
            assert cd.is_na(g[column]) == cd.is_na(w[column]), (column, dict(w))
            if not cd.is_na(w[column]):
                assert float(g[column]) == pytest.approx(float(w[column]), rel=1e-5, abs=1e-12)


@pytest.mark.parametrize("name", ["plain", "island", "no_features"])
def test_no_row_carries_the_dropped_uncalibrated_flag(runs: dict[str, Run], name: str) -> None:
    """ADR-0021 Amendment 1 / M5-06 decision 48."""
    for table in ("autocorrelation", "islands"):
        df = cd.read_table(getattr(runs[name], table))
        for cell in df["flags"]:
            assert "uncalibrated" not in cd.flag_set(cell)


@pytest.mark.parametrize("name", ["plain", "island", "no_features"])
def test_every_flag_is_in_the_vocabulary(runs: dict[str, Run], name: str) -> None:
    from giap.context.vocab import CONTEXT_FLAGS

    for table in ("autocorrelation", "islands"):
        df = cd.read_table(getattr(runs[name], table))
        for cell in df["flags"]:
            assert cd.flag_set(cell) <= {str(f) for f in CONTEXT_FLAGS}


# --------------------------------------------------------------------------- R4 and the rules


@pytest.mark.parametrize("name", ["plain", "island", "no_features"])
def test_r4_every_autocorrelation_row_is_labelled(runs: dict[str, Run], name: str) -> None:
    """R4: a q-value, or NA with ``not_estimable``, or ``descriptive`` (the ori-ter rows)."""
    df = cd.read_table(runs[name].autocorrelation)
    assert len(df) > 0
    for _, row in df.iterrows():
        flags = cd.flag_set(row["flags"])
        labelled = not cd.is_na(row["q_value"]) or bool(flags & {"not_estimable", "descriptive"})
        assert labelled, dict(row)


@pytest.mark.parametrize("name", ["plain", "island", "no_features"])
def test_r4_na_p_value_implies_na_q_value(runs: dict[str, Run], name: str) -> None:
    df = cd.read_table(runs[name].autocorrelation)
    for _, row in df.iterrows():
        if cd.is_na(row["p_value"]):
            assert cd.is_na(row["q_value"]), dict(row)


@pytest.mark.parametrize("name", ["plain", "island", "no_features"])
def test_oriter_rows_are_descriptive_only(runs: dict[str, Run], name: str) -> None:
    """Schema rule: ``statistic = oriter_slope`` ⇔ ``null_model = bootstrap`` ⇔ ``p_value`` NA
    with the flag ``descriptive``; ``lag`` NA iff ``statistic = oriter_slope``."""
    df = cd.read_table(runs[name].autocorrelation)
    oriter = df["statistic"] == "oriter_slope"
    assert oriter.any()
    for _, row in df[oriter].iterrows():
        assert row["null_model"] == "bootstrap"
        assert cd.is_na(row["p_value"]) and cd.is_na(row["q_value"])
        assert "descriptive" in cd.flag_set(row["flags"])
        assert cd.is_na(row["lag"])
    for _, row in df[~oriter].iterrows():
        assert row["null_model"] != "bootstrap"
        assert not cd.is_na(row["lag"])


@pytest.mark.parametrize("name", ["plain", "island", "no_features"])
def test_probabilities_are_in_the_unit_interval(runs: dict[str, Run], name: str) -> None:
    df = cd.read_table(runs[name].autocorrelation)
    for column in ("p_value", "q_value"):
        values = [float(v) for v in df[column] if not cd.is_na(v)]
        assert all(0.0 <= v <= 1.0 for v in values), column


def test_island_schema_rules_on_the_data(runs: dict[str, Run]) -> None:
    """``q_value`` always NA, ``p_value`` ∈ (0, alpha], ``method = scan_max_stat``,
    ``n_tus`` ≥ ``min_tus``, ``island_id`` matches ``ISL_[1-9][0-9]*``."""
    from giap.context.islands import METHOD

    run = runs["island"]
    alpha = float(run.cfg.stats.alpha.context)
    min_tus = int(run.cfg.context.islands.min_tus)
    df = cd.read_table(run.islands)
    assert len(df) >= 1
    for _, row in df.iterrows():
        assert cd.is_na(row["q_value"])
        p = float(row["p_value"])
        assert 0.0 < p <= alpha
        assert row["method"] == METHOD
        assert int(row["n_tus"]) >= min_tus
        assert re.fullmatch(r"ISL_[1-9][0-9]*", str(row["island_id"]))
        assert row["direction"] in ("high", "low")


def test_island_ids_are_numbered_from_one(runs: dict[str, Run]) -> None:
    df = cd.read_table(runs["island"].islands)
    assert list(df["island_id"]) == [f"ISL_{i + 1}" for i in range(len(df))]


# --------------------------------------------------------------------------- example 4


def test_example_4_no_island_writes_the_header_line_only(runs: dict[str, Run]) -> None:
    from giap.context.islands import ISLAND_COLUMNS

    run = runs["plain"]
    want, _ = cd.oracle_islands(run.fixture, run.cfg)
    assert len(want) == 0, "fixture guard: the 'plain' run must find no island"
    assert run.islands.read_text(encoding="utf-8") == "\t".join(ISLAND_COLUMNS) + "\n"


def test_example_4_every_gene_context_island_id_is_na(runs: dict[str, Run]) -> None:
    ctx = pd.read_parquet(runs["plain"].gene_context)
    assert len(ctx) > 0
    assert ctx["island_id"].isna().all()


# --------------------------------------------------------------------------- example 5


def test_example_5_without_features_the_run_succeeds(runs: dict[str, Run]) -> None:
    assert runs["no_features"].autocorrelation.exists()
    assert not (runs["no_features"].fixture.work / "features.parquet").exists()


def test_example_5_no_oric_row_for_each_circular_chromosome(runs: dict[str, Run]) -> None:
    df = cd.read_table(runs["no_features"].autocorrelation)
    oriter = df[df["statistic"] == "oriter_slope"]
    reps = runs["no_features"].fixture.replicons
    circular = reps[(reps["topology"] == "circular") & (reps["kind"] == "chromosome")]
    assert set(oriter["contig"]) == set(circular["contig"].astype(str))
    for _, row in oriter.iterrows():
        assert "no_oric" in cd.flag_set(row["flags"])
        assert cd.is_na(row["value"])


def test_example_5_the_missing_layer_is_logged(runs: dict[str, Run]) -> None:
    """Spec: "manifest/log records ``layer_missing:features``" (the run-flag pattern)."""
    messages = [r.getMessage() for r in runs["no_features"].records]
    assert any("layer_missing:features" in m for m in messages), messages


def test_with_features_no_row_is_flagged_no_oric(runs: dict[str, Run]) -> None:
    df = cd.read_table(runs["plain"].autocorrelation)
    oriter = df[df["statistic"] == "oriter_slope"]
    assert len(oriter) == 3  # pooled, right, left on the one circular chromosome
    assert not any("no_oric" in cd.flag_set(c) for c in oriter["flags"])


def test_without_features_the_gene_context_replichores_are_na(runs: dict[str, Run]) -> None:
    ctx = pd.read_parquet(runs["no_features"].gene_context)
    assert ctx["replichore"].isna().all()
    assert ctx["replichore_position"].isna().all()


# --------------------------------------------------------------------------- example 6


def test_gene_context_columns_and_one_row_per_gene(runs: dict[str, Run]) -> None:
    run = runs["island"]
    ctx = pd.read_parquet(run.gene_context)
    assert tuple(ctx.columns) == cd.GENE_CONTEXT_COLUMNS
    assert list(ctx["locus_tag"]) == list(run.fixture.genes["locus_tag"].astype(str))


def test_example_6_gene_context_matches_the_producers(runs: dict[str, Run]) -> None:
    """Example 6: a gene of an island carries that island's ID, its replichore and its position
    — exactly the M5-06 gene map and the M5-08 replichore columns, left-joined on
    ``locus_tag``."""
    run = runs["island"]
    islands, _ = cd.oracle_islands(run.fixture, run.cfg)
    want = cd.oracle_gene_context(run.fixture, islands)
    got = pd.read_parquet(run.gene_context)
    assert list(got["island_id"].astype(object)) == list(want["island_id"].astype(object))
    assert list(got["replichore"].astype(object)) == list(want["replichore"].astype(object))
    np.testing.assert_allclose(
        got["replichore_position"].to_numpy(float),
        want["replichore_position"].to_numpy(float),
        rtol=1e-12,
        equal_nan=True,
    )


def test_example_6_a_gene_in_an_island_has_all_three_columns(runs: dict[str, Run]) -> None:
    run = runs["island"]
    table = cd.read_table(run.islands)
    ids = set(table["island_id"])
    ctx = pd.read_parquet(run.gene_context).set_index("locus_tag")
    seen = set()
    for _, row in table.iterrows():
        for tag in str(row["genes"]).split(";"):
            got = ctx.loc[tag]
            assert got["island_id"] in ids  # a gene can sit in a high and a low island
            assert got["replichore"] in ("right", "left")
            assert 0.0 <= float(got["replichore_position"]) <= 1.0
            seen.add(str(got["island_id"]))
    assert seen == ids


def test_replichore_position_is_in_the_unit_interval(runs: dict[str, Run]) -> None:
    ctx = pd.read_parquet(runs["plain"].gene_context)
    values = ctx["replichore_position"].dropna().to_numpy(float)
    assert len(values) > 0
    assert values.min() >= 0.0 and values.max() <= 1.0


def test_genes_tsv_schema_has_the_context_group() -> None:
    """§9.3 genes.tsv context group: ``island_id``, ``replichore``, ``replichore_position``
    (the writer left-joins ``gene_context`` on ``locus_tag``; without it they stay NA)."""
    from giap.genetable.assemble import GENES_COLUMNS

    fields = [f["name"] for f in _schema("genes")["fields"]]
    for column in ("island_id", "replichore", "replichore_position"):
        assert column in fields
        assert column in GENES_COLUMNS


# --------------------------------------------------------------------------- example 7


@pytest.mark.parametrize("table", ["autocorrelation", "islands"])
def test_example_7_same_seed_gives_byte_identical_tables(runs: dict[str, Run],
                                                         table: str) -> None:  # fmt: skip
    a, b = runs["plain"], runs["again"]
    assert getattr(a, table).read_bytes() == getattr(b, table).read_bytes()


def test_example_7_gene_context_is_byte_identical(runs: dict[str, Run]) -> None:
    a, b = runs["plain"], runs["again"]
    assert a.gene_context.read_bytes() == b.gene_context.read_bytes()


# --------------------------------------------------------------------------- example 8


@pytest.mark.parametrize("table", ["regions_stats.tsv", "detection_summary.tsv"])
def test_example_8_m4_tables_do_not_move_with_the_m5_configuration(runs: dict[str, Run],
                                                                   table: str) -> None:  # fmt: skip
    """Spec: "M5 uses its own RNG streams" — changing only M5 inputs and M5 knobs (the
    island block length and span, the oriC file) leaves the M4 tables byte-identical."""
    plain = (runs["plain"].tables / table).read_bytes()
    assert (runs["m5knobs"].tables / table).read_bytes() == plain
    assert (runs["no_features"].tables / table).read_bytes() == plain


@pytest.mark.parametrize("table", ["regions_stats.tsv", "detection_summary.tsv"])
def test_example_8_m4_tables_are_still_deterministic(runs: dict[str, Run], table: str) -> None:
    plain = (runs["plain"].tables / table).read_bytes()
    assert (runs["again"].tables / table).read_bytes() == plain


# --------------------------------------------------------------------------- diagnostics


def test_island_diagnostics_are_persisted(runs: dict[str, Run]) -> None:
    """ADR-0021 Amendment 1: l-hat, l' and Î₁ per replicon survive the step."""
    from giap.context.diagnostics import DIAGNOSTICS_COLUMNS

    run = runs["plain"]
    df = pd.read_parquet(run.diagnostics)
    assert tuple(df.columns) == DIAGNOSTICS_COLUMNS
    assert list(df["contig"]) == list(run.fixture.replicons["contig"].astype(str))


def test_island_diagnostics_match_the_scan(runs: dict[str, Run]) -> None:
    from giap.context.diagnostics import island_diagnostics

    run = runs["plain"]
    _, scan = cd.oracle_islands(run.fixture, run.cfg)
    want = island_diagnostics(scan, run.cfg)
    got = pd.read_parquet(run.diagnostics)
    pd.testing.assert_frame_equal(got.reset_index(drop=True), want.reset_index(drop=True),
                                  check_dtype=False)  # fmt: skip


def test_island_diagnostics_never_report_a_zero_block_length(runs: dict[str, Run]) -> None:
    for name in ("plain", "island"):
        df = pd.read_parquet(runs[name].diagnostics)
        for column in ("l_hat", "block_len"):
            values = [v for v in df[column] if not pd.isna(v)]
            assert all(int(v) >= 1 for v in values), (name, column, list(df[column]))


def test_a_fixed_block_length_reports_no_l_hat(runs: dict[str, Run]) -> None:
    """The island run fixes ``block_len = 1``: l-hat was never estimated, so it is NA, not 0."""
    df = pd.read_parquet(runs["island"].diagnostics)
    scanned = df[df["scanned"].astype(bool)]
    assert len(scanned) >= 1
    assert (scanned["block_rule"] == "fixed").all()
    assert scanned["l_hat"].isna().all()
    assert (scanned["block_len"].astype(float) == 1.0).all()


# --------------------------------------------------------------------------- the step contract


def test_step_declares_the_new_inputs_and_outputs() -> None:
    """Spec "Public interface": inputs ``tus`` (C5, required) and ``features`` (C15, optional);
    outputs ``autocorrelation``, ``islands``, ``gene_context``."""
    import layers_chain as lc

    with pytest.MonkeyPatch.context() as mp:
        registry = lc.fresh_steps(mp)
        spec = registry.get("context_stats")
    assert "tus" in spec.inputs
    assert "features" in spec.optional_inputs
    assert {"autocorrelation", "islands", "gene_context"} <= set(spec.outputs)
    assert {"regions_stats", "detection_summary"} <= set(spec.outputs)


def test_step_leaves_its_inputs_untouched(tmp_path: Path) -> None:
    from giap.steps.context_stats import run_context_stats

    run = _run(tmp_path / "w", cd.config())
    first = run.autocorrelation.read_bytes()
    run_context_stats(run.fixture.work, run.cfg, None)
    assert run.autocorrelation.read_bytes() == first


@pytest.mark.parametrize("name", ["plain", "island", "no_features"])
def test_genes_tsv_is_the_gene_table_joined_with_gene_context(
    runs: dict[str, Run], name: str
) -> None:
    """Lead decision 3 (code review of #138): ``tables/genes.tsv`` has one row per gene in
    gene-table order, a 1-based ``start``, and the three context columns of ``gene_context``."""
    from giap.io.genes_tsv import read_genes_tsv

    run = runs[name]
    path = run.tables / "genes.tsv"
    assert path.exists()
    got = read_genes_tsv(path)
    genes = run.fixture.genes
    assert list(got["locus_tag"]) == list(genes["locus_tag"].astype(str))
    assert list(got["start"]) == [int(s) for s in genes["start"]]  # reader undoes the +1
    raw_start = [line.split("\t") for line in path.read_text().splitlines()]
    col = raw_start[0].index("start")
    assert [int(r[col]) for r in raw_start[1:]] == [int(s) + 1 for s in genes["start"]]
    ctx = pd.read_parquet(run.gene_context)
    for c in ("island_id", "replichore"):
        assert [None if pd.isna(v) else v for v in got[c]] == [
            None if pd.isna(v) else v for v in ctx[c]
        ], c
    # written %.6g (M1-04 write_tsv), so equal to six significant digits
    np.testing.assert_allclose(
        got["replichore_position"].astype(float),
        ctx["replichore_position"].astype(float),
        rtol=1e-5,
        equal_nan=True,
    )
