"""Preflight framework and genome rules (issue M1-09; GIAP_SPEC.md §7.0, §2.1 decision 12, §5.2,
§7.1.2; ADR-0005 item 3).

T-PRE-01 … T-PRE-10 map 1:1 to the genome ``rule_id`` values, each with a passing and a failing
example. Rules are pure functions of ``GenomeSummary`` + config, so most tests feed in-memory
``GenomeRecord`` objects; a few go through ``read_genome_fasta`` on tiny files in ``tmp_path``.
"""

from __future__ import annotations

import builtins
import dataclasses
import json
import pathlib
from collections.abc import Sequence
from pathlib import Path

import pytest

import giapsim
from giap.config import GiapConfig, resolve_config
from giap.io.fasta import GenomeRecord, read_genome_fasta, sanitize_id
from giap.model import Replicon
from giap.preflight.checks import Check, worst
from giap.preflight.genome import (
    GenomeSummary,
    classify_replicons,
    genome_checks,
    n50,
    resolve_topology,
    summarize_genome,
)

RULE_IDS = (
    "genome.parse",
    "genome.ids_unique",
    "genome.ids_renamed",
    "genome.contig_count",
    "genome.total_length",
    "genome.contiguity",
    "genome.gaps",
    "genome.gc",
    "genome.topology",
    "genome.small_replicons",
)
LEVELS = ("ok", "info", "warn", "error")

# Small thresholds so that sequences for rules other than the size rules stay tiny.
SMALL = [
    "preflight.min_total_bp_warn=1",
    "preflight.min_n50_warn=1",
    "preflight.min_largest_contig_warn=1",
    "preflight.small_replicon_bp=1",
]


# ---------------------------------------------------------------- helpers


def rec(
    original_id: str, seq: str, *, circular: bool = False, header: str | None = None
) -> GenomeRecord:
    """An in-memory record shaped like ``read_genome_fasta`` output (sequence upper-cased)."""
    if header is None:
        header = original_id + (" circular=true" if circular else "")
    return GenomeRecord(
        id=sanitize_id(original_id),
        original_id=original_id,
        header=header,
        sequence=seq.upper(),
        circular_tag=circular,
    )


def acgt(n: int) -> str:
    """``n`` bp with GC exactly 0.5 when ``n`` is a multiple of 4."""
    return ("ACGT" * (n // 4 + 1))[:n]


def cfg(*overrides: str) -> GiapConfig:
    return resolve_config(overrides=list(overrides))


def small_cfg(*overrides: str) -> GiapConfig:
    return cfg(*SMALL, *overrides)


def summary(records: Sequence[GenomeRecord], n_run_min_bp: int = 10) -> GenomeSummary:
    return summarize_genome(records, n_run_min_bp)


def checks_for(
    records: Sequence[GenomeRecord],
    config: GiapConfig | None = None,
    parse_error: str | None = None,
) -> dict[str, Check]:
    config = config if config is not None else small_cfg()
    s = summary(records, config.preflight.n_run_min_bp)
    out = genome_checks(s, config, parse_error)
    ids = [c.rule_id for c in out]
    assert len(ids) == len(set(ids)), f"rule reported twice: {ids}"
    unknown = set(ids) - set(RULE_IDS)
    assert not unknown, f"unknown rule ids: {unknown}"
    for c in out:
        assert c.input == "genome"
        assert c.level in LEVELS
        assert c.rule, f"{c.rule_id}: empty human rule"
        json.dumps(c.value)  # value is JSON-serialisable
    return {c.rule_id: c for c in out}


def check(records: Sequence[GenomeRecord], rule_id: str, config: GiapConfig | None = None) -> Check:
    return checks_for(records, config)[rule_id]


def write_fa(tmp_path: Path, text: str, name: str = "g.fa") -> Path:
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def from_file(tmp_path: Path, text: str, config: GiapConfig | None = None) -> dict[str, Check]:
    return checks_for(read_genome_fasta(write_fa(tmp_path, text)), config)


# ---------------------------------------------------------------- framework


def test_check_dataclass_defaults_and_frozen() -> None:
    c = Check(
        rule_id="genome.gc",
        input="genome",
        level="ok",
        rule="GC outside 25-75 % -> warn",
        value=0.5,
    )
    assert c.message == ""
    assert c.exit_code == 0
    with pytest.raises(dataclasses.FrozenInstanceError):
        c.level = "warn"  # type: ignore[misc]


def _c(level: str) -> Check:
    return Check(rule_id="genome.gc", input="genome", level=level, rule="r", value=None)  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("levels", "expected"),
    [
        (["ok"], "ok"),
        (["ok", "info"], "info"),
        (["info", "warn", "ok"], "warn"),
        (["warn", "error", "info", "ok"], "error"),
        (["error", "ok"], "error"),
    ],
)
def test_worst_orders_ok_info_warn_error(levels: list[str], expected: str) -> None:
    assert worst(_c(lv) for lv in levels) == expected


def test_every_rule_evaluated_and_recorded_ok_on_clean_genome() -> None:
    """ADR-0005 item 3: every evaluated rule is recorded, with level ``ok`` when it passes."""
    got = checks_for([rec("chr", acgt(4000), circular=True), rec("p1", acgt(400), circular=True)])
    assert set(got) == set(RULE_IDS)
    for rid, c in got.items():
        assert c.level == ("info" if rid == "genome.topology" else "ok"), (rid, c)


def test_rules_are_pure_no_file_access(monkeypatch: pytest.MonkeyPatch) -> None:
    config = small_cfg()  # reading defaults.yaml happens before the guard
    records = [rec("chr:1", "ACG" + "N" * 12 + "TACRT", circular=True)]

    def _no_io(*_a: object, **_k: object) -> None:
        raise AssertionError("preflight genome rules must not touch files")

    monkeypatch.setattr(builtins, "open", _no_io)
    monkeypatch.setattr(pathlib.Path, "open", _no_io)
    monkeypatch.setattr(pathlib.Path, "read_text", _no_io)
    s = summarize_genome(records, 10)
    out = genome_checks(s, config)
    topo = resolve_topology(s, "auto")
    classify_replicons(s, topo, small_bp=500_000)
    assert {c.rule_id for c in out} == set(RULE_IDS)


# ---------------------------------------------------------------- summary


def test_summarize_genome_fields() -> None:
    s = summary([rec("chr:1", "ACGTRRNN", circular=True), rec("p", "GGCCAATT")])
    assert s.ids == ("chr_1", "p")
    assert s.original_ids == ("chr:1", "p")
    assert s.lengths == (8, 8)
    assert s.circular_tags == (True, False)
    assert s.acgt == 12
    assert s.gc == 6
    # merged runs, 0-based half-open internally (§6.5; M1-09 review decision)
    assert s.ambiguous == (("chr_1", 4, 6, "R"),)
    assert s.n_runs == ()
    assert dict(s.invalid) == {}
    assert s.empty_records == ()


def test_gc_excludes_ambiguous_codes() -> None:
    c = check([rec("c", "ACGT" + "R" * 20 + "N" * 4)], "genome.gc")
    assert c.level == "ok"
    assert c.value == pytest.approx(0.5)


# ---------------------------------------------------------------- T-PRE-01 genome.parse


def test_T_PRE_01_pass(tmp_path: Path) -> None:
    c = from_file(tmp_path, ">c\nACGTN\n")["genome.parse"]
    assert c.level == "ok"


def test_T_PRE_01_fail(tmp_path: Path) -> None:
    # empty file: zero records → error, and no other genome rule is evaluated
    got = from_file(tmp_path, "")
    assert got["genome.parse"].level == "error"
    assert set(got) == {"genome.parse"}

    # non-IUPAC letter
    got = from_file(tmp_path, ">c\nACGU\n")
    assert got["genome.parse"].level == "error"
    assert got["genome.parse"].value == {"U": 1}

    # record with an empty sequence
    got = from_file(tmp_path, ">c\nACGT\n>e\n>d\nACGT\n")
    assert got["genome.parse"].level == "error"


def test_T_PRE_01_parse_error_from_reader() -> None:
    got = checks_for(
        [], small_cfg(), parse_error="sequence data before the first '>' header (line 1)"
    )
    assert got["genome.parse"].level == "error"
    assert set(got) == {"genome.parse"}


def test_parse_error_still_evaluates_other_rules() -> None:
    """When genome.parse is an error but records exist, the other rules are still reported."""
    got = checks_for([rec("c", "ACGU" + acgt(400))])
    assert got["genome.parse"].level == "error"
    assert set(got) == set(RULE_IDS)


def test_summary_invalid_and_empty_records() -> None:
    s = summary([rec("c", "ACGUU"), rec("e", ""), rec("d", "AC*G")])
    assert dict(s.invalid) == {"U": 2, "*": 1}
    assert s.empty_records == ("e",)


# ---------------------------------------------------------------- T-PRE-02 genome.ids_unique


def test_T_PRE_02_pass(tmp_path: Path) -> None:
    c = from_file(tmp_path, ">a\nACGT\n>b\nACGT\n")["genome.ids_unique"]
    assert c.level == "ok"


def test_T_PRE_02_fail(tmp_path: Path) -> None:
    c = from_file(tmp_path, ">a x\nACGT\n>a y\nACGT\n")["genome.ids_unique"]
    assert c.level == "error"
    # distinct original IDs that collide after safe mapping
    c = from_file(tmp_path, ">a|b\nACGT\n>a:b\nACGT\n")["genome.ids_unique"]
    assert c.level == "error"
    assert "a_b" in json.dumps(c.value)


# ---------------------------------------------------------------- T-PRE-03 genome.ids_renamed


def test_T_PRE_03_pass(tmp_path: Path) -> None:
    c = from_file(tmp_path, ">NZ_CP1.1\nACGT\n")["genome.ids_renamed"]
    assert c.level == "ok"


def test_T_PRE_03_fail(tmp_path: Path) -> None:
    c = from_file(tmp_path, ">chr:1\nACGT\n")["genome.ids_renamed"]
    assert c.level == "warn"
    assert c.value == [{"from": "chr:1", "to": "chr_1"}]


# ---------------------------------------------------------------- T-PRE-04 genome.contig_count


def test_T_PRE_04_pass() -> None:
    c = check([rec(f"c{i}", acgt(40)) for i in range(10)], "genome.contig_count")
    assert c.level == "ok"
    assert c.value == 10


def test_T_PRE_04_fail() -> None:
    c = check([rec(f"c{i}", acgt(40)) for i in range(11)], "genome.contig_count")
    assert c.level == "warn"
    assert c.value == 11


def test_contig_count_threshold_from_config() -> None:
    records = [rec(f"c{i}", acgt(40)) for i in range(3)]
    assert (
        check(records, "genome.contig_count", small_cfg("preflight.max_contigs_warn=3")).level
        == "ok"
    )
    assert (
        check(records, "genome.contig_count", small_cfg("preflight.max_contigs_warn=2")).level
        == "warn"
    )


# ---------------------------------------------------------------- T-PRE-05 genome.total_length


def test_T_PRE_05_pass() -> None:
    c = check([rec("c", acgt(100_000))], "genome.total_length", cfg())
    assert c.level == "ok"
    assert c.value == 100_000


def test_T_PRE_05_fail() -> None:
    c = check([rec("c", acgt(99_999))], "genome.total_length", cfg())
    assert c.level == "warn"
    assert c.value == 99_999


def test_total_length_threshold_from_config() -> None:
    records = [rec("c", acgt(400)), rec("d", acgt(100))]
    assert (
        check(records, "genome.total_length", cfg("preflight.min_total_bp_warn=500")).level == "ok"
    )
    assert (
        check(records, "genome.total_length", cfg("preflight.min_total_bp_warn=501")).level
        == "warn"
    )


# ---------------------------------------------------------------- T-PRE-06 genome.contiguity


def test_T_PRE_06_pass() -> None:
    """§7.0: a complete genome with two small plasmids raises no contiguity warning."""
    records = [
        rec("chr", acgt(4_600_000), circular=True),
        rec("p1", acgt(60_000), circular=True),
        rec("p2", acgt(5_000), circular=True),
    ]
    c = check(records, "genome.contiguity", cfg())
    assert c.level == "ok"
    assert c.value == {"n50": 4_600_000, "largest": 4_600_000, "median": 60_000}


def test_T_PRE_06_fail() -> None:
    # largest contig < 1 Mb (N50 fine)
    c = check([rec("a", acgt(900_000)), rec("b", acgt(800_000))], "genome.contiguity", cfg())
    assert c.level == "warn"
    assert c.value["largest"] == 900_000
    assert c.value["n50"] == 900_000
    # N50 < 100 kb (largest fine)
    records = [rec("big", acgt(1_200_000))] + [rec(f"s{i}", acgt(50_000)) for i in range(30)]
    c = check(records, "genome.contiguity", cfg())
    assert c.level == "warn"
    assert c.value["n50"] == 50_000
    assert c.value["largest"] == 1_200_000
    assert c.value["median"] == 50_000


def test_contiguity_thresholds_from_config() -> None:
    records = [
        rec("a", acgt(400)),
        rec("b", acgt(200)),
        rec("c", acgt(100)),
    ]  # N50 400, largest 400
    ok = cfg("preflight.min_n50_warn=400", "preflight.min_largest_contig_warn=400")
    assert check(records, "genome.contiguity", ok).level == "ok"
    assert (
        check(
            records,
            "genome.contiguity",
            cfg("preflight.min_n50_warn=401", "preflight.min_largest_contig_warn=1"),
        ).level
        == "warn"
    )
    assert (
        check(
            records,
            "genome.contiguity",
            cfg("preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=401"),
        ).level
        == "warn"
    )


def test_n50_edge_cases() -> None:
    with pytest.raises(ValueError):
        n50([])
    assert n50([5]) == 5
    assert n50([3, 3, 2]) == 3
    assert n50([2, 3, 3]) == 3  # order of input does not matter
    assert n50([4_600_000, 60_000, 5_000]) == 4_600_000
    assert n50([1, 1, 1, 1]) == 1
    assert n50([10, 10]) == 10  # cumulative sum reaches exactly half


# ---------------------------------------------------------------- T-PRE-07 genome.gaps


def test_T_PRE_07_pass() -> None:
    c = check([rec("c", "ACGT" + "N" * 9 + "ACGT")], "genome.gaps")
    assert c.level == "ok"
    assert summary([rec("c", "ACGT" + "N" * 9 + "ACGT")]).n_runs == ()


def test_T_PRE_07_fail() -> None:
    records = [rec("c", "ACG" + "N" * 10 + "T")]
    assert summary(records).n_runs == (("c", 3, 13),)  # 0-based half-open internally
    c = check(records, "genome.gaps")
    assert c.level == "warn"
    assert c.value == [{"contig": "c", "start": 4, "end": 13}]  # 1-based closed in results

    records = [rec("c", "ACRT")]
    assert summary(records).ambiguous == (("c", 2, 3, "R"),)
    c = check(records, "genome.gaps")
    assert c.level == "warn"
    assert c.value == [{"contig": "c", "start": 3, "end": 3, "code": "R"}]


def test_n_run_min_bp_from_config() -> None:
    records = [rec("c", "ACGT" + "N" * 5 + "ACGT")]
    assert check(records, "genome.gaps", small_cfg()).level == "ok"
    c = check(records, "genome.gaps", small_cfg("preflight.n_run_min_bp=5"))
    assert c.level == "warn"
    assert c.value == [{"contig": "c", "start": 5, "end": 9}]


def test_n_run_lowercase_and_across_line_break(tmp_path: Path) -> None:
    got = from_file(tmp_path, ">c\nACGTnnnnn\nNNNNNacgt\n")
    assert got["genome.gaps"].value == [{"contig": "c", "start": 5, "end": 14}]


def test_n_run_at_origin_of_circular_contig_reported_as_two_runs() -> None:
    records = [rec("c", "N" * 10 + "ACGTACGT" + "N" * 10, circular=True)]
    assert summary(records).n_runs == (("c", 0, 10), ("c", 18, 28))
    assert check(records, "genome.gaps").value == [
        {"contig": "c", "start": 1, "end": 10},
        {"contig": "c", "start": 19, "end": 28},
    ]


def test_gap_positions_use_safe_ids() -> None:
    c = check([rec("chr:1", "ACRT")], "genome.gaps")
    assert c.value == [{"contig": "chr_1", "start": 3, "end": 3, "code": "R"}]


# ---------------------------------------------------------------- T-PRE-08 genome.gc


@pytest.mark.parametrize(("seq", "gc"), [("ACGT", 0.5), ("AAAC", 0.25), ("ACCC", 0.75)])
def test_T_PRE_08_pass(seq: str, gc: float) -> None:
    """GC bounds are inclusive: 25 % and 75 % pass."""
    c = check([rec("c", seq)], "genome.gc")
    assert c.level == "ok"
    assert c.value == pytest.approx(gc)


@pytest.mark.parametrize(("seq", "gc"), [("AAAAAAAC", 0.125), ("GGGGGGGA", 0.875)])
def test_T_PRE_08_fail(seq: str, gc: float) -> None:
    c = check([rec("c", seq)], "genome.gc")
    assert c.level == "warn"
    assert c.value == pytest.approx(gc)


def test_gc_range_from_config() -> None:
    records = [rec("c", "AAAAAAAC")]  # 0.125
    assert (
        check(records, "genome.gc", small_cfg("preflight.gc_range_warn=[0.125, 0.9]")).level == "ok"
    )
    assert (
        check(records, "genome.gc", small_cfg("preflight.gc_range_warn=[0.13, 0.9]")).level
        == "warn"
    )


# ---------------------------------------------------------------- T-PRE-09 genome.topology


def test_T_PRE_09_pass() -> None:
    records = [rec("c", acgt(40), circular=True)]
    s = summary(records)
    assert resolve_topology(s, "auto") == ("circular",)
    c = check(records, "genome.topology")
    assert c.level == "info"
    assert "c" in json.dumps(c.value)


def test_T_PRE_09_fail() -> None:
    records = [rec("c", acgt(40)), rec("d", acgt(40), circular=True)]
    s = summary(records)
    assert resolve_topology(s, "auto") == ("linear", "circular")
    assert resolve_topology(s, True) == ("circular", "circular")  # complete: true → all circular
    c = check(records, "genome.topology")
    assert c.level == "info"


@pytest.mark.parametrize(
    "header", ["c circular=true", "c [topology=circular]", "c Escherichia coli circular=true"]
)
def test_topology_tags_from_file(tmp_path: Path, header: str) -> None:
    s = summary(read_genome_fasta(write_fa(tmp_path, f">{header}\nACGT\n")))
    assert s.circular_tags == (True,)
    assert resolve_topology(s, "auto") == ("circular",)


def test_topology_untagged_from_file_is_linear(tmp_path: Path) -> None:
    s = summary(read_genome_fasta(write_fa(tmp_path, ">c some description\nACGT\n")))
    assert resolve_topology(s, "auto") == ("linear",)


# ---------------------------------------------------------------- T-PRE-10 genome.small_replicons


def test_T_PRE_10_pass() -> None:
    c = check([rec("c", acgt(600_000))], "genome.small_replicons", cfg())
    assert c.level == "ok"
    assert c.value == []


def test_T_PRE_10_fail() -> None:
    records = [rec("chr", acgt(4_600_000)), rec("p1", acgt(60_000))]
    c = check(records, "genome.small_replicons", cfg())
    assert c.level == "info"
    assert c.value == ["p1"]


def test_small_replicon_bp_from_config() -> None:
    records = [rec("chr", acgt(400)), rec("p1", acgt(100))]
    c = check(records, "genome.small_replicons", small_cfg("preflight.small_replicon_bp=101"))
    assert (c.level, c.value) == ("info", ["p1"])
    c = check(records, "genome.small_replicons", small_cfg("preflight.small_replicon_bp=100"))
    assert (c.level, c.value) == ("ok", [])


# ---------------------------------------------------------------- replicon classification


def _classify(spec: list[tuple[str, int, bool]], **kw: int) -> list[Replicon]:
    # sequences of the stated length (in-memory; GC irrelevant)
    s = summary([rec(i, "A" * n, circular=circ) for i, n, circ in spec])
    return classify_replicons(
        s, resolve_topology(s, "auto"), small_bp=kw.pop("small_bp", 500_000), **kw
    )


def test_classify_chromosome_linear_large_plus_circular_plasmid() -> None:
    reps = _classify([("chr", 4_600_000, False), ("p1", 60_000, True)])
    assert [(r.id, r.length, r.topology, r.kind) for r in reps] == [
        ("chr", 4_600_000, "linear", "chromosome"),
        ("p1", 60_000, "circular", "plasmid"),
    ]


def test_classify_largest_circular_below_1mb_is_chromosome() -> None:
    reps = _classify([("c1", 900_000, True), ("c2", 50_000, False)])
    assert [r.kind for r in reps] == ["chromosome", "unknown"]


def test_classify_not_largest_and_not_small_is_unknown() -> None:
    reps = _classify([("c1", 900_000, False), ("c2", 600_000, True)])
    assert [r.kind for r in reps] == ["unknown", "unknown"]


def test_classify_tie_for_largest_goes_to_first_contig() -> None:
    reps = _classify(
        [("a", 1_000, True), ("b", 1_000, True)], small_bp=500, chromosome_min_bp=10_000
    )
    assert [r.kind for r in reps] == ["chromosome", "unknown"]


def test_classify_chromosome_min_bp_threshold() -> None:
    """The chromosome-size threshold is a parameter; preflight passes min_largest_contig_warn."""
    spec = [("chr", 2_000, False), ("p1", 100, True)]
    assert [r.kind for r in _classify(spec, small_bp=500, chromosome_min_bp=2_000)] == [
        "chromosome",
        "plasmid",
    ]
    assert [r.kind for r in _classify(spec, small_bp=500, chromosome_min_bp=2_001)] == [
        "unknown",
        "plasmid",
    ]


def test_classify_complete_true_makes_all_circular() -> None:
    s = summary([rec("c1", "A" * 900), rec("c2", "A" * 50)])
    reps = classify_replicons(s, resolve_topology(s, True), small_bp=500, chromosome_min_bp=10_000)
    assert [(r.topology, r.kind) for r in reps] == [
        ("circular", "chromosome"),
        ("circular", "plasmid"),
    ]


def test_topology_check_value_uses_min_largest_contig_warn_for_chromosome() -> None:
    """ADR-0005 item 3: the chromosome-size threshold reuses preflight.min_largest_contig_warn."""
    records = [rec("chr", acgt(2_000)), rec("p1", acgt(100), circular=True)]
    hi = json.dumps(
        check(records, "genome.topology", cfg("preflight.min_largest_contig_warn=2001")).value
    )
    lo = json.dumps(
        check(records, "genome.topology", cfg("preflight.min_largest_contig_warn=2000")).value
    )
    assert "chromosome" in lo
    assert "chromosome" not in hi


# ---------------------------------------------------------------- realistic genome (giapsim)


def test_giapsim_genome_passes(tmp_path: Path) -> None:
    sim = giapsim.simulate(
        giapsim.SimConfig(genome_length=120_000), seed=7, outdir=tmp_path / "sim"
    )
    records = read_genome_fasta(sim.genome_fasta)
    config = cfg("preflight.min_largest_contig_warn=100000", "preflight.small_replicon_bp=100000")
    got = checks_for(records, config)
    assert set(got) == set(RULE_IDS)
    levels = {rid: c.level for rid, c in got.items()}
    assert levels == {rid: ("info" if rid == "genome.topology" else "ok") for rid in RULE_IDS}
    assert got["genome.total_length"].value == 120_000
    s = summary(records)
    assert resolve_topology(s, "auto") == ("circular",)
    reps = classify_replicons(
        s, resolve_topology(s, "auto"), small_bp=100_000, chromosome_min_bp=100_000
    )
    assert [(r.id, r.kind) for r in reps] == [("chr1", "chromosome")]
    assert worst(got.values()) == "info"
