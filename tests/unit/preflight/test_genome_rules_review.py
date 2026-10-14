"""Review follow-ups for the genome preflight rules (issue M1-09; GIAP_SPEC.md §7.0, §7.1.2).

Pins lead decisions taken after code review: mutation gaps, list truncation, merged IUPAC
ambiguity runs with their own cap, a single source for the N-run threshold, and the
``genome.small_replicons`` exclusion of contigs classified as chromosome.
"""

from __future__ import annotations

import re
import time
from collections.abc import Sequence

import pytest

from giap.config import GiapConfig, resolve_config
from giap.io.fasta import GenomeRecord, sanitize_id
from giap.preflight.checks import Check
from giap.preflight.genome import (
    LIST_CAP,
    classify_replicons,
    genome_checks,
    n50,
    resolve_topology,
    summarize_genome,
)

SMALL = [
    "preflight.min_total_bp_warn=1",
    "preflight.min_n50_warn=1",
    "preflight.min_largest_contig_warn=1",
    "preflight.small_replicon_bp=1",
]


def rec(original_id: str, seq: str, *, circular: bool = False) -> GenomeRecord:
    return GenomeRecord(
        id=sanitize_id(original_id),
        original_id=original_id,
        header=original_id + (" circular=true" if circular else ""),
        sequence=seq.upper(),
        circular_tag=circular,
    )


def acgt(n: int) -> str:
    return ("ACGT" * (n // 4 + 1))[:n]


def cfg(*overrides: str) -> GiapConfig:
    return resolve_config(overrides=list(overrides))


def checks_for(records: Sequence[GenomeRecord], config: GiapConfig) -> dict[str, Check]:
    s = summarize_genome(records, config.preflight.n_run_min_bp)
    return {c.rule_id: c for c in genome_checks(s, config)}


def mentions(text: str, n: int) -> bool:
    """``text`` contains ``n`` as a whole number, with or without thousands separators."""
    return any(re.search(rf"(?<![\d,]){re.escape(f)}(?![\d,])", text) for f in {str(n), f"{n:,}"})


def ambiguity_entries(c: Check) -> list[dict[str, object]]:
    assert isinstance(c.value, list)
    return [e for e in c.value if isinstance(e, dict) and "code" in e]


def n_run_entries(c: Check) -> list[dict[str, object]]:
    assert isinstance(c.value, list)
    return [e for e in c.value if isinstance(e, dict) and "code" not in e]


# ---------------------------------------------------------------- 1. mutation gaps


def test_n50_odd_total_takes_first_length_reaching_half() -> None:
    assert n50([2, 1, 1]) == 2


def test_circular_contig_of_exactly_small_replicon_bp_is_not_plasmid() -> None:
    s = summarize_genome(
        [rec("chr", acgt(2000), circular=True), rec("p", acgt(1000), circular=True)], 10
    )
    reps = classify_replicons(s, resolve_topology(s, "auto"), small_bp=1000)
    assert [r.kind for r in reps] == ["chromosome", "unknown"]
    reps = classify_replicons(s, resolve_topology(s, "auto"), small_bp=1001)
    assert reps[1].kind == "plasmid"


def test_contiguity_all_empty_records_warns_with_zeros() -> None:
    c = checks_for([rec("a", ""), rec("b", "")], cfg(*SMALL))["genome.contiguity"]
    assert c.level == "warn"
    assert c.value == {"n50": 0, "largest": 0, "median": 0}


def test_gc_without_acgt_bases_warns_with_null_value() -> None:
    c = checks_for([rec("a", "N" * 50)], cfg(*SMALL))["genome.gc"]
    assert c.level == "warn"
    assert c.value is None
    assert c.message == "no A/C/G/T bases"


# ---------------------------------------------------------------- 2. truncation


def test_renamed_ids_truncated_to_cap_and_message_states_total() -> None:
    total = LIST_CAP + 234
    records = [rec(f"c|{i}", acgt(8)) for i in range(total)]
    c = checks_for(records, cfg(*SMALL))["genome.ids_renamed"]
    assert c.level == "warn"
    assert isinstance(c.value, list) and len(c.value) == LIST_CAP
    assert c.value[0] == {"from": "c|0", "to": "c_0"}
    assert mentions(c.message, total)


def test_n_runs_truncated_to_cap_and_message_states_total() -> None:
    total = LIST_CAP + 50
    seq = ("ACGT" + "N" * 10) * total + "ACGT"
    c = checks_for([rec("a", seq)], cfg(*SMALL))["genome.gaps"]
    assert c.level == "warn"
    assert isinstance(c.value, list) and len(c.value) == LIST_CAP
    assert c.value[0] == {"contig": "a", "start": 5, "end": 14}
    assert mentions(c.message, total)


# ---------------------------------------------------------------- 3. ambiguity runs


def test_adjacent_identical_iupac_codes_merge_into_one_run() -> None:
    c = checks_for([rec("a", "ACRRRGT")], cfg(*SMALL))["genome.gaps"]
    assert ambiguity_entries(c) == [{"contig": "a", "start": 3, "end": 5, "code": "R"}]


def test_different_adjacent_codes_and_single_codes_are_separate_runs() -> None:
    c = checks_for([rec("a", "RRYYAKA"), rec("b", "W")], cfg(*SMALL))["genome.gaps"]
    assert ambiguity_entries(c) == [
        {"contig": "a", "start": 1, "end": 2, "code": "R"},
        {"contig": "a", "start": 3, "end": 4, "code": "Y"},
        {"contig": "a", "start": 6, "end": 6, "code": "K"},
        {"contig": "b", "start": 1, "end": 1, "code": "W"},
    ]


def test_alternating_ry_5mb_is_fast_and_capped() -> None:
    n = 5_000_000
    config = cfg(*SMALL)
    t0 = time.perf_counter()
    s = summarize_genome([rec("a", "RY" * (n // 2))], config.preflight.n_run_min_bp)
    out = {c.rule_id: c for c in genome_checks(s, config)}
    elapsed = time.perf_counter() - t0
    # loose 10x guard; the 2 s budget is enforced in benchmarks/test_preflight_genome.py
    assert elapsed < 20.0, f"took {elapsed:.2f} s"
    assert len(s.ambiguous) <= LIST_CAP
    c = out["genome.gaps"]
    assert c.level == "warn"
    assert len(ambiguity_entries(c)) <= LIST_CAP
    assert mentions(c.message, n)


# ---------------------------------------------------------------- 4. separate caps


def test_iupac_entries_survive_n_run_overflow_and_both_totals_reported() -> None:
    n_runs = 1500
    seq = ("ACGT" + "N" * 10) * n_runs + "ACGT" + "RAYAK"
    c = checks_for([rec("a", seq)], cfg(*SMALL))["genome.gaps"]
    assert c.level == "warn"
    assert len(n_run_entries(c)) == LIST_CAP
    base = 14 * n_runs + 4
    assert ambiguity_entries(c) == [
        {"contig": "a", "start": base + 1, "end": base + 1, "code": "R"},
        {"contig": "a", "start": base + 3, "end": base + 3, "code": "Y"},
        {"contig": "a", "start": base + 5, "end": base + 5, "code": "K"},
    ]
    assert mentions(c.message, n_runs)
    assert mentions(c.message, 3)


# ---------------------------------------------------------------- 5. N-run threshold


@pytest.mark.parametrize("threshold", [1, 5, 25])
def test_summary_records_threshold_and_check_uses_it(threshold: int) -> None:
    config = cfg(*SMALL, f"preflight.n_run_min_bp={threshold}")
    seq = acgt(8) + "N" * threshold + acgt(8) + "N" * (threshold - 1) + acgt(8)
    s = summarize_genome([rec("a", seq)], config.preflight.n_run_min_bp)
    assert s.n_run_min_bp == threshold  # type: ignore[attr-defined]
    c = {x.rule_id: x for x in genome_checks(s, config)}["genome.gaps"]
    assert f"≥ {threshold} N" in c.rule
    assert n_run_entries(c) == [{"contig": "a", "start": 9, "end": 8 + threshold}]


def test_check_threshold_comes_from_summary_not_config() -> None:
    s = summarize_genome([rec("a", acgt(8) + "N" * 5 + acgt(8))], 5)
    c = {x.rule_id: x for x in genome_checks(s, cfg(*SMALL, "preflight.n_run_min_bp=20"))}
    gaps = c["genome.gaps"]
    assert "≥ 5 N" in gaps.rule
    assert "≥ 20 N" not in gaps.rule
    assert gaps.level == "warn"


def test_n_run_min_bp_zero_rejected_by_config() -> None:
    with pytest.raises(ValueError):
        cfg("preflight.n_run_min_bp=0")


# ---------------------------------------------------------------- 6. small replicons


def test_small_circular_chromosome_not_listed_as_small_replicon() -> None:
    c = checks_for([rec("chr", acgt(300_000), circular=True)], cfg())
    assert c["genome.topology"].value == [
        {"id": "chr", "length": 300_000, "topology": "circular", "kind": "chromosome"}
    ]
    small = c["genome.small_replicons"]
    assert small.value == []
    assert small.level == "ok"


def test_small_linear_sole_contig_still_listed_as_small_replicon() -> None:
    c = checks_for([rec("ctg", acgt(300_000))], cfg())
    assert c["genome.topology"].value == [
        {"id": "ctg", "length": 300_000, "topology": "linear", "kind": "unknown"}
    ]
    small = c["genome.small_replicons"]
    assert small.value == ["ctg"]
    assert small.level == "info"
