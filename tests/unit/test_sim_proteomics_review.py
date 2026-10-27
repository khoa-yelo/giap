"""gippasim proteomics: the lead decisions after code review of PR #55 (issue M3-19).

docs/specs/M3-19-sim-peptides-db-errors.md, "Lead decisions after code review", items 1-4:

1. A wrong-peptide shortfall is never silent: the truth metadata records ``pi_wrong_realised``
   next to ``pi_wrong``, a WARNING names both, and PEP/q use the realised share.
2. Decoy sampling honours ``missed_cleavages`` (``missed_cleavages=0, digest_len=(20, 30)``).
3. ``n_loci_genome`` is the count of GIPPA's merged exact matcher (six-frame, I = L; 0 allowed) for
   correct and wrong peptides; N-terminal peptides of GTG/TTG-start genes carry
   ``alt_start_nterm = true`` and are exactly the correct peptides without an exact hit.
4. Mutant killers: a tight wrong share (the 1/(1 - pi) factor), a peptide starting exactly at the
   origin, the charge-2 share.

Conventions pinned here where the decision is silent:

* "Truth metadata" is a JSON or YAML mapping under ``<out>/truth/`` (or ``<out>/``) that holds
  ``pi_wrong`` and ``pi_wrong_realised``; the realised share is counted like example 3 (wrong
  over distinct non-decoy peptides of ``truth/peptides.tsv``).
* The WARNING goes through ``logging`` (level WARNING) or ``warnings.warn``; its text contains
  ``pi_wrong`` and the requested and realised values (any float formatting, 3 decimals suffice).
* ``alt_start_nterm`` is a column of ``truth/peptides.tsv``, ``true``/``false`` on every correct
  row (empty or ``false`` elsewhere).
* The charge-2 share is the simulator's 0.6 (accepted with PR #55).
"""

from __future__ import annotations

import logging
import math
import re
import warnings
from pathlib import Path
from typing import Any

import pytest
from helpers.peptide_oracle import (
    CODE,
    STARTS,
    Run,
    empty,
    find_truth_meta,
    genome_slice,
    gippa_exact_counts,
    mixture_fdr,
    mixture_pep,
    truthy,
    trypsin_p_peptides,
)

import gippasim

MIX = (3.0, 1.0, 0.0, 1.0)
ALL_DETECTED = (50.0, 0.0, 0.0)
FEATURES: dict[str, Any] = {
    "genome_length": 100_000,
    "n_plasmids": 1,
    "plasmid_length": (6_000, 9_000),
    "linear_contig_length": 20_000,
    "origin_crossing_gene": True,
    "n_rrna_operons": 1,
    "n_trna": 4,
    "n_pseudogenes": 3,
    "n_sorfs": 4,
    "n_nterm_extensions": 4,
    "n_novel_orfs": (("intergenic", 3), ("cds_outframe", 3), ("antisense", 3)),
}
_FLOAT = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def proteomics(**kw: Any) -> Any:
    from gippasim.proteomics import ProteomicsConfig

    return ProteomicsConfig(**kw)


def simulate(out: Path, seed: int, prot: Any, **kw: Any) -> Run:
    cfg = gippasim.SimConfig(proteomics=prot, **kw)
    res = gippasim.simulate(cfg, seed, out)
    return Run(out, {"res": res, "cfg": cfg, "seed": seed})


def distinct_labels(run: Run) -> dict[str, bool]:
    return {r["peptide"]: truthy(r["is_wrong"]) for r in run.truth if r["origin"] != "decoy"}


def correct_rows(run: Run) -> list[dict[str, str]]:
    return [r for r in run.truth if r["origin"] == "gene"]


def wrong_rows(run: Run) -> list[dict[str, str]]:
    return [r for r in run.truth if r["origin"] in ("shuffled", "other_strain")]


# ---------------------------------------------------------------- 1: shortfall


@pytest.fixture(scope="module")
def shortfall(tmp_path_factory: pytest.TempPathFactory) -> tuple[Run, list[str]]:
    """The reviewer's case: a near-identical strain cannot supply pi_wrong = 0.5 wrong peptides."""
    prot = proteomics(pi_wrong=0.5, wrong_source="other_strain")
    st = gippasim.StrainConfig(snv_rate_per_kb=0.5)
    out = tmp_path_factory.mktemp("shortfall")
    logger = logging.getLogger()
    records: list[logging.LogRecord] = []

    class _Keep(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(record)

    h = _Keep(level=logging.WARNING)
    old = logger.level
    logger.addHandler(h)
    logger.setLevel(min(old, logging.WARNING) if old else logging.WARNING)
    try:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            run = simulate(out, 41, prot, **{**FEATURES, "strain": st})
    finally:
        logger.removeHandler(h)
        logger.setLevel(old)
    msgs = [r.getMessage() for r in records if r.levelno >= logging.WARNING]
    msgs += [str(w.message) for w in caught]
    return run, msgs


def test_shortfall_realised_share_recorded(shortfall: tuple[Run, list[str]]) -> None:
    run, _ = shortfall
    labels = distinct_labels(run)
    share = sum(labels.values()) / len(labels)
    assert 0.0 < share < 0.4, share  # the fixture really is a shortfall
    meta = find_truth_meta(run.out)
    assert float(meta["pi_wrong"]) == pytest.approx(0.5)
    assert float(meta["pi_wrong_realised"]) == pytest.approx(share, abs=1e-4)
    assert float(meta["pi_wrong_realised"]) < float(meta["pi_wrong"])


def test_shortfall_warning_names_both(shortfall: tuple[Run, list[str]]) -> None:
    run, msgs = shortfall
    labels = distinct_labels(run)
    share = sum(labels.values()) / len(labels)

    def names_both(m: str) -> bool:
        nums = [float(x) for x in _FLOAT.findall(m)]
        return (
            "pi_wrong" in m
            and any(abs(x - 0.5) < 1e-9 for x in nums)
            and any(abs(x - share) < 5e-4 for x in nums)
        )

    assert any(names_both(m) for m in msgs), msgs


def test_shortfall_pep_and_q_use_realised_share(shortfall: tuple[Run, list[str]]) -> None:
    run, _ = shortfall
    pi = float(find_truth_meta(run.out)["pi_wrong_realised"])
    for r in run.target_rows():
        s = float(r["score"])
        assert float(r["pep"]) == pytest.approx(mixture_pep(s, pi, MIX), rel=1e-3, abs=1e-9), r
        assert float(r["q_value"]) == pytest.approx(mixture_fdr(s, pi, MIX), rel=1e-3, abs=1e-9)


def test_shortfall_pep_calibrated_against_truth(shortfall: tuple[Run, list[str]]) -> None:
    """Sum of PEPs = expected number of wrong target rows (within 4 SE), q <= 0.05 controls FDP."""
    run, _ = shortfall
    rs = run.target_rows()
    peps = [float(r["pep"]) for r in rs]
    n_wrong = sum(run.is_wrong(r["peptide"]) for r in rs)
    se = math.sqrt(sum(p * (1 - p) for p in peps))
    assert abs(sum(peps) - n_wrong) <= 4 * se + 1, (sum(peps), n_wrong, se)
    acc = [r for r in rs if float(r["q_value"]) <= 0.05]
    fdp = sum(run.is_wrong(r["peptide"]) for r in acc) / len(acc)
    assert fdp <= 0.05 + 3 * math.sqrt(0.05 * 0.95 / len(acc)), fdp


# ---------------------------------------------------------------- 4a: tight wrong share


@pytest.fixture(scope="module")
def wrong30(tmp_path_factory: pytest.TempPathFactory) -> Run:
    prot = proteomics(pi_wrong=0.3, detect_logit=ALL_DETECTED)
    return simulate(tmp_path_factory.mktemp("wrong30"), 42, prot, **FEATURES)


def test_wrong_share_tight(wrong30: Run) -> None:
    """pi_wrong = 0.3 gives a share of 0.3 (dropping 1/(1 - pi) would give 0.23)."""
    labels = distinct_labels(wrong30)
    n = len(labels)
    share = sum(labels.values()) / n
    assert abs(share - 0.3) <= 3 * math.sqrt(0.3 * 0.7 / n), (share, n)


def test_no_shortfall_realised_equals_requested(wrong30: Run) -> None:
    meta = find_truth_meta(wrong30.out)
    labels = distinct_labels(wrong30)
    assert float(meta["pi_wrong"]) == pytest.approx(0.3)
    assert float(meta["pi_wrong_realised"]) == pytest.approx(
        sum(labels.values()) / len(labels), abs=1e-4
    )


# ---------------------------------------------------------------- 2: missed_cleavages = 0


def test_missed_cleavages_zero_long_peptides(tmp_path: Path) -> None:
    prot = proteomics(pi_wrong=0.05, missed_cleavages=0, digest_len=(20, 30))
    run = simulate(tmp_path, 43, prot, **FEATURES)
    targets, decoys = run.target_rows(), run.decoy_rows()
    assert targets and decoys
    for r in targets:
        p = r["peptide"]
        assert 20 <= len(p) <= 30 and not any(c in "KR" for c in p[:-1]), p
    products: set[str] = set()
    for prot_ in run.proteins.values():
        products |= trypsin_p_peptides(prot_.aa, (20, 30), 0)
    for r in decoys:
        assert r["peptide"][::-1] in products, r


# ---------------------------------------------------------------- 3: n_loci_genome, alt starts


@pytest.fixture(scope="module")
def complete(tmp_path_factory: pytest.TempPathFactory) -> Run:
    prot = proteomics(pi_wrong=0.05, detect_logit=ALL_DETECTED)
    return simulate(tmp_path_factory.mktemp("complete"), 44, prot, **FEATURES)


@pytest.fixture(scope="module")
def other(tmp_path_factory: pytest.TempPathFactory) -> Run:
    prot = proteomics(pi_wrong=0.05, wrong_source="other_strain")
    st = gippasim.StrainConfig(snv_rate_per_kb=10.0, indel_rate_per_kb=0.5)
    return simulate(tmp_path_factory.mktemp("other"), 45, prot, **{**FEATURES, "strain": st})


def _alt_start(run: Run, r: dict[str, str]) -> bool:
    """The row's first residue is an initiator M read from a GTG/TTG codon."""
    s0, e0 = int(r["start"]) - 1, int(r["end"])
    nt = genome_slice(run.genome[r["contig"]], s0, e0, r["strand"])
    return r["peptide"][0] == "M" and nt[:3] in STARTS and CODE[nt[:3]] != "M"


@pytest.mark.parametrize("name", ["complete", "other"])
def test_n_loci_genome_is_gippa_exact_count(name: str, request: pytest.FixtureRequest) -> None:
    run: Run = request.getfixturevalue(name)
    rows = correct_rows(run) + wrong_rows(run)
    assert wrong_rows(run)
    counts = gippa_exact_counts(run, [r["peptide"] for r in rows])
    bad = []
    for r in rows:
        assert not empty(r["n_loci_genome"]), r
        if int(r["n_loci_genome"]) != counts[r["peptide"]]:
            bad.append((r["peptide"], r["origin"], r["n_loci_genome"], counts[r["peptide"]]))
    assert not bad, (len(bad), bad[:5])


def test_alt_start_nterm_flag(complete: Run) -> None:
    rows = correct_rows(complete)
    assert "alt_start_nterm" in rows[0], list(rows[0])
    want = [_alt_start(complete, r) for r in rows]
    assert any(want), "fixture has no GTG/TTG-start N-terminal peptide"
    bad = [
        (r["peptide"], r["alt_start_nterm"], w)
        for r, w in zip(rows, want, strict=True)
        if empty(r["alt_start_nterm"]) or truthy(r["alt_start_nterm"]) != w
    ]
    assert not bad, bad[:5]


def test_alt_start_peptides_are_exactly_those_without_exact_hit(complete: Run) -> None:
    rows = correct_rows(complete)
    counts = gippa_exact_counts(complete, [r["peptide"] for r in rows])
    alt = {r["peptide"] for r in rows if truthy(r.get("alt_start_nterm", ""))}
    no_hit = {r["peptide"] for r in rows if counts[r["peptide"]] == 0}
    assert alt and alt == no_hit, (sorted(alt - no_hit)[:5], sorted(no_hit - alt)[:5])
    for r in rows:
        if r["peptide"] in alt:
            assert int(r["n_loci_genome"]) == 0, r


def test_wrong_rows_have_integer_n_loci(other: Run) -> None:
    counts = gippa_exact_counts(other, [r["peptide"] for r in wrong_rows(other)])
    for r in wrong_rows(other):
        assert int(r["n_loci_genome"]) == counts[r["peptide"]] >= 0, r


# ---------------------------------------------------------------- 4b: peptide at the origin


def test_peptide_starting_exactly_at_the_origin(tmp_path: Path) -> None:
    """Seed 32 (20 kb) has tryptic products of the origin-crossing gene whose unwrapped start is
    exactly the contig length: they must be reported at position 1, not n + 1."""
    prot = proteomics(pi_wrong=0.0, detect_logit=ALL_DETECTED)
    run = simulate(tmp_path, 32, prot, genome_length=20_000, origin_crossing_gene=True)
    rows = correct_rows(run)
    at_origin = []
    for r in rows:
        n = len(run.genome[r["contig"]])
        assert 1 <= int(r["start"]) <= n, r
        p = run.proteins[r["locus_tag"]]
        if int(r["start"]) == 1 and p.end > n:
            at_origin.append(r)
    assert at_origin, "fixture lost its origin-start peptide (seed 32)"
    for r in at_origin:
        nt = genome_slice(run.genome[r["contig"]], 0, int(r["end"]), r["strand"])
        aa = "".join(CODE[nt[i : i + 3]] for i in range(0, len(nt), 3))
        assert aa[1:] == r["peptide"][1:], r


# ---------------------------------------------------------------- 4c: charge-2 share


def test_charge_2_share(wrong30: Run) -> None:
    rs = wrong30.target_rows()
    n = len(rs)
    share = sum(r["charge"] == "2" for r in rs) / n
    assert abs(share - 0.6) <= 4 * math.sqrt(0.6 * 0.4 / n), (share, n)
