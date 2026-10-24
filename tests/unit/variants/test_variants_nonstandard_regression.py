"""M7-18 regression: ``run_variants`` and the ``variants`` step on peptides with B, J, O, U, X or
Z (docs/specs/M7-18-variants-ambiguous-residues.md: "The bug", rules 4-9, examples 1-4, 6 and 9
end to end, acceptance criteria 1-3 and the step event; lead decisions 1-3).

The M7-15 K-12 run failed in the ``variants`` step with ``ValueError: no table-11 codon encodes
'X'`` (``minimal_alt_codon``). These rows go through the real chain (``align_hits`` recomputes
each op from the genome) on the linear 66-nt ``c1`` of ``helpers.nonstandard_genome``.

The step test drives ``variants_step`` through ``registry.run_step`` with the genome FASTA and
the C9 ``peptides`` Parquet on disk; the ``peptide_hits.tsv`` reader and the gene-model reader
are replaced by the in-memory frames (their contracts are tested elsewhere).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.nonstandard_genome import P, linear, with_x

from giap.config import GiapConfig
from giap.steps.variants import run_variants

LIN = linear()
GENOME = LIN.genome()
EVENT = "variants.nonstandard_residues"


def _cfg(relation: str = "same") -> GiapConfig:
    return GiapConfig.model_validate({"genome_relation": relation, "sample_name": "S1"})


def _hit(pep: str, diff_types: str | None) -> dict[str, Any]:
    n = 0 if diff_types is None else len(diff_types.split(";"))
    return {
        "peptide": pep, "peptide_len": len(pep), "source": "denovo", "contig": "c1",
        "start": 0, "end": 36, "strand": "+", "n_diff": n,
        "diff_types": pd.NA if diff_types is None else diff_types,
        "placement_q": 0.0, "uniqueness": "unique",
    }  # fmt: skip


def _frames(rows: list[tuple[str, str | None]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    hits = pd.DataFrame([_hit(p, d) for p, d in rows])
    peptides = pd.DataFrame({"peptide": [p for p, _ in rows], "search_db": "same"})
    return hits, peptides


def _run(rows: list[tuple[str, str | None]], relation: str = "same") -> tuple[Any, Any, Any]:
    hits, peptides = _frames(rows)
    return run_variants(hits, None, peptides, LIN.gm(), GENOME, _cfg(relation))


def _codons(pile: pd.DataFrame) -> list[int]:
    return sorted(int(c) for c in pile["codon_start"])


SUB = "substitution"


@pytest.mark.parametrize(
    ("letter", "i"),
    [
        ("X", 2),  # over C (example 1)
        ("X", 0),  # first residue
        ("X", 11),  # last residue
        ("U", 2),  # over the C codon TGT (example 6)
        ("B", 11),  # over D
        ("Z", 4),  # over E
        ("O", 7),  # over the H codon (sense)
    ],
)
@pytest.mark.parametrize("relation", ["same", "related"])
def test_nonstandard_letter_does_not_crash_run_variants(letter: str, i: int, relation: str) -> None:
    """On main at 94bb0029 this raises ``ValueError: no table-11 codon encodes '<letter>'``."""
    pep = with_x(P, i, letter)
    pile, table, records = _run([(pep, SUB)], relation)
    assert _codons(pile) == [3 * k for k in range(12) if k != i]
    assert set(pile["call"]) == {"REF"}
    assert len(table) == 0 and records == []
    got = pile.attrs["nonstandard"]
    rule = "rare" if letter in "OU" else "ambiguous"
    assert got["skipped_sites"] == got["skipped_codons"] == got["skipped_peptides"] == 1
    assert got["by_residue"] == {letter: 1} and got["by_rule"][rule] == 1
    assert got["examples"] == [pep]


def test_ex4_literal_j_over_an_i_codon_gives_no_record_with_alt_equal_to_ref() -> None:
    """On main this gives a VCF record with ALT == REF (``minimal_alt_codon('ATT', 'J')``)."""
    pep = "AKCREFGHJLWD"
    pile, table, records = _run([(pep, SUB)])
    assert all(r.ref not in r.alts for r in records), "VCF record with ALT == REF"
    assert 24 not in _codons(pile)
    assert len(table) == 0 and records == []
    assert pile.attrs["nonstandard"]["by_residue"] == {"J": 1}


def test_ex2_ref_peptide_and_x_peptide() -> None:
    pile, table, records = _run([(P, None), ("AKXREFGHILWD", SUB)])
    r = pile[pile["codon_start"] == 6].iloc[0]
    assert (int(r["n_ref"]), int(r["n_il"]), r["call"]) == (1, 0, "REF")
    assert len(table) == 0 and records == []
    assert pile.attrs["nonstandard"]["skipped_sites"] == 1


def test_ex3_skipped_peptide_is_not_in_ad_or_dp() -> None:
    _, table, records = _run([("AKXREFGHILWD", SUB), ("AKDREFGHILWD", SUB)])
    (rec,) = records
    assert (rec.pos0, rec.ad, rec.dp) == (6, (0, 1), 1)
    (row,) = table.to_dict("records")
    assert (row["alt_aa"], row["call"], row["confidence"], int(row["n_alt"])) == (
        "D", "ALT", "candidate", 1,
    )  # fmt: skip
    _, alone_table, alone_records = _run([("AKDREFGHILWD", SUB)])
    assert alone_records == records
    assert alone_table.equals(table)


def test_ex9_swap_and_x() -> None:
    """On main this raises ``ValueError: no table-11 codon encodes 'X'``."""
    pep = "KAXREFGHILWD"
    pile, table, records = _run([(pep, "swap;substitution")])
    assert 6 not in _codons(pile) and 3 not in _codons(pile)
    (rec,) = records
    assert (rec.pos0, rec.ref, rec.alts) == (0, "GCTAAA", ("AAAGCA",))
    (row,) = table.to_dict("records")
    assert (row["alt_type"], row["confidence"]) == ("swap", "candidate")
    assert pile.attrs["nonstandard"]["skipped_sites"] == 1


# --------------------------------------------------------------------------- the step event


def _step(tmp: Path, rows: list[tuple[str, str | None]], monkeypatch: pytest.MonkeyPatch) -> Any:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    import giap.steps.variants as sv

    hits, peptides = _frames(rows)
    monkeypatch.setattr(sv, "_peptide_hits", lambda _p: hits.copy())
    monkeypatch.setattr(sv, "read_gene_model", lambda *_a, **_k: LIN.gm())
    (tmp / "mapping").mkdir()
    (tmp / "genes.parquet").write_bytes(b"")
    (tmp / "replicons.parquet").write_bytes(b"")
    (tmp / "genome.fna").write_text(f">c1\n{GENOME['c1']}\n")
    peptides.to_parquet(tmp / "peptides.parquet")
    inputs = {
        "mapping": tmp / "mapping", "genome": tmp / "genome.fna",
        "genes": tmp / "genes.parquet", "replicons": tmp / "replicons.parquet",
        "peptides": tmp / "peptides.parquet",
    }  # fmt: skip
    outputs = {
        "pileup": tmp / "out" / "pileup.tsv.gz", "variants": tmp / "out" / "variants.tsv",
        "variants_vcf": tmp / "out" / "variants.vcf",
        "variants_qc": tmp / "out" / "variants_qc.parquet",
    }  # fmt: skip
    records: list[logging.LogRecord] = []

    class _Rec(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(record)

    h = _Rec(level=logging.DEBUG)
    log = logging.getLogger("giap.steps.variants")
    old = log.level
    log.addHandler(h)
    log.setLevel(logging.DEBUG)
    try:
        reg.run_step("variants", _cfg(), inputs, outputs)
    finally:
        log.removeHandler(h)
        log.setLevel(old)
    return records


def _events(records: list[logging.LogRecord]) -> list[logging.LogRecord]:
    return [r for r in records if getattr(r, "event", None) == EVENT]


def test_step_emits_one_info_event_when_a_site_is_skipped(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """On main the step fails with ``ValueError: no table-11 codon encodes 'X'``."""
    pep = "AKXREFGHILWD"
    (ev,) = _events(_step(tmp_path, [(P, None), (pep, SUB)], monkeypatch))
    assert ev.levelno == logging.INFO
    assert ev.data == {  # type: ignore[attr-defined]
        "by_residue": {"X": 1},
        "by_rule": {"ambiguous": 1, "rare": 0},
        "examples": [pep],
        "skipped_codons": 1,
        "skipped_peptides": 1,
        "skipped_sites": 1,
    }
    assert (tmp_path / "out" / "variants.vcf").exists()


def test_step_emits_no_event_on_all_standard_input(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import giap.steps.variants as sv

    assert sv.NONSTANDARD_EVENT == EVENT
    assert _events(_step(tmp_path, [(P, None), ("AKDREFGHILWD", SUB)], monkeypatch)) == []
