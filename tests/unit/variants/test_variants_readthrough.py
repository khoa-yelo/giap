"""Variants on read-through rows (M7-16b rules 10-11; spec examples 1, 9, 12 and 19, edge cases
"Minus strand" and "TAA / TAG / TGA"; Q4 default).

* ``align_hit``/``align_hits`` gain ``stop_sites: frozenset[(contig, strand, codon_start)]``; a
  genome ``*`` is accepted only at a codon in that set and only as a substitution of cost 1; any
  other ``*`` (or ``X``) still raises; the recomputed cost must equal ``n_diff``.
* ``run_variants`` (the ``variants`` step's chain) builds that set itself from its gene model and
  genome (``site_keys(annotated_stops(gm, genome))``); its signature is unchanged.
* REF ``*`` with a standard ALT is an ordinary ALT with ``stop_codon``; ``nt_class`` from the stop
  codon. ``U`` over TGA and ``O`` over TAG are ALT ``U``/``O`` with ``stop_codon`` and
  ``vcf_unrepresentable``, neither ``snv`` nor ``multi_nt``, and no VCF record. Every other U/O
  stays ``rare``-skipped; B/J/X/Z stay ``ambiguous``-skipped (example 9: a patched ``X``).

World: the spec's toy gene ``G1`` = ``MSDKPLTEAGRWQLNEFK`` + stop (cell p = 18, codon
``[54, 57)``) + ``VPAGHTLDMRK`` + TAA (p = 30, *not* an annotated stop) on a linear 93-nt contig;
the stop codon is chosen per test. The minus world is its reverse complement (G1 ``-`` ``[36,
93)``, stop codon ``[36, 39)``). Input is de novo with ``genome_relation = same`` unless stated.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest
from helpers.readthrough_toys import back_translate, gene, model, rc

from gippa.config import GippaConfig
from gippa.model import Replicon
from gippa.steps.variants import run_variants
from gippa.variants.vcf_records import minimal_alt_codon

BODY, FLANK = "MSDKPLTEAGRWQLNEFK", "VPAGHTLDMRK"
N = 3 * (len(BODY) + 1 + len(FLANK) + 1)  # 93


def world(stop: str = "TGA", strand: str = "+") -> tuple[dict[str, str], Any]:
    nt = back_translate(BODY + "*" + FLANK + "*", {18: stop, 30: "TAA"})
    assert len(nt) == N
    reps = {"c1": Replicon("c1", N, "linear")}
    if strand == "+":
        return {"c1": nt}, model(reps, [gene("G1", "c1", 0, 57, "+")])
    return {"c1": rc(nt)}, model(reps, [gene("G1", "c1", N - 57, N, "-")])


def window(a: int, b: int, strand: str = "+") -> tuple[int, int]:
    """Genome interval of gene-strand aa ``[a, b)``."""
    s, e = 3 * a, 3 * b
    return (s, e) if strand == "+" else (N - e, N - s)


def _cfg(relation: str = "same") -> GippaConfig:
    return GippaConfig.model_validate({"genome_relation": relation, "sample_name": "S1"})


def _hit(
    pep: str,
    s: int,
    e: int,
    strand: str = "+",
    n_diff: int = 1,
    diff_types: str | None = "substitution",
    source: str = "denovo",
) -> dict[str, Any]:
    return {
        "peptide": pep,
        "peptide_len": len(pep),
        "source": source,
        "contig": "c1",
        "start": s,
        "end": e,
        "strand": strand,
        "n_diff": n_diff,
        "diff_types": pd.NA if diff_types is None else diff_types,
        "placement_q": 0.0,
        "uniqueness": "unique",
    }


def run(
    rows: list[dict[str, Any]], stop: str = "TGA", strand: str = "+", relation: str = "same"
) -> tuple[Any, Any, Any]:
    genome, gm = world(stop, strand)
    hits = pd.DataFrame(rows)
    peptides = pd.DataFrame({"peptide": hits["peptide"].unique(), "search_db": "same"})
    return run_variants(hits, None, peptides, gm, genome, _cfg(relation))


def ex1(alt: str = "W", strand: str = "+") -> dict[str, Any]:
    """``WQLNEFK`` + alt + ``VPAGHTL`` over aa ``[11, 26)`` (example 1's window)."""
    s, e = window(11, 26, strand)
    return _hit("WQLNEFK" + alt + "VPAGHTL", s, e, strand)


def only_site(table: pd.DataFrame) -> dict[str, Any]:
    assert len(table) == 1, table
    return table.to_dict("records")[0]


# ---------------------------------------------------------------- examples 1, 12, minus strand


@pytest.mark.parametrize(
    ("relation", "flags"),
    [
        ("same", {"snv", "stop_codon"}),
        ("related", {"snv", "strain_difference", "stop_codon"}),
    ],
)
def test_example_1_ref_stop_alt_w(relation: str, flags: set[str]) -> None:
    _pile, table, records = run([ex1()], relation=relation)
    r = only_site(table)
    assert (int(r["codon_start"]), r["strand"], r["gene"]) == (54, "+", "G1")
    assert (r["ref_codon"], r["ref_aa"], r["alt_aa"], r["alt_type"]) == (
        "TGA",
        "*",
        "W",
        "substitution",
    )
    assert (r["call"], r["confidence"], int(r["n_alt"])) == ("ALT", "candidate", 1)
    assert int(r["aa_pos"]) == 19
    assert set(r["flags"]) == flags
    (rec,) = records
    assert (rec.pos0 + 1, rec.ref, rec.alts) == (55, "TGA", ("TGG",))


def test_example_12_window_ending_on_the_stop() -> None:
    s, e = window(9, 19)
    assert (s, e) == (27, 57)
    _pile, table, records = run([_hit("GRWQLNEFKW", s, e)])
    r = only_site(table)
    assert (int(r["codon_start"]), r["alt_aa"], set(r["flags"])) == (54, "W", {"snv", "stop_codon"})
    assert [(x.pos0, x.ref, x.alts) for x in records] == [(54, "TGA", ("TGG",))]


def test_minus_strand_stop_is_the_lowest_codon() -> None:
    _pile, table, records = run([ex1(strand="-")], strand="-")
    r = only_site(table)
    assert (int(r["codon_start"]), r["strand"], r["ref_codon"], r["alt_aa"]) == (
        36,
        "-",
        "TGA",
        "W",
    )
    assert set(r["flags"]) == {"snv", "stop_codon"}
    (rec,) = records
    assert (rec.pos0, rec.ref, rec.alts) == (36, rc("TGA"), (rc("TGG"),))


# ---------------------------------------------------------------- rule 11 per codon


@pytest.mark.parametrize(
    ("stop", "alt", "cls"),
    [
        ("TGA", "W", "snv"),
        ("TAA", "Q", "snv"),
        ("TAG", "Q", "snv"),
        ("TAG", "Y", "snv"),
        ("TAA", "Y", "snv"),
        ("TGA", "R", "snv"),
        ("TAA", "R", "multi_nt"),
        ("TAG", "W", "snv"),
    ],
)
def test_standard_alt_over_each_stop_codon(stop: str, alt: str, cls: str) -> None:
    _pile, table, records = run([ex1(alt)], stop=stop)
    r = only_site(table)
    assert (r["ref_codon"], r["ref_aa"], r["alt_aa"], r["call"]) == (stop, "*", alt, "ALT")
    assert set(r["flags"]) == {cls, "stop_codon"}
    want, n = minimal_alt_codon(stop, alt)
    assert int(r["nt_changes"]) == n
    (rec,) = records
    assert (rec.pos0, rec.ref, rec.alts) == (54, stop, (want,))


@pytest.mark.parametrize(("stop", "alt"), [("TGA", "U"), ("TAG", "O")])
@pytest.mark.parametrize("strand", ["+", "-"])
def test_u_over_tga_and_o_over_tag_are_unrepresentable_alts(
    stop: str, alt: str, strand: str
) -> None:
    _pile, table, records = run([ex1(alt, strand)], stop=stop, strand=strand)
    r = only_site(table)
    assert (r["ref_codon"], r["ref_aa"], r["alt_aa"], r["call"], r["confidence"]) == (
        stop,
        "*",
        alt,
        "ALT",
        "candidate",
    )
    assert set(r["flags"]) == {"stop_codon", "vcf_unrepresentable"}
    assert records == []


@pytest.mark.parametrize(("stop", "alt"), [("TAA", "U"), ("TAG", "U"), ("TGA", "O"), ("TAA", "O")])
def test_other_u_and_o_over_a_stop_stay_rare_skipped(stop: str, alt: str) -> None:
    pile, table, records = run([ex1(alt)], stop=stop)
    assert len(table) == 0 and records == []
    ns = pile.attrs["nonstandard"]
    assert ns["by_residue"] == {alt: 1} and ns["by_rule"]["rare"] == 1
    assert 54 not in {int(c) for c in pile["codon_start"]}


@pytest.mark.parametrize("letter", ["B", "J", "X", "Z"])
def test_ambiguous_letters_over_a_stop_stay_skipped(letter: str) -> None:
    pile, table, records = run([ex1(letter)])
    assert len(table) == 0 and records == []
    ns = pile.attrs["nonstandard"]
    assert ns["by_residue"] == {letter: 1} and ns["by_rule"]["ambiguous"] == 1


def test_example_9_patched_x_gives_no_alt_and_the_other_codons_are_ref() -> None:
    """Second-pass peptide ``...NEFKXVPAG...`` over the patched ``X``: the X op is skipped
    (``ambiguous``, M7-18), so there is no ALT site; the other 14 codons are covered as REF."""
    pile, table, records = run([ex1("X")])
    assert len(table) == 0 and records == []
    covered = sorted(int(c) for c in pile["codon_start"])
    assert covered == [3 * k for k in range(11, 26) if k != 18]
    assert set(pile["call"]) == {"REF"}


# ---------------------------------------------------------------- rule 10: the re-alignment


def test_a_stop_outside_the_site_set_still_raises() -> None:
    """``PAGHTLDMRK`` + ``Q`` over the TAA at p = 30, which is not G1's stop."""
    s, e = window(20, 31)
    with pytest.raises(ValueError, match=r"\*"):
        run([_hit("PAGHTLDMRKQ", s, e)])


def test_recomputed_cost_must_equal_n_diff() -> None:
    row = ex1()
    row["n_diff"] = 2
    row["diff_types"] = "substitution;substitution"
    with pytest.raises(ValueError, match="recomputed 1"):
        run([row])


def test_stop_plus_one_mismatch_costs_two() -> None:
    """Example 4c's shape: a mismatch next to the read-through residue is a second substitution."""
    s, e = window(11, 26)
    _pile, table, _ = run(
        [_hit("WQLNEFKWVPAGHTA", s, e, n_diff=2, diff_types="substitution;substitution")]
    )
    got = {(int(r["codon_start"]), r["ref_aa"], r["alt_aa"]) for r in table.to_dict("records")}
    assert got == {(54, "*", "W"), (75, "L", "A")}


def test_align_hits_takes_the_site_set() -> None:
    from gippa.variants.align import align_hits

    genome, gm = world()
    hits = pd.DataFrame([ex1()])
    with pytest.raises(ValueError):
        align_hits(hits, None, genome, gm.replicons)  # no sites: the * raises as on main
    (a,) = align_hits(hits, None, genome, gm.replicons, stop_sites=frozenset({("c1", "+", 54)}))
    diffs = [o for o in a.ops if o.kind not in ("match", "il")]
    assert [(o.kind, o.codons, o.ref, o.alt) for o in diffs] == [("substitution", (54,), "*", "W")]
    with pytest.raises(ValueError):  # the wrong strand is not the site
        align_hits(hits, None, genome, gm.replicons, stop_sites=frozenset({("c1", "-", 54)}))


def test_align_hit_minus_strand_site_key_is_the_lowest_coordinate() -> None:
    from gippa.variants.align import align_hits

    genome, gm = world(strand="-")
    hits = pd.DataFrame([ex1(strand="-")])
    (a,) = align_hits(hits, None, genome, gm.replicons, stop_sites=frozenset({("c1", "-", 36)}))
    diffs = [o for o in a.ops if o.kind not in ("match", "il")]
    assert [(o.codons, o.ref, o.alt) for o in diffs] == [((36,), "*", "W")]


def test_example_19_strain_stop_loss_shape() -> None:
    """Rule 20's other-strain C-terminal peptide: TAA read as ``Y`` (reference TAC), window
    ending on the stop, ``related``: ALT ``Y``, ``snv;strain_difference;stop_codon``,
    ``candidate``."""
    s, e = window(9, 19)
    _pile, table, records = run(
        [_hit("GRWQLNEFKY", s, e, source="database")], stop="TAA", relation="related"
    )
    r = only_site(table)
    assert (r["alt_aa"], r["confidence"]) == ("Y", "candidate")
    assert set(r["flags"]) == {"snv", "strain_difference", "stop_codon"}
    assert [(x.pos0, x.ref, x.alts) for x in records] == [(54, "TAA", ("TAC",))]
