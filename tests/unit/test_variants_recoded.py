"""Variants over annotated recoded Sec/Pyl codons (M7-16c; spec rule 10 with rule 22
"``sixframe_index`` and ``variants_step`` both call ``annotated_stops(gm, genome,
recoded=True)``", rule 11, rule 26, example 2's variants row (toy shape), example 30's variants
part; test plan C7).

* M7-01 ``align_hits`` accepts a genome ``*`` at a codon in
  ``site_keys(annotated_stops(gm, genome, recoded=True)[0])`` and still raises for a ``*`` at an
  internal codon outside that set;
* ``run_variants`` builds that set with ``recoded=True``: ALT ``U`` over a recoded TGA is
  ``candidate`` (one distinct unique peptide) with flags ``stop_codon;vcf_unrepresentable`` and no
  VCF record; ``O`` over a recoded TAG likewise; a ``W`` over the recoded TGA is an ordinary ALT
  ``W`` with ``snv;stop_codon`` and a VCF record TGA -> TGG (rule 26).

World: the 16b toy ``G1`` (``MSDKPLTEAGRWQLNEFK`` + TGA + ``VPAGHTLDMRK`` + TAA, 93 nt, linear)
annotated as **one** CDS ``[0, 93)`` whose protein carries ``U`` at residue 18, so the TGA at
``[54, 57)`` is an internal recoded codon, not a terminal stop. The minus world is its reverse
complement (``G1`` ``-`` ``[0, 93)``, recoded codon ``[36, 39)``).
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import pandas as pd
import pytest
from helpers.readthrough_toys import back_translate, gene, model, rc

from gippa.config import GippaConfig
from gippa.model import Replicon
from gippa.steps.variants import run_variants

BODY, FLANK = "MSDKPLTEAGRWQLNEFK", "VPAGHTLDMRK"
N = 3 * (len(BODY) + 1 + len(FLANK) + 1)  # 93


def world(codon: str = "TGA", letter: str = "U", strand: str = "+") -> tuple[dict[str, str], Any]:
    nt = back_translate(BODY + "*" + FLANK + "*", {18: codon, 30: "TAA"})
    assert len(nt) == N
    reps = {"c1": Replicon("c1", N, "linear")}
    aa = BODY + letter + FLANK
    g1 = replace(gene("G1", "c1", 0, N, strand), aa=aa)
    return ({"c1": nt} if strand == "+" else {"c1": rc(nt)}), model(reps, [g1])


def window(a: int, b: int, strand: str = "+") -> tuple[int, int]:
    s, e = 3 * a, 3 * b
    return (s, e) if strand == "+" else (N - e, N - s)


def _hit(pep: str, s: int, e: int, strand: str = "+", n_diff: int = 1) -> dict[str, Any]:
    return {
        "peptide": pep,
        "peptide_len": len(pep),
        "source": "denovo",
        "contig": "c1",
        "start": s,
        "end": e,
        "strand": strand,
        "n_diff": n_diff,
        "diff_types": "substitution" if n_diff else pd.NA,
        "placement_q": 0.0,
        "uniqueness": "unique",
    }


def ex(alt: str, strand: str = "+") -> dict[str, Any]:
    s, e = window(11, 26, strand)
    return _hit("WQLNEFK" + alt + "VPAGHTL", s, e, strand)


def run(rows: list[dict[str, Any]], codon: str = "TGA", letter: str = "U", strand: str = "+"):  # type: ignore[no-untyped-def]
    genome, gm = world(codon, letter, strand)
    hits = pd.DataFrame(rows)
    peptides = pd.DataFrame({"peptide": hits["peptide"].unique(), "search_db": "same"})
    cfg = GippaConfig.model_validate({"genome_relation": "same", "sample_name": "S1"})
    return run_variants(hits, None, peptides, gm, genome, cfg)


def only(table: pd.DataFrame) -> dict[str, Any]:
    assert len(table) == 1, table
    return table.to_dict("records")[0]


# ---------------------------------------------------------------- the site set


def test_recoded_site_set_of_the_toy() -> None:
    from gippa.mapping.readthrough import annotated_stops, site_keys

    genome, gm = world()
    sites, _ = annotated_stops(gm, genome, recoded=True)
    assert site_keys(sites) == frozenset({("c1", "+", 54), ("c1", "+", 90)})
    gm_m = world(strand="-")
    sites_m, _ = annotated_stops(gm_m[1], gm_m[0], recoded=True)
    assert site_keys(sites_m) == frozenset({("c1", "-", 36), ("c1", "-", 0)})


def test_align_hits_accepts_a_star_at_a_recoded_codon() -> None:
    from gippa.mapping.readthrough import annotated_stops, site_keys
    from gippa.variants.align import align_hits

    genome, gm = world()
    hits = pd.DataFrame([ex("U")])
    keys = site_keys(annotated_stops(gm, genome, recoded=True)[0])
    (a,) = align_hits(hits, None, genome, gm.replicons, stop_sites=keys)
    diffs = [o for o in a.ops if o.kind not in ("match", "il")]
    assert [(o.kind, o.codons, o.ref, o.alt) for o in diffs] == [("substitution", (54,), "*", "U")]
    with pytest.raises(ValueError):  # the 16a/16b set (terminal only) does not hold it
        align_hits(
            hits, None, genome, gm.replicons, stop_sites=site_keys(annotated_stops(gm, genome)[0])
        )


def test_a_star_at_an_internal_codon_outside_the_set_still_raises() -> None:
    """The same TGA, but the protein carries ``A`` there: no recoded site, so ``*`` raises."""
    with pytest.raises(ValueError, match=r"\*"):
        run([ex("W")], letter="A")


def test_a_mismatched_letter_gives_no_site_and_raises() -> None:
    """Rule 21 (e): a ``U`` over TAA is skipped, so the ``*`` is not consumable."""
    with pytest.raises(ValueError, match=r"\*"):
        run([ex("U")], codon="TAA")


# ---------------------------------------------------------------- example 2's row (toy shape)


@pytest.mark.parametrize("strand", ["+", "-"])
def test_u_over_a_recoded_tga_is_an_unrepresentable_candidate_alt(strand: str) -> None:
    _pile, table, records = run([ex("U", strand)], strand=strand)
    r = only(table)
    cs = 54 if strand == "+" else 36
    assert (int(r["codon_start"]), r["strand"], r["gene"]) == (cs, strand, "G1")
    assert (r["ref_codon"], r["ref_aa"], r["alt_aa"], r["call"], r["confidence"]) == (
        "TGA",
        "*",
        "U",
        "ALT",
        "candidate",
    )
    assert int(r["n_alt"]) == 1 and int(r["aa_pos"]) == 19
    assert set(r["flags"]) == {"stop_codon", "vcf_unrepresentable"}
    assert records == []


def test_o_over_a_recoded_tag_is_an_unrepresentable_alt() -> None:
    _pile, table, records = run([ex("O")], codon="TAG", letter="O")
    r = only(table)
    assert (r["ref_codon"], r["ref_aa"], r["alt_aa"], r["call"]) == ("TAG", "*", "O", "ALT")
    assert set(r["flags"]) == {"stop_codon", "vcf_unrepresentable"}
    assert records == []


# ---------------------------------------------------------------- rule 26 / example 30


@pytest.mark.parametrize("strand", ["+", "-"])
def test_w_over_a_recoded_tga_is_an_ordinary_w_alt(strand: str) -> None:
    _pile, table, records = run([ex("W", strand)], strand=strand)
    r = only(table)
    assert (r["ref_codon"], r["ref_aa"], r["alt_aa"], r["call"]) == ("TGA", "*", "W", "ALT")
    assert set(r["flags"]) == {"snv", "stop_codon"}
    (rec,) = records
    if strand == "+":
        assert (rec.pos0 + 1, rec.ref, rec.alts) == (55, "TGA", ("TGG",))
    else:
        assert (rec.pos0, rec.ref, rec.alts) == (36, rc("TGA"), (rc("TGG"),))


def test_the_terminal_stop_of_the_recoded_gene_is_still_a_site() -> None:
    """``PAGHTLDMRK`` + ``Q`` over G1's own terminal TAA at p = 30 (codon ``[90, 93)``)."""
    s, e = window(20, 31)
    _pile, table, _ = run([_hit("PAGHTLDMRKQ", s, e)])
    r = only(table)
    assert (int(r["codon_start"]), r["ref_aa"], r["alt_aa"]) == (90, "*", "Q")
