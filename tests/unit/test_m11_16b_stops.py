"""M11-16b stop sets and recoded codons in ``mapping/readthrough.py`` and ``mapping/loci.py``
(docs/specs/M11-16-translation-tables-4-25.md items 4 (consumable-stop cells), 5.1, 5.2,
examples 5 and 6, decision D6, edge cases "minus strand", "selenoproteins under 4/25", "Pyl").

Interface read from the spec::

    annotated_stops(gm, genome, recoded=False, *, table=11) -> (sites, n_not_a_stop)
    recoded_skips(gm, genome, *, table=11) -> {"length_mismatch", "codon_mismatch"}
    LocusClassifier(gm, genome, *, table=11)
    annotate_peptide_hits(..., table=11)       # the spec's "annotate_loci"

Reading taken: the spec names ``annotate_loci``, which does not exist; the classifier entry point
that ``steps/mapping.loci`` calls is ``annotate_peptide_hits``, so it takes ``table=``. The
``recoded_skips`` keyword is implied by D6 (the step's ``recoded_skipped`` counts follow the table).
"""

from __future__ import annotations

from dataclasses import replace
from typing import Any

import pandas as pd
import pytest
from helpers.readthrough_toys import gene, model, rc
from test_mapping_outputs_readthrough import hit

from gippa.mapping import loci
from gippa.mapping.readthrough import (
    StopSite,
    annotated_stops,
    consumed_sites,
    frame_cells,
    recoded_skips,
)
from gippa.mapping.sixframe import build_sixframe
from gippa.model import Interval, Replicon

TABLES = (4, 11, 25)
EX5_GENE = "ATGGCTGCTGCTTGAGCTGCTGCTGCTTAA"  # [3, 33): M A A A (TGA) A A A A, TAA
ZERO = {"length_mismatch": 0, "codon_mismatch": 0}
SEC = StopSite("c", "+", 15, ("G1",), "TGA", "recoded")
PYL = StopSite("c", "+", 15, ("G1",), "TAG", "recoded")
TERM = StopSite("c", "+", 30, ("G1",), "TAA", "terminal")


def _put(seq: list[str], at: int, nt: str) -> None:
    for i, c in enumerate(nt):
        seq[at + i] = c


def ex5(aa: str = "MAAAUAAAA", codon4: str = "TGA") -> tuple[dict[str, str], Any]:
    """Example 5: ``c`` linear 60 nt (GCT filler), G1 cds ``+`` ``[3, 33)``, codon 4 at
    ``[15, 18)``."""
    seq = list("GCT" * 20)
    _put(seq, 3, EX5_GENE)
    _put(seq, 15, codon4)
    reps = {"c": Replicon("c", 60, "linear")}
    return {"c": "".join(seq)}, model(reps, [replace(gene("G1", "c", 3, 33, "+"), aa=aa)])


def _stops(gm: Any, genome: dict[str, str], table: int) -> tuple[tuple[StopSite, ...], int]:
    return annotated_stops(gm, genome, recoded=True, table=table)


# ---------------------------------------------------------------- example 5 (5.1)


def test_example_5_table_11_unchanged() -> None:
    genome, gm = ex5()
    assert _stops(gm, genome, 11) == ((SEC, TERM), 0)
    assert dict(recoded_skips(gm, genome, table=11)) == ZERO
    assert annotated_stops(gm, genome, recoded=True) == ((SEC, TERM), 0)  # the default is 11


@pytest.mark.parametrize("table", (4, 25))
def test_example_5_u_under_4_and_25_is_a_codon_mismatch(table: int) -> None:
    """D6: a ``U`` has no recoding codon under 4/25: no site, ``codon_mismatch`` once."""
    genome, gm = ex5()
    assert _stops(gm, genome, table) == ((TERM,), 0)
    assert dict(recoded_skips(gm, genome, table=table)) == {**ZERO, "codon_mismatch": 1}


def test_d6_codon_mismatch_counts_once_per_u_residue() -> None:
    seq = list("GCT" * 20)
    _put(seq, 3, "ATG" + "TGA" * 3 + "GCT" * 5 + "TAA")
    reps = {"c": Replicon("c", 60, "linear")}
    gm = model(reps, [replace(gene("G1", "c", 3, 33, "+"), aa="MUUUAAAAA")])
    genome = {"c": "".join(seq)}
    assert dict(recoded_skips(gm, genome, table=4)) == {**ZERO, "codon_mismatch": 3}
    assert dict(recoded_skips(gm, genome, table=11)) == ZERO


@pytest.mark.parametrize(("table", "aa"), [(4, "MAAAWAAAA"), (25, "MAAAGAAAA")])
def test_example_5_what_bakta_writes_at_4_and_25(table: int, aa: str) -> None:
    genome, gm = ex5(aa)
    assert _stops(gm, genome, table) == ((TERM,), 0)
    assert dict(recoded_skips(gm, genome, table=table)) == ZERO


@pytest.mark.parametrize("table", (4, 25))
def test_example_5_tga_terminal_codon_is_not_a_stop(table: int) -> None:
    """G2 ``[33, 45)`` = ATG GCT GCT TGA: no site, ``not_a_stop`` = 1 (rule 1)."""
    genome, _ = ex5()
    seq = list(genome["c"])
    _put(seq, 33, "ATGGCTGCTTGA")
    reps = {"c": Replicon("c", 60, "linear")}
    gm = model(reps, [gene("G2", "c", 33, 45, "+")])
    g = {"c": "".join(seq)}
    assert annotated_stops(gm, g, table=table) == ((), 1)
    assert annotated_stops(gm, g, recoded=True, table=table) == ((), 1)
    want = (StopSite("c", "+", 42, ("G2",), "TGA", "terminal"),)
    assert annotated_stops(gm, g, table=11) == (want, 0)


@pytest.mark.parametrize("table", (4, 25))
def test_minus_strand_tga_is_not_a_stop_under_4_and_25(table: int) -> None:
    """A ``-`` gene whose gene-strand stop is TGA (TCA on the plus string at ``[0, 3)``)."""
    nt = rc("ATG" + "GCT" * 8 + "TGA") + "GCT" * 10
    assert nt[0:3] == "TCA"
    reps = {"m": Replicon("m", len(nt), "linear")}
    gm = model(reps, [gene("M1", "m", 0, 30, "-")])
    assert annotated_stops(gm, {"m": nt}, table=table) == ((), 1)
    assert annotated_stops(gm, {"m": nt}, table=11) == (
        (StopSite("m", "-", 0, ("M1",), "TGA", "terminal"),),
        0,
    )


# ---------------------------------------------------------------- example 6: Pyl survives


@pytest.mark.parametrize("table", TABLES)
def test_example_6_pyl_over_tag_is_recoded_in_every_table(table: int) -> None:
    genome, gm = ex5("MAAAOAAAA", "TAG")
    assert _stops(gm, genome, table) == ((PYL, TERM), 0)
    assert dict(recoded_skips(gm, genome, table=table)) == ZERO


# ---------------------------------------------------------------- cells against the index


@pytest.mark.parametrize("table", TABLES)
def test_frame_cells_of_the_tables_sites_hold_stop(table: int) -> None:
    genome, gm = ex5("MAAAOAAAA", "TAG")
    index = build_sixframe(genome, dict(gm.replicons), table=table)
    sites, _ = _stops(gm, genome, table)
    cells = frame_cells(index, sites)
    assert [c.tolist() for c in cells] == [[5, 10], [], [], [], [], []]


@pytest.mark.parametrize("table", (4, 25))
def test_frame_cells_refuse_a_tga_site_under_4_and_25(table: int) -> None:
    """Item 4: ``frame_cells``' STOP check catches a site that cannot exist under the table."""
    genome, gm = ex5()
    index = build_sixframe(genome, dict(gm.replicons), table=table)
    with pytest.raises(ValueError, match="is not a stop"):
        frame_cells(index, (SEC, TERM))


# ---------------------------------------------------------------- 5.2: loci flags


def _peptides(peps: list[str]) -> pd.DataFrame:
    return pd.DataFrame({"peptide": peps, "search_db": "other", "protein_ids": [[] for _ in peps]})


def _annotate(rows: list[dict[str, Any]], world: Any, table: int, sites: Any) -> pd.DataFrame:
    genome, gm = world
    df = pd.DataFrame(rows)
    consumed = consumed_sites(df, sites, {"c": 60}, {"c": False})
    return loci.annotate_peptide_hits(
        df, gm, genome, _peptides(list(df["peptide"].unique())), "same",
        consumed=consumed, table=table,
    )  # fmt: skip


@pytest.mark.parametrize("table", TABLES)
def test_sec_site_only_under_11(table: int) -> None:
    """A row over codon 15 of example 5's G1 (``U`` annotated): ``sec_site`` under 11; under 4/25
    the row reads W/G there, consumes nothing and is never ``sec_site`` (5.2)."""
    world = ex5()
    sites, _ = _stops(world[1], world[0], table)
    letter = {11: "W", 4: "W", 25: "G"}[table]
    pep = "MAAA" + letter + "AAAA"
    route = "tolerant" if table == 11 else "exact"
    ann = _annotate([hit(pep, "c", 3, 30, route=route)], world, table, sites)
    flags = set(str(ann["flags"].iloc[0]).split(";")) - {""}
    assert ("sec_site" in flags) is (table == 11)
    assert "readthrough" not in flags and "pyl_site" not in flags
    assert ann["locus_class"].iloc[0] == "cds_inframe"


@pytest.mark.parametrize("table", TABLES)
def test_pyl_site_in_every_table(table: int) -> None:
    world = ex5("MAAAOAAAA", "TAG")
    sites, _ = _stops(world[1], world[0], table)
    ann = _annotate([hit("MAAAWAAAA", "c", 3, 30)], world, table, sites)
    assert ann["flags"].iloc[0] == "pyl_site"


# ---------------------------------------------------------------- 5.2: the N-terminal walk


def _plus_world(x: str) -> tuple[dict[str, str], Any]:
    """``n`` linear 90 nt: hit ``[0, 30) +``, in-phase codon ``x`` at ``[45, 48)``, CDS G
    ``+`` ``[60, 90)``."""
    nt = "GCT" * 15 + x + "GCT" * 4 + "ATG" + "GCT" * 8 + "TAA"
    assert len(nt) == 90 and nt[45:48] == x
    reps = {"n": Replicon("n", 90, "linear")}
    return {"n": nt}, model(reps, [gene("G", "n", 60, 90, "+")])


def _minus_world(x: str) -> tuple[dict[str, str], Any]:
    """``n`` linear 90 nt: CDS G ``-`` ``[0, 30)``, plus-string codon ``x`` at ``[45, 48)``,
    hit ``[60, 90) -`` (upstream of G's 5' end at 30)."""
    nt = rc("ATG" + "GCT" * 8 + "TAA") + "GCT" * 5 + x + "GCT" * 4 + "GCT" * 10
    assert len(nt) == 90 and nt[45:48] == x
    reps = {"n": Replicon("n", 90, "linear")}
    return {"n": nt}, model(reps, [gene("G", "n", 0, 30, "-")])


PLUS_HIT = Interval("n", 0, 30, "+")
MINUS_HIT = Interval("n", 60, 90, "-")


def _classify(world: Any, iv: Interval, table: int) -> tuple[str, str | None]:
    genome, gm = world
    return loci.LocusClassifier(gm, genome, table=table).classify(iv)


@pytest.mark.parametrize(
    ("x", "table", "want"),
    [
        ("TGA", 11, ("intergenic", None)),
        ("TGA", 4, ("nterm_extension", "G")),
        ("TGA", 25, ("nterm_extension", "G")),
        ("TAA", 4, ("intergenic", None)),
        ("TAG", 25, ("intergenic", None)),
        ("GCT", 11, ("nterm_extension", "G")),
        ("GCT", 4, ("nterm_extension", "G")),
    ],
)
def test_nterm_walk_plus_strand_stops_only_at_the_tables_stops(
    x: str, table: int, want: tuple[str, str | None]
) -> None:
    assert _classify(_plus_world(x), PLUS_HIT, table) == want


@pytest.mark.parametrize(
    ("x", "table", "want"),
    [
        ("TCA", 11, ("intergenic", None)),
        ("TCA", 4, ("nterm_extension", "G")),
        ("TCA", 25, ("nterm_extension", "G")),
        ("TTA", 4, ("intergenic", None)),
        ("CTA", 25, ("intergenic", None)),
        ("GCT", 4, ("nterm_extension", "G")),
    ],
)
def test_nterm_walk_minus_strand_uses_minus_stop_codons(
    x: str, table: int, want: tuple[str, str | None]
) -> None:
    """Edge case "minus strand": TCA on the plus string is a stop only under 11."""
    assert _classify(_minus_world(x), MINUS_HIT, table) == want


@pytest.mark.parametrize("table", TABLES)
def test_annotate_peptide_hits_passes_the_table_to_the_walk(table: int) -> None:
    genome, gm = _plus_world("TGA")
    df = pd.DataFrame([hit("AAAAAAAAAA", "n", 0, 30, route="exact")])
    ann = loci.annotate_peptide_hits(df, gm, genome, _peptides(["AAAAAAAAAA"]), "same", table=table)
    nterm = table != 11
    assert (ann["locus_class"].iloc[0] == "nterm_extension") is nterm
    assert ("nterm_extension_evidence" in str(ann["flags"].iloc[0])) is nterm


def test_classifier_default_table_is_11() -> None:
    genome, gm = _plus_world("TGA")
    assert loci.LocusClassifier(gm, genome).classify(PLUS_HIT) == ("intergenic", None)
    assert loci.LocusClassifier(gm, genome, table=11).classify(PLUS_HIT) == ("intergenic", None)
