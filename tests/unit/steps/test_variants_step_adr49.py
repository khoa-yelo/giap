"""ADR-0049 through the wired chain (``giap.steps.variants.run_variants``): the A6a left 2->1
block flag, swap flags, and flagged calls keeping their tier (local/handoff/adr0049_decision.md,
decisions 1-2 and the implementation rules; GIAP_SPEC.md §7.5.2 step 4, §7.5.3).

Decision 1: a mass-equivalent call that meets the call threshold (``variants.high_min_peptides``
distinct peptides, ``giap.variants.pileup.sites``) stays ``high`` and carries ``mass_equivalent``.
Implementation rule A6a: a substitution is flagged when an ALT peptide also aligns as the 2->1
block ``ref[p-1:p+1] -> alt`` with its prefix shifted by one (``q = 0`` or
``P[:q] == ref[p-1-q:p-1]``); no right-hand ``ref[p:p+2]`` form, no looser prefix.

Hand-built genome: linear ``c1`` with CDS ``G1`` + [30, ...) made of the segments below (fixed
codons), and linear ``c2`` with a minus-strand CDS ``G2`` holding the N-terminal segment again.
Every site has two distinct ALT peptides (``high``) unless a test drops one; residues of a
segment are named by their index in it.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest

from giap.config import GiapConfig
from giap.model import Replicon
from giap.model.coords import Interval
from giap.model.genes import Gene, GeneModel
from giap.steps.variants import run_variants

ME = "mass_equivalent"
CODON = {
    "A": "GCT", "C": "TGT", "D": "GAT", "E": "GAA", "F": "TTT", "G": "GGT", "H": "CAT",
    "I": "ATT", "K": "AAA", "L": "CTG", "M": "ATG", "N": "AAT", "P": "CCT", "Q": "CAA",
    "R": "CGT", "S": "TCT", "T": "ACT", "V": "GTT", "W": "TGG", "Y": "TAT",
}  # fmt: skip
_COMP = str.maketrans("ACGT", "TGCA")
SPACER = "SCS"

# name -> (segment, site index in it, ((peptide, first segment index), ...), alt_type, alt_aa)
SITES: dict[str, tuple[str, int, tuple[tuple[str, int], ...], str, str]] = {
    # N-terminal AG -> Q at ALT-peptide index 0
    "nterm": ("DWAGSVATEPIKR", 3, (("QSVATEPIK", 3), ("QSVATEPIKR", 3)), "substitution", "Q"),
    # internal: AQSVATAEPIK over ...AAG..., index 1, prefix A == ref[p-2]
    "internal": ("DKAAGSVATAEPIKR", 4, (("AQSVATAEPIK", 3), ("AQSVATAEPIKR", 3)),
                 "substitution", "Q"),
    # looser prefix: AG -> Q before the site, but ref[p-2] = K is not the peptide's A
    "loose": ("DKAGTWNHYFR", 3, (("AQTWNHYF", 2), ("AQTWNHYFR", 2)), "substitution", "Q"),
    # right-hand GA -> Q at the C-terminus (suffix shift trivially true)
    "right_c": ("PDFHYWGAMK", 6, (("DFHYWQ", 1), ("PDFHYWQ", 0)), "substitution", "Q"),
    # right-hand GA -> Q inside the peptide, suffix AA == ref[p+1:p+3] == ref[p+2:p+4]
    "right_in": ("PRHYWGAAAMK", 5, (("RHYWQAA", 1), ("HYWQAA", 2)), "substitution", "Q"),
    # adjacent swap SV -> VS
    "swap": ("KENSVLDRW", 3, (("ENVSLDR", 1), ("KENVSLDR", 0)), "swap", "VS"),
    # plain mass-equivalent substitution N -> D (deamidation): flagged today, stays high
    "deamid": ("YTKNWHRF", 3, (("YTKDWHR", 0), ("TKDWHRF", 1)), "substitution", "D"),
}  # fmt: skip

G1_START = 30
SEGS = list(SITES)
PROT = "M" + SPACER
OFFSET: dict[str, int] = {}
for _name in SEGS:
    OFFSET[_name] = len(PROT)
    PROT += SITES[_name][0] + SPACER
SEQ1 = "A" * G1_START + "".join(CODON[a] for a in PROT) + "A" * 60
L1 = len(SEQ1)

# c2: minus-strand G2 holding the N-terminal segment; residue i at [e-3(i+1), e-3i). The c2
# hits reuse the c1 peptides (two loci per peptide string; the chain keys support by peptide).
PROT2 = "M" + SPACER + SITES["nterm"][0] + SPACER
G2_START = 45
G2_END = G2_START + 3 * len(PROT2)
_G2_PLUS = "".join(CODON[a] for a in PROT2).translate(_COMP)[::-1]
SEQ2 = "C" * G2_START + _G2_PLUS + "C" * 45
L2 = len(SEQ2)
NTERM_OFF2 = 1 + len(SPACER)


def codon1(name: str) -> int:
    """0-based lowest coordinate of the site codon of ``name`` on c1 (+)."""
    return G1_START + 3 * (OFFSET[name] + SITES[name][1])


def codon2() -> int:
    k = NTERM_OFF2 + SITES["nterm"][1]
    return G2_END - 3 * (k + 1)


def _hit(pep: str, contig: str, start: int, end: int, strand: str, dt: str) -> dict[str, Any]:
    return {
        "peptide": pep, "peptide_len": len(pep), "source": "denovo", "contig": contig,
        "start": start, "end": end, "strand": strand, "n_diff": 1, "diff_types": dt,
        "placement_q": 0.0, "uniqueness": "unique",
    }  # fmt: skip


def _gm() -> GeneModel:
    reps = {"c1": Replicon("c1", L1, "linear"), "c2": Replicon("c2", L2, "linear")}
    g1 = Gene("G1", "cds", Interval("c1", G1_START, G1_START + 3 * len(PROT), "+"), PROT, None,
              "ATG")  # fmt: skip
    g2 = Gene("G2", "cds", Interval("c2", G2_START, G2_END, "-"), PROT2, None, "ATG")
    return GeneModel(reps, (g1, g2), ())


def _cfg() -> GiapConfig:
    return GiapConfig.model_validate({"genome_relation": "related", "sample_name": "S1"})


def _hits(drop: frozenset[str] = frozenset()) -> list[dict[str, Any]]:
    out = []
    for name in SEGS:
        _, _, peps, kind, _ = SITES[name]
        for pep, first in peps:
            if pep in drop:
                continue
            lo = G1_START + 3 * (OFFSET[name] + first)
            out.append(_hit(pep, "c1", lo, lo + 3 * len(pep), "+", kind))
    for pep, first in SITES["nterm"][2]:
        if "c2:" + pep in drop:
            continue
        a, b = NTERM_OFF2 + first, NTERM_OFF2 + first + len(pep)
        out.append(_hit(pep, "c2", G2_END - 3 * b, G2_END - 3 * a, "-", "substitution"))
    return out


def _run(drop: frozenset[str] = frozenset()) -> pd.DataFrame:
    hits = _hits(drop)
    peptides = pd.DataFrame({"peptide": sorted({h["peptide"] for h in hits}),
                             "search_db": "same"})  # fmt: skip
    _, table, _ = run_variants(pd.DataFrame(hits), None, peptides, _gm(),
                               {"c1": SEQ1, "c2": SEQ2}, _cfg())  # fmt: skip
    return table


def _flags(cell: object) -> list[str]:
    if isinstance(cell, str):
        return [f for f in cell.split(";") if f]
    return [str(f) for f in cell] if isinstance(cell, (tuple, list)) else []


def _row(table: pd.DataFrame, contig: str, codon: int, alt_type: str) -> dict[str, Any]:
    sel = table[(table["contig"] == contig) & (table["codon_start"] == codon)
                & (table["alt_type"] == alt_type)]  # fmt: skip
    assert len(sel) == 1, f"expected one {alt_type} site at {contig}:{codon}, got\n{table}"
    return sel.iloc[0].to_dict()


def site_row(table: pd.DataFrame, name: str) -> dict[str, Any]:
    seg, k, _, kind, alt = SITES[name]
    row = _row(table, "c1", codon1(name), kind)
    ref = seg[k : k + 2] if kind == "swap" else seg[k]
    assert row["ref_aa"] == ref[0] and row["alt_aa"] == alt, row
    return row


@pytest.fixture(scope="module")
def table() -> pd.DataFrame:
    return _run()


def test_fixture_builds_the_decision_examples() -> None:
    for name, (seg, k, peps, kind, alt) in SITES.items():
        for pep, first in peps:
            q = k - first
            want = seg[first : first + len(pep)]
            if kind == "swap":
                want = want[:q] + alt + want[q + 2 :]
            else:
                want = want[:q] + alt + want[q + 1 :]
            assert pep == want, (name, pep, want)
    assert SITES["nterm"][0][2:4] == "AG" and SITES["internal"][0][2:5] == "AAG"
    assert SITES["loose"][0][1:3] == "KA" and SITES["right_c"][0][6:8] == "GA"


def test_every_site_is_one_high_alt_call(table: pd.DataFrame) -> None:
    for name in SITES:
        row = site_row(table, name)
        assert row["call"] == "ALT" and row["n_alt"] == 2, (name, row)
        assert row["confidence"] == "high", (name, row)


# --------------------------------------------------------------------------- A6a positives


def test_n_terminal_ag_to_q_read_is_flagged_and_stays_high(table: pd.DataFrame) -> None:
    row = site_row(table, "nterm")
    assert ME in _flags(row["flags"]), row
    assert row["confidence"] == "high"


def test_n_terminal_ag_to_q_read_is_flagged_on_the_minus_strand(table: pd.DataFrame) -> None:
    """The same segment on a minus-strand gene: the residue before the site (transcription order)
    lies at the higher genome coordinate."""
    row = _row(table, "c2", codon2(), "substitution")
    assert (row["ref_aa"], row["alt_aa"], row["strand"]) == ("G", "Q", "-"), row
    assert ME in _flags(row["flags"]), row
    assert row["confidence"] == "high"


def test_internal_ag_to_q_read_with_shifted_prefix_is_flagged_and_stays_high(
    table: pd.DataFrame,
) -> None:
    row = site_row(table, "internal")
    assert ME in _flags(row["flags"]), row
    assert row["confidence"] == "high"


def test_one_block_reading_peptide_flags_the_site() -> None:
    """Any ALT peptide with the block reading flags the row (scope (d): ``nterm2_any``), and a
    single-peptide call is still flagged (the flag is independent of the support count)."""
    t = _run(drop=frozenset({"QSVATEPIKR"}))
    row = site_row(t, "nterm")
    assert row["n_alt"] == 1 and row["confidence"] == "candidate"
    assert ME in _flags(row["flags"]), row


# --------------------------------------------------------------------------- A6a negatives


@pytest.mark.parametrize("name", ["loose", "right_c", "right_in"])
def test_looser_and_right_hand_forms_are_not_flagged(table: pd.DataFrame, name: str) -> None:
    """These forms flagged true positives by chance in the statistician's data."""
    row = site_row(table, name)
    assert ME not in _flags(row["flags"]), (name, row)
    assert row["confidence"] == "high"


# --------------------------------------------------------------------------- swaps and tiers


def test_swap_call_is_flagged_and_stays_high(table: pd.DataFrame) -> None:
    row = site_row(table, "swap")
    assert ME in _flags(row["flags"]), row
    assert row["confidence"] == "high"


def test_mass_equivalent_substitution_stays_high(table: pd.DataFrame) -> None:
    """Decision 1 on a plain rule-1 flag (N -> D, deamidation)."""
    row = site_row(table, "deamid")
    assert ME in _flags(row["flags"]), row
    assert row["confidence"] == "high"


def test_flags_never_change_call_confidence_or_counts(table: pd.DataFrame) -> None:
    """Every row's tier is the §7.5.2 step 4 rule on ``n_alt`` alone, flagged or not."""
    k = _cfg().variants.high_min_peptides
    for r in table.to_dict("records"):
        assert r["confidence"] == ("high" if int(r["n_alt"]) >= k else "candidate"), r
        assert r["call"] == ("mixed" if int(r["n_ref"]) else "ALT"), r
