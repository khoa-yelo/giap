"""M7-14 rule 3: the ``id_mixups`` QC count (docs/specs/M7-14-te2-proxy-and-id-mixups.md examples
3-4, every ``id_mixups`` degenerate row, lead decision 1; GIAP_SPEC.md §7.5.1 row 1, §9.3).

§7.5.1: "Search against this strain's Bakta proteins (``search_db = same``) | cannot occur, except
at Sec sites and non-ATG starts, which appear as ALT calls | consistency check: ALT sites without
the ``start_codon`` or ``stop_codon`` flag are counted as ``id_mixups``".

Interface pinned here (test-author reading of the spec's signature
``id_mixups(sites, support) -> tuple[int, int]``)::

    giap.variants.qc.id_mixups(sites: pd.DataFrame, support: pd.DataFrame) -> tuple[int, int]

* ``sites``: the ``variants.tsv`` frame (M7-03 ``vcf_records`` table / ``read_variants``), one row
  per ALT site (``call`` ALT or mixed), keyed by ``variant_id``; ``flags`` is a tuple (or list) of
  §9.4 flag names. Only ``variant_id`` and ``flags`` decide the count; confidence never does
  (lead decision 1: any confidence).
* ``support``: one row per (site, peptide supporting that site's ALT): columns ``variant_id``,
  ``peptide``, ``source`` (``database`` | ``denovo``) and ``search_db`` (``same`` | ``other``,
  already resolved). A peptide listed twice for one site counts once; rows whose ``variant_id``
  is not in ``sites`` are ignored.
* Returns ``(count, n)`` as Python ``int``: ``n`` = ALT sites with at least one ``database``,
  ``search_db = same`` supporting peptide; ``count`` = those of them flagged neither
  ``start_codon`` nor ``stop_codon``.
"""

from __future__ import annotations

from typing import Any

import pandas as pd
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig

SUPPORT_COLUMNS = ["variant_id", "peptide", "source", "search_db"]
FLAGS = ("mass_equivalent", "snv", "multi_nt", "strain_difference", "database_derived",
         "start_codon", "stop_codon", "vcf_unrepresentable")  # fmt: skip


def _id_mixups(sites: pd.DataFrame, support: pd.DataFrame) -> tuple[int, int]:
    from giap.variants.qc import id_mixups

    return id_mixups(sites, support)


def _site(
    vid: str,
    *,
    ref: str = "N",
    alt: str = "D",
    flags: tuple[str, ...] = ("snv",),
    confidence: str = "high",
    call: str = "ALT",
    codon_start: int = 300,
    strand: str = "+",
    ref_codon: str = "AAC",
    aa_pos: int | None = 21,
) -> dict[str, Any]:
    return {
        "variant_id": vid, "contig": "chr1", "codon_start": codon_start, "strand": strand,
        "gene": "G1", "aa_pos": aa_pos, "ref_codon": ref_codon, "ref_aa": ref, "alt_aa": alt,
        "alt_type": "substitution", "call": call, "confidence": confidence, "n_ref": 0,
        "n_alt": 2, "source": "database", "flags": tuple(flags),
    }  # fmt: skip


def _sites(*rows: dict[str, Any]) -> pd.DataFrame:
    cols = list(_site("x"))
    return pd.DataFrame(list(rows), columns=cols)


def _support(*rows: tuple[str, str, str, str]) -> pd.DataFrame:
    return pd.DataFrame(list(rows), columns=SUPPORT_COLUMNS)


# ------------------------------------------------------------------ spec examples 3-4


def test_example_3_same_db_internal_n_to_d_is_one_mixup() -> None:
    sites = _sites(_site("VAR_1", ref="N", alt="D", flags=("snv",)))
    support = _support(("VAR_1", "PEPTIDEDK", "database", "same"),
                       ("VAR_1", "PEPTIDEDKR", "database", "same"))  # fmt: skip
    assert _id_mixups(sites, support) == (1, 1)


def test_example_4_same_alt_supported_only_by_de_novo_peptides_is_not_counted() -> None:
    sites = _sites(_site("VAR_1", ref="N", alt="D", flags=("snv",)))
    support = _support(("VAR_1", "PEPTIDEDK", "denovo", "same"),
                       ("VAR_1", "PEPTIDEDKR", "denovo", "same"))  # fmt: skip
    assert _id_mixups(sites, support) == (0, 0)


def test_return_values_are_python_ints() -> None:
    sites = _sites(_site("VAR_1"))
    count, n = _id_mixups(sites, _support(("VAR_1", "PEPTIDEDK", "database", "same")))
    assert type(count) is int and type(n) is int


# ------------------------------------------------------------------ degenerate rows


def test_no_same_db_peptides_gives_zero_with_n_zero() -> None:
    """Spec: "no ``search_db = same`` peptides | 0 with n 0 (not NA: the check ran)"."""
    sites = _sites(_site("VAR_1"), _site("VAR_2", codon_start=600))
    support = _support(("VAR_1", "PEPTIDEDK", "database", "other"),
                       ("VAR_2", "PEPTIDEAK", "denovo", "other"))  # fmt: skip
    assert _id_mixups(sites, support) == (0, 0)


def test_empty_inputs_give_zero_with_n_zero() -> None:
    assert _id_mixups(_sites(), _support()) == (0, 0)


def test_sites_without_any_support_rows_give_zero_with_n_zero() -> None:
    assert _id_mixups(_sites(_site("VAR_1")), _support()) == (0, 0)


def _flags(site: dict[str, Any], starts: set[tuple[str, str, int]], cfg: GiapConfig,
           sources: set[tuple[str, str]]) -> tuple[str, ...]:  # fmt: skip
    from giap.variants.flags import flags_for

    return tuple(flags_for(site, starts, cfg, sources))


def _cfg(relation: str = "same") -> GiapConfig:
    return GiapConfig.model_validate({"seed": 1, "genome_relation": relation})


def test_same_db_alt_m_at_a_gtg_start_is_not_a_mixup() -> None:
    """Spec: "a same-db ALT at a GTG start with ALT M | not a mix-up (flag ``start_codon``)",
    with the flag set by the M7-05 rule itself."""
    starts = {("chr1", "+", 300)}
    raw = _site("VAR_1", ref="V", alt="M", ref_codon="GTG", aa_pos=1, flags=())
    raw["flags"] = _flags(raw, starts, _cfg(), {("database", "same")})
    assert "start_codon" in raw["flags"]
    support = _support(("VAR_1", "MKLLVAGR", "database", "same"))
    assert _id_mixups(_sites(raw), support) == (0, 1)


def test_same_db_alt_other_than_m_at_a_start_codon_is_a_mixup() -> None:
    """Spec: "a same-db ALT with ALT != M at a start codon | a mix-up"."""
    starts = {("chr1", "+", 300)}
    raw = _site("VAR_1", ref="V", alt="L", ref_codon="GTG", aa_pos=1, flags=())
    raw["flags"] = _flags(raw, starts, _cfg(), {("database", "same")})
    assert "start_codon" not in raw["flags"]
    support = _support(("VAR_1", "LKLLVAGR", "database", "same"))
    assert _id_mixups(_sites(raw), support) == (1, 1)


def test_genome_relation_related_is_still_counted() -> None:
    """Spec: "``genome_relation = related`` | ``id_mixups`` still counted": the
    ``strain_difference`` flag the relation adds does not exempt a site."""
    raw = _site("VAR_1", flags=())
    raw["flags"] = _flags(raw, set(), _cfg("related"), {("database", "same")})
    assert "strain_difference" in raw["flags"]
    support = _support(("VAR_1", "PEPTIDEDK", "database", "same"))
    assert _id_mixups(_sites(raw), support) == (1, 1)


def test_stop_codon_flag_exempts_a_site() -> None:
    flags = ("multi_nt", "stop_codon")
    sites = _sites(_site("VAR_1", ref="*", alt="U", ref_codon="TGA", flags=flags))
    support = _support(("VAR_1", "PEPTUDEK", "database", "same"))
    assert _id_mixups(sites, support) == (0, 1)


@pytest.mark.parametrize("confidence", ["high", "candidate"])
def test_decision_1_any_confidence_is_counted(confidence: str) -> None:
    """Lead decision 1: contradicting sites of any confidence; TE2 is the ``high``-only test."""
    sites = _sites(_site("VAR_1", confidence=confidence))
    support = _support(("VAR_1", "PEPTIDEDK", "database", "same"))
    assert _id_mixups(sites, support) == (1, 1)


@pytest.mark.parametrize("call", ["ALT", "mixed"])
def test_mixed_codons_are_alt_sites_too(call: str) -> None:
    """An ALT site is any ``variants.tsv`` row (``call`` ALT or mixed, as in M9-08 rule 10)."""
    sites = _sites(_site("VAR_1", call=call))
    support = _support(("VAR_1", "PEPTIDEDK", "database", "same"))
    assert _id_mixups(sites, support) == (1, 1)


def test_one_same_db_peptide_among_others_makes_the_site_same_db() -> None:
    sites = _sites(_site("VAR_1"))
    support = _support(("VAR_1", "AAAAPEPK", "denovo", "same"),
                       ("VAR_1", "BBBBPEPK", "database", "other"),
                       ("VAR_1", "PEPTIDEDK", "database", "same"))  # fmt: skip
    assert _id_mixups(sites, support) == (1, 1)


def test_database_other_only_is_not_a_same_db_site() -> None:
    sites = _sites(_site("VAR_1", flags=("snv", "database_derived")))
    support = _support(("VAR_1", "PEPTIDEDK", "database", "other"))
    assert _id_mixups(sites, support) == (0, 0)


def test_a_site_counts_once_however_many_same_db_peptides_support_it() -> None:
    sites = _sites(_site("VAR_1"), _site("VAR_2", codon_start=900, flags=("start_codon",)))
    support = _support(*[("VAR_1", f"PEP{i}DK", "database", "same") for i in range(5)],
                       ("VAR_1", "PEP0DK", "database", "same"),
                       *[("VAR_2", f"MPEP{i}K", "database", "same") for i in range(3)])  # fmt: skip
    assert _id_mixups(sites, support) == (1, 2)


def test_support_rows_of_unknown_sites_are_ignored() -> None:
    sites = _sites(_site("VAR_1", flags=("start_codon",)))
    support = _support(("VAR_1", "MPEPK", "database", "same"),
                       ("VAR_99", "PEPTIDEDK", "database", "same"))  # fmt: skip
    assert _id_mixups(sites, support) == (0, 1)


def test_hand_built_mixture() -> None:
    """Five sites: two mix-ups (one candidate), one start, one stop, one de novo only."""
    sites = _sites(
        _site("VAR_1"),
        _site("VAR_2", codon_start=600, confidence="candidate", flags=("multi_nt",)),
        _site("VAR_3", codon_start=900, ref="V", alt="M", flags=("snv", "start_codon")),
        _site("VAR_4", codon_start=1200, ref="*", alt="U", flags=("stop_codon",)),
        _site("VAR_5", codon_start=1500),
    )
    support = _support(("VAR_1", "P1K", "database", "same"), ("VAR_2", "P2K", "database", "same"),
                       ("VAR_3", "P3K", "database", "same"), ("VAR_4", "P4K", "database", "same"),
                       ("VAR_5", "P5K", "denovo", "same"))  # fmt: skip
    assert _id_mixups(sites, support) == (2, 4)


# ------------------------------------------------------------------ properties


@st.composite
def _case(draw: st.DrawFn) -> tuple[pd.DataFrame, pd.DataFrame]:
    n_sites = draw(st.integers(0, 12))
    rows = []
    for i in range(n_sites):
        flags = draw(st.lists(st.sampled_from(FLAGS), unique=True, max_size=4))
        conf = draw(st.sampled_from(["high", "candidate"]))
        rows.append(_site(f"VAR_{i + 1}", codon_start=3 * i, confidence=conf, flags=tuple(flags)))
    sup = draw(st.lists(
        st.tuples(st.sampled_from([f"VAR_{i}" for i in range(1, n_sites + 2)]),
                  st.sampled_from(["PA", "PB", "PC", "PD"]),
                  st.sampled_from(["database", "denovo"]),
                  st.sampled_from(["same", "other"])),
        max_size=30,
    ))  # fmt: skip
    return _sites(*rows), _support(*sup)


def _oracle(sites: pd.DataFrame, support: pd.DataFrame) -> tuple[int, int]:
    same = set(support.loc[(support["source"] == "database") & (support["search_db"] == "same"),
                           "variant_id"])  # fmt: skip
    n = count = 0
    for vid, flags in zip(sites["variant_id"], sites["flags"], strict=True):
        if vid in same:
            n += 1
            count += not ({"start_codon", "stop_codon"} & set(flags))
    return count, n


@settings(max_examples=150, deadline=None)
@given(_case())
def test_property_equals_the_rule_3_oracle(case: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    sites, support = case
    got = _id_mixups(sites, support)
    assert got == _oracle(sites, support)
    count, n = got
    assert 0 <= count <= n <= len(sites)


@settings(max_examples=100, deadline=None)
@given(_case(), st.randoms(use_true_random=False))
def test_property_row_order_does_not_matter(
    case: tuple[pd.DataFrame, pd.DataFrame], rnd: Any
) -> None:
    sites, support = case
    want = _id_mixups(sites, support)
    s_idx, p_idx = list(range(len(sites))), list(range(len(support)))
    rnd.shuffle(s_idx)
    rnd.shuffle(p_idx)
    got = _id_mixups(sites.iloc[s_idx].reset_index(drop=True),
                     support.iloc[p_idx].reset_index(drop=True))  # fmt: skip
    assert got == want


@settings(max_examples=100, deadline=None)
@given(_case())
def test_property_de_novo_and_other_db_support_never_changes_the_result(
    case: tuple[pd.DataFrame, pd.DataFrame],
) -> None:
    sites, support = case
    want = _id_mixups(sites, support)
    extra = _support(*[(v, "XTRA", src, db) for v in sites["variant_id"]
                       for src, db in (("denovo", "same"), ("denovo", "other"),
                                       ("database", "other"))])  # fmt: skip
    assert _id_mixups(sites, pd.concat([support, extra], ignore_index=True)) == want


@settings(max_examples=100, deadline=None)
@given(_case(), st.sampled_from(["start_codon", "stop_codon"]))
def test_property_an_exemption_flag_never_raises_the_count(
    case: tuple[pd.DataFrame, pd.DataFrame], flag: str
) -> None:
    sites, support = case
    count, n = _id_mixups(sites, support)
    flagged = sites.copy()
    flagged["flags"] = [(*f, flag) if flag not in f else tuple(f) for f in sites["flags"]]
    assert _id_mixups(flagged, support) == (0, n)
    assert count <= n


@settings(max_examples=100, deadline=None)
@given(_case())
def test_property_confidence_is_ignored(case: tuple[pd.DataFrame, pd.DataFrame]) -> None:
    sites, support = case
    want = _id_mixups(sites, support)
    for conf in ("high", "candidate"):
        alt = sites.copy()
        alt["confidence"] = conf
        assert _id_mixups(alt, support) == want
