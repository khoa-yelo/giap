"""NCBI Datasets v2 taxonomy and genome reports; exact-taxon complete-genome choice (SPEC §4)."""

from __future__ import annotations

from typing import Any

import pytest
from giap_campaign.ncbi import (
    is_below_species,
    is_prokaryote,
    parse_assembly,
    parse_taxon,
    pick_assembly,
)
from hypothesis import given, settings
from hypothesis import strategies as st

ECOLI_LINEAGE = [1, 131567, 2, 3379134, 1224, 1236, 91347, 543, 561, 562]


def _taxon_report(
    tax_id: int,
    name: str,
    lineage: list[int],
    *,
    rank: str | None = None,
    species: tuple[str, int] | None = None,
    domain: tuple[str, int] | None = None,
) -> dict[str, Any]:
    classification: dict[str, Any] = {}
    if domain is not None:
        classification["domain"] = {"name": domain[0], "id": domain[1]}
    if species is not None:
        classification["species"] = {"name": species[0], "id": species[1]}
    taxonomy: dict[str, Any] = {
        "tax_id": tax_id,
        "organism_name": name,
        "lineage": lineage,
        "classification": classification,
    }
    if rank is not None:
        taxonomy["rank"] = rank
    return {"taxonomy_nodes": [{"taxonomy": taxonomy}]}


ECOLI_MG1655 = _taxon_report(
    511145,
    "Escherichia coli str. K-12 substr. MG1655",
    [*ECOLI_LINEAGE, 83333],
    rank="STRAIN",
    species=("Escherichia coli", 562),
    domain=("Bacteria", 2),
)
ECOLI = _taxon_report(
    562,
    "Escherichia coli",
    ECOLI_LINEAGE[:-1],
    rank="SPECIES",
    species=("Escherichia coli", 562),
    domain=("Bacteria", 2),
)
HUMAN = _taxon_report(
    9606,
    "Homo sapiens",
    [1, 131567, 2759, 33154, 33208, 6072, 33213, 33511, 7711, 89593, 7742, 9604, 9605],
    rank="SPECIES",
    species=("Homo sapiens", 9606),
    domain=("Eukaryota", 2759),
)
HVO = _taxon_report(
    309800,
    "Haloferax volcanii DS2",
    [1, 131567, 2157, 28890, 183963, 2235, 1644055, 2251, 2246],
    rank="STRAIN",
    species=("Haloferax volcanii", 2246),
    domain=("Archaea", 2157),
)


# --- parse_taxon --------------------------------------------------------------------------------


def test_parse_taxon_strain() -> None:
    t = parse_taxon(ECOLI_MG1655)
    assert t.tax_id == 511145
    assert t.name == "Escherichia coli str. K-12 substr. MG1655"
    assert t.rank == "STRAIN"
    assert t.lineage == (*ECOLI_LINEAGE, 83333)
    assert isinstance(t.lineage, tuple)
    assert t.species_id == 562
    assert t.domain == "Bacteria"


def test_parse_taxon_is_frozen() -> None:
    t = parse_taxon(ECOLI)
    with pytest.raises(AttributeError):
        t.tax_id = 1  # type: ignore[misc]


def test_parse_taxon_rank_upper_cased_or_empty() -> None:
    lower = _taxon_report(562, "Escherichia coli", ECOLI_LINEAGE[:-1], rank="species")
    assert parse_taxon(lower).rank == "SPECIES"
    absent = _taxon_report(562, "Escherichia coli", ECOLI_LINEAGE[:-1])
    assert parse_taxon(absent).rank == ""


def test_parse_taxon_without_classification() -> None:
    report = _taxon_report(2, "Bacteria", [1, 131567], rank="DOMAIN")
    report["taxonomy_nodes"][0]["taxonomy"].pop("classification")
    t = parse_taxon(report)
    assert t.species_id is None
    assert t.domain is None
    assert t.tax_id == 2


@pytest.mark.parametrize(
    "report",
    [
        {},
        {"taxonomy_nodes": []},
        {"taxonomy_nodes": [{}]},
        {
            "taxonomy_nodes": [
                {
                    "query": ["999999999"],
                    "errors": [{"reason": "Taxonomy ID not found", "message": "Invalid taxid"}],
                }
            ]
        },
        {
            "taxonomy_nodes": [
                {
                    "taxonomy": {"tax_id": 562, "organism_name": "Escherichia coli"},
                    "errors": [{"reason": "something went wrong"}],
                }
            ]
        },
    ],
    ids=["empty", "no-nodes", "empty-node", "errors-only", "errors-with-taxonomy"],
)
def test_parse_taxon_missing_or_error_raises_lookup_error(report: dict[str, Any]) -> None:
    with pytest.raises(LookupError):
        parse_taxon(report)


# --- is_prokaryote / is_below_species -----------------------------------------------------------


@pytest.mark.parametrize(
    ("report", "expected"),
    [
        (ECOLI_MG1655, True),
        (ECOLI, True),
        (HVO, True),
        (HUMAN, False),
        (_taxon_report(2, "Bacteria", [1, 131567], rank="DOMAIN"), True),
        (_taxon_report(2157, "Archaea", [1, 131567], rank="DOMAIN"), True),
        (_taxon_report(10239, "Viruses", [1], rank="SUPERKINGDOM"), False),
        (_taxon_report(12345, "unclassified", []), False),
    ],
    ids=["strain", "species", "archaeon", "human", "bacteria", "archaea", "viruses", "empty"],
)
def test_is_prokaryote(report: dict[str, Any], expected: bool) -> None:
    assert is_prokaryote(parse_taxon(report)) is expected


def test_is_prokaryote_uses_lineage_not_domain_label() -> None:
    # domain label disagrees with lineage: the lineage decides
    report = _taxon_report(1280, "Staphylococcus aureus", [1, 131567, 2, 1239], domain=("x", 0))
    assert is_prokaryote(parse_taxon(report)) is True


@pytest.mark.parametrize(
    ("report", "expected"),
    [
        (ECOLI_MG1655, True),
        (HVO, True),
        (ECOLI, False),
        (_taxon_report(83333, "Escherichia coli K-12", ECOLI_LINEAGE, rank="STRAIN"), False),
    ],
    ids=["strain", "archaeal-strain", "species", "no-species-entry"],
)
def test_is_below_species(report: dict[str, Any], expected: bool) -> None:
    assert is_below_species(parse_taxon(report)) is expected


# --- parse_assembly -----------------------------------------------------------------------------


def _assembly_report(
    accession: str,
    tax_id: int = 511145,
    *,
    level: str = "Complete Genome",
    category: str | None = None,
    date: str = "2020-01-01",
    strain: str | None = "K-12 substr. MG1655",
    organism: str = "Escherichia coli str. K-12 substr. MG1655",
) -> dict[str, Any]:
    org: dict[str, Any] = {"tax_id": tax_id, "organism_name": organism}
    if strain is not None:
        org["infraspecific_names"] = {"strain": strain}
    info: dict[str, Any] = {"assembly_level": level, "release_date": date}
    if category is not None:
        info["refseq_category"] = category
    return {"accession": accession, "organism": org, "assembly_info": info}


def test_parse_assembly() -> None:
    a = parse_assembly(
        _assembly_report("GCF_000005845.2", category="reference genome", date="2013-09-26")
    )
    assert a.accession == "GCF_000005845.2"
    assert a.tax_id == 511145
    assert a.organism == "Escherichia coli str. K-12 substr. MG1655"
    assert a.strain == "K-12 substr. MG1655"
    assert a.level == "Complete Genome"
    assert a.refseq_category == "reference genome"
    assert a.release_date == "2013-09-26"
    with pytest.raises(AttributeError):
        a.accession = "x"  # type: ignore[misc]


def test_parse_assembly_missing_optional_fields() -> None:
    a = parse_assembly(_assembly_report("GCA_000000001.1", strain=None, category=None))
    assert not a.strain
    assert not a.refseq_category
    assert a.accession == "GCA_000000001.1"


# --- pick_assembly ------------------------------------------------------------------------------


def _asm(accession: str, **kw: Any) -> Any:
    return parse_assembly(_assembly_report(accession, **kw))


def test_pick_assembly_single_candidate_is_unique() -> None:
    best, reason = pick_assembly([_asm("GCF_000005845.2")], 511145)
    assert best is not None
    assert best.accession == "GCF_000005845.2"
    assert reason == "unique"


def test_pick_assembly_none() -> None:
    assert pick_assembly([], 511145) == (None, "no_complete_genome")


def test_pick_assembly_filters_level_and_tax_id() -> None:
    assemblies = [
        _asm("GCF_000000001.1", tax_id=562),
        _asm("GCF_000000002.1", level="Chromosome"),
        _asm("GCF_000000003.1", level="Scaffold"),
        _asm("GCF_000000004.1", level="Contig"),
    ]
    assert pick_assembly(assemblies, 511145) == (None, "no_complete_genome")


def test_pick_assembly_n_counts_only_candidates() -> None:
    assemblies = [
        _asm("GCA_000000010.1"),
        _asm("GCA_000000011.1"),
        _asm("GCF_000000012.1", tax_id=562),
        _asm("GCF_000000013.1", level="Chromosome"),
    ]
    best, reason = pick_assembly(assemblies, 511145)
    assert reason == "chosen_of_2"
    assert best is not None and best.accession == "GCA_000000010.1"


@pytest.mark.parametrize(
    ("assemblies", "winner"),
    [
        # reference/representative category beats GCF_ prefix and a newer date
        (
            [
                ("GCF_000000002.1", None, "2024-01-01"),
                ("GCA_000000001.1", "reference genome", "2001-01-01"),
            ],
            "GCA_000000001.1",
        ),
        (
            [
                ("GCF_000000002.1", "na", "2024-01-01"),
                ("GCA_000000001.1", "representative genome", "2001-01-01"),
            ],
            "GCA_000000001.1",
        ),
        # reference and representative are one tier: then GCF_ before GCA_
        (
            [
                ("GCA_000000001.1", "reference genome", "2024-01-01"),
                ("GCF_000000002.1", "representative genome", "2001-01-01"),
            ],
            "GCF_000000002.1",
        ),
        # GCF_ beats GCA_ even when older
        (
            [("GCA_000000001.1", None, "2024-01-01"), ("GCF_000000002.1", None, "2001-01-01")],
            "GCF_000000002.1",
        ),
        # same tier and prefix: latest release date
        (
            [
                ("GCF_000000001.1", None, "2019-05-01"),
                ("GCF_000000003.1", None, "2021-12-31"),
                ("GCF_000000002.1", None, "2020-07-15"),
            ],
            "GCF_000000003.1",
        ),
        # all equal: accession ascending
        (
            [
                ("GCF_000000009.1", None, "2020-01-01"),
                ("GCF_000000002.2", None, "2020-01-01"),
                ("GCF_000000005.1", None, "2020-01-01"),
            ],
            "GCF_000000002.2",
        ),
    ],
    ids=[
        "reference-first",
        "representative-first",
        "category-tie-gcf",
        "gcf-over-gca",
        "newest",
        "accession",
    ],
)
def test_pick_assembly_ranking(assemblies: list[tuple[str, str | None, str]], winner: str) -> None:
    objs = [_asm(acc, category=cat, date=date) for acc, cat, date in assemblies]
    best, reason = pick_assembly(objs, 511145)
    assert best is not None
    assert best.accession == winner
    assert reason == f"chosen_of_{len(assemblies)}"


_ASM = st.tuples(
    st.sampled_from(["GCF_", "GCA_"]),
    st.integers(min_value=1, max_value=50),
    st.sampled_from([None, "reference genome", "representative genome", "na"]),
    st.sampled_from(["2015-01-01", "2020-06-30", "2020-07-01", "2024-02-29"]),
    st.sampled_from(["Complete Genome", "Chromosome"]),
    st.sampled_from([511145, 562]),
)


@settings(max_examples=100, deadline=None)
@given(specs=st.lists(_ASM, max_size=8, unique_by=lambda s: (s[0], s[1])), data=st.data())
def test_pick_assembly_property_order_independent_and_dominant(
    specs: list[tuple[str, int, str | None, str, str, int]], data: st.DataObject
) -> None:
    objs = [
        _asm(f"{pre}{num:09d}.1", category=cat, date=date, level=level, tax_id=tax)
        for pre, num, cat, date, level, tax in specs
    ]
    shuffled = data.draw(st.permutations(objs))
    best, reason = pick_assembly(objs, 511145)
    assert pick_assembly(shuffled, 511145) == (best, reason)

    candidates = [a for a in objs if a.level == "Complete Genome" and a.tax_id == 511145]
    if not candidates:
        assert (best, reason) == (None, "no_complete_genome")
        return
    assert reason == ("unique" if len(candidates) == 1 else f"chosen_of_{len(candidates)}")

    def key(a: Any) -> tuple[int, int, str]:
        top = a.refseq_category in {"reference genome", "representative genome"}
        return (int(top), int(a.accession.startswith("GCF_")), a.release_date)

    assert best in candidates
    assert key(best) == max(key(a) for a in candidates)
    ties = sorted(a.accession for a in candidates if key(a) == key(best))
    assert best.accession == ties[0]
