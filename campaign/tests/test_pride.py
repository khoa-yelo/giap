"""PRIDE Archive v3 records: projects, acquisition type and file records (SPEC §3)."""

from __future__ import annotations

from typing import Any

import pytest
from giap_campaign.pride import (
    acquisition,
    parse_file,
    parse_project,
)
from hypothesis import given, settings
from hypothesis import strategies as st

FTP = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD012345/"
HTTPS = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/PXD012345/"


def _cv(name: str, accession: str = "MS:0000000") -> dict[str, str]:
    return {"@type": "CvParam", "cvLabel": "MS", "accession": accession, "name": name}


def _project(**overrides: Any) -> dict[str, Any]:
    record: dict[str, Any] = {
        "accession": "PXD012345",
        "title": "Proteome of Escherichia coli K-12",
        "submissionType": "COMPLETE",
        "organisms": [
            {"@type": "CvParam", "cvLabel": "NEWT", "accession": "562", "name": "Escherichia coli"}
        ],
        "instruments": [_cv("Q Exactive", "MS:1001911"), _cv("Orbitrap Fusion", "MS:1002416")],
        "experimentTypes": [_cv("Bottom-up proteomics", "PRIDE:0000429")],
        "softwares": [_cv("MaxQuant", "MS:1001583")],
        "keywords": ["E. coli", "stress"],
        "references": [{"referenceLine": "A et al.", "pubmedID": 31234567, "doi": "10.1000/x1"}],
    }
    record.update(overrides)
    return record


def _file(name: str, *, url: str | None = "default", **overrides: Any) -> dict[str, Any]:
    record: dict[str, Any] = {
        "fileName": name,
        "fileCategory": {"@type": "CvParam", "cvLabel": "PRIDE", "value": "SEARCH"},
        "fileSizeBytes": 1000,
        "checksum": "0123456789abcdef0123456789abcdef01234567",
        "publicFileLocations": [
            {
                "@type": "CvParam",
                "name": "Aspera Protocol",
                "value": f"prd_ascp@fasp.ebi.ac.uk:{name}",
            },
        ],
    }
    if url == "default":
        url = FTP + name
    if url is not None:
        record["publicFileLocations"].append(
            {"@type": "CvParam", "name": "FTP Protocol", "value": url}
        )
    record.update(overrides)
    return record


# --- parse_project ------------------------------------------------------------------------------


def test_parse_project_full_record() -> None:
    p = parse_project(_project())
    assert p.accession == "PXD012345"
    assert p.title == "Proteome of Escherichia coli K-12"
    assert p.submission_type == "COMPLETE"
    assert list(p.taxids) == [562]
    assert list(p.organism_names) == ["Escherichia coli"]
    assert list(p.instruments) == ["Q Exactive", "Orbitrap Fusion"]
    assert list(p.experiment_types) == ["Bottom-up proteomics"]
    assert list(p.software) == ["MaxQuant"]
    assert list(p.keywords) == ["E. coli", "stress"]
    assert list(p.pubmed_ids) == ["31234567"]
    assert list(p.dois) == ["10.1000/x1"]


def test_parse_project_is_frozen() -> None:
    p = parse_project(_project())
    with pytest.raises(AttributeError):
        p.accession = "PXD999999"  # type: ignore[misc]


def test_parse_project_minimal_record_defaults() -> None:
    p = parse_project({"accession": "PXD000001"})
    assert p.accession == "PXD000001"
    assert p.title == ""
    assert p.submission_type == ""
    for field in (
        "taxids",
        "organism_names",
        "instruments",
        "experiment_types",
        "software",
        "keywords",
        "pubmed_ids",
        "dois",
    ):
        assert list(getattr(p, field)) == [], field


def test_parse_project_taxid_forms_dedup_and_junk() -> None:
    organisms = [
        {"accession": "562", "name": "Escherichia coli"},
        {"accession": "NEWT:562", "name": "Escherichia coli"},
        {"accession": "NCBITaxon:511145", "name": "Escherichia coli str. K-12 substr. MG1655"},
        "9606",
        "NEWT:1280",
        {"accession": "not applicable", "name": "Not available"},
        {"name": "no accession"},
        {"accession": ""},
        "",
        "unknown",
        {"accession": "NCBITaxon:9606", "name": "Homo sapiens"},
    ]
    p = parse_project(_project(organisms=organisms))
    assert list(p.taxids) == [562, 511145, 9606, 1280]
    assert all(isinstance(t, int) for t in p.taxids)


def test_parse_project_organism_names_from_dict_items() -> None:
    organisms = [
        {"accession": "NEWT:83333", "name": "Escherichia coli K-12"},
        {"accession": "9606", "name": "Homo sapiens (Human)"},
    ]
    p = parse_project(_project(organisms=organisms))
    assert list(p.organism_names) == ["Escherichia coli K-12", "Homo sapiens (Human)"]


def test_parse_project_pubmed_and_doi_filtering() -> None:
    refs = [
        {"pubmedID": 31234567, "doi": "10.1000/a"},
        {"pubmedID": "0", "doi": ""},
        {"pubmedID": 0},
        {"pubmedID": None, "doi": None},
        {"pubmedID": ""},
        {"pubmedID": "29876543", "doi": "10.1000/b"},
        {"referenceLine": "no ids"},
    ]
    p = parse_project(_project(references=refs))
    assert list(p.pubmed_ids) == ["31234567", "29876543"]
    assert all(isinstance(x, str) for x in p.pubmed_ids)
    assert list(p.dois) == ["10.1000/a", "10.1000/b"]


@settings(max_examples=60, deadline=None)
@given(
    items=st.lists(
        st.tuples(
            st.integers(min_value=1, max_value=3_000_000),
            st.sampled_from(["", "NEWT:", "NCBITaxon:"]),
            st.booleans(),
        ),
        max_size=8,
    )
)
def test_parse_project_taxids_property(items: list[tuple[int, str, bool]]) -> None:
    organisms: list[Any] = [
        {"accession": f"{prefix}{t}", "name": f"org {t}"} if as_dict else f"{prefix}{t}"
        for t, prefix, as_dict in items
    ]
    expected = list(dict.fromkeys(t for t, _, _ in items))
    assert list(parse_project(_project(organisms=organisms)).taxids) == expected


# --- acquisition --------------------------------------------------------------------------------


def _acq(title: str = "", keywords: list[str] | None = None, types: list[str] | None = None) -> str:
    record = _project(
        title=title,
        keywords=keywords or [],
        experimentTypes=[_cv(t) for t in (types or [])],
    )
    return acquisition(parse_project(record))


@pytest.mark.parametrize(
    ("title", "keywords", "types", "expected"),
    [
        ("E. coli DIA proteome", [], [], "DIA"),
        ("", ["dia"], [], "DIA"),
        ("", [], ["SWATH MS"], "DIA"),
        ("SWATH-MS of Bacillus subtilis", [], [], "DIA"),
        ("Data-independent acquisition of a biofilm", [], [], "DIA"),
        ("", ["diaPASEF"], [], "DIA"),
        ("", ["DIAPASEF", "timsTOF"], [], "DIA"),
        ("Comparison of DIA and DDA", [], [], "DIA"),
        ("", ["DDA"], ["SWATH MS"], "DIA"),
        ("E. coli DDA proteome", [], [], "DDA"),
        ("", [], ["Shotgun proteomics"], "DDA"),
        ("Data-dependent acquisition", [], [], "DDA"),
        ("", ["dda"], [], "DDA"),
        ("Proteome of E. coli", ["stress"], ["Bottom-up proteomics"], "unknown"),
        ("DIAMOND search of metaproteome", [], [], "unknown"),
        ("Diagnostic markers in media", ["DIAMOND"], [], "unknown"),
        ("", ["ADDA", "DDAX"], [], "unknown"),
        ("", [], [], "unknown"),
    ],
)
def test_acquisition(title: str, keywords: list[str], types: list[str], expected: str) -> None:
    assert _acq(title, keywords, types) == expected


# --- parse_file ---------------------------------------------------------------------------------


def test_parse_file_converts_pride_ftp_to_https() -> None:
    f = parse_file(_file("evidence.txt"))
    assert f.name == "evidence.txt"
    assert f.category == "SEARCH"
    assert f.size_bytes == 1000
    assert f.url == HTTPS + "evidence.txt"
    assert f.checksum == "0123456789abcdef0123456789abcdef01234567"


def test_parse_file_is_frozen() -> None:
    f = parse_file(_file("evidence.txt"))
    with pytest.raises(AttributeError):
        f.name = "x"  # type: ignore[misc]


def test_parse_file_category_upper_cased_and_defaults() -> None:
    record = _file("psm.tsv", fileCategory={"value": "result"})
    assert parse_file(record).category == "RESULT"
    record = _file("psm.tsv")
    del record["fileCategory"]
    del record["fileSizeBytes"]
    del record["checksum"]
    f = parse_file(record)
    assert f.category == ""
    assert f.size_bytes is None
    assert f.checksum is None


@pytest.mark.parametrize("checksum", ["", None, 12345])
def test_parse_file_checksum_only_non_empty_string(checksum: Any) -> None:
    assert parse_file(_file("psm.tsv", checksum=checksum)).checksum is None


@pytest.mark.parametrize(
    "url",
    [
        None,
        "prd_ascp@fasp.ebi.ac.uk:pride/data/archive/2021/03/PXD012345/psm.tsv",
        "s3://bucket/psm.tsv",
        "gsiftp://example.org/psm.tsv",
    ],
    ids=["no-ftp-location", "aspera-style", "s3", "gsiftp"],
)
def test_parse_file_other_location_gives_no_url(url: str | None) -> None:
    assert parse_file(_file("psm.tsv", url=url)).url is None


def test_parse_file_no_locations_at_all() -> None:
    record = _file("psm.tsv")
    del record["publicFileLocations"]
    assert parse_file(record).url is None
