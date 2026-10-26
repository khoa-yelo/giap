"""Sample sheet from PRIDE projects, NCBI taxa and assemblies, PRIDE file lists (SPEC §5).

Every lookup is an injected fake: no network. ``Project``, ``PrideFile``, ``Taxon`` and ``Assembly``
objects are built through the §3/§4 parsers so these tests do not depend on dataclass field order.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable
from typing import Any
from urllib.parse import parse_qs, urlsplit

import numpy as np
import pandas as pd
import pytest
from giap_campaign.catalog import build_sheet, iter_projects
from giap_campaign.ncbi import parse_assembly, parse_taxon
from giap_campaign.pride import parse_file, parse_project
from hypothesis import given, settings
from hypothesis import strategies as st

PRIDE = "https://www.ebi.ac.uk/pride/ws/archive/v3"
FTP = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/"
HTTPS = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/03/"

SHEET_COLUMNS = [
    "dataset",
    "tax_id",
    "organism",
    "below_species",
    "title",
    "submission_type",
    "acquisition",
    "instruments",
    "pubmed_ids",
    "dois",
    "genome_accession",
    "genome_choice",
    "raw_format",
    "n_raw",
    "raw_gb",
    "raw_files",
    "raw_urls",
    "raw_bytes",
    "raw_checksums",
    "formats_seen",
    "flags",
    "status",
]

BACTERIA = [1, 131567, 2]
ECOLI_LINEAGE = [*BACTERIA, 1224, 1236, 91347, 543, 561, 562]
MD5 = "0123456789abcdef0123456789abcdef"
SHA1 = "0123456789abcdef0123456789abcdef01234567"


# --- builders -----------------------------------------------------------------------------------


def _taxon_report(
    tax_id: int,
    name: str,
    lineage: list[int],
    *,
    rank: str = "SPECIES",
    species_id: int | None = None,
) -> dict[str, Any]:
    classification: dict[str, Any] = {}
    if species_id is not None:
        classification["species"] = {"name": name, "id": species_id}
    return {
        "taxonomy_nodes": [
            {
                "taxonomy": {
                    "tax_id": tax_id,
                    "organism_name": name,
                    "rank": rank,
                    "lineage": lineage,
                    "classification": classification,
                }
            }
        ]
    }


TAXA = {
    511145: _taxon_report(
        511145,
        "Escherichia coli str. K-12 substr. MG1655",
        [*ECOLI_LINEAGE, 83333],
        rank="STRAIN",
        species_id=562,
    ),
    562: _taxon_report(562, "Escherichia coli", ECOLI_LINEAGE[:-1], species_id=562),
    9606: _taxon_report(9606, "Homo sapiens", [1, 131567, 2759, 33154], species_id=9606),
    4932: _taxon_report(4932, "Saccharomyces cerevisiae", [1, 131567, 2759, 4751], species_id=4932),
    9: _taxon_report(9, "Buchnera aphidicola", [*BACTERIA, 1224, 1236, 32199], species_id=9),
    10: _taxon_report(10, "Cellvibrio", [*BACTERIA, 1224, 1236, 72274], rank="GENUS"),
    309800: _taxon_report(
        309800,
        "Haloferax volcanii DS2",
        [1, 131567, 2157, 28890, 2246],
        rank="STRAIN",
        species_id=2246,
    ),
}
PROKARYOTES = {511145, 562, 9, 10, 309800}


def _assembly(
    accession: str,
    tax_id: int = 511145,
    *,
    level: str = "Complete Genome",
    category: str | None = None,
    date: str = "2020-01-01",
) -> Any:
    info: dict[str, Any] = {"assembly_level": level, "release_date": date}
    if category is not None:
        info["refseq_category"] = category
    return parse_assembly(
        {
            "accession": accession,
            "organism": {"tax_id": tax_id, "organism_name": "x", "infraspecific_names": {}},
            "assembly_info": info,
        }
    )


def _project(
    accession: str = "PXD000001",
    taxids: list[int] | None = None,
    **overrides: Any,
) -> Any:
    record: dict[str, Any] = {
        "accession": accession,
        "title": "DIA proteome of Escherichia coli K-12",
        "submissionType": "COMPLETE",
        "organisms": [
            {"accession": f"NEWT:{t}", "name": f"taxon {t}"}
            for t in (taxids if taxids is not None else [511145])
        ],
        "instruments": [{"name": "Q Exactive"}, {"name": "Orbitrap Fusion"}],
        "experimentTypes": [{"name": "Bottom-up proteomics"}],
        "keywords": ["E. coli"],
        "references": [
            {"pubmedID": 31234567, "doi": "10.1000/x1"},
            {"pubmedID": "29876543", "doi": "10.1000/x2"},
        ],
    }
    record.update(overrides)
    return parse_project(record)


def _file(
    name: str,
    size: int | None = 1000,
    checksum: str | None = None,
    *,
    url: bool = True,
) -> Any:
    record: dict[str, Any] = {
        "fileName": name,
        "fileCategory": {"value": "RAW"},
        "publicFileLocations": [{"name": "Aspera Protocol", "value": f"prd_ascp@x:{name}"}],
    }
    if size is not None:
        record["fileSizeBytes"] = size
    if checksum is not None:
        record["checksum"] = checksum
    if url:
        record["publicFileLocations"].append({"name": "FTP Protocol", "value": FTP + name})
    return parse_file(record)


class Fakes:
    """Injected ``taxon``, ``assemblies`` and ``files`` with call recording."""

    def __init__(
        self,
        *,
        assemblies: dict[int, list[Any] | BaseException] | None = None,
        files: dict[str, list[Any] | BaseException] | None = None,
        taxon_errors: dict[int, BaseException] | None = None,
    ) -> None:
        self.assembly_map = assemblies or {}
        self.file_map = files or {}
        self.taxon_errors = taxon_errors or {}
        self.taxon_calls: Counter[int] = Counter()
        self.assembly_calls: Counter[int] = Counter()
        self.file_calls: Counter[str] = Counter()

    def taxon(self, taxid: int) -> Any:
        self.taxon_calls[taxid] += 1
        if taxid in self.taxon_errors:
            raise self.taxon_errors[taxid]
        if taxid not in TAXA:
            raise LookupError(f"taxid {taxid} not found")
        return parse_taxon(TAXA[taxid])

    def assemblies(self, taxid: int) -> list[Any]:
        self.assembly_calls[taxid] += 1
        value = self.assembly_map.get(taxid, [])
        if isinstance(value, BaseException):
            raise value
        return list(value)

    def files(self, accession: str) -> list[Any]:
        self.file_calls[accession] += 1
        value = self.file_map.get(accession, [_file("run1.raw")])
        if isinstance(value, BaseException):
            raise value
        return list(value)

    def build(self, projects: list[Any]) -> pd.DataFrame:
        return build_sheet(projects, taxon=self.taxon, assemblies=self.assemblies, files=self.files)


def _complete(taxid: int = 511145) -> dict[int, list[Any] | BaseException]:
    return {taxid: [_assembly("GCF_000005845.2", taxid, category="reference genome")]}


def _only_row(df: pd.DataFrame) -> pd.Series:
    assert len(df) == 1, df
    return df.iloc[0]


# --- rows per project and taxid -----------------------------------------------------------------


def test_empty_input_gives_empty_sheet_with_all_columns() -> None:
    df = Fakes().build([])
    assert len(df) == 0
    assert list(df.columns) == SHEET_COLUMNS


@pytest.mark.parametrize("taxids", [[9606], [9606, 4932], []], ids=["human", "eukaryotes", "none"])
def test_project_without_prokaryote_gives_no_row_and_no_file_listing(taxids: list[int]) -> None:
    fakes = Fakes()
    df = fakes.build([_project("PXD000001", taxids)])
    assert len(df) == 0
    assert list(df.columns) == SHEET_COLUMNS
    assert fakes.file_calls["PXD000001"] == 0


def test_multi_organism_project_one_row_for_the_prokaryote() -> None:
    fakes = Fakes(assemblies=_complete())
    row = _only_row(fakes.build([_project("PXD000001", [9606, 511145])]))
    assert int(row["tax_id"]) == 511145
    assert row["flags"] == "multi_organism"
    assert row["status"] == "ready"
    assert fakes.file_calls["PXD000001"] == 1


def test_single_organism_project_has_no_flags() -> None:
    row = _only_row(Fakes(assemblies=_complete()).build([_project("PXD000001", [511145])]))
    assert row["flags"] == ""


def test_partial_submission_flag_and_sorted_flags() -> None:
    fakes = Fakes(assemblies=_complete())
    row = _only_row(fakes.build([_project("PXD000001", [511145], submissionType="PARTIAL")]))
    assert row["flags"] == "partial_submission"
    assert row["submission_type"] == "PARTIAL"
    both = _only_row(fakes.build([_project("PXD000002", [511145, 9606], submissionType="PARTIAL")]))
    assert both["flags"] == "multi_organism;partial_submission"


def test_two_prokaryotes_two_rows_but_files_listed_once() -> None:
    fakes = Fakes(assemblies=_complete())
    df = fakes.build([_project("PXD000001", [511145, 309800, 9606])])
    assert [int(t) for t in df["tax_id"]] == [309800, 511145]
    assert set(df["flags"]) == {"multi_organism"}
    assert fakes.file_calls == Counter({"PXD000001": 1})


def test_files_called_once_per_project() -> None:
    fakes = Fakes(assemblies=_complete())
    fakes.build(
        [
            _project("PXD000001", [511145]),
            _project("PXD000002", [511145, 562]),
            _project("PXD000003", [9606]),
        ]
    )
    assert fakes.file_calls == Counter({"PXD000001": 1, "PXD000002": 1})


def test_taxon_memoised_across_projects() -> None:
    fakes = Fakes(assemblies=_complete())
    projects = [_project(f"PXD00000{i}", [511145, 9606]) for i in range(1, 5)]
    df = fakes.build(projects)
    assert len(df) == 4
    assert fakes.taxon_calls[511145] == 1
    assert fakes.taxon_calls[9606] == 1


def test_rows_sorted_by_dataset_then_numeric_tax_id() -> None:
    fakes = Fakes()
    df = fakes.build(
        [
            _project("PXD000002", [562]),
            _project("PXD000001", [10, 9, 562]),
        ]
    )
    assert [(d, int(t)) for d, t in zip(df["dataset"], df["tax_id"], strict=True)] == [
        ("PXD000001", 9),
        ("PXD000001", 10),
        ("PXD000001", 562),
        ("PXD000002", 562),
    ]


# --- genome choice and status -------------------------------------------------------------------


def test_species_level_taxid_needs_genome() -> None:
    fakes = Fakes(assemblies=_complete(562))
    row = _only_row(fakes.build([_project("PXD000001", [562])]))
    assert row["below_species"] == "false"
    assert row["organism"] == "Escherichia coli"
    assert row["genome_choice"] == "species_level"
    assert row["genome_accession"] == ""
    assert row["status"] == "needs_genome"
    # the raw files are still described
    assert row["raw_format"] == "thermo_raw"


@pytest.mark.parametrize(
    "assemblies",
    [
        [],
        [_assembly("GCF_000000001.1", level="Chromosome")],
        [_assembly("GCF_000000002.1", tax_id=562)],
    ],
    ids=["none", "chromosome-only", "other-taxon"],
)
def test_strain_without_complete_genome_needs_genome(assemblies: list[Any]) -> None:
    fakes = Fakes(assemblies={511145: assemblies})
    row = _only_row(fakes.build([_project("PXD000001", [511145])]))
    assert row["below_species"] == "true"
    assert row["genome_choice"] == "no_complete_genome"
    assert row["genome_accession"] == ""
    assert row["status"] == "needs_genome"


@pytest.mark.parametrize(
    "error", [LookupError("no reports"), OSError("HTTP 503")], ids=["lookup", "oserror"]
)
def test_assembly_lookup_failure_needs_genome(error: BaseException) -> None:
    fakes = Fakes(assemblies={511145: error})
    row = _only_row(fakes.build([_project("PXD000001", [511145])]))
    assert row["genome_accession"] == ""
    assert row["status"] == "needs_genome"
    assert int(row["tax_id"]) == 511145


@pytest.mark.parametrize(
    "error", [LookupError("taxid not found"), OSError("HTTP 503")], ids=["lookup", "oserror"]
)
def test_taxon_failure_gives_taxon_error_row(error: BaseException) -> None:
    """Assumed: the row keeps the project's taxid and project metadata; organism and the genome
    columns are empty (the taxon is unknown)."""
    fakes = Fakes(assemblies=_complete(), taxon_errors={424242: error})
    row = _only_row(fakes.build([_project("PXD000001", [424242])]))
    assert row["dataset"] == "PXD000001"
    assert int(row["tax_id"]) == 424242
    assert row["organism"] == ""
    assert row["genome_accession"] == ""
    assert row["title"] == "DIA proteome of Escherichia coli K-12"
    assert row["status"] == "taxon_error"


def test_taxon_error_beside_a_good_taxid() -> None:
    fakes = Fakes(assemblies=_complete(), taxon_errors={424242: LookupError("nope")})
    df = fakes.build([_project("PXD000001", [511145, 424242])])
    by_tax = {int(t): s for t, s in zip(df["tax_id"], df["status"], strict=True)}
    assert by_tax == {511145: "ready", 424242: "taxon_error"}
    assert set(df["flags"]) == {"multi_organism"}
    assert fakes.file_calls["PXD000001"] == 1


def test_only_wiff_files_gives_no_raw_files() -> None:
    files = [_file("x.wiff"), _file("x.wiff.scan"), _file("y.wiff2"), _file("evidence.txt")]
    fakes = Fakes(assemblies=_complete(), files={"PXD000001": files})
    row = _only_row(fakes.build([_project("PXD000001", [511145])]))
    assert row["raw_format"] == ""
    assert int(row["n_raw"]) == 0
    assert row["raw_files"] == ""
    assert "sciex_wiff" in row["formats_seen"].split(";")
    assert row["formats_seen"] == "sciex_wiff"
    assert row["status"] == "no_raw_files"


def test_files_lookup_failure_gives_no_raw_files() -> None:
    """Assumed: an OSError/LookupError from ``files`` leaves the row without raw files."""
    fakes = Fakes(assemblies=_complete(), files={"PXD000001": OSError("HTTP 500")})
    row = _only_row(fakes.build([_project("PXD000001", [511145])]))
    assert row["status"] == "no_raw_files"
    assert int(row["n_raw"]) == 0


@pytest.mark.parametrize(
    ("taxid", "taxon_error", "files", "expected"),
    [
        (424242, True, [_file("x.wiff")], "taxon_error"),
        (424242, True, [_file("a.raw")], "taxon_error"),
        (562, False, [_file("x.wiff")], "no_raw_files"),
        (511145, False, [_file("x.wiff")], "no_raw_files"),
        (562, False, [_file("a.raw")], "needs_genome"),
        (511145, False, [_file("a.raw")], "ready"),
    ],
    ids=[
        "taxon-error-over-no-raw",
        "taxon-error-alone",
        "no-raw-over-species-level",
        "no-raw-over-ready",
        "needs-genome",
        "ready",
    ],
)
def test_status_precedence(taxid: int, taxon_error: bool, files: list[Any], expected: str) -> None:
    errors: dict[int, BaseException] = {taxid: LookupError("x")} if taxon_error else {}
    assemblies = {**_complete(511145), **_complete(562)}
    fakes = Fakes(assemblies=assemblies, files={"PXD000001": files}, taxon_errors=errors)
    row = _only_row(fakes.build([_project("PXD000001", [taxid])]))
    assert row["status"] == expected


@pytest.mark.parametrize(
    ("names", "fmt"),
    [
        (["b.d.zip", "a.d.tar.gz"], "bruker_d"),
        (["c.mzML.gz", "d.mzml"], "mzml"),
        (["a.raw", "x.wiff", "y.wiff", "z.wiff"], "thermo_raw"),
    ],
)
def test_raw_format_from_choose_raw(names: list[str], fmt: str) -> None:
    fakes = Fakes(assemblies=_complete(), files={"PXD000001": [_file(n) for n in names]})
    row = _only_row(fakes.build([_project("PXD000001", [511145])]))
    assert row["raw_format"] == fmt
    assert row["status"] == "ready"


# --- the ready row, column by column ------------------------------------------------------------


def test_ready_row_every_column() -> None:
    files = [
        _file("c_run3.raw", None, SHA1),
        _file("b_run2.raw", 2_500_000_000, MD5),
        _file("a_run1.raw", 1_234_567_890, None),
        _file("z_nourl.raw", 5, None, url=False),
        _file("b.d.zip", 10),
        _file("c.mzML.gz", 10),
        _file("x.wiff", 10),
        _file("x.wiff.scan", 10),
        _file("results.txt", 10),
    ]
    assemblies = {
        511145: [
            _assembly("GCA_000801205.1", date="2024-01-01"),
            _assembly("GCF_000005845.2", category="reference genome", date="2013-09-26"),
            _assembly("GCF_000000099.1", level="Contig"),
        ]
    }
    fakes = Fakes(assemblies=assemblies, files={"PXD012345": files})
    df = fakes.build([_project("PXD012345", [511145])])
    assert list(df.columns) == SHEET_COLUMNS
    row = _only_row(df)

    assert row["dataset"] == "PXD012345"
    assert isinstance(row["tax_id"], int | np.integer)
    assert row["tax_id"] == 511145
    assert row["organism"] == "Escherichia coli str. K-12 substr. MG1655"
    assert row["below_species"] == "true"
    assert row["title"] == "DIA proteome of Escherichia coli K-12"
    assert row["submission_type"] == "COMPLETE"
    assert row["acquisition"] == "DIA"
    assert row["instruments"] == "Q Exactive;Orbitrap Fusion"
    assert row["pubmed_ids"] == "31234567;29876543"
    assert row["dois"] == "10.1000/x1;10.1000/x2"
    assert row["genome_accession"] == "GCF_000005845.2"
    assert row["genome_choice"] == "chosen_of_2"
    assert row["raw_format"] == "thermo_raw"
    assert int(row["n_raw"]) == 3
    assert float(row["raw_gb"]) == 3.73  # (1_234_567_890 + 2_500_000_000) / 1e9, 2 decimals
    assert row["raw_files"] == "a_run1.raw;b_run2.raw;c_run3.raw"
    assert row["raw_urls"] == ";".join(
        HTTPS + n for n in ("a_run1.raw", "b_run2.raw", "c_run3.raw")
    )
    assert row["raw_bytes"] == "1234567890;2500000000;"
    assert row["raw_checksums"] == f";{MD5};{SHA1}"
    assert row["formats_seen"] == "bruker_d;mzml;sciex_wiff;thermo_raw"
    assert row["flags"] == ""
    assert row["status"] == "ready"


def test_ready_row_unique_genome_and_list_columns_empty() -> None:
    fakes = Fakes(
        assemblies=_complete(),
        files={"PXD000001": [_file("only.raw", 999_999)]},
    )
    project = _project("PXD000001", [511145], instruments=[], references=[], title="")
    row = _only_row(fakes.build([project]))
    assert row["genome_choice"] == "unique"
    assert row["instruments"] == ""
    assert row["pubmed_ids"] == ""
    assert row["dois"] == ""
    assert row["title"] == ""
    assert row["acquisition"] == "unknown"
    assert float(row["raw_gb"]) == 0.0  # 999_999 bytes rounds to 0.00 GB
    assert row["raw_bytes"] == "999999"
    assert row["raw_checksums"] == ""


def test_no_nan_or_none_cells() -> None:
    fakes = Fakes(
        assemblies={511145: LookupError("x")},
        files={"PXD000001": [_file("a.raw", None, None)], "PXD000002": [_file("x.wiff")]},
        taxon_errors={424242: LookupError("x")},
    )
    df = fakes.build(
        [_project("PXD000001", [511145, 562]), _project("PXD000002", [424242, 309800])]
    )
    assert len(df) == 4
    for column in SHEET_COLUMNS:
        for value in df[column]:
            assert value is not None, column
            assert not (isinstance(value, float) and np.isnan(value)), column
            assert str(value) not in {"nan", "NaN", "NA", "None", "<NA>"}, column
    known = df[df["status"] != "taxon_error"]
    assert set(known["below_species"]) == {"true", "false"}


_POOL = [511145, 562, 9, 10, 309800, 9606, 4932]


@settings(max_examples=60, deadline=None)
@given(
    projects=st.lists(
        st.tuples(
            st.integers(min_value=1, max_value=60),
            st.lists(st.sampled_from(_POOL), max_size=4, unique=True),
        ),
        max_size=6,
        unique_by=lambda p: p[0],
    )
)
def test_property_rows_order_and_calls(projects: list[tuple[int, list[int]]]) -> None:
    fakes = Fakes(assemblies={t: [_assembly(f"GCF_{t:09d}.1", t)] for t in _POOL})
    objs = [_project(f"PXD{num:06d}", taxids) for num, taxids in projects]
    df = fakes.build(objs)

    expected = sorted(
        (f"PXD{num:06d}", t) for num, taxids in projects for t in taxids if t in PROKARYOTES
    )
    got = [(d, int(t)) for d, t in zip(df["dataset"], df["tax_id"], strict=True)]
    assert got == expected
    with_rows = {f"PXD{num:06d}" for num, taxids in projects if set(taxids) & PROKARYOTES}
    assert fakes.file_calls == Counter(dict.fromkeys(with_rows, 1))
    assert all(n == 1 for n in fakes.taxon_calls.values())
    flags = dict(zip(df["dataset"], df["flags"], strict=True))
    for num, taxids in projects:
        acc = f"PXD{num:06d}"
        if acc in flags:
            assert flags[acc] == ("multi_organism" if len(taxids) > 1 else "")


# --- iter_projects ------------------------------------------------------------------------------


class FakeApi:
    """``get_json`` fake serving ``/projects`` pages and ``/projects/{acc}`` records."""

    def __init__(
        self,
        records: list[dict[str, Any]],
        wrap: Callable[[list[dict[str, Any]]], Any] = lambda page: page,
    ) -> None:
        self.records = records
        self.wrap = wrap
        self.urls: list[str] = []

    def __call__(self, url: str, **kwargs: Any) -> Any:
        self.urls.append(url)
        parts = urlsplit(url)
        assert parts.scheme == "https"
        path = f"{parts.scheme}://{parts.netloc}{parts.path}"
        if path == f"{PRIDE}/projects":
            query = parse_qs(parts.query)
            size = int(query["pageSize"][0])
            page = int(query["page"][0])
            return self.wrap(self.records[page * size : (page + 1) * size])
        prefix = f"{PRIDE}/projects/"
        assert path.startswith(prefix), url
        acc = path[len(prefix) :]
        for r in self.records:
            if r["accession"] == acc:
                return r
        raise LookupError(acc)

    def page_params(self) -> list[tuple[int, int]]:
        out = []
        for url in self.urls:
            query = parse_qs(urlsplit(url).query)
            out.append((int(query["pageSize"][0]), int(query["page"][0])))
        return out


def _records(n: int) -> list[dict[str, Any]]:
    return [{"accession": f"PXD{i:06d}", "title": f"project {i}"} for i in range(1, n + 1)]


def test_iter_projects_pages_until_empty_page() -> None:
    api = FakeApi(_records(7))
    got = list(iter_projects(api, page_size=3))
    assert [r["accession"] for r in got] == [f"PXD{i:06d}" for i in range(1, 8)]
    assert all(urlsplit(u).path.endswith("/pride/ws/archive/v3/projects") for u in api.urls)
    assert all(u.startswith(f"{PRIDE}/projects?") for u in api.urls)
    # pages 0, 1, 2 hold 3 + 3 + 1 records; page 3 is empty and ends paging
    assert api.page_params() == [(3, 0), (3, 1), (3, 2), (3, 3)]


def test_iter_projects_exact_multiple_and_empty_archive() -> None:
    api = FakeApi(_records(4))
    assert len(list(iter_projects(api, page_size=2))) == 4
    assert api.page_params() == [(2, 0), (2, 1), (2, 2)]
    empty = FakeApi([])
    assert list(iter_projects(empty, page_size=5)) == []
    assert empty.page_params() == [(5, 0)]


def test_iter_projects_default_page_size() -> None:
    api = FakeApi(_records(3))
    assert len(list(iter_projects(api))) == 3
    assert api.page_params()[0] == (100, 0)


def test_iter_projects_max_projects_truncates_mid_page() -> None:
    api = FakeApi(_records(10))
    got = list(iter_projects(api, page_size=3, max_projects=4))
    assert [r["accession"] for r in got] == [f"PXD{i:06d}" for i in range(1, 5)]
    assert api.page_params() == [(3, 0), (3, 1)]


def test_iter_projects_max_projects_at_page_boundary_stops_fetching() -> None:
    api = FakeApi(_records(10))
    got = list(iter_projects(api, page_size=3, max_projects=3))
    assert len(got) == 3
    assert api.page_params() == [(3, 0)]


def test_iter_projects_is_lazy() -> None:
    api = FakeApi(_records(10))
    it = iter_projects(api, page_size=3)
    assert api.urls == []
    first = next(it)
    assert first["accession"] == "PXD000001"
    assert len(api.urls) == 1


@pytest.mark.parametrize(
    "wrap",
    [
        lambda page: {"_embedded": {"projects": page}, "page": {"size": 3}},
        lambda page: {"list": page},
        lambda page: {"content": page, "totalElements": 99},
    ],
    ids=["_embedded.projects", "list", "content"],
)
def test_iter_projects_accepts_dict_wrapped_pages(
    wrap: Callable[[list[dict[str, Any]]], Any],
) -> None:
    api = FakeApi(_records(5), wrap=wrap)
    got = list(iter_projects(api, page_size=3))
    assert [r["accession"] for r in got] == [f"PXD{i:06d}" for i in range(1, 6)]
    assert api.page_params() == [(3, 0), (3, 1), (3, 2)]


def test_iter_projects_accessions_in_given_order() -> None:
    api = FakeApi(_records(5))
    got = list(iter_projects(api, accessions=["PXD000004", "PXD000002", "PXD000005"]))
    assert [r["accession"] for r in got] == ["PXD000004", "PXD000002", "PXD000005"]
    assert api.urls == [
        f"{PRIDE}/projects/PXD000004",
        f"{PRIDE}/projects/PXD000002",
        f"{PRIDE}/projects/PXD000005",
    ]
