"""PRIDE Archive REST v3 records: projects, acquisition guess and file records (campaign SPEC §3).

Pure parsing of the JSON records; fetching is the caller's (``http.get_json``).
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any, Literal

__all__ = [
    "PRIDE_API",
    "Acquisition",
    "PrideFile",
    "Project",
    "acquisition",
    "parse_file",
    "parse_project",
]

PRIDE_API = "https://www.ebi.ac.uk/pride/ws/archive/v3"
_PRIDE_FTP = "ftp://ftp.pride.ebi.ac.uk/"
_PRIDE_HTTPS = "https://ftp.pride.ebi.ac.uk/"

Acquisition = Literal["DIA", "DDA", "unknown"]

_TRAILING_INT = re.compile(r"(\d+)$")


def _words(terms: Iterable[str]) -> re.Pattern[str]:
    alternatives = "|".join(re.escape(t) for t in terms)
    return re.compile(rf"\b(?:{alternatives})\b", re.IGNORECASE)


_DIA = _words(["DIA", "SWATH", "SWATH-MS", "data-independent", "diaPASEF"])
_DDA = _words(["DDA", "data-dependent", "ddaPASEF", "shotgun"])


@dataclass(frozen=True, slots=True)
class Project:
    accession: str
    title: str
    submission_type: str
    taxids: tuple[int, ...]
    organism_names: tuple[str, ...]
    instruments: tuple[str, ...]
    experiment_types: tuple[str, ...]
    software: tuple[str, ...]
    keywords: tuple[str, ...]
    pubmed_ids: tuple[str, ...]
    dois: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PrideFile:
    name: str
    category: str
    size_bytes: int | None
    url: str | None
    checksum: str | None


def _str(value: Any) -> str:
    return value if isinstance(value, str) else ""


def _list(value: Any) -> list[Any]:
    return value if isinstance(value, list) else []


def _names(items: Any) -> tuple[str, ...]:
    """``name`` of dict items (CvParams) or the strings themselves; empty ones dropped."""
    out = [_str(i.get("name")) if isinstance(i, dict) else _str(i) for i in _list(items)]
    return tuple(n for n in out if n)


def _taxid(item: Any) -> int | None:
    """Trailing integer of ``"562"``, ``"NEWT:562"``, ``"NCBITaxon:562"`` or a dict's accession."""
    value = item.get("accession") if isinstance(item, dict) else item
    if isinstance(value, int) and not isinstance(value, bool):
        return value
    match = _TRAILING_INT.search(value) if isinstance(value, str) else None
    return int(match.group(1)) if match else None


def _pubmed_id(value: Any) -> str | None:
    if value is None or isinstance(value, bool):
        return None
    text = str(value).strip()
    return None if text in ("", "0") else text


def parse_project(record: dict[str, Any]) -> Project:
    """A PRIDE v3 project record (``/projects`` item or ``/projects/{acc}``)."""
    organisms = _list(record.get("organisms"))
    taxids = (_taxid(o) for o in organisms)
    references = [r for r in _list(record.get("references")) if isinstance(r, dict)]
    pubmed = (_pubmed_id(r.get("pubmedID")) for r in references)
    return Project(
        accession=_str(record.get("accession")),
        title=_str(record.get("title")),
        submission_type=_str(record.get("submissionType")),
        taxids=tuple(dict.fromkeys(t for t in taxids if t is not None)),
        organism_names=_names([o for o in organisms if isinstance(o, dict)]),
        instruments=_names(record.get("instruments")),
        experiment_types=_names(record.get("experimentTypes")),
        software=_names(record.get("softwares")),
        keywords=_names(record.get("keywords")),
        pubmed_ids=tuple(p for p in pubmed if p is not None),
        dois=tuple(d for d in (_str(r.get("doi")) for r in references) if d),
    )


def acquisition(project: Project) -> Acquisition:
    """Metadata guess: whole-word, case-insensitive search of experiment types, keywords, title."""
    text = "\n".join([*project.experiment_types, *project.keywords, project.title])
    if _DIA.search(text):
        return "DIA"
    if _DDA.search(text):
        return "DDA"
    return "unknown"


def _https_url(record: dict[str, Any]) -> str | None:
    for location in _list(record.get("publicFileLocations")):
        if not isinstance(location, dict) or location.get("name") != "FTP Protocol":
            continue
        value = _str(location.get("value"))
        if value.startswith(_PRIDE_FTP):
            return _PRIDE_HTTPS + value[len(_PRIDE_FTP) :]
        if value.startswith("https://"):
            return value
    return None


def parse_file(record: dict[str, Any]) -> PrideFile:
    """A ``/projects/{acc}/files/all`` record; the PRIDE FTP location rewritten to https."""
    category = record.get("fileCategory")
    size = record.get("fileSizeBytes")
    checksum = record.get("checksum")
    return PrideFile(
        name=_str(record.get("fileName")),
        category=_str(category.get("value")).upper() if isinstance(category, dict) else "",
        size_bytes=size if isinstance(size, int) and not isinstance(size, bool) else None,
        url=_https_url(record),
        checksum=checksum if isinstance(checksum, str) and checksum else None,
    )
