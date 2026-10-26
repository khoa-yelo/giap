"""NCBI Datasets v2 taxonomy and genome reports; exact-taxon complete genome (campaign SPEC §4).

Pure parsing of the JSON reports; fetching is the caller's (``http.get_json``).
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import Any

__all__ = [
    "NCBI_API",
    "Assembly",
    "Taxon",
    "is_below_species",
    "is_prokaryote",
    "parse_assembly",
    "parse_taxon",
    "pick_assembly",
]

NCBI_API = "https://api.ncbi.nlm.nih.gov/datasets/v2"
COMPLETE_GENOME = "Complete Genome"
_PROKARYOTE_ROOTS = frozenset({2, 2157})  # Bacteria, Archaea
_TOP_CATEGORIES = frozenset({"reference genome", "representative genome"})


@dataclass(frozen=True, slots=True)
class Taxon:
    tax_id: int
    name: str
    rank: str
    lineage: tuple[int, ...]
    species_id: int | None
    domain: str | None


@dataclass(frozen=True, slots=True)
class Assembly:
    accession: str
    tax_id: int
    organism: str
    strain: str
    level: str
    refseq_category: str
    release_date: str


def _dict(value: Any) -> dict[str, Any]:
    return value if isinstance(value, dict) else {}


def _str(value: Any) -> str:
    return value if isinstance(value, str) else ""


def _int(value: Any) -> int | None:
    if isinstance(value, bool):
        return None
    if isinstance(value, int):
        return value
    if isinstance(value, str) and value.strip().isdigit():
        return int(value)
    return None


def parse_taxon(report: dict[str, Any]) -> Taxon:
    """First node of a ``/taxonomy/taxon/{taxid}`` report; no node or ``errors`` → LookupError."""
    nodes = report.get("taxonomy_nodes")
    if not isinstance(nodes, list) or not nodes or not isinstance(nodes[0], dict):
        raise LookupError("taxonomy report has no node")
    node = nodes[0]
    if node.get("errors"):
        raise LookupError(f"taxonomy report error: {node['errors']}")
    taxonomy = _dict(node.get("taxonomy"))
    tax_id = _int(taxonomy.get("tax_id"))
    if tax_id is None:
        raise LookupError("taxonomy report has no tax_id")
    lineage = taxonomy.get("lineage")
    classification = _dict(taxonomy.get("classification"))
    domain = _str(_dict(classification.get("domain")).get("name"))
    return Taxon(
        tax_id=tax_id,
        name=_str(taxonomy.get("organism_name")),
        rank=_str(taxonomy.get("rank")).upper(),
        lineage=tuple(
            i for i in map(_int, lineage if isinstance(lineage, list) else []) if i is not None
        ),
        species_id=_int(_dict(classification.get("species")).get("id")),
        domain=domain or None,
    )


def is_prokaryote(t: Taxon) -> bool:
    """Bacteria (2) or Archaea (2157) is the taxon itself or in its lineage."""
    return t.tax_id in _PROKARYOTE_ROOTS or not _PROKARYOTE_ROOTS.isdisjoint(t.lineage)


def is_below_species(t: Taxon) -> bool:
    return t.species_id is not None and t.species_id != t.tax_id


def parse_assembly(report: dict[str, Any]) -> Assembly:
    """One ``reports[]`` item of a ``/genome/taxon/{taxid}/dataset_report`` response."""
    organism = _dict(report.get("organism"))
    info = _dict(report.get("assembly_info"))
    return Assembly(
        accession=_str(report.get("accession")),
        tax_id=_int(organism.get("tax_id")) or 0,
        organism=_str(organism.get("organism_name")),
        strain=_str(_dict(organism.get("infraspecific_names")).get("strain")),
        level=_str(info.get("assembly_level")),
        refseq_category=_str(info.get("refseq_category")),
        release_date=_str(info.get("release_date")),
    )


def _rank(a: Assembly) -> tuple[bool, bool, str]:
    return (a.refseq_category in _TOP_CATEGORIES, a.accession.startswith("GCF_"), a.release_date)


def pick_assembly(assemblies: Iterable[Assembly], tax_id: int) -> tuple[Assembly | None, str]:
    """Best complete genome of exactly ``tax_id``: top category, GCF_, newest, accession."""
    candidates = [a for a in assemblies if a.level == COMPLETE_GENOME and a.tax_id == tax_id]
    if not candidates:
        return None, "no_complete_genome"
    top = max(_rank(a) for a in candidates)
    best = min((a for a in candidates if _rank(a) == top), key=lambda a: a.accession)
    return best, "unique" if len(candidates) == 1 else f"chosen_of_{len(candidates)}"
