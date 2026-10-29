"""The ``peptide_hits`` schema's ``flags`` pattern (M7-16b rule 18; lead decision Q16; acceptance
"the ``peptide_hits`` schema ``flags`` pattern accepts every committed ``peptide_hits.tsv`` and the
new flag, and rejects an unknown value").

The pattern has the form ``variants.schema.json`` uses, ``(?:a|b|...)(?:;(?:a|b|...))*``, over every
value the code writes today plus ``readthrough``: ``engine_protein_mismatch``, ``low_complexity``,
``multi``, ``multi_identical``, ``nterm_extension_evidence``, ``origin_spanning``, ``readthrough``,
``same_db``, ``strain_difference``. Its description names ``readthrough``. 16c adds ``sec_site``
and ``pyl_site``; this file does not require or forbid them. M11-20 (ADR-0068) adds
``initiator_m`` (spec § 5), now required.
"""

from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCHEMA = REPO / "src" / "gippa" / "schemas" / "tables" / "peptide_hits.schema.json"
VARIANTS = REPO / "src" / "gippa" / "schemas" / "tables" / "variants.schema.json"
VALUES = (
    "engine_protein_mismatch",
    "initiator_m",
    "low_complexity",
    "multi",
    "multi_identical",
    "nterm_extension_evidence",
    "origin_spanning",
    "readthrough",
    "same_db",
    "strain_difference",
)
COMMITTED = sorted(
    p
    for p in (
        REPO / "tests/data/golden/sim50k/mapping/peptide_hits.tsv",
        REPO / "tests/data/snapshot/sim_empty/OUT/mapping/peptide_hits.tsv",
        REPO / "tests/data/snapshot/sim_plots/OUT/mapping/peptide_hits.tsv",
    )
    if p.exists()
)  # workflow/tests/data/genetable/mapping/peptide_hits.tsv is a one-line stub, not a table


def _text(path: Path) -> str:
    """Committed TSVs, some stored gzip-compressed under a ``.tsv`` name (the snapshot)."""
    data = path.read_bytes()
    return (gzip.decompress(data) if data[:2] == b"\x1f\x8b" else data).decode("utf-8")


def _field(path: Path, name: str) -> dict:  # type: ignore[type-arg]
    doc = json.loads(path.read_text(encoding="utf-8"))
    (f,) = [f for f in doc["fields"] if f["name"] == name]
    return f  # type: ignore[no-any-return]


def _pattern() -> re.Pattern[str]:
    f = _field(SCHEMA, "flags")
    pat = f.get("constraints", {}).get("pattern")
    assert pat, "peptide_hits.schema.json: field 'flags' has no constraints.pattern (rule 18)"
    return re.compile(pat)


def test_pattern_has_the_variants_schema_form() -> None:
    pat = _pattern().pattern
    m = re.fullmatch(r"\(\?:([a-z_|]+)\)\(\?:;\(\?:([a-z_|]+)\)\)\*", pat)
    assert m, pat
    assert m.group(1) == m.group(2)
    listed = set(m.group(1).split("|"))
    assert set(VALUES) <= listed, sorted(set(VALUES) - listed)
    assert listed - set(VALUES) <= {"sec_site", "pyl_site"}, sorted(listed - set(VALUES))
    # the same form as the variants schema
    vpat = _field(VARIANTS, "flags")["constraints"]["pattern"]
    assert re.fullmatch(r"\(\?:[a-z_|]+\)\(\?:;\(\?:[a-z_|]+\)\)\*", vpat)


def test_description_names_readthrough() -> None:
    assert "readthrough" in _field(SCHEMA, "flags")["description"]


@pytest.mark.parametrize(
    "value",
    [
        *VALUES,
        "readthrough;strain_difference",
        "multi;multi_identical;origin_spanning",
        "low_complexity;same_db",
        "nterm_extension_evidence;readthrough;strain_difference",
    ],
)
def test_pattern_accepts_written_values(value: str) -> None:
    assert _pattern().fullmatch(value), value


@pytest.mark.parametrize(
    "value",
    [
        "bogus",
        "readthrough;bogus",
        ";readthrough",
        "readthrough;",
        "readthrough,strain_difference",
        "Readthrough",
        "stop_codon",
        "readthrough;;same_db",
    ],
)
def test_pattern_rejects_unknown_values(value: str) -> None:
    assert not _pattern().fullmatch(value), value


def test_every_committed_peptide_hits_tsv_matches_the_pattern() -> None:
    assert COMMITTED, "no committed peptide_hits.tsv found"
    pat = _pattern()
    for path in COMMITTED:
        lines = _text(path).splitlines()
        header = lines[0].split("\t")
        k = header.index("flags")
        bad = sorted(
            {
                v
                for ln in lines[1:]
                if ln
                for v in [ln.split("\t")[k]]
                if v not in ("", "NA") and not pat.fullmatch(v)
            }
        )
        assert not bad, f"{path.relative_to(REPO)}: flags {bad[:5]} do not match the pattern"


@pytest.mark.parametrize("path", COMMITTED, ids=lambda p: str(p.relative_to(REPO)))
def test_committed_peptide_hits_validate_against_the_schema(path: Path, tmp_path: Path) -> None:
    frictionless = pytest.importorskip("frictionless")
    plain = tmp_path / "peptide_hits.tsv"
    plain.write_text(_text(path), encoding="utf-8")
    descriptor = json.loads(SCHEMA.read_text(encoding="utf-8"))
    descriptor = {k: v for k, v in descriptor.items() if k != "coordinates"}
    _pattern()  # the check is about the new pattern; fail first without it
    resource = frictionless.Resource(
        path=plain.name,
        basepath=str(plain.parent),
        schema=frictionless.Schema.from_descriptor(descriptor),
    )
    report = resource.validate()
    errors = [e for t in report.tasks for e in t.errors if "flags" in str(e.message)]
    assert not errors, [e.message for e in errors[:5]]
