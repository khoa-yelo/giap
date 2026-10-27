"""M11-02: CITATION.cff complete, schema-valid, in agreement with the version and the changelog.

Spec: docs/specs/M11-02-changelog-citation.md (Definition items 3-5, examples 3-6, acceptance
criteria 3-5). The CFF 1.2.0 JSON Schema is vendored under ``tests/data/schemas`` (see its
PROVENANCE.md); nothing here touches the network.
"""

from __future__ import annotations

import datetime as dt
import hashlib
import json
import re
from pathlib import Path
from typing import Any

import jsonschema
import pytest
import yaml

import gippa

ROOT = Path(__file__).resolve().parents[2]
CFF = ROOT / "CITATION.cff"
CHANGELOG = ROOT / "CHANGELOG.md"
SCHEMA = ROOT / "tests" / "data" / "schemas" / "cff-1.2.0.schema.json"
PROVENANCE = SCHEMA.parent / "PROVENANCE.md"
SCHEMA_SHA256 = "0b8d22140da702d766df318dcff3a91af2f39521298dcf36d76315fd99cc169b"

KEYWORDS = {"proteogenomics", "bacteria", "genome context", "operons", "target-decoy", "Nextflow"}
EMAIL_RE = re.compile(r"[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}")


class _NoDateLoader(yaml.SafeLoader):
    """SafeLoader that keeps ISO dates as strings, as the CFF schema expects."""


_NoDateLoader.yaml_implicit_resolvers = {
    k: [(tag, rx) for tag, rx in v if tag != "tag:yaml.org,2002:timestamp"]
    for k, v in yaml.SafeLoader.yaml_implicit_resolvers.items()
}


def load_cff(text: str) -> dict[str, Any]:
    data = yaml.load(text, Loader=_NoDateLoader)
    assert isinstance(data, dict)
    return data


def agreement_errors(cff: dict[str, Any], version: str, changelog: str) -> list[str]:
    """Item 4: CFF version == gippa.__version__; a date-released has a matching changelog heading."""
    errors = []
    cff_version = str(cff.get("version"))
    if cff_version != version:
        errors.append(f"CITATION.cff version {cff_version} != gippa {version}")
    if "date-released" in cff:
        date = str(cff["date-released"])
        try:
            dt.date.fromisoformat(date)
        except ValueError:
            errors.append(f"CITATION.cff date-released {date!r} is not an ISO date")
        else:
            heading = f"## GIPPA {cff_version} ({date})"
            if not any(line.strip() == heading for line in changelog.splitlines()):
                errors.append(f"CHANGELOG.md has no heading {heading!r}")
    return errors


def email_errors(cff: dict[str, Any], text: str) -> list[str]:
    """R12: no personal email address anywhere in the citation file."""
    errors = [
        f"{key} {i} has an email key"
        for key in ("authors", "contact")
        for i, a in enumerate(cff.get(key) or [])
        if isinstance(a, dict) and "email" in a
    ]
    errors += [f"email address in CITATION.cff: {m}" for m in EMAIL_RE.findall(text)]
    return errors


@pytest.fixture(scope="module")
def cff_text() -> str:
    return CFF.read_text(encoding="utf-8")


@pytest.fixture(scope="module")
def cff(cff_text: str) -> dict[str, Any]:
    return load_cff(cff_text)


@pytest.fixture(scope="module")
def schema() -> dict[str, Any]:
    return json.loads(SCHEMA.read_text(encoding="utf-8"))


# ------------------------------------------------------------------ vendored schema (item 5)


def test_vendored_schema_matches_provenance() -> None:
    assert SCHEMA.is_file(), "tests/data/schemas/cff-1.2.0.schema.json is missing"
    assert hashlib.sha256(SCHEMA.read_bytes()).hexdigest() == SCHEMA_SHA256
    prov = PROVENANCE.read_text(encoding="utf-8")
    for needle in (SCHEMA_SHA256, "1.2.0", "CC-BY-4.0", "github.com/citation-file-format"):
        assert needle in prov


def test_schema_is_the_cff_1_2_0_schema(schema: dict[str, Any]) -> None:
    assert schema["$id"] == "https://citation-file-format.github.io/1.2.0/schema.json"
    jsonschema.Draft7Validator.check_schema(schema)


def test_citation_validates_against_schema(cff: dict[str, Any], schema: dict[str, Any]) -> None:
    validator = jsonschema.Draft7Validator(schema, format_checker=jsonschema.FormatChecker())
    errors = sorted(validator.iter_errors(cff), key=lambda e: list(e.path))
    assert errors == [], [f"{list(e.path)}: {e.message}" for e in errors]


def test_schema_rejects_a_broken_citation(schema: dict[str, Any]) -> None:
    validator = jsonschema.Draft7Validator(schema, format_checker=jsonschema.FormatChecker())
    bad = {"cff-version": "1.2.0", "title": "x", "message": "m", "authors": [], "type": "paper"}
    assert list(validator.iter_errors(bad))


# ------------------------------------------------------------------ item 3: keys due now


def test_cff_version_and_type(cff: dict[str, Any]) -> None:
    assert cff.get("cff-version") == "1.2.0"
    assert cff.get("type") == "software"


def test_title_and_message(cff: dict[str, Any]) -> None:
    assert str(cff.get("title", "")).startswith("GIPPA")
    assert cff.get("message")


def test_abstract_is_present_and_at_most_120_words(cff: dict[str, Any]) -> None:
    abstract = str(cff.get("abstract", ""))
    words = abstract.split()
    assert words, "CITATION.cff has no abstract"
    assert len(words) <= 120, f"abstract has {len(words)} words"


def test_authors_have_given_and_family_names(cff: dict[str, Any]) -> None:
    authors = cff.get("authors")
    assert isinstance(authors, list) and authors
    for a in authors:
        assert a.get("given-names") and a.get("family-names"), a
    assert {"given-names": "Khoa", "family-names": "Hoang"}.items() <= authors[0].items()


def test_license_and_repository(cff: dict[str, Any]) -> None:
    assert cff.get("license") == "MIT"
    assert cff.get("repository-code") == "https://github.com/khoa-yelo/gippa"


def test_keywords(cff: dict[str, Any]) -> None:
    keywords = cff.get("keywords")
    assert isinstance(keywords, list)
    missing = KEYWORDS - set(keywords)
    assert not missing, f"missing keywords: {sorted(missing)}"


def test_optional_keys_when_present_are_well_formed(cff: dict[str, Any]) -> None:
    # date-released, identifiers, url become due at M11-13; if present they must be sane now.
    if "url" in cff:
        assert str(cff["url"]).startswith("https://")
    if "identifiers" in cff:
        assert any(i.get("type") == "doi" for i in cff["identifiers"])
    if "date-released" in cff:
        dt.date.fromisoformat(str(cff["date-released"]))


def test_no_email_in_citation(cff: dict[str, Any], cff_text: str) -> None:
    assert email_errors(cff, cff_text) == []


# ------------------------------------------------------------------ item 4: agreement


def test_citation_agrees_with_version_and_changelog(cff: dict[str, Any]) -> None:
    errors = agreement_errors(cff, gippa.__version__, CHANGELOG.read_text(encoding="utf-8"))
    assert errors == []


def test_changelog_has_no_release_section_before_release_prep(cff: dict[str, Any]) -> None:
    # Acceptance 5: CHANGELOG.md is written once, by the release prep (M11-13), which also adds
    # date-released; until then it carries no "## GIPPA" release heading.
    if "date-released" not in cff:
        text = CHANGELOG.read_text(encoding="utf-8")
        assert "<!-- towncrier release notes start -->" in text
        assert not re.search(r"^## GIPPA ", text, flags=re.MULTILINE)


# examples 3-6 on synthetic inputs (self-tests of the checks above)

BASE = """\
cff-version: 1.2.0
message: "If you use GIPPA, please cite it as below."
title: "GIPPA"
type: software
authors:
  - family-names: Hoang
    given-names: Khoa
version: 0.0.1
"""


def test_example3_matching_version_passes() -> None:
    assert agreement_errors(load_cff(BASE), "0.0.1", "# Changelog\n") == []


def test_example4_bumped_version_fails() -> None:
    errors = agreement_errors(load_cff(BASE), "1.0.0", "# Changelog\n")
    assert errors == ["CITATION.cff version 0.0.1 != gippa 1.0.0"]


def test_example5_date_without_changelog_heading_fails() -> None:
    cff = load_cff(BASE.replace("0.0.1", "1.0.0") + "date-released: 2026-11-01\n")
    assert cff["date-released"] == "2026-11-01"  # kept as a string for the schema
    errors = agreement_errors(
        cff, "1.0.0", "# Changelog\n\n<!-- towncrier release notes start -->\n"
    )
    assert len(errors) == 1
    assert "## GIPPA 1.0.0 (2026-11-01)" in errors[0]
    ok = "# Changelog\n\n## GIPPA 1.0.0 (2026-11-01)\n\n### Features\n"
    assert agreement_errors(cff, "1.0.0", ok) == []


def test_example5_invalid_date_fails() -> None:
    cff = load_cff(BASE + "date-released: 2026-13-40\n")
    assert any("not an ISO date" in e for e in agreement_errors(cff, "0.0.1", ""))


def test_example6_author_email_fails() -> None:
    text = BASE.replace(
        "given-names: Khoa\n", "given-names: Khoa\n    email: someone@example.org\n"
    )
    errors = email_errors(load_cff(text), text)
    assert any("email key" in e for e in errors)
    assert any("someone@example.org" in e for e in errors)
    assert email_errors(load_cff(BASE), BASE) == []
