"""M11-11: the preprint outline cites only evidence that exists and invents no numbers."""

from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUTLINE = ROOT / "docs" / "preprint" / "outline.md"
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")  # M10-05 guard
ESTIMATE = re.compile(r"\b(q|FDR|FWER|CI|FDP|F1|KS)\b|rejection rate|^unestimated|^n/a")
REQUIRED_SOURCES = (
    [f"§13.2 b{i}" for i in range(1, 8)]
    + [f"§15.6 b{i}" for i in range(1, 6)]
    + [f"§15.7 b{i}" for i in range(1, 5)]  # pre-registered hypotheses (review S8)
)


def _rows() -> list[dict[str, str]]:
    text = OUTLINE.read_text(encoding="utf-8")
    table = text.split("## 4. Claim → evidence map", 1)[1].split("\n## ", 1)[0]
    lines = [ln for ln in table.splitlines() if ln.startswith("| C")]
    keys = ["id", "source", "claim", "number", "estimate", "evidence", "status"]
    return [
        dict(zip(keys, (c.strip() for c in ln.strip("|").split(" | ")), strict=True))
        for ln in lines
    ]


def _paths(cell: str) -> list[str]:
    spans = re.findall(r"`([^`]+)`", cell)
    return [s.split()[0] for s in spans if "/" in s or s.endswith((".md", ".json"))]


def test_every_spec_bullet_is_a_row() -> None:
    found = set()
    for r in _rows():
        head = r["source"].split(" ", 1)[0]
        found |= {f"{head} b{b}" for b in re.findall(r"\bb(\d+)\b", r["source"])}
    assert not set(REQUIRED_SOURCES) - found, sorted(set(REQUIRED_SOURCES) - found)


def test_error_estimate_rule() -> None:  # R4
    assert not ESTIMATE.search("TBD") and not ESTIMATE.search("median of 3")  # example 3
    for r in _rows():
        assert ESTIMATE.search(r["estimate"]), (r["id"], r["estimate"])


def test_supported_rows_cite_existing_paths() -> None:
    for r in _rows():
        paths = _paths(r["evidence"])
        assert paths, r["id"]
        assert r["status"].startswith(("supported", "pending")), r["id"]
        if r["status"].startswith("supported"):
            missing = [p for p in paths if not (ROOT / p).exists()]
            assert not missing, (r["id"], missing)
        else:
            assert re.search(r"M\d+-\d+|#\d+|Khoa|spec", r["status"]), r["id"]


def _numbers(text: str) -> set[str]:
    return {n.replace(",", "") for n in re.findall(r"\d[\d,]*", text)}


def test_no_invented_numbers() -> None:  # M10-15 notes rule
    text = OUTLINE.read_text(encoding="utf-8")
    cited = {p for p in _paths(text) if (ROOT / p).is_file()}
    corpus = (ROOT / "GIPPA_SPEC.md").read_text(encoding="utf-8") + " ".join(cited)
    for p in cited:
        corpus += (ROOT / p).read_text(encoding="utf-8", errors="replace")
    known = _numbers(corpus)
    unknown = sorted(n for n in _numbers(text) if len(n) >= 4 and n not in known)
    assert not unknown, unknown


def test_no_row_level_sequences() -> None:
    m = PEPTIDE_LIKE.search(OUTLINE.read_text(encoding="utf-8"))
    assert m is None, m.group(0) if m else None


def test_uncalibrated_flags_are_named() -> None:  # review B4, B5: never imply calibration
    text = OUTLINE.read_text(encoding="utf-8")
    for flag in ("ora_uncalibrated", "null_uncalibrated", "ADR-0049", "ADR-0055"):
        assert flag in text.split("## 6. Limitations", 1)[1], flag
