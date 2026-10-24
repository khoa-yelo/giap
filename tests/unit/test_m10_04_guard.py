"""Guard: M10-04 commits no run output and no biological number (acceptance criterion 7).

Spec: docs/specs/M10-04-pilot-runs.md, "What is committed (and what is not)"; Khoa's rule of
2026-10-23 (nothing derived from the DIA-NN reports is committed); Q-K1 parked with its default
(lead decision 2): the README section holds provenance only. Extends M10-03's guard.

The README rule is an allow-list, as the spec asks. A number in the section "Pass-1 runs
(M10-04)" is allowed only as (test-author reading):

* part of an identifier or path (a token with a letter in it: ``out_a1``, ``CP171242.1``,
  ``M10-04``, ``wheat-07``, ``.r2``), a ``§`` reference, a version (``1.12.2``, ``DB 6.0``);
* a hexadecimal digest or git sha (7-64 hex characters), a Slurm job ID (6-9 digits);
* a byte size or memory (``2,145,332 bytes``, ``54 GB``, ``2104832K``, ``54.GB``);
* a duration (``05:59:12``, ``2-00:00:00``, ``5 h 59 min``, ``12 s``), a date or timestamp;
* an exit code (``0:0``, ``exit 1``, ``exit code 4``, ``exit_code 0``);
* a CPU count (``14 CPUs``, ``--threads 14``, ``AllocCPUS 14``), an attempt (``attempt 2``);
* the seed ``20261014``; a lone digit (an ordinal, an attempt or an exit code in a table cell);
* a Markdown list ordinal at the start of a line.

Fenced code blocks (the commands) are not read. Anything else with a digit fails: counts of
genes, percentages, rho, Jaccards.
"""

from __future__ import annotations

import fnmatch
import re
import subprocess
from pathlib import PurePosixPath

from helpers.slurm_sandbox import REPO

README = REPO / "validation" / "pilot" / "README.md"
SECTION = "Pass-1 runs (M10-04)"
FORBIDDEN = ("out_a*", "cmp_*", "ATTEMPTS.tsv", "report.html", "compare.html")
DELIVERABLES = ("validation/pilot/run_evidence.py", "scripts/slurm/pilot_compare.sbatch")

_ALLOWED = [
    r"```.*?```",  # fenced code blocks
    r"\b\d{4}-\d{2}-\d{2}(?:[T ][0-9:.]+(?:Z|[+-]\d{2}:?\d{2})?)?",  # dates, timestamps
    r"\b(?:\d+-)?\d{1,3}:\d{2}:\d{2}(?:\.\d+)?\b",  # Slurm durations
    r"(?<![\w.])\d[\d,]*(?:\.\d+)?\s*(?:bytes?|[KMGT]i?B|[KMGT])\b",  # sizes, memory
    r"(?<![\w.])\d+(?:\.\d+)?\.(?:[KMGT]B)\b",  # Nextflow memory (54.GB)
    r"(?<![\w.])\d+(?:\.\d+)?\s*(?:ms|s|sec|min|minutes?|h|hours?|d|days?)\b",  # durations
    r"(?i)\bexit(?:[ _]?code)?[\s:=]*\d+\b",  # exit codes
    r"\b\d+:\d+\b",  # sacct ExitCode
    r"(?i)(?<![\w.])\d+\s*CPUs?\b|--threads[ =]\d+|AllocCPUS[\s:=]*\d+",  # CPU counts
    r"(?i)\battempts?\s+\d+\b",
    r"\b20261014\b",  # the seed
    r"(?i)\b(?:version|v|db|bakta)\s*\d+(?:\.\d+)+\b|\b\d+\.\d+\.\d+\b",  # versions
    r"§\d+(?:\.\d+)*",
    r"\b(?=[0-9a-f]*\d)(?=[0-9a-f]*[a-f])[0-9a-f]{7,64}\b|\b[0-9a-f]{64}\b",  # digests, shas
    r"\b\d{6,9}\b",  # Slurm job IDs
    r"[\w.\-/]*[A-Za-z_][\w.\-/]*",  # identifiers and paths with a letter in them
    r"(?m)^\s*\d+\.\s",  # list ordinals
    r"(?<![\w.,])\d(?![\w.,%])",  # a lone digit
]


def non_provenance_numbers(text: str) -> list[str]:
    """Digits left after removing every allowed provenance form (empty = the text passes)."""
    for pat in _ALLOWED:
        text = re.sub(pat, " ", text, flags=re.S)
    return re.findall(r"\S*\d\S*", text)


def section(md: str, title: str) -> str | None:
    m = re.search(rf"^(#+)\s+{re.escape(title)}\s*$", md, re.M)
    if m is None:
        return None
    level = len(m.group(1))
    rest = md[m.end() :]
    nxt = re.search(rf"^#{{1,{level}}}\s", rest, re.M)
    return rest[: nxt.start()] if nxt else rest


def _tracked(*paths: str) -> list[str]:
    res = subprocess.run(["git", "-C", str(REPO), "ls-files", "-z", "--", *paths],
                         capture_output=True, text=True, check=True)  # fmt: skip
    return [p for p in res.stdout.split("\0") if p]


def _forbidden(path: str) -> bool:
    p = PurePosixPath(path)
    if any(fnmatch.fnmatchcase(part, pat) for part in p.parts for pat in FORBIDDEN):
        return True
    return "evidence" in p.parts[:-1] and p.suffix == ".json"


def test_no_run_output_is_tracked_under_validation_pilot() -> None:
    # the patterns catch every name the runs and the compare job write (self-check)
    for name in ("validation/pilot/out_a1/run_manifest.json", "validation/pilot/out_a2.r2/x",
                 "validation/pilot/cmp_a1_a2/tables/compare_detection.tsv",
                 "validation/pilot/evidence/a1.json", "validation/pilot/evidence/compare.json",
                 "validation/pilot/ATTEMPTS.tsv", "validation/pilot/x/report/report.html",
                 "validation/pilot/report/compare.html"):  # fmt: skip
        assert _forbidden(name), name
    assert not _forbidden("validation/pilot/run_evidence.py")
    assert not _forbidden("validation/pilot/db_a2/makedb_manifest.json")
    tracked = _tracked("validation/pilot/", "scripts/slurm/")
    missing = [d for d in DELIVERABLES if d not in tracked]
    assert not missing, f"M10-04 deliverables not tracked: {missing}"
    bad = [p for p in tracked if p.startswith("validation/pilot/") and _forbidden(p)]
    assert not bad, f"run outputs are never committed (Khoa 2026-10-23): {bad}"


def test_readme_pass1_section_holds_provenance_numbers_only() -> None:
    # the allow-list keeps provenance and rejects biological numbers (self-check)
    ok = (
        "| A2 | 1873779 | COMPLETED | 0:0 | 05:59:12 | 2104832K | wheat-07 | `0123abc` |\n"
        "Attempt 2 of `out_a1.r2` (`resume-infra`), exit code 4, 14 CPUs, 54 GB, seed 20261014.\n"
        "`report.html`: sha256 `" + "ab12" * 16 + "`, 2,145,332 bytes; on 2026-10-24.\n"
        "Bakta 1.12.2, DB 6.0, §15.2; 5 h 59 min; `--threads 14 --memory 54.GB`.\n"
    )
    assert non_provenance_numbers(ok) == [], non_provenance_numbers(ok)
    for bio in ("A2 detected 1,234 genes.", "Spearman rho = 0.82 between arms.",
                "island Jaccard 0.41", "62% of peptides mapped", "A1: 812 genes"):  # fmt: skip
        assert non_provenance_numbers(bio), bio
    body = section(README.read_text(encoding="utf-8"), SECTION)
    assert body is not None, f"validation/pilot/README.md has no section {SECTION!r} (M10-04)"
    assert body.strip(), "the section is empty"
    bad = non_provenance_numbers(body)
    assert not bad, f"numbers outside the provenance fields (Q-K1): {bad}"
