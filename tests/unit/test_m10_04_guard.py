"""Guard: M10-04 commits no run output and no biological number (acceptance criterion 7).

Spec: docs/specs/M10-04-pilot-runs.md, "What is committed (and what is not)"; Khoa's rule of
2026-10-23 (nothing derived from the DIA-NN reports is committed); Q-K1 parked with its default
(lead decision 2): the README section holds provenance only. Extends M10-03's guard.

The README rule is an allow-list, as the spec asks (tightened by the lead's decision on review
finding B3). A number in the section "Pass-1 runs (M10-04)" is allowed only as:

* part of an identifier or path (a token with a letter in it: ``out_a1``, ``CP171242.1``,
  ``M10-04``, ``ADR-0046``, ``wheat-07``, ``.r2``, ``2nd``), a ``§`` reference, a version
  (``1.12.2``, ``DB 6.0``);
* a hexadecimal digest or git sha (7-64 hex characters with a letter);
* a **Slurm job ID** (6-9 digits) only in context: ``Slurm`` or ``job``/``jobs``/``JobID`` in the
  same clause (clauses end at ``.``/``;``/``!``/``?`` before a space, a line end or a table
  ``|``), or in a table column whose header has ``job``. When a committed ``ATTEMPTS``-style list
  exists under ``validation/pilot/``, the ID must also appear in it;
* a byte size or memory (``2,145,332 bytes``, ``54 GB``, ``2104832K``, ``54.GB``); a bare unit
  letter K/M/G/T only with no space (``3012 M`` is a count, not a size);
* a duration (``05:59:12``, ``2-00:00:00``, ``5 h 59 min``, ``12 s``), a date or timestamp;
* an exit code (``0:0``, ``exit 1``, ``exit code 4``, ``exit_code 0``);
* a CPU count (``14 CPUs``, ``--threads 14``, ``AllocCPUS 14``); the seed ``20261014``;
* a lone digit only after ``attempt``, ``exit``, ``step``, ``§`` or ``ADR-``, as an ordinal, in a
  table column headed attempt/exit/step, or as a Markdown list ordinal at the start of a line.

Fenced code blocks are read too. A code line (backslash continuations joined) that starts with
an allow-listed command word (``sbatch gippa python nextflow sacct bash squeue``) has its numeric
flag values (``--threads 14``, ``-j 1873779``) allowed; the rest of that line, and every other
code line, is scanned like prose. Anything else with a digit fails: counts of genes, peptides or
precursors, percentages, rho, Jaccards.
"""

from __future__ import annotations

import fnmatch
import re
import subprocess
from pathlib import PurePosixPath

from helpers.slurm_sandbox import REPO
from hypothesis import given
from hypothesis import strategies as st

README = REPO / "validation" / "pilot" / "README.md"
SECTION = "Pass-1 runs (M10-04)"
FORBIDDEN = ("out_a*", "cmp_*", "ATTEMPTS.tsv", "report.html", "compare.html")
DELIVERABLES = ("validation/pilot/run_evidence.py", "scripts/slurm/pilot_compare.sbatch")

_ALLOWED = [
    r"\b\d{4}-\d{2}-\d{2}(?:[T ][0-9:.]+(?:Z|[+-]\d{2}:?\d{2})?)?",  # dates, timestamps
    r"\b(?:\d+-)?\d{1,3}:\d{2}:\d{2}(?:\.\d+)?\b",  # Slurm durations
    r"(?<![\w.])\d[\d,]*(?:\.\d+)?(?:\s*(?:bytes?|[KMGT]i?B)|[KMGT])\b",  # sizes, memory
    r"(?<![\w.])\d+(?:\.\d+)?\.(?:[KMGT]B)\b",  # Nextflow memory (54.GB)
    r"(?<![\w.])\d+(?:\.\d+)?\s*(?:ms|s|sec|min|minutes?|h|hours?|d|days?)\b",  # durations
    r"(?i)\bexit(?:[ _]?code)?[\s:=]*\d+\b",  # exit codes
    r"\b\d+:\d+\b",  # sacct ExitCode
    r"(?i)(?<![\w.])\d+\s*CPUs?\b|--threads[ =]\d+|AllocCPUS[\s:=]*\d+",  # CPU counts
    r"(?i)\b(?:attempts?|steps?)\s+\d+\b",
    r"\b20261014\b",  # the seed
    r"(?i)\b(?:version|v|db|bakta)\s*\d+(?:\.\d+)+\b|\b\d+\.\d+\.\d+\b",  # versions
    r"§\d+(?:\.\d+)*",
    r"\bADR-\d+\b",
    r"\b(?=[0-9a-f]*\d)(?=[0-9a-f]*[a-f])[0-9a-f]{7,64}\b|\b[0-9a-f]{64}\b",  # digests, shas
    r"[\w.\-/]*[A-Za-z_][\w.\-/]*",  # identifiers, paths, ordinals (2nd), .r2
    r"(?m)^\s*\d+\.\s",  # list ordinals
]
_CMD_WORD = re.compile(r"(?:sbatch|gippa|python3?|nextflow|sacct|bash|squeue)")
_FLAG_VALUE = re.compile(r"(?<!\S)(--?[A-Za-z][\w.-]*)(?:=|\s+)\d[\w.:,-]*")
_FENCE = re.compile(r"^[ \t]*(```|~~~)[^\n]*\n(.*?)^[ \t]*\1[ \t]*$", re.M | re.S)
_JOB_WORD = re.compile(r"(?i)\b(?:slurm|job(?:s|ids?)?)\b")
_JOB_ID = re.compile(r"(?<![\w.,:-])\d{6,9}(?![\w,:-]|\.\d)")
_CLAUSE_END = re.compile(r"([.;!?](?=\s|$)|\n|\|)")
_LONE = re.compile(r"(?<![\w.,])\d(?![\w.,%])")


def _code_block(body: str) -> str:
    """Each logical line of a fenced block: flag values dropped if it starts with a command."""
    out = []
    for line in re.sub(r"\\\n", " ", body).splitlines():
        words = re.sub(r"^\s*\$\s+", "", line).split()
        first = words[0].strip("\"'").rsplit("/", 1)[-1] if words else ""
        if _CMD_WORD.fullmatch(first):
            line = _FLAG_VALUE.sub(r"\1", line)
        out.append(line)
    return "\n".join(out) + "\n"


def _tables(text: str, known: set[str] | None) -> str:
    """Cells under a ``job`` header lose their job IDs; under attempt/exit/step, lone digits."""
    lines, out, header = text.split("\n"), [], None
    for i, line in enumerate(lines):
        if not line.lstrip().startswith("|"):
            header = None
            out.append(line)
            continue
        cells = line.strip().strip("|").split("|")
        nxt = lines[i + 1] if i + 1 < len(lines) else ""
        if header is None and re.fullmatch(r"\s*\|?(?:\s*:?-+:?\s*\|)+\s*:?-*:?\s*", nxt):
            header = cells
        elif header is not None and not re.fullmatch(r"[\s|:-]*", line):
            for j, cell in enumerate(cells):
                h = header[j] if j < len(header) else ""
                if _JOB_WORD.search(h):
                    cell = _drop_job_ids(cell, known, force=True)
                if re.search(r"(?i)attempt|exit|step", h):
                    cell = _LONE.sub(" ", cell)
                cells[j] = cell
            line = "|" + "|".join(cells) + "|"
        out.append(line)
    return "\n".join(out)


def _drop_job_ids(clause: str, known: set[str] | None, *, force: bool = False) -> str:
    if not force and not _JOB_WORD.search(clause):
        return clause
    return _JOB_ID.sub(lambda m: " " if known is None or m.group() in known else m.group(), clause)


def non_provenance_numbers(text: str, job_ids: set[str] | None = None) -> list[str]:
    """Digits left after removing every allowed provenance form (empty = the text passes).

    ``job_ids``: the committed ``ATTEMPTS``-style job IDs, when such a list exists (else None:
    the context rule alone).
    """
    text = _FENCE.sub(lambda m: _code_block(m.group(2)), text)
    text = _tables(text, job_ids)
    parts = _CLAUSE_END.split(text)
    text = "".join(p if i % 2 else _drop_job_ids(p, job_ids) for i, p in enumerate(parts))
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


#: provenance the README section legitimately holds (must pass)
_LEGIT = (
    "| Arm | Job | State | ExitCode | Elapsed | MaxRSS | Node | Commit |\n"
    "|---|---|---|---|---|---|---|---|\n"
    "| A2 | 1873779 | COMPLETED | 0:0 | 05:59:12 | 2104832K | wheat-07 | `0123abc` |\n",
    "Attempt 2 of `out_a1.r2` (`resume-infra`), exit code 4, 14 CPUs, 54 GB, seed 20261014.\n",
    "`report.html`: sha256 `" + "ab12" * 16 + "`, 2,145,332 bytes; on 2026-10-24.\n",
    "Bakta 1.12.2, DB 6.0, §15.2; 5 h 59 min; `--threads 14 --memory 54.GB`.\n",
    "The Slurm job 1873779 ended COMPLETED; the compare job was 1873801.\n",
    "Head jobs 1873779 and 1873790 (A2, A1), JobID 1873801.\n",
    "| Arm | Attempt | Head job | Class | Exit |\n|---|---|---|---|---|\n"
    "| A2 | 1 | 1873779 | fresh | 0 |\n| A2 | 2 | 1873790 | resume-infra | 0 |\n",
    "1. A2 first, then A1 (step 8, ADR-0046, the 2nd attempt, `out_a2.r2`).\n",
    "```bash\nsbatch -J gippa-run-pilot-a2 scripts/slurm/gippa_run.sbatch \\\n"
    "  --genome $G --threads 14 --memory 54.GB --outdir $R/out_a2\n"
    "sacct -P -n -j 1873779 -o JobID,State,ExitCode\n"
    "squeue -h -u $USER -n gippa-run-pilot-a2\n```\n",
    "MaxRSS 2104832K, 1.5G and 3012M; 54 GB requested.\n",
)

#: biological numbers (must be rejected), including the reviewer's leaked strings (B3)
_BIO = (
    "A2 detected 1,234 genes.",
    "Spearman rho = 0.82 between arms.",
    "island Jaccard 0.41",
    "62% of peptides mapped",
    "A1: 812 genes",
    "1234567 peptides mapped",
    "123456 precursors",
    "A1: 3 genes detected",
    "3012 M peptides",
    "2104 K genes",
    "```\nwc -l genes.tsv\n4983\n```",
    "~~~\n4983 genes\n~~~",
    "```bash\nsbatch --threads 14 x.sbatch 4983 genes\n```",
    "```bash\ngzip -dc a2.tsv.gz | wc -l\n123456\n```",
    "Slurm job 1873779; 1234567 peptides mapped.",
    "Slurm job 1873779.\n1234567 peptides mapped.",
    "| Arm | Peptides |\n|---|---|\n| A2 | 1234567 |\n",
    "| Arm | Genes |\n|---|---|\n| A2 | 3 |\n",
    "    4983 genes (indented code)",
)


def test_guard_allow_list_self_check() -> None:
    """The allow-list keeps provenance and rejects biological numbers (B3's leaks included)."""
    for ok in _LEGIT:
        assert non_provenance_numbers(ok) == [], (ok, non_provenance_numbers(ok))
    for bio in _BIO:
        assert non_provenance_numbers(bio), bio


def test_guard_job_ids_cross_checked_against_a_committed_list() -> None:
    text = "The Slurm job 1873779 ended COMPLETED.\n"
    assert non_provenance_numbers(text, job_ids={"1873779"}) == []
    assert non_provenance_numbers(text, job_ids={"1873790"}) == ["1873779"]
    table = "| Arm | Job |\n|---|---|\n| A2 | 1873779 |\n"
    assert non_provenance_numbers(table, job_ids={"1873790"}) == ["1873779"]


@given(st.integers(min_value=100_000, max_value=999_999_999).filter(lambda n: n != 20261014))
def test_guard_job_id_needs_its_context(n: int) -> None:
    """Property: a 6-9 digit number passes next to ``Slurm job`` and fails next to a noun (the
    calibration seed 20261014 is allow-listed everywhere, so it is left out)."""
    assert non_provenance_numbers(f"Slurm job {n} ended COMPLETED.") == []
    assert non_provenance_numbers(f"| Head job |\n|---|\n| {n} |\n") == []
    for noun in ("peptides", "precursors", "genes mapped"):
        assert non_provenance_numbers(f"{n} {noun}") == [str(n)]


@given(st.integers(min_value=0, max_value=99_999), st.sampled_from(["K", "M", "G", "T"]))
def test_guard_bare_unit_letters_need_no_space(n: int, unit: str) -> None:
    """Property: ``<n><unit>`` is a size; ``<n> <unit> genes`` is a count."""
    assert non_provenance_numbers(f"MaxRSS {n}{unit}") == []
    assert non_provenance_numbers(f"{n} {unit} genes") == [str(n)]


@given(st.integers(min_value=0, max_value=9))
def test_guard_lone_digits_need_their_word(d: int) -> None:
    for ok in (f"attempt {d}", f"exit {d}", f"step {d}", f"§{d}", f"ADR-{d:04d}", f"`out_a1.r{d}`"):
        assert non_provenance_numbers(ok) == [], ok
    for bio in (f"A1: {d} genes detected", f"{d} islands", f"| A2 | {d} |"):
        assert non_provenance_numbers(bio) == [str(d)], bio


def _committed_job_ids() -> set[str] | None:
    """Job IDs of a committed ``ATTEMPTS``-style list under ``validation/pilot/`` (None if none)."""
    files = [p for p in _tracked("validation/pilot/")
             if re.search(r"(?i)attempts", PurePosixPath(p).name)]  # fmt: skip
    if not files:
        return None
    text = "\n".join((REPO / p).read_text(encoding="utf-8") for p in files)
    return set(re.findall(r"\b\d{6,9}\b", text))


def test_readme_pass1_section_holds_provenance_numbers_only() -> None:
    body = section(README.read_text(encoding="utf-8"), SECTION)
    assert body is not None, f"validation/pilot/README.md has no section {SECTION!r} (M10-04)"
    assert body.strip(), "the section is empty"
    bad = non_provenance_numbers(body, _committed_job_ids())
    assert not bad, f"numbers outside the provenance fields (Q-K1): {bad}"
