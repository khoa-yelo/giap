"""Nothing of the A3 run is committed (M10-10 criterion 11; Khoa's data rules, Q10).

Spec: docs/specs/M10-10-a3-run.md, "Deliverables" (nothing else is committed) and criterion 11.
A regression guard: it passes today and must keep passing after the runs. The committed M10-08
files ``validation/pilot/metrics/a3_feasibility.json`` and ``validation/pilot/a3_assets.tsv``
(which the spec's ``a3*.tsv*`` pattern would otherwise catch) are exempt by path.
"""

from __future__ import annotations

import fnmatch
import re
import subprocess
from pathlib import PurePosixPath

from helpers.slurm_sandbox import REPO

FORBIDDEN = (
    "*.ssl*", "*.mzML", "a3*.tsv*", "a3.convert.*", "*feasibility.json", "out_a3*",
    "cmp_a1_a2_a3*", "compare3.json", "ssl.sha256",
)  # fmt: skip
EXEMPT = {"validation/pilot/metrics/a3_feasibility.json", "validation/pilot/a3_assets.tsv"}
PEPTIDE_LIKE = re.compile(r"[ACDEFGHIKLMNPQRSTVWY]{7,}")
README = REPO / "validation" / "pilot" / "README.md"


def _files() -> list[str]:
    out = subprocess.run(["git", "-C", str(REPO), "ls-files"], capture_output=True, text=True,
                         check=True).stdout  # fmt: skip
    return out.splitlines()


def test_no_a3_run_output_is_committed() -> None:
    bad = []
    for f in _files():
        if f in EXEMPT:
            continue
        parts = PurePosixPath(f).parts
        if any(fnmatch.fnmatchcase(p, pat) for p in parts for pat in FORBIDDEN):
            bad.append(f)
    assert bad == [], f"A3 outputs are never committed: {bad}"


def _a3_section(text: str) -> str | None:
    m = re.search(r"^(#+)\s*A3 run \(M10-10\).*$", text, re.M)
    if m is None:
        return None
    level = len(m.group(1))
    rest = text[m.end() :]
    end = re.search(rf"^#{{1,{level}}}\s", rest, re.M)
    return rest[: end.start()] if end else rest


def test_readme_a3_section_has_no_peptide_like_string() -> None:
    """Checked once the section exists (it is written after the runs)."""
    section = _a3_section(README.read_text(encoding="utf-8"))
    if section is None:
        return
    hits = PEPTIDE_LIKE.findall(section)
    assert hits == [], f"peptide-like strings in README > A3 run (M10-10): {hits}"


def test_section_finder_and_pattern() -> None:
    """The guard itself: it finds a section and flags a sequence in it."""
    text = "# Pilot\n## A3 run (M10-10)\nrows 31204, PEPTLDEKR\n## Next\nSAMPLLERK\n"
    section = _a3_section(text)
    assert section is not None and "SAMPLLERK" not in section
    assert PEPTIDE_LIKE.findall(section) == ["PEPTLDEKR"]
    assert fnmatch.fnmatchcase("a3.tsv.gz", "a3*.tsv*")
    assert fnmatch.fnmatchcase("WT_1.ssl", "*.ssl*")
