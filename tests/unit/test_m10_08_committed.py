"""What M10-08a commits, and what it must never commit (acceptance criteria 11-12, deliverable 11).

Spec: docs/specs/M10-08-a3-feasibility.md. R5 spirit: Cascadia, TIMSCONVERT and pyTDFSDK are never
imported from ``src/`` or ``validation/`` (only ``scripts/pilot/cascadia_offline.py`` imports
Cascadia, inside its env). R11: no weights, image, mzML, ``.ssl``, raw data, TDF-SDK library or
Unimod copy is tracked; nothing above 5 MB under ``validation/`` or the recorded Cascadia
directory, which holds only ``PROVENANCE.md``. ADR-M10-B is drafted in 08a (Status: proposed,
Needs Khoa: yes). These read only tracked files. No statistical claim, so no calibration test.

Note: the import and ``git ls-files`` checks are regression guards; they already hold before
08a lands and must keep holding after it.
"""

from __future__ import annotations

import ast
import re
import subprocess
import sys
from pathlib import Path

from helpers.a3_kit import LAUNCHER, require
from helpers.slurm_sandbox import REPO

R11_MAX_BYTES = 5 * 1024 * 1024
FORBIDDEN_IMPORTS = {"cascadia", "timsconvert", "pytdfsdk"}


def _tracked(*prefixes: str) -> list[str]:
    res = subprocess.run(["git", "-C", str(REPO), "ls-files", "--", *prefixes],
                         capture_output=True, text=True, check=True)  # fmt: skip
    return [ln for ln in res.stdout.splitlines() if ln]


def _imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(), filename=str(path))
    out: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            out |= {a.name.split(".")[0].lower() for a in node.names}
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            out.add(node.module.split(".")[0].lower())
    return out


def test_check_imports_passes() -> None:
    cmd = [sys.executable, str(REPO / "scripts" / "ci" / "check_imports.py"), str(REPO / "src")]
    res = subprocess.run(cmd, capture_output=True, text=True, check=False)
    assert res.returncode == 0, res.stdout


def test_src_and_validation_never_import_pilot_tools() -> None:
    bad = []
    for top in ("src", "validation"):
        for path in sorted((REPO / top).rglob("*.py")):
            hit = _imports(path) & FORBIDDEN_IMPORTS
            if hit:
                bad.append(f"{path.relative_to(REPO)}: {sorted(hit)}")
    assert not bad, bad


def test_only_the_launcher_imports_cascadia() -> None:
    assert "cascadia" in _imports(require(LAUNCHER, "deliverable 4"))
    others = [p for p in _tracked("scripts") if p.endswith(".py")
              and (REPO / p).resolve() != LAUNCHER.resolve()
              and "cascadia" in _imports(REPO / p)]  # fmt: skip
    assert not others, others


_NEVER = re.compile(
    r"(\.(ckpt|sif|mzml|ssl|d|tdf|tdf_bin)$)|((^|/)(?!changelog\.d/)[^/]+\.d/)"
    r"|((^|/)libtimsdata[^/]*$)|((^|/)unimod_tables\.xml)",
    re.I,
)


def test_no_binary_a3_artefact_is_tracked() -> None:
    bad = [p for p in _tracked() if _NEVER.search(p)]
    assert not bad, bad


def test_no_large_file_under_validation_or_recorded_cascadia() -> None:
    big = [p for p in _tracked("validation", "tests/data/recorded/cascadia")
           if (REPO / p).stat().st_size > R11_MAX_BYTES]  # fmt: skip
    assert not big, big


def test_recorded_cascadia_holds_only_provenance() -> None:
    files = _tracked("tests/data/recorded/cascadia")
    assert all(Path(p).name == "PROVENANCE.md" for p in files), files
    assert all(p.startswith("tests/data/recorded/cascadia/0.0.7/") for p in files), files


def test_adr_m10_b_drafted() -> None:
    adrs = sorted((REPO / "docs" / "adr").glob("ADR-*-a3-external-artefacts.md"))
    assert len(adrs) == 1, "deliverable 11: docs/adr/ADR-NNNN-a3-external-artefacts.md"
    text = adrs[0].read_text()
    assert re.match(r"ADR-\d{4}-", adrs[0].name)
    assert re.search(r"^- Status: proposed", text, re.M)
    assert re.search(r"^- Needs Khoa: yes", text, re.M)
    for needle in ("TDF-SDK", "drive.google.com", "psims", "B2", "B1"):
        assert needle in text, needle
