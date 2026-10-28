"""``preflight.write_outputs`` keeps an unchanged report file (M7-15 item 3; review #278 B2, the
"preflight keep-unchanged call" survivor): writing the same report twice leaves
``OUT/preflight/preflight.{json,md}`` as the earlier files (same inode and mtime), so a
``--resume`` run whose consumers hash path and mtime caches; a changed report replaces them.
"""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from gippa.config import resolve_config
from gippa.preflight.environment import GIB, EnvProbe, requirements
from gippa.preflight.run import PreflightInputs, run_preflight, write_outputs

LOW = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]
PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


@pytest.fixture(autouse=True)
def _isolated(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for name, sub in (("GIPPA_DB_DIR", "db"), ("GIPPA_REPO", "repo"), ("GIPPA_ROOT", "repo")):
        (tmp_path / sub).mkdir(exist_ok=True)
        monkeypatch.setenv(name, str(tmp_path / sub))


def _probe() -> EnvProbe:
    found: dict[str, str | None] = {
        r.name: r.min_version or "1.0"
        for r in requirements(["core", "mobile", "bgc", "function"], "genome")
    }
    return EnvProbe(found=found, free_disk_bytes=500 * GIB, total_ram_bytes=64 * GIB)


def _report(tmp_path: Path, peptide: str = "PEPTIDEK"):  # type: ignore[no-untyped-def]
    g = tmp_path / "g.fa"
    g.write_text(">chr1 circular=true\n" + ("ACGGTCAT" * 18_750) + "\n")
    p = tmp_path / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    pep = tmp_path / "peptides.tsv"
    pep.write_text(f"peptide\tprotein_ids\n{peptide}\tP1\n")
    inputs = PreflightInputs(genome=g, proteins=p, peptides=pep)
    return run_preflight(inputs, resolve_config(overrides=LOW), _probe(), "preflight")


def _stamp(path: Path) -> tuple[int, int]:
    st = path.stat()
    return st.st_ino, st.st_mtime_ns


def test_unchanged_report_files_are_the_earlier_files(tmp_path: Path) -> None:
    report, norm = _report(tmp_path)
    out = tmp_path / "out"
    write_outputs(out, report, norm)
    rep = out / "preflight"
    names = ("preflight.json", "preflight.md")
    old = {n: rep / n for n in names}
    past = 1_000_000_000_000_000_000  # 2001-09-09: any rewrite would get a newer mtime
    for p in old.values():
        os.utime(p, ns=(past, past))
    first = {n: _stamp(p) for n, p in old.items()}
    content = {n: p.read_bytes() for n, p in old.items()}
    write_outputs(out, report, norm)
    for n in names:
        assert (rep / n).read_bytes() == content[n], n
        assert _stamp(rep / n) == first[n], f"{n} rewritten although unchanged"
    assert sorted(p.name for p in out.iterdir()) == ["preflight", "work"], "temp dir left"


def test_changed_report_replaces_the_files(tmp_path: Path) -> None:
    report, norm = _report(tmp_path)
    out = tmp_path / "out"
    write_outputs(out, report, norm)
    js = out / "preflight" / "preflight.json"
    past = 1_000_000_000_000_000_000
    os.utime(js, ns=(past, past))
    before = js.read_bytes()
    js.write_bytes(before + b" ")  # a different earlier report
    os.utime(js, ns=(past, past))
    write_outputs(out, report, norm)
    assert js.read_bytes() == before
    assert js.stat().st_mtime_ns != past
