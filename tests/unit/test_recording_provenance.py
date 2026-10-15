"""Recorded-fixture layout and provenance validator (M2-10; GIAP_SPEC.md §11.4, R11).

§11.4: "Recorded native outputs live in `tests/data/recorded/<tool>/<version>/<input>/`, with a
`PROVENANCE.md` giving the command, environment lock or digest, date and host." M2-10's
`record_fixture.py` writes one ``PROVENANCE.md`` per ``<input>`` directory with the command line,
env lock sha256 or image digest, DB name+version, input SHA-256s, date (UTC) and host.

Machine-checkable convention pinned here (test author): a wave-A ``PROVENANCE.md`` lists every
committed file of its directory (recursively, except itself) in ``sha256sum`` format,
``<64 hex>  <path relative to the directory>``, and those sums must match the committed bytes.
Every ``sha256sum``-format line in any recorded ``PROVENANCE.md`` whose path exists must match.
"""

from __future__ import annotations

import gzip
import hashlib
import re
from dataclasses import dataclass
from pathlib import Path

import pytest

DATA = Path(__file__).resolve().parents[1] / "data"
RECORDED = DATA / "recorded"
MAX_FILE = 5 * 1024 * 1024
MAX_TOTAL = 50 * 1024 * 1024

SHA_LINE = re.compile(r"^\s*([0-9a-f]{64})\s+\*?(\S+)\s*$", re.MULTILINE)
HEX64 = re.compile(r"\b[0-9a-f]{64}\b")
DATE = re.compile(r"\b20\d\d-\d\d-\d\d(?!\d)")
HOST = re.compile(r"(?im)^[\s>*|-]*\**host\**\b\W*\S")
COMMAND = re.compile(r"(?im)^[\s>*|#-]*\**command")
DIGEST = re.compile(r"sha256:[0-9a-f]{64}")
LOCK = re.compile(r"locks/\S+\.explicit\.txt")


@dataclass(frozen=True)
class Leaf:
    tool: str
    version: str  # regex for the <version> directory
    inputs: tuple[str, ...]  # each a regex for the <input> directory
    executable: str  # must appear in the recorded command


# Wave A (spec "Fixture layout produced"). geNOMAD/antiSMASH accept the sample labels already
# committed by M2-18 (k12_mg1655, nissle1917). BacTermFinder (opt-in, env not built) stays with its
# pending_sample contract test.
WAVE_A = (
    Leaf("bakta", r"1\.12\.2-light", ("sim100k",), "bakta"),
    Leaf("bakta", r"1\.12\.2-full", ("sim100k", "k12", "nissle"), "bakta"),
    Leaf("terminatornet", r"[0-9a-f]{12}", ("sim100k", "k12"), "TerminatorNet.py"),
    Leaf(
        "genomad", r"1\.12\.0(-db1\.9)?", ("sim100k", "k12(_mg1655)?", "nissle(1917)?"), "genomad"
    ),
    Leaf("antismash", r"8\.0\.4", ("sim100k", "k12(_mg1655)?", "nissle(1917)?"), "antismash"),
)


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _files(d: Path) -> list[Path]:
    return sorted(p for p in d.rglob("*") if p.is_file() and p.name != "PROVENANCE.md")


def _match(parent: Path, pattern: str) -> list[Path]:
    rx = re.compile(pattern)
    if not parent.is_dir():
        return []
    return sorted(p for p in parent.iterdir() if p.is_dir() and rx.fullmatch(p.name))


def _label(pattern: str) -> str:
    """Readable test id for a directory regex (`1\\.12\\.0(-db1\\.9)?` -> `1.12.0`)."""
    return re.sub(r"\(.*?\)\?", "", pattern).replace("\\", "").replace("[0-9a-f]{12}", "git12")


def _wave_a_dirs() -> list[tuple[str, str, Path | None]]:
    """(id, executable, dir or None when missing) for every expected wave-A input directory."""
    out: list[tuple[str, str, Path | None]] = []
    for leaf in WAVE_A:
        versions = _match(RECORDED / leaf.tool, leaf.version)
        for inp in leaf.inputs:
            hits = [d for v in versions for d in _match(v, inp)]
            out.append(
                (
                    f"{leaf.tool}/{_label(leaf.version)}/{_label(inp)}",
                    leaf.executable,
                    hits[0] if hits else None,
                )
            )
    return out


def provenance_problems(d: Path, executable: str) -> list[str]:
    """Everything wrong with ``d/PROVENANCE.md`` under the §11.4 rule and the sha256 convention."""
    prov = d / "PROVENANCE.md"
    if not prov.is_file():
        return [f"{d}: no PROVENANCE.md"]
    text = prov.read_text()
    bad: list[str] = []
    if not COMMAND.search(text) or executable not in text:
        bad.append(f"{prov}: no command line with `{executable}`")
    if not (DIGEST.search(text) or (LOCK.search(text) and HEX64.search(text))):
        bad.append(
            f"{prov}: no image digest (sha256:...) or env lock (locks/<tool>.explicit.txt + sha256)"
        )
    if not DATE.search(text):
        bad.append(f"{prov}: no date (YYYY-MM-DD)")
    if not HOST.search(text):
        bad.append(f"{prov}: no `host` line")
    if not re.search(r"(?i)\binputs?\b", text):
        bad.append(f"{prov}: no input SHA-256s")
    listed = {p: h for h, p in SHA_LINE.findall(text)}
    for f in _files(d):
        rel = f.relative_to(d).as_posix()
        if rel not in listed:
            bad.append(f"{prov}: `{rel}` has no sha256 line")
        elif listed[rel] != _sha(f):
            bad.append(f"{prov}: sha256 of `{rel}` does not match")
    # sums of uncommitted sources (inputs, full recordings) may be listed too; they are not checked
    return bad


# ----------------------------------------------------------------- the validator itself


def _toy(tmp: Path, prov: str, files: dict[str, bytes]) -> Path:
    d = tmp / "tool" / "1.0" / "sim100k"
    d.mkdir(parents=True)
    for name, data in files.items():
        (d / name).parent.mkdir(parents=True, exist_ok=True)
        (d / name).write_bytes(data)
    (d / "PROVENANCE.md").write_text(prov)
    return d


def _good_prov(files: dict[str, bytes]) -> str:
    sums = "\n".join(f"{hashlib.sha256(b).hexdigest()}  {n}" for n, b in files.items())
    return (
        "# tool 1.0 on sim100k\n"
        "- Command: `tool --in genome.fa --out out`\n"
        f"- Image: docker://x/tool@sha256:{'a' * 64}\n"
        "- Date: 2026-10-15T03:00:00Z\n"
        "- Host: iron-01 (node iron-01, SLURM_JOB_ID 1)\n"
        f"- Inputs: genome.fa sha256 {'b' * 64}\n\n"
        f"## Files\n```\n{sums}\n```\n"
    )


def test_validator_accepts_a_complete_provenance(tmp_path: Path) -> None:
    files = {"out.tsv.gz": gzip.compress(b"x\n", mtime=0), "sub/a.tsv": b"a\n"}
    d = _toy(tmp_path, _good_prov(files), files)
    assert provenance_problems(d, "tool") == []


def test_validator_rejects_a_changed_file(tmp_path: Path) -> None:
    files = {"out.tsv": b"x\n"}
    d = _toy(tmp_path, _good_prov(files), files)
    (d / "out.tsv").write_bytes(b"y\n")
    assert any("does not match" in p for p in provenance_problems(d, "tool"))


def test_validator_rejects_an_unlisted_file(tmp_path: Path) -> None:
    files = {"out.tsv": b"x\n"}
    d = _toy(tmp_path, _good_prov(files), files)
    (d / "extra.tsv").write_bytes(b"z\n")
    assert any("`extra.tsv` has no sha256 line" in p for p in provenance_problems(d, "tool"))


@pytest.mark.parametrize("drop", ["Command", "Image", "Date", "Host", "Inputs"])
def test_validator_rejects_a_missing_field(tmp_path: Path, drop: str) -> None:
    files = {"out.tsv": b"x\n"}
    prov = "\n".join(ln for ln in _good_prov(files).splitlines() if not ln.startswith(f"- {drop}"))
    d = _toy(tmp_path, prov, files)
    assert provenance_problems(d, "tool"), drop


# ----------------------------------------------------------------- the committed fixtures

WAVE_A_DIRS = _wave_a_dirs()


@pytest.mark.parametrize(("ident", "exe", "d"), WAVE_A_DIRS, ids=[w[0] for w in WAVE_A_DIRS])
def test_wave_a_dir_exists_with_complete_provenance(ident: str, exe: str, d: Path | None) -> None:
    assert d is not None, f"missing tests/data/recorded/{ident}/ (M2-10 wave A, §11.4 layout)"
    assert _files(d), f"{d} holds no recorded files"
    assert provenance_problems(d, exe) == []


def test_terminatornet_version_dir_is_the_pinned_commit() -> None:
    """`terminatornet/<git12>/`: the first 12 hex digits of the pinned commit (ADR, M2-01)."""
    assert (RECORDED / "terminatornet" / "721803095459").is_dir()


def test_every_recorded_file_has_a_dated_provenance() -> None:
    """Older fixtures may keep one PROVENANCE.md per tool directory: the nearest one up to
    `recorded/<tool>/` counts, and it must at least carry a date."""
    missing: list[str] = []
    for f in sorted(RECORDED.rglob("*")):
        if not f.is_file() or f.name == "PROVENANCE.md":
            continue
        tool_dir = RECORDED / f.relative_to(RECORDED).parts[0]
        d = f.parent
        while not (d / "PROVENANCE.md").is_file() and d != tool_dir:
            d = d.parent
        prov = d / "PROVENANCE.md"
        if not prov.is_file() or not DATE.search(prov.read_text()):
            missing.append(str(f.relative_to(RECORDED)))
    assert missing == []


def test_listed_sha256_sums_match_committed_files() -> None:
    wrong: list[str] = []
    for prov in [*sorted(RECORDED.rglob("PROVENANCE.md")), DATA / "genomes" / "PROVENANCE.md"]:
        if not prov.is_file():
            continue
        for h, rel in SHA_LINE.findall(prov.read_text()):
            f = prov.parent / rel
            if f.is_file() and _sha(f) != h:
                wrong.append(f"{prov.relative_to(DATA)}: {rel}")
    assert wrong == []


def test_recorded_gz_are_reproducible() -> None:
    """Example 3: `gzip -n` (no name, mtime 0) so re-recording gives identical bytes."""
    bad = []
    for f in sorted(DATA.rglob("*.gz")):
        raw = f.read_bytes()
        if raw[:2] != b"\x1f\x8b" or raw[3] & 0x08 or raw[4:8] != b"\0\0\0\0":
            bad.append(str(f.relative_to(DATA)))
    assert bad == []


def test_r11_size_limits() -> None:
    """R11: no file > 5 MB; tests/data < 50 MB in total (example 5)."""
    files = [f for f in DATA.rglob("*") if f.is_file()]
    big = [str(f.relative_to(DATA)) for f in files if f.stat().st_size > MAX_FILE]
    assert big == []
    assert sum(f.stat().st_size for f in files) < MAX_TOTAL
