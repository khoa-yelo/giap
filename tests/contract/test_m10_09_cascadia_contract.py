"""Contract test: Cascadia 0.0.7 ``.ssl`` sample through ``convert_cascadia`` (M10-09; §5.4, §11.4).

Spec: docs/specs/M10-09-cascadia-convert.md, "The contract test" (checks C1-C8), criteria 10
and 12; facts W1-W7 from the pinned 0.0.7 source.

Data live on FarmShare only (Khoa, 2026-10-23: recorded samples and converted tables are never
committed). The tests read

* ``SAMPLE = $GIPPA_ROOT/pilot/recorded/cascadia/0.0.7/WT_1.ssl.gz`` (the header plus the first
  1,000 data rows of ``WT_1.ssl``, M10-08 deliverable 13), and, for the ``slow`` variant,
* ``FULL = $GIPPA_ROOT/pilot/A3_cascadia/WT_1/WT_1.ssl``.

They **skip, never fail**, when ``GIPPA_ROOT`` is unset (CI) or the file is missing, and read no
environment at import time. The raw file is read with the stdlib ``csv`` module, independently of
the converter. No golden output file is read or written. Assertion messages give counts, row
positions and modification tokens only, never a row: criterion 12 pastes the output into the PR.
The count of scores in (0.8, 0.8000005) (C4) is printed and attached with ``record_property``;
run with ``-s`` or ``-rA`` to see it.
"""

from __future__ import annotations

import csv
import gzip
import io
import math
import os
import re
import subprocess
from collections import Counter
from functools import cache
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
W1 = [
    "file",
    "scan",
    "charge",
    "sequence",
    "score-type",
    "score",
    "retention-time",
    "start-time",
    "end-time",
]
W5_TOKENS = frozenset(
    {"[+15.99]", "[+0.98]", "[+57.02]", "[+42.01]", "[+43.01]", "[-17.03]", "[+25.980265]"}
)
CUTOFF = 0.8
NO_ROOT = "GIPPA_ROOT is not set: the Cascadia sample lives on FarmShare only (Khoa 2026-10-23)"
MISSING = "M10-08 smoke run not done or verdict fail"
KINDS = [
    pytest.param("sample", id="sample"),
    pytest.param("full", id="full", marks=pytest.mark.slow),
]


# ------------------------------------------------------------------ data (skip when absent)


def _root() -> Path:
    root = os.environ.get("GIPPA_ROOT")
    if not root:
        pytest.skip(NO_ROOT)
    return Path(root)


def _git_ok(*args: str) -> bool:
    res = subprocess.run(["git", "-C", str(REPO), *args], capture_output=True, check=False)
    return res.returncode == 0


def _outside_checkout(path: Path) -> Path:
    """``path`` is outside the checkout, or git-ignored and untracked in it (Khoa's rule)."""
    resolved, repo = path.resolve(), REPO.resolve()
    if repo not in resolved.parents:
        return path
    rel = str(resolved.relative_to(repo))
    ignored = _git_ok("check-ignore", "-q", "--", rel)
    tracked = _git_ok("ls-files", "--error-unmatch", "--", rel)
    assert ignored and not tracked, (
        f"{path} lies inside the checkout and is not git-ignored and untracked (Khoa's rule)"
    )
    return path


def data_path(kind: str) -> Path:
    root = _root()
    if kind == "sample":
        p = root / "pilot" / "recorded" / "cascadia" / "0.0.7" / "WT_1.ssl.gz"
    else:
        p = root / "pilot" / "A3_cascadia" / "WT_1" / "WT_1.ssl"
    if not p.is_file():
        pytest.skip(f"{p} is missing: {MISSING}")
    return _outside_checkout(p)


@cache
def _raw(path: Path) -> tuple[tuple[str, ...], tuple[dict[str, str], ...]]:
    """Header and rows by the stdlib ``csv`` module (tab, no quoting); gzip by magic bytes."""
    data = path.read_bytes()
    if data[:2] == b"\x1f\x8b":
        data = gzip.decompress(data)
    text = data.decode("utf-8")
    reader = csv.reader(io.StringIO(text, newline=""), delimiter="\t", quoting=csv.QUOTE_NONE)
    header = next(reader)
    rows = []
    for lineno, cells in enumerate(reader, start=2):
        assert len(cells) == len(header), f"line {lineno}: {len(cells)} fields"
        rows.append(dict(zip(header, cells, strict=True)))
    return tuple(header), tuple(rows)


def raw(kind: str) -> tuple[tuple[str, ...], tuple[dict[str, str], ...]]:
    return _raw(data_path(kind))


def _stem(file: str) -> str:
    """D1: after the last ``/``, drop one ``.gz``, then one ``.mzML`` (case-insensitive)."""
    name = file.rsplit("/", 1)[-1]
    name = re.sub(r"\.gz$", "", name, flags=re.IGNORECASE)
    return re.sub(r"\.mzml$", "", name, flags=re.IGNORECASE)


def _positions(bad: list[int]) -> str:
    return f"{len(bad)} rows, e.g. data rows {bad[:5]}"


# ------------------------------------------------------------------ C1-C7 on the raw file


@pytest.mark.parametrize("kind", KINDS)
def test_c1_header(kind: str) -> None:
    header, _ = raw(kind)
    assert list(header) == W1, f"header differs from W1: {list(header)}"


@pytest.mark.parametrize("kind", KINDS)
def test_c2_score_type(kind: str) -> None:
    _, rows = raw(kind)
    assert rows, "no data rows"
    bad = [i for i, r in enumerate(rows, 1) if r["score-type"] != "UNKNOWN"]
    assert not bad, f"score-type is not UNKNOWN: {_positions(bad)}"


def _score_ok(cell: str) -> bool:
    try:
        s = float(cell)
    except ValueError:
        return False
    return math.isfinite(s) and CUTOFF < s <= 1.0


@pytest.mark.parametrize("kind", KINDS)
def test_c3_score(kind: str) -> None:
    _, rows = raw(kind)
    bad = [i for i, r in enumerate(rows, 1) if not _score_ok(r["score"])]
    assert not bad, (
        f"score not a finite float in (0.8, 1] (threshold changed? M10-08 finding): "
        f"{_positions(bad)}"
    )


@pytest.mark.parametrize("kind", KINDS)
def test_c4_scan(kind: str, record_property: pytest.RecordProperty) -> None:
    _, rows = raw(kind)
    bad = [i for i, r in enumerate(rows, 1) if not re.fullmatch(r"[0-9]+", r["scan"])]
    assert not bad, f"scan is not a non-negative integer: {_positions(bad)}"
    scans = [int(r["scan"]) for r in rows]
    not_increasing = [i for i in range(2, len(scans) + 1) if scans[i - 1] <= scans[i - 2]]
    assert not not_increasing, f"scan not strictly increasing: {_positions(not_increasing)}"
    boundary = sum(1 for r in rows if _score_ok(r["score"]) and float(r["score"]) < 0.8000005)
    record_property("scores_in_0.8_0.8000005", boundary)
    print(f"C4 {kind}: {boundary} of {len(rows)} scores in (0.8, 0.8000005)")


@pytest.mark.parametrize("kind", KINDS)
def test_c5_file(kind: str) -> None:
    _, rows = raw(kind)
    files = {r["file"] for r in rows}
    assert len(files) == 1, f"{len(files)} distinct file values"
    assert _stem(next(iter(files))) == "WT_1"


@pytest.mark.parametrize("kind", KINDS)
def test_c6_sequence(kind: str) -> None:
    from gippa.io.converters.strip import strip_modifications

    _, rows = raw(kind)
    not_plain, with_i, bad_dash = [], [], []
    tokens: Counter[str] = Counter()
    for i, r in enumerate(rows, 1):
        seq = r["sequence"]
        pep = strip_modifications(seq)
        if not re.fullmatch(r"[A-Z]+", pep):
            not_plain.append(i)
        if "I" in pep:
            with_i.append(i)
        tokens.update(re.findall(r"\[[^\]]*\]", seq))
        lead = re.match(r"\[[^\]]*\]-", seq)
        rest = seq[lead.end() :] if lead else seq
        if "-" in re.sub(r"\[[^\]]*\]", "", rest):
            bad_dash.append(i)
    assert not not_plain, f"stripped sequence is not [A-Z]+: {_positions(not_plain)}"
    assert not with_i, (
        f"an I in 0.0.7 output (tokenizer I->L setting changed, or not 0.0.7?): "
        f"{_positions(with_i)}"
    )
    unknown = sorted(set(tokens) - W5_TOKENS)
    assert not unknown, (
        f"modification tokens outside the W5 mass-shift set: {unknown}; name tokens mean "
        "Cascadia ran without M10-08's offline Unimod launcher (its edge case 'modified "
        "sequences without mass shifts')"
    )
    assert not bad_dash, f"'-' not right after a leading N-terminal token: {_positions(bad_dash)}"


@pytest.mark.parametrize("kind", KINDS)
def test_c7_charge(kind: str) -> None:
    _, rows = raw(kind)
    bad = [i for i, r in enumerate(rows, 1) if not re.fullmatch(r"[1-9][0-9]*", r["charge"])]
    assert not bad, f"charge is not an integer >= 1: {_positions(bad)}"


# ------------------------------------------------------------------ C8 the conversion


@pytest.mark.parametrize("kind", KINDS)
def test_c8_conversion(kind: str) -> None:
    path = data_path(kind)  # skip first: no import before the data check (criterion 10)
    _, rows = raw(kind)
    from gippa.config.model import InputConfig
    from gippa.io.converters.cascadia import convert_cascadia

    res = convert_cascadia([path], search_db=None, cfg=InputConfig())
    errs = sorted({i.code for i in res.issues if i.severity == "error"})
    assert not errs, f"error issues: {errs}"
    df = res.peptides
    assert len(df) == len(rows), f"{len(df)} output rows for {len(rows)} input rows"
    got = Counter(
        zip(df["run"], df["spectrum_id"], df["modified_peptide"], df["charge"].astype(int),
            df["score"], strict=True)
    )  # fmt: skip
    want = Counter(
        (_stem(r["file"]), f"cascadia:{r['scan']}", r["sequence"], int(r["charge"]),
         float(r["score"]))
        for r in rows
    )  # fmt: skip
    assert got == want, f"{sum((got - want).values())} rows differ from the raw rows"
    assert set(df["source"]) == {"denovo"}
    assert "q_value" not in df.columns
    assert "pep" not in df.columns


@pytest.mark.parametrize("kind", KINDS)
def test_c8_cli_is_deterministic_and_valid(kind: str, tmp_path: Path) -> None:
    from click.testing import CliRunner

    from gippa.cli import cli
    from gippa.config.model import InputConfig
    from gippa.exitcodes import ExitCode
    from gippa.io.peptides import read_peptide_table, validate_peptides

    path = data_path(kind)
    (tmp_path / "t").mkdir()
    blobs = []
    for i in range(2):
        out = tmp_path / "t" / f"a3_{i}.tsv.gz"
        argv = ["convert", "--from", "cascadia", str(path), "--out", str(out)]
        result = CliRunner().invoke(cli, argv)
        assert result.exit_code == ExitCode.SUCCESS, f"exit {result.exit_code}"
        blobs.append(out.read_bytes())
    assert blobs[0] == blobs[1], "two CLI runs differ"
    df, found = read_peptide_table(tmp_path / "t" / "a3_0.tsv.gz")
    _, found2 = validate_peptides(df, InputConfig())
    errs = sorted({i.code for i in [*found, *found2] if i.severity == "error"})
    assert not errs, f"error issues on the written table: {errs}"
