"""The committed exact Wallenius reference, ``tests/data/golden/wallenius/`` (docs/specs/
M8-12b-golden-median-exact.md B6, B8, example 9, §4.2; ADR-0044 decisions 1 and 4).

Layout (the M8-13 contract, ``tests/helpers/golden_r_gate.py::read_exact_tails`` on
``m8/golden-equivalence``): a gzip TSV with header ``category tail p_exact``, 300 rows,
150 categories in ``ora_categories`` order, ``over`` before ``under``, 17-digit values.

These tests fail until the lead commits the output of ``scripts/slurm/golden_wallenius.sbatch``
(lead decision 3).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import math
import re
from fractions import Fraction
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
ROOT = REPO / "tests" / "data" / "golden" / "wallenius"
EXACT = ROOT / "exact_tails.tsv.gz"
PROVENANCE = ROOT / "PROVENANCE.md"
INPUTS = REPO / "tests" / "data" / "golden" / "r" / "inputs"
ORA = ("ora_genes.tsv.gz", "ora_categories.tsv.gz")
FIVE_MB = 5 * 1024 * 1024
SPOT = {  # spec example 9 (ADR-0044)
    ("CAT_010", "over"): "6.0681368749300294e-46",
    ("CAT_010", "under"): "1",
    ("CAT_012", "over"): "2.2283895161517117e-10",
    ("CAT_013", "under"): "1.0597943893784124e-17",
    ("CAT_004", "over"): "4.1695578881156134e-07",
}


def _read(path: Path) -> list[list[str]]:
    assert path.is_file(), f"missing {path.relative_to(REPO)}"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as fh:
        return list(csv.reader(fh, delimiter="\t"))


@pytest.fixture(scope="module")
def rows() -> list[list[str]]:
    return _read(EXACT)


@pytest.fixture(scope="module")
def tails(rows: list[list[str]]) -> dict[tuple[str, str], float]:
    """The ``read_exact_tails`` contract of M8-13, restated (that helper is not on this branch)."""
    header, body = rows[0], rows[1:]
    assert len(header) == 3 and header[:2] == ["category", "tail"], header
    out: dict[tuple[str, str], float] = {}
    for cat, tail, value in body:
        assert tail in ("over", "under"), tail
        assert (cat, tail) not in out, (cat, tail)
        out[(cat, tail)] = float(value)
    return out


def _categories() -> dict[str, list[str]]:
    cats: dict[str, list[str]] = {}
    for cat, gene in _read(INPUTS / "ora_categories.tsv.gz")[1:]:
        cats.setdefault(cat, []).append(gene)
    return cats


# --- B6: the file ------------------------------------------------------------------------------


def test_file_is_a_small_house_gzip() -> None:
    assert EXACT.is_file(), f"missing {EXACT.relative_to(REPO)}"
    raw = EXACT.read_bytes()
    assert len(raw) <= FIVE_MB
    assert raw[:3] == b"\x1f\x8b\x08" and raw[3] == 0 and raw[4:8] == b"\0\0\0\0"
    text = gzip.decompress(raw).decode("utf-8")
    assert "\r" not in text and text.endswith("\n")
    assert len(text.splitlines()) == 301


def test_rows_are_150_categories_in_file_order(rows: list[list[str]]) -> None:
    assert rows[0] == ["category", "tail", "p_exact"]
    want = [[c, t] for c in _categories() for t in ("over", "under")]
    assert len(want) == 300
    assert [r[:2] for r in rows[1:]] == want


def test_values_are_17_digit_probabilities(
    rows: list[list[str]], tails: dict[tuple[str, str], float]
) -> None:
    for _, _, value in rows[1:]:
        p = float(value)
        assert value == format(p, ".17g"), value
        assert 1e-300 <= p <= 1.0, value
    for cat in {c for c, _ in tails}:
        assert tails[(cat, "over")] + tails[(cat, "under")] >= 1.0, cat


def test_adr_0044_spot_values(rows: list[list[str]]) -> None:
    by_key = {(c, t): v for c, t, v in rows[1:]}
    for key, value in SPOT.items():
        assert by_key[key] == value, key


def test_vacuity_floor(tails: dict[tuple[str, str], float]) -> None:
    """ADR-0044 decision 4: >= 2 tails in (0, 1e-20) and >= 1 in [1e-20, 1e-8)."""
    values = list(tails.values())
    assert sum(0 < p < 1e-20 for p in values) >= 2
    assert sum(1e-20 <= p < 1e-8 for p in values) >= 1


def test_exactly_six_tails_below_1e_3(tails: dict[tuple[str, str], float]) -> None:
    small = sorted(k for k, p in tails.items() if p < 1e-3)
    assert small == [
        ("CAT_003", "over"),
        ("CAT_004", "over"),
        ("CAT_010", "over"),
        ("CAT_011", "over"),
        ("CAT_012", "over"),
        ("CAT_013", "under"),
    ]


# --- B8: PROVENANCE ----------------------------------------------------------------------------


def test_provenance_records_the_run() -> None:
    assert PROVENANCE.is_file(), f"missing {PROVENANCE.relative_to(REPO)}"
    text = PROVENANCE.read_text()
    assert "golden_wallenius.sbatch" in text and "make_wallenius_exact.py" in text, "command"
    assert re.search(r"[Jj]ob\D{0,20}\d{5,}", text), "Slurm job ID"
    assert re.search(r"commit\W{0,3}[0-9a-f]{40}\b", text), "git commit of the script"
    assert "(scripts modified after this commit)" not in text, "generated from committed code"
    assert re.search(r"\b20\d\d-\d\d-\d\d\b", text), "date"
    assert re.search(r"nmant\W{0,3}63", text), "the long double description"
    assert re.search(r"numpy\W{0,3}\d+\.\d+", text) and re.search(r"Python\W{0,3}3\.\d+", text)
    assert "ADR-0044" in text
    assert re.search(r"byte-identical", text), "two runs in the job"
    assert re.search(r"[Rr]eason", text), "regeneration only with a reason"


def test_provenance_ties_the_reference_to_the_committed_ora_bytes() -> None:
    text = PROVENANCE.read_text() if PROVENANCE.is_file() else ""
    for name in ORA:
        digest = hashlib.sha256((INPUTS / name).read_bytes()).hexdigest()
        assert digest in text, f"{name} sha256 {digest} not in PROVENANCE.md"


def test_provenance_summarises_the_checks() -> None:
    text = PROVENANCE.read_text() if PROVENANCE.is_file() else ""
    named = set(re.findall(r"CAT_\d{3}\W{1,3}(?:over|under)", text))
    assert len(named) >= 10, "the 10 spot-check tails are named"
    assert re.search(r"\d(\.\d+)?e-\d+", text), "the max relative differences are recorded"


# --- calibration: an independent float64 urn recursion agrees on all 300 tails ------------------


def _omega(p_all: list[Fraction], members: list[int]) -> float:
    s_in, s_all, n, M = sum(p_all[i] for i in members), sum(p_all), len(members), len(p_all)
    return float((s_in / n) / ((s_all - s_in) / (M - n)))


def _urn_pmf(M: int, n: int, N: int, omega: float) -> np.ndarray:
    m2 = M - n
    probs = np.ones(1)
    for t in range(N):
        r = np.arange(probs.size, dtype=np.float64)
        red = omega * (n - r)
        white = np.maximum(m2 - (t - r), 0.0)
        den = red + white
        p_red = np.divide(red, den, out=np.zeros_like(den), where=den > 0)
        p_white = np.divide(white, den, out=np.zeros_like(den), where=den > 0)
        new = np.zeros(min(probs.size + 1, n + 1))
        new[: probs.size] += probs * p_white
        new[1 : probs.size + 1] += (probs * p_red)[: new.size - 1]
        probs = new
    return probs


@pytest.mark.slow
def test_reference_agrees_with_a_float64_urn_on_every_tail(
    tails: dict[tuple[str, str], float],
) -> None:
    """§3.4: R1 80-bit vs float64 within 1e-13 (ADR-0044 measured 3.8e-15)."""
    genes = _read(INPUTS / "ora_genes.tsv.gz")[1:]
    index = {g: i for i, (g, _, _) in enumerate(genes)}
    p_all = [Fraction(p) for _, p, _ in genes]
    det = [int(d) for _, _, d in genes]
    M, N = len(genes), sum(det)
    worst = 0.0
    for cat, names in _categories().items():
        members = [index[g] for g in names]
        k = sum(det[i] for i in members)
        pmf = _urn_pmf(M, len(members), N, _omega(p_all, members))
        over, under = math.fsum(pmf[k:].tolist()), math.fsum(pmf[: k + 1].tolist())
        for tail, mine in (("over", over), ("under", under)):
            ref = tails[(cat, tail)]
            worst = max(worst, abs(mine - ref) / ref)
    assert worst <= 1e-13, worst
