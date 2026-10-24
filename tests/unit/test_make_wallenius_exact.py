"""The exact Wallenius reference generator (docs/specs/M8-12b-golden-median-exact.md §3, B6, B7,
B10, examples 7-9; ADR-0044 decision 1).

Interface fixed here (the spec names the functions; these tests pin the details it leaves open):

- ``scripts/golden/make_wallenius_exact.py``, importable by path, with
  ``main(argv: list[str] | None = None) -> int`` and the CLI of spec §3.1
  (``--inputs DIR --out FILE --workers INT --categories A,B --log FILE``).
- ``r1_pmf(M, n, N, omega) -> np.ndarray`` (dtype ``np.longdouble``): ``pmf[x - lo]`` for
  ``x = lo..hi``, ``lo = max(0, N - (M - n))``, ``hi = min(n, N)``.
- ``r2_pmf(x, M, n, N, omega, *, p_start: int | None = None) -> Fraction``: R2 at one point. It
  calls the module-level ``r2_sum(x, M, n, N, omega, P) -> int`` (the fixed-point alternating
  sum S at precision P, spec §3.3 step 3; looked up as a module global, so a test can spy on it)
  starting at ``P = p_start`` (default ``264 + 2x``) and raising P until
  ``S > 0`` and ``bitlen(S) >= 200 + bitlen(x + 1) + 1``.
- ``exact_tails(k, M, n, N, omega) -> (over, under, method)``: ``over = P(X >= k)``,
  ``under = P(X <= k)`` (``Fraction`` or float; compared through ``Fraction(...)``) and
  ``method = {"over": str, "under": str}``, ``"R1"`` for an 80-bit urn tail and ``"R2"`` for the
  exact-integral small tail.
- With ``--categories`` the count-based checks of §3.4 (≥ 10 spot checks, ≥ 3 Σf categories,
  the vacuity floor) are not enforced on the subset; a subset run on valid inputs exits 0.
"""

from __future__ import annotations

import ast
import gzip
import importlib.util
import itertools
import json
import math
import subprocess
import sys
from collections.abc import Callable
from fractions import Fraction
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "golden" / "make_wallenius_exact.py"
GOLDEN_INPUTS = REPO / "tests" / "data" / "golden" / "r" / "inputs"
REL_R2 = Fraction(1, 2**190)


def _load() -> ModuleType:
    assert SCRIPT.is_file(), f"missing {SCRIPT.relative_to(REPO)}"
    spec = importlib.util.spec_from_file_location("make_wallenius_exact", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = mod  # multiprocessing pickles functions by module name
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def mod() -> ModuleType:
    return _load()


# --- the test's own exact reference: the urn recursion in rationals -----------------------------


def exact_pmf(M: int, n: int, N: int, omega: float) -> dict[int, Fraction]:
    """P(X = x) exactly, by the urn definition with the dyadic ω (independent of R1 and R2)."""
    w, m2 = Fraction(omega), M - n
    probs = {0: Fraction(1)}
    for t in range(N):
        new: dict[int, Fraction] = {}
        for r, p in probs.items():
            red, white = w * (n - r), Fraction(m2 - (t - r))
            if red > 0:
                new[r + 1] = new.get(r + 1, Fraction(0)) + p * red / (red + white)
            if white > 0:
                new[r] = new.get(r, Fraction(0)) + p * white / (red + white)
        probs = new
    return probs


def _rel(a: Fraction, b: Fraction) -> Fraction:
    return abs(a - b) / abs(b)


def _code(fn: Callable[[], Any]) -> Any:
    try:
        return fn()
    except SystemExit as exc:
        return exc.code


# --- examples 7-8 -------------------------------------------------------------------------------


def test_example_7_r2_pmf(mod: ModuleType) -> None:
    want = {0: Fraction(1, 15), 1: Fraction(3, 5), 2: Fraction(1, 3)}
    assert exact_pmf(4, 2, 2, 2.0) == want  # the urn: (4/6)(2/4) = 1/3, (2/6)(1/5) = 1/15
    for x, f in want.items():
        got = Fraction(mod.r2_pmf(x, 4, 2, 2, 2.0))
        assert _rel(got, f) < REL_R2, x
        assert float(got) == float(f), x


def test_example_7_exact_tails(mod: ModuleType) -> None:
    over, under, method = mod.exact_tails(2, 4, 2, 2, 2.0)
    assert float(over) == float(Fraction(1, 3))
    assert Fraction(under) == 1, "k = hi: under is exactly 1"
    assert method == {"over": "R1", "under": "R1"}, "both tails >= 1e-3 are R1"


def test_example_8_single_draw(mod: ModuleType) -> None:
    assert _rel(Fraction(mod.r2_pmf(1, 3, 1, 1, 2.0)), Fraction(1, 2)) < REL_R2


@pytest.mark.parametrize(("M", "n", "N"), [(5, 2, 3), (10, 4, 6), (20, 7, 9), (30, 25, 12)])
def test_omega_one_is_the_central_hypergeometric(mod: ModuleType, M: int, n: int, N: int) -> None:
    lo, hi = max(0, N - (M - n)), min(n, N)
    for x in range(lo, hi + 1):
        hyper = Fraction(math.comb(n, x) * math.comb(M - n, N - x), math.comb(M, N))
        assert _rel(Fraction(mod.r2_pmf(x, M, n, N, 1.0)), hyper) < REL_R2, x


@pytest.mark.parametrize(
    ("M", "n", "N", "omega"), [(8, 3, 4, 2.0), (12, 5, 7, 0.5), (15, 9, 6, 4.0), (20, 6, 15, 0.3)]
)
def test_r2_matches_the_exact_urn_and_sums_to_one(
    mod: ModuleType, M: int, n: int, N: int, omega: float
) -> None:
    exact = exact_pmf(M, n, N, omega)
    r2 = {x: Fraction(mod.r2_pmf(x, M, n, N, omega)) for x in exact}
    for x, f in exact.items():
        assert _rel(r2[x], f) < REL_R2, x
    assert abs(sum(r2.values()) - 1) < REL_R2, "Σf = 1"


@pytest.mark.parametrize("omega", [2.0, 4.0, 0.5])
def test_colour_symmetry_with_dyadic_inverse(mod: ModuleType, omega: float) -> None:
    M, n, N = 14, 5, 8
    for x in range(max(0, N - (M - n)), min(n, N) + 1):
        a = Fraction(mod.r2_pmf(x, M, n, N, omega))
        b = Fraction(mod.r2_pmf(N - x, M, M - n, N, 1 / omega))
        assert _rel(a, b) < REL_R2, x


# --- R1 ----------------------------------------------------------------------------------------


def test_r1_is_long_double_over_the_support(mod: ModuleType) -> None:
    M, n, N, omega = 10, 6, 7, 1.5  # lo = 3, hi = 6
    pmf = mod.r1_pmf(M, n, N, omega)
    assert pmf.dtype == np.longdouble
    assert len(pmf) == 4
    exact = exact_pmf(M, n, N, omega)
    for x in range(3, 7):
        assert abs(float(pmf[x - 3]) - float(exact[x])) <= 1e-15 * float(exact[x]), x


@settings(max_examples=25, deadline=None)
@given(
    M=st.integers(2, 60),
    data=st.data(),
    omega=st.floats(0.2, 5.0, allow_nan=False, allow_infinity=False),
)
def test_r1_and_r2_agree_on_small_urns(
    mod: ModuleType, M: int, data: st.DataObject, omega: float
) -> None:
    n = data.draw(st.integers(1, M - 1))
    N = data.draw(st.integers(1, min(M - 1, 25)))
    lo, hi = max(0, N - (M - n)), min(n, N)
    r1 = mod.r1_pmf(M, n, N, omega)
    assert len(r1) == hi - lo + 1
    for x in range(lo, hi + 1):
        r2 = float(mod.r2_pmf(x, M, n, N, omega))
        if r2 > 1e-300:
            assert abs(float(r1[x - lo]) - r2) <= 1e-14 * r2, (x, float(r1[x - lo]), r2)


# --- regimes and exact tails --------------------------------------------------------------------


def test_small_tail_is_r2_and_the_large_one_its_complement(mod: ModuleType) -> None:
    M, n, N, omega, k = 40, 10, 10, 0.25, 5
    exact = exact_pmf(M, n, N, omega)
    over_x = sum(v for x, v in exact.items() if x >= k)
    under_x = sum(v for x, v in exact.items() if x <= k)
    assert over_x < Fraction(1, 1000)  # precondition: over is the small tail
    over, under, method = mod.exact_tails(k, M, n, N, omega)
    assert method["over"] == "R2"
    assert set(method) == {"over", "under"}
    assert _rel(Fraction(over), over_x) < REL_R2
    assert _rel(Fraction(under), under_x) < REL_R2, "1 - small + pmf(k)"
    assert Fraction(over) + Fraction(under) >= 1


def test_small_lower_tail_is_r2(mod: ModuleType) -> None:
    M, n, N, omega, k = 40, 10, 10, 4.0, 0
    exact = exact_pmf(M, n, N, omega)
    assert exact[0] < Fraction(1, 1000)
    over, under, method = mod.exact_tails(k, M, n, N, omega)
    assert method["under"] == "R2"
    assert _rel(Fraction(under), exact[0]) < REL_R2
    assert Fraction(over) == 1, "k = lo: over is exactly 1"


@pytest.mark.parametrize("k", [1, 2, 3])
def test_tails_at_or_above_the_threshold_are_r1(mod: ModuleType, k: int) -> None:
    M, n, N, omega = 40, 10, 10, 4.0
    exact = exact_pmf(M, n, N, omega)
    over, under, method = mod.exact_tails(k, M, n, N, omega)
    assert method == {"over": "R1", "under": "R1"}
    for got, tail in ((over, [x for x in exact if x >= k]), (under, [x for x in exact if x <= k])):
        want = sum(exact[x] for x in tail)
        assert abs(float(got) - float(want)) <= 1e-15 * float(want)


# --- the precision rule (self-verifying, spec §3.3 step 3) --------------------------------------


def test_precision_rises_from_a_tiny_start(
    mod: ModuleType, monkeypatch: pytest.MonkeyPatch
) -> None:
    real = mod.r2_sum
    calls: list[tuple[int, int, int]] = []  # (x passed, P, S)

    def spy(*args: Any, **kwargs: Any) -> int:
        s = int(real(*args, **kwargs))
        p = kwargs.get("P", args[5] if len(args) > 5 else None)
        calls.append((int(args[0]), int(p), s))
        return s

    monkeypatch.setattr(mod, "r2_sum", spy)
    M, n, N, omega, x = 400, 20, 210, 3.0, 10
    got = Fraction(mod.r2_pmf(x, M, n, N, omega, p_start=4))
    ps = [p for _, p, _ in calls]
    assert ps[0] == 4
    assert len(ps) >= 2 and all(b > a for a, b in itertools.pairwise(ps)), ps
    x_used, _, s = calls[-1]
    assert s > 0 and s.bit_length() >= 200 + (x_used + 1).bit_length() + 1
    for x_used, _, s in calls[:-1]:
        assert s <= 0 or s.bit_length() < 200 + (x_used + 1).bit_length() + 1, "stopped late"
    calls.clear()
    default = Fraction(mod.r2_pmf(x, M, n, N, omega))
    assert calls[0][1] == 264 + 2 * calls[0][0], "default start P = 264 + 2x"
    assert _rel(got, default) < REL_R2


def test_precision_result_matches_the_exact_urn(mod: ModuleType) -> None:
    M, n, N, omega = 60, 12, 40, 3.0
    exact = exact_pmf(M, n, N, omega)
    for x in (0, 3, 12):
        assert _rel(Fraction(mod.r2_pmf(x, M, n, N, omega, p_start=2)), exact[x]) < REL_R2, x


# --- CLI on tiny inputs: layout, ω, error exits -------------------------------------------------


def _write_gz(path: Path, text: str) -> None:
    path.write_bytes(gzip.compress(text.encode(), mtime=0))


P_DETECT = [
    "0.900000", "0.800000", "0.750000", "0.700000", "0.650000", "0.600000",
    "0.400000", "0.350000", "0.300000", "0.250000", "0.200000", "0.123457",
]  # fmt: skip
DETECTED = [1, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0]
CATS = {"CAT_A": [0, 1, 2, 3], "CAT_B": [4, 6, 8, 9, 11], "CAT_C": [5, 7, 10]}


def _tiny_inputs(root: Path, cats: dict[str, list[int]] = CATS, p: list[str] = P_DETECT) -> Path:
    root.mkdir(parents=True, exist_ok=True)
    genes = "".join(
        f"GENE_{i + 1:04d}\t{q}\t{d}\n" for i, (q, d) in enumerate(zip(p, DETECTED, strict=True))
    )
    _write_gz(root / "ora_genes.tsv.gz", "gene\tp_detect\tdetected\n" + genes)
    body = "".join(f"{c}\tGENE_{g + 1:04d}\n" for c, members in cats.items() for g in members)
    _write_gz(root / "ora_categories.tsv.gz", "category\tgene\n" + body)
    return root


def _run(*args: str) -> subprocess.CompletedProcess[str]:
    assert SCRIPT.is_file(), f"missing {SCRIPT.relative_to(REPO)}"
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True, cwd=REPO, timeout=600
    )


def _rows(path: Path) -> list[list[str]]:
    text = gzip.decompress(path.read_bytes()).decode("utf-8")
    assert "\r" not in text and text.endswith("\n")
    return [line.split("\t") for line in text.splitlines()]


def _exact_omega(members: list[int]) -> float:
    p = [Fraction(q) for q in P_DETECT]
    s_in, s_all, n, M = sum(p[i] for i in members), sum(p), len(members), len(p)
    return float((s_in / n) / ((s_all - s_in) / (M - n)))


def _assert_gzip_house_format(raw: bytes) -> None:
    assert raw[:3] == b"\x1f\x8b\x08"
    assert raw[3] == 0, "no file name (FLG = 0)"
    assert raw[4:8] == b"\x00\x00\x00\x00", "mtime 0"
    assert raw[8] == 2, "compression level 9"


def test_cli_tiny_inputs_layout_and_values(tmp_path: Path) -> None:
    inputs = _tiny_inputs(tmp_path / "in")
    out, log = tmp_path / "exact.tsv.gz", tmp_path / "log.json"
    proc = _run(
        "--inputs", str(inputs), "--out", str(out), "--log", str(log),
        "--categories", "CAT_A,CAT_B,CAT_C",
    )  # fmt: skip
    assert proc.returncode == 0, proc.stderr
    _assert_gzip_house_format(out.read_bytes())
    rows = _rows(out)
    assert rows[0] == ["category", "tail", "p_exact"]
    assert [r[:2] for r in rows[1:]] == [[c, t] for c in CATS for t in ("over", "under")]
    M, N = len(P_DETECT), sum(DETECTED)
    for i, (cat, members) in enumerate(CATS.items()):
        k = sum(DETECTED[g] for g in members)
        exact = exact_pmf(M, len(members), N, _exact_omega(members))
        want = {
            "over": float(sum(v for x, v in exact.items() if x >= k)),
            "under": float(sum(v for x, v in exact.items() if x <= k)),
        }
        for row in rows[1 + 2 * i : 3 + 2 * i]:
            got = float(row[2])
            assert row[2] == format(got, ".17g"), "17 significant digits, round-trip exact"
            assert abs(got - want[row[1]]) <= 1e-15 * want[row[1]], (cat, row)
            assert 0 < got <= 1
        assert float(rows[1 + 2 * i][2]) + float(rows[2 + 2 * i][2]) >= 1
    assert isinstance(json.loads(log.read_text()), dict)


def test_cli_rows_do_not_depend_on_the_other_categories(tmp_path: Path) -> None:
    inputs = _tiny_inputs(tmp_path / "in")
    a, b = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    assert _run("--inputs", str(inputs), "--out", str(a), "--categories", "CAT_B").returncode == 0
    proc = _run("--inputs", str(inputs), "--out", str(b), "--categories", "CAT_A,CAT_B,CAT_C")
    assert proc.returncode == 0, proc.stderr
    assert _rows(a)[1:] == [r for r in _rows(b)[1:] if r[0] == "CAT_B"]


def test_cli_is_deterministic(tmp_path: Path) -> None:
    inputs = _tiny_inputs(tmp_path / "in")
    outs = [tmp_path / f"{i}.tsv.gz" for i in range(2)]
    for out, workers in zip(outs, ("1", "2"), strict=True):
        proc = _run("--inputs", str(inputs), "--out", str(out), "--workers", workers,
                    "--categories", "CAT_A,CAT_B,CAT_C")  # fmt: skip
        assert proc.returncode == 0, proc.stderr
    assert outs[0].read_bytes() == outs[1].read_bytes()


def test_cli_empty_categories_exit_without_writing(tmp_path: Path) -> None:
    inputs = _tiny_inputs(tmp_path / "in", cats={})
    out = tmp_path / "exact.tsv.gz"
    assert _run("--inputs", str(inputs), "--out", str(out)).returncode != 0
    assert not out.exists()
    assert not list(tmp_path.glob("exact*")), "no partial or temp file left"


def test_cli_undefined_omega_names_the_category(tmp_path: Path) -> None:
    p = ["0.500000"] * 4 + ["0.000000"] * 8  # mean_out = 0 for CAT_Z
    inputs = _tiny_inputs(tmp_path / "in", cats={"CAT_Z": [0, 1, 2, 3]}, p=p)
    out = tmp_path / "exact.tsv.gz"
    proc = _run("--inputs", str(inputs), "--out", str(out), "--categories", "CAT_Z")
    assert proc.returncode != 0
    assert "CAT_Z" in proc.stdout + proc.stderr
    assert not out.exists()


def test_refuses_without_80_bit_long_double(
    mod: ModuleType, monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    inputs = _tiny_inputs(tmp_path / "in")
    out = tmp_path / "exact.tsv.gz"
    real = np.finfo

    class Float64Info:
        nmant = 52
        bits = 64
        dtype = np.dtype(np.float64)

    def fake(t: Any) -> Any:
        return Float64Info() if np.dtype(t) == np.dtype(np.longdouble) else real(t)

    monkeypatch.setattr(np, "finfo", fake)
    argv = ["--inputs", str(inputs), "--out", str(out), "--categories", "CAT_A"]
    assert _code(lambda: mod.main(argv)) not in (0, None)
    assert not out.exists()


# --- B10: imports -------------------------------------------------------------------------------


def test_stdlib_and_numpy_only() -> None:
    assert SCRIPT.is_file(), f"missing {SCRIPT.relative_to(REPO)}"
    tree = ast.parse(SCRIPT.read_text())
    allowed = set(sys.stdlib_module_names) | {"numpy", "__future__"}
    for node in ast.walk(tree):
        names: list[str] = []
        if isinstance(node, ast.Import):
            names = [a.name for a in node.names]
        elif isinstance(node, ast.ImportFrom):
            assert node.level == 0, "no relative imports"
            names = [node.module or ""]
        elif isinstance(node, ast.Call) and getattr(node.func, "id", "") == "__import__":
            pytest.fail("no dynamic imports")
        for name in names:
            assert name.split(".")[0] in allowed, f"import of {name}"


# --- example 9 on the golden ORA inputs (slow: R2 at about 9 s per support point) ---------------


@pytest.fixture(scope="module")
def golden_runs(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    tmp = tmp_path_factory.mktemp("wallenius_golden")
    runs = {
        "w1": ("CAT_010", "1"),
        "w2": ("CAT_010", "2"),
        "pair": ("CAT_010,CAT_013", "2"),
    }
    out: dict[str, Path] = {}
    for key, (cats, workers) in runs.items():
        path = tmp / f"{key}.tsv.gz"
        proc = _run("--inputs", str(GOLDEN_INPUTS), "--out", str(path),
                    "--categories", cats, "--workers", workers)  # fmt: skip
        assert proc.returncode == 0, proc.stderr
        out[key] = path
    return out


@pytest.mark.slow
def test_cat_010_bytes_do_not_depend_on_workers(golden_runs: dict[str, Path]) -> None:
    assert golden_runs["w1"].read_bytes() == golden_runs["w2"].read_bytes()
    _assert_gzip_house_format(golden_runs["w1"].read_bytes())
    assert _rows(golden_runs["w1"]) == [
        ["category", "tail", "p_exact"],
        ["CAT_010", "over", "6.0681368749300294e-46"],
        ["CAT_010", "under", "1"],
    ]


@pytest.mark.slow
def test_adr_0044_spot_values(golden_runs: dict[str, Path]) -> None:
    rows = {(c, t): v for c, t, v in _rows(golden_runs["pair"])[1:]}
    assert list(rows) == [
        ("CAT_010", "over"),
        ("CAT_010", "under"),
        ("CAT_013", "over"),
        ("CAT_013", "under"),
    ]
    assert rows[("CAT_010", "over")] == "6.0681368749300294e-46"
    assert rows[("CAT_010", "under")] == "1"
    assert rows[("CAT_013", "under")] == "1.0597943893784124e-17"
    assert rows[("CAT_013", "over")] == "1"
