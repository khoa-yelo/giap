"""Median-exact golden GSEA inputs (docs/specs/M8-12b-golden-median-exact.md §1, B1-B3,
examples 1-5; M8-13 lead decision 1).

Interface fixed here (spec §1.3):

- ``scripts/golden/make_golden_inputs.py`` gains a pure function
  ``_balance_signs(micro: np.ndarray, frozen: np.ndarray, taken: set[int]) -> np.ndarray``,
  called by ``gsea_inputs`` (looked up as a module global, so a test can patch it) after
  ``_make_unique`` and before ``values = micro / SCALE``. ``micro`` is the int64 micro-unit
  array, ``frozen`` the boolean mask of the special TUs, ``taken`` the reserved values
  (``{0, 612345, -734567}`` in the generator). It returns a new array; it raises ``SystemExit``
  when ``N - Z`` is odd or when ``|d|`` exceeds the free TUs on the majority side.
- ``gsea_inputs`` raises ``SystemExit`` unless the median is exactly 0.0 with 997 / 6 / 997
  signs (post-conditions, example 5).

The committed-file tests are red until the lead commits the regenerated data (lead decision 3).
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import subprocess
import sys
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "golden" / "make_golden_inputs.py"
COMMITTED = REPO / "tests" / "data" / "golden" / "r" / "inputs"
FILES = ("gsea_stats.tsv.gz", "gsea_sets.tsv.gz", "ora_genes.tsv.gz", "ora_categories.tsv.gz")

# Spec §2.3: the ORA inputs keep their M8-12 bytes.
ORA_SHA256 = {
    "ora_genes.tsv.gz": "e94df51857577faf669c1e59e16454eb911ffdc7149d0c38b8956683254b1023",
    "ora_categories.tsv.gz": "114fcef698799d06250a8561e14022d9ea83af0cc6602193b48c69b35f9fef34",
}
# Spec §1.3 table at GOLDEN_SEED = 50000: TU -> (before, after) in micro-units.
FLIPPED = {
    "TU_0807": (-554, 555),  # +554 is taken: one step away from zero
    "TU_1362": (-1995, 1995),
    "TU_0809": (-5135, 5135),
    "TU_1686": (-5627, 5627),
    "TU_1539": (-5635, 5635),
    "TU_0883": (-5775, 5775),
    "TU_1339": (-7387, 7387),
    "TU_0108": (-8749, 8749),
}
# The special TUs of the M8-12 file (the special draw is unchanged, so they stay the same TUs).
ZERO_TUS = ["TU_0193", "TU_0446", "TU_1046", "TU_1530", "TU_1658", "TU_1672"]
TIE_IDS_SHA256 = {  # sha256 of the ","-joined TU IDs (file order) holding each tie value
    "0.612345": "74681ab0744de15fd386e4539a03b3834444eb2a386a94cc05523f2f81adff63",
    "-0.734567": "8bebecafc5738cf3150d8f17ed9622236a4b9928185437ba3c8aa6a9d203e759",
}
TAKEN = {0, 612345, -734567}


def _load() -> ModuleType:
    assert SCRIPT.is_file(), f"missing {SCRIPT.relative_to(REPO)}"
    spec = importlib.util.spec_from_file_location("make_golden_inputs_m8_12b", SCRIPT)
    assert spec is not None and spec.loader is not None
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _generate(out: Path, seed: int) -> None:
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--out", str(out), "--seed", str(seed)],
        capture_output=True,
        text=True,
        cwd=REPO,
    )
    assert proc.returncode == 0, proc.stderr


def _scores(path: Path) -> tuple[list[str], list[str]]:
    with gzip.open(path, "rt", newline="") as fh:
        rows = list(csv.reader(fh, delimiter="\t"))[1:]
    return [r[0] for r in rows], [r[1] for r in rows]


def _micro(strings: list[str]) -> np.ndarray:
    return np.array([round(float(s) * 10**6) for s in strings], dtype=np.int64)


def _check_median_exact(strings: list[str]) -> None:
    scores = np.array([float(s) for s in strings])
    assert float(np.median(scores)) == 0.0  # exact ==, spec §1.2
    assert int((scores < 0).sum()) == 997
    assert int((scores == 0).sum()) == 6
    assert int((scores > 0).sum()) == 997
    assert all(s != "-0.000000" for s in strings), "never -0.0"
    order = np.sort(scores, kind="stable")
    assert (order[997:1003] == 0.0).all(), "zeros at sorted positions 997-1002"
    assert order[996] < 0.0 < order[1003]


@pytest.fixture(scope="module")
def fresh(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("median_exact")
    _generate(out, 50000)
    return out


# --- B1: the median is exactly 0.0 (example 3) -----------------------------------------------


def test_golden_seed_is_50000(golden_seed: int) -> None:
    assert golden_seed == 50000  # the pinned tables below are for this seed


def test_median_is_exactly_zero_with_997_6_997_signs(fresh: Path) -> None:
    _check_median_exact(_scores(fresh / "gsea_stats.tsv.gz")[1])


def test_median_is_exactly_zero_for_another_seed(tmp_path: Path, golden_seed: int) -> None:
    _generate(tmp_path, golden_seed + 1)
    _check_median_exact(_scores(tmp_path / "gsea_stats.tsv.gz")[1])


def test_ties_and_zeros_are_kept_on_the_same_tus(fresh: Path) -> None:
    ids, strings = _scores(fresh / "gsea_stats.tsv.gz")
    assert [i for i, s in zip(ids, strings, strict=True) if s == "0.000000"] == ZERO_TUS
    for value, digest in TIE_IDS_SHA256.items():
        tied = [i for i, s in zip(ids, strings, strict=True) if s == value]
        assert len(tied) == 50, value
        assert hashlib.sha256(",".join(tied).encode()).hexdigest() == digest, value
    free = [s for s in strings if s not in ("0.000000", *TIE_IDS_SHA256)]
    assert len(free) == 1894 == len(set(free)), "free scores distinct, non-zero, not tie values"


def test_exactly_the_eight_pinned_rows_flip(monkeypatch: pytest.MonkeyPatch) -> None:
    mod = _load()
    real = mod._balance_signs
    seen: dict[str, Any] = {}

    def spy(micro: np.ndarray, frozen: np.ndarray, taken: set[int]) -> np.ndarray:
        seen["before"], seen["frozen"] = micro.copy(), frozen.copy()
        seen["after"] = real(micro, frozen, taken)
        return seen["after"]

    monkeypatch.setattr(mod, "_balance_signs", spy)
    rng = np.random.default_rng(np.random.SeedSequence(50000).spawn(2)[0])
    stats, _ = mod.gsea_inputs(rng)
    before, after = seen["before"], seen["after"]
    changed = np.flatnonzero(before != after)
    ids = [f"TU_{i + 1:04d}" for i in changed]
    assert sorted(ids) == sorted(FLIPPED)
    for i, tu in zip(changed, ids, strict=True):
        assert (int(before[i]), int(after[i])) == FLIPPED[tu], tu
    assert not seen["frozen"][changed].any(), "frozen TUs are never touched"
    written = _micro([s for _, s in stats])
    assert np.array_equal(written, after), "the written scores are the balanced ones"


def test_balancing_removed_fails_loudly(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Example 5: the mutation (no balancing) exits non-zero naming the failing condition."""
    mod = _load()
    assert hasattr(mod, "_balance_signs"), "_balance_signs is missing"
    monkeypatch.setattr(mod, "_balance_signs", lambda micro, frozen, taken: micro.copy())
    with pytest.raises(SystemExit) as exc:
        mod.main(["--out", str(tmp_path), "--seed", "50000"])
    msg = str(exc.value.code)
    assert "0.0" in msg or "997" in msg, msg
    assert not (tmp_path / "gsea_stats.tsv.gz").exists()


# --- B2: ORA bytes, determinism ----------------------------------------------------------------


@pytest.mark.parametrize("name", sorted(ORA_SHA256))
def test_ora_inputs_keep_their_pinned_bytes(fresh: Path, name: str) -> None:
    assert hashlib.sha256((fresh / name).read_bytes()).hexdigest() == ORA_SHA256[name]


def test_ora_bytes_do_not_depend_on_the_gsea_stream(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    mod = _load()
    assert hasattr(mod, "_balance_signs"), "_balance_signs is missing"

    def other_gsea(rng: np.random.Generator) -> tuple[list[tuple[str, ...]], ...]:
        rng.standard_normal(12345)  # a GSEA-stream change: different draws
        return [("TU_0001", "0.000000")], [("SET_001", "TU_0001")]

    monkeypatch.setattr(mod, "gsea_inputs", other_gsea)
    assert mod.main(["--out", str(tmp_path), "--seed", "50000"]) == 0
    for name, digest in ORA_SHA256.items():
        assert hashlib.sha256((tmp_path / name).read_bytes()).hexdigest() == digest, name


def test_two_runs_are_byte_identical(fresh: Path, tmp_path: Path) -> None:
    _generate(tmp_path, 50000)
    _check_median_exact(_scores(fresh / "gsea_stats.tsv.gz")[1])
    for name in FILES:
        assert (tmp_path / name).read_bytes() == (fresh / name).read_bytes(), name


# --- B3: the committed inputs (red until the data commit) --------------------------------------


def test_committed_scores_have_an_exactly_zero_median() -> None:
    _check_median_exact(_scores(COMMITTED / "gsea_stats.tsv.gz")[1])


def test_committed_inputs_equal_a_fresh_median_exact_generation(fresh: Path) -> None:
    _check_median_exact(_scores(fresh / "gsea_stats.tsv.gz")[1])
    for name in FILES:
        assert (COMMITTED / name).read_bytes() == (fresh / name).read_bytes(), name


# --- _balance_signs unit cases (examples 1-2, edge cases) --------------------------------------


def _balance(micro: list[int], frozen: list[bool], taken: set[int] | None = None) -> np.ndarray:
    fn = _load()._balance_signs
    out = fn(np.array(micro, dtype=np.int64), np.array(frozen, dtype=bool), taken or {0})
    return np.asarray(out)


def test_example_1_flip_with_a_collision_step() -> None:
    # free [-6, -3, -2, -1, 1, 5] and two frozen zeros: d = 1, -1 -> +1 is taken -> +2.
    micro = [-6, -3, -2, -1, 1, 5, 0, 0]
    frozen = [False] * 6 + [True, True]
    out = _balance(micro, frozen)
    assert out.tolist() == [-6, -3, -2, 2, 1, 5, 0, 0]
    assert sorted(out.tolist()) == [-6, -3, -2, 0, 0, 1, 2, 5]
    assert float(np.median(out / 10**6)) == 0.0


def test_input_array_is_not_modified() -> None:
    micro = np.array([-6, -3, -2, -1, 1, 5, 0, 0], dtype=np.int64)
    frozen = np.array([False] * 6 + [True, True])
    _load()._balance_signs(micro, frozen, {0})
    assert micro.tolist() == [-6, -3, -2, -1, 1, 5, 0, 0]


def test_collision_steps_past_frozen_and_free_values() -> None:
    # -5 -> 5 is a frozen tie value, 6 is a free value, so it lands on 7.
    micro = [-5, -8, -9, -10, 5, 6, 0, 0]
    frozen = [False, False, False, False, True, False, True, True]
    out = _balance(micro, frozen, {0, 5})
    assert out.tolist() == [7, -8, -9, -10, 5, 6, 0, 0]


def test_positive_majority_flips_the_smallest_positive_away_from_zero() -> None:
    micro = [3, 8, 9, -3, 10, 0, 0, -20]
    frozen = [False] * 5 + [True, True, False]
    out = _balance(micro, frozen)
    assert out.tolist() == [-4, 8, 9, -3, 10, 0, 0, -20]


def test_ties_in_magnitude_are_broken_by_index() -> None:
    # d = 1; candidates -2 (index 3) and -2 (index 1) tie on |micro|: the lower index flips.
    micro = [-9, -2, 4, -2, -7, 0, 0, 0, 0, 11]
    frozen = [False] * 5 + [True] * 4 + [False]
    out = _balance(micro, frozen)
    assert out.tolist() == [-9, 2, 4, -2, -7, 0, 0, 0, 0, 11]


def test_balanced_input_is_unchanged() -> None:
    micro = [-4, -1, 0, 0, 2, 3]
    frozen = [False, False, True, True, False, False]
    assert _balance(micro, frozen).tolist() == micro


def test_example_2_odd_free_count_exits() -> None:
    with pytest.raises(SystemExit):
        _balance([-2, -1, 0, 1, 2, 3], [False, False, True, False, False, False])


def test_too_few_free_values_on_the_majority_side_exits() -> None:
    # Three frozen negatives, one free positive: d = 1 but no free negative to flip.
    with pytest.raises(SystemExit):
        _balance([-7, -7, -7, 1, 0, 0], [True, True, True, False, True, True], {0, -7})


@settings(max_examples=60, deadline=None)
@given(
    free=st.lists(
        st.integers(-400, 400).filter(lambda v: v != 0), min_size=2, max_size=40, unique=True
    ),
    n_zero=st.integers(2, 6),
    n_tie=st.integers(0, 3),
)
def test_balance_property(free: list[int], n_zero: int, n_tie: int) -> None:
    """Any feasible input: signs balance, frozen untouched, |d| sign flips, values stay legal."""
    tie = -1000  # a frozen negative tie value, never a free value
    if (len(free) + n_tie) % 2:  # N - Z even (example 2 covers the odd case)
        free = free[:-1]
    micro = [*free, *([tie] * n_tie), *([0] * n_zero)]
    frozen = [False] * len(free) + [True] * (n_tie + n_zero)
    arr = np.array(micro)
    d = (int((arr < 0).sum()) - int((arr > 0).sum())) // 2
    side = [v for v in free if (v < 0 if d > 0 else v > 0)]
    if abs(d) > len(side):
        return
    out = _balance(micro, frozen, {0, tie}).tolist()
    assert out[len(free) :] == micro[len(free) :], "frozen untouched"
    new_free = out[: len(free)]
    assert sum(v < 0 for v in out) == sum(v > 0 for v in out)
    assert len(set(new_free)) == len(new_free)
    assert all(v not in (0, tie) for v in new_free)
    flips = [(a, b) for a, b in zip(free, new_free, strict=True) if a != b]
    assert len(flips) == abs(d)
    assert all((a < 0) != (b < 0) and abs(b) >= abs(a) for a, b in flips)
    flipped = sorted(abs(a) for a, _ in flips)
    kept = [abs(v) for v, w in zip(free, new_free, strict=True) if v == w and v in side]
    assert not kept or not flipped or max(flipped) <= min(kept), "the closest to zero flip"
    assert float(np.median(np.array(out) / 10**6)) == 0.0
