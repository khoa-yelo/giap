"""Unit tests of the calibration harness ``tests/calibration/_criteria.py`` (issue M4-05).

GIAP_SPEC.md §8.4 criteria, with k = ``se_multiplier`` and R the number of independent replicates:

- level:   rate of p <= alpha  <=  alpha + k * sqrt(alpha (1 - alpha) / R)
- uniform: KS p-value against U(0, 1)  >  ``ks_min_p`` (continuous p only; > 5 % ties refused)
- FDR:     mean FDP  <=  q + k * sd(FDP, ddof=1) / sqrt(R)
- FWER:    share of replicates with >= 1 false discovery <= alpha + k * sqrt(alpha (1 - alpha) / R)

Seeds: ``suite_seed(suite, rep)`` = ``SeedSequence([calibration_base, crc32(suite), rep])`` drawn
as one uint64 (the ``giap.rng.for_step`` key [seed, crc32(name)] extended by the replicate index);
replicate ``rep`` uses ``default_rng(suite_seed(suite, rep))`` so any shard reproduces it alone.

Nightly job array (§8.4, §12 ``scripts/slurm/nightly.sbatch``): one suite per task, task i = the
i-th key of ``replicates_full`` in file order; ``--shard K/N`` splits a suite's replicates into N
contiguous blocks. Results are JSON lines in ``local/logs/calibration/<suite>.jsonl``
(``GIAP_CALIBRATION_LOG_DIR`` overrides the directory).

Unmarked tests run on every commit; the operating-characteristic tests at the end are marked
``calibration``.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import math
import zlib
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import _criteria as C
import numpy as np
import pytest
import yaml
from hypothesis import given, settings
from hypothesis import strategies as st
from scipy import stats as sps

HERE = Path(__file__).resolve().parent
TOL_PATH = HERE / "tolerances.yaml"
FAST = settings(max_examples=200, deadline=None)


def _raw() -> dict[str, object]:
    with TOL_PATH.open() as fh:
        data: dict[str, object] = yaml.safe_load(fh)
    return data


# --------------------------------------------------------------------------- tolerances


def test_load_tolerances_matches_protected_file() -> None:
    raw = _raw()
    tol = C.load_tolerances()
    assert isinstance(tol, C.Tolerances)
    assert tuple(tol.levels) == tuple(float(x) for x in raw["levels"])  # type: ignore[union-attr]
    assert tol.se_multiplier == float(raw["se_multiplier"])  # type: ignore[arg-type]
    assert tol.ks_min_p == float(raw["ks_min_p"])  # type: ignore[arg-type]
    assert tol.calibration_base == int(raw["seeds"]["calibration_base"])  # type: ignore[index]
    assert dict(tol.replicates_fast) == raw["replicates_fast"]
    assert dict(tol.replicates_full) == raw["replicates_full"]


def test_tolerances_are_frozen() -> None:
    tol = C.load_tolerances()
    with pytest.raises(dataclasses.FrozenInstanceError):
        tol.se_multiplier = 3.0  # type: ignore[misc]


def test_load_tolerances_reads_the_given_file(tmp_path: Path) -> None:
    raw = _raw()
    raw["se_multiplier"] = 3.0
    raw["seeds"] = {"GOLDEN_SEED": 1, "calibration_base": 7}  # type: ignore[assignment]
    raw["replicates_fast"] = {"regions": 11}
    raw["replicates_full"] = {"regions": 22}
    p = tmp_path / "tol.yaml"
    p.write_text(yaml.safe_dump(raw))
    tol = C.load_tolerances(p)
    assert tol.se_multiplier == 3.0
    assert tol.calibration_base == 7
    assert dict(tol.replicates_fast) == {"regions": 11}
    assert dict(tol.replicates_full) == {"regions": 22}


def test_harness_matches_existing_calibration_settings(tolerances: dict[str, object]) -> None:
    """The existing suites read the same values through the session ``tolerances`` fixture."""
    fast, full = tolerances["replicates_fast"], tolerances["replicates_full"]
    assert isinstance(fast, dict) and isinstance(full, dict)
    for suite in full:
        assert C.replicates(suite, full=False) == int(fast[suite])
        assert C.replicates(suite, full=True) == int(full[suite])
    k = float(tolerances["se_multiplier"])  # type: ignore[arg-type]
    for a in tolerances["levels"]:  # type: ignore[union-attr]
        for r in (int(fast["regions"]), int(full["regions"])):
            # inline formula in test_stats_permutation_calibration.py::_check_level
            assert C.level_bound(float(a), r) == pytest.approx(
                float(a) + k * np.sqrt(float(a) * (1 - float(a)) / r), rel=1e-12
            )


# --------------------------------------------------------------------------- examples 1-8


def test_example1_level_bound_005_100() -> None:
    assert C.level_bound(0.05, 100) == pytest.approx(0.05 + 2 * math.sqrt(0.0475 / 100), abs=1e-12)
    assert C.level_bound(0.05, 100) == pytest.approx(0.093589, abs=5e-7)


def test_example2_level_bound_001_1000() -> None:
    assert C.level_bound(0.01, 1000) == pytest.approx(
        0.01 + 2 * math.sqrt(0.0099 / 1000), abs=1e-12
    )
    assert C.level_bound(0.01, 1000) == pytest.approx(0.016293, abs=5e-7)


def test_example3_check_level_fails() -> None:
    res = C.check_level([0.001] * 3 + [0.5] * 97, 0.01)
    assert isinstance(res, C.CheckResult)
    assert res.ok is False
    assert res.observed == pytest.approx(0.03)
    assert res.bound == pytest.approx(0.0299, abs=1e-4)
    assert res.bound == pytest.approx(C.level_bound(0.01, 100), rel=1e-12)
    assert res.R == 100
    assert isinstance(res.detail, str) and res.detail


def test_example4_check_level_passes() -> None:
    res = C.check_level([0.001] * 2 + [0.5] * 98, 0.01)
    assert res.ok is True
    assert res.observed == pytest.approx(0.02)
    assert res.bound == pytest.approx(0.0299, abs=1e-4)
    assert res.R == 100


def test_check_level_counts_p_equal_to_alpha_as_rejection() -> None:
    res = C.check_level([0.05] * 10 + [0.9] * 90, 0.05)
    assert res.observed == pytest.approx(0.10)
    assert res.ok is False  # 0.10 > 0.093589


def test_example5_check_fdr() -> None:
    res = C.check_fdr([0.0] * 49 + [0.5], 0.01)
    sd = math.sqrt((49 * 0.01**2 + 0.49**2) / 49)  # ddof = 1
    assert sd == pytest.approx(0.0707, abs=1e-4)
    assert res.observed == pytest.approx(0.01, abs=1e-12)
    assert res.bound == pytest.approx(0.01 + 2 * sd / math.sqrt(50), abs=1e-12)
    assert res.bound == pytest.approx(0.0300, abs=1e-4)
    assert res.R == 50
    assert res.ok is True


def test_example6_check_fwer() -> None:
    res = C.check_fwer([True] * 6 + [False] * 94)
    assert res.observed == pytest.approx(0.06)
    assert res.bound == pytest.approx(0.05 + 2 * math.sqrt(0.0475 / 100), abs=1e-12)
    assert res.bound == pytest.approx(0.0936, abs=1e-4)
    assert res.R == 100
    assert res.ok is True


def test_check_fwer_other_alpha_and_failure() -> None:
    res = C.check_fwer([True] * 10 + [False] * 90, alpha=0.05)
    assert res.ok is False  # 0.10 > 0.0936
    res = C.check_fwer([True] * 10 + [False] * 90, alpha=0.10)
    assert res.bound == pytest.approx(0.10 + 2 * math.sqrt(0.09 / 100), abs=1e-12)
    assert res.ok is True


def test_example7_replicates_regions() -> None:
    assert C.replicates("regions", full=False) == 100
    assert C.replicates("regions", full=True) == 1000
    assert C.replicates("gene_fdr", full=False) == 50
    assert C.replicates("gene_fdr", full=True) == 300


def test_example8_suite_seed_deterministic() -> None:
    a, b = C.suite_seed("regions", 0), C.suite_seed("regions", 0)
    assert isinstance(a, int) and a == b
    assert a != C.suite_seed("regions", 1)
    assert a != C.suite_seed("gene_fdr", 0)


# --------------------------------------------------------------------------- edge cases


def test_check_fdr_all_zero_has_bound_q() -> None:
    res = C.check_fdr([0.0] * 50, 0.05)
    assert res.observed == 0.0
    assert res.bound == 0.05
    assert res.ok is True
    assert res.R == 50


def test_check_fdr_fails_when_mean_exceeds_band() -> None:
    res = C.check_fdr([0.2] * 25 + [0.0] * 25, 0.05)  # mean 0.1, sd 0.10102, SE 0.014286
    assert res.observed == pytest.approx(0.1)
    assert res.bound == pytest.approx(
        0.05 + 2 * np.std([0.2] * 25 + [0.0] * 25, ddof=1) / math.sqrt(50)
    )
    assert res.ok is False


def test_unknown_suite_raises_keyerror() -> None:
    with pytest.raises(KeyError):
        C.replicates("no_such_suite", full=False)
    with pytest.raises(KeyError):
        C.replicates("no_such_suite", full=True)


def test_check_uniform_refuses_discrete_p() -> None:
    p = (1 + np.arange(1000) % 200) / 200.0  # permutation-like grid: 80 % ties
    with pytest.raises(ValueError, match=r"(?i)discrete|ties"):
        C.check_uniform(p)


def test_check_uniform_continuous_null_passes() -> None:
    p = np.random.default_rng(1).uniform(size=1000)
    res = C.check_uniform(p)
    assert res.observed == pytest.approx(sps.kstest(p, "uniform").pvalue, rel=1e-9)
    assert res.bound == C.load_tolerances().ks_min_p
    assert res.R == 1000
    assert res.ok is True


def test_check_uniform_skewed_p_fails() -> None:
    p = np.random.default_rng(2).beta(0.3, 1.0, size=1000)
    res = C.check_uniform(p)
    assert res.ok is False
    assert res.observed <= res.bound


def test_empty_inputs_raise() -> None:
    with pytest.raises(ValueError):
        C.check_level([], 0.05)
    with pytest.raises(ValueError):
        C.check_fdr([], 0.05)
    with pytest.raises(ValueError):
        C.check_fwer([])


# --------------------------------------------------------------------------- seeds and runner


def _spec_seed(suite: str, rep: int) -> int:
    base = int(_raw()["seeds"]["calibration_base"])  # type: ignore[index]
    ss = np.random.SeedSequence([base, zlib.crc32(suite.encode()), rep])
    return int(ss.generate_state(1, np.uint64)[0])


@pytest.mark.parametrize("suite", ["regions", "gene_fdr", "regions.perm_p"])
@pytest.mark.parametrize("rep", [0, 1, 999])
def test_suite_seed_formula(suite: str, rep: int) -> None:
    s = C.suite_seed(suite, rep)
    assert s == _spec_seed(suite, rep)
    assert 0 <= s < 2**64


def test_replicate_rng_reproduces_from_logged_seed() -> None:
    got = C.replicate_rng("regions", 7).normal(size=5)
    want = np.random.default_rng(C.suite_seed("regions", 7)).normal(size=5)
    np.testing.assert_array_equal(got, want)


def test_run_replicates_is_deterministic_and_shard_safe() -> None:
    def draw(rng: np.random.Generator) -> float:
        return float(rng.normal())

    full = C.run_replicates("regions.unit", draw, 10)
    assert full == C.run_replicates("regions.unit", draw, 10)
    left = C.run_replicates("regions.unit", draw, 4, start=0)
    right = C.run_replicates("regions.unit", draw, 6, start=4)
    assert left + right == full
    assert full == [float(C.replicate_rng("regions.unit", i).normal()) for i in range(10)]
    assert len(set(full)) == 10


def test_replicate_streams_are_thread_safe() -> None:
    def one(i: int) -> float:
        return float(C.replicate_rng("regions.threads", i).uniform())

    serial = [one(i) for i in range(32)]
    with ThreadPoolExecutor(max_workers=4) as ex:
        parallel = list(ex.map(one, reversed(range(32))))
    assert parallel[::-1] == serial


@FAST
@given(st.text(min_size=1, max_size=20), st.integers(0, 10**6))
def test_suite_seed_property(suite: str, rep: int) -> None:
    assert C.suite_seed(suite, rep) == C.suite_seed(suite, rep)
    assert C.suite_seed(suite, rep) != C.suite_seed(suite, rep + 1)


# --------------------------------------------------------------------------- band invariants


@FAST
@given(st.sampled_from([0.01, 0.05, 0.10]), st.integers(1, 10**5))
def test_level_bound_properties(alpha: float, r: int) -> None:
    b = C.level_bound(alpha, r)
    assert b > alpha
    assert C.level_bound(alpha, r + 1) < b
    k = C.load_tolerances().se_multiplier
    assert b == pytest.approx(alpha + k * math.sqrt(alpha * (1 - alpha) / r), rel=1e-12)


@FAST
@given(
    st.lists(st.floats(0.0, 1.0, allow_nan=False), min_size=2, max_size=200),
    st.sampled_from([0.01, 0.05, 0.10]),
    st.randoms(use_true_random=False),
)
def test_checks_are_order_invariant(p: list[float], a: float, rnd: object) -> None:
    q = list(p)
    rnd.shuffle(q)  # type: ignore[attr-defined]
    assert C.check_level(p, a) == C.check_level(q, a)
    r1, r2 = C.check_fdr(p, a), C.check_fdr(q, a)
    assert r1.R == r2.R
    assert r1.observed == pytest.approx(r2.observed, abs=1e-12)
    assert r1.bound == pytest.approx(r2.bound, abs=1e-12)


@FAST
@given(
    st.lists(st.floats(0.0, 1.0, allow_nan=False), min_size=2, max_size=200),
    st.sampled_from([0.01, 0.05, 0.10]),
)
def test_check_result_consistency(p: list[float], a: float) -> None:
    for res in (C.check_level(p, a), C.check_fdr(p, a), C.check_fwer([x <= a for x in p], a)):
        assert len(p) == res.R
        assert res.ok == (res.observed <= res.bound)


# --------------------------------------------------------------------------- job array


def test_array_suites_are_the_full_replicate_keys_in_file_order() -> None:
    keys = list(_raw()["replicates_full"])  # type: ignore[call-overload]
    assert list(C.array_suites()) == keys
    assert [C.task_suite(i) for i in range(len(keys))] == keys


@pytest.mark.parametrize("bad", [-1, 10**6])
def test_task_suite_out_of_range(bad: int) -> None:
    with pytest.raises(ValueError):
        C.task_suite(bad)


@FAST
@given(st.integers(0, 5000), st.integers(1, 64))
def test_shards_partition_replicates(r: int, n: int) -> None:
    blocks = [C.shard_replicates(r, n, k) for k in range(n)]
    assert all(isinstance(b, range) for b in blocks)
    flat = [i for b in blocks for i in b]
    assert flat == list(range(r))  # contiguous, ordered, exactly once
    sizes = [len(b) for b in blocks]
    assert max(sizes) - min(sizes) <= 1


def test_shard_out_of_range() -> None:
    with pytest.raises(ValueError):
        C.shard_replicates(100, 4, 4)
    with pytest.raises(ValueError):
        C.shard_replicates(100, 0, 0)


def test_parse_args_task_id_and_shard() -> None:
    suites = list(C.array_suites())
    i = suites.index("regions")
    ns = C.parse_args(["--task-id", str(i), "--shard", "1/4"])
    assert ns.suite == "regions"
    assert ns.full is True  # the job array runs the full suite by default
    assert ns.replicates == range(250, 500)


def test_parse_args_reads_slurm_array_task_id(monkeypatch: pytest.MonkeyPatch) -> None:
    suites = list(C.array_suites())
    monkeypatch.setenv("SLURM_ARRAY_TASK_ID", str(suites.index("gene_fdr")))
    ns = C.parse_args([])
    assert ns.suite == "gene_fdr"
    assert ns.replicates == range(0, 300)


def test_parse_args_fast_and_suite_name(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
    ns = C.parse_args(["--suite", "regions", "--fast"])
    assert ns.suite == "regions" and ns.full is False
    assert ns.replicates == range(0, 100)


def test_parse_args_requires_a_task(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
    with pytest.raises(SystemExit):
        C.parse_args([])
    with pytest.raises(SystemExit):
        C.parse_args(["--suite", "no_such_suite"])


def test_main_dry_run_prints_pytest_command_without_submitting(
    capsys: pytest.CaptureFixture[str], monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.delenv("SLURM_ARRAY_TASK_ID", raising=False)
    i = list(C.array_suites()).index("regions")
    assert C.main(["--task-id", str(i), "--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "regions" in out and "calibration_full" in out and "pytest" in out


def test_pytest_args_select_one_full_suite() -> None:
    args = C.pytest_args("regions", full=True)
    assert "tests/calibration" in " ".join(args)
    assert "-m" in args and args[args.index("-m") + 1] == "calibration_full"
    assert "regions" in " ".join(args)
    fast = C.pytest_args("regions", full=False)
    assert fast[fast.index("-m") + 1] == "calibration and not calibration_full"


# --------------------------------------------------------------------------- records and report


def test_log_result_appends_jsonl(tmp_path: Path) -> None:
    res = C.check_level([0.001] * 2 + [0.5] * 98, 0.01)
    p1 = C.log_result("regions", "level", 0.01, res, seed=123, duration=1.5, log_dir=tmp_path)
    p2 = C.log_result(
        "regions", "level", 0.05, C.check_level([0.5] * 100, 0.05), 123, 2.0, tmp_path
    )
    assert p1 == p2 == tmp_path / "regions.jsonl"
    lines = p1.read_text().splitlines()
    assert len(lines) == 2
    rec = json.loads(lines[0])
    assert rec == {
        "suite": "regions",
        "criterion": "level",
        "level": 0.01,
        "observed": pytest.approx(0.02),
        "bound": pytest.approx(C.level_bound(0.01, 100)),
        "R": 100,
        "seed": 123,
        "duration": 1.5,
        "ok": True,
    }
    assert C.read_results(p1)[1]["level"] == 0.05


def test_log_dir_env_override(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setenv("GIAP_CALIBRATION_LOG_DIR", str(tmp_path / "logs"))
    path = C.log_result("gene_fdr", "fdr", 0.05, C.check_fdr([0.0] * 50, 0.05), 1, 0.1)
    assert path == tmp_path / "logs" / "gene_fdr.jsonl"
    assert path.exists()


def test_default_log_dir_is_local_logs_calibration(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GIAP_CALIBRATION_LOG_DIR", raising=False)
    d = C.log_dir()
    assert d.parts[-3:] == ("local", "logs", "calibration")


def test_summary_markdown(tmp_path: Path) -> None:
    C.log_result(
        "regions", "level", 0.01, C.check_level([0.001] * 3 + [0.5] * 97, 0.01), 1, 0.2, tmp_path
    )
    C.log_result("regions", "level", 0.05, C.check_level([0.5] * 100, 0.05), 1, 0.3, tmp_path)
    md = C.summary_markdown(C.read_results(tmp_path / "regions.jsonl"))
    lines = [ln for ln in md.strip().splitlines() if ln.startswith("|")]
    header = lines[0].lower()
    for col in ("suite", "criterion", "level", "observed", "bound", "r", "ok"):
        assert col in header
    body = lines[2:]
    assert len(body) == 2
    assert "FAIL" in body[0] and "pass" in body[1]


# --------------------------------------------------------------------------- protection


def test_harness_never_writes_tolerances(tmp_path: Path) -> None:
    before = (hashlib.sha256(TOL_PATH.read_bytes()).hexdigest(), TOL_PATH.stat().st_mtime_ns)
    C.load_tolerances()
    C.replicates("regions", full=True)
    C.suite_seed("regions", 0)
    C.check_level([0.5] * 10, 0.05)
    C.check_fdr([0.0] * 10, 0.05)
    C.check_fwer([False] * 10)
    C.check_uniform(np.random.default_rng(0).uniform(size=100))
    C.log_result("regions", "level", 0.05, C.check_level([0.5] * 10, 0.05), 0, 0.0, tmp_path)
    after = (hashlib.sha256(TOL_PATH.read_bytes()).hexdigest(), TOL_PATH.stat().st_mtime_ns)
    assert before == after


def test_no_tolerance_constants_duplicated_in_harness() -> None:
    src = (HERE / "_criteria.py").read_text()
    assert str(_raw()["seeds"]["calibration_base"]) not in src  # type: ignore[index]
    assert "tolerances.yaml" in src  # constants come from the protected file


# --------------------------------------------------------------------------- `full` fixture


def test_full_fixture_false_when_unmarked(full: bool) -> None:
    assert full is False


@pytest.mark.calibration
def test_full_fixture_false_under_calibration(full: bool) -> None:
    assert full is False


@pytest.mark.calibration_full
def test_full_fixture_true_under_calibration_full(full: bool) -> None:
    assert full is True


# ------------------------------------------------------------------ operating characteristics
# Statistical claims of the methods note (docs/methods/calibration.md): on an exact null the
# level check fails with probability P(Bin(R, alpha) / R > level_bound) (same float comparison as
# the check), and the uniformity check fails U(0, 1) samples at rate <= ks_min_p. Over M
# meta-replicates the failure count is compared with that probability by an exact binomial test
# passing at p > ks_min_p (one threshold for every level, no new constants).

M_META = 2000


@pytest.mark.calibration
def test_level_check_false_alarm_rate_matches_binomial_tail() -> None:
    tol = C.load_tolerances()
    r = C.replicates("regions", full=False)
    for a in tol.levels:
        fails = C.run_replicates(
            f"harness.level_oc.{a}",
            lambda rng, a=a: not C.check_level(rng.uniform(size=r), a).ok,
            M_META,
        )
        bound = C.level_bound(a, r)
        counts = np.arange(r + 1)
        exact = float(sps.binom.pmf(counts, r, a)[counts / r > bound].sum())
        n_fail = int(np.sum(fails))
        test = sps.binomtest(n_fail, M_META, exact)
        assert test.pvalue > tol.ks_min_p, (a, n_fail / M_META, exact, test.pvalue)


@pytest.mark.calibration
def test_uniform_check_false_alarm_rate_is_ks_min_p() -> None:
    tol = C.load_tolerances()
    fails = C.run_replicates(
        "harness.uniform_oc", lambda rng: not C.check_uniform(rng.uniform(size=1000)).ok, M_META
    )
    rate = float(np.mean(fails))
    a = tol.ks_min_p
    assert rate <= a + tol.se_multiplier * math.sqrt(a * (1 - a) / M_META) + 1.0 / M_META, rate


@pytest.mark.calibration
def test_fdr_check_detects_inflated_fdp() -> None:
    """Power side: FDP with true mean 4q fails in >= 95 % of meta-replicates (gene_fdr fast)."""
    r = C.replicates("gene_fdr", full=False)
    for q in C.load_tolerances().levels:
        flagged = C.run_replicates(
            f"harness.fdr_power.{q}",
            lambda rng, q=q: not C.check_fdr(rng.beta(4 * q * 20, 20 - 4 * q * 20, size=r), q).ok,
            200,
        )
        assert float(np.mean(flagged)) >= 0.95, q
