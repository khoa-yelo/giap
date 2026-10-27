"""Attack nulls for the M7-10 class false-ALT rate (GIPPA_SPEC.md §7.5.4, §8.1 row "ALT call is
real", §8.5; docs/specs/M7-10-class-error-rates.md rules 1-5; statistician review of PR #239).

The M7-10 calibration gate (tests/calibration/test_variant_class_error.py) first ran only on a
scenario without planted variants, where truth and estimate both sit at 0.97-1.0, so a constant
``rate = 1`` estimator passed it. The planted arm now gates too (lead decision 9). These attacks
work on the pure functions with synthetic paired tables whose per-class error rate is planted,
so a broken estimator fails:

1. ``known_rate_mix``: planted rates 0, 0.05, 0.3 and 0.8 with a *variant-aware* reference (the
   reference search can return the ALT peptide). The estimate must equal the brute-force paired
   oracle exactly and be unbiased for the planted rate
   (|bias| <= ``denovo_precision_max_abs_bias``), and the Wilson 95% interval from (rate, n) must
   cover the planted rate in at least 0.95 - 3 SE(R) of the replicates.
   Traps: correct reads that differ from the reference only by I/L; spectrum IDs reused across runs
   with different peptides; q = 0.02 and confident decoy rows that disagree. A constant estimator,
   a missing I/L fold, a (spectrum_id)-only key, or a missing q / decoy filter fails here.
2. ``reference_db_degenerate`` (characterisation, rule 5): when the reference search uses a
   database without the variant, a correct ALT read can never agree with its spectrum's reference
   identification, so the estimate is exactly 1 whatever the true error rate. The "upper bound" of
   rule 5 is then vacuous; the methods note must say so.
3. ``quality_selection`` (characterisation): when de novo errors and confident reference IDs both
   depend on spectrum quality, the paired rate is biased *low* for the all-reads error rate. The
   rate is therefore not an upper bound in general; it estimates the error rate on spectra a
   database search identifies confidently.
4. ``duplicate_charge_rows``: preflight keeps one spectrum's de novo read at two charges
   (``validate_peptides`` dedups on (run, spectrum_id, peptide, charge)) and C17 drops ``charge``,
   so one spectrum yields two identical C17 rows. Rule 2 counts one observation per spectrum;
   counting both is pseudo-replication (n inflated, interval too narrow).
5. ``reference_false_matches`` (characterisation): a confident reference ID that is wrong turns a
   correct de novo read into an "error", so the estimate is e + (1 - e) * phi for false-match
   share phi; this explains the small upward bias of the planted calibration arm.

R = ``replicates_fast.variants`` (``replicates_full.variants`` under ``calibration_full``).
Seeds: ``replicate_rng("variants.class_error.attack.<name>", i)``.
"""

from __future__ import annotations

import string
from collections.abc import Iterator

import _criteria as C
import numpy as np
import pandas as pd
import pytest
import yaml

from gippa.variants.error_rates import class_rates, paired_observations

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("variants")]

SUITE = "variants"
BASE = "variants.class_error.attack"
Q_MAX = 0.01
# I and L excluded so a planted error never vanishes under folding.
_AA = np.array(list(set(string.ascii_uppercase) - set("BJOUXZIL")))
PLANTED = {
    ("substitution", False): 0.05,
    ("substitution", True): 0.30,
    ("block", True): 0.80,
    ("swap", False): 0.0,
}


def _max_bias() -> float:
    raw = yaml.safe_load(C.TOLERANCES_PATH.read_text())
    return float(raw["denovo_precision_max_abs_bias"])


def _seeds(name: str, full: bool) -> Iterator[np.random.Generator]:
    R = C.replicates(SUITE, full)
    for i in C.replicate_range(R):
        yield C.replicate_rng(f"{BASE}.{name}", i)


def _peptide(rng: np.random.Generator, length: int = 12) -> str:
    return "".join(rng.choice(_AA, size=length))


def _mutate(rng: np.random.Generator, pep: str) -> str:
    """One residue changed to a different non-I/L residue (differs after folding)."""
    i = int(rng.integers(0, len(pep)))
    alt = rng.choice([a for a in _AA if a != pep[i]])
    return pep[:i] + str(alt) + pep[i + 1 :]


def _with_il(rng: np.random.Generator, pep: str) -> tuple[str, str]:
    """``pep`` with an L planted, and the same with that position as I (equal after folding)."""
    i = int(rng.integers(0, len(pep)))
    return pep[:i] + "L" + pep[i + 1 :], pep[:i] + "I" + pep[i + 1 :]


def _wilson(e: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    p = e / n
    den = 1 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return mid - half, mid + half


def _support(rows: list[tuple[str, str, bool]]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=["peptide", "alt_type", "mass_equivalent"])
    df.insert(0, "variant_id", [f"v{i}" for i in range(len(df))])
    return df.astype({"mass_equivalent": bool})


def _ref(rows: list[tuple[str, str, str, float, bool]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["run", "spectrum_id", "peptide", "q_value", "is_decoy"])


def _known_rate_replicate(
    rng: np.random.Generator, n_spec: int = 300, p_pair: float = 0.6
) -> tuple[dict, dict, dict]:
    """Returns ``(estimate {class: (rate, n, errors)}, paired oracle {class: (n, e)},
    all-reads truth {class: (n, e)})``."""
    dn, ref, sup = [], [], []
    truth: dict = {}
    k = 0
    for cls, p in PLANTED.items():
        for _ in range(n_spec):
            run = str(rng.choice(["a", "b"]))
            sid = f"s{k % (n_spec * 3)}"  # spectrum IDs reused across classes and runs
            k += 1
            true_l, true_i = _with_il(rng, _peptide(rng))
            err = bool(rng.random() < p)
            read = _mutate(rng, true_l) if err else true_i  # correct reads differ only by I/L
            dn.append((run, sid, read))
            sup.append((read, *cls))
            t = truth.setdefault(cls, [0, 0])
            t[0] += 1
            t[1] += int(err)
            u = rng.random()
            if u < p_pair:  # confident target: the variant-aware search finds the truth
                ref.append((run, sid, true_l, float(rng.uniform(0, Q_MAX)), False))
            elif u < p_pair + 0.2:  # not confident: must be ignored
                ref.append((run, sid, _peptide(rng), 0.02, False))
            elif u < p_pair + 0.3:  # confident decoy: must be ignored
                ref.append((run, sid, _peptide(rng), 0.001, True))
    # Run/spectrum-ID collisions between classes are possible; keep one reference per spectrum
    # (preflight rejects conflicts, ADR-0005) and recompute the oracle from what remains.
    dnf = pd.DataFrame(dn, columns=["run", "spectrum_id", "peptide"])
    reff = _ref(ref).drop_duplicates(["run", "spectrum_id"], keep=False)
    dnf = dnf.drop_duplicates(["run", "spectrum_id"], keep=False)
    supf = _support(sup)
    obs = paired_observations(dnf, reff, supf, q_max=Q_MAX)
    est = {(r.alt_type, r.mass_equivalent): (r.rate, r.n, r.errors) for r in class_rates(obs)}
    # brute-force paired oracle on the de-collided tables
    conf = reff[(reff["q_value"] <= Q_MAX) & ~reff["is_decoy"]]
    cols = zip(conf["run"], conf["spectrum_id"], conf["peptide"], strict=True)
    best = {(a, b): c for a, b, c in cols}
    cls_of: dict[str, set] = {}
    for pep, t, me in sup:
        cls_of.setdefault(pep, set()).add((t, me))
    oracle: dict = {}
    for run, sid, read in dnf.itertuples(index=False):
        hit = best.get((run, sid))
        if hit is None:
            continue
        e = read.replace("I", "L") != hit.replace("I", "L")
        for c in cls_of[read]:
            o = oracle.setdefault(c, [0, 0])
            o[0] += 1
            o[1] += int(e)
    return est, {c: tuple(v) for c, v in oracle.items()}, {c: tuple(v) for c, v in truth.items()}


def _known_rate_mix(full: bool) -> None:
    max_bias = _max_bias()
    diffs: dict = {c: [] for c in PLANTED}
    cover: dict = {c: [] for c in PLANTED}
    for rng in _seeds("known_rate_mix", full):
        est, oracle, _truth = _known_rate_replicate(rng)
        assert set(est) == set(oracle), (est, oracle)
        for c, (n, e) in oracle.items():
            rate, n_hat, e_hat = est[c]
            assert (n_hat, e_hat) == (n, e), f"class {c}: estimate {(n_hat, e_hat)} != oracle"
            assert rate == pytest.approx(e / n)
            diffs[c].append(rate - PLANTED[c])
            lo, hi = _wilson(e, n)
            cover[c].append(lo <= PLANTED[c] <= hi)
    for c, d in diffs.items():
        bias = float(np.mean(d))
        assert abs(bias) <= max_bias, f"class {c}: bias {bias:.4f} > {max_bias}"
        if PLANTED[c] > 0:  # Wilson at p = 0 covers trivially
            share, R = float(np.mean(cover[c])), len(cover[c])
            floor = 0.95 - 3 * np.sqrt(0.95 * 0.05 / R)  # 3 binomial SEs below nominal
            assert share >= floor, f"class {c}: Wilson coverage {share:.2f} < {floor:.2f}"


def test_known_rate_mix() -> None:
    _known_rate_mix(False)


@pytest.mark.calibration_full
def test_known_rate_mix_full() -> None:
    _known_rate_mix(True)


def test_reference_db_degenerate() -> None:
    """Reference search against a database without the variant: every paired ALT-supporting read
    disagrees, so the estimate is 1 while the true error rate is 0.2 (rule 5 is vacuous here)."""
    for rng in _seeds("reference_db_degenerate", False):
        dn, ref, sup = [], [], []
        n_err = n_all = 0
        for k in range(400):
            sid = f"s{k}"
            ref_pep = _peptide(rng)  # the reference-proteome peptide at the locus
            err = bool(rng.random() < 0.2)
            if err:  # a reference-peptide spectrum misread as an ALT; the search IDs it right
                read = _mutate(rng, ref_pep)
                if rng.random() < 0.6:
                    ref.append(("r", sid, ref_pep, 0.001, False))
            else:  # a true variant spectrum read correctly; the ALT peptide is not searchable
                read = _mutate(rng, ref_pep)
                if rng.random() < 0.3:  # sometimes the search still IDs the REF peptide
                    ref.append(("r", sid, ref_pep, 0.005, False))
            n_all += 1
            n_err += int(err)
            dn.append(("r", sid, read))
            sup.append((read, "substitution", False))
        obs = paired_observations(
            pd.DataFrame(dn, columns=["run", "spectrum_id", "peptide"]), _ref(ref), _support(sup)
        )
        (r,) = class_rates(obs)
        assert r.n > 0
        assert r.rate == 1.0, "the estimator is not degenerate here: revisit the methods note"
        assert n_err / n_all < 0.3


def test_quality_selection() -> None:
    """De novo error prob 0.6 (1 - s) and pairing prob s, s ~ U(0, 1): all-reads truth 0.30,
    paired truth 0.20. The estimate tracks the paired truth and is biased low for all reads."""
    max_bias = _max_bias()
    d_all, d_paired = [], []
    for rng in _seeds("quality_selection", False):
        dn, ref, sup = [], [], []
        n_err = 0
        n = 2000
        for k in range(n):
            s = float(rng.random())
            true = _peptide(rng)
            err = bool(rng.random() < 0.6 * (1 - s))
            read = _mutate(rng, true) if err else true
            n_err += int(err)
            dn.append(("r", f"s{k}", read))
            sup.append((read, "substitution", False))
            if rng.random() < s:
                ref.append(("r", f"s{k}", true, 0.001, False))
        obs = paired_observations(
            pd.DataFrame(dn, columns=["run", "spectrum_id", "peptide"]), _ref(ref), _support(sup)
        )
        (r,) = class_rates(obs)
        d_all.append(r.rate - n_err / n)
        d_paired.append(r.rate - 0.2)
    assert float(np.mean(d_all)) < -0.05, "expected a downward selection bias vs all reads"
    assert abs(float(np.mean(d_paired))) <= max_bias


def test_duplicate_charge_rows_count_once() -> None:
    """One spectrum read at charges 2 and 3 becomes two identical C17 rows; rule 2 counts one
    observation per de novo spectrum, so n must be 1, not 2."""
    dn = pd.DataFrame(
        [("r", "s1", "PEPTMDEKAA"), ("r", "s1", "PEPTMDEKAA")],
        columns=["run", "spectrum_id", "peptide"],
    )
    ref = _ref([("r", "s1", "PEPTWDEKAA", 0.001, False)])
    obs = paired_observations(dn, ref, _support([("PEPTMDEKAA", "substitution", False)]))
    (r,) = class_rates(obs)
    assert (r.n, r.errors) == (1, 1), f"one spectrum counted {r.n} times"


def test_reference_false_matches() -> None:
    """Confident reference IDs are wrong at share phi (q <= 0.01 bounds the FDP only on average):
    a correct de novo read on such a spectrum counts as an error. The estimate is then
    e + (1 - e) * phi for de novo error rate e on paired spectra, biased up by (1 - e) * phi, and
    estimate - paired truth equals (spurious - missed) / n exactly. This is the +0.004 to +0.009
    bias of the planted calibration arm (methods note, "Reference false matches")."""
    e_true, phi, n = 0.2, 0.05, 2000
    d_formula = []
    for rng in _seeds("reference_false_matches", False):
        dn, ref, sup = [], [], []
        n_err = spurious = missed = 0
        for k in range(n):
            true = _peptide(rng)
            err = bool(rng.random() < e_true)
            read = _mutate(rng, true) if err else true
            ref_wrong = bool(rng.random() < phi)
            hit = (read if rng.random() < 0.1 else _mutate(rng, true)) if ref_wrong else true
            n_err += int(err)
            spurious += int(not err and hit != true)
            missed += int(err and hit == read)
            dn.append(("r", f"s{k}", read))
            sup.append((read, "substitution", False))
            ref.append(("r", f"s{k}", hit, 0.001, False))
        obs = paired_observations(
            pd.DataFrame(dn, columns=["run", "spectrum_id", "peptide"]), _ref(ref), _support(sup)
        )
        (r,) = class_rates(obs)
        assert r.n == n
        assert r.rate == pytest.approx(n_err / n + (spurious - missed) / n)
        d_formula.append(r.rate - (e_true + (1 - e_true) * phi))
    assert abs(float(np.mean(d_formula))) <= 0.01, "estimate != e + (1 - e) * phi on average"
