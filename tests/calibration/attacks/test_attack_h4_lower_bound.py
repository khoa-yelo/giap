"""Attack (GIPPA_SPEC.md §8.4, §8.5, statistician): hypothesis 4's verdict rule (M10-10b;
docs/specs/M10-10-a3-run.md, "Statistician review (2026-10-24)" > "Attack plan", S2 and S3).

H4 compares GIPPA's genome-anchored de novo precision per score bin with the A2-confirmed share
``L_b``, a lower bound on the bin's realised precision. The verdict is ``not_supported`` when any
reportable bin has ``relation_adj = below_lower_bound`` (U_b, the genome CI's upper half-width
rescaled to 1 - 0.025/B, below the Wilson lower limit at 1 - 0.05/B). Under H4 (the genome
estimator is unbiased) a ``not_supported`` is a false refutation; its family-wise rate over the B
bins must stay at or below 0.05 + k sqrt(0.05 · 0.95 / R) (k = ``se_multiplier`` = 2).

Generator (one replicate): 10 latent score groups with true per-peptide correctness
``P_BIN`` = 0.30 … 0.99, each group's scores inside its own slice of Cascadia's (0.8, 1] (so the
step's decile bins recover the groups unless scores tie). Peptides are random strings without I
(9 or 12 residues: the length classes "9" and "12-14"). A correct peptide matches the target
exactly with probability c_true, a wrong one with r; every peptide matches the shuffled decoy with
r. The reference c_l comes from ``c_by_class`` on n_ref = 2,000 reference peptides (half of each
length) matching at c_ref. A2 (D2) confirms each correct peptide with probability φ and each wrong
one with probability ε; half of D2's strings are written with I for L (DIA-NN), and D2 also holds
1,000 unrelated peptides. GIPPA's side is the real ``denovo_precision`` (1,000 bootstrap
resamples, values rounded to ``%.6g`` as the TSV writer does) on the aggregated (score = max)
queries; H4's side is the real ``hypothesis4`` (and, in ``ties_and_reps``, the real
``denovo_queries`` on a written Cascadia-style table).

Scenarios (plan table): gated at FWER — ``worst_case``, ``small_bins``, ``low_denominator``,
``partial_confirm``, ``a2_false_confirm``, ``ties_and_reps`` (plus: every peptide's bin equals
GIPPA's, S1), ``contaminants``; reported, not gated — ``clustered`` (the documented CI limit, S3:
runs of 5 dependent, nested peptides); power reported — ``power_easy_reference``,
``power_one_bin``. Every scenario also asserts that the verdict is never ``supported`` and that
``relation_95 ⊇ relation_adj``. Reported values go to ``compare_attack_h4_recorded.jsonl`` and the
gates to ``compare.jsonl`` (``_criteria.log_result``).

R = 200 (marker ``calibration``) or 1,000 (``calibration_full``, Slurm). Seeds:
``gippa.rng.for_step("attack_h4", seed)`` with seed = ``seeds.calibration_base`` + 1010 + 10,000 ·
(scenario index) + replicate, one generator per replicate for the simulation and the bootstrap,
so the fast replicates are the first 200 of the full run. If a gated scenario fails, the merge
blocks (§8.5); the fallback, fixed in the spec, is ``relation_adj`` in two adjacent bins.

Measured at 1e06c770 (statistician review of PR #366): every gated scenario passes except
``a2_false_confirm`` (0.085 at R = 200, bound 0.081; 0.095 at R = 1,000, bound 0.064). Its
refutations sit in the top bins, where 1 % false confirmations of k_b means ε ≈ 1 (p = 0.99) and
``L_b`` exceeds the realised precision: the lower bound itself fails by A2's FDR, not the U_b
rescaling. The test is left failing (R3: no bound is widened); the lead decides the fix.
"""

from __future__ import annotations

import functools
import gzip
import json
import math
import sys
import tempfile
import time
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest

from gippa.config import GippaConfig
from gippa.mapping.denovo_precision import c_by_class, denovo_precision, score_bins
from gippa.rng import for_step
from gippa.steps.mapping import length_label

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from helpers.pilot_metrics_m10 import load_metrics

SUITE = "compare"
pytestmark = [pytest.mark.calsuite(SUITE)]
MODES = [pytest.param(False, marks=pytest.mark.calibration, id="fast"),
         pytest.param(True, marks=pytest.mark.calibration_full, id="full")]  # fmt: skip

R_FAST = 200
R_FULL = 1_000
SEED_OFFSET = 1010
ALPHA = 0.05  # the FWER level of the verdict (S2)
P_BIN = np.array([0.30, 0.45, 0.60, 0.70, 0.80, 0.85, 0.90, 0.95, 0.98, 0.99])
AA = np.array(list("ACDEFGHKLMNPQRSTVWY"))  # 19 residues, no I (folding is the identity)
LENGTHS = (9, 12)
N_REF = 2_000
N_D2_OTHER = 1_000
N_BOOT = 1_000
CFG = GippaConfig(seed=7)
CLASSES = CFG.mapping.length_classes


@dataclass(frozen=True)
class Scenario:
    n: int = 2_000  # peptides per latent group
    c_true: float = 0.95
    c_ref: float = 0.95
    r: float = 0.08
    phi: float = 1.0  # A2 confirms a correct peptide
    false_confirm: float = 0.0  # share of k_b that are false confirmations (sets ε per group)
    reps: bool = False  # 1-3 rows per peptide, first row lowest, top groups tied at 1.0
    contam: float = 0.0  # share of de novo peptides that are correct contaminants
    ref_contam: float = 0.0  # share of reference rows that are contaminants (never match)
    cluster: int = 1  # runs of peptides sharing correctness and match status
    bias_group: int | None = None  # group whose genome precision is biased down by 0.2
    role: str = "gated"  # gated | reported | power
    note: str = field(default="", compare=False)


SCENARIOS: dict[str, Scenario] = {
    "worst_case": Scenario(note="phi 1, eps 0, n 2000/bin, c 0.95, r 0.08"),
    "small_bins": Scenario(n=25, note="n 25/bin (just over n >= 20)"),
    "low_denominator": Scenario(c_true=0.6, c_ref=0.6, r=0.4, note="c - r = 0.2"),
    "partial_confirm": Scenario(phi=0.6, n=500, note="phi 0.6, n 500"),
    "a2_false_confirm": Scenario(false_confirm=0.01, note="1 % of k_b false confirmations"),
    "ties_and_reps": Scenario(n=500, reps=True, note="1-3 rows, first lowest; top tied at 1.0"),
    "contaminants": Scenario(contam=0.05, ref_contam=0.03, note="5 % de novo, 3 % reference"),
    "clustered": Scenario(cluster=5, role="reported", note="runs of 5 nested peptides"),
    "power_easy_reference": Scenario(n=500, c_true=0.85, c_ref=0.99, role="power",
                                     note="c_ref 0.99 > c_true 0.85"),
    "power_one_bin": Scenario(n=500, bias_group=6, role="power",
                              note="group 6 (p 0.90) biased down by 0.2"),
}  # fmt: skip
ORDER = list(SCENARIOS)


def _base_seed() -> int:
    return C.load_tolerances().calibration_base + SEED_OFFSET


def _k() -> float:
    return C.load_tolerances().se_multiplier


def _bound(R: int) -> float:
    return ALPHA + _k() * math.sqrt(ALPHA * (1.0 - ALPHA) / R)


# ------------------------------------------------------------------ one replicate


def _strings(rng: np.random.Generator, n: int, length: int) -> np.ndarray:
    codes = rng.integers(0, AA.size, size=(n, length))
    return np.array(["".join(row) for row in AA[codes]], dtype=object)


def _group(rng: np.random.Generator, s: Scenario, b: int) -> pd.DataFrame:
    """Peptides of latent group ``b``: strings, max score, correctness, matches, D2 status."""
    p = float(P_BIN[b])
    n, k = s.n, s.cluster
    if k == 1:
        lengths = rng.choice(LENGTHS, size=n)
        peps = np.empty(n, dtype=object)
        for L in LENGTHS:
            sel = lengths == L
            peps[sel] = _strings(rng, int(sel.sum()), L)
        lead = np.arange(n)
    else:  # one 12-mer and its 9-mers at offsets 0..k-2: nested peptides of one locus
        n_cl = -(-n // k)
        parents = _strings(rng, n_cl, 12)
        peps = np.array([q if j == 0 else q[j - 1 : j + 8] for q in parents for j in range(k)],
                        dtype=object)[:n]  # fmt: skip
        lead = np.repeat(np.arange(n_cl), k)[:n]
    m = int(lead.max()) + 1
    correct = (rng.random(m) < p)[lead]
    u_t, u_d = rng.random(m)[lead], rng.random(m)[lead]
    c_b = s.c_true
    if s.bias_group == b:  # p (c_b - r) / (c_ref - r) = p - 0.2
        c_b = s.r + (s.c_ref - s.r) * (1.0 - 0.2 / p)
    contam = rng.random(n) < s.contam
    correct = correct | contam
    target = np.where(correct, u_t < c_b, u_t < s.r) & ~contam
    decoy = u_d < s.r
    eps = 0.0
    if s.false_confirm and p < 1.0:  # n(1-p) eps = f (n p phi + n(1-p) eps)
        f = s.false_confirm
        eps = min(1.0, f * p * s.phi / ((1.0 - f) * (1.0 - p)))
    u_c = rng.random(n)
    in_d2 = np.where(correct & ~contam, u_c < s.phi, np.where(contam, False, u_c < eps))
    lo, hi = 0.8 + 0.02 * b, 0.8 + 0.02 * (b + 1)
    score = rng.uniform(lo, hi, n)
    score[score <= lo] = hi  # (lo, hi]
    if s.reps and b >= 7:
        score[:] = 1.0  # Cascadia's top ties: the step's edges coincide and bins merge
    score = np.array([float(f"{x:.6g}") for x in score])  # as the converter writes them
    return pd.DataFrame({"peptide": peps, "score": score, "correct": correct, "target": target,
                         "decoy": decoy, "in_d2": in_d2})  # fmt: skip


def _simulate(rng: np.random.Generator, s: Scenario) -> pd.DataFrame:
    df = pd.concat([_group(rng, s, b) for b in range(P_BIN.size)], ignore_index=True)
    return df.drop_duplicates("peptide").reset_index(drop=True)  # chance collisions: one query


def _reference(rng: np.random.Generator, s: Scenario) -> dict[str, float]:
    lengths = np.repeat(LENGTHS, N_REF // len(LENGTHS))
    cont = rng.random(lengths.size) < s.ref_contam
    ref = pd.DataFrame({"peptide": ["A" * int(L) for L in lengths], "q_value": 0.001,
                        "is_decoy": False,
                        "target_exact": (rng.random(lengths.size) < s.c_ref) & ~cont})  # fmt: skip
    c, assumed = c_by_class(ref, CLASSES)
    assert not assumed
    return c


def _g6(v: float | None) -> float | None:
    return None if v is None else float(f"{v:.6g}")


def _mapping_fdr(
    df: pd.DataFrame, c: Mapping[str, float], rng: np.random.Generator
) -> pd.DataFrame:
    """GIPPA's ``denovo_precision`` rows as ``mapping_fdr.tsv`` carries them."""
    frame = pd.DataFrame({"peptide": df["peptide"], "score": df["score"],
                          "length_class": [length_label(len(p), CLASSES) for p in df["peptide"]],
                          "target_exact": df["target"], "decoy_exact": df["decoy"]})  # fmt: skip
    rows = denovo_precision(frame, c, False, rng, n_boot=N_BOOT)
    return pd.DataFrame(
        [{"section": "denovo_precision", "decoy_type": "shuffle", "score_bin": b.score_bin,
          "precision": _g6(b.precision), "ci_low": _g6(b.ci_low), "ci_high": _g6(b.ci_high),
          "c_value": _g6(b.c_value), "flags": ";".join(b.flags)} for b in rows]
    )  # fmt: skip


def _d2(rng: np.random.Generator, df: pd.DataFrame) -> set[str]:
    conf = df.loc[df["in_d2"].to_numpy(dtype=bool), "peptide"].astype(str).tolist()
    flip = rng.random(len(conf)) < 0.5  # DIA-NN may write I where Cascadia writes L
    out = {q.replace("L", "I") if f else q for q, f in zip(conf, flip, strict=True)}
    return out | set(_strings(rng, N_D2_OTHER, 15).tolist())


def _cascadia_table(rng: np.random.Generator, df: pd.DataFrame, path: Path) -> Path:
    """1-3 rows per peptide; the planted score is the maximum and the first row the lowest."""
    lines = ["peptide\tmodified_peptide\tscore\trun\tcharge\tspectrum_id\tsource"]
    n_rows = rng.integers(1, 4, size=len(df))
    scan = 0
    for pep, sc, nr in zip(df["peptide"], df["score"], n_rows, strict=True):
        d = np.sort(rng.uniform(0.0, sc - 0.8, size=int(nr) - 1))[::-1]  # first row lowest
        scores = [float(sc - d[0])] if nr > 1 else []
        scores += [float(sc)] + [float(sc - x) for x in d[1:]]
        for j, x in enumerate(scores):
            scan += 1
            run = ("WT_1", "WT_2", "WT_3")[j]
            lines.append(f"{pep}\t{pep}\t{x:.6g}\t{run}\t2\t{run}:scan={scan}\tdenovo")
    with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fh:
        fh.write(("\n".join(lines) + "\n").encode("utf-8"))
    return path


def _queries(rng: np.random.Generator, s: Scenario, df: pd.DataFrame) -> pd.DataFrame:
    M = load_metrics()
    if not s.reps:
        return pd.DataFrame({"peptide": df["peptide"], "score": df["score"], "source": "denovo"})
    from gippa.io.peptides import read_peptide_table

    with tempfile.TemporaryDirectory() as tmp:
        P, _ = read_peptide_table(_cascadia_table(rng, df, Path(tmp) / "a3.tsv.gz"))
    q = M.denovo_queries(P, CFG)
    # S1: every query in the bin GIPPA's step gives it (score = max over the peptide's rows)
    gippa_bin = dict(zip(df["peptide"], score_bins(df["score"].to_numpy(dtype=float)), strict=True))
    ours = score_bins(pd.to_numeric(q["score"]).to_numpy(dtype=float))
    assert len(q) == len(df), (len(q), len(df))
    bad = sum(gippa_bin[p] != b for p, b in zip(q["peptide"].astype(str), ours, strict=True))
    assert bad == 0, f"{bad} queries in another bin than GIPPA's (S1)"
    return q


@dataclass(frozen=True)
class Rep:
    verdict: str
    any95: bool  # some reportable bin below at 95 % (the unadjusted rule)
    B: int
    n_bins: int
    nested: float  # mean per-bin nested_share


def _replicate(s: Scenario, seed: int) -> Rep:
    M = load_metrics()
    rng = for_step("attack_h4", seed)
    df = _simulate(rng, s)
    c = _reference(rng, s)
    fdr = _mapping_fdr(df, c, rng)
    d2 = _d2(rng, df)
    q = _queries(rng, s, df)
    out = M.hypothesis4(q, d2, fdr, n_denovo_logged=len(q))
    assert out["verdict"] in M.H4_VERDICTS and out["verdict"] != "supported", out["verdict"]
    rep = [b for b in out["bins"] if b["reportable"]]
    for b in rep:  # relation_95 ⊇ relation_adj
        if b["relation_adj"] == "below_lower_bound":
            assert b["relation_95"] == "below_lower_bound", b
    nested = [b["nested_share"]["value"] for b in out["bins"] if b["n_b"]]
    return Rep(out["verdict"], any(b["relation_95"] == "below_lower_bound" for b in rep),
               out["B"], len(out["bins"]), float(np.mean(nested)))  # fmt: skip


@dataclass(frozen=True)
class Run:
    name: str
    R: int
    not_supported: int
    not_evaluable: int
    any95: int
    mean_B: float
    mean_bins: float
    mean_nested: float
    seconds: float

    @property
    def rate(self) -> float:
        return self.not_supported / self.R

    @property
    def rate95(self) -> float:
        return self.any95 / self.R


@functools.cache
def run(name: str, R: int) -> Run:
    s = SCENARIOS[name]
    seed0 = _base_seed() + 10_000 * ORDER.index(name)
    t0 = time.perf_counter()
    reps = [_replicate(s, seed0 + i) for i in range(R)]
    return Run(name, R, sum(x.verdict == "not_supported" for x in reps),
               sum(x.verdict == "not_evaluable" for x in reps), sum(x.any95 for x in reps),
               float(np.mean([x.B for x in reps])), float(np.mean([x.n_bins for x in reps])),
               float(np.mean([x.nested for x in reps])), time.perf_counter() - t0)  # fmt: skip


def _log_recorded(record: dict[str, Any]) -> None:
    path = C.log_dir() / f"{SUITE}_attack_h4_recorded.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as fh:
        fh.write(json.dumps(record, sort_keys=True) + "\n")


def _record(r: Run) -> dict[str, Any]:
    s = SCENARIOS[r.name]
    rec = {"scenario": r.name, "role": s.role, "plant": s.note, "R": r.R,
           "seed": _base_seed() + 10_000 * ORDER.index(r.name),
           "not_supported_rate": r.rate, "relation_95_rate": r.rate95,
           "not_evaluable": r.not_evaluable, "mean_B": r.mean_B, "mean_bins": r.mean_bins,
           "mean_nested_share": r.mean_nested, "seconds": round(r.seconds, 2)}  # fmt: skip
    if s.role == "gated":
        rec["bound"] = _bound(r.R)
    _log_recorded(rec)
    print(json.dumps(rec, sort_keys=True))
    return rec


def _R(full: bool) -> int:
    return R_FULL if full else R_FAST


# ------------------------------------------------------------------ tests


GATED = [n for n, s in SCENARIOS.items() if s.role == "gated"]
REPORTED = [n for n, s in SCENARIOS.items() if s.role != "gated"]


@pytest.mark.parametrize("name", GATED)
@pytest.mark.parametrize("full", MODES)
def test_h4_false_refutation_fwer(full: bool, name: str) -> None:
    """Under H4 (unbiased genome estimator) the share of replicates with verdict
    ``not_supported`` is at most 0.05 + k sqrt(0.05 · 0.95 / R) (§8.4 FWER)."""
    R = _R(full)
    r = run(name, R)
    _record(r)
    bound = _bound(R)
    ok = r.rate <= bound
    detail = (f"{name}: false not_supported {r.not_supported}/{R} = {r.rate:.4f} "
              f"(bound {bound:.4f}; unadjusted relation_95 rule {r.rate95:.4f})")  # fmt: skip
    res = C.CheckResult(ok, r.rate, bound, R, detail)
    seed = _base_seed() + 10_000 * ORDER.index(name)
    C.log_result(SUITE, f"h4_attack.{name}.fwer", ALPHA, res, seed, r.seconds)
    assert r.mean_B > 0, f"{name}: no reportable bin, the check is vacuous"
    assert ok, detail


@pytest.mark.parametrize("name", REPORTED)
@pytest.mark.parametrize("full", MODES)
def test_h4_reported_rates(full: bool, name: str) -> None:
    """``clustered``: the false-refutation rate under dependent rows (S3), reported, not gated.
    ``power_*``: the share of replicates refuting a biased estimator, reported (the plan expects
    > 0.8 for ``power_easy_reference`` at n = 500). The per-replicate invariants (never
    ``supported``, ``relation_95 ⊇ relation_adj``, a reportable bin exists) still hold."""
    r = run(name, _R(full))
    _record(r)
    assert r.mean_B > 0, f"{name}: no reportable bin"
