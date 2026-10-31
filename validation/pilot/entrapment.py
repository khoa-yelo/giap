"""A1's free entrapment estimate (GIPPA_SPEC.md §15.5 metric 8, §8.7); spec
docs/specs/M10-07-free-entrapment.md with its binding statistician decisions.

K-12 peptides absent from Nissle's translation are the entrapment (E), K-12 peptides present in
it the original targets (T). Wen et al. 2025 (doi:10.1038/s41592-025-02719-x), equation (1):
``fdp_combined = N_E (1 + 1/r) / (N_T + N_E)``; equation (2): ``fdp_lower = N_E / (N_T + N_E)``.
Rule 1 of ``present`` is GIPPA's exact route, so it credits an M at the non-ATG starts that
``--annotation`` (a Bakta JSON, possibly gzipped) annotates on ``--genome`` (ADR-0068). Pure
functions plus a thin CLI::

    python validation/pilot/entrapment.py compute --pep-a1 PEP_A1 --k12-fasta K12_FASTA \\
        --contaminants CONTAMINANTS --genome GENOME --strain-db STRAIN_DB \\
        --annotation ANNOTATION --runs-manifest RUNS_MANIFEST --out METRICS_DIR \\
        [--out-a1 OUT_A1] [--diann-cmd command_step1.txt] [--provenance provenance.json]

Exit codes: 0 written; 1 cross-check 1 failed (two exact matchers disagree); 2 a precondition
(a missing, unreadable or malformed input, non-wild-type run, FASTA hash, topology tags, DIA-NN
settings, git unavailable, ``--out`` inside the checkout unless ``git check-ignore`` says it is
ignored); 4 an unexpected exception (its traceback goes to stderr). Nothing is written unless
the exit code is 0. File helpers, ``Precondition`` and ``Inconsistent`` are M10-05's
(``metrics.py``).
"""

from __future__ import annotations

import argparse
import sys
import traceback
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import asdict, dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

import ahocorasick
import numpy as np
import pandas as pd
from scipy.stats import binom
from validation.pilot import metrics as m5
from validation.pilot.metrics import Inconsistent, Precondition, sha256

from gippa.annotate.bakta import parse_bakta_json
from gippa.config import resolve_config
from gippa.config.model import Enzyme, TolerantConfig
from gippa.io.fasta import read_genome_fasta
from gippa.io.peptides import read_peptide_table
from gippa.io.protein_fasta import read_protein_fasta
from gippa.io.tsv import write_tsv
from gippa.mapping.decoys import DecoyType
from gippa.mapping.digest import digest
from gippa.mapping.exact import match_exact
from gippa.mapping.frameshift import match_frameshift
from gippa.mapping.initiator import StartCells, annotated_starts, start_frame_cells
from gippa.mapping.masses import MOD_DELTA, RESIDUE_MASS
from gippa.mapping.sixframe import Frame, SixFrameIndex, build_sixframe
from gippa.mapping.tolerant import match_tolerant
from gippa.model.coords import Replicon
from gippa.model.genes import GeneModel

SPEC = "docs/specs/M10-07-free-entrapment.md"
SCHEMA = "gippa-pilot-entrapment/1"
ALPHA = Fraction(1, 100)
WATER, PROTON = 18.010565, 1.007276
STANDARD = frozenset("ACDEFGHIKLMNPQRSTVWY")
BROAD = frozenset({"block", "swap", "near_block"})  # S2 sensitivity ``present_broad``
CAVEATS = (
    "near-homolog entrapment: the estimate errs high (§15.5)",
    "an estimated upper bound in expectation (Wen 2025, Suppl. Note 2.2), not a bound on this "
    "sample",
    "assumes the Nissle genome is complete, is the sample's strain, and that no K-12 material is "
    "in the sample",
    "assessment is made at precursor level, the unit DIA-NN's global q controls",
)
NOT_ASSESSED = {
    "peptide": "peptide-level FDR is not controlled by the search (q_value is DIA-NN's "
    "precursor-level Global.Q.Value)",
    "per_run": "runs were selected by the global and the run-level q-value together; no FDR "
    "claim attaches to a run's intersection",
}
CI_METHOD = "binomial_inversion_equal_chance"
CI_ASSUMES = (
    "equal chance and independent false matches; not valid under clustered (protein-level) errors"
)
# The entrapment needs target rows only: an empty decoy space makes the matchers' decoy pass a
# no-op (every query still needs a decoy type).
_NO_DECOY: DecoyType = "shuffle"
_NO_DECOYS: dict[DecoyType, tuple[Frame, ...]] = {_NO_DECOY: ()}


# ---------------------------------------------------------------- the searchable space


@dataclass(frozen=True)
class Settings:
    """The A1 digest and precursor space (``diann_A1.cfg``); ``charges``/``mz`` inclusive."""

    enzyme: Enzyme = "trypsin_p"
    missed: int = 1
    min_len: int = 7
    max_len: int = 30
    charges: tuple[int, int] = (1, 4)
    mz: tuple[float, float] = (300.0, 1800.0)
    var_mods: int = 1


DEFAULT_SETTINGS = Settings()


def norm(s: str) -> str:
    return s.replace("I", "L")


def _digest(proteins: Iterable[str], settings: Settings) -> tuple[set[str], set[str]]:
    """(standard, non-standard) normalised digest products."""
    std: set[str] = set()
    other: set[str] = set()
    st = settings
    for p in (p for f in proteins
              for _, _, p in digest(f, st.enzyme, st.missed, st.min_len, st.max_len)):  # fmt: skip
        (std if STANDARD.issuperset(p) else other).add(norm(p))
    return std, other


def digest_set(proteins: Iterable[str], settings: Settings = DEFAULT_SETTINGS) -> set[str]:
    """D(F): the normalised products of standard residues only."""
    return _digest(proteins, settings)[0]


def _mass(peptide: str, oxidised: Iterable[int] = ()) -> float:
    ox = list(oxidised)
    if any(peptide[i] != "M" for i in ox):
        raise ValueError(f"oxidised position {ox} is not an M")
    m = sum(RESIDUE_MASS[c] for c in peptide) + WATER
    return m + peptide.count("C") * MOD_DELTA["Carbamidomethyl"] + len(ox) * MOD_DELTA["Oxidation"]


def mz(peptide: str, charge: int, oxidised_at: int | None = None) -> float:
    """m/z with carbamidomethyl on every C and Met oxidation at ``oxidised_at`` (0-based)."""
    m = _mass(peptide, () if oxidised_at is None else (oxidised_at,))
    return (m + charge * PROTON) / charge


def n_prec(peptide: str, settings: Settings = DEFAULT_SETTINGS) -> int:
    """Precursors of ``peptide``: every form with at most ``var_mods`` oxidised M, at every
    charge whose m/z lies inside the bounds."""
    ms = [i for i, c in enumerate(peptide) if c == "M"]
    (zlo, zhi), (lo, hi) = settings.charges, settings.mz
    n = 0
    for j in range(min(settings.var_mods, len(ms)) + 1):
        for ox in combinations(ms, j):
            m = _mass(peptide, ox)
            n += sum(lo <= (m + z * PROTON) / z <= hi for z in range(zlo, zhi + 1))
    return n


# ---------------------------------------------------------------- present (one function)


@dataclass(frozen=True)
class Index:
    """Nissle's six-frame translation, the ``strain_db`` text, the tolerant settings and the
    annotated non-ATG start cells (no decoy cells: the decoy space is empty)."""

    sixframe: SixFrameIndex
    strain_text: str | None
    tolerant: TolerantConfig
    starts: StartCells


@dataclass(frozen=True)
class Presence:
    """``rule`` (exact, strain_db, mass_equivalent or None) and the ``diff_types`` of every
    tolerant target row (searched only for peptides without an exact rule-1 match)."""

    rule: str | None
    tolerant: tuple[frozenset[str], ...] = ()

    @property
    def broad(self) -> bool:
        return self.rule is not None or any(t and t <= BROAD for t in self.tolerant)


def build_index(genome: Mapping[str, str], replicons: Mapping[str, Replicon],
                strain_db: Sequence[str] | None = None, tolerant: TolerantConfig | None = None,
                *, annotation: GeneModel) -> Index:  # fmt: skip
    """``strain_db`` None: L2 answered "no" (six-frame translation only). ``annotation`` gives
    the start cells, with the config's ``start_codons`` and ``translation_table``, as GIPPA's
    ``sixframe_index`` step builds them."""
    text = None if strain_db is None else "*".join(norm(p) for p in strain_db)
    config = resolve_config()
    cfg = tolerant if tolerant is not None else config.mapping.tolerant
    ann = config.annotation
    six = build_sixframe(genome, replicons, table=ann.translation_table)
    sites, _ = annotated_starts(annotation, genome, ann.start_codons, table=ann.translation_table)
    starts = StartCells(sites, start_frame_cells(six, sites), {_NO_DECOY: ()})
    return Index(six, text, cfg, starts)


def _substrings(peptides: Sequence[str], text: str) -> set[str]:
    if not peptides or not text:
        return set()
    auto = ahocorasick.Automaton(ahocorasick.STORE_ANY)
    for p in peptides:
        auto.add_word(p, p)
    auto.make_automaton()
    return {p for _, p in auto.iter(text)}


def present(peptides: Iterable[str], index: Index) -> dict[str, Presence]:
    """Rule 1 (exact six-frame match, with an M credited at an annotated non-ATG start, or a
    ``strain_db`` substring; I -> L) else rule 2 (a tolerant target row whose ``diff_types`` is
    non-empty and only ``block``)."""
    peps = sorted(set(peptides))
    q = {p: norm(p) for p in peps}
    texts = sorted(set(q.values()))
    decoy_of = dict.fromkeys(texts, _NO_DECOY)
    hx = match_exact(texts, index.sixframe, _NO_DECOYS, decoy_of, min_len=1, max_len=10**6,
                     starts=index.starts)  # fmt: skip
    exact = set(hx.loc[hx["search"] == "target", "peptide"])
    rest = [t for t in texts if t not in exact]
    sdb = _substrings(rest, index.strain_text or "")
    left = [t for t in rest if t not in sdb]
    rows: dict[str, list[frozenset[str]]] = {}
    # frameshift rows too (GIPPA's tolerant route); they only move the ``far`` stratum (N1)
    for match in (match_tolerant, match_frameshift) if left and index.tolerant.enabled else ():
        part = match(left, index.sixframe, _NO_DECOYS, dict.fromkeys(left, _NO_DECOY), set(left),
                     index.tolerant)  # fmt: skip
        tgt = part[part["search"] == "target"]
        for p, d in zip(tgt["peptide"], tgt["diff_types"], strict=True):
            rows.setdefault(str(p), []).append(frozenset(x for x in str(d).split(";") if x))

    def presence(t: str) -> Presence:
        if t in exact or t in sdb:
            return Presence("exact" if t in exact else "strain_db")
        types = tuple(sorted(set(rows.get(t, [])), key=sorted))
        return Presence("mass_equivalent" if any(x and x <= {"block"} for x in types) else None,
                        types)  # fmt: skip

    return {p: presence(t) for p, t in q.items()}


# ---------------------------------------------------------------- the estimator (pure API)


def _defined(n_t: int, n_e: int, k: int, s: int) -> bool:
    if n_e and s and not k:
        raise ValueError("an entrapment discovery with no entrapment space (|k12_only| = 0)")
    return n_t + n_e > 0 and s > 0 and k > 0


def combined(n_t: int, n_e: int, n_k12_only: int, n_shared: int) -> Fraction | None:
    """Equation (1), N_E (1 + 1/r) / (N_T + N_E) with r = k/s, exactly; unclipped."""
    if not _defined(n_t, n_e, n_k12_only, n_shared):
        return None
    return Fraction(n_e * (n_k12_only + n_shared), n_k12_only * (n_t + n_e))


def lower(n_t: int, n_e: int) -> Fraction | None:
    """Equation (2), N_E / (N_T + N_E)."""
    return None if n_t + n_e == 0 else Fraction(n_e, n_t + n_e)


def ucb95(n_t: int, n_e: int, n_k12_only: int, n_shared: int, level: float = 0.05) -> float | None:
    """S3: ``max{M in N_E..N : BinomCDF(N_E; M, pi) > level} / N``, pi = k/(k+s)."""
    if not _defined(n_t, n_e, n_k12_only, n_shared):
        return None
    n = n_t + n_e
    m = np.arange(n_e, n + 1)
    ok = m[binom.cdf(n_e, m, n_k12_only / (n_k12_only + n_shared)) > level]
    return float(ok.max()) / n


def assess(combined: Fraction | None, lower: Fraction | None, alpha: Fraction = ALPHA) -> str:
    """The paper's readings (S6): lower > alpha first, then combined <= alpha."""
    if lower is not None and lower > alpha:
        return "evidence_of_failure"
    if combined is None:
        return "not_evaluable"
    return "consistent_with_control" if combined <= alpha else "inconclusive"


Fracs = tuple[dict[str, Any], Fraction | None, Fraction | None]


def _estimates(n_t: int, n_e: int, k: int, s: int) -> Fracs:
    """(the float block, combined, lower); the reason precedence is no shared peptide, no
    discovery, no entrapment space."""
    reason = "no_shared_peptides" if s == 0 else "empty_discoveries" if n_t + n_e == 0 else None
    c = lo = None
    if reason is None:
        c, lo = combined(n_t, n_e, k, s), lower(n_t, n_e)
        reason = "no_entrapment_space" if c is None else None
    out: dict[str, Any] = {"fdp_combined": None if c is None else float(c),
                           "fdp_lower": None if lo is None else float(lo),
                           "ucb95": ucb95(n_t, n_e, k, s) if s and n_t + n_e else None}  # fmt: skip
    return out | ({"reason": reason} if reason else {}), c, lo


def _ratio(a: int | float, b: int | float) -> float | None:
    return None if b == 0 else a / b


def _level(n_t: int, n_e: int, k: int, s: int, other: Mapping[str, int],
           scope: str) -> dict[str, Any]:  # fmt: skip
    """One level's block; ``scope`` is ``precursor`` (assessed), ``peptide`` or ``per_run``."""
    est, c, lo = _estimates(n_t, n_e, k, s)
    ci = {"ucb95": est.pop("ucb95"), "method": CI_METHOD, "assumes": CI_ASSUMES}
    b: dict[str, Any] = {"N_T": n_t, "N_E": n_e, **other, **est, "ci": ci}
    b["above_one"] = c is not None and c > 1
    if scope == "precursor":
        b["assessment"] = assess(c, lo)
    else:
        b["assessment"], b["assessment_reason"] = "not_assessed", NOT_ASSESSED[scope]
    return b


# ---------------------------------------------------------------- evaluate


def _accepted(rows: pd.DataFrame) -> pd.DataFrame:
    """M10-05's rule: a target row (NA ``is_decoy`` reads as target) with q <= Q_ACCEPT."""
    q = m5._q(rows)
    return rows[m5._target(rows) & ~np.isnan(q) & (q <= m5.Q_ACCEPT)]


def evaluate(k12: Sequence[str], contaminants: Sequence[str], index: Index, rows: pd.DataFrame,
             settings: Settings = DEFAULT_SETTINGS) -> dict[str, Any]:  # fmt: skip
    """The ``entrapment.json`` blocks (``universe``, ``peptide``, ``precursor``, ``per_run``,
    ``warnings``) and ``discoveries`` (peptide, bucket, present_rule) for the TSV."""
    k12_all, k12_ns = _digest(k12, settings)
    cont = digest_set(contaminants, settings)
    nprec = {p: n_prec(p, settings) for p in sorted(k12_all - cont)}
    uni = sorted(p for p, n in nprec.items() if n > 0)
    pres_u = present(uni, index)
    shared = [p for p in uni if pres_u[p].rule]
    only = [p for p in uni if not pres_u[p].rule]
    S, K = sum(nprec[p] for p in shared), sum(nprec[p] for p in only)
    rules = [pres_u[p].rule for p in shared]
    universe = {"k12_peptides": len(k12_all), "excluded_contaminant": len(k12_all & cont),
                "non_standard": len(k12_ns), "shared": len(shared), "k12_only": len(only),
                "shared_exact": rules.count("exact"),
                "shared_strain_db_only": rules.count("strain_db"),
                "shared_mass_equivalent": rules.count("mass_equivalent"),
                "r": _ratio(len(only), len(shared)), "shared_precursors": S,
                "k12_only_precursors": K, "r_prec": _ratio(K, S)}  # fmt: skip

    acc = _accepted(rows)
    peps = [norm(str(p)) for p in acc["peptide"]]
    in_u = set(uni)
    pres_d = present(sorted({p for p in peps if p in in_u and p not in cont}), index)

    def bucket(p: str) -> str:
        if p in cont:
            return "excluded_contaminant"
        if p not in in_u:
            return "not_in_universe"
        return "target" if pres_d[p].rule else "entrapment"

    def counts(keys: Iterable[str]) -> tuple[int, int, dict[str, int]]:
        bs = [bucket(p) for p in keys]
        other = {b: bs.count(b) for b in ("excluded_contaminant", "not_in_universe")}
        return bs.count("target"), bs.count("entrapment"), other

    k, s = len(only), len(shared)
    n_t, n_e, oth = counts(set(peps))
    peptide = _level(n_t, n_e, k, s, oth, "peptide")
    mods = zip(acc["modified_peptide"], acc["charge"], peps, strict=True)
    pt, pe, po = counts(p for _, _, p in {(norm(str(m)), int(z), p) for m, z, p in mods})
    precursor = _level(pt, pe, K, S, po, "precursor")
    runs = [str(r) for r in acc["run"]]
    per_run = {}
    for run in sorted({str(r) for r in rows["run"]}):
        rt, re_, ro = counts({p for p, r in zip(peps, runs, strict=True) if r == run})
        per_run[run] = _level(rt, re_, k, s, ro, "per_run")

    disc = sorted(set(peps))
    found = {p for p in disc if bucket(p) in ("target", "entrapment")}
    # S2: present_broad
    bs = [p for p in uni if pres_u[p].broad]
    bt = sum(pres_d[p].broad for p in found)
    bk, be = len(uni) - len(bs), len(found) - bt
    broad = {"shared": len(bs), "k12_only": bk, "r": _ratio(bk, len(bs)), "N_T": bt, "N_E": be,
             **_estimates(bt, be, bk, len(bs))[0]}  # fmt: skip
    # S4: r_prec weighted
    rp = _estimates(n_t, n_e, K, S)[0]
    weighted = {"fdp_combined_rprec": rp["fdp_combined"], "mean_nprec_shared": _ratio(S, s),
                "mean_nprec_k12_only": _ratio(K, k)} | _reason(rp)  # fmt: skip
    # S4: far stratum
    min10 = index.tolerant.min_len

    def stratum(p: str) -> str:
        return "lt10" if len(p) < min10 else "near_homolog" if pres_u[p].tolerant else "far"

    e_found = [p for p in found if not pres_d[p].rule]
    e_far = [p for p in only if stratum(p) == "far"]
    t10 = [p for p in shared if len(p) >= min10]
    nt10 = sum(len(p) >= min10 for p in found if pres_d[p].rule)
    nef = sum(stratum(p) == "far" for p in e_found)
    est = _estimates(nt10, nef, len(e_far), len(t10))[0]
    breakdown = {"k12_only": _strata([stratum(p) for p in only]),
                 "N_E": _strata([stratum(p) for p in e_found])}  # fmt: skip
    far = {"r_far": _ratio(len(e_far), len(t10)), "N_T_10": nt10, "N_E_far": nef,
           "fdp_combined_far": est["fdp_combined"], "fdp_lower_far": est["fdp_lower"],
           "ucb95_far": est["ucb95"], "breakdown": breakdown} | _reason(est)  # fmt: skip
    peptide["sensitivity"] = {"present_broad": broad, "r_prec_weighted": weighted, "far": far}

    warnings: list[Any] = []
    if oth["not_in_universe"]:
        warnings.append({"code": "digest_mismatch", "count": oth["not_in_universe"]})
    if n_t + n_e and not n_e:
        warnings.append({"code": "n_entrapment_zero",
                         "message": "no K-12-only peptide was accepted"})  # fmt: skip
    discoveries = [(p, bucket(p), (pres_d[p].rule or "") if p in found else "") for p in disc]
    return {"universe": universe, "peptide": peptide, "precursor": precursor, "per_run": per_run,
            "warnings": warnings, "discoveries": discoveries}  # fmt: skip


def _strata(labels: Sequence[str]) -> dict[str, int]:
    return {x: labels.count(x) for x in ("lt10", "near_homolog", "far")}


def _reason(est: Mapping[str, Any]) -> dict[str, Any]:
    return {"reason": est["reason"]} if "reason" in est else {}


# ---------------------------------------------------------------- cross-checks and settings


def parse_diann_cmd(tokens: Sequence[str]) -> dict[str, list[str]]:
    """Flag -> values of a one-token-per-line DIA-NN command (a bare flag has ``[]``)."""
    out: dict[str, list[str]] = {}
    flag = None
    for t in tokens:
        if t.startswith("--"):
            flag = t
            out.setdefault(flag, [])
        elif flag is not None:
            out[flag].append(t)
            flag = None
    return out


def settings_mismatch(flags: Mapping[str, list[str]], settings: Settings) -> list[str]:
    """Differences between the search's flags and ``settings`` (cross-check 4)."""
    cut = {"trypsin_p": "K*,R*", "trypsin": "K*,R*,!*P"}[settings.enzyme]
    want: dict[str, str | float] = {
        "--cut": cut, "--missed-cleavages": settings.missed, "--min-pep-len": settings.min_len,
        "--max-pep-len": settings.max_len, "--min-pr-charge": settings.charges[0],
        "--max-pr-charge": settings.charges[1], "--min-pr-mz": settings.mz[0],
        "--max-pr-mz": settings.mz[1], "--var-mods": settings.var_mods,
        "--var-mod": "UniMod:35,15.994915,M"}  # fmt: skip
    bad = []
    for f, v in want.items():
        got = flags.get(f, [])
        same = len(got) == 1 and (got[0] == v if isinstance(v, str) else _num(got[0]) == float(v))
        if not same:
            bad.append(f"{f} {got} != {v}")
    if "--unimod4" not in flags:
        bad.append("--unimod4 missing (carbamidomethyl C)")
    if "--met-excision" in flags:
        bad.append("--met-excision set (GIPPA's digest keeps the initiator M)")
    return bad


def _num(s: str) -> float | None:
    try:
        return float(s)
    except ValueError:
        return None


CLASSIFIED = ("target", "entrapment")


Discovery = tuple[str, str, str]  # (peptide, bucket, present_rule)


def cross_check_run(discoveries: Sequence[Discovery], out_a1: Path, n_e: int,
                    ) -> tuple[dict[str, Any], list[dict[str, Any]]]:  # fmt: skip
    """Cross-checks 1 and 2 against the A1 run directory. 1 looks at the T and E buckets only
    and raises ``Inconsistent`` when an A1 ``route = exact`` peptide is not ``exact`` here (the
    exact route is six-frame only, so ``strain_db`` is no agreement)."""
    hits = m5._run_table(out_a1, "peptide_hits.tsv")
    if hits is None:
        raise Precondition(f"{out_a1}: no peptide_hits.tsv under tables/ or mapping/")
    routes = m5._routes(m5._table(hits, "peptide_hits"))
    exact = {norm(p) for p, r in routes.items() if "exact" in r}
    rule = {p: r for p, b, r in discoveries if b in CLASSIFIED}
    if bad := sorted(p for p in exact & set(rule) if rule[p] != "exact"):
        raise Inconsistent(f"{len(bad)} discovered peptide(s) have a route = exact row in "
                           f"{hits} but no exact match here")  # fmt: skip
    ds = m5._run_table(out_a1, "detection_summary.tsv")
    gi = None if ds is None else m5.genome_inconsistent(m5._table(ds, "detection_summary"))
    count = None if gi is None else gi["count"]
    targets = [(p, r) for p, b, r in discoveries if b == "target"]
    unclassified = sum(p in exact and b not in CLASSIFIED for p, b, _ in discoveries)
    check = {"exact_route_agrees": True, "n_exact_route": len(exact & set(rule)),
             "n_exact_route_unclassified": unclassified, "genome_inconsistent_count": count,
             "n_entrapment": n_e,
             "present_not_exact": sum(r in ("strain_db", "mass_equivalent") for _, r in targets),
             "exact_rejected": sum(r == "exact" and p not in exact
                                   for p, r in targets)}  # fmt: skip
    if count is None:
        return check, [{"code": "genome_inconsistent_unavailable", "out_a1": str(out_a1)}]
    if n_e > count:
        return check, [{"code": "n_entrapment_above_genome_inconsistent", "n_entrapment": n_e,
                        "genome_inconsistent_count": count}]  # fmt: skip
    return check, []


# ---------------------------------------------------------------- CLI (file I/O)

INPUTS = ("pep_a1", "k12_fasta", "contaminants", "genome", "strain_db", "annotation",
          "runs_manifest")  # fmt: skip
OPTIONAL = ("out_a1", "diann_cmd", "provenance")
#: provenance.json ``inputs.fasta[]`` basename -> the argument it vouches for (S2: by role)
PROVENANCE_ROLES = {"UP000000625.fasta": "k12_fasta", "contaminants.fasta": "contaminants"}
PEP_COLUMNS = ("peptide", "modified_peptide", "charge", "run", "q_value")


def check_provenance(prov: Any, path: Path, files: Mapping[str, Path]) -> None:
    """Each role's record (by basename) holds exactly the sha256 of its argument."""
    inputs = prov.get("inputs") if isinstance(prov, dict) else None
    recs = inputs.get("fasta") if isinstance(inputs, dict) else None
    if not isinstance(recs, list):
        raise Precondition(f"{path}: no inputs.fasta list")
    by: dict[str, set[Any]] = {}
    for r in recs:
        if isinstance(r, dict):
            by.setdefault(Path(str(r.get("path", ""))).name, set()).add(r.get("sha256"))
    for name, dest in PROVENANCE_ROLES.items():
        if by.get(name) != {sha256(files[dest])}:
            raise Precondition(f"{files[dest]}: sha256 is not the one {path} records for {name}")


def _proteins(path: Path) -> list[str]:
    return m5._load(path, lambda p: [r.sequence for r in read_protein_fasta(p)])


def _index(genome: Path, strain_db: Path, annotation: Path) -> Index:
    recs = m5._load(genome, read_genome_fasta)
    if not any(r.circular_tag for r in recs):
        raise Precondition(f"{genome}: no topology tag on any contig")
    reps = {r.id: Replicon(r.id, len(r.sequence), "circular" if r.circular_tag else "linear")
            for r in recs}  # fmt: skip
    seqs = {r.id: r.sequence for r in recs}
    table = resolve_config().annotation.translation_table

    def model(path: Path) -> GeneModel:
        return parse_bakta_json(path, reps, seqs, translation_table=table).model

    return build_index(seqs, reps, _proteins(strain_db), annotation=m5._load(annotation, model))


def _pep_rows(path: Path, manifest: Path) -> pd.DataFrame:
    rows = m5._load(path, lambda p: read_peptide_table(p)[0])
    if missing := [c for c in PEP_COLUMNS if c not in rows.columns]:
        raise Precondition(f"{path}: missing column(s) {', '.join(missing)}")
    if _accepted(rows)["charge"].isna().any():
        raise Precondition(f"{path}: NA charge on an accepted row")
    wt = m5._rows(manifest)
    if not wt or "run" not in wt[0]:
        raise Precondition(f"{manifest}: no run column")
    if foreign := sorted({str(r) for r in rows["run"]} - {str(r["run"]) for r in wt}):
        raise Precondition(f"{path}: run(s) not in {manifest}: {', '.join(foreign)}")
    return rows


def compute(a: argparse.Namespace) -> tuple[dict[str, Any], list[Discovery]]:
    """``entrapment.json`` and the discovery rows from the parsed arguments (reads only)."""
    m5._outside_repo(a.out, "--out", allow_ignored=True)  # B2: ``git check-ignore``
    files = {k: getattr(a, k) for k in (*INPUTS, *OPTIONAL) if getattr(a, k) is not None}
    for p in files.values():
        if not p.exists():
            raise Precondition(f"missing input: {p}")
    files.pop("out_a1", None)
    if a.provenance is not None:
        check_provenance(m5._json(a.provenance), a.provenance, files)
    settings = DEFAULT_SETTINGS
    if a.diann_cmd is not None:
        lines = m5._load(a.diann_cmd, lambda p: p.read_text(encoding="utf-8").split())
        if bad := settings_mismatch(parse_diann_cmd(lines), settings):
            raise Precondition(f"{a.diann_cmd}: settings differ: {'; '.join(bad)}")
    rows = _pep_rows(a.pep_a1, a.runs_manifest)
    if not (k12 := _proteins(a.k12_fasta)):
        raise Precondition(f"{a.k12_fasta}: no protein")
    index = _index(a.genome, a.strain_db, a.annotation)
    doc = evaluate(k12, _proteins(a.contaminants), index, rows, settings)
    discoveries = doc.pop("discoveries")
    nu = doc["peptide"]["not_in_universe"]
    checks: dict[str, Any] = {"digest": {"not_in_universe": nu, "agrees": nu == 0}}
    if a.diann_cmd is not None:
        checks["settings"] = {"agrees": True}
    if a.out_a1 is not None:
        out_a1 = m5.resolve_attempt(a.out_a1)[0]
        checks["run"], warn = cross_check_run(discoveries, out_a1, doc["peptide"]["N_E"])
        doc["warnings"] += warn
    doc |= {"schema": SCHEMA, "spec": m5._spec_sha(SPEC), "cross_checks": checks,
            "caveats": list(CAVEATS),
            "inputs": {k: {"path": str(p.resolve()), "sha256": sha256(p)}
                       for k, p in files.items()},
            "settings": asdict(settings) | {"tolerant": index.tolerant.model_dump()}}  # fmt: skip
    return doc, discoveries


def write_outputs(out: Path, doc: Mapping[str, Any], discoveries: Sequence[Sequence[str]]) -> None:
    """The TSV first, then the JSON (written last, so a failed TSV leaves no JSON behind)."""
    out.mkdir(parents=True, exist_ok=True)
    write_tsv(out / "entrapment_discoveries.tsv", ["peptide", "bucket", "present_rule"],
              discoveries)  # fmt: skip
    (out / "entrapment.json").write_text(m5._dump(doc), encoding="utf-8")


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute")
    for dest in (*INPUTS, "out", *OPTIONAL):
        c.add_argument(f"--{dest.replace('_', '-')}", type=Path, required=dest not in OPTIONAL)
    a = ap.parse_args(argv)
    try:
        doc, discoveries = compute(a)
        write_outputs(a.out, doc, discoveries)
    except Precondition as err:
        print(f"precondition: {err}", file=sys.stderr)
        return 2
    except Inconsistent as err:
        print(f"inconsistent: {err}", file=sys.stderr)
        return 1
    except Exception:  # exit 1 is cross-check 1's: a bug exits 4 with its traceback
        traceback.print_exc()
        return 4
    p = doc["precursor"]
    print(f"wrote {a.out}/entrapment.json; precursor {p['assessment']}; "
          f"{len(doc['warnings'])} warning(s)")  # fmt: skip
    return 0


if __name__ == "__main__":
    sys.exit(main())
