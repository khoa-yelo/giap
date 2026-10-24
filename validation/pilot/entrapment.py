"""A1's free entrapment estimate (GIAP_SPEC.md §15.5 metric 8, §8.7); spec
docs/specs/M10-07-free-entrapment.md with its binding statistician decisions.

K-12 peptides absent from Nissle's translation are the entrapment (E), K-12 peptides present in
it the original targets (T). Wen et al. 2025 (doi:10.1038/s41592-025-02719-x), equation (1):
``fdp_combined = N_E (1 + 1/r) / (N_T + N_E)``; equation (2): ``fdp_lower = N_E / (N_T + N_E)``.
Pure functions plus a thin CLI::

    python validation/pilot/entrapment.py compute --pep-a1 PEP_A1 --k12-fasta K12_FASTA \\
        --contaminants CONTAMINANTS --genome GENOME --strain-db STRAIN_DB \\
        --runs-manifest RUNS_MANIFEST --out METRICS_DIR \\
        [--out-a1 OUT_A1] [--diann-cmd command_step1.txt] [--provenance provenance.json]

Exit codes: 0 written; 1 cross-check 1 failed (two exact matchers disagree); 2 a precondition
(missing input, non-wild-type run, FASTA hash, topology tags, DIA-NN settings, ``--out`` inside
the checkout). Nothing is written unless the exit code is 0.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Any

import ahocorasick
import numpy as np
import pandas as pd
from scipy.stats import binom

from giap.config import resolve_config
from giap.config.model import Enzyme, TolerantConfig
from giap.mapping.decoys import DecoyType
from giap.mapping.digest import digest
from giap.mapping.exact import match_exact
from giap.mapping.frameshift import match_frameshift
from giap.mapping.masses import MOD_DELTA, RESIDUE_MASS
from giap.mapping.sixframe import Frame, SixFrameIndex, build_sixframe
from giap.mapping.tolerant import match_tolerant
from giap.model.coords import Replicon

REPO = Path(__file__).resolve().parents[2]
SPEC = "docs/specs/M10-07-free-entrapment.md"
SCHEMA = "giap-pilot-entrapment/1"
ALPHA = Fraction(1, 100)
Q_ACCEPT = 0.01  # accepted: q_value <= 0.01 (inclusive, the written value); NA not accepted
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
    "per_run": "runs were selected by the global q-value, not a run-level q-value",
}
CI_METHOD = "binomial_inversion_equal_chance"
CI_ASSUMES = (
    "equal chance and independent false matches; not valid under clustered (protein-level) errors"
)
# The entrapment needs target rows only: an empty decoy space makes the matchers' decoy pass a
# no-op (every query still needs a decoy type).
_NO_DECOY: DecoyType = "shuffle"
_NO_DECOYS: dict[DecoyType, tuple[Frame, ...]] = {_NO_DECOY: ()}


class Precondition(Exception):
    """Exit 2: nothing is written."""


class Inconsistent(Exception):
    """Exit 1: two exact matchers disagree on one genome (cross-check 1)."""


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
    for prot in proteins:
        for _, _, p in digest(
            prot, settings.enzyme, settings.missed, settings.min_len, settings.max_len
        ):
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
    """Nissle's six-frame translation, the ``strain_db`` text and the tolerant settings."""

    sixframe: SixFrameIndex
    strain_text: str | None
    tolerant: TolerantConfig


@dataclass(frozen=True)
class Presence:
    """``rule`` (exact, strain_db, mass_equivalent or None) and the ``diff_types`` of every
    tolerant target row (searched only for peptides without an exact rule-1 match)."""

    rule: str | None
    tolerant: tuple[frozenset[str], ...] = ()

    @property
    def broad(self) -> bool:
        return self.rule is not None or any(t and t <= BROAD for t in self.tolerant)


def build_index(
    genome: Mapping[str, str],
    replicons: Mapping[str, Replicon],
    strain_db: Sequence[str] | None = None,
    tolerant: TolerantConfig | None = None,
) -> Index:
    """``strain_db`` None: L2 answered "no" (six-frame translation only)."""
    text = None if strain_db is None else "*".join(norm(p) for p in strain_db)
    cfg = tolerant if tolerant is not None else resolve_config().mapping.tolerant
    return Index(build_sixframe(genome, replicons), text, cfg)


def _substrings(peptides: Sequence[str], text: str) -> set[str]:
    if not peptides or not text:
        return set()
    auto = ahocorasick.Automaton(ahocorasick.STORE_ANY)
    for p in peptides:
        auto.add_word(p, p)
    auto.make_automaton()
    return {p for _, p in auto.iter(text)}


def present(peptides: Iterable[str], index: Index) -> dict[str, Presence]:
    """Rule 1 (exact six-frame match, or a ``strain_db`` substring; I -> L) else rule 2 (a
    tolerant target row whose ``diff_types`` is non-empty and only ``block``)."""
    peps = sorted(set(peptides))
    q = {p: norm(p) for p in peps}
    texts = sorted(set(q.values()))
    decoy_of = dict.fromkeys(texts, _NO_DECOY)
    hx = match_exact(texts, index.sixframe, _NO_DECOYS, decoy_of, min_len=1, max_len=10**6)
    exact = set(hx.loc[hx["search"] == "target", "peptide"])
    rest = [t for t in texts if t not in exact]
    sdb = _substrings(rest, index.strain_text or "")
    left = [t for t in rest if t not in sdb]
    rows: dict[str, list[frozenset[str]]] = {}
    if left and index.tolerant.enabled:
        dec = {t: _NO_DECOY for t in left}
        cfg = index.tolerant
        parts = [
            match_tolerant(left, index.sixframe, _NO_DECOYS, dec, set(left), cfg),
            match_frameshift(left, index.sixframe, _NO_DECOYS, dec, set(left), cfg),
        ]
        for part in parts:
            tgt = part[part["search"] == "target"]
            for p, d in zip(tgt["peptide"], tgt["diff_types"], strict=True):
                rows.setdefault(str(p), []).append(frozenset(x for x in str(d).split(";") if x))
    out: dict[str, Presence] = {}
    for p, t in q.items():
        if t in exact:
            out[p] = Presence("exact")
        elif t in sdb:
            out[p] = Presence("strain_db")
        else:
            types = tuple(sorted(set(rows.get(t, [])), key=sorted))
            me = any(x and x <= {"block"} for x in types)
            out[p] = Presence("mass_equivalent" if me else None, types)
    return out


# ---------------------------------------------------------------- the estimator (pure API)


def _defined(n_t: int, n_e: int, k: int, s: int) -> bool:
    if n_t + n_e == 0 or s == 0:
        return False
    if k == 0:
        if n_e:
            raise ValueError("an entrapment discovery with no entrapment space (|k12_only| = 0)")
        return False
    return True


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


def _estimates(n_t: int, n_e: int, k: int, s: int) -> dict[str, Any]:
    if s == 0 or n_t + n_e == 0:
        reason = "no_shared_peptides" if s == 0 else "empty_discoveries"
        return {"fdp_combined": None, "fdp_lower": None, "ucb95": None, "reason": reason}
    c, lo = combined(n_t, n_e, k, s), lower(n_t, n_e)
    out: dict[str, Any] = {
        "fdp_combined": None if c is None else float(c),
        "fdp_lower": float(lo or 0),
        "ucb95": ucb95(n_t, n_e, k, s),
    }
    if c is None:
        out["reason"] = "no_entrapment_space"
    return out


def _ratio(a: int | float, b: int | float) -> float | None:
    return None if b == 0 else a / b


def _level(
    n_t: int, n_e: int, k: int, s: int, other: Mapping[str, int], scope: str
) -> dict[str, Any]:
    """One level's block; ``scope`` is ``precursor`` (assessed), ``peptide`` or ``per_run``."""
    est = _estimates(n_t, n_e, k, s)
    ub = est.pop("ucb95")
    b: dict[str, Any] = {
        "N_T": n_t,
        "N_E": n_e,
        **other,
        **est,
        "ci": {"ucb95": ub, "method": CI_METHOD, "assumes": CI_ASSUMES},
    }
    b["above_one"] = est["fdp_combined"] is not None and est["fdp_combined"] > 1
    if scope == "precursor":
        defined = est["fdp_combined"] is not None
        b["assessment"] = assess(
            combined(n_t, n_e, k, s) if defined else None,
            lower(n_t, n_e) if est["fdp_lower"] is not None else None,
        )
    else:
        b["assessment"], b["assessment_reason"] = "not_assessed", NOT_ASSESSED[scope]
    return b


# ---------------------------------------------------------------- evaluate


def _accepted(rows: pd.DataFrame) -> pd.DataFrame:
    q = pd.to_numeric(rows["q_value"], errors="coerce")
    decoy = rows["is_decoy"].astype("boolean").fillna(False).astype(bool)
    return rows[(~decoy & q.notna() & (q <= Q_ACCEPT)).to_numpy()]


def evaluate(
    k12: Sequence[str],
    contaminants: Sequence[str],
    index: Index,
    rows: pd.DataFrame,
    settings: Settings = DEFAULT_SETTINGS,
) -> dict[str, Any]:
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
    universe = {
        "k12_peptides": len(k12_all),
        "excluded_contaminant": len(k12_all & cont),
        "non_standard": len(k12_ns),
        "shared": len(shared),
        "k12_only": len(only),
        "shared_exact": rules.count("exact"),
        "shared_strain_db_only": rules.count("strain_db"),
        "shared_mass_equivalent": rules.count("mass_equivalent"),
        "r": _ratio(len(only), len(shared)),
        "shared_precursors": S,
        "k12_only_precursors": K,
        "r_prec": _ratio(K, S),
    }

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
    prec = {
        (norm(str(m)), int(z), p)
        for m, z, p in zip(acc["modified_peptide"], acc["charge"], peps, strict=True)
    }
    pt, pe, po = counts(p for _, _, p in prec)
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
    broad = {
        "shared": len(bs),
        "k12_only": len(uni) - len(bs),
        "r": _ratio(len(uni) - len(bs), len(bs)),
        "N_T": bt,
        "N_E": len(found) - bt,
        **_estimates(bt, len(found) - bt, len(uni) - len(bs), len(bs)),
    }
    # S4: r_prec weighted
    rp = _estimates(n_t, n_e, K, S)
    weighted = {
        "fdp_combined_rprec": rp["fdp_combined"],
        "mean_nprec_shared": _ratio(S, s),
        "mean_nprec_k12_only": _ratio(K, k),
    }
    if "reason" in rp:
        weighted["reason"] = rp["reason"]
    # S4: far stratum
    min10 = index.tolerant.min_len

    def stratum(p: str) -> str:
        return "lt10" if len(p) < min10 else "near_homolog" if pres_u[p].tolerant else "far"

    e_found = [p for p in found if not pres_d[p].rule]
    e_far = [p for p in only if stratum(p) == "far"]
    t10 = [p for p in shared if len(p) >= min10]
    nt10 = sum(len(p) >= min10 for p in found if pres_d[p].rule)
    nef = sum(stratum(p) == "far" for p in e_found)
    est = _estimates(nt10, nef, len(e_far), len(t10))
    far = {
        "r_far": _ratio(len(e_far), len(t10)),
        "N_T_10": nt10,
        "N_E_far": nef,
        "fdp_combined_far": est["fdp_combined"],
        "fdp_lower_far": est["fdp_lower"],
        "ucb95_far": est["ucb95"],
        "breakdown": {
            "k12_only": _strata([stratum(p) for p in only]),
            "N_E": _strata([stratum(p) for p in e_found]),
        },
    }
    if "reason" in est:
        far["reason"] = est["reason"]
    peptide["sensitivity"] = {"present_broad": broad, "r_prec_weighted": weighted, "far": far}

    warnings: list[Any] = []
    if oth["not_in_universe"]:
        warnings.append({"code": "digest_mismatch", "count": oth["not_in_universe"]})
    if n_t + n_e and not n_e:
        warnings.append(
            {"code": "n_entrapment_zero", "message": "no K-12-only peptide was accepted"}
        )
    discoveries = [(p, bucket(p), (pres_d[p].rule or "") if p in found else "") for p in disc]
    return {
        "universe": universe,
        "peptide": peptide,
        "precursor": precursor,
        "per_run": per_run,
        "warnings": warnings,
        "discoveries": discoveries,
    }


def _strata(labels: Sequence[str]) -> dict[str, int]:
    return {x: labels.count(x) for x in ("lt10", "near_homolog", "far")}


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
        "--cut": cut,
        "--missed-cleavages": settings.missed,
        "--min-pep-len": settings.min_len,
        "--max-pep-len": settings.max_len,
        "--min-pr-charge": settings.charges[0],
        "--max-pr-charge": settings.charges[1],
        "--min-pr-mz": settings.mz[0],
        "--max-pr-mz": settings.mz[1],
        "--var-mods": settings.var_mods,
        "--var-mod": "UniMod:35,15.994915,M",
    }
    bad = []
    for f, v in want.items():
        got = flags.get(f, [])
        same = len(got) == 1 and (got[0] == v if isinstance(v, str) else _num(got[0]) == float(v))
        if not same:
            bad.append(f"{f} {got} != {v}")
    if "--unimod4" not in flags:
        bad.append("--unimod4 missing (carbamidomethyl C)")
    if "--met-excision" in flags:
        bad.append("--met-excision set (GIAP's digest keeps the initiator M)")
    return bad


def _num(s: str) -> float | None:
    try:
        return float(s)
    except ValueError:
        return None


def _table(path: Path) -> list[dict[str, str | None]]:
    from giap.io.tsv import read_tsv

    header, rows = read_tsv(path)
    return [dict(zip(header, r, strict=True)) for r in rows]


def cross_check_run(
    discoveries: Sequence[tuple[str, str, str]], out_a1: Path, n_e: int
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    """Cross-checks 1 and 2 against the A1 run directory (raises ``Inconsistent`` for 1)."""
    hits = next(
        (
            out_a1 / d / "peptide_hits.tsv"
            for d in ("tables", "mapping")
            if (out_a1 / d / "peptide_hits.tsv").is_file()
        ),
        None,
    )
    if hits is None:
        raise Precondition(f"{out_a1}: no tables/peptide_hits.tsv")
    exact = {norm(str(r["peptide"])) for r in _table(hits) if r.get("route") == "exact"}
    rule = {p: r for p, _, r in discoveries}
    bad = sorted(p for p in exact & set(rule) if rule[p] not in ("exact", "strain_db"))
    if bad:
        raise Inconsistent(
            f"{len(bad)} discovered peptide(s) have a route = exact row in "
            f"{hits} but are not present here"
        )
    ds = out_a1 / "tables" / "detection_summary.tsv"
    gi = [
        r
        for r in (_table(ds) if ds.is_file() else [])
        if r.get("section") == "qc" and r.get("metric") == "genome_inconsistent_rate"
    ]
    value, n = (gi[0]["value"], gi[0]["n"]) if gi else (None, None)
    count = None if value is None or n is None else round(float(value) * int(n))
    targets = [(p, r) for p, b, r in discoveries if b == "target"]
    check = {
        "exact_route_agrees": True,
        "n_exact_route": len(exact & set(rule)),
        "genome_inconsistent_count": count,
        "n_entrapment": n_e,
        "present_not_exact": sum(r in ("strain_db", "mass_equivalent") for _, r in targets),
        "exact_rejected": sum(r == "exact" and p not in exact for p, r in targets),
    }
    warn = []
    if count is None or n_e > count:
        warn.append(
            {
                "code": "n_entrapment_above_genome_inconsistent",
                "n_entrapment": n_e,
                "genome_inconsistent_count": count,
            }
        )
    return check, warn


# ---------------------------------------------------------------- CLI (file I/O)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _spec_sha() -> str:
    res = subprocess.run(
        ["git", "-C", str(REPO), "log", "-1", "--format=%H", "--", SPEC],
        capture_output=True,
        text=True,
        check=False,
    )
    sha = res.stdout.strip()
    return sha if re.fullmatch(r"[0-9a-f]{40}", sha) else "0" * 40


def compute(a: argparse.Namespace) -> tuple[dict[str, Any], list[tuple[str, str, str]]]:
    """``entrapment.json`` and the discovery rows from the parsed arguments (reads only)."""
    from giap.io.fasta import read_genome_fasta
    from giap.io.peptides import read_peptide_table
    from giap.io.protein_fasta import read_protein_fasta

    out = a.out.resolve()
    if out == REPO or REPO in out.parents:
        raise Precondition(f"--out {a.out} is inside the repository")
    files = {
        "pep_a1": a.pep_a1,
        "k12_fasta": a.k12_fasta,
        "contaminants": a.contaminants,
        "genome": a.genome,
        "strain_db": a.strain_db,
        "runs_manifest": a.runs_manifest,
        "diann_cmd": a.diann_cmd,
        "provenance": a.provenance,
    }
    files = {k: v for k, v in files.items() if v is not None}
    for p in [*files.values(), *([a.out_a1] if a.out_a1 else [])]:
        if not p.exists():
            raise Precondition(f"missing input: {p}")
    if a.provenance is not None:
        prov = json.loads(a.provenance.read_text(encoding="utf-8"))
        recorded = {f.get("sha256") for f in prov.get("inputs", {}).get("fasta", [])}
        for p in (a.k12_fasta, a.contaminants):
            if sha256(p) not in recorded:
                raise Precondition(f"{p}: sha256 is not the one {a.provenance} records")
    settings = DEFAULT_SETTINGS
    if a.diann_cmd is not None:
        tokens = [t.strip() for t in a.diann_cmd.read_text(encoding="utf-8").splitlines()]
        if bad := settings_mismatch(parse_diann_cmd([t for t in tokens if t]), settings):
            raise Precondition(f"{a.diann_cmd}: settings differ: {'; '.join(bad)}")
    rows, _ = read_peptide_table(a.pep_a1)
    wt = {r.get("run") for r in _table(a.runs_manifest)}
    if foreign := sorted({str(r) for r in rows["run"]} - {str(r) for r in wt}):
        raise Precondition(f"{a.pep_a1}: run(s) not in {a.runs_manifest}: {', '.join(foreign)}")
    k12 = [r.sequence for r in read_protein_fasta(a.k12_fasta)]
    if not k12:
        raise Precondition(f"{a.k12_fasta}: no protein")
    recs = read_genome_fasta(a.genome)
    if not any(r.circular_tag for r in recs):
        raise Precondition(f"{a.genome}: no topology tag on any contig")
    reps = {
        r.id: Replicon(r.id, len(r.sequence), "circular" if r.circular_tag else "linear")
        for r in recs
    }
    index = build_index(
        {r.id: r.sequence for r in recs},
        reps,
        [r.sequence for r in read_protein_fasta(a.strain_db)],
    )
    cont = [r.sequence for r in read_protein_fasta(a.contaminants)]
    doc = evaluate(k12, cont, index, rows, settings)
    discoveries = doc.pop("discoveries")
    checks: dict[str, Any] = {
        "digest": {
            "not_in_universe": doc["peptide"]["not_in_universe"],
            "agrees": doc["peptide"]["not_in_universe"] == 0,
        }
    }
    if a.diann_cmd is not None:
        checks["settings"] = {"agrees": True}
    if a.out_a1 is not None:
        checks["run"], warn = cross_check_run(discoveries, a.out_a1, doc["peptide"]["N_E"])
        doc["warnings"] += warn
    doc |= {
        "schema": SCHEMA,
        "spec": _spec_sha(),
        "cross_checks": checks,
        "caveats": list(CAVEATS),
        "inputs": {k: {"path": str(p.resolve()), "sha256": sha256(p)} for k, p in files.items()},
        "settings": {
            "enzyme": settings.enzyme,
            "missed": settings.missed,
            "min_len": settings.min_len,
            "max_len": settings.max_len,
            "charges": list(settings.charges),
            "mz": list(settings.mz),
            "var_mods": settings.var_mods,
            "tolerant": index.tolerant.model_dump(),
        },
    }
    return doc, discoveries


def write_outputs(out: Path, doc: Mapping[str, Any], discoveries: Sequence[Sequence[str]]) -> None:
    from giap.io.tsv import write_tsv

    out.mkdir(parents=True, exist_ok=True)
    text = json.dumps(doc, sort_keys=True, indent=2) + "\n"
    (out / "entrapment.json").write_text(text, encoding="utf-8")
    write_tsv(
        out / "entrapment_discoveries.tsv", ["peptide", "bucket", "present_rule"], discoveries
    )


def main(argv: Sequence[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n", 1)[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("compute")
    for name in (
        "pep-a1",
        "k12-fasta",
        "contaminants",
        "genome",
        "strain-db",
        "runs-manifest",
        "out",
    ):
        c.add_argument(f"--{name}", type=Path, required=True)
    for name in ("out-a1", "diann-cmd", "provenance"):
        c.add_argument(f"--{name}", type=Path, default=None)
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
    p = doc["precursor"]
    print(
        f"wrote {a.out}/entrapment.json; precursor {p['assessment']}; "
        f"{len(doc['warnings'])} warning(s)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
