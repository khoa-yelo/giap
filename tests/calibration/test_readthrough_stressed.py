"""Stressed read-through arm (M7-16a; statistician review change 5; ADR-0058 Evidence; spec test
plan "Stressed arm"). Suite ``placement``: ``replicates_fast.placement`` /
``replicates_full.placement``, existing keys ``levels`` and ``se_multiplier`` only (R3).

A sampled arm in which the read-through mass is not vacuous, through the real chain:
``sixframe_index`` (with ``genes``) -> ``map_exact`` -> ``map_tolerant`` -> ``placement``
(``best_costs`` + placement q-values), via ``gippa.steps.registry.run_step``.

**Build** (one linear contig of ``LENGTH`` nt per replicate):

- dense short annotated genes: runs of 5-15 same-strand CDS, each ``ATG`` + 28-38 law codons +
  a stop (TAA/TAG/TGA), 0-2 law codons between genes of a run (so the 3' flank of a stop is the
  next gene's coding frame), runs on alternating strands;
- a reduced-alphabet coding law: ``K``, ``P``, ``Y``, ``F`` equally likely (s = sum f^2 = 0.25;
  none of them is in a BLOCKS / NEAR_BLOCKS pair);
- ``N_WRONG`` wrong 10-12-mers drawn i.i.d. from the law (no origin; 75 % 12-mers, which carry
  the read-through rows: measured locally ~0.2 read-through rows per wrong 12-mer at 1 Mb), and
  ``N_TRUE`` correct peptides: 13-16-aa windows of the genes with 1 or 2 planted substitutions
  (strata 1 and 2), origin = the window. Their lower E puts them at the head of each m >= 1
  list, so every list is accepted past the correct block into the wrong peptides until
  (D + 1) / T reaches q: the FDP criterion binds where the read-through mass sits. Database
  input (``reverse``) or de novo input (``shuffle``).

The density, in chance-mass terms: about one marked cell per 37 coding positions (G/N_coding
~ 2.7e-2), so the expected read-through share of the m = 2 chance mass of a law peptide is
~ 2G/(N_coding (1 - s)) ~ 7 % (Nissle: ~ 0.08 %). Over all six frames G/N ~ 4.5e-3, ~ 10x the
Appendix B density (the statistician's "50x" holds for the coding-frame mass the law peptides
see, not for G/N over six frames; reported to the lead).

**Per replicate** (wrong peptides only, rows at the peptide's best cost in each search, no
frameshift rows): T_rt = target rows whose window holds a target cell; D_rt = decoy rows whose
window holds a mirrored cell; T = all target rows; FDP per level q = false / accepted among the
accepted peptides (placement_q <= q) with m_target >= 1 (strata 1 and 2 pooled); false = wrong,
or correct with no best-cost target locus overlapping its origin.

**Criteria, per arm:**

0. realised mean T_rt >= 100 per replicate (``T_RT_MIN``), asserted first: otherwise the arm is
   mis-built (a failure, never "vacuous");
1. FDP in the m >= 1 strata <= q + se_multiplier * SE at every q in ``levels`` (harness
   ``calibrate_levels``, ADR-0015 confirmation);
2. mean(D_rt) >= mean(T_rt) - se_multiplier * SE(D_rt - T_rt) (one-sided, ADR-0019 form);
3. **negative control**: the same replicates with the decoy cells removed (decoy searches
   without marks, ``match_tolerant`` with no target queries) must **fail** criterion 2; a
   passing control fails the test.

Frameshift splits are off in this arm (``mapping.tolerant.frameshift.enabled=false``): a split
side never consumes a cell (rule 7), and the frameshift decoy marks are covered by
tests/unit/test_frameshift_readthrough.py and property (e). Diagnostics per arm (mean T_rt,
D_rt, T, read-through share T_rt / T, FDP and accepted per level, per-replicate accepted target
rows that consumed a cell, seconds per replicate) go to
``<log_dir>/placement_readthrough_stressed.jsonl`` for the statistician.
"""

from __future__ import annotations

import json
import shutil
import time
from dataclasses import dataclass, field
from pathlib import Path

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _gate import calibrate_levels
from helpers.readthrough_counts import best_cost_rows, holds_cell
from helpers.readthrough_toys import gene, model, rc, write_annotation, write_genome

pytestmark = [
    pytest.mark.calibration,
    pytest.mark.calsuite("placement"),
    pytest.mark.calibration_unshardable,  # direct checks over every replicate (shard 0, whole)
]

SUITE = "placement"
LENGTH = 1_000_000
N_WRONG = 1_000
WRONG_LEN = ((10, 0.125), (11, 0.125), (12, 0.75))  # 12-mers carry the read-through rows
N_TRUE = 2_000
TRUE_LEN = (13, 16)
T_RT_MIN = 100.0
LAW = "KPYF"
LAW_CODONS = {"K": ("AAA", "AAG"), "P": ("CCT", "CCC", "CCA", "CCG"), "Y": ("TAT", "TAC"),
              "F": ("TTT", "TTC")}  # fmt: skip
STOPS = ("TAA", "TAG", "TGA")
OVERRIDES = ["mapping.tolerant.frameshift.enabled=false"]
LEVELS = C.load_tolerances().levels


@dataclass
class Rep:
    t_rt: float
    d_rt: float
    t_all: float
    d_rt_control: float
    fdp: np.ndarray
    accepted: np.ndarray
    accepted_rt: np.ndarray  # accepted target rows (best cost) that consumed a cell, per level
    secs: float
    extra: dict[str, float] = field(default_factory=dict)


# ---------------------------------------------------------------- the synthetic genome


def _law(rng: np.random.Generator, n: int) -> str:
    return "".join(rng.choice(list(LAW), size=n))


def _codons(aa: str, rng: np.random.Generator) -> str:
    return "".join(LAW_CODONS[a][int(rng.integers(len(LAW_CODONS[a])))] for a in aa)


def build_genome(seed: int, length: int = LENGTH):  # type: ignore[no-untyped-def]
    """(sequence, genes [(tag, start, end, strand)], protein per gene) of one replicate."""
    from gippa.model import Replicon

    rng = np.random.default_rng([seed, 1616])
    parts: list[str] = []
    genes: list[tuple[str, int, int, str]] = []
    prot: dict[str, str] = {}
    pos = 0
    strand = "+"
    k = 0
    while pos < length - 2_000:
        run: list[tuple[str, str]] = []  # (gene-strand nt, protein) incl. trailing spacer
        for _ in range(int(rng.integers(5, 16))):
            body = _law(rng, int(rng.integers(28, 39)))
            nt = "ATG" + _codons(body, rng) + STOPS[int(rng.integers(3))]
            spacer = _codons(_law(rng, int(rng.integers(0, 3))), rng)
            run.append((nt + spacer, "M" + body))
        block = "".join(x for x, _ in run)
        if strand == "-":
            block = rc(block)
        off = 0
        lens = [len(x) for x, _ in run]
        for (_, p), ln in zip(run, lens, strict=True):
            glen = 3 * (len(p) + 1)
            s = pos + off if strand == "+" else pos + len(block) - off - glen
            tag = f"R{k:05d}"
            k += 1
            genes.append((tag, s, s + glen, strand))
            prot[tag] = p
            off += ln
        parts.append(block)
        pos += len(block)
        strand = "-" if strand == "+" else "+"
    tail = length - pos
    parts.append(_codons(_law(rng, tail // 3), rng) + "A" * (tail % 3))
    seq = "".join(parts)
    assert len(seq) == length
    reps = {"s": Replicon("s", length, "linear")}
    gm = model(reps, [gene(t, "s", a, b, st) for t, a, b, st in genes])
    return seq, reps, gm, genes, prot


def _peptides(seed: int, genes, prot, source: str):  # type: ignore[no-untyped-def]
    """(peptide table, wrong set, {peptide: origin (start, end, strand)})."""
    rng = np.random.default_rng([seed, 2016])
    wrong: set[str] = set()
    lens = [n for n, _ in WRONG_LEN]
    probs = [p for _, p in WRONG_LEN]
    while len(wrong) < N_WRONG:
        wrong.add(_law(rng, int(rng.choice(lens, p=probs))))
    origin: dict[str, tuple[int, int, str]] = {}
    while len(origin) < N_TRUE:
        tag, s, e, strand = genes[int(rng.integers(len(genes)))]
        p = prot[tag]
        L = int(rng.integers(TRUE_LEN[0], TRUE_LEN[1] + 1))
        a = int(rng.integers(1, len(p) - L + 1))  # skip the start M
        w = list(p[a : a + L])
        n_sub = int(rng.integers(1, 3))
        for j in rng.choice(L, size=n_sub, replace=False):
            w[int(j)] = str(rng.choice([x for x in LAW if x != w[int(j)]]))
        pep = "".join(w)
        if pep in wrong or pep in origin:
            continue
        if strand == "+":
            iv = (s + 3 * a, s + 3 * (a + L), "+")
        else:
            iv = (e - 3 * (a + L), e - 3 * a, "-")
        origin[pep] = iv
    peps = sorted(wrong) + sorted(origin)
    n = len(peps)
    df = pd.DataFrame(
        {
            "peptide": peps,
            "source": [source] * n,
            "q_value": [0.001] * n if source == "database" else [np.nan] * n,
            "pep": [0.01] * n if source == "database" else [np.nan] * n,
            "score": [np.nan] * n if source == "database" else list(rng.uniform(0.5, 1.0, n)),
            "is_decoy": [False] * n,
            "search_db": ["other"] * n,
        }
    )
    return df, wrong, origin


# ---------------------------------------------------------------- one replicate


def run_replicate(seed: int, kind: str, tmp: Path) -> Rep:
    from gippa.config import resolve_config
    from gippa.mapping.sixframe import load_index, load_stop_cells
    from gippa.mapping.tolerant import match_tolerant
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg
    from gippa.steps.mapping import decoy_types, select_queries

    t0 = time.perf_counter()
    load_step_modules()
    tmp.mkdir(parents=True, exist_ok=True)
    seq, _reps, gm, genes, prot = build_genome(seed)
    source = "database" if kind == "reverse" else "denovo"
    table, wrong, origin = _peptides(seed, genes, prot, source)
    ann = write_annotation(gm, tmp / "ann")
    files = {
        "genome": write_genome({"s": seq}, tmp / "genome.fa"),
        "peptides": tmp / "peptides.parquet",
        **ann,
    }
    table.to_parquet(files["peptides"], index=False)
    cfg = resolve_config(overrides=[f"seed={seed}", *OVERRIDES])
    for name in ("sixframe_index", "map_exact", "map_tolerant", "placement"):
        spec = reg.get(name)
        ins = {k: files[k] for k in (*spec.inputs, *spec.optional_inputs) if k in files}
        outs = {k: tmp / f"{k}.out" for k in spec.outputs}
        reg.run_step(name, cfg, ins, outs)
        files.update(outs)
    cells = load_stop_cells(files["sixframe"])
    hits = pd.concat(
        [pd.read_parquet(files["hits_exact"]), pd.read_parquet(files["hits_tolerant"])],
        ignore_index=True,
    )
    best = best_cost_rows(hits)
    best = best.assign(rt=holds_cell(best, cells))
    is_wrong = best["peptide"].isin(wrong)
    tgt = best[(best["search"] == "target") & is_wrong]
    dec = best[(best["search"] == "decoy") & is_wrong & (best["decoy_type"] == kind)]
    # negative control: the same decoy search with the decoy cells removed
    index, decoys = load_index(files["sixframe"])
    import gippa.mapping.readthrough as rt

    none = tuple(np.zeros(0, dtype=np.int64) for _ in index.frames)
    broken = rt.StopCells(
        sites=cells.sites, target=cells.target, decoys=dict.fromkeys(decoys, none)
    )
    q = select_queries(pd.read_parquet(files["peptides"]), cfg)
    peps = sorted(set(q["peptide"]) & wrong)
    ctl = match_tolerant(peps, index, decoys, decoy_types(q, cfg), set(), cfg.mapping.tolerant,
                         cells=broken)  # fmt: skip
    ctl = best_cost_rows(ctl)
    d_ctl = float(holds_cell(ctl, broken).sum()) if len(ctl) else 0.0
    # FDP over the accepted m >= 1 peptides at every level
    per = pd.read_parquet(files["placement"])
    per = per[per["decoy_type"] == kind]
    tbest = best[best["search"] == "target"]
    loci: dict[str, list[tuple[int, int, str, bool]]] = {}
    for r in tbest.itertuples(index=False):
        loci.setdefault(str(r.peptide), []).append((int(r.start), int(r.end), str(r.strand),
                                                    bool(r.rt)))  # fmt: skip

    def is_true(p: str) -> bool:
        o = origin.get(p)
        if o is None:
            return False
        return any(st == o[2] and s < o[1] and o[0] < e for s, e, st, _ in loci.get(p, ()))

    fdp, acc, acc_rt = [], [], []
    m1 = per[per["m_target"].fillna(-1) >= 1]
    for level in LEVELS:
        a = m1[m1["placement_q"].fillna(np.inf) <= level]["peptide"].astype(str).tolist()
        n_false = sum(not is_true(p) for p in a)
        fdp.append(n_false / len(a) if a else 0.0)
        acc.append(float(len(a)))
        acc_rt.append(float(sum(x[3] for p in a for x in loci.get(p, ()))))
    secs = time.perf_counter() - t0
    shutil.rmtree(tmp, ignore_errors=True)
    return Rep(
        t_rt=float(tgt["rt"].sum()), d_rt=float(dec["rt"].sum()), t_all=float(len(tgt)),
        d_rt_control=d_ctl, fdp=np.array(fdp), accepted=np.array(acc),
        accepted_rt=np.array(acc_rt), secs=secs,
        extra={"n_cells": float(sum(len(c) for c in cells.target)), "n_sites": len(cells.sites)},
    )  # fmt: skip


# ---------------------------------------------------------------- the arm


def balance_fails(t: np.ndarray, d: np.ndarray, k: float) -> tuple[bool, float]:
    """(criterion 2 fails, SE): mean(D) < mean(T) - k SE(D - T)."""
    diff = d - t
    se = float(diff.std(ddof=1) / np.sqrt(len(diff))) if len(diff) > 1 else 0.0
    return bool(d.mean() < t.mean() - k * se), se


def _diagnostics(record: dict[str, object]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "placement_readthrough_stressed.jsonl").open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


def run_arm(kind: str, full: bool, tmp: Path) -> dict[int, Rep]:
    stream = f"placement.readthrough_stressed.{kind}"
    R = C.replicates(SUITE, full)
    k = C.load_tolerances().se_multiplier
    seeds = {C.suite_seed(stream, i): i for i in range(2 * R)}
    memo: dict[int, Rep] = {}

    def rep(i: int) -> Rep:
        if i not in memo:
            seed = int(C.replicate_rng(stream, i).integers(0, 2**31 - 1))
            memo[i] = run_replicate(seed, kind, tmp / f"r{i}")
        return memo[i]

    for i in range(R):  # block 1 first: the build check precedes every criterion
        rep(i)
    t_rt = np.array([memo[i].t_rt for i in range(R)])
    assert t_rt.mean() >= T_RT_MIN, (
        f"{kind} arm mis-built: realised mean T_rt {t_rt.mean():.1f} < {T_RT_MIN} per replicate"
    )
    gate = calibrate_levels(
        SUITE, "fdr", LEVELS, lambda g: rep(seeds[int(g.bit_generator.seed_seq.entropy)]).fdp, R,
        stream=stream,
    )  # fmt: skip
    reps = [memo[i] for i in sorted(memo)]
    t = np.array([r.t_rt for r in reps])
    d = np.array([r.d_rt for r in reps])
    dc = np.array([r.d_rt_control for r in reps])
    fails, se = balance_fails(t, d, k)
    ctl_fails, ctl_se = balance_fails(t, dc, k)
    _diagnostics({
        "suite": stream, "R_run": len(reps), "full": full, "mean_T_rt": float(t.mean()),
        "mean_D_rt": float(d.mean()), "se_D_minus_T": se,
        "mean_T": float(np.mean([r.t_all for r in reps])),
        "rt_share": float(t.sum() / max(1.0, sum(r.t_all for r in reps))),
        "mean_D_rt_control": float(dc.mean()), "control_se": ctl_se,
        "levels": list(LEVELS), "mean_fdp": np.mean([r.fdp for r in reps], axis=0).tolist(),
        "mean_accepted": np.mean([r.accepted for r in reps], axis=0).tolist(),
        "accepted_rt_rows_per_replicate": [r.accepted_rt.tolist() for r in reps],
        "secs_per_replicate": float(np.mean([r.secs for r in reps])),
        "mean_cells": float(np.mean([r.extra["n_cells"] for r in reps])),
        "gate": gate.detail(),
    })  # fmt: skip
    gate.assert_ok()
    assert not fails, (
        f"{kind}: mean D_rt {d.mean():.2f} < mean T_rt {t.mean():.2f} - {k} SE ({se:.3f})"
    )
    assert ctl_fails, (
        f"{kind} negative control passed: without decoy cells mean D_rt {dc.mean():.2f} vs "
        f"mean T_rt {t.mean():.2f} - {k} SE ({ctl_se:.3f}); the balance check has no teeth"
    )
    return memo


@pytest.mark.parametrize("kind", ["reverse", "shuffle"])
def test_readthrough_stressed_arm(tmp_path: Path, kind: str) -> None:
    run_arm(kind, False, tmp_path)


@pytest.mark.calibration_full
@pytest.mark.parametrize("kind", ["reverse", "shuffle"])
def test_readthrough_stressed_arm_full(tmp_path: Path, kind: str) -> None:
    run_arm(kind, True, tmp_path)


# ---------------------------------------------------------------- the build itself (no matcher)


def test_build_is_dense_and_annotated_stops_are_stops() -> None:
    """Unmarked sanity check of the generator (runs with the calibration marker of the module):
    every annotated gene ends in a stop on its strand, genes are 30-40 codons, and the coding
    frame of a gene's 3' flank inside a run is the next gene."""
    seq, _reps, _gm, genes, prot = build_genome(5, length=60_000)
    assert len(genes) >= 400
    for tag, s, e, strand in genes:
        codon = seq[e - 3 : e] if strand == "+" else rc(seq[s : s + 3])
        assert codon in STOPS, tag
        assert 30 <= (e - s) // 3 <= 40 and (e - s) % 3 == 0
        assert len(prot[tag]) == (e - s) // 3 - 1
