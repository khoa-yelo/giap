"""Attack nulls for the annotated-stop decoy cells (M7-16a, ADR-0058 items 1-2; statistician review
change 6; GIPPA_SPEC.md §7.2.1, §8.2, §8.5).

Placement FDR counts decoy hits to estimate chance target hits. ADR-0058 lets the tolerant
matcher consume a ``*`` at an annotated terminal stop (a *cell*) at cost 1, and gives each decoy
type the same rule at the *mirrored* cells (``reverse``: n - 1 - p; ``shuffle``: p). These attacks
check that the decoy read-through chance mass tracks the target's under the null, on the cells
the registered ``sixframe_index`` step writes into the index (format 2; ``load_stop_cells``), and
that a deliberately broken cell set fails.

**Exact expected mass, no peptide sampling.** For a wrong peptide of length L with budget
m = k(L) (``mapping.tolerant.max_diff``) drawn from a position-specific law, the chance mass of a
search space is the sum over its L-windows of P(window within m substitutions). A window is
searchable when every position is a residue or a marked cell; a marked cell never matches (it
costs exactly one substitution, any query letter; spec rule 4). T_rt / D_rt is the mass of the
windows holding >= 1 cell (target / decoy), T_reg the mass of cell-free windows. Multiplicity is
counted (expected number of windows hit, the quantity the placement decoy count estimates).
Swaps and blocks are left out: they never read a cell (rule 4), so the read-through part is
exactly the substitution mass.

Two worlds:

- ``gippasim``: gippasim genomes of ``GENOME_NT`` (circular chromosome with its origin-crossing
  gene, one circular plasmid, one linear contig), one per replicate, the biological gene density.
  Laws: i.i.d. from the target six-frame composition, and tryptic (C-terminal K/R, no internal
  K/R).
  (a) reverse, i.i.d.: D_rt = T_rt per replicate and L to rtol 1e-9 (a theorem; also the regular
  mass); (b) shuffle, i.i.d.: mean D_rt/T_rt >= 1 - EQUIV - se_multiplier * SE per L (the EQUIV
  margin of test_attack_genome_decoys.py); (c) tryptic, both decoys: the net shift of the replicate
  means, (mean D_rt - mean T_rt) / mean T_reg >= -EQUIV per L. Structural: per frame the decoy has
  the target's cell count, and the number of consumable windows (searchable, >= 1 cell) is equal
  in all three spaces for every L (so the shuffle identity mirror cannot mark more consumable
  positions than the target).
- ``stressed``: the dense build of tests/calibration/test_readthrough_stressed.py
  (``build_genome``, ~1 cell per 37 coding positions) at ``STRESSED_NT``, the i.i.d. law of that
  arm (K, P, Y, F uniform). (a) and the structural check as above; (b') shuffle one-sided:
  mean D_rt/T_rt >= 1 - EQUIV - se_multiplier * SE; (c') net shift >= -EQUIV for both decoys.
  Recorded (not a criterion): the same build with every gene's start ATG replaced by a law codon,
  which locates the shuffle excess of the stressed arm (Slurm 1889791: shuffle D_rt ~ 3 x T_rt).

**Teeth** (a passing broken build fails the test): no decoy cells (fails (a), (b), and (c') in the
dense build), the target cells reused for the reverse decoy (fails (a) and the structural check),
and an off-by-one reverse mirror n - p (fails (a)).

(d) ``recoded`` cells (ADR-0058 item 6, M7-16c) have their own attack, d1-d4, in
``test_attack_readthrough_recoded.py``.

R = ``replicates_fast.placement`` (marker ``calibration``) or ``replicates_full.placement``
(``calibration_full``); the stressed world uses R // 5 (each replicate carries ~1,800 cells).
Seeds: ``_criteria.replicate_rng`` streams. Per-L ratios go to
``<log_dir>/placement_readthrough_attack.jsonl``. No tolerance, bound, alpha or R is new (R3).
"""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping
from pathlib import Path

import _criteria as C
import numpy as np
import pytest
from helpers.readthrough_index import peptides_table
from helpers.readthrough_toys import rc, write_annotation, write_genome

from gippa.mapping.sixframe import READTHROUGH, STOP, Frame, encode_target

SUITE = "placement"
GENOME_NT = 100_000
STRESSED_NT = 200_000
EQUIV = 0.01  # the equivalence margin of test_attack_genome_decoys.py
LENGTHS = (10, 11, 12, 16)
KINDS = ("reverse", "shuffle")
K_R = [int(encode_target(a)[0]) for a in "KR"]
STRESSED_LAW = [int(encode_target(a)[0]) for a in "KPYF"]

pytestmark = [pytest.mark.calsuite(SUITE)]

SUITES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]

Cells = tuple[np.ndarray, ...]
Law = Callable[[int], np.ndarray]  # L -> (L, 19) per-position residue probabilities


# ---------------------------------------------------------------- exact mass


def _budgets() -> dict[int, int]:
    from gippa.config import resolve_config
    from gippa.mapping.tolerant import budget, parse_max_diff

    tcfg = resolve_config().mapping.tolerant
    ranges = parse_max_diff(tcfg.max_diff)
    out = {L: budget(L, tcfg, ranges) for L in LENGTHS}
    assert all(k is not None and k >= 1 for k in out.values()), out
    return out  # type: ignore[return-value]


def _p_within(qs: list[np.ndarray], m: int) -> np.ndarray:
    """P(<= m mismatches) for independent positions with match probabilities ``qs``."""
    c = [np.ones_like(qs[0])] + [np.zeros_like(qs[0]) for _ in range(m)]
    for q in qs:
        c = [c[0] * q] + [c[j] * q + c[j - 1] * (1 - q) for j in range(1, m + 1)]
    return np.sum(c, axis=0)


def masses(
    frames: tuple[Frame, ...], cells: Cells, law: np.ndarray, m: int
) -> tuple[float, float, int]:
    """(regular mass, read-through mass, consumable windows) of one space at L = len(law).

    The cells are written as ``READTHROUGH`` (whatever residue they cover, as the search text
    would hold them); a cell never matches, so it is one forced substitution.
    """
    L = law.shape[0]
    Q = np.zeros((L, READTHROUGH + 1))
    Q[:, :STOP] = law
    reg = rt = 0.0
    n_rt = 0
    for fr, cl in zip(frames, cells, strict=True):
        if len(fr.aa) < L:
            continue
        aa = fr.aa.astype(np.int64)
        aa[np.asarray(cl, dtype=np.int64)] = READTHROUGH
        W = np.lib.stride_tricks.sliding_window_view(aa, L)
        mark = W == READTHROUGH
        ok = ((W < STOP) | mark).all(axis=1)
        nm = mark.sum(axis=1)[ok]
        Wv = W[ok]
        p = _p_within([Q[j, Wv[:, j]] for j in range(L)], m)
        reg += float(p[nm == 0].sum())
        rt += float(p[nm > 0].sum())
        n_rt += int(((nm > 0) & (nm <= m)).sum())
    return reg, rt, n_rt


def _iid(f: np.ndarray) -> Law:
    return lambda L: np.tile(f, (L, 1))


def _tryptic(f: np.ndarray) -> Law:
    inner = f.copy()
    inner[K_R] = 0.0
    inner /= inner.sum()
    last = np.zeros_like(f)
    last[K_R] = f[K_R]
    last /= last.sum()
    return lambda L: np.vstack([np.tile(inner, (L - 1, 1)), last[None, :]])


def _composition(frames: tuple[Frame, ...]) -> np.ndarray:
    aa = np.concatenate([f.aa for f in frames])
    return np.bincount(aa[aa < STOP], minlength=STOP) / int((aa < STOP).sum())


# ---------------------------------------------------------------- worlds (the real step)


def _index(genome: Mapping[str, str], gm, tmp: Path, seed: int):  # type: ignore[no-untyped-def]
    """``sixframe_index`` (registered step, both decoy types) -> (index, decoys, StopCells)."""
    from gippa.config import resolve_config
    from gippa.mapping.sixframe import load_index, load_stop_cells
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    tmp.mkdir(parents=True, exist_ok=True)
    ann = write_annotation(gm, tmp / "ann")
    ins = {
        "genome": write_genome(genome, tmp / "genome.fa"),
        "replicons": ann["replicons"],
        "genes": ann["genes"],
        "peptides": peptides_table(tmp / "peptides.parquet"),  # one database, one de novo query
    }
    out = tmp / "sixframe.npz"
    reg.run_step("sixframe_index", resolve_config(overrides=[f"seed={seed}"]), ins,
                 {"sixframe": out})  # fmt: skip
    index, decoys = load_index(out, table=11)
    assert set(decoys) == set(KINDS), decoys.keys()
    return index, decoys, load_stop_cells(out, table=11)


def gippasim_world(seed: int, tmp: Path):  # type: ignore[no-untyped-def]
    import gippasim
    from gippa.annotate.bakta import parse_bakta_json

    cfg = gippasim.SimConfig(genome_length=GENOME_NT, n_plasmids=1, plasmid_length=(8_000, 20_000),
                            linear_contig_length=10_000)  # fmt: skip
    res = gippasim.simulate(cfg, seed=seed, outdir=tmp / "sim")
    ann = parse_bakta_json(res.bakta_json)
    return _index(ann.sequences, ann.model, tmp / "idx", seed)


def stressed_world(seed: int, tmp: Path, no_start_met: bool = False):  # type: ignore[no-untyped-def]
    from test_readthrough_stressed import build_genome

    seq, _reps, gm, genes, _prot = build_genome(seed, length=STRESSED_NT)
    if no_start_met:  # every gene's ATG -> a law codon (AAA = K) on its strand
        s = list(seq)
        for _tag, a, b, strand in genes:
            at, codon = (a, "AAA") if strand == "+" else (b - 3, rc("AAA"))
            s[at : at + 3] = codon
        seq = "".join(s)
    return _index({"s": seq}, gm, tmp / "idx", seed)


# ---------------------------------------------------------------- broken cell sets (teeth)


def _none(cells: Cells) -> Cells:
    return tuple(np.zeros(0, dtype=np.int64) for _ in cells)


def broken_sets(index, cells):  # type: ignore[no-untyped-def]
    """name -> {kind: decoy cells} of deliberately wrong builds."""
    n = [len(f.aa) for f in index.frames]
    tgt = cells.target
    return {
        "no_decoy_cells": {k: _none(tgt) for k in KINDS},
        "target_cells_on_reverse": {"reverse": tgt, "shuffle": cells.decoys["shuffle"]},
        "reverse_off_by_one": {
            "reverse": tuple(np.sort(nf - c[c > 0]) for nf, c in zip(n, tgt, strict=True)),
            "shuffle": cells.decoys["shuffle"],
        },
    }


# ---------------------------------------------------------------- one replicate


def replicate(index, decoys, cells, laws: Mapping[str, Law], sets=None):  # type: ignore[no-untyped-def]
    """{(set, law, L, space): (reg, rt, n_rt)} for the built cells (set "built") and ``sets``."""
    budgets = _budgets()
    spaces: dict[str, Mapping[str, Cells]] = {"built": dict(cells.decoys), **(sets or {})}
    out: dict[tuple[str, str, int, str], tuple[float, float, int]] = {}
    for law_name, law in laws.items():
        for L, m in budgets.items():
            q = law(L)
            out[("built", law_name, L, "target")] = masses(index.frames, cells.target, q, m)
            for name, dc in spaces.items():
                for kind in KINDS:
                    out[(name, law_name, L, kind)] = masses(decoys[kind], dc[kind], q, m)
    return out


def _structural(index, decoys, cells, res) -> None:  # type: ignore[no-untyped-def]
    for kind in KINDS:
        for fid, (t, d) in enumerate(zip(cells.target, cells.decoys[kind], strict=True)):
            assert len(t) == len(d), (kind, fid)
            assert np.all(decoys[kind][fid].aa[d] == STOP), (kind, fid)
    for (name, law, L, space), (_reg, _rt, n_rt) in res.items():
        if name == "built" and space in KINDS:
            n_t = res[("built", law, L, "target")][2]
            assert n_rt == n_t, f"{space} L={L}: {n_rt} consumable windows vs target {n_t}"


# ---------------------------------------------------------------- criteria


def _se(x: np.ndarray) -> float:
    return float(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0


def crit_reverse_exact(reps, name: str, law: str) -> list[str]:  # type: ignore[no-untyped-def]
    """(a): reverse D_rt = T_rt per replicate and L (rtol 1e-9); failures listed."""
    bad = []
    for i, r in enumerate(reps):
        for L in LENGTHS:
            t = r[("built", law, L, "target")]
            d = r[(name, law, L, "reverse")]
            if not np.isclose(d[1], t[1], rtol=1e-9, atol=0.0) or d[2] != t[2]:
                bad.append(f"rep {i} L={L}: D_rt {d[1]:.6g} vs T_rt {t[1]:.6g}")
            if name == "built" and not np.isclose(d[0], t[0], rtol=1e-9, atol=0.0):
                bad.append(f"rep {i} L={L}: regular D {d[0]:.6g} vs T {t[0]:.6g}")
    return bad


def ratios(reps, name: str, law: str, kind: str, L: int) -> np.ndarray:  # type: ignore[no-untyped-def]
    return np.array([r[(name, law, L, kind)][1] / r[("built", law, L, "target")][1] for r in reps])


def crit_ratio(reps, name: str, law: str, kind: str, k: float) -> list[str]:  # type: ignore[no-untyped-def]
    """(b): mean D_rt/T_rt >= 1 - EQUIV - k SE, per L (one-sided)."""
    bad = []
    for L in LENGTHS:
        x = ratios(reps, name, law, kind, L)
        if not x.mean() >= 1.0 - EQUIV - k * _se(x):
            bad.append(f"{kind} L={L}: mean D_rt/T_rt {x.mean():.4f} (SE {_se(x):.4f})")
    return bad


def net_shift(reps, name: str, law: str, kind: str, L: int) -> float:  # type: ignore[no-untyped-def]
    t = np.mean([r[("built", law, L, "target")][1] for r in reps])
    d = np.mean([r[(name, law, L, kind)][1] for r in reps])
    reg = np.mean([r[("built", law, L, "target")][0] for r in reps])
    return float((d - t) / reg)


def crit_net(reps, name: str, law: str, kinds=KINDS) -> list[str]:  # type: ignore[no-untyped-def]
    """(c): (mean D_rt - mean T_rt) / mean T_reg >= -EQUIV, per L and decoy type."""
    return [
        f"{kind} L={L}: net shift {s:+.2e}"
        for kind in kinds
        for L in LENGTHS
        if not (s := net_shift(reps, name, law, kind, L)) >= -EQUIV
    ]


def _log(record: dict[str, object]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "placement_readthrough_attack.jsonl").open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


def _summary(reps, names, laws, full: bool, world: str) -> dict[str, object]:  # type: ignore[no-untyped-def]
    rec: dict[str, object] = {"suite": f"placement.readthrough_attack.{world}", "full": full,
                              "R": len(reps)}  # fmt: skip
    for name in names:
        for law in laws:
            for kind in KINDS:
                for L in LENGTHS:
                    x = ratios(reps, name, law, kind, L)
                    rec[f"{name}.{law}.{kind}.L{L}"] = {
                        "mean_ratio": float(x.mean()), "se": _se(x),
                        "min": float(x.min()), "max": float(x.max()),
                        "net_shift": net_shift(reps, name, law, kind, L),
                    }  # fmt: skip
    for law in laws:
        for L in LENGTHS:
            share = [r[("built", law, L, "target")][1] / r[("built", law, L, "target")][0]
                     for r in reps]  # fmt: skip
            rec[f"rt_share.{law}.L{L}"] = float(np.mean(share))
    return rec


# ---------------------------------------------------------------- tests


def _stream(world: str, i: int) -> int:
    return int(C.replicate_rng(f"placement.readthrough_attack.{world}", i).integers(0, 2**31 - 1))


@pytest.mark.parametrize("full", SUITES)
def test_attack_readthrough_decoys_gippasim(full: bool, tmp_path: Path) -> None:
    R = C.replicates(SUITE, full)
    k = C.load_tolerances().se_multiplier
    reps = []
    n_cells = []
    for i in range(R):
        index, decoys, cells = gippasim_world(_stream("gippasim", i), tmp_path / f"g{i}")
        f = _composition(index.frames)
        res = replicate(index, decoys, cells, {"iid": _iid(f), "tryptic": _tryptic(f)},
                        broken_sets(index, cells))  # fmt: skip
        _structural(index, decoys, cells, res)
        n_cells.append(sum(len(c) for c in cells.target))
        reps.append(res)
    rec = _summary(reps, ["built", "no_decoy_cells"], ["iid", "tryptic"], full, "gippasim")
    rec["mean_cells"] = float(np.mean(n_cells))
    _log(rec)
    assert min(n_cells) >= 50, f"too few cells per genome: {min(n_cells)}"
    # (a) reverse exact, (b) shuffle within EQUIV + k SE, (c) tryptic net shift
    assert not (bad := crit_reverse_exact(reps, "built", "iid")), bad
    assert not (bad := crit_ratio(reps, "built", "iid", "shuffle", k)), bad
    assert not (bad := crit_net(reps, "built", "tryptic")), bad
    # teeth: every broken set fails (a) for reverse; no decoy cells also fails (b)
    for name in ("no_decoy_cells", "target_cells_on_reverse", "reverse_off_by_one"):
        assert crit_reverse_exact(reps, name, "iid"), f"broken set {name} passed (a)"
    assert crit_ratio(reps, "no_decoy_cells", "iid", "shuffle", k), "no decoy cells passed (b)"


@pytest.mark.parametrize("full", SUITES)
def test_attack_readthrough_decoys_stressed(full: bool, tmp_path: Path) -> None:
    R = max(2, C.replicates(SUITE, full) // 5)
    k = C.load_tolerances().se_multiplier
    f = np.zeros(STOP)
    f[STRESSED_LAW] = 0.25
    laws = {"iid": _iid(f)}
    reps, reps_nomet = [], []
    for i in range(R):
        seed = _stream("stressed", i)
        index, decoys, cells = stressed_world(seed, tmp_path / f"s{i}")
        res = replicate(index, decoys, cells, laws, {"no_decoy_cells": {
            kind: _none(cells.target) for kind in KINDS}})  # fmt: skip
        _structural(index, decoys, cells, res)
        reps.append(res)
        index, decoys, cells = stressed_world(seed, tmp_path / f"m{i}", no_start_met=True)
        reps_nomet.append(replicate(index, decoys, cells, laws))
    rec = _summary(reps, ["built", "no_decoy_cells"], ["iid"], full, "stressed")
    rec["no_start_met"] = _summary(reps_nomet, ["built"], ["iid"], full, "stressed_no_start_met")
    _log(rec)
    assert not (bad := crit_reverse_exact(reps, "built", "iid")), bad
    assert not (bad := crit_ratio(reps, "built", "iid", "shuffle", k)), bad
    assert not (bad := crit_net(reps, "built", "iid")), bad
    # teeth at the dense build: no decoy cells fails (a), (b) and the net-shift criterion (c')
    assert crit_reverse_exact(reps, "no_decoy_cells", "iid"), "no decoy cells passed (a)"
    assert crit_ratio(reps, "no_decoy_cells", "iid", "shuffle", k), "no decoy cells passed (b)"
    assert crit_net(reps, "no_decoy_cells", "iid"), "no decoy cells passed (c')"
