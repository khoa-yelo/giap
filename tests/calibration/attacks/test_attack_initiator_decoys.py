"""Attack nulls for the initiator-M decoy cells (M11-20 § Statistics S1, § Tests "Attack", required
before merge; GIPPA_SPEC.md §7.2.5 "decoys are searched identically", §8.5). Modelled on
test_attack_readthrough_decoys.py.

A query whose first residue is M is credited at an annotated non-ATG start *cell*: target and
``shuffle`` decoy at cell p, window [p, p + k), the query's residues 1..k-1 against frame
p+1..p+k-1; ``reverse`` decoy at cell c = n - 1 - p with the **mirrored** window [c - k + 1, c + 1)
(cell last), residues 1..k-1 against decoy c-k+1..c-1 (rule 3). The cells are the ones the
registered ``sixframe_index`` step writes (index format 4, ``load_start_cells``).

**Exact expected mass, no peptide sampling.** For an M-initial wrong peptide of length k from a
position-specific law (position 0 = M), the credited chance mass of a space is the sum over its
credited windows of P(residues 1..k-1 match); multiplicity counted (the expected number of
credited windows hit). ``plain`` is the target's M-initial plain mass (windows starting at an M).

Worlds: ``gippasim`` (GENOME_NT circular chromosome with its origin-crossing gene, a plasmid and a
linear contig; starts ATG 0.5 / GTG 0.35 / TTG 0.15) and ``stressed`` (every gene GTG/TTG, short
sORFs, genes at both linear contig ends, a GTG across the origin as example 7b, a minus gene
across the origin of a second circular contig). Checks per world:

- (a) ``reverse``, i.i.d. law: credited decoy mass = target mass per replicate and k, rtol 1e-9;
- (b) ``shuffle``, i.i.d.: mean D/T >= 1 - EQUIV - se_multiplier * SE per k (one-sided);
- (c) tryptic law (inner residues without K/R, last K/R): (mean D - mean T) / mean plain >=
  -EQUIV - se_multiplier * SE per k and decoy type (SE of the per-replicate (D - T) / mean plain;
  one-sided, the form of (b): few credited windows are tryptic-compatible, so the replicate
  shifts are wide; at R = 20 on a correct build the plain -EQUIV bound already fails at k = 12);
- (d) structural: per frame the decoy cell count equals the target's, reverse cells hold the
  target's start residue, the usable credited windows are equal for every k in 7..40, and the
  implementation (``match_exact(starts=...)``) returns every model window: target rows for
  ``M`` + frame[p+1:p+k], reverse rows for ``M`` + decoy[c-k+1:c] (mirrored onto [p, p+k)),
  shuffle rows for ``M`` + decoy[p+1:p+k];
- (e) recorded only: the N-terminome law (M + the position-specific profile of the world's own
  annotated N-termini), net shift (D - T) / plain per decoy type.

**Teeth** (each fails): no decoy cells (a, b, d), the target cells reused as reverse cells (a),
the off-by-one mirror n - p (a), and the forward-reading reverse window [c, c + k) (a, d).

R = ``replicates_fast.placement`` (``calibration``) or ``replicates_full.placement``
(``calibration_full``); the stressed world R // 5. Seeds: ``_criteria.replicate_rng`` streams. Per-k
records go to ``<log_dir>/placement_initiator_attack.jsonl``. No tolerance, bound or R is new (R3).
"""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
from helpers.initiator_toys import cds, rest_of, usable, windows
from helpers.readthrough_index import peptides_table
from helpers.readthrough_toys import model, rc, write_annotation, write_genome

from gippa.mapping.decoys import mirror_aa
from gippa.mapping.sixframe import STOP, Frame, decode, encode_target, translate_text

SUITE = "placement"
GENOME_NT = 100_000
EQUIV = 0.01  # the equivalence margin of test_attack_genome_decoys.py
LENGTHS = (7, 10, 12)  # the spec's k (S2, (e))
KS = range(7, 41)
KINDS = ("reverse", "shuffle")
M = int(encode_target("M")[0])
K_R = [int(encode_target(a)[0]) for a in "KR"]
STARTS = (("ATG", 0.5), ("GTG", 0.35), ("TTG", 0.15))

pytestmark = [pytest.mark.calsuite(SUITE)]
SUITES = [
    pytest.param(False, marks=pytest.mark.calibration, id="fast"),
    pytest.param(True, marks=pytest.mark.calibration_full, id="full"),
]
Cells = tuple[np.ndarray, ...]
Law = Callable[[int], np.ndarray]  # k -> (k, 19); row 0 is the initiator (M)


# ---------------------------------------------------------------- exact mass


def credited(frames: tuple[Frame, ...], cells: Cells, law: np.ndarray, rule: str) -> float:
    """Credited mass of one space; ``rule``: target / shuffle / reverse / forward (a tooth)."""
    k, out = law.shape[0], 0.0
    for fr, cl in zip(frames, cells, strict=True):
        rr = "shuffle" if rule == "forward" else rule  # forward: [c, c + k) on the reverse frame
        for _, a, b in windows(rr, len(fr.aa), np.asarray(cl), k):
            res = rest_of(rr, fr.aa, a, b).astype(np.int64)
            if np.all(res < STOP):
                out += float(law[0, M] * np.prod(law[np.arange(1, k), res]))
    return out


def plain(frames: tuple[Frame, ...], law: np.ndarray) -> float:
    k, out = law.shape[0], 0.0
    for fr in frames:
        if len(fr.aa) < k:
            continue
        W = np.lib.stride_tricks.sliding_window_view(fr.aa.astype(np.int64), k)
        W = W[(W[:, 0] == M) & (W < STOP).all(axis=1)]
        out += float(law[0, M] * np.prod(law[np.arange(1, k), W[:, 1:]], axis=1).sum())
    return out


def _with_m(rest: Callable[[int], np.ndarray]) -> Law:
    def law(k: int) -> np.ndarray:
        first = np.zeros((1, STOP))
        first[0, M] = 1.0
        return np.vstack([first, rest(k)])

    return law


def iid(f: np.ndarray) -> Law:
    return _with_m(lambda k: np.tile(f, (k - 1, 1)))


def tryptic(f: np.ndarray) -> Law:
    inner, last = f.copy(), np.zeros_like(f)
    inner[K_R] = 0.0
    inner /= inner.sum()
    last[K_R] = f[K_R]
    last /= last.sum()
    return _with_m(lambda k: np.vstack([np.tile(inner, (k - 2, 1)), last[None, :]]))


def nterminome(gm: Any) -> Law:
    """Position-specific profile of the annotated proteins' residues 1..k-1 (I -> L)."""
    aas = [encode_target(g.aa) for g in gm.genes if g.aa and g.aa[0] == "M"]

    def rest(k: int) -> np.ndarray:
        out = np.empty((k - 1, STOP))
        for j in range(1, k):
            col = np.array([a[j] for a in aas if len(a) > j and a[j] < STOP], dtype=np.int64)
            out[j - 1] = np.bincount(col, minlength=STOP)[:STOP] + 0.5  # Jeffreys: no zero row
            out[j - 1] /= out[j - 1].sum()
        return out

    return _with_m(rest)


def composition(frames: tuple[Frame, ...]) -> np.ndarray:
    aa = np.concatenate([f.aa for f in frames])
    return np.bincount(aa[aa < STOP], minlength=STOP) / int((aa < STOP).sum())


# ---------------------------------------------------------------- worlds (the real step)


def _index(genome: Mapping[str, str], gm: Any, tmp: Path, seed: int) -> tuple[Any, Any, Any]:
    from gippa.config import resolve_config
    from gippa.mapping.sixframe import load_index, load_start_cells
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    tmp.mkdir(parents=True, exist_ok=True)
    ann = write_annotation(gm, tmp / "ann")
    ins = {
        "genome": write_genome(genome, tmp / "genome.fa"),
        "replicons": ann["replicons"],
        "genes": ann["genes"],
        "peptides": peptides_table(tmp / "peptides.parquet"),
    }
    out = tmp / "sixframe.npz"
    reg.run_step(
        "sixframe_index", resolve_config(overrides=[f"seed={seed}"]), ins, {"sixframe": out}
    )
    index, decoys = load_index(out, table=11)
    assert set(decoys) == set(KINDS), decoys.keys()
    return index, decoys, load_start_cells(out, table=11)


def gippasim_world(seed: int, tmp: Path) -> tuple[Any, Any, Any, Any]:
    import gippasim
    from gippa.annotate.bakta import parse_bakta_json

    cfg = gippasim.SimConfig(
        genome_length=GENOME_NT,
        n_plasmids=1,
        plasmid_length=(8_000, 20_000),
        linear_contig_length=10_000,
        origin_crossing_gene=True,
        start_codons=STARTS,
    )
    res = gippasim.simulate(cfg, seed=seed, outdir=tmp / "sim")
    ann = parse_bakta_json(res.bakta_json)
    return (*_index(ann.sequences, ann.model, tmp / "idx", seed), ann.model)


_SENSE = [
    a + b + c
    for a in "ACGT"
    for b in "ACGT"
    for c in "ACGT"
    if a + b + c not in ("TAA", "TAG", "TGA")
]


def stressed_world(seed: int, tmp: Path) -> tuple[Any, Any, Any, Any]:
    """Circular 30 kb and 5 kb + linear 6 kb, every gene GTG/TTG, ~1/3 short sORFs (8-30 aa)."""
    from gippa.model import Replicon

    rng = np.random.default_rng(seed)
    sizes = {"s": 30_000, "p": 5_000, "l": 6_000}
    seqs = {c: list("".join(rng.choice(list("ACGT"), size=n))) for c, n in sizes.items()}
    genes: list[Any] = []

    def plant(contig: str, start: int, strand: str, n_aa: int) -> int:
        seq, n = seqs[contig], len(seqs[contig])
        first = str(rng.choice(["GTG", "TTG"]))
        body = "".join(rng.choice(_SENSE, size=n_aa - 1))
        gs = first + body + "TAA"
        for j, b in enumerate(gs if strand == "+" else rc(gs)):
            seq[(start + j) % n] = b
        kind = "sorf" if n_aa <= 30 else "cds"
        genes.append(
            cds(
                f"{contig}{len(genes)}",
                contig,
                start,
                start + len(gs),
                strand,
                "M" + translate_text(body, 11),
                first,
                kind,
            )
        )
        return start + len(gs)

    plant("s", 30_000 - 1, "+", 20)  # example 7b: the GTG at bases L-1, 0, 1
    plant("p", 5_000 - 30, "-", 25)  # a minus gene across the origin, first codon [45, 48)
    plant("l", 0, "+", 30)  # first codon [0, 3)
    plant("l", 6_000 - 93, "-", 30)  # first codon [L - 3, L)
    for contig, lo, hi in (
        ("s", 120, 30_000 - 200),
        ("p", 120, 5_000 - 200),
        ("l", 120, 6_000 - 200),
    ):
        pos = lo
        while True:
            n_aa = int(rng.integers(8, 31)) if rng.random() < 1 / 3 else int(rng.integers(31, 200))
            if pos + 3 * (n_aa + 1) > hi:
                break
            pos = plant(contig, pos, str(rng.choice(["+", "-"])), n_aa) + int(rng.integers(3, 40))
    reps = {c: Replicon(c, n, "linear" if c == "l" else "circular") for c, n in sizes.items()}
    genome = {c: "".join(s) for c, s in seqs.items()}
    gm = model(reps, genes)
    return (*_index(genome, gm, tmp / "idx", seed), gm)


# ---------------------------------------------------------------- teeth


def broken_sets(index: Any, cells: Any) -> dict[str, tuple[dict[str, Cells], str]]:
    """name -> ({kind: decoy cells}, reverse rule) of deliberately wrong builds."""
    n = [len(f.aa) for f in index.frames]
    tgt, sh, rv = cells.target, cells.decoys["shuffle"], cells.decoys["reverse"]
    none = tuple(np.zeros(0, dtype=np.int64) for _ in tgt)
    return {
        "no_decoy_cells": ({"reverse": none, "shuffle": none}, "reverse"),
        "target_cells_on_reverse": ({"reverse": tgt, "shuffle": sh}, "reverse"),
        "reverse_off_by_one": (
            {
                "reverse": tuple(np.sort(nf - c[c > 0]) for nf, c in zip(n, tgt, strict=True)),
                "shuffle": sh,
            },
            "reverse",
        ),
        "forward_reverse_window": ({"reverse": rv, "shuffle": sh}, "forward"),
    }


# ---------------------------------------------------------------- one replicate


def replicate(index: Any, decoys: Any, cells: Any, laws: Mapping[str, Law]) -> dict[tuple, float]:
    """{(set, law, k, space): mass} of the built cells and every broken set ("plain": plain)."""
    spaces = {
        "built": ({k: cells.decoys[k] for k in KINDS}, "reverse"),
        **broken_sets(index, cells),
    }
    out: dict[tuple, float] = {}
    for law_name, law in laws.items():
        for k in LENGTHS:
            q = law(k)
            out[("built", law_name, k, "target")] = credited(
                index.frames, cells.target, q, "target"
            )
            out[("built", law_name, k, "plain")] = plain(index.frames, q)
            for name, (dc, rrule) in spaces.items():
                out[(name, law_name, k, "reverse")] = credited(
                    decoys["reverse"], dc["reverse"], q, rrule
                )
                out[(name, law_name, k, "shuffle")] = credited(
                    decoys["shuffle"], dc["shuffle"], q, "shuffle"
                )
    return out


def structural(
    index: Any, decoys: Any, dc: Mapping[str, Cells], tgt: Cells, rrule: str = "reverse"
) -> list[str]:
    """(d) without the implementation check: counts, reverse residues, usable windows per k."""
    bad = []
    for fid, f in enumerate(index.frames):
        rep_len = index.replicons[f.contig].length
        for kind in KINDS:
            d = np.asarray(dc[kind][fid], dtype=np.int64)
            if len(d) != len(tgt[fid]):
                bad.append(f"{kind} frame {fid}: {len(d)} cells vs {len(tgt[fid])}")
                continue
            if (
                kind == "reverse"
                and len(d)
                and not np.array_equal(np.sort(decoys[kind][fid].aa[d]), np.sort(f.aa[tgt[fid]]))
            ):
                bad.append(f"reverse frame {fid}: cells do not hold the start residues")
            rule = ("shuffle" if rrule == "forward" else "reverse") if kind == "reverse" else kind
            for k in KS:
                u = usable(rule, decoys[kind][fid], d, k, rep_len)
                t = usable("target", f, tgt[fid], k, rep_len)
                if u != t:
                    bad.append(f"{kind} frame {fid} k={k}: {u} usable vs target {t}")
    return bad


def implementation_agrees(index: Any, decoys: Any, cells: Any) -> list[str]:
    """(d): ``match_exact(starts=...)`` returns every model credited window (k in LENGTHS + 40)."""
    from gippa.mapping.exact import match_exact

    want: dict[str, set[tuple]] = {"reverse": set(), "shuffle": set()}
    for kind in KINDS:
        for fid, f in enumerate(index.frames):
            n = len(f.aa)
            if index.replicons[f.contig].length < 3 * 40:
                continue
            for k in (*LENGTHS, 40):
                for _, a, b in windows("target", n, cells.target[fid], k):
                    r = rest_of("target", f.aa, a, b)
                    if np.all(r < STOP):
                        iv = index.to_interval(fid, a, b)
                        want[kind].add(("M" + decode(r), "target", iv.start, iv.end, iv.strand))
                for _, a, b in windows(kind, n, cells.decoys[kind][fid], k):
                    r = rest_of(kind, decoys[kind][fid].aa, a, b)
                    if np.all(r < STOP):
                        iv = index.to_interval(fid, *mirror_aa(kind, n, a, b))
                        want[kind].add(("M" + decode(r), "decoy", iv.start, iv.end, iv.strand))
    bad = []
    for kind, rows in want.items():
        qs = sorted({r[0] for r in rows})
        df = match_exact(qs, index, decoys, dict.fromkeys(qs, kind), starts=cells)
        got = set(
            zip(df["peptide"], df["search"], df["start"], df["end"], df["strand"], strict=True)
        )
        miss = rows - got
        if miss:
            bad.append(
                f"{kind}: {len(miss)} of {len(rows)} model windows missing, e.g. {sorted(miss)[:3]}"
            )
    return bad


# ---------------------------------------------------------------- criteria


def _se(x: np.ndarray) -> float:
    x = x[~np.isnan(x)]
    return float(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else 0.0


def crit_reverse_exact(reps: list[dict], name: str, law: str = "iid") -> list[str]:
    bad = []
    for i, r in enumerate(reps):
        for k in LENGTHS:
            t, d = r[("built", law, k, "target")], r[(name, law, k, "reverse")]
            if not np.isclose(d, t, rtol=1e-9, atol=0.0):
                bad.append(f"rep {i} k={k}: reverse D {d:.6g} vs T {t:.6g}")
    return bad


def ratios(reps: list[dict], name: str, law: str, kind: str, k: int) -> np.ndarray:
    """D/T per replicate; NaN where T = 0 (e.g. no K/R-free credited 16-mer: recorded laws)."""
    t = np.array([r[("built", law, k, "target")] for r in reps])
    d = np.array([r[(name, law, k, kind)] for r in reps])
    return np.divide(d, t, out=np.full(len(t), np.nan), where=t > 0)


def crit_ratio(reps: list[dict], name: str, kind: str, mult: float, law: str = "iid") -> list[str]:
    bad = []
    for k in LENGTHS:
        x = ratios(reps, name, law, kind, k)
        if not x.mean() >= 1.0 - EQUIV - mult * _se(x):
            bad.append(f"{kind} k={k}: mean D/T {x.mean():.4f} (SE {_se(x):.4f})")
    return bad


def net_shift(reps: list[dict], name: str, law: str, kind: str, k: int) -> float:
    t = np.mean([r[("built", law, k, "target")] for r in reps])
    d = np.mean([r[(name, law, k, kind)] for r in reps])
    return float((d - t) / np.mean([r[("built", law, k, "plain")] for r in reps]))


def crit_net(reps: list[dict], name: str, mult: float, law: str = "tryptic") -> list[str]:
    """(c): net shift >= -EQUIV - mult * SE, SE of the per-replicate (D - T) / mean plain."""
    bad = []
    for kind in KINDS:
        for k in LENGTHS:
            s = net_shift(reps, name, law, kind, k)
            pl = np.mean([r[("built", law, k, "plain")] for r in reps])
            se = _se(
                np.array([r[(name, law, k, kind)] - r[("built", law, k, "target")] for r in reps])
                / pl
            )
            if not s >= -EQUIV - mult * se:
                bad.append(f"{kind} k={k}: net shift {s:+.2e} (SE {se:.2e})")
    return bad


def _log(rec: dict[str, object]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "placement_initiator_attack.jsonl").open("a") as fh:
        fh.write(json.dumps(rec, default=float) + "\n")


def _summary(reps: list[dict], full: bool, world: str, extra: dict[str, object]) -> dict:
    rec: dict[str, object] = {
        "suite": f"placement.initiator_attack.{world}",
        "full": full,
        "R": len(reps),
        **extra,
    }
    for name in ("built", "no_decoy_cells"):
        for law in ("iid", "tryptic", "nterminome"):
            for kind in KINDS:
                for k in LENGTHS:
                    x = ratios(reps, name, law, kind, k)
                    rec[f"{name}.{law}.{kind}.k{k}"] = {
                        "mean_ratio": float(np.nanmean(x)) if np.isfinite(x).any() else None,
                        "se": _se(x),
                        "net_shift": net_shift(reps, name, law, kind, k),
                    }
    for law in ("iid", "nterminome"):
        for k in LENGTHS:
            rec[f"credited_share.{law}.k{k}"] = float(
                np.mean(
                    [r[("built", law, k, "target")] / r[("built", law, k, "plain")] for r in reps]
                )
            )
    return rec


# ---------------------------------------------------------------- tests


def _stream(world: str, i: int) -> int:
    return int(C.replicate_rng(f"placement.initiator_attack.{world}", i).integers(0, 2**31 - 1))


def run_attack(world: str, R: int, full: bool, tmp: Path) -> None:
    build = gippasim_world if world == "gippasim" else stressed_world
    mult = C.load_tolerances().se_multiplier
    reps, n_sites, struct, teeth_struct = [], [], [], {}
    for i in range(R):
        index, decoys, cells, gm = build(_stream(world, i), tmp / f"{world}{i}")
        f = composition(index.frames)
        reps.append(
            replicate(
                index,
                decoys,
                cells,
                {"iid": iid(f), "tryptic": tryptic(f), "nterminome": nterminome(gm)},
            )
        )
        n_sites.append(len(cells.sites))
        struct += structural(index, decoys, cells.decoys, cells.target)
        struct += implementation_agrees(index, decoys, cells)
        for name, (dc, rrule) in broken_sets(index, cells).items():
            teeth_struct[name] = teeth_struct.get(name, 0) + len(
                structural(index, decoys, dc, cells.target, rrule)
            )
    _log(
        _summary(
            reps,
            full,
            world,
            {"mean_sites": float(np.mean(n_sites)), "teeth_structural_failures": teeth_struct},
        )
    )
    assert min(n_sites) >= 20, f"too few start sites per genome: {min(n_sites)}"
    assert not struct, struct[:10]  # (d)
    assert not (bad := crit_reverse_exact(reps, "built")), bad  # (a)
    assert not (bad := crit_ratio(reps, "built", "shuffle", mult)), bad  # (b)
    assert not (bad := crit_net(reps, "built", mult)), bad  # (c)
    for name in (
        "no_decoy_cells",
        "target_cells_on_reverse",
        "reverse_off_by_one",
        "forward_reverse_window",
    ):
        assert crit_reverse_exact(reps, name), f"broken build {name} passed (a)"
    assert crit_ratio(reps, "no_decoy_cells", "shuffle", mult), "no decoy cells passed (b)"
    assert teeth_struct["no_decoy_cells"] > 0, "no decoy cells passed (d)"
    assert teeth_struct["forward_reverse_window"] > 0, "forward reverse window passed (d)"


@pytest.mark.parametrize("full", SUITES)
def test_attack_initiator_decoys_gippasim(full: bool, tmp_path: Path) -> None:
    run_attack("gippasim", C.replicates(SUITE, full), full, tmp_path)


@pytest.mark.parametrize("full", SUITES)
def test_attack_initiator_decoys_stressed(full: bool, tmp_path: Path) -> None:
    run_attack("stressed", max(2, C.replicates(SUITE, full) // 5), full, tmp_path)
