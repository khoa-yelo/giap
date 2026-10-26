"""Attack (d): annotated recoded Sec/Pyl cells in the placement decoys (M7-16c, ADR-0058 item 6;
spec "16c amendment" > "Re-runs (16c)" item 2; statistician pre-code check Q2; GIAP_SPEC.md §7.2.1,
§8.2, §8.5).

16c adds a cell at every annotated internal ``U``-over-TGA / ``O``-over-TAG codon (site kind
``recoded``), mirrored into each decoy like a terminal cell. At biological density these are 3 of
~4,700 cells and ~1e-7 of the chance mass, so the pooled 16a criteria
(``test_attack_readthrough_decoys.py``) cannot see them: a build that drops recoded cells from
the decoys only would pass (a) to rtol 1e-9. This file therefore splits every mass by kind and
judges the recoded part on its own. Masses are the exact expected chance masses of the 16a attack
(``masses`` there; same laws, budgets, LENGTHS, EQUIV, se_multiplier): a window is *recoded* when
it holds >= 1 recoded cell, *terminal* when it holds only terminal cells, *regular* when it holds
none.

- **d1** recoded-dense giapsim world: the 16a ``giapsim`` world (GENOME_NT, circular chromosome,
  one circular plasmid, one linear contig, biological terminal density, ~120 terminal sites) with
  ``n_selenoproteins = N_SEC`` (60 of ~106-126 eligible genes), so every replicate holds >= 50
  recoded cells. On the recoded part: (a) reverse, i.i.d.: D_rec = T_rec per replicate and L to
  rtol 1e-9 with equal consumable windows; (b) shuffle, i.i.d.: mean D_rec/T_rec >= 1 - EQUIV -
  se_multiplier * SE per L; (c) tryptic, both decoys: (mean D_rec - mean T_rec) / mean T_reg >=
  -EQUIV per L. Structural, per kind (terminal, recoded): per frame the decoy has the target's
  cell count, every decoy cell holds ``*``, and the consumable windows are equal in all three
  spaces for every L. The recoded sites of the index are exactly ``truth/selenoproteins.tsv``.
  The decoy recoded cells are derived independently of ``mirror_cells`` (n - 1 - p / p) and must
  be a subset of the built decoy cells.
- **d2** teeth (same replicates): ``decoys_terminal_only`` (decoys mirrored from the terminal
  cells only; the target keeps its recoded cells) fails (a) and the structural check;
  ``recoded_off_by_one`` (reverse recoded cells at n - p) fails (a); ``target_without_recoded``
  (target drops its recoded cells, decoys keep theirs) fails the structural check.
- **d3** (``calibration_full`` only) the Nissle and K-12 recordings through ``sixframe_index``:
  3 recoded sites, structural equality per kind, reverse i.i.d. exactness of both parts; shuffle
  and tryptic ratios and net shifts logged, only net shift >= -EQUIV asserted.
- **d4** the sites ``load_stop_cells`` returns equal the ``site_keys`` the variants step passes
  to M7-01 ``align_hits`` (``run_variants_qc`` on the same annotation files the index was built
  from), on d4 giapsim worlds and (full) the two recordings.

R = ``replicates_fast.placement`` / ``replicates_full.placement``; d4 uses R // 10 worlds. Seeds:
``_criteria.replicate_rng`` streams ``placement.readthrough_attack.recoded*``. Per-L ratios go to
``<log_dir>/placement_readthrough_attack.jsonl``. No tolerance, bound, alpha or R is new (R3).
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from helpers.readthrough_index import peptides_table
from helpers.readthrough_toys import write_annotation, write_genome
from test_attack_readthrough_decoys import (
    EQUIV,
    GENOME_NT,
    KINDS,
    LENGTHS,
    SUITE,
    SUITES,
    Law,
    _budgets,
    _composition,
    _iid,
    _p_within,
    _se,
    _tryptic,
)

from giap.mapping.sixframe import READTHROUGH, STOP, Frame

N_SEC = 60
MIN_REC_CELLS = 50
REPO = Path(__file__).resolve().parents[3]
RECORDINGS = {
    "nissle": (REPO / "tests/data/genomes/nissle1917.fna.gz",
               REPO / "tests/data/recorded/bakta/1.12.2-full/nissle/nissle1917.json.gz"),
    "k12": (REPO / "tests/data/genomes/k12_mg1655.fna.gz",
            REPO / "tests/data/recorded/bakta/1.12.2-full/k12/k12_mg1655.json.gz"),
}  # fmt: skip
PARTS = ("terminal", "recoded")

pytestmark = [pytest.mark.calsuite(SUITE)]

Cells = tuple[np.ndarray, ...]
Split = tuple[Cells, Cells]  # (terminal cells, recoded cells) per frame id
_E = np.zeros(0, dtype=np.int64)


# ---------------------------------------------------------------- kind-split exact mass


def split_masses(
    frames: tuple[Frame, ...], split: Split, law: np.ndarray, m: int
) -> tuple[float, float, float, int, int]:
    """(regular, terminal-only, recoded) mass and consumable windows (terminal-only, recoded) of
    one space at L = len(law). Same mass as the 16a ``masses``; a cell never matches."""
    L = law.shape[0]
    Q = np.zeros((L, READTHROUGH + 1))
    Q[:, :STOP] = law
    reg = rt_t = rt_r = 0.0
    n_t = n_r = 0
    for fr, term, rec in zip(frames, *split, strict=True):
        if len(fr.aa) < L:
            continue
        aa = fr.aa.astype(np.int64)
        cells = np.concatenate([np.asarray(term, np.int64), np.asarray(rec, np.int64)])
        aa[cells] = READTHROUGH
        is_rec = np.zeros(len(aa), dtype=bool)
        is_rec[np.asarray(rec, np.int64)] = True
        W = np.lib.stride_tricks.sliding_window_view(aa, L)
        mark = W == READTHROUGH
        ok = ((W < STOP) | mark).all(axis=1)
        nm = mark.sum(axis=1)[ok]
        nr = np.lib.stride_tricks.sliding_window_view(is_rec, L).sum(axis=1)[ok]
        Wv = W[ok]
        p = _p_within([Q[j, Wv[:, j]] for j in range(L)], m)
        reg += float(p[nm == 0].sum())
        rt_t += float(p[(nm > 0) & (nr == 0)].sum())
        rt_r += float(p[nr > 0].sum())
        cons = (nm > 0) & (nm <= m)
        n_t += int((cons & (nr == 0)).sum())
        n_r += int((cons & (nr > 0)).sum())
    return reg, rt_t, rt_r, n_t, n_r


# ---------------------------------------------------------------- worlds (the real step)


def build_index(genome: Mapping[str, str], gm: Any, tmp: Path, seed: int) -> dict[str, Any]:
    """``sixframe_index`` (registered step, both decoy types) on ``genome``/``gm``; the loaded
    index, decoys, StopCells and the input files (for d4)."""
    from giap.config import resolve_config
    from giap.mapping.sixframe import load_index, load_stop_cells
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

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
    reg.run_step("sixframe_index", resolve_config(overrides=[f"seed={seed}"]), ins,
                 {"sixframe": out})  # fmt: skip
    index, decoys = load_index(out)
    assert set(decoys) == set(KINDS), decoys.keys()
    return {"index": index, "decoys": decoys, "cells": load_stop_cells(out), "inputs": ins}


def recoded_world(seed: int, tmp: Path) -> tuple[dict[str, Any], set[tuple[str, str, int]]]:
    import giapsim
    from giap.annotate.bakta import parse_bakta_json

    cfg = giapsim.SimConfig(genome_length=GENOME_NT, n_plasmids=1, plasmid_length=(8_000, 20_000),
                            linear_contig_length=10_000, n_selenoproteins=N_SEC)  # fmt: skip
    res = giapsim.simulate(cfg, seed=seed, outdir=tmp / "sim")
    ann = parse_bakta_json(res.bakta_json)
    truth = pd.read_csv(Path(res.outdir) / "truth" / "selenoproteins.tsv", sep="\t")
    keys = {(str(r.contig), str(r.strand), int(r.codon_start)) for r in truth.itertuples()}
    return build_index(ann.sequences, ann.model, tmp / "idx", seed), keys


def recording_world(name: str, tmp: Path) -> dict[str, Any]:
    from giap.annotate.bakta import parse_bakta_json
    from giap.io.fasta import read_genome_fasta

    fna, js = RECORDINGS[name]
    genome = {r.id: r.sequence for r in read_genome_fasta(fna)}
    gm = parse_bakta_json(js, genome=genome).model
    return build_index(genome, gm, tmp / name, 1)


# ---------------------------------------------------------------- the kind split


def mirror(kind: str, frames: tuple[Frame, ...], cells: Cells) -> Cells:
    """Independent of ``mirror_cells``: n - 1 - p (reverse), p (shuffle)."""
    if kind == "reverse":
        return tuple(np.sort(len(f.aa) - 1 - np.asarray(c, np.int64)) for f, c in
                     zip(frames, cells, strict=True))  # fmt: skip
    return tuple(np.asarray(c, np.int64).copy() for c in cells)


def built_split(w: Mapping[str, Any]) -> tuple[dict[str, Split], list[str]]:
    """{space: (terminal, recoded)} of the built index, and the containment failures."""
    from giap.mapping.readthrough import frame_cells

    index, decoys, cells = w["index"], w["decoys"], w["cells"]
    rec_t = frame_cells(index, [s for s in cells.sites if s.kind == "recoded"])
    bad: list[str] = []
    out: dict[str, Split] = {}
    for space, built, rec in [("target", cells.target, rec_t)] + [
        (k, cells.decoys[k], mirror(k, decoys[k], rec_t)) for k in KINDS
    ]:
        term = []
        for fid, (b, r) in enumerate(zip(built, rec, strict=True)):
            inside = np.isin(r, b)
            if not inside.all():
                bad.append(f"{space} frame {fid}: {int((~inside).sum())} recoded cells not built")
            term.append(np.setdiff1d(b, r).astype(np.int64))
        out[space] = (tuple(term), tuple(np.asarray(r, np.int64) for r in rec))
    return out, bad


def broken_splits(w: Mapping[str, Any], built: Mapping[str, Split]) -> dict[str, dict[str, Split]]:
    """d2: name -> {space: split} of deliberately wrong builds."""
    decoys = w["decoys"]
    t_term, t_rec = built["target"]
    none = tuple(_E for _ in t_rec)
    return {
        "decoys_terminal_only": {"target": built["target"],
                                 **{k: (mirror(k, decoys[k], t_term), none) for k in KINDS}},
        "recoded_off_by_one": {"target": built["target"], "shuffle": built["shuffle"],
                               "reverse": (built["reverse"][0], tuple(
                                   np.sort(len(f.aa) - c[c > 0]) for f, c in
                                   zip(decoys["reverse"], t_rec, strict=True)))},
        "target_without_recoded": {"target": (t_term, none),
                                   "reverse": built["reverse"], "shuffle": built["shuffle"]},
    }  # fmt: skip


# ---------------------------------------------------------------- one replicate


Res = dict[tuple[str, str, int, str], tuple[float, float, float, int, int]]


def replicate(w: Mapping[str, Any], sets: Mapping[str, Mapping[str, Split]],
              laws: Mapping[str, Law], broken_laws: tuple[str, ...] = ("iid",)) -> Res:  # fmt: skip
    """{(set, law, L, space): split_masses} for ``sets`` ("built" with every law, the others with
    ``broken_laws``)."""
    frames = {"target": w["index"].frames, **{k: w["decoys"][k] for k in KINDS}}
    out: Res = {}
    for law_name, law in laws.items():
        for L, m in _budgets().items():
            q = law(L)
            for name, sp in sets.items():
                if name != "built" and law_name not in broken_laws:
                    continue
                for space, split in sp.items():
                    out[(name, law_name, L, space)] = split_masses(frames[space], split, q, m)
    return out


def structural(w: Mapping[str, Any], sp: Mapping[str, Split], res: Res, name: str) -> list[str]:
    """Per kind: equal per-frame cell counts, decoy cells on ``*``, equal consumable windows."""
    bad = []
    tgt = sp["target"]
    for kind in KINDS:
        for j, part in enumerate(PARTS):
            for fid, (t, d) in enumerate(zip(tgt[j], sp[kind][j], strict=True)):
                if len(t) != len(d):
                    bad.append(f"{kind} {part} frame {fid}: {len(d)} cells vs target {len(t)}")
                aa = w["decoys"][kind][fid].aa
                if len(d) and (d.min() < 0 or d.max() >= len(aa) or np.any(aa[d] != STOP)):
                    bad.append(f"{kind} {part} frame {fid}: a cell not on '*'")
    for (nm, law, L, space), v in res.items():
        if nm == name and space in KINDS:
            t = res[(nm, law, L, "target")]
            for j, part in enumerate(PARTS):
                if v[3 + j] != t[3 + j]:
                    bad.append(f"{space} {part} L={L}: {v[3 + j]} consumable windows vs {t[3 + j]}")
    return bad


# ---------------------------------------------------------------- criteria (recoded part)

_REC, _TERM = 2, 1  # index of the part's mass in split_masses


def crit_exact(reps: list[Res], name: str, law: str, part: int = _REC) -> list[str]:
    """(a): reverse D = T for ``part`` per replicate and L (rtol 1e-9), equal windows."""
    bad = []
    for i, r in enumerate(reps):
        for L in LENGTHS:
            t, d = r[(name, law, L, "target")], r[(name, law, L, "reverse")]
            if not np.isclose(d[part], t[part], rtol=1e-9, atol=0.0) or d[part + 2] != t[part + 2]:
                bad.append(f"rep {i} L={L}: D {d[part]:.6g} vs T {t[part]:.6g}")
    return bad


def ratios(reps: list[Res], name: str, law: str, kind: str, L: int, part: int = _REC) -> np.ndarray:
    return np.array([r[(name, law, L, kind)][part] / r[(name, law, L, "target")][part]
                     for r in reps])  # fmt: skip


def crit_ratio(reps: list[Res], law: str, kind: str, k: float) -> list[str]:
    """(b): mean D_rec/T_rec >= 1 - EQUIV - k SE, per L (one-sided)."""
    bad = []
    for L in LENGTHS:
        x = ratios(reps, "built", law, kind, L)
        if not x.mean() >= 1.0 - EQUIV - k * _se(x):
            bad.append(f"{kind} L={L}: mean D_rec/T_rec {x.mean():.4f} (SE {_se(x):.4f})")
    return bad


def net_shift(reps: list[Res], law: str, kind: str, L: int, part: int = _REC) -> float:
    t = np.mean([r[("built", law, L, "target")][part] for r in reps])
    d = np.mean([r[("built", law, L, kind)][part] for r in reps])
    reg = np.mean([r[("built", law, L, "target")][0] for r in reps])
    return float((d - t) / reg)


def crit_net(reps: list[Res], law: str) -> list[str]:
    """(c): (mean D_rec - mean T_rec) / mean T_reg >= -EQUIV, per L and decoy type."""
    return [f"{kind} L={L}: net shift {s:+.2e}" for kind in KINDS for L in LENGTHS
            if not (s := net_shift(reps, law, kind, L)) >= -EQUIV]  # fmt: skip


def _summary(reps: list[Res], laws: tuple[str, ...], full: bool, world: str) -> dict[str, Any]:
    rec: dict[str, Any] = {"suite": f"placement.readthrough_attack.{world}", "full": full,
                           "R": len(reps)}  # fmt: skip
    for law in laws:
        for kind in KINDS:
            for L in LENGTHS:
                for pname, part in (("recoded", _REC), ("terminal", _TERM)):
                    x = ratios(reps, "built", law, kind, L, part)
                    rec[f"{pname}.{law}.{kind}.L{L}"] = {
                        "mean_ratio": float(x.mean()), "se": _se(x),
                        "min": float(x.min()), "max": float(x.max()),
                        "net_shift": net_shift(reps, law, kind, L, part),
                    }  # fmt: skip
        for L in LENGTHS:
            rec[f"rec_share.{law}.L{L}"] = float(np.mean(
                [r[("built", law, L, "target")][_REC] / r[("built", law, L, "target")][0]
                 for r in reps]))  # fmt: skip
    return rec


def _log(record: dict[str, Any]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "placement_readthrough_attack.jsonl").open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


def _stream(world: str, i: int) -> int:
    return int(C.replicate_rng(f"placement.readthrough_attack.{world}", i).integers(0, 2**31 - 1))


# ---------------------------------------------------------------- d4 helper


class _Seen(Exception):
    pass


def variants_site_keys(w: Mapping[str, Any]) -> frozenset[tuple[str, str, int]]:
    """The ``stop_sites`` the variants step's chain hands to ``align_hits``, read the way the
    step reads its inputs (``read_gene_model`` + ``read_genome_fasta`` on the same files)."""
    import giap.steps.variants as V
    from giap.config import resolve_config
    from giap.io.fasta import read_genome_fasta
    from giap.steps.mapping import read_gene_model

    ins = w["inputs"]
    genome = {r.id: r.sequence for r in read_genome_fasta(ins["genome"])}
    gm = read_gene_model(ins["genes"], ins["replicons"])
    seen: list[frozenset[tuple[str, str, int]]] = []

    def spy(*a: Any, **k: Any) -> Any:
        seen.append(frozenset(k["stop_sites"]))
        raise _Seen

    ph = pd.DataFrame({c: [] for c in ("peptide", "source", "placement_q", "uniqueness")})
    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(V, "align_hits", spy)
        with pytest.raises(_Seen):
            V.run_variants_qc(ph, None, pd.DataFrame({"peptide": [], "search_db": []}), gm,
                              genome, resolve_config())  # fmt: skip
    return seen[0]


def index_site_keys(w: Mapping[str, Any]) -> frozenset[tuple[str, str, int]]:
    return frozenset((s.contig, s.strand, s.codon_start) for s in w["cells"].sites)


# ---------------------------------------------------------------- tests


@pytest.mark.parametrize("full", SUITES)
def test_attack_d1_d2_recoded_giapsim(full: bool, tmp_path: Path) -> None:
    R = C.replicates(SUITE, full)
    k = C.load_tolerances().se_multiplier
    reps: list[Res] = []
    n_rec, n_term, truth_bad, struct_bad = [], [], [], []
    teeth_struct: dict[str, int] = {"decoys_terminal_only": 0, "target_without_recoded": 0}
    for i in range(R):
        w, truth = recoded_world(_stream("recoded", i), tmp_path / f"g{i}")
        got = {(s.contig, s.strand, s.codon_start) for s in w["cells"].sites if s.kind == "recoded"}
        if got != truth or len(truth) != N_SEC:
            truth_bad.append(f"rep {i}: {len(got)} recoded sites vs {len(truth)} truth")
        built, contain = built_split(w)
        sets = {"built": built, **broken_splits(w, built)}
        f = _composition(w["index"].frames)
        res = replicate(w, sets, {"iid": _iid(f), "tryptic": _tryptic(f)})
        struct_bad += [f"rep {i}: {b}" for b in contain + structural(w, built, res, "built")]
        for name in teeth_struct:
            teeth_struct[name] += bool(structural(w, sets[name], res, name))
        n_rec.append(sum(len(c) for c in built["target"][1]))
        n_term.append(sum(len(c) for c in built["target"][0]))
        reps.append(res)
        if not full and i == 0:
            print(f"rep 0: {n_rec[0]} recoded / {n_term[0]} terminal target cells")
    rec = _summary(reps, ("iid", "tryptic"), full, "recoded_giapsim")
    rec |= {"N_SEC": N_SEC, "mean_rec_cells": float(np.mean(n_rec)), "min_rec_cells": min(n_rec),
            "mean_term_cells": float(np.mean(n_term)),
            "teeth_struct_fail_reps": teeth_struct}  # fmt: skip
    for name in ("decoys_terminal_only", "recoded_off_by_one"):
        rec[f"teeth.{name}.a_fail_reps"] = len({b.split(" ")[1] for b in
                                                crit_exact(reps, name, "iid")})  # fmt: skip
    _log(rec)
    # d1 preconditions: the world is recoded-dense and the index holds exactly the truth codons
    assert min(n_rec) >= MIN_REC_CELLS, f"too few recoded cells per replicate: {min(n_rec)}"
    assert not truth_bad, truth_bad[:5]
    # d1 structural per kind, (a) reverse exact on both parts, (b) shuffle, (c) tryptic net shift
    assert not struct_bad, struct_bad[:5]
    assert not (bad := crit_exact(reps, "built", "iid")), bad[:5]
    assert not (bad := crit_exact(reps, "built", "iid", _TERM)), bad[:5]
    assert not (bad := crit_ratio(reps, "iid", "shuffle", k)), bad
    assert not (bad := crit_net(reps, "tryptic")), bad
    # d2 teeth
    assert crit_exact(reps, "decoys_terminal_only", "iid"), "decoys_terminal_only passed (a)"
    assert crit_exact(reps, "recoded_off_by_one", "iid"), "recoded_off_by_one passed (a)"
    assert teeth_struct["decoys_terminal_only"] == R, teeth_struct
    assert teeth_struct["target_without_recoded"] == R, teeth_struct


@pytest.mark.parametrize("full", SUITES)
def test_attack_d4_index_and_variants_site_sets(full: bool, tmp_path: Path) -> None:
    n = max(2, C.replicates(SUITE, full) // 10)
    for i in range(n):
        w, truth = recoded_world(_stream("recoded_d4", i), tmp_path / f"v{i}")
        keys = index_site_keys(w)
        assert truth <= keys
        assert variants_site_keys(w) == keys, f"world {i}: variants and index site sets differ"


@pytest.mark.calibration_full
@pytest.mark.parametrize("name", sorted(RECORDINGS))
def test_attack_d3_recordings(name: str, tmp_path: Path) -> None:
    w = recording_world(name, tmp_path)
    rec_sites = [s for s in w["cells"].sites if s.kind == "recoded"]
    assert len(rec_sites) == 3 and {s.codon for s in rec_sites} == {"TGA"}, rec_sites
    built, contain = built_split(w)
    f = _composition(w["index"].frames)
    res = replicate(w, {"built": built}, {"iid": _iid(f), "tryptic": _tryptic(f)})
    n_cells = {sp: [sum(len(c) for c in part) for part in split] for sp, split in built.items()}
    out: dict[str, Any] = {"suite": f"placement.readthrough_attack.recoded_{name}", "full": True,
                           "cells": n_cells}  # fmt: skip
    for (_n, law, L, space), v in res.items():
        t = res[("built", law, L, "target")]
        out[f"{law}.{space}.L{L}"] = {
            "reg": v[0], "rt_term": v[1], "rt_rec": v[2], "windows": [v[3], v[4]],
            "rec_ratio": v[2] / t[2] if t[2] else None, "rec_net_shift": (v[2] - t[2]) / t[0],
        }  # fmt: skip
    _log(out)
    assert not contain, contain
    assert not (bad := structural(w, built, res, "built")), bad
    assert n_cells["target"][1] == 3 and all(v == n_cells["target"] for v in n_cells.values())
    assert not (bad := crit_exact([res], "built", "iid")), bad
    assert not (bad := crit_exact([res], "built", "iid", _TERM)), bad
    bad = [f"{law} {kind} L={L}: {s:+.2e}" for law in ("iid", "tryptic") for kind in KINDS
           for L in LENGTHS if not (s := net_shift([res], law, kind, L)) >= -EQUIV]  # fmt: skip
    assert not bad, bad
    # d4 on the recording
    assert variants_site_keys(w) == index_site_keys(w)
