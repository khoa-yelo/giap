"""Calibration arms of M11-20 (docs/specs/M11-20-initiator-methionine.md § Tests "Calibration",
S3-S6, AC 10, AC 12): fast tier, R = ``replicates_fast.<suite>`` at 100 kb; no full-tier arm
(statistician Q-S1). Existing tolerances only (R3).

- **I1, placement FDR level with credit** (suite ``placement``): gippasim 100 kb with start codons
  ATG 0.5 / GTG 0.35 / TTG 0.15 and the M3-19 wrong-peptide model (``pi_wrong`` 0.05,
  ``search_db = other``); half of the wrong peptides (or wrong de novo reads) not already M-initial
  are rewritten ``M + p[1:]`` when that string occurs in no sample protein (I = L) and not in the
  input. Database (``reverse``) and de novo (``shuffle``) arms through ``run_fdr_suite``: mean FDP
  <= q + 2 SE at every level and list (an unexercised list is a smoke check, as in
  test_placement_fdr_sim.py). Recorded per replicate (``<log_dir>/placement_initiator_i1.jsonl``):
  T0 / D0 (peptides with an exact target / decoy row) with credit (the pipeline's
  ``hits_exact``) and without (``match_exact(starts=None)`` on the same index), the credited
  target and decoy row counts, and the alt-start N-terminal truth peptides placed exactly; and over
  the replicates the exact (conditional binomial, Clopper-Pearson 95 %) CI of the credited
  decoy/target row ratio. Asserted: the arm loads cells (credited target rows > 0).
- **I2, novel** (suite ``novel_orfs``): the M11-16f code-25 G4 arm (``test_m11_16f_novel_arms``,
  PR #481) at its 100-kb size under the default ``annotation.start_codons`` and gippasim's default
  annotation: ``nonatg_nterm_recurrence`` = 0 in every replicate (S5) and the per-class FDP at
  level as before. Skips until #481 is on main.
- **AC 10, non-ATG N-terminal recall**: every 7-40-aa N-terminal tryptic peptide of a non-ATG gene
  of a gippasim genome maps exactly at its own gene (100 %).
"""

from __future__ import annotations

import csv
import dataclasses
import json
import re
import shutil
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _placement_suite import run_fdr_suite
from helpers.initiator_toys import row_set

# I1 and AC 10 run in the placement lane (as test_exact_recall.py); the I2 tests carry NOVEL
pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("placement")]
NOVEL = pytest.mark.calsuite("novel_orfs")
UNSHARDABLE = pytest.mark.calibration_unshardable
LEVELS = C.load_tolerances().levels
STARTS = (("ATG", 0.5), ("GTG", 0.35), ("TTG", 0.15))


# ---------------------------------------------------------------- I1


def _hay(res: Any) -> str:
    with Path(res.proteins_faa).open() as fh:
        text = "".join("|" if line.startswith(">") else line.strip() for line in fh)
    return text.replace("I", "L")


def force_m(src: Path, dst: Path, wrong: set[str], hay: str, seed: int) -> dict[str, str]:
    """Rewrite half of the wrong non-M peptides of ``src`` to ``M + p[1:]``; {old: new}."""
    with src.open(newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    present = {r["peptide"] for r in rows}
    cand = sorted(p for p in wrong & present if p[0] != "M")
    coin = np.random.default_rng([seed, 2020]).random(len(cand)) < 0.5
    sub: dict[str, str] = {}
    for p, take in zip(cand, coin, strict=True):
        new = "M" + p[1:]
        if (
            take
            and new not in present
            and new not in sub.values()
            and new.replace("I", "L") not in hay
        ):
            sub[p] = new
    for r in rows:
        if r["is_decoy"] != "true" and r["peptide"] in sub:
            r["peptide"] = sub[r["peptide"]]
    with dst.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    return sub


def _chain(kind: str, seed: int, tmp: Path) -> Any:
    from helpers.mapping_steps import drive, fresh_steps, preflight, roots_of
    from helpers.placement_sims import Chain
    from helpers.placement_truth import database_truth, denovo_truth

    from gippa.config import resolve_config
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SimConfig, simulate

    denovo = None
    if kind == "denovo":
        from gippasim.denovo import DenovoConfig

        denovo = DenovoConfig()
    pro = ProteomicsConfig(search_db="other", pi_wrong=0.05)  # type: ignore[arg-type]
    res = simulate(
        SimConfig(
            genome_length=100_000,
            start_codons=STARTS,
            proteomics=pro,
            denovo=denovo,
            write_precomputed=True,
        ),
        seed,
        tmp / "sim",
    )
    truth = denovo_truth(res) if kind == "denovo" else database_truth(res)
    src = res.outdir / "peptides_denovo.tsv" if kind == "denovo" else Path(res.peptides)
    tmp.mkdir(parents=True, exist_ok=True)
    sub = force_m(src, tmp / "peptides_m.tsv", truth.wrong, _hay(res), seed)
    for old, new in sub.items():
        truth.add(new, truth.origins.get(old, []), truth.genes.get(old, ()))
        truth.wrong = (truth.wrong - {old}) | {new}
    run = preflight(dataclasses.replace(res, peptides=tmp / "peptides_m.tsv"), tmp / "pf")
    cfg = resolve_config(run.config_yaml, [f"seed={seed}"])
    roots = roots_of(run)
    if kind == "denovo":
        roots["reference_peptides"] = res.outdir / "reference_peptides.tsv"
    with pytest.MonkeyPatch.context() as mp:
        avail = drive(fresh_steps(mp), cfg, tmp / "work", roots)
    return Chain(res, avail, truth, cfg), len(sub)


def credit_record(ch: Any) -> dict[str, float]:
    """T0 / D0 with and without credit; credited target / decoy rows; alt-start N-termini.
    Invariants (statistician): the plain rows are a subset of the credited run's rows, so
    T0_credit >= T0_none and D0_credit >= D0_none."""
    from gippa.mapping.exact import match_exact
    from gippa.mapping.sixframe import load_index
    from gippa.steps.mapping import decoy_types, select_queries

    q = select_queries(pd.read_parquet(ch.avail["peptides"]), ch.cfg)
    index, decoys = load_index(ch.avail["sixframe"], table=11)
    ex = ch.cfg.mapping.exact
    old = match_exact(
        q["peptide"].tolist(),
        index,
        decoys,
        decoy_types(q, ch.cfg),
        min_len=ex.min_len,
        max_len=ex.max_len,
    )
    new = pd.read_parquet(ch.avail["hits_exact"])
    lost = row_set(old) - row_set(new)
    assert not lost, f"credit removed plain exact rows: {sorted(lost)[:5]}"  # rule 3: additive
    extra = row_set(new) - row_set(old)
    tr = pd.read_csv(ch.res.truth["peptides"], sep="\t", dtype=str)
    nterm = set(tr.loc[(tr["alt_start_nterm"] == "true") & (tr["is_wrong"] == "false"), "peptide"])

    def n(df: pd.DataFrame, s: str) -> float:
        return float(df.loc[df["search"] == s, "peptide"].nunique())

    placed = set(new.loc[new["search"] == "target", "peptide"])
    for s in ("target", "decoy"):  # T0 and D0 can only rise with the credit (S3)
        assert n(new, s) >= n(old, s), (s, n(new, s), n(old, s))
    return {
        "T0_credit": n(new, "target"),
        "D0_credit": n(new, "decoy"),
        "T0_none": n(old, "target"),
        "D0_none": n(old, "decoy"),
        "credited_target_rows": float(sum(r[1] == "target" for r in extra)),
        "credited_decoy_rows": float(sum(r[1] == "decoy" for r in extra)),
        "nterm_truth": float(len(nterm & set(q["peptide"]))),
        "nterm_exact": float(len(nterm & placed)),
    }


def poisson_ratio_ci(d: int, t: int, level: float = 0.95) -> tuple[float, float]:
    """Exact CI of the rate ratio D/T given D + T (Clopper-Pearson on D / (D + T))."""
    from scipy.stats import beta

    n, a = d + t, (1 - level) / 2
    lo = 0.0 if d == 0 else float(beta.ppf(a, d, n - d + 1))
    hi = 1.0 if d == n else float(beta.ppf(1 - a, d + 1, n - d))
    return lo / (1 - lo), (hi / (1 - hi) if hi < 1 else float("inf"))


def _i1(kind: str, tmp: Path) -> None:
    from helpers.placement_sims import score_chain

    from gippa.mapping.sixframe import load_start_cells

    def replicate(seed: int, d: Path) -> Any:
        ch, n_forced = _chain(kind, seed, d)
        rep = score_chain(ch, LEVELS)
        rec = {
            "stream": f"placement.m11_20_i1.{kind}",
            "seed": seed,
            "forced_m": n_forced,
            "sites": len(load_start_cells(ch.avail["sixframe"], table=11).sites),
            **credit_record(ch),
        }
        rep.extra["m11_20"] = rec  # type: ignore[assignment]
        out = C.log_dir()
        out.mkdir(parents=True, exist_ok=True)
        with (out / "placement_initiator_i1.jsonl").open("a") as fh:
            fh.write(json.dumps(rec, default=float) + "\n")
        shutil.rmtree(d, ignore_errors=True)
        return rep

    memo = run_fdr_suite(f"placement.m11_20_i1.{kind}", replicate, tmp, ())
    recs = [r.extra["m11_20"] for r in memo.values()]
    if not recs:
        return
    t = int(sum(r["credited_target_rows"] for r in recs))  # type: ignore[index]
    d = int(sum(r["credited_decoy_rows"] for r in recs))  # type: ignore[index]
    summary = {
        "stream": f"placement.m11_20_i1.{kind}.summary",
        "R": len(recs),
        "credited_target_rows": t,
        "credited_decoy_rows": d,
        "ratio_ci95": poisson_ratio_ci(d, t) if t + d else None,
        "mean": {
            k: float(np.mean([r[k] for r in recs]))  # type: ignore[index]
            for k in recs[0]
            if isinstance(recs[0][k], float)
        },
    }  # type: ignore[index]
    with (C.log_dir() / "placement_initiator_i1.jsonl").open("a") as fh:
        fh.write(json.dumps(summary, default=float) + "\n")
    assert t > 0, f"I1 {kind}: no credited target row in {len(recs)} replicates (cells not loaded)"


@UNSHARDABLE
def test_i1_database_placement_fdr_with_credit(tmp_path: Path) -> None:
    _i1("database", tmp_path)


@UNSHARDABLE
def test_i1_denovo_placement_fdr_with_credit(tmp_path: Path) -> None:
    _i1("denovo", tmp_path)


def test_poisson_ratio_ci_brackets_the_ratio() -> None:
    lo, hi = poisson_ratio_ci(10, 20)
    assert lo < 0.5 < hi
    assert poisson_ratio_ci(0, 5)[0] == 0.0


# ---------------------------------------------------------------- I2


def _i2_arms() -> Any:
    from gippa.mapping.initiator import annotated_starts  # noqa: F401  (M11-20: fails until 20i)

    return pytest.importorskip("test_m11_16f_novel_arms", reason="M11-16f arms (#481) not on main")


@NOVEL
def test_i2_code25_nonatg_nterm_recurrence_is_zero() -> None:
    arms = _i2_arms()
    M = arms.M
    R = C.replicates("novel_orfs", False)
    per = []
    for i in range(R):
        rec = arms.record(25, i, "100k")
        per.append(
            sum(
                int(v)
                for reg in M.REGIMES
                for v in rec["diag"][reg].get(M.NONATG_NTERM, {}).values()
            )
        )
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "novel_orfs_initiator_i2.jsonl").open("a") as fh:
        fh.write(
            json.dumps({"stream": "novel_orfs.m11_20_i2.code25", "R": R, M.NONATG_NTERM: per})
            + "\n"
        )
    assert not any(per), {i: n for i, n in enumerate(per) if n}


@NOVEL
@pytest.mark.parametrize("regime", ["picked_decoy", "denovo_subgroup"])
@pytest.mark.parametrize("cls", ["intergenic", "cds_outframe", "antisense", "nterm_extension"])
@pytest.mark.parametrize("q", LEVELS)
def test_i2_code25_class_fdr_at_level(regime: str, cls: str, q: float) -> None:
    arms = _i2_arms()
    R, key = C.replicates("novel_orfs", False), arms.M.g4_key(regime, cls, q)
    rec = C.calibrate_fdr_counts(
        "novel_orfs",
        f"novel_orfs.m11_20_i2.{regime}.{cls}.code25",
        key,
        q,
        lambda i: arms.record(25, i, "100k")["g4"],
        R,
        vacuous={"n_targets": 1.0 / q},
    )
    arms.verdict(rec)


# ---------------------------------------------------------------- AC 10


def test_ac10_nonatg_nterminal_peptides_map_at_their_gene(tmp_path: Path) -> None:
    from gippa.io.fasta import read_genome_fasta
    from gippa.mapping import initiator as ini
    from gippa.mapping.decoys import build_decoys
    from gippa.mapping.exact import match_exact
    from gippa.mapping.sixframe import build_sixframe
    from gippa.model import Replicon
    from gippa.model.genes import gene_model_from_bakta_json
    from gippa.rng import for_step
    from gippasim.simulate import SimConfig, simulate

    found = total = 0
    for seed in (101, 102):
        res = simulate(
            SimConfig(genome_length=100_000, start_codons=STARTS),
            seed=seed,
            outdir=tmp_path / str(seed),
        )
        (rec,) = read_genome_fasta(res.genome_fasta)
        rep = Replicon(rec.id, len(rec.sequence), "circular")
        gm = gene_model_from_bakta_json(res.bakta_json, {rec.id: rep})
        genome = {rec.id: rec.sequence}
        idx = build_sixframe(genome, {rec.id: rep})
        decoys = {"reverse": build_decoys(idx, "reverse", 50, for_step("sixframe_index", seed))}
        sites, _ = ini.annotated_starts(gm, genome, ("ATG", "GTG", "TTG"), table=11)
        tc = ini.start_frame_cells(idx, sites)
        cells = ini.StartCells(
            sites, tc, {"reverse": ini.mirror_start_cells("reverse", decoys["reverse"], tc)}
        )
        want = {}
        for g in gm.genes:
            m = re.search(r"[KR](?!P)", g.aa)
            pep = g.aa[: m.end()] if m else g.aa
            if (
                g.feature_type not in ("cds", "sorf")
                or not g.aa.startswith("M")
                or g.start_codon == "ATG"
                or not 7 <= len(pep) <= 40
            ):
                continue
            iv = g.interval
            s = iv.start if iv.strand == "+" else iv.end - 3 * len(pep)
            want.setdefault(pep, set()).add((s % rep.length, iv.strand))
        df = match_exact(sorted(want), idx, decoys, dict.fromkeys(want, "reverse"), starts=cells)
        tgt = df[df["search"] == "target"]
        got = {
            (p, int(s), st)
            for p, s, st in zip(tgt["peptide"], tgt["start"], tgt["strand"], strict=True)
        }
        total += sum(len(v) for v in want.values())
        found += sum((p, s, st) in got for p, v in want.items() for s, st in v)
    assert total > 50
    assert found == total, (found, total)
