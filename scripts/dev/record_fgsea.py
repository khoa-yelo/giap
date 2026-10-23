"""Record real ``giap_fgsea.R`` calls for tests/data/recorded/fgsea/1.36.2/ (M8-17a §6; §11.4).

Usage: python scripts/dev/record_fgsea.py <outdir> [name ...]   (PYTHONPATH=src:tests/sim:tests)
Runs only where the giap-fgsea env is built (``scripts/slurm/record_fgsea.sbatch``). Every
recording is one call of the script through ``giap.io.fgsea_r`` (``write_call``, ``command``,
``child_env``) in ``<outdir>/<name>/``, which keeps ``stats.tsv``, ``jobs.tsv``, ``sets.tsv``,
``command.txt``, ``results.tsv``, ``versions.txt`` and ``PROVENANCE.md``. All inputs are synthetic:

- ``sim_main``: the main call of the M8-06 step fixture (``tests/unit/steps/test_function_step.py``
  ``World``: ``world_pre`` + its function DB, the detectability step, cfg seed 7 with
  ``function.gsea.n_perm_simple = 200``), built as M8-17 §5 and §8 describe: the ranking in
  ``for_step("function.gsea", 7).permutation(n)`` order, the seed the next
  ``integers(1, 2**31 - 1)`` draw, one ``main:<collection>`` job per collection with in-range
  sets (``function.gsea_size``, fewer than n TUs), members in natural TU order;
- ``sim_retest``: that run's GO re-test call shape (§8: a term re-tested on its set minus the
  genes of significant descendants, ``[1, n]``). No fixture GO term reaches q <= alpha, so the
  significance is forced: each tested go_bp term with tested descendants loses the genes of its
  lowest-p tested descendant;
- ``ties``: the T5 ranking (``tests/contract/test_fgsea_r_contract.py`` ``_ties_call``);
- ``unbalanced``: 499 steeply decaying positive scores and one negligible negative one; a
  set at the bottom gets ``pval`` NA;
- ``log2err_na`` (extra): 1 negative and 499 uniform positive scores; the bottom set gets a tiny
  finite p with ``log2err`` NA (``isCpGeHalf``);
- ``tiny_p``: 1000 TUs with a 25-TU set packed at the top, so p < 1e-20 with a finite log2err.
"""

from __future__ import annotations

import hashlib
import os
import platform
import shlex
import subprocess
import sys
import tempfile
from collections.abc import Sequence
from datetime import UTC, datetime
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
SIM_COLLECTIONS = ("cog", "kegg", "go_bp", "go_slim")
NAMES = ("sim_main", "sim_retest", "ties", "unbalanced", "tiny_p", "log2err_na")


def _call(ids: Sequence[str], scores: Sequence[float], jobs: list, seed: int, params=None):
    from giap.config.model import GseaConfig
    from giap.function.gsea import EngineCall

    return EngineCall(tuple(ids), np.asarray(scores, dtype=np.float64), tuple(jobs),
                      params or GseaConfig(), seed)  # fmt: skip


def run(call, d: Path, note: str) -> None:
    """One real call in ``d``; fails loudly on a non-zero exit."""
    from giap.io import fgsea_r as fr

    launcher = fr.resolve_launcher(os.environ)
    d.mkdir(parents=True, exist_ok=True)
    for p in d.iterdir():
        p.unlink()
    fr.write_call(call, d)
    cmd = fr.command(call)
    proc = subprocess.run(cmd, cwd=d, env=fr.child_env(launcher, os.environ), capture_output=True,
                          text=True, timeout=1800, check=False)  # fmt: skip
    if proc.returncode != 0:
        raise SystemExit(f"{d.name}: exit {proc.returncode}\n{proc.stderr}")
    (d / "command.txt").write_text(shlex.join([Path(cmd[0]).name, *cmd[1:]]) + "\n")
    assert fr.read_versions(d / "versions.txt").fgsea == fr.FGSEA_VERSION
    fr.read_results(d / "results.tsv", call)  # the recording parses
    lock = REPO / "workflow" / "envs" / "locks" / "fgsea.explicit.txt"
    head = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"],
                          capture_output=True, text=True, check=False).stdout.strip()  # fmt: skip
    (d / "PROVENANCE.md").write_text(
        f"# {d.name}: real fgsea {fr.FGSEA_VERSION} output (giap_fgsea.R)\n\n"
        f"- Input: {note}\n"
        f"- Command (in this directory, giap-fgsea env at {launcher.path}): `"
        f"{(d / 'command.txt').read_text().strip()}`\n"
        f"- Env lock workflow/envs/locks/fgsea.explicit.txt sha256 "
        f"{hashlib.sha256(lock.read_bytes()).hexdigest()}\n"
        f"- Date {datetime.now(UTC).strftime('%Y-%m-%dT%H:%M:%SZ')}, host {platform.node()}, "
        f"Slurm job {os.environ.get('SLURM_JOB_ID', 'none')}, giap {head}\n"
        f"- Recorded by scripts/dev/record_fgsea.py (scripts/slurm/record_fgsea.sbatch).\n",
        encoding="utf-8",
    )
    print(f"{d.name}: {len((d / 'results.tsv').read_text().splitlines()) - 1} rows")


# ---------------------------------------------------------------- synthetic rankings


def ties_call():
    from giap.function.gsea import GseaJob

    tied = [f"T{i:02d}" for i in range(20)]
    rest = [f"U{i:02d}" for i in range(30)]
    scores = [2.0] * 20 + [1.5 - 0.1 * i for i in range(30)]
    members = tuple(f"T{i:02d}" for i in range(0, 20, 2))
    other = tuple(f"U{i:02d}" for i in range(0, 30, 3))
    return _call(tied + rest, scores, [GseaJob("main:ties", {"S": members, "O": other}, 5, 49)],
                 20261022)  # fmt: skip


def _natural(ids: Sequence[str], order: Sequence[str]) -> tuple[str, ...]:
    keep = set(ids)
    return tuple(t for t in order if t in keep)


def _bottom_top_call(scores: Sequence[float], job: str, seed: int):
    """500 TUs (``scores`` in rank order, the last one negative); sets BOTTOM (the last 15),
    TOP (the first 15) and two random 20-sets, members in natural order."""
    from giap.function.gsea import GseaJob

    r = np.random.default_rng(20261023)
    natural = [f"TU_{i:03d}" for i in range(1, len(scores) + 1)]
    ids = list(r.permutation(natural))
    bottom = _natural(ids[-15:], natural)
    top = _natural(ids[:15], natural)
    rand = {f"R{k}": _natural(r.choice(ids, 20, replace=False).tolist(), natural) for k in (1, 2)}
    sets = {"BOTTOM": bottom, "TOP": top, **rand}
    return _call(ids, scores, [GseaJob(job, sets, 10, 500)], seed)


def unbalanced_call():
    """Steeply decaying positive scores and a negligible negative one: a random set's first hit
    dominates its ES, so fewer than 10 of the nPermSimple random sets have ES <= 0 (fgsea's
    ``modeFraction < 10``) and BOTTOM (ES < 0) gets pval, NES and log2err NA."""
    scores = [*(4.0 * 0.9 ** np.arange(499)).tolist(), -1e-25]
    return _bottom_top_call(scores, "main:unbalanced", 424242)


def log2err_na_call():
    """Uniform positive scores and one negative: BOTTOM gets a tiny finite p whose log2err is
    NA (fgsea's ``isCpGeHalf``, "likely overestimated"; M8-17 example 4)."""
    r = np.random.default_rng(20261025)
    scores = [*np.sort(r.uniform(0.05, 4.0, 499))[::-1].tolist(), -0.5]
    return _bottom_top_call(scores, "main:log2err_na", 424242)


def tiny_p_call():
    from giap.function.gsea import GseaJob

    r = np.random.default_rng(20261024)
    n = 1000
    natural = [f"TU_{i:04d}" for i in range(1, n + 1)]
    ids = list(r.permutation(natural))
    scores = np.sort(r.normal(0.0, 1.0, n))[::-1] + np.where(np.arange(n) < 25, 3.0, 0.0)
    sets = {"PACKED": _natural(ids[:25], natural)}
    sets |= {f"R{k}": _natural(r.choice(ids, 40, replace=False).tolist(), natural)
             for k in (1, 2, 3)}  # fmt: skip
    return _call(ids, scores.tolist(), [GseaJob("main:tiny_p", sets, 10, 500)], 777)


# ---------------------------------------------------------------- the M8-06 step fixture


def _world(root: Path):
    """The step fixture's gene table (``World`` in tests/unit/steps/test_function_step.py)."""
    from helpers.function_db import write_function_db
    from helpers.function_step_data import world_pre, write_gene_table

    from giap.config import GiapConfig
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    cfg = GiapConfig.model_validate({"seed": 7, "stats": {"bootstrap": 100},
                                     "function": {"gsea": {"n_perm_simple": 200}}})  # fmt: skip
    _, xrefs = world_pre()
    db = write_function_db(root / "db", xrefs)
    pre, _ = world_pre(root / "db")
    pre_path = write_gene_table(root / "genetable" / "genes_pre.parquet", pre)
    gene_table = root / "genetable" / "genes.parquet"
    reg.run_step("detectability", cfg, {"genes_pre": pre_path},
                 {"gene_table": gene_table,
                  "detectability_summary": root / "genetable" / "det.parquet"})  # fmt: skip
    return cfg, gene_table, db


def sim_calls(root: Path, out: Path, names: Sequence[str]) -> None:
    from giap.function.collections import build_collections
    from giap.function.gsea import GO_NAMESPACES, GseaJob, _closure, tu_ranking, tu_sets
    from giap.genetable.eligibility import eligible
    from giap.io.fgsea_r import read_results
    from giap.io.function_resources import read_resources
    from giap.io.genes_tsv import read_genes_parquet
    from giap.rng import for_step

    cfg, gene_table, db = _world(root)
    genes = read_genes_parquet(gene_table)
    colls, _ = build_collections(genes, eligible(genes), read_resources(dict(db.inputs())))
    rank = tu_ranking(genes)
    n = len(rank.tu_ids)
    rng = for_step("function.gsea", cfg.seed)
    perm = rng.permutation(n)
    seed = int(rng.integers(1, 2**31 - 1))
    ids = [rank.tu_ids[int(i)] for i in perm]
    scores = np.asarray(rank.scores, dtype=np.float64)[perm]
    lo, hi = cfg.function.gsea_size
    colls = [c for c in colls if c.collection in SIM_COLLECTIONS]
    by_coll = {c.collection: c for c in colls}
    jobs = []
    for coll in colls:
        tus = tu_sets(rank, coll)
        sets = {t: _natural(s, rank.tu_ids) for t, s in tus.items()
                if lo <= len(s) <= hi and len(s) < n}  # fmt: skip
        if sets:
            jobs.append(GseaJob(f"main:{coll.collection}", sets, lo, hi))
    main = _call(ids, scores, jobs, seed, cfg.function.gsea)
    note = (
        "M8-06 step fixture (World: world_pre, function DB, detectability step; seed 7, "
        "function.gsea.n_perm_simple 200), main call per M8-17 §5 and §8, collections "
        f"{', '.join(SIM_COLLECTIONS)} (§6; the fixture's go_mf and go_cc are left out)"
    )
    run(main, out / "sim_main", note)
    if "sim_retest" not in names:
        return
    result = read_results(out / "sim_main" / "results.tsv", main)
    tu_of = {g: tu for tu, gs in rank.tu_genes.items() for g in gs}
    retests = []
    for job in jobs:
        coll = by_coll[job.job_id.removeprefix("main:")]
        if coll.collection not in GO_NAMESPACES:
            continue
        rows = result[job.job_id]
        p = {r.term: r.p for r in rows}
        children: dict[str, set[str]] = {}
        for t, ps in coll.parents.items():
            for parent in ps:
                children.setdefault(parent, set()).add(t)
        below = _closure(children)
        for r in rows:
            # Forced significance (no fixture GO term has q <= alpha): each tested term with
            # tested descendants is re-tested without its lowest-p tested descendant's genes.
            tested = sorted(below(r.term) & set(p), key=lambda d: (p[d], d))
            if not tested:
                continue
            desc = tested[:1]
            removed = frozenset().union(*(coll.sets[d] for d in desc))
            tus = {tu_of[g] for g in coll.sets[r.term] - removed if g in tu_of}
            if tus:
                key = f"retest:{coll.collection}:{r.term}"
                retests.append(GseaJob(key, {r.term: _natural(tus, rank.tu_ids)}, 1, n))
    print(f"sim_main: GO re-tests {[j.job_id for j in retests]}")
    if not retests:
        raise SystemExit("sim_retest: no GO re-test job")
    retest_note = note.replace("main call", "GO re-test call") + (
        "; sim_main has no GO term with q <= alpha, so each tested go_bp term with tested "
        "descendants is re-tested without the genes of its lowest-p tested descendant "
        "(forced significance), on [1, n] with the same ranking and seed"
    )
    run(_call(ids, scores, retests, seed, cfg.function.gsea), out / "sim_retest", retest_note)


def main(argv: list[str]) -> int:
    out = Path(argv[1])
    names = argv[2:] or list(NAMES)
    bad = sorted(set(names) - set(NAMES))
    if bad:
        raise SystemExit(f"unknown recording {bad}; valid: {NAMES}")
    makers = {
        "ties": (ties_call, "the T5 tie ranking: 50 TUs, 20 tied at the top"),
        "unbalanced": (unbalanced_call, "1 negative and 499 steeply decaying positive scores"),
        "tiny_p": (tiny_p_call, "1000 TUs, 25-TU set packed at the top"),
        "log2err_na": (log2err_na_call, "1 negative and 499 uniform positive scores"),
    }
    for name in names:
        if name in makers:
            make, note = makers[name]
            run(make(), out / name, note)
    if {"sim_main", "sim_retest"} & set(names):
        with tempfile.TemporaryDirectory() as tmp:
            sim_calls(Path(tmp), out, names)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
