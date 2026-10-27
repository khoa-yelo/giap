"""Test doubles and fixtures for the M8-17b engine switch (docs/specs/M8-17-r-fgsea-engine.md
§7-§11; ADR-0038).

- :class:`RecordingEngine` (T6): stores every ``EngineCall`` and returns scripted rows. By default
  each sent set gets ``EngineRow(size = #members, es 0.5, nes 1.5, p 0.01, log2err 0.25,
  leading_edge = first two members)``; ``rows[(job_id, term)]`` overrides single fields and
  ``drop`` removes rows. It never draws from any generator.
- :class:`ReplayEngine` (T7): serves recorded ``gippa_fgsea.R`` calls from
  ``tests/data/recorded/fgsea/1.36.2/<input>/``. It asserts that ``write_call(call)`` equals the
  recorded input TSVs byte for byte and that the call's options equal ``command.txt`` (so the
  recorded ``n_perm_simple`` 200 is passed, lead note 3), then returns ``read_results`` of the
  recorded ``results.tsv``. ``override`` replaces single fields of the replayed rows.
- :func:`sim_world` builds the M8-06 step fixture exactly as ``scripts/dev/record_fgsea.py``
  ``_world`` does (``world_pre`` + its function DB + the detectability step, cfg seed 7 with
  ``function.gsea.n_perm_simple`` 200) and returns what ``gsea_table`` needs, restricted to the
  four recorded collections (lead note 2: sim_main was recorded for cog, kegg, go_bp, go_slim).
- :func:`expected_sim_main` is the test-author's oracle for A6: the §8 row mapping applied to the
  recorded sim_main output, written with the M8-06 writer. ``python tests/helpers/
  fgsea_switch_data.py`` (PYTHONPATH=src:tests/sim:tests) regenerates the committed file
  ``tests/data/expected/fgsea_replay/sim_main_enrichment_gsea.tsv``.
"""

from __future__ import annotations

import dataclasses
import math
import re
import shlex
import sys
import tempfile
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[2]
RECORDED = REPO / "tests" / "data" / "recorded" / "fgsea" / "1.36.2"
EXPECTED = REPO / "tests" / "data" / "expected" / "fgsea_replay"
EXPECTED_SIM_MAIN = EXPECTED / "sim_main_enrichment_gsea.tsv"
SIM_COLLECTIONS = ("cog", "kegg", "go_bp", "go_slim")
SIM_SEED = 7
SIM_N_PERM = 200
MULTILEVEL = "fgsea_multilevel"


# ---------------------------------------------------------------- RecordingEngine (T6)


@dataclass
class RecordingEngine:
    """A ``GseaEngine`` fake: records calls, returns scripted rows."""

    method: str = MULTILEVEL
    rows: dict[tuple[str, str], dict[str, Any]] = field(default_factory=dict)
    drop: set[tuple[str, str]] = field(default_factory=set)
    calls: list[Any] = field(default_factory=list)
    provenance: Any = None

    def __call__(self, call: Any) -> Any:
        from gippa.function.gsea import EngineRow

        self.calls.append(call)
        out: dict[str, tuple[Any, ...]] = {}
        for job in call.jobs:
            got = []
            for term, members in job.sets.items():
                if (job.job_id, term) in self.drop:
                    continue
                base: dict[str, Any] = {
                    "term": term,
                    "size": len(set(members)),
                    "es": 0.5,
                    "nes": 1.5,
                    "p": 0.01,
                    "log2err": 0.25,
                    "leading_edge": tuple(members[:2]),
                }
                base.update(self.rows.get((job.job_id, term), {}))
                got.append(EngineRow(**base))
            out[job.job_id] = tuple(got)
        return out


class ExplodingEngine:
    """An engine that must never be called."""

    method = MULTILEVEL

    def __call__(self, call: Any) -> Any:
        raise AssertionError(f"the engine was called with jobs {[j.job_id for j in call.jobs]}")


# ---------------------------------------------------------------- ReplayEngine (T7)


def _options(command_txt: str) -> list[str]:
    """``command.txt`` without the script name: the option words."""
    words = shlex.split(command_txt.strip())
    return words[1:]


@dataclass
class ReplayEngine:
    """Replays recorded calls in order; ``dirs[i]`` serves the i-th call."""

    dirs: Sequence[Path]
    override: Mapping[tuple[str, str], Mapping[str, Any]] = field(default_factory=dict)
    method: str = MULTILEVEL
    calls: list[Any] = field(default_factory=list)
    provenance: Any = None

    def __call__(self, call: Any) -> Any:
        from gippa.io.fgsea_r import command, read_results, read_versions, write_call

        i = len(self.calls)
        self.calls.append(call)
        assert i < len(self.dirs), f"unexpected engine call #{i + 1}"
        d = self.dirs[i]
        with tempfile.TemporaryDirectory() as tmp:
            write_call(call, Path(tmp))
            for name in ("stats.tsv", "jobs.tsv", "sets.tsv"):
                got = (Path(tmp) / name).read_bytes()
                want = (d / name).read_bytes()
                assert got == want, f"call #{i + 1}: {name} differs from {d / name}"
        assert command(call)[1:] == _options((d / "command.txt").read_text(encoding="utf-8")), (
            f"call #{i + 1}: options differ from {d / 'command.txt'}"
        )
        self.provenance = read_versions(d / "versions.txt")
        result = read_results(d / "results.tsv", call)
        if not self.override:
            return result
        return {
            job: tuple(
                dataclasses.replace(r, **dict(self.override.get((job, r.term), {}))) for r in rows
            )
            for job, rows in result.items()
        }


# ---------------------------------------------------------------- the M8-06 step fixture


@dataclass(frozen=True)
class SimWorld:
    cfg: Any
    genes: pd.DataFrame
    colls: list[Any]
    replicons: pd.DataFrame
    gene_table: Path
    replicons_path: Path
    db: Any


def sim_cfg() -> Any:
    from gippa.config import GippaConfig

    return GippaConfig.model_validate(
        {
            "seed": SIM_SEED,
            "stats": {"bootstrap": 100},
            "function": {"gsea": {"n_perm_simple": SIM_N_PERM}},
        }
    )


def sim_world(root: Path) -> SimWorld:
    """``record_fgsea.py``'s ``_world`` and its collections, restricted to SIM_COLLECTIONS."""
    from gippa.function.collections import build_collections
    from gippa.genetable.eligibility import eligible
    from gippa.io.function_resources import read_resources
    from gippa.io.genes_tsv import read_genes_parquet
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg
    from helpers.function_db import write_function_db
    from helpers.function_step_data import world_pre, write_gene_table, write_replicons

    load_step_modules()
    cfg = sim_cfg()
    _, xrefs = world_pre()
    db = write_function_db(root / "db", xrefs)
    pre, _ = world_pre(root / "db")
    pre_path = write_gene_table(root / "genetable" / "genes_pre.parquet", pre)
    gene_table = root / "genetable" / "genes.parquet"
    reg.run_step(
        "detectability",
        cfg,
        {"genes_pre": pre_path},
        {"gene_table": gene_table, "detectability_summary": root / "genetable" / "det.parquet"},
    )
    genes = read_genes_parquet(gene_table)
    colls, _ = build_collections(genes, eligible(genes), read_resources(dict(db.inputs())))
    colls = [c for c in colls if c.collection in SIM_COLLECTIONS]
    assert [c.collection for c in colls] == list(SIM_COLLECTIONS), [c.collection for c in colls]
    rep_path = write_replicons(root / "replicons.parquet")
    return SimWorld(cfg, genes, colls, pd.read_parquet(rep_path), gene_table, rep_path, db)


def run_sim_gsea(w: SimWorld, engine: Any) -> pd.DataFrame:
    from gippa.function.gsea import gsea_table
    from gippa.rng import for_step

    fn = w.cfg.function
    return gsea_table(
        w.genes,
        w.colls,
        (fn.gsea_size[0], fn.gsea_size[1]),
        fn.gsea,
        w.cfg.stats.alpha.enrichment,
        for_step("function.gsea", w.cfg.seed),
        w.replicons,
        engine=engine,
    )


# ---------------------------------------------------------------- the A6 oracle


def _recorded_rows(path: Path) -> dict[str, list[dict[str, str]]]:
    lines = path.read_text(encoding="utf-8").split("\n")
    assert lines[-1] == ""
    header = lines[0].split("\t")
    out: dict[str, list[dict[str, str]]] = {}
    for line in lines[1:-1]:
        row = dict(zip(header, line.split("\t"), strict=True))
        out.setdefault(row["job_id"], []).append(row)
    return out


def _num(cell: str) -> float:
    return math.nan if cell == "NA" else float(cell)


def expected_sim_main(w: SimWorld) -> pd.DataFrame:
    """§8 row mapping of the recorded sim_main output, independent of ``gsea_table``.

    No GO term of sim_main reaches q <= alpha (recordings README), so there is no re-test call,
    ``redundant`` is False and ``parent_group`` NA on every row."""
    from gippa.function import gsea as g
    from gippa.stats.multitest import bh

    rank = g.tu_ranking(w.genes)
    det = g._Detection.of(w.genes)
    alpha = w.cfg.stats.alpha.enrichment
    rec = _recorded_rows(RECORDED / "sim_main" / "results.tsv")
    rows: list[dict[str, Any]] = []
    for coll in w.colls:
        got = rec.get(f"main:{coll.collection}", [])
        if not got:
            continue
        p = np.array([_num(r["pval"]) for r in got], dtype=np.float64)
        q = bh(p)
        if coll.collection in g.GO_NAMESPACES:
            assert not np.any(q <= alpha), "sim_main now has a significant GO term"
        tus = g.tu_sets(rank, coll)
        deltas = g.cluster_deltas(
            rank, w.genes, w.replicons, {r["pathway"]: tus[r["pathway"]] for r in got}
        )
        for i, r in enumerate(got):
            term = r["pathway"]
            es, nes, err = _num(r["ES"]), _num(r["NES"]), _num(r["log2err"])
            pv = float(p[i])
            assert int(r["size"]) == len(tus[term])
            edge_tus = [t for t in r["leadingEdge"].split(";") if t] if es != 0 else []
            gs = coll.sets[term]
            edge = [gene for tu in edge_tus for gene in rank.tu_genes[tu] if gene in gs]
            flags = ["ci_unestimated"]
            if math.isnan(pv):
                flags.append("unestimated")
            elif math.isnan(err):
                flags.append("log2err_unestimated")
            if g.is_genomic_cluster(deltas[term]):
                flags.append("genomic_cluster")
            if g.is_detection_confounded(det.z(gs), es):
                flags.append("detection_confounded")
            rows.append(
                {
                    "collection": coll.collection,
                    "term": term,
                    "name": coll.names.get(term, "NA"),
                    "size_tus": len(tus[term]),
                    "es": es,
                    "nes": nes,
                    "log2err": err,
                    "leading_edge": edge,
                    "redundant": False,
                    "parent_group": None,
                    "coverage": float(coll.coverage),
                    "stat": es,
                    "effect": nes,
                    "ci_low": math.nan,
                    "ci_high": math.nan,
                    "p_value": pv,
                    "q_value": float(q[i]),
                    "method": MULTILEVEL,
                    "null_model": "tu_set_permutation",
                    "n_perm": SIM_N_PERM,
                    "flags": flags,
                }
            )
    df = pd.DataFrame(rows, columns=list(g.COLUMNS))
    for c in (
        "es",
        "nes",
        "log2err",
        "coverage",
        "stat",
        "effect",
        "ci_low",
        "ci_high",
        "p_value",
        "q_value",
    ):
        df[c] = df[c].astype("float64")
    for c in ("size_tus", "n_perm"):
        df[c] = df[c].astype("int64")
    df["redundant"] = df["redundant"].astype("boolean")
    return df


def write_expected(path: Path = EXPECTED_SIM_MAIN) -> Path:
    from gippa.io.function_tables import write_enrichment_gsea_tsv

    with tempfile.TemporaryDirectory() as tmp:
        w = sim_world(Path(tmp))
        df = expected_sim_main(w)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_enrichment_gsea_tsv(path, df)
    return path


SCHEMA = REPO / "src" / "gippa" / "schemas" / "tables" / "enrichment_gsea.schema.json"


def validate_gsea(path: Path) -> None:
    """``path`` validates against ``enrichment_gsea.schema.json`` (frictionless, as M8-06a)."""
    import json

    import pytest

    frictionless = pytest.importorskip("frictionless")
    d = json.loads(SCHEMA.read_text(encoding="utf-8"))
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(
        path=path.name, basepath=str(path.parent), schema=frictionless.Schema.from_descriptor(d)
    )
    report = res.validate()
    assert report.valid, report.flatten(["rowNumber", "fieldName", "type", "note"])


def flag_set(cell: Any) -> set[str]:
    if isinstance(cell, str):
        return set() if cell in ("", "NA") else set(cell.split(";"))
    return {str(x) for x in cell}


_JOB = re.compile(r"^(main):([a-z_]+)$|^(retest):([a-z_]+):(.+)$")


def parse_job_id(job_id: str) -> tuple[str, ...]:
    m = _JOB.match(job_id)
    assert m is not None, job_id
    return tuple(x for x in m.groups() if x is not None)


if __name__ == "__main__":
    print(write_expected(Path(sys.argv[1]) if len(sys.argv) > 1 else EXPECTED_SIM_MAIN))
