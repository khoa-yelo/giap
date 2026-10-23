"""``denovo_only_unsupported`` through the ``novel_orfs`` step
(docs/specs/M7-17-denovo-support-flag.md acceptance criterion 2): examples 6, 12 and 13, the
side-table column (name, dtype, position, null without ``peptides``), the
``novel_orfs.denovo_support`` INFO event with its exact ``data`` and message, and no new WARNING.

Inputs: the giapsim SIM100K local run of ``helpers.novel_steps_fixture`` (as
tests/unit/steps/test_novel_class_fdr_step.py), its step rerun on hand-built ``peptides``: each
candidate ORF's own peptides as ``database`` or ``denovo`` C9 rows.

Readings pinned by the test author:

* the oracle of the flag reads the written ``novel_orfs.tsv`` (``class``, ``rbs_supported``) and
  the sources the test wrote; it is the spec's rule 4 (``intergenic``/``nterm_extension`` never
  exempt);
* the side-table cell is ``False`` (not null) on every row of a run with ``peptides``;
* a run with ``peptides`` but no de novo inputs (``hits_exact`` ...) keeps its M7-11b
  ``layer_missing`` WARNINGs; "no WARNING is added" is checked as: no WARNING record carries the
  new event or names the flag.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import novel_steps_fixture as nf
from helpers.novel_fdr_fixture import db_id, pep_row, peptides_frame

FLAG = "denovo_only_unsupported"
EVENT = "novel_orfs.denovo_support"
RBS_CLASSES = frozenset({"intergenic", "nterm_extension"})
NEW_INPUTS = ("peptides", "features", "hits_exact", "hits_tolerant", "placement_costs")
OLD_SIDE = ("orf_id", "class", "fdr_class", "regime", "q", "class_fdr", "reason",
            "class_mismatch", "class_rule_met", "called")  # fmt: skip
SIDE_COLUMNS = (*OLD_SIDE, FLAG)
TABLE_KEYS = ("novel_orfs", "novel_orfs_gff3", "novel_orfs_faa", "annotation_patch")
DATA_KEYS = {"n_orfs", "n_denovo_only", "n_flagged", "n_called", "n_flagged_called"}


# ---------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def run1(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, Path, Path]:
    root = tmp_path_factory.mktemp("m7_17_sim100k")
    sim = nf.simulate(root / "sim")
    out = root / "OUT"
    with nf.db_env(root):
        code, log = nf.run_giap(sim, out)
    assert code == 0, log
    return sim, out, out / "work"


@pytest.fixture(scope="module")
def cfg1(run1: tuple[Any, Path, Path]) -> Any:
    from giap.config import resolve_config

    return resolve_config(Path(run1[0].config_yaml), [f"seed={nf.SEED}"])


@pytest.fixture(scope="module")
def orfs(run1: tuple[Any, Path, Path]) -> list[Any]:
    from giap.steps.novel import read_candidates

    got, _ = read_candidates(run1[2] / "novel_candidates")
    assert len(got) >= 4, got
    return sorted(got, key=lambda o: o.orf_id)


def _layer_path(sim: Any, cfg: Any, key: str) -> Path:
    p = Path(getattr(cfg.annotation.precomputed, key))
    return p if p.is_absolute() else Path(sim.config_yaml).parent / p


def _base(run1: tuple[Any, Path, Path], cfg: Any, **extra: Any) -> dict[str, Any]:
    sim, _, work = run1
    prom = _layer_path(sim, cfg, "promoters")
    ins = nf.step_inputs("novel_orfs", work, {
        "rbs_novel": _layer_path(sim, cfg, "rbs_novel"), "novel_promoter_windows": prom,
        "promoter_windows": prom, "terminators": _layer_path(sim, cfg, "terminators"),
    }, skip=NEW_INPUTS)  # fmt: skip
    ins.update(extra)
    return ins


class _Rec(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _run(cfg: Any, ins: Mapping[str, Any], out: Path) -> tuple[dict, list[logging.LogRecord]]:
    rec = _Rec()
    log = logging.getLogger("giap.steps")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        outs = nf.run("novel_orfs", cfg, ins, out)
    finally:
        log.removeHandler(rec)
        log.setLevel(old)
    return outs, rec.records


def _events(records: list[logging.LogRecord]) -> list[logging.LogRecord]:
    return [r for r in records if getattr(r, "event", None) == EVENT]


def _flags(row: Mapping[str, str | None]) -> list[str]:
    f = row["flags"]
    return [] if f in (None, "", "NA") else str(f).split(";")


def _write(df: pd.DataFrame, path: Path) -> Path:
    df.to_parquet(path, index=False)
    return path


def _peptides(orfs: list[Any], src: Mapping[str, Mapping[str, str]], tmp: Path,
              extra: list[dict[str, Any]] = ()) -> Path:  # type: ignore[assignment]  # fmt: skip
    """C9 rows: ORF ``o``'s peptide ``p`` with source ``src[o.orf_id][p]`` (default database);
    a peptide shared by two ORFs keeps its first source."""
    rows, seen, k = [], set(), 0
    for o in orfs:
        for p in o.peptides:
            if p in seen:
                continue
            seen.add(p)
            k += 1
            s = src.get(o.orf_id, {}).get(p, "database")
            ids = [o.overlapping_gene or "b9999"]
            rows.append(pep_row(p, ids, k, is_decoy=False, source=s,
                                q=math.nan if s == "denovo" else None))  # fmt: skip
    return _write(peptides_frame([*rows, *extra]), tmp / "peptides.parquet")


def _expected(cls: str, has_dn: bool, db: bool, rbs: str | None) -> bool:
    return has_dn and not db and not (rbs == "true" and cls not in RBS_CLASSES)


def _pair(o: Any, src: Mapping[str, Mapping[str, str]], first: Mapping[str, str]) -> tuple:
    """(has_denovo, db_support) from the sources the test wrote (``first`` = peptide -> source
    actually written)."""
    s = {first[p] for p in o.peptides if p in first}
    return "denovo" in s, "database" in s


def _written_sources(orfs: list[Any], src: Mapping[str, Mapping[str, str]]) -> dict[str, str]:
    out: dict[str, str] = {}
    for o in orfs:
        for p in o.peptides:
            out.setdefault(p, src.get(o.orf_id, {}).get(p, "database"))
    return out


# ---------------------------------------------------------------- example 12 (no peptides)


def test_example_12_without_peptides_no_flag_null_cell_no_event(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    _, out, _ = run1
    o, rec = _run(cfg1, _base(run1, cfg1), tmp_path / "a")
    schema = pq.read_schema(o["novel_class_fdr"])
    assert tuple(schema.names) == SIDE_COLUMNS
    assert pa.types.is_boolean(schema.field(FLAG).type)
    side = pd.read_parquet(o["novel_class_fdr"])
    assert side[FLAG].isna().all(), side[FLAG].tolist()
    assert not _events(rec)
    for r in nf.tsv_dicts(o["novel_orfs"]):
        assert FLAG not in _flags(r), r
    # the published table (database peptides only) carries no flag either and is unchanged
    assert o["novel_orfs"].read_bytes() == (out / "tables" / "novel_orfs.tsv").read_bytes()


# ---------------------------------------------------------------- database peptides


def test_run_own_database_peptides_never_flag(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    _, _, work = run1
    a, _ = _run(cfg1, _base(run1, cfg1), tmp_path / "a")
    o, rec = _run(cfg1, _base(run1, cfg1, peptides=work / "peptides",
                              features=work / "features"), tmp_path / "o")  # fmt: skip
    for key in TABLE_KEYS:
        assert a[key].read_bytes() == o[key].read_bytes(), key
    schema = pq.read_schema(o["novel_class_fdr"])
    assert tuple(schema.names) == SIDE_COLUMNS
    assert pa.types.is_boolean(schema.field(FLAG).type)
    side = pd.read_parquet(o["novel_class_fdr"])
    tsv = nf.tsv_dicts(o["novel_orfs"])
    assert list(side["orf_id"]) == [r["orf_id"] for r in tsv]
    assert side[FLAG].notna().all() and not side[FLAG].astype(bool).any()
    (ev,) = _events(rec)
    assert ev.levelno == logging.INFO
    n_called = sum(r["called"] == "true" for r in tsv)
    assert dict(ev.data) == {"n_orfs": len(tsv), "n_denovo_only": 0, "n_flagged": 0,  # type: ignore[attr-defined]
                             "n_called": n_called, "n_flagged_called": 0}  # fmt: skip


def test_example_6_picked_decoy_orf_with_database_peptides_not_flagged(
    run1: tuple[Any, Path, Path], cfg1: Any, orfs: list[Any], tmp_path: Path
) -> None:
    """A makedb run: an ORF on its six-frame entry with a losing decoy is ``picked_decoy``; its
    database peptides clear the flag (side cell false)."""
    _, _, work = run1
    six = [o for o in orfs if o.cls in ("intergenic", "cds_outframe", "antisense")]
    assert six, [o.cls for o in orfs]
    o0 = six[0]
    entry = (db_id(o0.cls, o0.contig, o0.stop - 300, o0.stop, "+") if o0.strand == "+"
             else db_id(o0.cls, o0.contig, o0.stop + 3, o0.stop + 303, "-"))  # fmt: skip
    rows = [pep_row(p, [entry], k, is_decoy=False) for k, p in enumerate(o0.peptides, 1)]
    rows.append(pep_row("DECYAAAAK", ["rev_" + entry], 900, is_decoy=True))
    others = [o for o in orfs if o is not o0]
    peps = _peptides(others, {}, tmp_path, extra=rows)
    out, rec = _run(cfg1, _base(run1, cfg1, peptides=peps, features=work / "features"),
                    tmp_path / "o")  # fmt: skip
    side = pd.read_parquet(out["novel_class_fdr"]).set_index("orf_id")
    assert side.loc[o0.orf_id, "regime"] == "picked_decoy"
    assert side.loc[o0.orf_id, FLAG] == False  # noqa: E712
    assert not side[FLAG].astype(bool).any()
    assert len(_events(rec)) == 1


# ---------------------------------------------------------------- example 13 (mixed sources)


def _mixed_sources(orfs: list[Any]) -> dict[str, dict[str, str]]:
    """Every ORF de novo except: the first ORF with >= 2 peptides gets one database peptide
    (example 5), the last ORF all database."""
    src = {o.orf_id: dict.fromkeys(o.peptides, "denovo") for o in orfs}
    multi = [o for o in orfs if len(o.peptides) >= 2]
    assert multi, "no ORF with two peptides"
    src[multi[0].orf_id][multi[0].peptides[0]] = "database"
    src[orfs[-1].orf_id] = dict.fromkeys(orfs[-1].peptides, "database")
    return src


def test_example_13_mixed_run_flag_side_cell_and_event(
    run1: tuple[Any, Path, Path], cfg1: Any, orfs: list[Any], tmp_path: Path
) -> None:
    _, _, work = run1
    src = _mixed_sources(orfs)
    first = _written_sources(orfs, src)
    peps = _peptides(orfs, src, tmp_path)
    out, rec = _run(cfg1, _base(run1, cfg1, peptides=peps, features=work / "features"),
                    tmp_path / "o")  # fmt: skip
    tsv = nf.tsv_dicts(out["novel_orfs"])
    side = pd.read_parquet(out["novel_class_fdr"])
    assert list(side["orf_id"]) == [r["orf_id"] for r in tsv]
    by_id = {o.orf_id: o for o in orfs}
    n = {k: 0 for k in DATA_KEYS}
    for r, cell in zip(tsv, side[FLAG], strict=True):
        o = by_id[str(r["orf_id"])]
        dn, db = _pair(o, src, first)
        want = _expected(str(r["class"]), dn, db, r["rbs_supported"])
        assert (FLAG in _flags(r)) is want, (r, dn, db)
        assert bool(cell) is want and not pd.isna(cell), (r["orf_id"], cell)
        if FLAG in _flags(r):
            assert _flags(r)[-1] == FLAG, r
        called = r["called"] == "true"
        n["n_orfs"] += 1
        n["n_denovo_only"] += int(dn and not db)
        n["n_flagged"] += int(want)
        n["n_called"] += int(called)
        n["n_flagged_called"] += int(want and called)
    assert n["n_flagged"] >= 1 and n["n_denovo_only"] >= n["n_flagged"], n
    (ev,) = _events(rec)
    assert ev.levelno == logging.INFO
    assert set(ev.data) == DATA_KEYS  # type: ignore[attr-defined]
    assert dict(ev.data) == n  # type: ignore[attr-defined]
    assert ev.getMessage() == (
        f"novel_orfs: {n['n_flagged_called']} of {n['n_called']} calls are de novo-only without "
        "database or RBS support (ADR-0050: exact de novo errors are outside the class FDR)"
    )
    for r in rec:
        if r.levelno >= logging.WARNING:
            assert getattr(r, "event", None) != EVENT, r.getMessage()
            assert FLAG not in r.getMessage(), r.getMessage()


def test_no_candidate_orf_event_counts_zero(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    from giap.io.parquet import write_parquet_deterministic
    from giap.steps.novel import _CAND_COLS

    _, _, work = run1
    empty = tmp_path / "novel_candidates.parquet"
    write_parquet_deterministic(empty, _CAND_COLS, [])
    o, rec = _run(cfg1, _base(run1, cfg1, novel_candidates=empty, peptides=work / "peptides",
                              features=work / "features"), tmp_path / "o")  # fmt: skip
    side = pd.read_parquet(o["novel_class_fdr"])
    assert len(side) == 0 and tuple(side.columns) == SIDE_COLUMNS
    (ev,) = _events(rec)
    assert dict(ev.data) == dict.fromkeys(DATA_KEYS, 0)  # type: ignore[attr-defined]
