"""The wired ``novel_orfs`` class FDR (docs/specs/M7-11b-class-fdr-wiring.md acceptance
criterion 2): the new optional inputs (``peptides``, ``features``, ``hits_exact``,
``hits_tolerant``, ``placement_costs``) and output (``novel_class_fdr``), M7-11 decisions 7-9
through the step (pinned by M7-11 decision 10), examples 5, 6, 10 and 12, the three events and
two WARNINGs with exact ``data``, the missing-input degenerate rows, the side-table columns,
dtypes and row order, ``placement_costs`` from ``placement``, and the published side table of a
``giap run``.

Inputs: one giapsim SIM100K local run (``helpers.novel_steps_fixture``; database peptides with
``rev_<locus_tag>`` decoys, nine candidate ORFs). The step is rerun on its work files with
hand-built ``peptides`` / hit tables / ``placement_costs``: the candidate ORFs' own peptides as
makedb-tagged target rows (entry = the six-frame ORF with the same ``(contig, strand, stop)``,
``nterm_extension`` = the run's own ``nterm_entries`` record), or as de novo rows.

Readings pinned by the test author:

* example 6's "byte-identical to main" is checked as: the four files are the same bytes with no
  new input (the pass-1 path), with ``peptides`` without decoy rows (``no_decoys``) and with the
  run's own ``rev_<locus_tag>`` peptides (``decoys_not_makedb``); and the run's published
  ``tables/novel_orfs.tsv`` equals the no-new-input rerun;
* ``class_rule_met`` is checked against the pass-1 ``called`` (= rule and >= 2 peptides) for
  ORFs with >= 2 peptides;
* NaN fields of ``novel_orfs.class_fdr`` data (``c`` of a picked class) may be None or NaN in
  the log record (``events.jsonl`` writes null either way);
* the floor WARNING text is checked for its fixed parts only (number formats are open).
"""

from __future__ import annotations

import logging
import math
from collections.abc import Iterator, Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import novel_steps_fixture as nf
from helpers.class_call_fixture import COST_COLUMNS, costs_frame, decoy_hits, target_hits
from helpers.novel_fdr_fixture import db_id, pep_row, peptides_frame

from giap.model.coords import Interval

NEW_INPUTS = ("peptides", "features", "hits_exact", "hits_tolerant", "placement_costs")
SIDE_COLUMNS = ("orf_id", "class", "fdr_class", "regime", "q", "class_fdr", "reason",
                "class_mismatch", "class_rule_met", "called")  # fmt: skip
TABLE_KEYS = ("novel_orfs", "novel_orfs_gff3", "novel_orfs_faa", "annotation_patch")
CLASS_EVENT = "novel_orfs.class_fdr"
FLOOR_EVENT = "novel_orfs.class_q_floor_above_alpha"
MAKEDB_EVENT = "novel_orfs.decoys_not_makedb"
SUMMARY_KEYS = {"regime", "class", "estimated", "reason", "n_targets", "n_target_wins",
                "n_decoys", "paired_share", "c", "q_floor", "floor_above_alpha", "alpha",
                "n_called"}  # fmt: skip
ALPHA = 0.01


# ---------------------------------------------------------------- fixtures


@pytest.fixture(scope="module")
def run1(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, Path, Path]:
    root = tmp_path_factory.mktemp("m7_11b_sim100k")
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
def world(run1: tuple[Any, Path, Path]) -> dict[str, Any]:
    from giap.io.fasta import read_genome_fasta
    from giap.makedb.build import _nterm_entries
    from giap.steps.mapping import read_gene_model
    from giap.steps.novel import read_candidates

    _, _, work = run1
    orfs, _ = read_candidates(work / "novel_candidates")
    gm = read_gene_model(work / "genes", work / "replicons", work / "features")
    genome = {r.id: r.sequence for r in read_genome_fasta(work / "genome")}
    nterm = {e.locus_tag: e.db_id for e in _nterm_entries(gm, genome)}
    by_cls: dict[str, list[Any]] = {}
    for o in sorted(orfs, key=lambda o: o.orf_id):
        by_cls.setdefault(o.cls, []).append(o)
    nts = [o for o in by_cls.get("nterm_extension", []) if o.overlapping_gene in nterm]
    assert by_cls.get("cds_outframe") and by_cls.get("intergenic"), sorted(by_cls)
    assert len(by_cls.get("antisense", [])) >= 2 and len(nts) >= 2, sorted(by_cls)
    return {"orfs": orfs, "by_cls": by_cls, "nterm": nterm, "nts": nts, "gm": gm,
            "genome": genome}  # fmt: skip


def _layer_path(sim: Any, cfg: Any, key: str) -> Path:
    p = Path(getattr(cfg.annotation.precomputed, key))
    return p if p.is_absolute() else Path(sim.config_yaml).parent / p


def _base(run1: tuple[Any, Path, Path], cfg: Any, **extra: Any) -> dict[str, Any]:
    """The inputs the local executor gives ``novel_orfs`` (precomputed layers), without the new
    keys, plus ``extra``."""
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


def _run(name: str, cfg: Any, ins: Mapping[str, Any], out: Path) -> tuple[dict, list]:
    rec = _Rec()
    log = logging.getLogger("giap.steps")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        outs = nf.run(name, cfg, ins, out)
    finally:
        log.removeHandler(rec)
        log.setLevel(old)
    return outs, rec.records


def _events(records: list[logging.LogRecord], name: str) -> list[logging.LogRecord]:
    return [r for r in records if getattr(r, "event", None) == name]


def _side(outs: Mapping[str, Path]) -> pd.DataFrame:
    return pd.read_parquet(outs["novel_class_fdr"])


def _na(v: Any) -> Any:
    """A missing scalar (None, NaN, NA) as None; ``q`` keeps its NaN."""
    return None if v is None or v is pd.NA or (isinstance(v, float) and math.isnan(v)) else v


def _rows(outs: Mapping[str, Path]) -> dict[str, dict[str, Any]]:
    """Side-table rows by ``orf_id``; missing ``reason`` as None."""
    out = {}
    for r in _side(outs).to_dict("records"):
        out[r["orf_id"]] = {k: (v if k == "q" else _na(v)) for k, v in r.items()}
    return out


def _tsv(outs: Mapping[str, Path]) -> list[dict[str, str | None]]:
    return nf.tsv_dicts(outs["novel_orfs"])


def _write(df: pd.DataFrame, path: Path) -> Path:
    df.to_parquet(path, index=False)
    return path


def _entry(o: Any, tag: str | None = None) -> str:
    """The makedb six-frame entry closed by ``o``'s stop codon."""
    cls = tag or o.cls
    if o.strand == "+":
        return db_id(cls, o.contig, o.stop - 300, o.stop, "+")
    return db_id(cls, o.contig, o.stop + 3, o.stop + 303, "-")


def _features(work: Path) -> Path:
    return work / "features"


def _ranks() -> Iterator[int]:
    k = 0
    while True:
        k += 1
        yield k


def _flags(row: Mapping[str, str | None]) -> list[str]:
    f = row["flags"]
    return [] if f in (None, "", "NA") else str(f).split(";")


def _close(a: Any, b: float) -> bool:
    return a is not None and not (isinstance(a, float) and math.isnan(a)) and math.isclose(a, b)


def _nan_or_none(v: Any) -> bool:
    return v is None or (isinstance(v, float) and math.isnan(v))


# ---------------------------------------------------------------- registration


def test_registration_adds_optional_inputs_and_the_side_table() -> None:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    load_step_modules()
    spec = reg.get("novel_orfs")
    assert set(spec.inputs) == {"novel_candidates", "gene_table", "genes", "replicons", "genome"}
    assert set(spec.optional_inputs) == {
        "novel_promoter_windows",
        "rbs_novel",
        "promoter_windows",
        "terminators",
        *NEW_INPUTS,
    }
    assert set(spec.outputs) == {*TABLE_KEYS, "novel_class_fdr"}
    assert "placement_costs" in reg.get("placement").outputs


# ---------------------------------------------------------------- example 6 + side-table shape


def test_example_6_no_decoys_byte_identical_and_side_table(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    _, out, work = run1
    pep = pd.read_parquet(work / "peptides")
    clean = _write(pep[~pep["is_decoy"].astype(bool)].reset_index(drop=True),
                   tmp_path / "peptides.parquet")  # fmt: skip
    a, ra = _run("novel_orfs", cfg1, _base(run1, cfg1), tmp_path / "a")
    b, rb = _run("novel_orfs", cfg1, _base(run1, cfg1, peptides=clean,
                                            features=_features(work)), tmp_path / "b")  # fmt: skip
    for key in TABLE_KEYS:
        assert a[key].read_bytes() == b[key].read_bytes(), key
    assert a["novel_orfs"].read_bytes() == (out / "tables" / "novel_orfs.tsv").read_bytes()
    for side in (_side(a), _side(b)):
        assert set(side["regime"]) <= {"unestimated"} and set(side["reason"]) <= {"no_decoys"}
        assert len(side) == len(_tsv(a))
    warn = {getattr(r, "event", None) for r in ra if r.levelno >= logging.WARNING}
    warn_b = {getattr(r, "event", None) for r in rb if r.levelno >= logging.WARNING}
    assert warn_b == warn  # no new WARNING
    assert not _events(rb, CLASS_EVENT)


def test_side_table_columns_dtypes_and_row_order(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    import pyarrow as pa
    import pyarrow.parquet as pq

    _, _, work = run1
    o, _ = _run("novel_orfs", cfg1, _base(run1, cfg1, peptides=work / "peptides",
                                          features=_features(work)), tmp_path / "o")  # fmt: skip
    schema = pq.read_schema(o["novel_class_fdr"])
    assert tuple(schema.names) == SIDE_COLUMNS
    t = {f.name: f.type for f in schema}
    for c in ("orf_id", "class", "fdr_class", "regime", "class_fdr", "reason"):
        assert pa.types.is_string(t[c]) or pa.types.is_large_string(t[c]), (c, t[c])
    assert pa.types.is_float64(t["q"])
    for c in ("class_mismatch", "class_rule_met", "called"):
        assert pa.types.is_boolean(t[c]), (c, t[c])
    side = _side(o)
    tsv = _tsv(o)
    assert list(side["orf_id"]) == [r["orf_id"] for r in tsv]
    assert list(side["class"]) == [r["class"] for r in tsv]
    assert list(side["class_fdr"]) == [r["class_fdr"] for r in tsv]
    assert [bool(x) for x in side["called"]] == [r["called"] == "true" for r in tsv]


# ---------------------------------------------------------------- example 5 (run's own decoys)


def test_example_5_rev_locus_decoys_are_not_makedb(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    _, _, work = run1
    pep = pd.read_parquet(work / "peptides")
    n_dec = int(pep["is_decoy"].astype(bool).sum())
    assert n_dec > 0
    a, _ = _run("novel_orfs", cfg1, _base(run1, cfg1), tmp_path / "a")
    c, rc = _run("novel_orfs", cfg1, _base(run1, cfg1, peptides=work / "peptides",
                                            features=_features(work)), tmp_path / "c")  # fmt: skip
    for key in TABLE_KEYS:
        assert a[key].read_bytes() == c[key].read_bytes(), key
    side = _side(c)
    assert set(side["regime"]) == {"unestimated"}
    assert set(side["reason"]) == {"decoys_not_makedb"}
    assert set(side["class_fdr"]) == {"unestimated"}
    ev = _events(rc, MAKEDB_EVENT)
    assert len(ev) == 1 and ev[0].levelno == logging.WARNING
    assert dict(ev[0].data) == {"n_decoy_rows": n_dec}  # type: ignore[attr-defined]
    assert not _events(rc, CLASS_EVENT)
    assert not [r for r in rc if "layer_missing:placement_costs" in r.getMessage()]


def test_no_candidate_orf_empty_side_table_and_makedb_warning(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    from helpers.novel_outputs_fixture import NOVEL_ORFS_COLUMNS

    from giap.io.parquet import write_parquet_deterministic
    from giap.steps.novel import _CAND_COLS

    _, _, work = run1
    empty = tmp_path / "novel_candidates.parquet"
    write_parquet_deterministic(empty, _CAND_COLS, [])
    o, rec = _run("novel_orfs", cfg1, _base(run1, cfg1, novel_candidates=empty,
                                            peptides=work / "peptides",
                                            features=_features(work)), tmp_path / "o")  # fmt: skip
    assert o["novel_orfs"].read_text() == "\t".join(NOVEL_ORFS_COLUMNS) + "\n"
    side = _side(o)
    assert len(side) == 0 and tuple(side.columns) == SIDE_COLUMNS
    assert len(_events(rec, MAKEDB_EVENT)) == 1
    assert not _events(rec, CLASS_EVENT)


# ---------------------------------------------------------------- decisions 7-9 (picked)


def _makedb_peptides(w: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Target rows on makedb entries: cds_outframe (its decoy loses), the intergenic ORF on a
    ``cds_outframe``-tagged entry (decision 9), antisense[0] beaten by its decoy, antisense[1]
    without decoy, nterm[0] on its gene's ``nterm_entries`` record (paired by locus tag),
    nterm[1] on an ID with no such record (decision 7), the rest on locus tags."""
    bc, rk = w["by_cls"], _ranks()
    rows: list[dict[str, Any]] = []
    role: dict[str, Any] = {}

    def targets(o: Any, pid: str) -> None:
        for p in o.peptides:
            rows.append(pep_row(p, [pid], next(rk), is_decoy=False))

    cds = bc["cds_outframe"]
    for o in cds:
        targets(o, _entry(o))
    ig = bc["intergenic"][0]
    targets(ig, _entry(ig, "cds_outframe"))
    a0, a1 = bc["antisense"][:2]
    rows.append(pep_row("DECYAAAAK", ["rev_" + _entry(a0)], 0.5, is_decoy=True))
    targets(a0, _entry(a0))
    targets(a1, _entry(a1))
    for o in cds:
        rows.append(pep_row(f"DECY{o.stop}K", ["rev_" + _entry(o)], 900, is_decoy=True))
    n0, n1 = w["nts"][:2]
    targets(n0, w["nterm"][n0.overlapping_gene])
    targets(n1, f"GIAPORF_nterm_extension_{n1.contig}_1_99_p")
    done = {o.orf_id for o in [*cds, ig, a0, a1, n0, n1]}
    for o in w["orfs"]:
        if o.orf_id not in done:
            targets(o, o.overlapping_gene or "b9999")
    role.update(cds=cds, ig=ig, a0=a0, a1=a1, n0=n0, n1=n1, done=done)
    return peptides_frame(rows), role


def _check_call_rule(
    side: list[dict[str, Any]], n_pep: Mapping[str, int], alpha: float = ALPHA
) -> None:
    for r in side:
        if r["reason"] is None:
            assert math.isfinite(r["q"]), r
            want = r["class_rule_met"] and r["q"] <= alpha
        else:
            assert math.isnan(r["q"]), r
            want = False if r["reason"] == "lost_pick" else (
                r["class_rule_met"] and n_pep[r["orf_id"]] >= 2)  # fmt: skip
        assert bool(r["called"]) is bool(want), r


def test_decisions_7_8_9_through_the_step(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path
) -> None:
    _, _, work = run1
    peps, role = _makedb_peptides(world)
    p = _write(peps, tmp_path / "peptides.parquet")
    o, rec = _run("novel_orfs", cfg1, _base(run1, cfg1, peptides=p, features=_features(work)),
                  tmp_path / "o")  # fmt: skip
    side = _rows(o)
    tsv = {r["orf_id"]: r for r in _tsv(o)}
    assert set(side) == {x.orf_id for x in world["orfs"]}
    assert {r["regime"] for r in side.values()} == {"picked_decoy"}
    t_cds = len(role["cds"]) + 1
    for x in role["cds"]:
        r = side[x.orf_id]
        assert (r["fdr_class"], r["reason"], bool(r["class_mismatch"])) == (
            "cds_outframe", None, False)  # fmt: skip
        assert r["q"] == pytest.approx(1 / t_cds)
    ig = side[role["ig"].orf_id]  # decision 9
    assert (ig["class"], ig["fdr_class"], bool(ig["class_mismatch"])) == (
        "intergenic", "cds_outframe", True)  # fmt: skip
    assert ig["q"] == pytest.approx(1 / t_cds)  # from the cds_outframe list
    assert tsv[role["ig"].orf_id]["class"] == "intergenic"
    assert "class_mismatch" not in _flags(tsv[role["ig"].orf_id])  # lead decision 3
    a0 = side[role["a0"].orf_id]
    assert (a0["reason"], a0["class_fdr"], bool(a0["called"])) == ("lost_pick", "1", False)
    assert "class_fdr_unestimated" not in _flags(tsv[role["a0"].orf_id])
    assert side[role["a1"].orf_id]["reason"] is None
    n0, n1 = side[role["n0"].orf_id], side[role["n1"].orf_id]  # decision 7
    assert (n0["fdr_class"], n0["reason"]) == ("nterm_extension", None)
    assert math.isfinite(n0["q"])
    assert (n1["reason"], n1["class_fdr"]) == ("no_makedb_pair", "unestimated")
    assert "class_fdr_unestimated" in _flags(tsv[role["n1"].orf_id])
    for oid in set(side) - role["done"]:
        assert side[oid]["reason"] == "untagged_class", oid
    for oid, r in side.items():
        assert r["class_fdr"] == tsv[oid]["class_fdr"], oid
        if r["reason"] is None:
            assert tsv[oid]["class_fdr"] == f"{r['q']:.6g}", oid
    _check_call_rule(list(side.values()), {x.orf_id: len(x.peptides) for x in world["orfs"]})
    assert not _events(rec, MAKEDB_EVENT)


def test_picked_events_and_floor_warning_data(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path
) -> None:
    _, _, work = run1
    peps, role = _makedb_peptides(world)
    p = _write(peps, tmp_path / "peptides.parquet")
    _, rec = _run("novel_orfs", cfg1, _base(run1, cfg1, peptides=p, features=_features(work)),
                  tmp_path / "o")  # fmt: skip
    evs = {(e.data["regime"], e.data["class"]): e for e in _events(rec, CLASS_EVENT)}  # type: ignore[attr-defined]
    assert all(e.levelno == logging.INFO for e in evs.values())
    assert len(evs) == len(_events(rec, CLASS_EVENT))  # one per summary
    cds = evs[("picked_decoy", "cds_outframe")].data  # type: ignore[attr-defined]
    assert set(cds) == SUMMARY_KEYS
    t = len(role["cds"]) + 1
    want = {"regime": "picked_decoy", "class": "cds_outframe", "estimated": True, "reason": None,
            "n_targets": t, "n_target_wins": t, "n_decoys": 0, "paired_share": 1.0,
            "floor_above_alpha": True, "alpha": ALPHA, "n_called": 0}  # fmt: skip
    assert {k: cds[k] for k in want} == want
    assert _close(cds["q_floor"], 1 / t) and _nan_or_none(cds["c"])
    anti = evs[("picked_decoy", "antisense")].data  # type: ignore[attr-defined]
    assert anti["n_targets"] == 2 and anti["n_target_wins"] == 1 and anti["n_decoys"] == 1
    floors = {(e.data["regime"], e.data["class"]): e for e in _events(rec, FLOOR_EVENT)}  # type: ignore[attr-defined]
    assert set(floors) == {k for k, e in evs.items() if e.data["floor_above_alpha"]}  # type: ignore[attr-defined]
    f = floors[("picked_decoy", "cds_outframe")]
    assert f.levelno == logging.WARNING
    assert set(f.data) == {"regime", "class", "n_targets", "q_floor", "alpha"}  # type: ignore[attr-defined]
    assert {k: f.data[k] for k in ("regime", "class", "n_targets", "alpha")} == {  # type: ignore[attr-defined]
        "regime": "picked_decoy",
        "class": "cds_outframe",
        "n_targets": t,
        "alpha": ALPHA,
    }
    assert _close(f.data["q_floor"], 1 / t)  # type: ignore[attr-defined]
    msg = f.getMessage()
    assert msg.startswith(f"novel_orfs: class cds_outframe (picked_decoy) has {t} targets; "
                          "its q floor "), msg  # fmt: skip
    assert "> alpha" in msg and "so none can be called" in msg, msg
    assert "a class with fewer than about 1/alpha = " in msg, msg
    assert "candidates cannot be called at alpha = " in msg, msg


# ---------------------------------------------------------------- de novo (examples 10, 12)


def _denovo_inputs(
    w: dict[str, Any],
    tmp: Path,
    *,
    dn: list[Any],
    pi: float = math.nan,
    e_target: Mapping[str, float] | None = None,
    decoy_rows: list[tuple[str, Interval]] = (),  # type: ignore[assignment]
    e_decoy: Mapping[str, float] | None = None,
    db_rows: list[dict[str, Any]] = (),  # type: ignore[assignment]
) -> dict[str, Path]:
    """``dn`` ORFs' peptides as de novo C9 rows with target hits and ``placement_costs`` rows
    (E 1e-3 unless ``e_target`` says otherwise, π̂ ``pi``); every other ORF's peptides as
    database rows on its locus tag (no decoy rows unless ``db_rows`` adds some)."""
    dn_ids = {o.orf_id for o in dn}
    dn_peps = sorted({p for o in dn for p in o.peptides})
    rows = [pep_row(p, [], 0, is_decoy=False, source="denovo", q=math.nan, score=40.0)
            for p in dn_peps]  # fmt: skip
    rk = _ranks()
    for o in w["orfs"]:
        if o.orf_id not in dn_ids:
            for p in o.peptides:
                rows.append(pep_row(p, [o.overlapping_gene or "b9999"], next(rk), is_decoy=False))
    rows += list(db_rows)
    et = [(e_target or {}).get(p, 1e-3) for p in dn_peps]
    d_peps = [p for p, _ in decoy_rows]
    ed = [(e_decoy or {}).get(p, 1e-2) for p in dn_peps]
    costs = costs_frame(
        dn_peps,
        e_target=et,
        e_decoy=[ed[i] if p in d_peps else math.nan for i, p in enumerate(dn_peps)],
        pi=pi,
    )
    for p in d_peps:
        costs.loc[costs["peptide"] == p, "m_decoy"] = 0.0
    hx = target_hits(
        [p for p, e in zip(dn_peps, et, strict=True) if math.isfinite(e)],
        [Interval("chr1", 100, 130, "+")] * sum(math.isfinite(e) for e in et),
    )
    ht = decoy_hits(d_peps, [iv for _, iv in decoy_rows])
    return {
        "peptides": _write(peptides_frame(rows), tmp / "peptides.parquet"),
        "hits_exact": _write(hx, tmp / "hits_exact.parquet"),
        "hits_tolerant": _write(ht, tmp / "hits_tolerant.parquet"),
        "placement_costs": _write(costs, tmp / "placement_costs.parquet"),
    }


def _dn_ins(run1: tuple[Any, Path, Path], cfg: Any, files: dict[str, Path],
            drop: tuple[str, ...] = ()) -> dict[str, Any]:  # fmt: skip
    files = {**files, "features": _features(run1[2])}
    return _base(run1, cfg, **{k: v for k, v in files.items() if k not in drop})


def test_example_10_denovo_floor_warning(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path
) -> None:
    anti = world["by_cls"]["antisense"]
    n_t = len({p for o in anti for p in o.peptides})
    files = _denovo_inputs(world, tmp_path, dn=anti)  # π̂ missing: c = 1
    o, rec = _run("novel_orfs", cfg1, _dn_ins(run1, cfg1, files), tmp_path / "o")
    side = _rows(o)
    for x in anti:
        r = side[x.orf_id]
        assert (r["regime"], r["fdr_class"], r["reason"]) == ("denovo_subgroup", "antisense",
                                                              None)  # fmt: skip
        assert not bool(r["called"])
    ev = {(e.data["regime"], e.data["class"]): e.data for e in _events(rec, CLASS_EVENT)}  # type: ignore[attr-defined]
    d = ev[("denovo_subgroup", "antisense")]
    assert set(d) == SUMMARY_KEYS
    assert (d["n_targets"], d["n_decoys"], d["estimated"]) == (n_t, 0, True)
    assert _close(d["c"], 1.0) and _close(d["q_floor"], 1 / n_t)
    assert _nan_or_none(d["n_target_wins"]) and _nan_or_none(d["paired_share"])
    floors = [
        e
        for e in _events(rec, FLOOR_EVENT)
        if e.data["regime"] == "denovo_subgroup" and e.data["class"] == "antisense"
    ]  # type: ignore[attr-defined]
    assert len(floors) == 1
    f = floors[0].data  # type: ignore[attr-defined]
    assert {k: f[k] for k in ("regime", "class", "n_targets", "alpha")} == {
        "regime": "denovo_subgroup",
        "class": "antisense",
        "n_targets": n_t,
        "alpha": ALPHA,
    }
    assert _close(f["q_floor"], 1 / n_t)
    for x in world["orfs"]:
        if x not in anti:
            assert side[x.orf_id]["reason"] == "no_decoys", x.orf_id


@pytest.mark.parametrize(("pi_per_peptide", "called"), [(0.0100000001, False), (0.0099999, True)])
def test_call_compares_q_as_float_before_formatting(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path,
    pi_per_peptide: float, called: bool,
) -> None:  # fmt: skip
    """q = c / T with no decoys; c = π̂ = ``pi_per_peptide`` * T, so q formats as ``0.01`` (or
    ``0.0099999``) and only the float comparison decides (rule 6)."""
    cds = world["by_cls"]["cds_outframe"]
    n_t = len({p for o in cds for p in o.peptides})
    assert all(len(o.peptides) >= 2 for o in cds)  # the cds_outframe rule holds
    files = _denovo_inputs(world, tmp_path, dn=cds, pi=pi_per_peptide * n_t)
    o, _ = _run("novel_orfs", cfg1, _dn_ins(run1, cfg1, files), tmp_path / "o")
    side = _rows(o)
    tsv = {r["orf_id"]: r for r in _tsv(o)}
    for x in cds:
        r = side[x.orf_id]
        assert r["q"] == pytest.approx(pi_per_peptide, rel=1e-6)
        assert bool(r["class_rule_met"]) is True
        assert bool(r["called"]) is called, (r["q"], r["class_fdr"])
        assert tsv[x.orf_id]["called"] == ("true" if called else "false")
        assert tsv[x.orf_id]["class_fdr"] == ("0.01" if not called else "0.0099999")


def test_example_12_mixed_run_one_event_per_regime(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path
) -> None:
    cds = world["by_cls"]["cds_outframe"]
    anti = world["by_cls"]["antisense"]
    db = [pep_row("DECYAAAAK", ["rev_" + _entry(cds[0])], 900, is_decoy=True)]
    files = _denovo_inputs(world, tmp_path, dn=anti, pi=0.5, db_rows=db)
    # cds_outframe ORFs on their makedb entries
    peps = pd.read_parquet(files["peptides"])
    for x in cds:
        for p in x.peptides:
            i = peps.index[peps["peptide"] == p]
            peps.loc[i, "protein_ids"] = pd.Series([[_entry(x)]] * len(i), index=i)
    files["peptides"] = _write(peps, tmp_path / "peptides_mixed.parquet")
    o, rec = _run("novel_orfs", cfg1, _dn_ins(run1, cfg1, files), tmp_path / "o")
    side = _rows(o)
    for x in cds:
        assert side[x.orf_id]["regime"] == "picked_decoy"
    for x in anti:
        assert side[x.orf_id]["regime"] == "denovo_subgroup"
    keys = {(e.data["regime"], e.data["class"]) for e in _events(rec, CLASS_EVENT)}  # type: ignore[attr-defined]
    assert ("picked_decoy", "cds_outframe") in keys
    assert ("denovo_subgroup", "antisense") in keys
    assert not any(r == "denovo_subgroup" and c == "cds_outframe" for r, c in keys)
    d = next(
        e.data
        for e in _events(rec, CLASS_EVENT)  # type: ignore[attr-defined]
        if e.data["regime"] == "denovo_subgroup" and e.data["class"] == "antisense"
    )  # type: ignore[attr-defined]
    assert d["n_targets"] == len({p for x in anti for p in x.peptides})


@pytest.mark.parametrize("missing", ["placement_costs", "hits_exact", "hits_tolerant"])
def test_denovo_inputs_missing_warns_layer_missing(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path, missing: str
) -> None:
    anti = world["by_cls"]["antisense"]
    files = _denovo_inputs(world, tmp_path, dn=anti, pi=0.5)
    o, rec = _run("novel_orfs", cfg1, _dn_ins(run1, cfg1, files, drop=(missing,)),
                  tmp_path / "o")  # fmt: skip
    side = _rows(o)
    for x in anti:
        r = side[x.orf_id]
        assert (r["regime"], r["reason"], r["class_fdr"]) == (
            "denovo_subgroup", "denovo_inputs_missing", "unestimated")  # fmt: skip
    lm = [
        r
        for r in _events(rec, "novel_orfs.layer_missing")
        if r.data.get("flag") == f"layer_missing:{missing}"
    ]  # type: ignore[attr-defined]
    assert len(lm) == 1 and lm[0].levelno == logging.WARNING
    assert not [e for e in _events(rec, CLASS_EVENT) if e.data["regime"] == "denovo_subgroup"]  # type: ignore[attr-defined]


def test_denovo_peptides_without_target_e_are_no_target_hit(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path
) -> None:
    a0, a1 = world["by_cls"]["antisense"][:2]
    files = _denovo_inputs(world, tmp_path, dn=[a0, a1], pi=0.5,
                           e_target={p: math.nan for p in a0.peptides})  # fmt: skip
    o, _ = _run("novel_orfs", cfg1, _dn_ins(run1, cfg1, files), tmp_path / "o")
    side = _rows(o)
    assert side[a0.orf_id]["reason"] == "no_target_hit"
    assert side[a0.orf_id]["class_fdr"] == "unestimated"
    assert side[a1.orf_id]["reason"] is None


def test_nan_e_decoy_fails_the_step(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path
) -> None:
    from giap.steps.registry import StepError

    a0 = world["by_cls"]["antisense"][0]
    p0 = a0.peptides[0]
    files = _denovo_inputs(world, tmp_path, dn=[a0], pi=0.5,
                           decoy_rows=[(p0, Interval("chr1", 5000, 5030, "-"))],
                           e_decoy={p0: math.nan})  # fmt: skip
    with pytest.raises(StepError) as ei:
        nf.run("novel_orfs", cfg1, _dn_ins(run1, cfg1, files), tmp_path / "o")
    assert isinstance(ei.value.cause, ValueError)


# ---------------------------------------------------------------- placement_costs


def test_placement_writes_its_best_costs(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    import pyarrow.parquet as pq

    _, _, work = run1
    outs = nf.run("placement", cfg1, nf.step_inputs("placement", work), tmp_path / "p")
    path = outs["placement_costs"]
    assert tuple(pq.read_schema(path).names) == COST_COLUMNS
    costs = pd.read_parquet(path)
    assert list(costs["peptide"]) == sorted(costs["peptide"])
    assert costs["peptide"].is_unique
    pep = pd.read_parquet(work / "peptides")
    assert set(costs["peptide"]) <= set(pep.loc[~pep["is_decoy"].astype(bool), "peptide"])
    per = pd.read_parquet(outs["placement"]).set_index("peptide")
    c = costs.set_index("peptide")
    assert set(per.index) <= set(c.index)
    for col in ("m_target", "e_target", "length"):
        pd.testing.assert_series_equal(
            c.loc[per.index, col].astype("float64"), per[col].astype("float64"), check_names=False
        )
    assert c["pi_denovo"].isna().all()  # a database run: no de novo π̂
    for col in ("m_target", "m_decoy", "e_target", "e_decoy", "pep", "q_value", "pi_denovo"):
        assert costs[col].dtype == "float64", col


# ---------------------------------------------------------------- the giap run


def test_giap_run_publishes_the_side_table_and_the_makedb_warning(
    run1: tuple[Any, Path, Path],
) -> None:
    import json

    _, out, work = run1
    assert (work / "placement_costs").exists()
    side = pd.read_parquet(out / "novel" / "novel_class_fdr.parquet")
    tsv = nf.tsv_dicts(out / "tables" / "novel_orfs.tsv")
    assert list(side["orf_id"]) == [r["orf_id"] for r in tsv]
    assert set(side["reason"]) == {"decoys_not_makedb"}
    pep = pd.read_parquet(work / "peptides")
    evs = [json.loads(ln) for ln in (out / "logs" / "events.jsonl").read_text().splitlines()
           if ln.strip()]  # fmt: skip
    mk = [e for e in evs if e.get("event") == MAKEDB_EVENT]
    assert len(mk) == 1 and mk[0]["level"] == "WARNING"
    assert mk[0]["data"] == {"n_decoy_rows": int(pep["is_decoy"].astype(bool).sum())}
    man = json.loads((out / "run_manifest.json").read_text())
    assert f"novel_orfs:{MAKEDB_EVENT}" in man["warnings"]
    assert man["novel_orfs"] is None  # no class event in a decoys_not_makedb run
    assert not any("decoys_not_makedb" in f for f in man["flags"])


# ---------------------------------------------------------------- review round (#256 B1, B2)
#
# B1: ``n_called`` of a ``novel_orfs.class_fdr`` event on a class with calls (mutant M11). With
# nine candidate ORFs a picked list has q >= 1 / 9, so the picked test runs at
# ``stats.alpha.novel_orfs = 0.5``; the de novo test keeps alpha 0.01 (q = c / T < alpha).
# B2: picked ``n_targets`` = won + lost_pick (mutant M8). ``orf_regimes`` and
# ``picked_class_fdr`` look an ORF's peptides up in the same target rows, so a picked-regime ORF
# always has a target score through the step; the test forces one ORF without target rows into
# the picked regime (``giap.novel.fdr.orf_regimes`` wrapped), with its makedb entry kept in the
# run by an unrelated target peptide on the same ID.


def _n_called_by_list(side: Mapping[str, Mapping[str, Any]]) -> dict[tuple[str, str], int]:
    out: dict[tuple[str, str], int] = {}
    for r in side.values():
        if bool(r["called"]):
            k = (str(r["regime"]), str(r["fdr_class"]))
            out[k] = out.get(k, 0) + 1
    return out


def _mismatch_peptides(w: dict[str, Any]) -> tuple[pd.DataFrame, dict[str, Any]]:
    """``_makedb_peptides`` with antisense[1] (its class rule holds) on the ``cds_outframe``-tagged
    entry closed by its stop instead of its own entry (decision 9)."""
    peps, role = _makedb_peptides(w)
    a1 = role["a1"]
    own, tagged = _entry(a1), _entry(a1, "cds_outframe")
    peps["protein_ids"] = [[tagged if i == own else i for i in ids] for ids in peps["protein_ids"]]
    return peps, role


def test_b1_picked_n_called_counts_calls_under_the_entry_class(
    run1: tuple[Any, Path, Path], world: dict[str, Any], tmp_path: Path
) -> None:
    from giap.config import resolve_config

    sim, _, work = run1
    cfg = resolve_config(Path(sim.config_yaml), [f"seed={nf.SEED}", "stats.alpha.novel_orfs=0.75"])
    assert cfg.stats.alpha.novel_orfs == 0.75
    peps, role = _mismatch_peptides(world)
    p = _write(peps, tmp_path / "peptides.parquet")
    o, rec = _run("novel_orfs", cfg, _base(run1, cfg, peptides=p, features=_features(work)),
                  tmp_path / "o")  # fmt: skip
    side = _rows(o)
    a1 = side[role["a1"].orf_id]
    assert (a1["class"], a1["fdr_class"], bool(a1["class_mismatch"])) == (
        "antisense", "cds_outframe", True)  # fmt: skip
    assert bool(a1["class_rule_met"]) and bool(a1["called"]), a1  # called from cds_outframe
    by_list = _n_called_by_list(side)
    t_cds = len(role["cds"]) + 2  # the cds_outframe ORFs, ig and a1 (both on tagged entries)
    called_cds = sum(bool(side[x.orf_id]["called"]) for x in [*role["cds"], role["ig"]])
    assert called_cds >= 1
    assert by_list[("picked_decoy", "cds_outframe")] == called_cds + 1
    evs = {(e.data["regime"], e.data["class"]): e.data for e in _events(rec, CLASS_EVENT)}  # type: ignore[attr-defined]
    cds = evs[("picked_decoy", "cds_outframe")]
    assert (cds["n_targets"], cds["alpha"], cds["floor_above_alpha"]) == (t_cds, 0.75, False)
    assert cds["n_called"] == called_cds + 1 >= 2
    anti = evs[("picked_decoy", "antisense")]
    assert (anti["n_targets"], anti["n_called"]) == (1, 0)  # a0 lost_pick; a1 is not counted here
    for k, d in evs.items():
        assert d["n_called"] == by_list.get(k, 0), k
    tsv_called = sum(r["called"] == "true" for r in _tsv(o))
    untagged = sum(bool(r["called"]) for r in side.values() if r["reason"] == "untagged_class")
    assert sum(d["n_called"] for d in evs.values()) == tsv_called - untagged


def test_b1_denovo_n_called_counts_the_called_class(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path
) -> None:
    cds = world["by_cls"]["cds_outframe"]
    n_t = len({p for o in cds for p in o.peptides})
    files = _denovo_inputs(world, tmp_path, dn=cds, pi=0.0099999 * n_t)  # q = 0.0099999
    o, rec = _run("novel_orfs", cfg1, _dn_ins(run1, cfg1, files), tmp_path / "o")
    side = _rows(o)
    assert all(bool(side[x.orf_id]["called"]) for x in cds)
    d = next(e.data for e in _events(rec, CLASS_EVENT)  # type: ignore[attr-defined]
             if (e.data["regime"], e.data["class"]) == ("denovo_subgroup", "cds_outframe"))  # type: ignore[attr-defined]  # fmt: skip
    assert d["n_called"] == len(cds) >= 1
    assert d["n_called"] == _n_called_by_list(side)[("denovo_subgroup", "cds_outframe")]


def test_b2_picked_n_targets_excludes_no_target_score(
    run1: tuple[Any, Path, Path], cfg1: Any, world: dict[str, Any], tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:  # fmt: skip
    from giap.novel import fdr as F

    _, _, work = run1
    peps, role = _makedb_peptides(world)
    ig, a1 = role["ig"], role["a1"]
    gone = set(ig.peptides) | set(a1.peptides)
    peps = peps[~peps["peptide"].isin(gone)].reset_index(drop=True)
    # keep both entries in the run (an unrelated target peptide on each ID), no ORF score on them
    extra = [pep_row("GGGGGGGGK", [_entry(ig, "cds_outframe")], 950, is_decoy=False),
             pep_row("GGGGGGGGGK", [_entry(a1)], 951, is_decoy=False)]  # fmt: skip
    peps = pd.concat([peps, peptides_frame(extra)], ignore_index=True)
    p = _write(peps, tmp_path / "peptides.parquet")
    forced = {ig.orf_id, a1.orf_id}
    orig = F.orf_regimes

    def regimes(*a: Any, **k: Any) -> dict[str, Any]:
        out = orig(*a, **k)
        assert {out[i] for i in forced} == {"unestimated"}  # no target row of their own
        return {i: ("picked_decoy" if i in forced else r) for i, r in out.items()}

    monkeypatch.setattr(F, "orf_regimes", regimes)
    o, rec = _run("novel_orfs", cfg1, _base(run1, cfg1, peptides=p, features=_features(work)),
                  tmp_path / "o")  # fmt: skip
    side = _rows(o)
    for oid, cls in ((ig.orf_id, "cds_outframe"), (a1.orf_id, "antisense")):
        r = side[oid]
        assert (r["regime"], r["fdr_class"], r["reason"], r["class_fdr"]) == (
            "picked_decoy", cls, "no_target_score", "unestimated"), r  # fmt: skip
    a0 = side[role["a0"].orf_id]
    assert a0["reason"] == "lost_pick"
    evs = {(e.data["regime"], e.data["class"]): e.data for e in _events(rec, CLASS_EVENT)}  # type: ignore[attr-defined]
    for cls in ("cds_outframe", "antisense"):
        d = evs[("picked_decoy", cls)]
        rows = [r for r in side.values() if (r["regime"], r["fdr_class"]) == ("picked_decoy", cls)]
        won = sum(1 for r in rows if r["reason"] is None)
        lost = sum(1 for r in rows if r["reason"] == "lost_pick")
        assert d["estimated"] is True, (cls, d)
        assert d["n_targets"] == won + lost, (cls, d, won, lost)
    t = len(role["cds"])  # the cds_outframe ORFs win; ig has no target score
    cds = evs[("picked_decoy", "cds_outframe")]
    assert (cds["n_targets"], cds["n_target_wins"], cds["floor_above_alpha"]) == (t, t, True)
    anti = evs[("picked_decoy", "antisense")]
    assert (anti["n_targets"], anti["n_target_wins"], anti["n_decoys"]) == (1, 0, 1)
    floors = {(e.data["regime"], e.data["class"]): e for e in _events(rec, FLOOR_EVENT)}  # type: ignore[attr-defined]
    f = floors[("picked_decoy", "cds_outframe")]
    assert f.data["n_targets"] == t  # type: ignore[attr-defined]
    assert f.getMessage().startswith(
        f"novel_orfs: class cds_outframe (picked_decoy) has {t} targets; "
    ), f.getMessage()
    assert ("picked_decoy", "antisense") not in floors  # no target win: no floor
