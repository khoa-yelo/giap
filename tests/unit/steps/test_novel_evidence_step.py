"""The wired ``novel_orfs`` step writes the M7-09 evidence (issue M7-09 acceptance criterion 4;
docs/specs/M7-09-novel-evidence-and-ranking.md rules 2-10; ADR-0012): ``rbs_supported``,
``promoter_supported``, the terminator columns, ``called`` and ``rank``, with flags in §9.4
order. The M7-08 executor-equality and TE7 tests stay in tests/unit/steps/test_novel_steps.py
and tests/integration (unchanged).

The expectations are computed here from the run's own work files, independently of
``giap.novel.evidence``: the sim100k run (tests/helpers/novel_steps_fixture.py) has no detected
gene with an OSTIR score, so its published table exercises the "threshold unestimated" branch;
the second test reruns the step on a gene table with every scored non-partial gene detected.

Readings pinned by the test author: the step reads the terminator layer it is given (under the
local executor the precomputed ``terminators`` TSV); ``promoter_supported`` is NA when
``promoter_score`` is NA.
"""

from __future__ import annotations

import itertools
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers import novel_steps_fixture as nf
from helpers.novel_evidence_fixture import FLAG_ORDER, flags_of, is_na

RBS_CLASSES = ("nterm_extension", "intergenic")


@pytest.fixture(scope="module")
def run1(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, Path, Path]:
    root = tmp_path_factory.mktemp("m7_09_sim100k")
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


def _layer_path(sim: Any, cfg: Any, key: str) -> Path:
    p = Path(getattr(cfg.annotation.precomputed, key))
    return p if p.is_absolute() else Path(sim.config_yaml).parent / p


def _threshold(gt: pd.DataFrame, q: float) -> float | None:
    part = gt["partial"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    det = gt["detected"].astype("boolean").fillna(False).to_numpy(dtype=bool)
    tir = pd.to_numeric(gt["tir_ostir"], errors="coerce").to_numpy(dtype=float)
    vals = tir[~part & det & ~np.isnan(tir)]
    return None if len(vals) < 20 else float(np.quantile(vals, q))


def _terminator_oracle(orf: Any, layer: pd.DataFrame, n: int) -> tuple[float | None, int | None]:
    """ADR-0012 d1 on a circular replicon: same strand, 5' end 0..300 nt after the stop's last
    base; the best confidence, ties -> the nearer."""
    last = orf.stop + 2 if orf.strand == "+" else orf.stop
    best: tuple[float, int] | None = None
    for r in layer.to_dict("records"):
        if r["contig"] != orf.contig or r["strand"] != orf.strand:
            continue
        five = int(r["start"]) if orf.strand == "+" else int(r["end"]) - 1
        d = ((five - last) if orf.strand == "+" else (last - five)) % n
        if d <= 300:
            c = float(r["confidence"])
            if best is None or c > best[0] or (c == best[0] and d < best[1]):
                best = (c, d)
    return (None, None) if best is None else best


def _check(table: pd.DataFrame, work: Path, sim: Any, cfg: Any, thr: float | None) -> None:
    from giap.annotate.precomputed import load_precomputed
    from giap.steps.novel import read_candidates

    orfs = {o.orf_id: o for o in read_candidates(work / "novel_candidates")[0]}
    ann = nf.annotation_of(sim)
    reps = dict(ann.model.replicons)
    terms = load_precomputed("terminators", _layer_path(sim, cfg, "terminators"), reps)
    pyes = cfg.annotation.promoters.promotech_yes
    tyes = cfg.annotation.terminators.terminatornet_yes
    recs = table.to_dict("records")
    assert recs
    for r in recs:
        oid, cls, f = r["orf_id"], r["class"], flags_of(r)
        assert f == sorted(f, key=FLAG_ORDER.index), (oid, f)
        tir = r["ostir_tir"]
        if thr is None or is_na(tir):
            assert is_na(r["rbs_supported"]), oid
            rbs = False
        else:
            rbs = bool(float(tir) >= thr)
            assert r["rbs_supported"] == rbs, oid
        assert ("rbs_supported" in f) is rbs, oid
        assert ("rbs_threshold_unestimated" in f) is (thr is None and cls in RBS_CLASSES), oid
        ps = r["promoter_score"]
        if is_na(ps):
            assert is_na(r["promoter_supported"]), oid
        else:
            assert r["promoter_supported"] == (float(ps) >= pyes), oid
        assert ("promoter_supported" in f) is (not is_na(ps) and float(ps) >= pyes), oid
        o = orfs[oid]
        conf, dist = _terminator_oracle(o, terms, reps[o.contig].length)
        if conf is None:
            assert is_na(r["terminator_conf"]) and is_na(r["terminator_distance_nt"]), oid
            assert r["terminator_downstream"] == False  # noqa: E712
        else:
            assert r["terminator_conf"] == pytest.approx(conf), oid
            assert r["terminator_distance_nt"] == dist, oid
            assert r["terminator_downstream"] == (conf >= tyes), oid
        assert ("terminator_supported" in f) is (conf is not None and conf >= tyes), oid
        two = int(r["n_peptides"]) >= cfg.novel_orfs.min_peptides_without_decoys
        if cls in ("cds_outframe", "antisense"):
            want = int(r["n_peptides"]) >= 2 and two
        elif cls == "intergenic":
            want = (two and rbs and not is_na(r["start_longest"])
                    and int(r["length_aa_best"]) >= cfg.novel_orfs.min_codons)  # fmt: skip
        elif cls == "nterm_extension":
            want = two and rbs  # the M7-06 walk only yields starts upstream of the gene
        elif cls == "cterm_extension":
            want = False
        else:
            want = two
        assert r["called"] == want, (oid, cls)
        assert r["class_fdr"] == "unestimated"
    called = [r for r in recs if r["called"]]
    assert sorted(int(r["rank"]) for r in called) == list(range(1, len(called) + 1))
    assert all(is_na(r["rank"]) for r in recs if not r["called"])
    by_rank = sorted(called, key=lambda r: int(r["rank"]))
    for a, b in itertools.pairwise(by_rank):
        assert int(a["n_peptides"]) >= int(b["n_peptides"]), (a["orf_id"], b["orf_id"])


def test_published_table_threshold(run1: tuple[Any, Path, Path], cfg1: Any) -> None:
    from giap.io.novel_tables import read_novel_orfs

    sim, out, work = run1
    gt = pd.read_parquet(work / "gene_table")
    thr = _threshold(gt, cfg1.novel_orfs.rbs_min_quantile)
    assert thr is not None, "M4-24 decision 8: detected genes give an estimated RBS threshold"
    table = read_novel_orfs(out / "tables" / "novel_orfs.tsv")
    assert set(table["class"]) & set(RBS_CLASSES), "the run plants nterm/intergenic ORFs"
    _check(table, work, sim, cfg1, thr)
    nf.assert_valid_novel_orfs(out / "tables" / "novel_orfs.tsv")


def test_step_with_detected_genes_applies_the_rbs_threshold(
    run1: tuple[Any, Path, Path], cfg1: Any, tmp_path: Path
) -> None:
    from giap.io.novel_tables import read_novel_orfs

    sim, _, work = run1
    gt = pd.read_parquet(work / "gene_table")
    scored = gt["tir_ostir"].notna() & ~gt["partial"].astype(bool)
    gt["detected"] = scored.to_numpy()
    gt_path = tmp_path / "gene_table.parquet"
    gt.to_parquet(gt_path, index=False)
    thr = _threshold(gt, cfg1.novel_orfs.rbs_min_quantile)
    assert thr is not None
    prom = _layer_path(sim, cfg1, "promoters")
    ins = nf.step_inputs("novel_orfs", work, {
        "gene_table": gt_path, "rbs_novel": _layer_path(sim, cfg1, "rbs_novel"),
        "novel_promoter_windows": prom, "promoter_windows": prom,
        "terminators": _layer_path(sim, cfg1, "terminators"),
    })  # fmt: skip
    o = nf.run("novel_orfs", cfg1, ins, tmp_path / "o")
    table = read_novel_orfs(o["novel_orfs"])
    _check(table, work, sim, cfg1, thr)
    rbs_rows = table[table["class"].isin(RBS_CLASSES) & table["ostir_tir"].notna()]
    assert len(rbs_rows) and rbs_rows["rbs_supported"].notna().all()
    assert table["called"].any() and table["rank"].notna().any()
    nf.assert_valid_novel_orfs(o["novel_orfs"])


def test_gff3_carries_rank_and_called_of_the_table(run1: tuple[Any, Path, Path]) -> None:
    from giap.io.novel_tables import read_novel_orfs

    _, out, _ = run1
    table = read_novel_orfs(out / "tables" / "novel_orfs.tsv")
    want = {r["orf_id"]: (str(bool(r["called"])).lower(),
                          "NA" if is_na(r["rank"]) else str(int(r["rank"])))
            for r in table.to_dict("records")}  # fmt: skip
    got = {}
    for line in (out / "tables" / "novel_orfs.gff3").read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or "\tCDS\t" not in line:
            continue
        attrs = dict(kv.split("=", 1) for kv in line.split("\t")[8].split(";"))
        got[attrs["ID"]] = (attrs["called"], attrs.get("rank", "NA"))
    assert got == want
    assert any(r == "1" for _, r in got.values())
