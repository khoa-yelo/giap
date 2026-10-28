"""Review round 1 of #366 (M10-10b: metric 6 and hypothesis 4); docs/specs/M10-10-a3-run.md,
the code review (B1, S1-S3) and the statistician review (the ``missing:<table>`` point), with
the lead decisions of this round.

Everything is synthetic (``helpers.pilot_metrics_m10``, ``helpers.a3_metrics_m10_10b``);
nothing comes from ``$GIPPA_ROOT``.

Readings fixed here (test-author choices, reversible):

* **B1, where the placement count lives.** The line ``placement: N queries (M de novo)`` is a
  plain ``ctx.log.info`` of the placement step, so it reaches only that task's
  ``OUT/work/<hash>/.command.log``. ``metrics.py`` finds the task through ``OUT/logs/trace.tsv``
  (the last ``COMPLETED`` row whose step is ``placement``) and
  ``gippa.runner.nextflow._task_dir`` (the ``workdir`` column, else the abbreviated ``hash``
  under ``OUT/work``). No trace, no task directory or no line: the warning
  ``denovo_count_unchecked`` (M10-14 cleans ``work/``). A differing M: ``Inconsistent``, exit 1.
* **B1, what M is compared with (like with like).** S6 says "check |S|" on the premise "with
  the default 7-40 window every de novo query has a length class, so M = |S|". The log's M is
  ``is_dn.sum()`` over the placement step's queries: the de novo queries after
  ``select_queries`` (the ``exact.min_len..max_len`` window, decoys, q) and **before** the
  ``mapping.length_classes`` filter. These tests compare M with that same count, which equals
  |S| under the default config, so that a run whose length classes leave a gap in the window
  is not a false ``Inconsistent`` (``test_b1_m_is_compared_before_the_length_class_filter``).
* **S2.** A de novo query with a NaN score forms no bin (``score_bins`` gives -1; the step's
  ``denovo_precision`` drops it) and is not an ``Inconsistent``; it still counts in the log's M.
* **S3.** ``render --final`` refuses a ``missing:`` reason on an A2 or A3 table of metric 6
  (``m6.A2.*``, ``m6.A3.*``, ``m6.overlap_a2_a3.*``); an A1 ``missing:`` stays a plain value.
  Section 7 (``not_evaluated``) lists those reasons at depth three as well.
* **S1.** One test per surviving mutant of the code review.

The H4 verdict rule is not pinned here beyond what the existing tests pin.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.a3_metrics_m10_10b import (
    A3,
    ALT_TYPES,
    BUCKETS,
    CONFIDENCE,
    NOVEL_CLASSES,
    PLACEMENT_ROW,
    build_a3,
    dn_row,
    example9_novel,
    example10_variants,
    example_stability,
    fdr_row,
    frame,
    queries_frame,
    stab,
    synthetic_peptides,
    write_denovo_table,
)
from helpers.pilot_gate_m10_15 import (
    FIXTURE_TEMPLATE,
    NOTES_OK,
    final_metrics,
    render_final,
    runbook_text,
    write_json,
    write_text,
)
from helpers.pilot_metrics_m10 import (
    Pilot,
    World,
    build_pilot,
    load_metrics,
    pep_row,
    run_cli,
    sha256,
    sim_world,
    write_peptides,
)


@pytest.fixture(scope="module")
def M() -> Any:
    return load_metrics()


def codes(warnings: Any) -> list[str]:
    return [w if isinstance(w, str) else str(w.get("code")) for w in warnings or []]


def _h4(
    M: Any, tmp: Path, q: pd.DataFrame, d2: Any, fdr: list[dict[str, object]], **kw: Any
) -> Any:
    kw.setdefault("n_denovo_logged", len(q))
    return M.hypothesis4(q, d2, frame("mapping_fdr", fdr, tmp), **kw)


def _tied_bins(groups: list[tuple[float, int]]) -> tuple[list[str], list[float], list[int]]:
    """Peptides in groups of tied scores, so ``score_bins`` gives one bin per group with the
    group's size (unequal n per bin); returns peptides, scores and each group's bin label."""
    from gippa.mapping.denovo_precision import score_bins

    n = sum(k for _, k in groups)
    peps = synthetic_peptides(n)
    scores = [s for s, k in groups for _ in range(k)]
    bins = score_bins(np.asarray(scores, dtype=float))
    labels, i = [], 0
    for _, k in groups:
        assert len(set(bins[i : i + k].tolist())) == 1, "a group split across bins"
        labels.append(int(bins[i]))
        i += k
    assert len(set(labels)) == len(groups), labels
    return peps, scores, labels


# ================================================================== S2: NaN de novo scores


def test_s2_nan_scores_are_dropped_not_inconsistent(M: Any, tmp_path: Path) -> None:
    """Two unscored queries: no bin -1, no ``Inconsistent``, ten bins of the 40 scored ones."""
    peps = synthetic_peptides(42)
    scores = [0.8 + 0.2 * (i + 1) / 40 for i in range(40)] + [float("nan")] * 2
    q = queries_frame(peps, scores)
    fdr = [PLACEMENT_ROW] + [fdr_row(b, 0.6, 0.4, 0.8) for b in range(10)]
    out = _h4(M, tmp_path, q, set(peps[:5]) | {peps[40]}, fdr, n_denovo_logged=None)
    assert [b["score_bin"] for b in out["bins"]] == list(range(10))
    assert sum(b["n_b"] for b in out["bins"]) == 40
    assert sum(b["k_b"] for b in out["bins"]) == 5  # the unscored confirmed peptide is in no bin


def test_s2_unscored_queries_still_count_in_the_log_check(M: Any, tmp_path: Path) -> None:
    """The step's "(M de novo)" counts every de novo query, scored or not."""
    peps = synthetic_peptides(42)
    scores = [0.8 + 0.2 * (i + 1) / 40 for i in range(40)] + [float("nan")] * 2
    q = queries_frame(peps, scores)
    fdr = [PLACEMENT_ROW] + [fdr_row(b, 0.6, 0.4, 0.8) for b in range(10)]
    out = _h4(M, tmp_path, q, set(), fdr, n_denovo_logged=42)
    assert "denovo_count_unchecked" not in codes(out["warnings"])
    with pytest.raises(M.Inconsistent):
        _h4(M, tmp_path, q, set(), fdr, n_denovo_logged=41)


# ================================================================== S1: surviving mutants


def _read_peptides(path: Path) -> pd.DataFrame:
    from gippa.io.peptides import read_peptide_table

    return read_peptide_table(path)[0]


def test_s1_contaminant_rows_drops_decoys_and_q(M: Any, tmp_path: Path) -> None:
    """Accepted target rows only: a decoy contaminant row and a q > 0.01 contaminant row are not
    counted (mutant: the ``is_decoy`` filter dropped)."""
    from gippa.config import GippaConfig

    rows = [
        pep_row("GADEFGHNK", "P1"),
        pep_row("GADEFGHQK", "P2"),
        pep_row("WWWGGGWWK", "Cont_P99"),
        pep_row("WWWGGGWYK", "Cont_P98", decoy=True),
        pep_row("WWWGGGWHK", "Cont_P97", decoy=True),
        pep_row("WWWGGGWNK", "Cont_P96", q=0.05),
    ]
    P = _read_peptides(write_peptides(tmp_path / "a2.tsv.gz", rows))
    out = M.contaminant_rows(P, GippaConfig())
    assert (out["k"], out["n"]) == (1, 3)


def test_s1_mean_genome_precision_is_n_weighted(M: Any, tmp_path: Path) -> None:
    """Unequal n per bin (40 and 20): the mean is (40 * 0.9 + 20 * 0.5) / 60, not 0.7."""
    peps, scores, labels = _tied_bins([(0.85, 40), (0.95, 20)])
    g = {labels[0]: (0.9, 0.85, 0.95), labels[1]: (0.5, 0.4, 0.6)}
    fdr = [PLACEMENT_ROW] + [fdr_row(b, *g[b]) for b in labels]
    out = _h4(M, tmp_path, queries_frame(peps, scores), set(peps[:10]), fdr)
    assert out["B"] == 2
    mean = out["mean_genome_precision"]
    assert mean["value"] == pytest.approx((40 * 0.9 + 20 * 0.5) / 60)
    assert mean["n"] == 60


def test_s1_pooled_l_over_reportable_bins_only(M: Any, tmp_path: Path) -> None:
    """A non-reportable bin (n = 10 < 20) with k = 5 stays out of ``pooled_L``."""
    peps, scores, labels = _tied_bins([(0.85, 30), (0.9, 25), (0.95, 10)])
    fdr = [PLACEMENT_ROW] + [fdr_row(b, 0.6, 0.4, 0.8) for b in labels]
    d2 = set(peps[:15]) | set(peps[30:40]) | set(peps[55:60])
    out = _h4(M, tmp_path, queries_frame(peps, scores), d2, fdr)
    by = {b["score_bin"]: b for b in out["bins"]}
    assert [by[b]["n_b"] for b in labels] == [30, 25, 10]
    assert [by[b]["k_b"] for b in labels] == [15, 10, 5]
    assert by[labels[2]]["reportable"] is False and out["B"] == 2
    pooled = out["pooled_L"]
    assert (pooled["k"], pooled["n"]) == (25, 55)
    assert pooled["value"] == pytest.approx(25 / 55)


def test_s1_null_score_bin_row_is_inconsistent(M: Any, tmp_path: Path) -> None:
    """A ``denovo_precision`` row without a ``score_bin`` is a malformed ``mapping_fdr.tsv``."""
    peps = synthetic_peptides(40)
    q = queries_frame(peps, [0.8 + 0.2 * (i + 1) / 40 for i in range(40)])
    fdr = [PLACEMENT_ROW] + [fdr_row(b, 0.6, 0.4, 0.8) for b in range(10)]
    _h4(M, tmp_path, q, set(), fdr)  # the well-formed table passes
    bad = [*fdr, fdr_row(None, 0.6, 0.4, 0.8)]  # type: ignore[arg-type]
    with pytest.raises(M.Inconsistent):
        _h4(M, tmp_path, q, set(), bad)


def test_s1_l_b_searchable_counts_only_searchable_confirmed(M: Any, tmp_path: Path) -> None:
    """A confirmed peptide outside the searchable set does not count in the numerator."""
    peps = synthetic_peptides(30)
    q = queries_frame(peps, [0.9] * 30)
    d2 = set(peps[:6]) | {peps[20], peps[25]}  # two confirmed, not searchable
    out = _h4(M, tmp_path, q, d2, [fdr_row(0, 0.6, 0.4, 0.8)], searchable=set(peps[:10]))
    s = out["bins"][0]["L_b_searchable"]
    assert (s["k"], s["n"]) == (6, 10)
    assert out["bins"][0]["k_b"] == 8


def test_s1_compare_row_with_a3_listed_first(M: Any, tmp_path: Path) -> None:
    """compare may write the pair as (A3, A2): its counts are swapped back before the ALT
    cross-check (example 10: A2 1 site, A3 2 sites, 1 shared)."""
    n2, n3 = example9_novel()
    v2, v3 = example10_variants()
    rows = [
        r
        for r in example_stability()
        if not (
            {r["run_a"], r["run_b"]} == {"A2", "A3"} and r["family"] in ("alt_sites", "novel_orfs")
        )
    ]
    rows.append(stab("novel_orfs", "A3", "A2", 3, 2, 0, 0.0))
    rows.append(stab("alt_sites", "A3", "A2", 2, 1, 1, 0.5))
    tables = {
        "A1": {
            "novel_orfs": frame("novel_orfs", [], tmp_path),
            "variants": frame("variants", [], tmp_path),
        },
        "A2": {
            "novel_orfs": frame("novel_orfs", n2, tmp_path),
            "variants": frame("variants", v2, tmp_path),
        },
        "A3": {
            "novel_orfs": frame("novel_orfs", n3, tmp_path),
            "variants": frame("variants", v3, tmp_path),
        },
    }
    out = M.metric6(tables, frame("compare_stability", rows, tmp_path))
    alt = out["overlap_a2_a3"]["alt_sites"]["alt"]
    assert (alt["n_a2"], alt["n_a3"], alt["n_shared"]) == (1, 2, 1)
    assert "duplicate_orf_key" not in codes(out["warnings"])


# ------------------------------------------------------------------ CLI fixtures


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> World:
    return sim_world(tmp_path_factory.mktemp("m10_10b_review_sim"))


@pytest.fixture
def pilot(tmp_path: Path, world: World) -> Pilot:
    return build_pilot(tmp_path / "pilot", world)


def _compute(p: Pilot, out: Path, *extra: str) -> Any:
    return run_cli(*p.compute_args(out, *extra), cwd=p.root)


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _metrics(out: Path) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads((out / "metrics.json").read_text(encoding="utf-8"))
    return doc


def _nothing_written(out: Path) -> bool:
    return not out.exists() or not any(out.iterdir())


def _patch_json(path: Path, fn: Any) -> None:
    doc = json.loads(path.read_text(encoding="utf-8"))
    fn(doc)
    path.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def test_s1_gippa_commit_differs_warns_for_a3(pilot: Pilot, world: World, tmp_path: Path) -> None:
    a3 = build_a3(pilot, world)
    other = "fedcba9876543210fedcba9876543210fedcba98"
    _patch_json(a3.out_a3 / "run_manifest.json", lambda d: d.__setitem__("git_commit", other))
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out, *a3.args()))
    ws = [w for w in _metrics(out)["warnings"] if isinstance(w, dict)]
    hit = [w for w in ws if w.get("code") == "gippa_commit_differs" and "A3" in w]
    assert hit and hit[0]["A3"] == other, ws


def test_s1_missing_mapping_fdr_is_a_precondition(
    pilot: Pilot, world: World, tmp_path: Path
) -> None:
    a3 = build_a3(pilot, world)
    (a3.out_a3 / "mapping" / "mapping_fdr.tsv").unlink()
    out = tmp_path / "metrics"
    res = _compute(pilot, out, *a3.args())
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert res.stderr.startswith("precondition:") and "mapping_fdr" in res.stderr, res.stderr
    assert _nothing_written(out)


def test_s1_strain_db_without_a3_is_a_precondition(
    pilot: Pilot, world: World, tmp_path: Path
) -> None:
    a3 = build_a3(pilot, world)
    out = tmp_path / "metrics"
    res = _compute(pilot, out, "--strain-db", str(a3.strain_db))
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert res.stderr.startswith("precondition:") and "--strain-db" in res.stderr, res.stderr
    assert _nothing_written(out)


# ================================================================== B1: the placement count

PLACEMENT_TASK = "GIPPA:MAPPING:PLACEMENT (placement)"
PLACEMENT_HASH = "3f/9a1b2c"
PLACEMENT_DIR = "3f/9a1b2c4d5e6f708192a3b4c5d6e7f80"


def placement_line(n: int, m: int) -> str:
    return (
        "2026-10-24 12:00:00,000 INFO gippa.steps.mapping: "
        f"placement: {n + 1000} queries ({m} de novo), 900 with a target hit, 850 accepted\n"
    )


def write_task_logs(
    out: Path, m: int | None, *, workdir_column: bool = True, failed_m: int | None = None
) -> None:
    """``OUT/logs/trace.tsv`` (Nextflow's columns) with an earlier task, an optional FAILED
    placement attempt and the COMPLETED placement task; the completed task's
    ``.command.log`` holds the placement line with ``(m de novo)`` (no line when ``m`` is None).
    ``workdir_column=False`` leaves ``workdir`` as ``-`` so the directory comes from ``hash``."""
    from helpers.manifest_data import trace_row

    work = out / "work"
    rows = []
    exact_dir = work / "1a" / "2b3c4d5e6f"
    exact_dir.mkdir(parents=True)
    (exact_dir / ".command.log").write_text("map_exact: 30 queries, 12 hits\n", encoding="utf-8")
    rows.append(
        trace_row("GIPPA:MAPPING:MAP_EXACT (map_exact)", hash="1a/2b3c4d", workdir=str(exact_dir))
    )
    if failed_m is not None:
        fdir = work / "77" / "aa00bb11cc"
        fdir.mkdir(parents=True)
        (fdir / ".command.log").write_text(placement_line(0, failed_m), encoding="utf-8")
        rows.append(trace_row(PLACEMENT_TASK, status="FAILED", hash="77/aa00bb", workdir=str(fdir)))
    tdir = work / PLACEMENT_DIR
    tdir.mkdir(parents=True)
    text = "placement: loading the index\n" + ("" if m is None else placement_line(30, m))
    (tdir / ".command.log").write_text(text, encoding="utf-8")
    rows.append(
        trace_row(
            PLACEMENT_TASK,
            hash=PLACEMENT_HASH,
            workdir=str(tdir) if workdir_column else "-",
        )
    )
    path = out / "logs" / "trace.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0])
    lines = ["\t".join(cols)] + ["\t".join(str(r[c]) for c in cols) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _with_logs(p: Pilot, world: World, tmp: Path, m: int | None, **kw: Any) -> tuple[Any, Path, A3]:
    a3 = build_a3(p, world)
    write_task_logs(a3.out_a3, m, **kw)
    out = tmp / "metrics"
    return _compute(p, out, *a3.args()), out, a3


@pytest.mark.parametrize("workdir_column", [True, False], ids=["workdir", "hash"])
def test_b1_matching_count_is_checked(
    pilot: Pilot, world: World, tmp_path: Path, workdir_column: bool
) -> None:
    res, out, a3 = _with_logs(pilot, world, tmp_path, 30, workdir_column=workdir_column)
    _ok(res)
    assert len(a3.denovo) == 30
    doc = _metrics(out)
    assert "denovo_count_unchecked" not in codes(doc["h4"]["warnings"])
    assert "denovo_count_unchecked" not in codes(doc["warnings"])


def test_b1_last_completed_placement_task_is_read(
    pilot: Pilot, world: World, tmp_path: Path
) -> None:
    """An earlier FAILED attempt's log (a different count) is not the run's count."""
    res, out, _ = _with_logs(pilot, world, tmp_path, 30, failed_m=29)
    _ok(res)
    assert "denovo_count_unchecked" not in codes(_metrics(out)["h4"]["warnings"])


def test_b1_differing_count_exit_1(pilot: Pilot, world: World, tmp_path: Path) -> None:
    res, out, _ = _with_logs(pilot, world, tmp_path, 29)
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert res.stderr.startswith("inconsistent:"), res.stderr
    assert _nothing_written(out)


@pytest.mark.parametrize("case", ["no_trace", "no_task_dir", "no_line"])
def test_b1_no_log_warns_denovo_count_unchecked(
    pilot: Pilot, world: World, tmp_path: Path, case: str
) -> None:
    a3 = build_a3(pilot, world)
    if case != "no_trace":
        write_task_logs(a3.out_a3, None if case == "no_line" else 29)
    if case == "no_task_dir":
        import shutil

        shutil.rmtree(a3.out_a3 / "work")  # M10-14's cleanup removes work/
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out, *a3.args()))
    doc = _metrics(out)
    assert "denovo_count_unchecked" in codes(doc["h4"]["warnings"])
    assert "denovo_count_unchecked" in codes(doc["warnings"])


LONG_PEPTIDE = "GADEFGHNQSTVWYGADEFGHNQSTVWYGAK"  # 31 residues, no I/L, one K at the end


def _gap_in_length_classes(a3: A3) -> None:
    """A3's resolved config with no length class for 31 residues and one de novo peptide of
    that length added to ``--peptides-a3`` (shas updated): |S| stays 30, the step's M is 31."""
    assert all(len(p) != 31 for p in a3.denovo)
    with gzip.open(a3.pep_a3, "rt", encoding="utf-8") as fh:
        text = fh.read()
    lines = [ln for ln in text.splitlines() if ln]
    cols = lines[0].split("\t")
    rows = [dict(zip(cols, ln.split("\t"), strict=True)) for ln in lines[1:]]
    rows.append(dn_row(LONG_PEPTIDE, 0.85, scan=5000))
    write_denovo_table(a3.pep_a3, rows)
    new = sha256(a3.pep_a3)

    def manifest(d: dict[str, Any]) -> None:
        for i in d["inputs"]:
            if i.get("role") == "peptides":
                i["sha256"] = new
        d["config"]["mapping"]["length_classes"] = [
            [7, 7], [8, 8], [9, 9], [10, 11], [12, 14], [15, 19], [20, 30], [32, 40],
        ]  # fmt: skip

    def conv(d: dict[str, Any]) -> None:
        d["outputs"]["peptides"]["sha256"] = new

    _patch_json(a3.out_a3 / "run_manifest.json", manifest)
    _patch_json(a3.conv_a3, conv)


def test_b1_m_is_compared_before_the_length_class_filter(
    pilot: Pilot, world: World, tmp_path: Path
) -> None:
    """Like with like: M counts the de novo queries before the length-class filter, so a
    31-residue query without a class is in M (31) but not in |S| (30)."""
    a3 = build_a3(pilot, world)
    _gap_in_length_classes(a3)
    write_task_logs(a3.out_a3, 31)
    out = tmp_path / "metrics"
    _ok(_compute(pilot, out, *a3.args()))
    doc = _metrics(out)
    assert doc["h4"]["n_queries"] == 30 or sum(b["n_b"] for b in doc["h4"]["bins"]) == 30
    assert "denovo_count_unchecked" not in codes(doc["h4"]["warnings"])
    import shutil

    shutil.rmtree(a3.out_a3 / "work")
    (a3.out_a3 / "logs" / "trace.tsv").unlink()
    write_task_logs(a3.out_a3, 30)
    out2 = tmp_path / "metrics2"
    res = _compute(pilot, out2, *a3.args())
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


# ================================================================== S3: missing: under --final


def _summary_novel() -> dict[str, Any]:
    return {
        "n_rows": 0,
        "by_class": {c: dict.fromkeys(BUCKETS, 0) for c in NOVEL_CLASSES},
        "totals": dict.fromkeys(BUCKETS, 0),
    }


def _summary_variants() -> dict[str, Any]:
    return {
        "alt": {t: dict.fromkeys(CONFIDENCE, 0) for t in ALT_TYPES},
        "mixed": 0,
        "class_error": {"numeric": 0, "unestimated": 0},
    }


def _overlap_set() -> dict[str, Any]:
    empty = {
        "value": None,
        "ci_low": None,
        "ci_high": None,
        "k": 0,
        "n": 0,
        "reason": "empty_denominator",
    }
    return {
        "n_a2": 0,
        "n_a3": 0,
        "n_shared": 0,
        "union": 0,
        "jaccard": None,
        "reason": "both_empty",
        "a3_confirmed_by_a2": dict(empty),
        "a2_found_by_a3": dict(empty),
    }


def m6_doc(missing: tuple[str, str] | None) -> dict[str, Any]:
    """A computed ``m6`` block; ``missing = (arm, table)`` replaces one table's summary (or one
    overlap family) by ``{"reason": "missing:<file>"}``."""
    m6: dict[str, Any] = {
        a: {"novel_orfs": _summary_novel(), "variants": _summary_variants()}
        for a in ("A1", "A2", "A3")
    }
    m6["overlap_a2_a3"] = {
        "novel_orfs": {"called": _overlap_set(), "any_row": _overlap_set(), "orf_id_jaccard": None},
        "alt_sites": {"alt": _overlap_set(), "high": _overlap_set()},
    }
    m6["warnings"] = []
    if missing is not None:
        arm, table = missing
        file = "novel_orfs.tsv" if table == "novel_orfs" else "variants.tsv"
        m6[arm][table] = {"reason": f"missing:{file}"}
    return m6


M6_TEMPLATE = FIXTURE_TEMPLATE.replace(
    "## 5. TE2-TE5 and TE7",
    "### Metric 6\n\n<!-- m6 class rows -->\n\n<!-- m6 alt rows -->\n\n## 5. TE2-TE5 and TE7",
)


def _render(tmp: Path, m6: dict[str, Any], *, final: bool) -> tuple[Any, Path]:
    doc = final_metrics()
    doc["m6"] = m6
    m = write_json(tmp / "in" / "metrics.json", doc)
    n = write_text(tmp / "in" / "report_notes.md", NOTES_OK)
    t = write_text(tmp / "in" / "report_template.md", M6_TEMPLATE)
    r = write_text(tmp / "in" / "reproduce.md", runbook_text()) if final else None
    out = tmp / "out" / "REPORT.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    return render_final(m, out, notes=n, template=t, runbook=r, final=final), out


def test_s3_baseline_computed_m6_renders_final(tmp_path: Path) -> None:
    res, out = _render(tmp_path, m6_doc(None), final=True)
    _ok(res)
    assert "**missing:" not in out.read_text(encoding="utf-8")


@pytest.mark.parametrize("table", ["novel_orfs", "variants"])
def test_s3_a1_missing_table_stays_a_plain_value_under_final(tmp_path: Path, table: str) -> None:
    res, out = _render(tmp_path, m6_doc(("A1", table)), final=True)
    _ok(res)
    text = out.read_text(encoding="utf-8")
    assert "missing:" in text and "**missing:" not in text


A2_A3_MISSING = [
    ("A2", "novel_orfs"),
    ("A2", "variants"),
    ("A3", "novel_orfs"),
    ("A3", "variants"),
    ("overlap_a2_a3", "novel_orfs"),
    ("overlap_a2_a3", "alt_sites"),
]


@pytest.mark.parametrize("where", A2_A3_MISSING, ids=lambda w: ".".join(w))
def test_s3_final_refuses_a2_a3_missing_table(tmp_path: Path, where: tuple[str, str]) -> None:
    arm, table = where
    res, out = _render(tmp_path, m6_doc(where), final=True)
    assert res.returncode == 1, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert not out.exists(), "nothing is written on a refused render"
    msg = res.stdout + res.stderr
    assert (
        f"m6.{arm}" in msg or ("novel_orfs.tsv" if table == "novel_orfs" else "variants.tsv") in msg
    )
    draft, _ = _render(tmp_path / "draft", m6_doc(where), final=False)
    _ok(draft)  # the draft keeps rendering it


@pytest.mark.parametrize("where", A2_A3_MISSING, ids=lambda w: ".".join(w))
def test_s3_section_7_lists_depth_3_reasons(tmp_path: Path, where: tuple[str, str]) -> None:
    arm, table = where
    res, out = _render(tmp_path, m6_doc(where), final=False)
    _ok(res)
    sec7 = out.read_text(encoding="utf-8").split("## 7. Not evaluated", 1)[1].split("## 10.", 1)[0]
    file = "novel_orfs.tsv" if table == "novel_orfs" else "variants.tsv"
    assert f"`m6.{arm}.{table}`" in sec7 and f"missing:{file}" in sec7, sec7


def test_s3_not_evaluated_walks_depth_3(M: Any) -> None:
    lines = M.not_evaluated({"m6": m6_doc(("A3", "novel_orfs"))})
    assert any("`m6.A3.novel_orfs`" in ln and "missing:novel_orfs.tsv" in ln for ln in lines), lines
