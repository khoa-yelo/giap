"""TE7 of ``validation/pilot/te_real.py``: ``classify_tree``, ``diann_layer`` and the ``te7`` CLI
(M10-06 acceptance criteria 6, 7, 8 (pair half) and 10 (``te7`` half)).

Spec: docs/specs/M10-06-te-real-and-te7.md ("TE7: a second A2 run, byte for byte", examples 7-10,
edge cases "TE7 reference resumed or on another commit", "TE7 with a module switched off", ".gz
gate files", "DIA-NN replicate fails or lands on another class"). Trees, ``outputs.sha256`` lists
and converted tables are hand-made; the code-equivalence check runs on a temporary git
repository.

Interface fixed here (test-author reading of the spec "Interface"):

* ``classify_tree(files_a, files_b, known_columns, producers, labels=("ref", "rep1"))``:
  ``files_*`` map an OUT-relative POSIX path to a zero-argument callable returning the bytes;
  ``known_columns`` maps a table file name (``regions_stats.tsv``) to its LAPACK-dependent
  column names; ``producers`` maps a path to ``{"producer": step, "cpu_classes": {label: [class,
  ...]}}`` (a path absent from it has ``producer = "unknown"``); ``labels`` names the two runs
  (the spec's signature has no labels, example 8 needs ``only_in:rep1``). Returns te7.json's
  ``giap`` object: ``status``, ``n_gate``, ``n_external``, ``n_figures``, ``n_excluded``,
  ``differing`` (dicts with ``path``, ``class``, ``sub_class``, ``columns``, ``n_cells``,
  ``max_abs``, ``max_rel``, ``producer``, ``cpu_classes``, in sorted path order),
  ``external_differing`` and ``figures_differing`` (paths, or dicts with ``path``). Classes:
  ``only_in:<label>``, ``gzip_header_only``, ``row_order_only``, ``numeric_only`` (``sub_class``
  ``known_lapack`` or None), ``content``. ``n_*`` count distinct paths over both runs;
* ``diann_layer(sums_a, sums_b, conv_a, conv_b, *, conv_sha)``: ``sums_*`` map an output file
  name to its sha256 (a parsed ``outputs.sha256``); ``conv_*`` are the converted tables (frames,
  read only when L2 differs; None allowed then); ``conv_sha`` = (sha256 of A's ``a2.tsv.gz``,
  of B's). The spec's signature has no place for L2's two hashes, hence the keyword. Returns
  ``status`` (``identical_report`` / ``identical_converted`` / ``nondeterministic``), ``l0``
  {``files_equal``, ``files_differ``: sorted names}, ``l3`` (None, or {``only_a``, ``only_b``,
  ``both``, ``jaccard_precursors``, ``jaccard_peptides``, ``q_value_differs``,
  ``acceptance_flips``, ``intensity_differs``, ``max_abs_log2_intensity_ratio``});
* ``te_real.py te7`` writes ``te7.json`` (schema ``giap-pilot-te7/1``) and ``te7_files.tsv``
  (``path``, ``class``, ``sha256_a``, ``sha256_b``, ``equal``).
"""

from __future__ import annotations

import csv
import json
import math
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.te_real_m10 import (
    GENES_TSV,
    KNOWN,
    KNOWN_COLUMNS,
    REPO,
    base_tree,
    conv_frame,
    edit_manifest,
    ex7_replicate,
    ex8_replicate,
    load_te_real,
    make_repo,
    producers_for,
    readers,
    run_te_real,
    sha256_bytes,
    sums_text,
    write_conv,
    write_run,
)
from hypothesis import given, settings
from hypothesis import strategies as st

LABELS = ("ref", "rep1")


@pytest.fixture(scope="module")
def tr() -> Any:
    return load_te_real()


def classify(tr: Any, a: dict[str, bytes], b: dict[str, bytes], producers: Any = None) -> Any:
    return tr.classify_tree(readers(a), readers(b), KNOWN, producers or {}, labels=LABELS)


def _paths(xs: Any) -> list[str]:
    return [x if isinstance(x, str) else x["path"] for x in xs]


def _by_path(giap: Any) -> dict[str, Any]:
    return {d["path"]: d for d in giap["differing"]}


# ------------------------------------------------------------------ criterion 6: example 7


def test_example_7_pass(tr: Any) -> None:
    g = classify(tr, base_tree(), ex7_replicate())
    assert g["status"] == "pass"
    assert g["differing"] == []
    assert _paths(g["external_differing"]) == ["annotation/tools/promotech/promotech.csv"]
    assert _paths(g["figures_differing"]) == []
    # gate: tables x4, mapping x1, annotation/genes.gff3, plots/data x1, schemas/datapackage.json
    assert g["n_gate"] == 8
    assert g["n_external"] == 1 and g["n_figures"] == 1
    # work/**, logs/**, .nextflow*, run_manifest.json, report/**, preflight/**
    assert g["n_excluded"] == 6


def test_figures_are_reported_not_judged(tr: Any) -> None:
    rep = base_tree()
    rep["plots/landscape.svg"] = b"<svg>other</svg>\n"
    g = classify(tr, base_tree(), rep)
    assert g["status"] == "pass"
    assert _paths(g["figures_differing"]) == ["plots/landscape.svg"]


def test_plots_data_is_gate(tr: Any) -> None:
    rep = base_tree()
    rep["plots/data/landscape.tsv"] = b"x\ty\n1\t3\n"
    g = classify(tr, base_tree(), rep)
    assert g["status"] == "fail"
    assert list(_by_path(g)) == ["plots/data/landscape.tsv"]


def test_schemas_datapackage_is_gate(tr: Any) -> None:
    a, b = base_tree(), base_tree()
    a["schemas/datapackage.json"] = b"{}\n"
    g = classify(tr, a, b)
    assert g["status"] == "fail" and list(_by_path(g)) == ["schemas/datapackage.json"]


# ------------------------------------------------------------------ criterion 6: example 8


def test_example_8_classes(tr: Any) -> None:
    prod = producers_for(
        "tables/regions_stats.tsv", "context_stats", ["CPU_SKU:6148"], ["CPU_SKU:9754"]
    )
    g = classify(tr, base_tree(), ex8_replicate(), prod)
    assert g["status"] == "fail"
    d = _by_path(g)
    assert list(d) == sorted(d)  # sorted path order
    assert set(d) == {
        "tables/genes.tsv",
        "tables/pileup.tsv.gz",
        "tables/islands.tsv",
        "tables/regions_stats.tsv",
        "mapping/extra.tsv",
    }
    genes = d["tables/genes.tsv"]
    assert genes["class"] == "numeric_only" and genes["sub_class"] is None
    assert genes["columns"] == ["abundance_log2"] and genes["n_cells"] == 3
    expected_abs = max(
        abs(10.000000000000002 - 10.0),
        abs(11.500000000000002 - 11.5),
        abs(12.250000000000002 - 12.25),
    )
    assert genes["max_abs"] == pytest.approx(expected_abs, rel=1e-9) and genes["max_abs"] <= 2e-15
    assert 0 < genes["max_rel"] < 1e-15 and math.isfinite(genes["max_rel"])
    assert d["tables/pileup.tsv.gz"]["class"] == "gzip_header_only"
    assert d["tables/islands.tsv"]["class"] == "row_order_only"
    rs = d["tables/regions_stats.tsv"]
    assert rs["class"] == "numeric_only" and rs["sub_class"] == "known_lapack"
    assert sorted(rs["columns"]) == ["effect", "p_value"] and rs["n_cells"] == 2
    assert rs["producer"] == "context_stats"
    assert d["mapping/extra.tsv"]["class"] == "only_in:rep1"
    assert genes["producer"] == "unknown"


def test_known_lapack_needs_different_cpu_classes(tr: Any) -> None:
    rep = base_tree()
    rep["tables/regions_stats.tsv"] = ex8_replicate()["tables/regions_stats.tsv"]
    same = producers_for(
        "tables/regions_stats.tsv", "context_stats", ["CPU_SKU:9754"], ["CPU_SKU:9754"]
    )
    rs = _by_path(classify(tr, base_tree(), rep, same))["tables/regions_stats.tsv"]
    assert rs["class"] == "numeric_only" and rs["sub_class"] is None
    rs = _by_path(classify(tr, base_tree(), rep, {}))["tables/regions_stats.tsv"]
    assert rs["sub_class"] is None and rs["producer"] == "unknown"


def test_known_lapack_needs_listed_columns(tr: Any) -> None:
    a, b = base_tree(), base_tree()
    b["tables/regions_stats.tsv"] = a["tables/regions_stats.tsv"].replace(b"\t12\t", b"\t12.0\t")
    prod = producers_for(
        "tables/regions_stats.tsv", "context_stats", ["CPU_SKU:6148"], ["CPU_SKU:9754"]
    )
    rs = _by_path(classify(tr, a, b, prod))["tables/regions_stats.tsv"]
    assert rs["class"] == "numeric_only" and rs["columns"] == ["n_genes"]
    assert rs["sub_class"] is None  # n_genes is not a listed LAPACK column


def test_known_lapack_only_is_still_fail(tr: Any) -> None:
    rep = base_tree()
    rep["tables/regions_stats.tsv"] = ex8_replicate()["tables/regions_stats.tsv"]
    prod = producers_for(
        "tables/regions_stats.tsv", "context_stats", ["CPU_SKU:6148"], ["CPU_SKU:9754"]
    )
    g = classify(tr, base_tree(), rep, prod)
    assert [d["sub_class"] for d in g["differing"]] == ["known_lapack"]
    assert g["status"] == "fail"


def test_content_class_names_no_content(tr: Any) -> None:
    rep = base_tree()
    rep["tables/genes.tsv"] = GENES_TSV.replace(b"thrB", b"thrX")
    g = classify(tr, base_tree(), rep)
    d = _by_path(g)["tables/genes.tsv"]
    assert d["class"] == "content"
    text = json.dumps(g)
    assert "thrB" not in text and "thrX" not in text


def test_header_change_is_content(tr: Any) -> None:
    rep = base_tree()
    rep["tables/genes.tsv"] = GENES_TSV.replace(b"abundance_log2", b"abundance_log10")
    assert _by_path(classify(tr, base_tree(), rep))["tables/genes.tsv"]["class"] == "content"


def test_extra_row_is_content_not_numeric(tr: Any) -> None:
    rep = base_tree()
    rep["tables/genes.tsv"] = GENES_TSV + b"T_00050\tthrD\t1.0\t1\n"
    assert _by_path(classify(tr, base_tree(), rep))["tables/genes.tsv"]["class"] == "content"


def test_only_in_ref(tr: Any) -> None:
    a = base_tree()
    a["tables/extra_ref.tsv"] = b"a\n"
    g = classify(tr, a, base_tree())
    assert _by_path(g)["tables/extra_ref.tsv"]["class"] == "only_in:ref"


def test_gzip_with_different_payload_is_not_header_only(tr: Any) -> None:
    from helpers.te_real_m10 import gz

    rep = base_tree()
    rep["tables/pileup.tsv.gz"] = gz(b"contig\tpos\tdepth\nchr1\t1\t3\nchr1\t2\t5\n", 0)
    c = _by_path(classify(tr, base_tree(), rep))["tables/pileup.tsv.gz"]["class"]
    assert c != "gzip_header_only"


@settings(max_examples=40, deadline=None)
@given(
    files=st.dictionaries(
        st.sampled_from(
            [
                "tables/a.tsv",
                "tables/b.tsv",
                "mapping/c.tsv",
                "annotation/d.gff3",
                "plots/data/e.tsv",
                "logs/x.log",
                "work/y",
                "plots/f.svg",
                "annotation/tools/bakta/g.json",
            ]
        ),
        st.binary(min_size=1, max_size=40),
        min_size=1,
    ),
    data=st.data(),
)
def test_property_identical_trees_pass_one_gate_byte_fails(
    files: dict[str, bytes], data: Any
) -> None:
    tr = load_te_real()
    g = tr.classify_tree(readers(files), readers(dict(files)), KNOWN, {}, labels=LABELS)
    assert g["status"] == "pass" and g["differing"] == []
    gate = sorted(
        p
        for p in files
        if p.startswith(("tables/", "mapping/", "plots/data/")) or p == "annotation/d.gff3"
    )
    if gate:
        p = data.draw(st.sampled_from(gate))
        b = dict(files)
        b[p] = files[p] + b"\x00"
        g = tr.classify_tree(readers(files), readers(b), KNOWN, {}, labels=LABELS)
        assert g["status"] == "fail" and [d["path"] for d in g["differing"]] == [p]


# ------------------------------------------------------------------ criterion 7: DIA-NN layer

H = {k: sha256_bytes(k.encode()) for k in ("rep", "rep2", "lib", "log1", "log2", "conv", "conv2")}


def test_example_9_identical_after_conversion(tr: Any) -> None:
    sa = {"A2_report.parquet": H["rep"], "A2_lib.parquet": H["lib"], "A2.log.txt": H["log1"]}
    sb = {"A2_report.parquet": H["rep2"], "A2_lib.parquet": H["lib"], "A2.log.txt": H["log2"]}
    r = tr.diann_layer(sa, sb, None, None, conv_sha=(H["conv"], H["conv"]))
    assert r["status"] == "identical_converted"
    assert r["l3"] is None
    assert r["l0"]["files_equal"] == ["A2_lib.parquet"]
    assert r["l0"]["files_differ"] == ["A2.log.txt", "A2_report.parquet"]


def test_identical_report(tr: Any) -> None:
    s = {"A2_report.parquet": H["rep"]}
    r = tr.diann_layer(s, dict(s), None, None, conv_sha=(H["conv"], H["conv"]))
    assert r["status"] == "identical_report" and r["l3"] is None


def _ex10_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Example 10: 12 rows only in A, 9 only in B, 3 acceptance flips, max |log2 ratio| 0.004."""
    common = [("WT_1", f"PEPTIDE{i:02d}K", 2, False, 0.001, 1000.0) for i in range(30)]
    a = list(common)
    b = list(common)
    for i in range(3):  # acceptance flips
        a[i] = (a[i][0], a[i][1], 2, False, 0.009, 1000.0)
        b[i] = (b[i][0], b[i][1], 2, False, 0.011, 1000.0)
    a[5] = ("WT_1", a[5][1], 2, False, 0.001, 1000.0)
    b[5] = ("WT_1", b[5][1], 2, False, 0.001, 1000.0 * 2**0.004)
    b[6] = ("WT_1", b[6][1], 2, False, 0.002, 1000.0)  # q differs, both accepted
    a += [("WT_2", f"ONLYA{i:02d}K", 3, False, 0.001, 50.0) for i in range(11)]
    a += [("WT_2", "DECOYAK", 2, True, 0.001, 5.0)]
    b += [("WT_2", f"ONLYB{i:02d}K", 3, False, 0.001, 50.0) for i in range(9)]
    return conv_frame(a), conv_frame(b)


def _oracle_l3(a: pd.DataFrame, b: pd.DataFrame) -> dict[str, Any]:
    key = ["run", "modified_peptide", "charge", "is_decoy"]
    ka = {tuple(r[k] for k in key): r for _, r in a.iterrows()}
    kb = {tuple(r[k] for k in key): r for _, r in b.iterrows()}
    both = ka.keys() & kb.keys()

    def acc(r: Any) -> bool:
        return (not r["is_decoy"]) and r["q_value"] <= 0.01

    pa = {k for k, r in ka.items() if acc(r)}
    pb = {k for k, r in kb.items() if acc(r)}
    peps_a = {r["peptide"] for r in ka.values() if acc(r)}
    peps_b = {r["peptide"] for r in kb.values() if acc(r)}
    ratios = [abs(math.log2(kb[k]["intensity"] / ka[k]["intensity"])) for k in both]
    return {
        "only_a": len(ka.keys() - kb.keys()),
        "only_b": len(kb.keys() - ka.keys()),
        "both": len(both),
        "jaccard_precursors": len(pa & pb) / len(pa | pb),
        "jaccard_peptides": len(peps_a & peps_b) / len(peps_a | peps_b),
        "q_value_differs": sum(ka[k]["q_value"] != kb[k]["q_value"] for k in both),
        "acceptance_flips": sum(acc(ka[k]) != acc(kb[k]) for k in both),
        "intensity_differs": sum(ka[k]["intensity"] != kb[k]["intensity"] for k in both),
        "max_abs_log2_intensity_ratio": max(ratios),
    }


def test_example_10_nondeterministic(tr: Any) -> None:
    a, b = _ex10_tables()
    s = {"A2_report.parquet": H["rep"]}
    sb = {"A2_report.parquet": H["rep2"]}
    r = tr.diann_layer(s, sb, a, b, conv_sha=(H["conv"], H["conv2"]))
    assert r["status"] == "nondeterministic"
    l3 = r["l3"]
    want = _oracle_l3(a, b)
    assert (want["only_a"], want["only_b"], want["acceptance_flips"]) == (12, 9, 3)
    for k, v in want.items():
        assert l3[k] == pytest.approx(v, abs=1e-12), k
    assert l3["max_abs_log2_intensity_ratio"] == pytest.approx(0.004, abs=1e-9)


# ------------------------------------------------------------------ te7 CLI fixtures


@pytest.fixture(scope="module")
def repo(tmp_path_factory: pytest.TempPathFactory) -> tuple[Path, dict[str, str]]:
    root = tmp_path_factory.mktemp("m10_06_repo") / "repo"
    return root, make_repo(root)


def _pair(
    tmp: Path,
    repo: tuple[Path, dict[str, str]],
    rep_files: dict[str, bytes] | None = None,
    rep_commit: str = "base",
    **rep_over: Any,
) -> tuple[Path, Path]:
    _, shas = repo
    ref = write_run(
        tmp / "out_a2", base_tree(), commit=shas["base"], outdir_label="out_a2", hosts=("rye-02",)
    )
    rep = write_run(
        tmp / "out_a2_rep1",
        rep_files if rep_files is not None else ex7_replicate(),
        commit=shas[rep_commit],
        outdir_label="out_a2_rep1",
        hosts=("wheat-03",),
        **rep_over,
    )
    return ref, rep


def _te7(ref: Path, rep: Path, repo: tuple[Path, dict[str, str]], out: Path, *extra: str) -> Any:
    return run_te_real(
        "te7",
        "--ref",
        str(ref),
        "--rep",
        str(rep),
        "--repo",
        str(repo[0]),
        *extra,
        "--out",
        str(out),
    )


def _ok(res: Any) -> None:
    assert res.returncode == 0, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"


def _nothing_written(out: Path) -> bool:
    return not out.exists() or not any(out.iterdir())


def _doc(out: Path) -> dict[str, Any]:
    doc: dict[str, Any] = json.loads((out / "te7.json").read_text(encoding="utf-8"))
    return doc


def test_te7_cli_pass_and_layout(tmp_path: Path, repo: Any) -> None:
    ref, rep = _pair(tmp_path, repo)
    out = tmp_path / "metrics"
    _ok(_te7(ref, rep, repo, out))
    raw = (out / "te7.json").read_text(encoding="utf-8")
    doc = json.loads(raw)
    assert raw == json.dumps(doc, sort_keys=True, indent=2) + "\n"
    assert doc["schema"] == "giap-pilot-te7/1"
    assert doc["giap"]["status"] == "pass"
    assert _paths(doc["giap"]["external_differing"]) == ["annotation/tools/promotech/promotech.csv"]
    assert set(doc) >= {"pair", "giap", "diann", "end_to_end"}
    with (out / "te7_files.tsv").open(newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    assert rows and set(rows[0]) == {"path", "class", "sha256_a", "sha256_b", "equal"}
    assert [r["path"] for r in rows] == sorted(r["path"] for r in rows)
    # the DIA-NN layer was not given: never part of the verdict, no end-to-end claim
    assert doc["end_to_end"] is False


def test_te7_cli_fail_is_exit_0_numeric_only(tmp_path: Path, repo: Any) -> None:
    rep_files = ex7_replicate()
    rep_files["tables/genes.tsv"] = ex8_replicate()["tables/genes.tsv"]
    ref, rep = _pair(tmp_path, repo, rep_files)
    out = tmp_path / "metrics"
    _ok(_te7(ref, rep, repo, out))
    g = _doc(out)["giap"]
    assert g["status"] == "fail"
    d = _by_path(g)["tables/genes.tsv"]
    assert d["class"] == "numeric_only" and d["columns"] == ["abundance_log2"]


def test_te7_cli_docs_only_change_keeps_the_pair_eligible(tmp_path: Path, repo: Any) -> None:
    ref, rep = _pair(tmp_path, repo, rep_commit="docs")
    _ok(_te7(ref, rep, repo, tmp_path / "metrics"))


def _manifest_case(name: str) -> Any:
    def inputs(m: dict[str, Any]) -> None:
        m["inputs"][0]["sha256"] = "f" * 64

    def seeds(m: dict[str, Any]) -> None:
        m["seeds"]["master"] = m["seeds"]["master"] + 1

    def databases(m: dict[str, Any]) -> None:
        m["databases"] = [
            {
                "name": "contaminants",
                "version": "v2",
                "checksum": "sha256:" + "e" * 64,
                "retrieved": "2026-10-01",
                "size_bytes": 10,
            }
        ]

    def config(m: dict[str, Any]) -> None:
        # a module switched off in one run only (edge case)
        m["config"]["modules"] = [x for x in m["config"]["modules"] if x != "function"]

    def resume(m: dict[str, Any]) -> None:
        m["resume"] = True

    def cached(m: dict[str, Any]) -> None:
        m["steps"][0]["status"] = "cached"

    def dirty(m: dict[str, Any]) -> None:
        m["git_dirty"] = True

    def failed(m: dict[str, Any]) -> None:
        m["exit_code"] = 1
        m["status"] = "failed"

    return {
        "inputs": inputs,
        "seeds": seeds,
        "databases": databases,
        "config": config,
        "resume": resume,
        "cached": cached,
        "dirty": dirty,
        "failed": failed,
    }[name]


@pytest.mark.parametrize(
    ("which", "case"),
    [
        ("rep", "inputs"),
        ("rep", "seeds"),
        ("rep", "databases"),
        ("rep", "config"),
        ("ref", "resume"),
        ("ref", "cached"),
        ("ref", "dirty"),
        ("rep", "failed"),
        ("ref", "failed"),
    ],
)
def test_te7_pair_preconditions_exit_2(tmp_path: Path, repo: Any, which: str, case: str) -> None:
    ref, rep = _pair(tmp_path, repo)
    edit_manifest(ref if which == "ref" else rep, _manifest_case(case))
    out = tmp_path / "metrics"
    res = _te7(ref, rep, repo, out)
    assert res.returncode == 2, f"{case}: {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert _nothing_written(out)


def test_te7_non_code_equivalent_commits_exit_2(tmp_path: Path, repo: Any) -> None:
    ref, rep = _pair(tmp_path, repo, rep_commit="src")
    out = tmp_path / "metrics"
    res = _te7(ref, rep, repo, out)
    assert res.returncode == 2, res.stdout + res.stderr
    assert _nothing_written(out)


def test_te7_outdir_and_command_line_differences_are_allowed(tmp_path: Path, repo: Any) -> None:
    """The manifests differ in outdir, command line and timestamps (write_run): still a pair."""
    ref, rep = _pair(tmp_path, repo)
    edit_manifest(
        rep,
        lambda m: m.update(
            started_at="2026-10-24T01:00:00+00:00", finished_at="2026-10-24T01:30:00+00:00"
        ),
    )
    _ok(_te7(ref, rep, repo, tmp_path / "metrics"))


def test_te7_out_inside_repository_exit_2(tmp_path: Path, repo: Any) -> None:
    ref, rep = _pair(tmp_path, repo)
    out = REPO / ".m10_06_te7_out_test"
    try:
        res = _te7(ref, rep, repo, out)
        assert res.returncode == 2, res.stdout + res.stderr
        assert _nothing_written(out)
    finally:
        shutil.rmtree(out, ignore_errors=True)


def test_te7_different_hosts_no_difference_pass_hosts_recorded(tmp_path: Path, repo: Any) -> None:
    ref, rep = _pair(tmp_path, repo, base_tree())
    out = tmp_path / "metrics"
    _ok(_te7(ref, rep, repo, out))
    doc = _doc(out)
    assert doc["giap"]["status"] == "pass"
    text = json.dumps(doc)
    assert "rye" in text and "wheat" in text  # CPU classes / hosts recorded


# ------------------------------------------------------------------ DIA-NN layer through the CLI


def _diann_dirs(tmp: Path, report_same: bool, conv_same: bool) -> list[str]:
    da, db = tmp / "diann_a" / "A2", tmp / "diann_b" / "A2"
    for d, rep, node in (
        (da, H["rep"], "rye-02"),
        (db, H["rep"] if report_same else H["rep2"], "rye-05"),
    ):
        d.mkdir(parents=True)
        (d / "outputs.sha256").write_text(
            sums_text({"A2_report.parquet": rep, "A2.log.txt": sha256_bytes(node.encode())}),
            encoding="utf-8",
        )
        (d / "provenance.json").write_text(json.dumps({"node": node}) + "\n", encoding="utf-8")
    a, b = _ex10_tables()
    pa = write_conv(tmp / "conv_a" / "a2.tsv.gz", a)
    pb = write_conv(tmp / "conv_b" / "a2.tsv.gz", a if conv_same else b)
    return [
        "--diann-a",
        str(da),
        "--diann-b",
        str(db),
        "--peptides-a",
        str(pa),
        "--peptides-b",
        str(pb),
    ]


def test_te7_cli_example_9_end_to_end(tmp_path: Path, repo: Any) -> None:
    ref, rep = _pair(tmp_path, repo)
    out = tmp_path / "metrics"
    _ok(_te7(ref, rep, repo, out, *_diann_dirs(tmp_path, report_same=False, conv_same=True)))
    doc = _doc(out)
    assert doc["diann"]["status"] == "identical_converted"
    assert doc["diann"]["l3"] is None
    assert doc["end_to_end"] is True
    nodes = json.dumps(doc["diann"]["nodes"])
    assert "rye-02" in nodes and "rye-05" in nodes


def test_te7_cli_example_10_nondeterministic_is_not_the_verdict(tmp_path: Path, repo: Any) -> None:
    ref, rep = _pair(tmp_path, repo)
    out = tmp_path / "metrics"
    _ok(_te7(ref, rep, repo, out, *_diann_dirs(tmp_path, report_same=False, conv_same=False)))
    doc = _doc(out)
    assert doc["diann"]["status"] == "nondeterministic"
    assert doc["diann"]["l3"]["only_a"] == 12 and doc["diann"]["l3"]["only_b"] == 9
    assert doc["giap"]["status"] == "pass"  # never folded into the verdict
    assert doc["end_to_end"] is False
    assert "PEPTIDE" not in json.dumps(doc)  # aggregates only


def test_te7_end_to_end_false_when_giap_fails(tmp_path: Path, repo: Any) -> None:
    rep_files = ex7_replicate()
    rep_files["tables/islands.tsv"] = ex8_replicate()["tables/islands.tsv"]
    ref, rep = _pair(tmp_path, repo, rep_files)
    out = tmp_path / "metrics"
    _ok(_te7(ref, rep, repo, out, *_diann_dirs(tmp_path, report_same=True, conv_same=True)))
    doc = _doc(out)
    assert doc["diann"]["status"] == "identical_report"
    assert doc["giap"]["status"] == "fail" and doc["end_to_end"] is False


# ------------------------------------------------------------------ criterion 10 (te7 half)


def test_te7_twice_is_byte_identical(tmp_path: Path, repo: Any) -> None:
    rep_files = ex8_replicate()
    ref, rep = _pair(tmp_path, repo, rep_files)
    o1, o2 = tmp_path / "m1", tmp_path / "m2"
    extra = _diann_dirs(tmp_path, report_same=False, conv_same=False)
    _ok(_te7(ref, rep, repo, o1, *extra))
    _ok(_te7(ref, rep, repo, o2, *extra))
    for name in ("te7.json", "te7_files.tsv"):
        assert (o1 / name).read_bytes() == (o2 / name).read_bytes(), name


# ------------------------------------------------------------------ the committed column list


def test_known_columns_file_lists_numeric_schema_columns() -> None:
    assert KNOWN_COLUMNS.is_file(), f"{KNOWN_COLUMNS.relative_to(REPO)} is not committed"
    with KNOWN_COLUMNS.open(newline="", encoding="utf-8") as fh:
        rows = [
            r for r in csv.DictReader((ln for ln in fh if not ln.startswith("#")), delimiter="\t")
        ]
    assert rows and {"table", "column"} <= set(rows[0])
    tables = {Path(r["table"]).name for r in rows}
    assert "regions_stats.tsv" in tables
    for r in rows:
        stem = Path(r["table"]).name.removesuffix(".tsv")
        schema = REPO / "src" / "giap" / "schemas" / "tables" / f"{stem}.schema.json"
        assert schema.is_file(), r
        types = {
            f["name"]: f.get("type")
            for f in json.loads(schema.read_text(encoding="utf-8"))["fields"]
        }
        assert r["column"] in types, r
        assert types[r["column"]] in ("number", "integer"), r


def test_numeric_cells_numpy_float_parsing_sanity() -> None:
    """Self-check of the example 8 fixture: every differing genes.tsv cell parses as a float."""
    a = [ln.split("\t") for ln in GENES_TSV.decode().splitlines()[1:]]
    b = [ln.split("\t") for ln in ex8_replicate()["tables/genes.tsv"].decode().splitlines()[1:]]
    diff = [
        (x, y) for ra, rb in zip(a, b, strict=True) for x, y in zip(ra, rb, strict=True) if x != y
    ]
    assert len(diff) == 3
    assert all(np.isfinite(float(x)) and np.isfinite(float(y)) for x, y in diff)
