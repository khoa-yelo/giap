"""M10-06b, TE7: ``validation/pilot/te7.py`` (``classify_tree``, ``diann_layer``) and the ``te7``
subcommand that ``te_real.py`` keeps and delegates (M10-06 acceptance criteria 6, 7, 8 (pair half)
and 10 (``te7`` half); review round 1 S1-S4, N3 and test gaps 1-3 and 5-7).

Spec: docs/specs/M10-06-te-real-and-te7.md ("TE7: a second A2 run, byte for byte", examples 7-10,
edge cases "TE7 reference resumed or on another commit", "TE7 with a module switched off", ".gz
gate files", "DIA-NN replicate fails or lands on another class"). Trees, ``outputs.sha256`` lists
and converted tables are hand-made; the code-equivalence check runs on a temporary git
repository.

Interface fixed here (test-author reading of the spec "Interface"; lead decisions of review
round 1: TE7 lives in ``validation/pilot/te7.py``, the CLI ``te_real.py te7`` is unchanged):

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
  (``path``, ``class``, ``sha256_a``, ``sha256_b``, ``equal``). Run labels: the ``rep\\d+`` of
  each run directory's name, else ``ref`` for the reference and ``rep1`` for the replicate (N3).
  The producing step of a gate file and its CPU classes come from ``logs/trace.tsv``'s ``name``
  (Nextflow aliases such as ``VARIANTS_STEP`` map to the registry step ``variants``; S3);
  ``known_lapack`` needs both runs' class lists non-empty and different (S1); ``--peptides-a``
  must have the sha256 of both manifests' ``peptides`` input, else exit 2 (S2).
"""

from __future__ import annotations

import csv
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.te7_m10 import (
    GENES_TSV,
    KNOWN,
    KNOWN_COLUMNS,
    REGIONS_TSV_LAPACK,
    REPO,
    base_tree,
    conv_frame,
    edit_manifest,
    ex7_replicate,
    ex8_replicate,
    gz,
    load_te7,
    make_repo,
    producers_for,
    readers,
    set_peptides_sha,
    sha256_bytes,
    sums_text,
    write_conv,
    write_named_trace,
    write_run,
)
from helpers.te_real_m10 import load_te_real, run_te_real
from hypothesis import given, settings
from hypothesis import strategies as st

LABELS = ("ref", "rep1")


@pytest.fixture(scope="module")
def tr() -> Any:
    return load_te7()


def test_te7_module_is_the_one_copy(tr: Any) -> None:
    """TE7 moved to ``te7.py`` (L14 split); ``te_real.py`` keeps no second copy of it."""
    te_real = load_te_real()
    for name in ("classify_tree", "diann_layer", "path_class", "pair_problems"):
        fn = getattr(tr, name)
        assert callable(fn), name
        assert getattr(te_real, name, fn) is fn, f"te_real.{name} is a second copy"


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
    rep = base_tree()
    rep["tables/pileup.tsv.gz"] = gz(b"contig\tpos\tdepth\nchr1\t1\t3\nchr1\t2\t5\n", 0)
    c = _by_path(classify(tr, base_tree(), rep))["tables/pileup.tsv.gz"]["class"]
    assert c != "gzip_header_only"


def test_gz_gate_file_that_is_not_gzip_is_content(tr: Any) -> None:
    """Test gap 7: a ``.gz`` path whose bytes are no gzip stream falls back to ``content``."""
    a, b = base_tree(), base_tree()
    a["tables/pileup.tsv.gz"] = b"contig\tpos\tdepth\nchr1\t1\t3\n"
    b["tables/pileup.tsv.gz"] = b"contig\tpos\tdepth\nchr1\t1\t4\n"
    g = classify(tr, a, b)
    assert g["status"] == "fail"
    assert _by_path(g)["tables/pileup.tsv.gz"]["class"] == "content"


def test_non_utf8_gate_file_is_content(tr: Any) -> None:
    """Test gap 7: bytes that do not decode as UTF-8 are ``content`` (no crash, no exit 1)."""
    a, b = base_tree(), base_tree()
    a["annotation/genes.gff3"] = b"##gff-version 3\n\xff\xfe\x00\x01\n"
    b["annotation/genes.gff3"] = b"##gff-version 3\n\xff\xfe\x00\x02\n"
    d = _by_path(classify(tr, a, b))["annotation/genes.gff3"]
    assert d["class"] == "content" and d["sub_class"] is None


@pytest.mark.parametrize(
    ("ref_classes", "rep_classes"),
    [(["rye"], []), ([], ["wheat"]), ([], [])],
    ids=["rep_has_none", "ref_has_none", "neither_has_any"],
)
def test_known_lapack_needs_cpu_classes_in_both_runs(
    tr: Any, ref_classes: list[str], rep_classes: list[str]
) -> None:
    """S1: a step absent from one run's trace (an empty class list) is no evidence of different
    CPU classes: listed columns only, still ``numeric_only`` without ``known_lapack``."""
    rep = base_tree()
    rep["tables/regions_stats.tsv"] = REGIONS_TSV_LAPACK
    prod = producers_for("tables/regions_stats.tsv", "context_stats", ref_classes, rep_classes)
    rs = _by_path(classify(tr, base_tree(), rep, prod))["tables/regions_stats.tsv"]
    assert rs["class"] == "numeric_only"
    assert rs["sub_class"] is None
    assert rs["producer"] == "context_stats"


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
    tr = load_te7()
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


@pytest.mark.parametrize(
    ("ia", "ib"),
    [(0.0, 1000.0), (1000.0, 0.0), (float("nan"), 1000.0), (1000.0, float("nan"))],
    ids=["zero_in_a", "zero_in_b", "na_in_a", "na_in_b"],
)
def test_l3_max_log2_ratio_skips_zero_and_na_intensities(tr: Any, ia: float, ib: float) -> None:
    """Surviving mutant 19: the maximum |log2 ratio| is over the pairs with both intensities
    finite and positive; a zero or NA intensity is counted in ``intensity_differs`` only."""
    rows = [("WT_1", f"PEPTIDE{i:02d}K", 2, False, 0.001, 1000.0) for i in range(4)]
    a, b = list(rows), list(rows)
    b[1] = ("WT_1", rows[1][1], 2, False, 0.001, 1000.0 * 2**0.25)  # the one finite ratio
    a[2] = ("WT_1", rows[2][1], 2, False, 0.001, ia)
    b[2] = ("WT_1", rows[2][1], 2, False, 0.001, ib)
    r = tr.diann_layer(
        {"A2_report.parquet": H["rep"]},
        {"A2_report.parquet": H["rep2"]},
        conv_frame(a),
        conv_frame(b),
        conv_sha=(H["conv"], H["conv2"]),
    )
    l3 = r["l3"]
    assert l3["max_abs_log2_intensity_ratio"] == pytest.approx(0.25, abs=1e-12)
    assert l3["intensity_differs"] == 2
    assert l3["both"] == 4 and l3["only_a"] == l3["only_b"] == 0


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
    *,
    ref_files: dict[str, bytes] | None = None,
    names: tuple[str, str] = ("out_a2", "out_a2_rep1"),
    **rep_over: Any,
) -> tuple[Path, Path]:
    _, shas = repo
    ref = write_run(
        tmp / names[0],
        ref_files if ref_files is not None else base_tree(),
        commit=shas["base"],
        outdir_label=names[0],
        hosts=("rye-02",),
    )
    rep = write_run(
        tmp / names[1],
        rep_files if rep_files is not None else ex7_replicate(),
        commit=shas[rep_commit],
        outdir_label=names[1],
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
    assert doc["diann"]["status"] == "not_run"  # test gap 5
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
        ("rep", "resume"),  # test gap 6: both runs must be fresh, not only the reference
        ("rep", "cached"),
        ("rep", "dirty"),
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


def _diann_dirs(
    tmp: Path,
    report_same: bool,
    conv_same: bool,
    pair: tuple[Path, Path] | None = None,
) -> list[str]:
    """The four DIA-NN flags. With ``pair``, both manifests' ``peptides`` input gets the sha256 of
    ``--peptides-a`` (the GIAP pair ran on the first search's converted table; S2)."""
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
    for out in pair or ():
        set_peptides_sha(out, sha256_bytes(pa.read_bytes()))
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
    _ok(_te7(ref, rep, repo, out, *_diann_dirs(tmp_path, False, True, (ref, rep))))
    doc = _doc(out)
    assert doc["diann"]["status"] == "identical_converted"
    assert doc["diann"]["l3"] is None
    assert doc["end_to_end"] is True
    nodes = json.dumps(doc["diann"]["nodes"])
    assert "rye-02" in nodes and "rye-05" in nodes


def test_te7_cli_example_10_nondeterministic_is_not_the_verdict(tmp_path: Path, repo: Any) -> None:
    ref, rep = _pair(tmp_path, repo)
    out = tmp_path / "metrics"
    _ok(_te7(ref, rep, repo, out, *_diann_dirs(tmp_path, False, False, (ref, rep))))
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
    _ok(_te7(ref, rep, repo, out, *_diann_dirs(tmp_path, True, True, (ref, rep))))
    doc = _doc(out)
    assert doc["diann"]["status"] == "identical_report"
    assert doc["giap"]["status"] == "fail" and doc["end_to_end"] is False


def test_te7_peptides_a_not_the_pairs_input_exit_2(tmp_path: Path, repo: Any) -> None:
    """S2: ``end_to_end`` rests on ``--peptides-a`` being the table the GIAP pair read; a
    ``--peptides-a`` whose sha256 differs from the manifests' ``peptides`` input is exit 2."""
    ref, rep = _pair(tmp_path, repo)
    out = tmp_path / "metrics"
    extra = _diann_dirs(tmp_path, False, True)  # manifests keep their own peptides sha256
    res = _te7(ref, rep, repo, out, *extra)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert _nothing_written(out)


@pytest.mark.parametrize("drop", ["--diann-a", "--diann-b", "--peptides-a", "--peptides-b"])
def test_te7_some_diann_flags_only_exit_2(tmp_path: Path, repo: Any, drop: str) -> None:
    """Test gap 5: the four DIA-NN flags go together; three of them is a precondition."""
    ref, rep = _pair(tmp_path, repo)
    extra = _diann_dirs(tmp_path, False, True, (ref, rep))
    i = extra.index(drop)
    del extra[i : i + 2]
    out = tmp_path / "metrics"
    res = _te7(ref, rep, repo, out, *extra)
    assert res.returncode == 2, f"exit {res.returncode}\n{res.stdout}\n{res.stderr}"
    assert _nothing_written(out)


# ------------------------------------------------------------------ run labels (N3)


def test_te7_two_fresh_replicates_are_labelled_rep1_and_rep2(tmp_path: Path, repo: Any) -> None:
    """N3: with two fresh replicates (``out_a2_rep1``, ``out_a2_rep2``) both labels are the
    ``rep\\d+`` of the directory; a file only in the first is ``only_in:rep1``, not ``ref``."""
    ref_files = base_tree()
    ref_files["tables/extra_first.tsv"] = b"a\n1\n"
    ref, rep = _pair(
        tmp_path, repo, base_tree(), ref_files=ref_files, names=("out_a2_rep1", "out_a2_rep2")
    )
    out = tmp_path / "metrics"
    _ok(_te7(ref, rep, repo, out))
    doc = _doc(out)
    assert set(doc["pair"]) == {"rep1", "rep2"}
    assert _by_path(doc["giap"])["tables/extra_first.tsv"]["class"] == "only_in:rep1"


def test_te7_reference_without_rep_number_is_ref(tmp_path: Path, repo: Any) -> None:
    ref, rep = _pair(tmp_path, repo, base_tree())
    out = tmp_path / "metrics"
    _ok(_te7(ref, rep, repo, out))
    assert set(_doc(out)["pair"]) == {"ref", "rep1"}


# ------------------------------------------------------- producers from real traces (S1, S3, S4)

#: real trace ``name`` values: the Nextflow alias after the last ``:`` (S3)
TRACE_NAMES = (
    "GIAP:CONTEXT:CONTEXT_STATS (chr1)",
    "GIAP:VARIANTS:VARIANTS_STEP (chr1)",
    "GIAP:FUNCTION:FUNCTION_STEP",
    "GIAP:OPERONS:OPERONS_STEP",
    "GIAP:DETERMINANTS:DETERMINANTS_STEP",
    "GIAP:REPORT:PLOTS",
    "GIAP:REPORT:REPORT_HTML",
)
REGIONS_A = (
    b"region_id\tmetric\tn_genes\thl_shift_log2\teffect\tci_low\tci_high\tp_value\n"
    b"prophage_1\tdetected\t12\t0.5\t0.25\t0.1\t0.4\t0.031\n"
)
#: only M4-20 decision 9 columns move (effect, ci_low); the p-value does not (decision 9)
REGIONS_B = REGIONS_A.replace(b"\t0.25\t0.1\t", b"\t0.25000000000000006\t0.10000000000000002\t")
AUTOCORR_HEADER = (
    b"contig\tstatistic\tlag\tnull_model\tvalue\tnull_mean\tnull_lo\tnull_hi\tstat\teffect"
    b"\tci_low\tci_high\tp_value\tq_value\tmethod\tn_perm\tflags\n"
)
AUTOCORR_A = AUTOCORR_HEADER + (
    b"chr1\tjoin_count_bb\t1\tpermutation\t1382\t760.125\t700.0\t820.0\tjoin_count_bb"
    b"\t621.875\tNA\tNA\t0.0009\t0.0018\tpermutation\t1000\tNA\n"
)
#: M5-11 decisions 44/46/48: null_mean, effect (decision 46's drifting stat), p_value, q_value
AUTOCORR_B = (
    AUTOCORR_A.replace(b"\t760.125\t", b"\t760.1250000000001\t")
    .replace(b"\tjoin_count_bb\t621.875\t", b"\tjoin_count_bb\t621.8749999999999\t")
    .replace(b"\t0.0009\t0.0018\t", b"\t0.001\t0.002\t")
)
STEP_FILES = {
    "tables/variants.tsv": "variants",
    "tables/category_summary.tsv": "function",
    "tables/operons.tsv": "operons",
    "tables/determinants_features.tsv": "determinants",
}


def _producer_pair(
    tmp: Path, repo: Any, rep_trace: tuple[str, ...] | None = TRACE_NAMES
) -> dict[str, Any]:
    """Ref on rye, replicate on wheat; every STEP_FILES table and the two LAPACK tables differ.
    ``rep_trace`` None: the replicate has no ``logs/trace.tsv`` at all."""
    a, b = base_tree(), base_tree()
    a["tables/regions_stats.tsv"], b["tables/regions_stats.tsv"] = REGIONS_A, REGIONS_B
    a["tables/autocorrelation.tsv"], b["tables/autocorrelation.tsv"] = AUTOCORR_A, AUTOCORR_B
    for path in STEP_FILES:
        a[path] = b"id\tscore\nx\t1.0\n"
        b[path] = b"id\tscore\nx\t1.5\n"
    ref, rep = _pair(tmp, repo, b, ref_files=a)
    write_named_trace(ref, [(n, "rye-02") for n in TRACE_NAMES])
    if rep_trace is None:
        (rep / "logs" / "trace.tsv").unlink()
    else:
        write_named_trace(rep, [(n, "wheat-03") for n in rep_trace])
    out = tmp / "metrics"
    _ok(_te7(ref, rep, repo, out))
    return _by_path(_doc(out)["giap"])


def test_producers_from_real_trace_names_and_aliases(tmp_path: Path, repo: Any) -> None:
    """S3: ``*_STEP`` aliases map to their registry steps; each table reports its producer's
    CPU classes per run."""
    d = _producer_pair(tmp_path, repo)
    both = {"ref": ["rye"], "rep1": ["wheat"]}
    for path, step in STEP_FILES.items():
        assert d[path]["producer"] == step, path
        assert d[path]["cpu_classes"] == both, path
        assert d[path]["sub_class"] is None, path  # score is not a listed LAPACK column
    rs = d["tables/regions_stats.tsv"]
    assert rs["producer"] == "context_stats" and rs["cpu_classes"] == both
    assert rs["class"] == "numeric_only" and rs["sub_class"] == "known_lapack"
    assert rs["columns"] == ["effect", "ci_low"]


def test_autocorrelation_permutation_p_values_are_known_lapack(tmp_path: Path, repo: Any) -> None:
    """S4 (M5-11 decision 48): the join-count p_value and q_value move with the BLAS-drifted
    observed value; with null_mean and effect they are all in te7_known_columns.tsv."""
    ac = _producer_pair(tmp_path, repo)["tables/autocorrelation.tsv"]
    assert ac["producer"] == "context_stats"
    assert ac["class"] == "numeric_only"
    assert ac["columns"] == ["null_mean", "effect", "p_value", "q_value"]
    assert ac["sub_class"] == "known_lapack"


def test_step_missing_from_one_trace_is_not_known_lapack(tmp_path: Path, repo: Any) -> None:
    """S1 / test gap 1: the replicate's trace has no CONTEXT_STATS task."""
    rep_trace = tuple(n for n in TRACE_NAMES if "CONTEXT_STATS" not in n)
    rs = _producer_pair(tmp_path, repo, rep_trace)["tables/regions_stats.tsv"]
    assert rs["class"] == "numeric_only" and rs["producer"] == "context_stats"
    assert rs["cpu_classes"] == {"ref": ["rye"], "rep1": []}
    assert rs["sub_class"] is None


def test_missing_trace_file_is_not_known_lapack(tmp_path: Path, repo: Any) -> None:
    """S1 / test gap 1: the replicate has no logs/trace.tsv."""
    d = _producer_pair(tmp_path, repo, None)
    for path in ("tables/regions_stats.tsv", "tables/autocorrelation.tsv"):
        assert d[path]["class"] == "numeric_only", path
        assert d[path]["cpu_classes"] == {"ref": ["rye"], "rep1": []}, path
        assert d[path]["sub_class"] is None, path


# ------------------------------------------------------------------ criterion 10 (te7 half)


def test_te7_twice_is_byte_identical(tmp_path: Path, repo: Any) -> None:
    rep_files = ex8_replicate()
    ref, rep = _pair(tmp_path, repo, rep_files)
    o1, o2 = tmp_path / "m1", tmp_path / "m2"
    extra = _diann_dirs(tmp_path, False, False, (ref, rep))
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


def _known_rows() -> tuple[list[str], list[dict[str, str]]]:
    lines = KNOWN_COLUMNS.read_text(encoding="utf-8").splitlines()
    comments = [ln for ln in lines if ln.startswith("#")]
    rows = list(csv.DictReader((ln for ln in lines if not ln.startswith("#")), delimiter="\t"))
    return comments, rows


def test_known_columns_exact_sets_per_table() -> None:
    """S4 (lead decision): regions_stats keeps its M4-20 decision 9 lstsq columns; autocorrelation
    lists null_mean, null_lo, null_hi, effect (decisions 44, 46, 47) and p_value, q_value
    (decision 48), and drops ``value`` (an integer join count)."""
    _, rows = _known_rows()
    by: dict[str, set[str]] = {}
    for r in rows:
        by.setdefault(r["table"], set()).add(r["column"])
    assert by == {
        "tables/regions_stats.tsv": {"hl_shift_log2", "effect", "ci_low", "ci_high"},
        "tables/autocorrelation.tsv": {
            "null_mean",
            "null_lo",
            "null_hi",
            "effect",
            "p_value",
            "q_value",
        },
    }
    for r in rows:
        if r["table"] == "tables/autocorrelation.tsv" and r["column"] in ("p_value", "q_value"):
            assert "48" in r["source"], r
    assert len(rows) == len({(r["table"], r["column"]) for r in rows})  # no duplicate row


def test_known_columns_comment_places_decision_46_in_effect() -> None:
    """S4: the comment says decision 46's drifting value is today's ``effect`` column (not a
    string column)."""
    comments, _ = _known_rows()
    text = " ".join(c.lstrip("# ") for c in comments)
    assert re.search(r"46.{0,160}effect|effect.{0,160}46", text), comments
    assert "string column" not in text


def test_numeric_cells_numpy_float_parsing_sanity() -> None:
    """Self-check of the example 8 fixture: every differing genes.tsv cell parses as a float."""
    a = [ln.split("\t") for ln in GENES_TSV.decode().splitlines()[1:]]
    b = [ln.split("\t") for ln in ex8_replicate()["tables/genes.tsv"].decode().splitlines()[1:]]
    diff = [
        (x, y) for ra, rb in zip(a, b, strict=True) for x, y in zip(ra, rb, strict=True) if x != y
    ]
    assert len(diff) == 3
    assert all(np.isfinite(float(x)) and np.isfinite(float(y)) for x, y in diff)
