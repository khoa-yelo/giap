"""The K-12 calibration driver ``validation/calibrate_cutoffs.py`` (issue M6-12, part B).

Spec ``docs/specs/M6-12-uniop-terminator-cutoffs.md`` (examples 6-9, the degenerate table, the
acceptance criteria, lead decisions 1-7) and methods ``docs/methods/operons.md`` "K-12 cut-off
calibration (M6-12)" (the authority where it is more specific). Offline: every input is a toy
written under ``tmp_path`` by ``tests/helpers/k12_cutoffs_world.py``; no RegulonDB bytes, no
``$GIAP_ROOT/recordings``, no ``$GIAP_DB_DIR`` (the env vars point at a missing directory here).

Interface pinned here (``validation.calibrate_cutoffs`` imported as ``cc``):

- ``STEP == "m6_12_cutoffs"``; ``UNIOP_PIN`` (40-hex ADR-0008 commit); ``K12_PINS: Pins``;
- ``Pins(terminatornet_commit, genome_sha256, uniop_pred_sha256, uniop_gff_sha256)``, keyword
  constructible, the expected recording identities (``K12_PINS`` holds the real K-12 values);
- ``normalized_genome_sha256(path) -> str``: SHA-256 of the genome as preflight's
  ``genome.norm.fa`` (``write_fasta`` of ``read_genome_fasta``; gz / line width do not matter);
- ``terminatornet_commit(terms_path) -> str``: the 40-hex ``tool commit `<sha>``` in the
  ``PROVENANCE.md`` beside the terms file; ``ValueError`` when the file or the line is missing;
- ``coordinate_id_map(gff, genes) -> dict[str, str]``: GFF ID -> locus tag by
  ``(start1, end, strand, contig) == (interval.start + 1, interval.end, strand, contig)``; a row
  matching no gene or two genes is a ``ValueError`` naming that row's ID (example 8);
- ``directon_ids(genes, replicon) -> dict[str, int]``: locus tag -> directon, over the whole gene
  model; maximal same-strand runs in genome order, last and first joined on a circular replicon
  when they share a strand;
- ``pair_table(genes, replicon, truth, probs) -> DataFrame``: the labelled rows of
  ``pair_truth.tsv`` (positive / negative) in file order with columns ``gene_a, gene_b, label,
  y, uniop_p`` (NaN without a UniOP p), ``d`` (``intergenic_distance``) and ``directon_id``;
  ``ValueError`` on a locus tag absent from the model or a pair that is not adjacent same-strand;
- ``truth_recall(calls, confidence, truth, windows, replicon, threshold, *, tol_bp=50) ->
  (k, n)``: n = truth intervals inside a same-strand window (containment, any shift); k = those
  within ``tol_bp`` (interval gap) of at least one *eligible* call with confidence >= threshold;
  a terminator counts once (decision 7, gap 10);
- ``run(regulondb_dir, bakta_json, genome, uniop_dir, terminatornet, out_dir, *, seed=20261014,
  B=2000, R=200, pins=K12_PINS, contig="U00096.3", min_labelled=500) -> dict``: verifies every
  identity, then writes ``cutoffs.json``, ``uniop_curve.tsv``, ``terminatornet_curve.tsv`` and
  ``terminatornet_calls.tsv`` into ``out_dir`` and returns the ``cutoffs.json`` content. Raises
  ``ValueError`` (writing nothing) on an identity mismatch, an unmatched GFF row, the < 500 /
  one-class guard, an empty ``uniop.pred`` or more than 1 % of labelled pairs without a UniOP p;
  reads truth from ``regulondb_dir/truth/``; decompresses ``terms.tsv.gz`` itself;
- ``main(argv) -> int``: flags ``--regulondb-dir --bakta-json --genome --uniop-dir
  --terminatornet --out-dir`` (required), ``--seed``, ``--B``; uses ``K12_PINS`` looked up at
  call time; 0 on success, 1 on any of ``run``'s ``ValueError`` (message on stderr, no
  ``cutoffs.json``); argparse usage errors exit 2.

Bootstrap stream order (methods "Bootstrap streams"): one ``rng = for_step(STEP, seed)``; UniOP
AUROC, UniOP precision then recall at the fixed t*, distance AUROC (-d), the paired AUROC
difference, d0, s, then the same seven with ``row_id`` clusters (iid). Anything else (the t*
quantiles, the coincidence null, the CV folds) draws after these.
"""

from __future__ import annotations

import gzip
import importlib
import json
import math
import re
import sys
from collections.abc import Sequence
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pytest
from helpers.k12_cutoffs_world import (
    CIRC10K,
    CONTIG,
    EX6_CALLS,
    EX6_GENES,
    EX6_LABELS,
    EX6_PROBS,
    EX6_TRUTH,
    LIN10K,
    TN_COMMIT,
    UNIOP_PIN,
    BigSpec,
    Call,
    G,
    World,
    big_world,
    build_world,
    example6,
    example7,
    gene_obj,
    refresh_pins,
    sha256,
    write_uniop,
)
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.annotate.operons import intergenic_distance
from giap.annotate.uniop_adapter import PairProb
from giap.io.fasta import write_fasta
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene
from giap.rng import for_step
from giap.stats.cutoffs import (
    auroc,
    auroc_ci,
    auroc_diff_ci,
    choose_cutoff,
    distance_fit_ci,
    distance_log_loss,
    downstream_windows,
    fit_distance_model,
    fixed_cutoff_ci,
    interval_gap,
    precision_lower_bound,
    terminator_match,
    threshold_curve,
    wilson_ci,
    window_of,
)

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
cc: Any = importlib.import_module("validation.calibrate_cutoffs")

FAST = {"B": 30, "R": 5}


@pytest.fixture(autouse=True)
def _no_real_data(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    """Point every data root at a missing directory: the driver must only read its arguments."""
    missing = tmp_path / "no-such-root"
    for var in ("GIAP_ROOT", "GIAP_DB_DIR", "GIAP_DATA_DIR"):
        monkeypatch.setenv(var, str(missing))


def run(w: World, out: Path, **kw: Any) -> dict[str, Any]:
    kw.setdefault("pins", cc.Pins(**w.pins))
    for k, v in FAST.items():
        kw.setdefault(k, v)
    result: dict[str, Any] = cc.run(*w.run_args(out), **kw)
    return result


def _strict(text: str) -> Any:
    def bad(c: str) -> None:
        raise AssertionError(f"non-finite JSON constant {c}")

    return json.loads(text, parse_constant=bad)


def cutoffs(out: Path) -> dict[str, Any]:
    doc: dict[str, Any] = _strict((out / "cutoffs.json").read_text(encoding="utf-8"))
    return doc


def tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def ci(b: Any) -> list[float]:
    return [b.ci_low, b.ci_high]


def approx_ci(got: Any, want: list[float]) -> None:
    assert isinstance(got, list) and len(got) == 2
    for g, w in zip(got, want, strict=True):
        if math.isnan(w):
            assert g is None
        else:
            assert g == pytest.approx(w, rel=1e-9, abs=1e-12)


def gene(tag: str, start: int, end: int, strand: str, contig: str = CONTIG) -> Gene:
    return gene_obj(G(tag, start, end, strand, contig))


def iv(start: int, end: int, strand: str) -> Interval:
    return Interval(CONTIG, start, end, strand)  # type: ignore[arg-type]


# ===================================================================== pins and identities


def test_k12_pins_are_the_recorded_identities() -> None:
    p = cc.K12_PINS
    assert p.terminatornet_commit == "72180309545988dba63b8bd5db5d0b09bf4dba70"
    assert p.genome_sha256 == "918f35c7bfc8a3c8d824ed281a4185903025a0361242219b6b0e0ea6ecaa14f2"
    assert p.uniop_pred_sha256 == (
        "7d852a93a436a9bcef2cf4cf38d608bcd5d59054d98bd1b1ddc21196690b5ee0"
    )
    assert p.uniop_gff_sha256 == "ea5870894e26c649b4e7520c65cf57e1138e1786688db7663a7b3e3eea062906"
    assert cc.UNIOP_PIN == UNIOP_PIN
    assert cc.STEP == "m6_12_cutoffs"


def test_normalized_genome_sha_of_the_committed_k12_genome_is_the_pin() -> None:
    """The pin is preflight's genome.norm.fa, not the committed .fna.gz bytes."""
    genome = REPO / "tests/data/genomes/k12_mg1655.fna.gz"
    assert cc.normalized_genome_sha256(genome) == cc.K12_PINS.genome_sha256


def test_normalized_genome_sha_ignores_gzip_and_line_width(tmp_path: Path) -> None:
    seq = "ACGT" * 300
    a = tmp_path / "a.fa"
    write_fasta(a, [("U00096.3", seq)], width=80)
    b = tmp_path / "b.fa.gz"
    with gzip.open(b, "wt") as fh:
        fh.write(">U00096.3 some description\n")
        fh.writelines(seq[i : i + 60].lower() + "\n" for i in range(0, len(seq), 60))
    assert cc.normalized_genome_sha256(a) == sha256(a)
    assert cc.normalized_genome_sha256(b) == sha256(a)


def test_terminatornet_commit_from_the_committed_provenance() -> None:
    terms = REPO / "tests/data/recorded/terminatornet/721803095459/k12/terms.tsv.gz"
    assert cc.terminatornet_commit(terms) == cc.K12_PINS.terminatornet_commit


def test_terminatornet_commit_missing_provenance_raises(tmp_path: Path) -> None:
    terms = tmp_path / "terms.tsv.gz"
    terms.write_bytes(gzip.compress(b""))
    with pytest.raises(ValueError, match="PROVENANCE"):
        cc.terminatornet_commit(terms)


# ================================================================ example 8: coordinate map


def _gff(path: Path, rows: Sequence[tuple[str, int, int, str]]) -> Path:
    lines = ["##gff-version 3"]
    lines += [f"{CONTIG}\tGIAP\tCDS\t{s}\t{e}\t.\t{st}\t0\tID={i}" for i, s, e, st in rows]
    path.write_text("\n".join(lines) + "\n")
    return path


def test_example8_coordinate_map_ignores_the_neighbours_tag(tmp_path: Path) -> None:
    """Map row c_1 -> T1 is the shifted-tag defect; coordinates give c_1 -> T2."""
    genes = [gene("T1", 2, 98, "+"), gene("T2", 189, 255, "+")]
    gff = _gff(tmp_path / "u.gff", [("c_1", 190, 255, "+")])
    assert cc.coordinate_id_map(gff, genes) == {"c_1": "T2"}  # T1 has no row: allowed


@pytest.mark.parametrize(
    ("row", "genes"),
    [
        (("c_1", 191, 255, "+"), [gene("T2", 189, 255, "+")]),  # shifted by 1 bp
        (("c_1", 190, 255, "-"), [gene("T2", 189, 255, "+")]),  # strand differs
        (("c_1", 190, 255, "+"), [gene("T2", 189, 255, "+"), gene("T3", 189, 255, "+")]),
    ],
    ids=["no-match", "wrong-strand", "two-genes"],
)
def test_example8_unmatched_or_ambiguous_row_raises_naming_it(
    tmp_path: Path, row: tuple[str, int, int, str], genes: list[Gene]
) -> None:
    gff = _gff(tmp_path / "u.gff", [row])
    with pytest.raises(ValueError, match="c_1"):
        cc.coordinate_id_map(gff, genes)


def test_first_failing_row_is_named(tmp_path: Path) -> None:
    genes = [gene("T1", 9, 98, "+"), gene("T2", 189, 255, "+")]
    gff = _gff(tmp_path / "u.gff", [("c_1", 10, 98, "+"), ("c_2", 7, 9, "+"), ("c_3", 1, 5, "+")])
    with pytest.raises(ValueError, match=r"\bc_2\b"):
        cc.coordinate_id_map(gff, genes)


def test_run_stops_on_a_shifted_gff_and_writes_nothing(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    gff = w.uniop_dir / "uniop_cds.gff"
    lines = gff.read_text().splitlines()
    f = lines[2].split("\t")  # second gene, g2
    f[3] = str(int(f[3]) + 1)
    lines[2] = "\t".join(f)
    gff.write_text("\n".join(lines) + "\n")
    refresh_pins(w)
    out = tmp_path / "out"
    with pytest.raises(ValueError, match=f"{CONTIG}_2"):
        run(w, out, min_labelled=2)
    assert not (out / "cutoffs.json").exists()


def test_run_never_reads_uniop_map(tmp_path: Path) -> None:
    """The recorded map is shifted by one; deleting it changes nothing."""
    w = example6(tmp_path / "w")
    a = run(w, tmp_path / "a", min_labelled=2)
    (w.uniop_dir / "uniop_map.tsv").unlink()
    b = run(w, tmp_path / "b", min_labelled=2)
    assert a["uniop"]["auroc"] == b["uniop"]["auroc"] == 1.0
    assert a["uniop"]["threshold"] == b["uniop"]["threshold"] == 0.9


# ======================================================================= directons and pairs


def test_directons_linear_split_at_strand_changes() -> None:
    genes = [
        gene("a", 0, 100, "+"),
        gene("b", 200, 300, "+"),
        gene("c", 400, 500, "-"),
        gene("d", 600, 700, "-"),
        gene("e", 800, 900, "+"),
    ]
    got = cc.directon_ids(genes, LIN10K)
    assert got["a"] == got["b"] != got["c"] == got["d"] != got["e"]
    assert got["a"] != got["e"]


def test_directons_circular_join_the_last_and_first_runs() -> None:
    genes = [
        gene("a", 0, 100, "+"),
        gene("b", 200, 300, "+"),
        gene("c", 400, 500, "-"),
        gene("e", 800, 900, "+"),
    ]
    got = cc.directon_ids(genes, CIRC10K)
    assert got["e"] == got["a"] == got["b"] != got["c"]
    assert len(set(got.values())) == 2


def test_directons_input_order_does_not_matter() -> None:
    genes = [gene("e", 800, 900, "+"), gene("c", 400, 500, "-"), gene("a", 0, 100, "+")]
    got = cc.directon_ids(genes, CIRC10K)
    assert got["a"] == got["e"] != got["c"]


@settings(max_examples=60, deadline=None)
@given(strands=st.lists(st.sampled_from("+-"), min_size=1, max_size=25), circ=st.booleans())
def test_directon_count_property(strands: list[str], circ: bool) -> None:
    """Linear: changes + 1. Circular: the ring's changes (>= 1 run); each pair in one directon."""
    rep = CIRC10K if circ else LIN10K
    genes = [gene(f"x{i:02d}", 100 * i, 100 * i + 50, s) for i, s in enumerate(strands)]
    got = cc.directon_ids(genes, rep)
    assert set(got) == {g.locus_tag for g in genes}
    n = len(strands)
    lin_changes = sum(strands[i] != strands[i + 1] for i in range(n - 1))
    if circ:
        ring = lin_changes + (strands[-1] != strands[0])
        want = max(1, ring)
    else:
        want = lin_changes + 1
    assert len(set(got.values())) == want
    for i in range(n - 1):
        same = strands[i] == strands[i + 1]
        assert (got[f"x{i:02d}"] == got[f"x{i + 1:02d}"]) is same


def _truth_df(rows: Sequence[tuple[str, str, str]]) -> pd.DataFrame:
    return pd.DataFrame(
        [{"gene_a": a, "gene_b": b, "label": lab, "tu_ids": ""} for a, b, lab in rows]
    )


def test_example7_labelled_wrap_pair_distance_crosses_the_origin() -> None:
    """gB -> gA has d = 200 - 9900 + 10000 = 300 and shares gA's directon; UniOP never scores
    the wrap pair, so its p is NaN."""
    genes = [gene("gA", 200, 800, "+"), gene("gB", 9500, 9900, "+")]
    truth = _truth_df([("gA", "gB", "negative"), ("gB", "gA", "positive")])
    probs = [PairProb("gA", "gB", 0.4, "uniop")]
    tab = cc.pair_table(genes, CIRC10K, truth, probs)
    assert list(tab["gene_a"]) == ["gA", "gB"]
    assert list(tab["y"]) == [0, 1]
    assert list(tab["d"]) == [8700, 300]
    assert tab["uniop_p"].iloc[0] == 0.4
    assert math.isnan(tab["uniop_p"].iloc[1])
    assert tab["directon_id"].iloc[0] == tab["directon_id"].iloc[1]


def test_pair_table_drops_ambiguous_but_directons_ignore_labels() -> None:
    """Directons come from the whole model *before* ambiguous pairs are dropped: an ambiguous
    middle pair does not split a + + + + run."""
    genes = [gene(t, 1000 * i, 1000 * i + 600, "+") for i, t in enumerate("abcd")]
    truth = _truth_df([("a", "b", "positive"), ("b", "c", "ambiguous"), ("c", "d", "negative")])
    probs = [PairProb(x, y, p, "uniop") for (x, y), p in {("a", "b"): 0.9, ("c", "d"): 0.1}.items()]
    tab = cc.pair_table(genes, LIN10K, truth, probs)
    assert list(zip(tab["gene_a"], tab["gene_b"], strict=True)) == [("a", "b"), ("c", "d")]
    assert set(tab.columns) >= {"gene_a", "gene_b", "label", "y", "uniop_p", "d", "directon_id"}
    assert tab["directon_id"].nunique() == 1
    assert list(tab["d"]) == [400, 400]


def test_pair_table_minus_strand_distance_is_genome_order() -> None:
    genes = [gene("m1", 100, 700, "-"), gene("m2", 650, 1300, "-")]
    tab = cc.pair_table(genes, LIN10K, _truth_df([("m1", "m2", "positive")]), [])
    assert list(tab["d"]) == [-50]
    assert math.isnan(tab["uniop_p"].iloc[0])


@pytest.mark.parametrize(
    "row",
    [("a", "zz", "positive"), ("a", "c", "positive")],
    ids=["unknown-locus", "not-adjacent"],
)
def test_pair_table_rejects_truth_from_another_model(row: tuple[str, str, str]) -> None:
    genes = [gene(t, 1000 * i, 1000 * i + 600, "+") for i, t in enumerate("abc")]
    with pytest.raises(ValueError):
        cc.pair_table(genes, LIN10K, _truth_df([row]), [])


# ========================================================= truth-level recall (decision 7)


def test_truth_recall_example6() -> None:
    genes = [gene_obj(g) for g in EX6_GENES]
    windows = downstream_windows(genes, LIN10K)
    calls = [iv(c.start, c.end, c.strand) for c in EX6_CALLS]
    conf = [c.conf for c in EX6_CALLS]
    truth = [iv(t.start, t.end, t.strand) for t in EX6_TRUTH]
    assert cc.truth_recall(calls, conf, truth, windows, LIN10K, 0.95) == (1, 3)
    assert cc.truth_recall(calls, conf, truth, windows, LIN10K, 0.85) == (2, 3)
    assert cc.truth_recall(calls, conf, truth, windows, LIN10K, 0.0) == (2, 3)  # c2: gap 70


def test_truth_recall_counts_a_terminator_once_and_ignores_ineligible_calls() -> None:
    genes = [gene_obj(g) for g in EX6_GENES]
    windows = downstream_windows(genes, LIN10K)
    truth = [iv(1120, 1150, "+"), iv(1470, 1490, "+")]  # both inside g1's window [950, 1500)
    calls = [
        iv(1100, 1140, "+"),  # matches truth 0
        iv(1150, 1180, "+"),  # matches truth 0 again
        iv(1460, 1520, "+"),  # overlaps truth 1 but crosses g2's start: ineligible
    ]
    conf = [0.9, 0.95, 0.99]
    assert cc.truth_recall(calls, conf, truth, windows, LIN10K, 0.5) == (1, 2)


def test_truth_recall_any_call_within_50_not_only_the_nearest() -> None:
    """A call 10 bp from one terminator and 40 bp from another recalls both."""
    genes = [gene("p", 100, 1000, "+"), gene("q", 3000, 4000, "+")]
    windows = downstream_windows(genes, LIN10K)
    truth = [iv(1210, 1220, "+"), iv(1100, 1160, "+")]
    calls = [iv(1180, 1200, "+")]  # gap 10 to truth 0 and 20 to truth 1
    assert cc.truth_recall(calls, [0.7], truth, windows, LIN10K, 0.5) == (2, 2)


def test_truth_recall_circular_shift() -> None:
    genes = [gene("gA", 200, 800, "+"), gene("gB", 9500, 9900, "+")]
    windows = downstream_windows(genes, CIRC10K)
    assert cc.truth_recall([iv(20, 60, "+")], [0.8], [iv(80, 100, "+")], windows, CIRC10K, 0.5) == (
        1,
        1,
    )


# ============================================================ example 6 through the driver


def test_example6_terminatornet_through_the_driver(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    out = tmp_path / "out"
    res = run(w, out, min_labelled=2)
    doc = cutoffs(out)
    assert res == doc
    tn = doc["terminatornet"]
    assert tn["commit"] == TN_COMMIT
    assert tn["terms_sha256"] == sha256(w.terminatornet)
    assert (tn["n_calls"], tn["n_eligible"], tn["n_truth"], tn["n_truth_eligible"]) == (6, 4, 4, 3)
    assert tn["threshold"] == 0.95
    assert tn["rule"] == "precision_0.90"
    assert (tn["tp"], tn["fp"]) == (1, 0)
    assert tn["precision"] == 1.0
    assert tn["precision_is_lower_bound"] is True
    assert tn["recall"] == pytest.approx(1 / 3, abs=1e-15)
    assert tn["recall_all"] == pytest.approx(1 / 4, abs=1e-15)
    approx_ci(tn["precision_ci"], list(wilson_ci(1, 1)))
    approx_ci(tn["recall_ci"], list(wilson_ci(1, 3)))
    assert tn["auroc_diagnostic"] == pytest.approx(0.75, abs=1e-15)
    assert tn["calibrated"] is True
    # recall by RegulonDB type, over eligible truth (T1, T2 rho-independent; T3 rho-dependent)
    assert tn["recall_by_type"]["rho-independent"]["recall"] == pytest.approx(0.5)
    assert tn["recall_by_type"]["rho-dependent"]["recall"] == 0.0
    assert "coincidence_rate" in tn
    if tn["coincidence_rate"] is not None:
        want = precision_lower_bound(tn["precision"], tn["coincidence_rate"])
        if math.isnan(want):  # c = 1: undefined, written as null
            assert tn["precision_lower_bound"] is None
        else:
            assert tn["precision_lower_bound"] == pytest.approx(want)


def test_example6_calls_tsv(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    out = tmp_path / "out"
    run(w, out, min_labelled=2)
    calls = tsv(out / "terminatornet_calls.tsv")
    assert list(calls.columns) == [
        "contig",
        "start",
        "end",
        "strand",
        "confidence",
        "gene",
        "matched_source_id",
    ]
    got = [
        (
            r.contig,
            int(r.start),
            int(r.end),
            r.strand,
            float(r.confidence),
            r.gene,
            r.matched_source_id,
        )
        for r in calls.itertuples()
    ]
    assert got == [
        (CONTIG, 1100, 1140, "+", 0.95, "g1", "T1"),
        (CONTIG, 2460, 2490, "+", 0.90, "g2", ""),
        (CONTIG, 2600, 2640, "-", 0.85, "g3", "T3"),
        (CONTIG, 4100, 4130, "-", 0.60, "g4", ""),
    ]


def test_example6_terminatornet_curve_tsv(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    out = tmp_path / "out"
    run(w, out, min_labelled=2)
    curve = tsv(out / "terminatornet_curve.tsv")
    assert list(curve.columns) == ["threshold", "tp", "fp", "precision", "recall", "f1"]
    assert [float(x) for x in curve["threshold"]] == [0.95, 0.90, 0.85, 0.60]
    assert [int(x) for x in curve["tp"]] == [1, 1, 2, 2]
    assert [int(x) for x in curve["fp"]] == [0, 1, 1, 2]
    assert [float(x) for x in curve["precision"]] == pytest.approx([1.0, 0.5, 2 / 3, 0.5])


def test_example6_uniop_and_refit_halves(tmp_path: Path) -> None:
    """Two labelled pairs (min_labelled lowered): d = 500 (positive) and 1000 (negative) are
    separated, so the refit fails and keeps (50, 20)."""
    w = example6(tmp_path / "w")
    out = tmp_path / "out"
    doc = run(w, out, min_labelled=2)
    u = doc["uniop"]
    assert (u["n_labelled"], u["n_positive"], u["n_negative"], u["n_without_p"]) == (2, 1, 1, 0)
    assert u["pin"] == UNIOP_PIN
    assert u["pred_sha256"] == sha256(w.uniop_dir / "uniop.pred")
    assert u["gff_sha256"] == sha256(w.uniop_dir / "uniop_cds.gff")
    assert u["auroc"] == 1.0
    assert (u["threshold"], u["rule"], u["tp"], u["fp"]) == (0.9, "precision_0.90", 1, 0)
    assert (u["precision"], u["recall"]) == (1.0, 1.0)
    assert u["precision_cv"] is None  # 2 directons < 10 folds: not computed, not a crash
    assert doc["n_directons"] == 2
    dm = doc["distance_model"]
    assert dm["converged"] is False
    assert (dm["d0"], dm["s"]) == (50.0, 20.0)
    assert dm["reason"] == "separation"
    assert dm["logloss_placeholder"] == pytest.approx(
        distance_log_loss([500, 1000], [1, 0], 50, 20)
    )
    assert dm["logloss"] == dm["logloss_placeholder"]
    assert dm["auroc"] == 1.0  # AUROC of -d


# ============================================================ example 7 through the driver


def test_example7_wrap_window_eligible_and_matched(tmp_path: Path) -> None:
    w = example7(tmp_path / "w")
    out = tmp_path / "out"
    doc = run(w, out, min_labelled=3)
    tn = doc["terminatornet"]
    assert (tn["n_calls"], tn["n_eligible"], tn["n_truth"], tn["n_truth_eligible"]) == (1, 1, 1, 1)
    assert (tn["threshold"], tn["tp"], tn["fp"], tn["recall"]) == (0.8, 1, 0, 1.0)
    calls = tsv(out / "terminatornet_calls.tsv")
    row = calls.iloc[0]
    assert (int(row.start), int(row.end), row.gene, row.matched_source_id) == (20, 60, "gB", "TW")
    # gB joins gA's directon through the origin: {gB, gA, gC, gD} and {gE, gF}
    assert doc["n_directons"] == 2
    assert doc["uniop"]["n_labelled"] == 3


# ============================================================== degenerate rows (part B)


def test_fewer_than_500_labelled_pairs_raises(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run(w, out)  # default min_labelled = 500
    assert not (out / "cutoffs.json").exists()


def test_one_class_raises(tmp_path: Path) -> None:
    labels = {k: "negative" for k in EX6_LABELS}
    w = build_world(tmp_path / "w", EX6_GENES, [LIN10K], labels, EX6_PROBS, EX6_CALLS, EX6_TRUTH)
    with pytest.raises(ValueError):
        run(w, tmp_path / "out", min_labelled=2)


def test_empty_uniop_pred_escalates(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    (w.uniop_dir / "uniop.pred").write_text(
        "operon prediction scores\nGene A\tGene B\tPrediction\n"
    )
    refresh_pins(w)
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run(w, out, min_labelled=2)
    assert not (out / "cutoffs.json").exists()


def _drop_probs(w: World, spec: BigSpec, n: int) -> list[tuple[str, str]]:
    labelled = [(a, b) for a, b, lab in w.pairs if lab != "ambiguous" and (a, b) in spec.probs]
    gone = labelled[5 : 5 + n]
    probs = {k: v for k, v in spec.probs.items() if k not in gone}
    for f in ("uniop.pred", "uniop_cds.gff", "uniop_map.tsv"):
        (w.uniop_dir / f).unlink()
    write_uniop(w.uniop_dir, w.genes, w.replicons, probs)
    refresh_pins(w)
    return gone


def test_more_than_one_percent_without_p_escalates(tmp_path: Path) -> None:
    w, spec = big_world(tmp_path / "w")
    n_lab = sum(1 for p in w.pairs if p[2] != "ambiguous")
    _drop_probs(w, spec, math.floor(0.01 * n_lab) + 1)
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run(w, out)
    assert not (out / "cutoffs.json").exists()


def test_one_pair_without_p_is_excluded_and_counted(tmp_path: Path) -> None:
    w, spec = big_world(tmp_path / "w")
    n_lab = sum(1 for p in w.pairs if p[2] != "ambiguous")
    _drop_probs(w, spec, 1)
    u = run(w, tmp_path / "out")["uniop"]
    assert u["n_without_p"] == 1
    assert u["n_labelled"] == n_lab - 1 == u["n_positive"] + u["n_negative"]


def test_no_eligible_calls_leaves_terminatornet_uncalibrated(tmp_path: Path) -> None:
    calls = [EX6_CALLS[3]]  # c4: in no + window
    w = build_world(tmp_path / "w", EX6_GENES, [LIN10K], EX6_LABELS, EX6_PROBS, calls, EX6_TRUTH)
    doc = run(w, tmp_path / "out", min_labelled=2)
    assert doc["terminatornet"]["calibrated"] is False
    assert doc["terminatornet"]["n_eligible"] == 0
    assert doc["uniop"]["threshold"] == 0.9  # the UniOP half still lands


def test_zero_matches_leaves_terminatornet_uncalibrated(tmp_path: Path) -> None:
    truth = [EX6_TRUTH[3]]  # T4 only: on the other strand from c5, in no + window
    w = build_world(tmp_path / "w", EX6_GENES, [LIN10K], EX6_LABELS, EX6_PROBS, EX6_CALLS, truth)
    doc = run(w, tmp_path / "out", min_labelled=2)
    tn = doc["terminatornet"]
    assert tn["calibrated"] is False
    assert tn["n_eligible"] == 4
    assert doc["uniop"]["threshold"] == 0.9


def test_other_contigs_are_ignored(tmp_path: Path) -> None:
    plasmid = Replicon("pX", 3000, "linear")
    genes = [*EX6_GENES, G("p1", 100, 900, "+", "pX"), G("p2", 1000, 1600, "+", "pX")]
    calls = [*EX6_CALLS, Call(950, 990, "+", 0.99, "pX")]
    truth = [*EX6_TRUTH]
    w = build_world(tmp_path / "w", genes, [LIN10K, plasmid], EX6_LABELS, EX6_PROBS, calls, truth)
    doc = run(w, tmp_path / "out", min_labelled=2)
    tn = doc["terminatornet"]
    assert (tn["n_calls"], tn["n_eligible"], tn["threshold"]) == (6, 4, 0.95)
    assert doc["uniop"]["n_labelled"] == 2


# ===================================================================== identity mismatches


def _break_tn_commit(w: World) -> None:
    prov = w.terminatornet.parent / "PROVENANCE.md"
    prov.write_text(prov.read_text().replace(TN_COMMIT, "0" * 40))


def _drop_tn_provenance(w: World) -> None:
    (w.terminatornet.parent / "PROVENANCE.md").unlink()


def _break_genome(w: World) -> None:
    w.pins["genome_sha256"] = "0" * 64


def _break_pred(w: World) -> None:
    w.pins["uniop_pred_sha256"] = "0" * 64


def _break_gff(w: World) -> None:
    w.pins["uniop_gff_sha256"] = "0" * 64


def _break_truth_summary(w: World) -> None:
    path = w.regulondb_dir / "truth" / "truth_summary.json"
    doc = json.loads(path.read_text())
    doc["sha256"]["TerminatorSet"] = "0" * 64
    path.write_text(json.dumps(doc))


def _break_release(w: World) -> None:
    path = w.regulondb_dir / "TerminatorSet.txt"
    path.write_text(path.read_text() + "TMX\t1\t2\tforward\ta\t\trho-independent\n")


BREAKERS = {
    "terminatornet-commit": _break_tn_commit,
    "terminatornet-provenance-missing": _drop_tn_provenance,
    "genome-sha": _break_genome,
    "uniop-pred-sha": _break_pred,
    "uniop-gff-sha": _break_gff,
    "truth-summary-sha": _break_truth_summary,
    "release-provenance": _break_release,
}


@pytest.mark.parametrize("name", list(BREAKERS))
def test_identity_mismatch_raises_and_writes_nothing(tmp_path: Path, name: str) -> None:
    w = example6(tmp_path / "w")
    BREAKERS[name](w)
    out = tmp_path / "out"
    with pytest.raises(ValueError):
        run(w, out, min_labelled=2)
    assert not (out / "cutoffs.json").exists()


def test_truth_summary_mismatch_names_the_dataset(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    _break_truth_summary(w)
    with pytest.raises(ValueError, match="TerminatorSet"):
        run(w, tmp_path / "out", min_labelled=2)


# ============================================== the K-12-sized-down world: oracle and streams


@pytest.fixture(scope="module")
def big(tmp_path_factory: pytest.TempPathFactory) -> tuple[World, BigSpec, dict[str, Any], Path]:
    root = tmp_path_factory.mktemp("big")
    w, spec = big_world(root / "w")
    out = root / "out"
    doc = cc.run(*w.run_args(out), pins=cc.Pins(**w.pins), seed=20261014, **FAST)
    return w, spec, doc, out


def _oracle_L(w: World, spec: BigSpec) -> pd.DataFrame:
    by = {g.locus_tag: g for g in w.genes}
    rep = w.replicons[CONTIG]
    rows = [
        {
            "a": a,
            "b": b,
            "y": int(lab == "positive"),
            "p": spec.probs[(a, b)],
            "d": float(intergenic_distance(by[a], by[b], rep)),
            "dir": spec.directon_of[a],
        }
        for a, b, lab in w.pairs
        if lab != "ambiguous" and (a, b) in spec.probs
    ]
    return pd.DataFrame(rows)


def test_big_world_uniop_point_estimates(big: Any) -> None:
    w, spec, doc, _ = big
    L = _oracle_L(w, spec)
    u = doc["uniop"]
    assert u["n_labelled"] == len(L) >= 500
    assert (u["n_positive"], u["n_negative"], u["n_without_p"]) == (
        int(L.y.sum()),
        int((1 - L.y).sum()),
        0,
    )
    assert u["auroc"] == pytest.approx(auroc(L.p, L.y), rel=1e-12)
    cut = choose_cutoff(threshold_curve(L.p, L.y))
    assert u["threshold"] == cut.threshold  # exact: the observed score itself, never rounded
    assert u["threshold"] in set(L.p)
    assert (u["rule"], u["tp"], u["fp"]) == (cut.rule, cut.tp, cut.fp)
    assert u["literal_threshold"] == cut.literal_threshold
    assert (u["precision"], u["recall"], u["f1"]) == pytest.approx(
        (cut.precision, cut.recall, cut.f1)
    )
    assert doc["n_directons"] == L["dir"].nunique()


def test_big_world_stream_order(big: Any) -> None:
    """Streams 1-8 of the methods note, recomputed with part A's functions from one rng."""
    w, spec, doc, _ = big
    L = _oracle_L(w, spec)
    B = FAST["B"]
    rng = for_step("m6_12_cutoffs", 20261014)
    t = choose_cutoff(threshold_curve(L.p, L.y)).threshold
    p, y, d, cl = L.p.to_numpy(), L.y.to_numpy(), L.d.to_numpy(), L["dir"].to_numpy()
    rid = np.arange(len(L))
    want: dict[str, list[float]] = {}
    for suffix, c in (("", cl), ("_iid", rid)):
        want["uniop.auroc_ci" + suffix] = ci(auroc_ci(p, y, c, rng=rng, B=B))
        prec, rec = fixed_cutoff_ci(p, y, c, t, rng=rng, B=B)
        want["uniop.precision_ci" + suffix] = ci(prec)
        want["uniop.recall_ci" + suffix] = ci(rec)
        want["distance_model.auroc_ci" + suffix] = ci(auroc_ci(-d, y, c, rng=rng, B=B))
        want["distance_model.auroc_diff_ci" + suffix] = ci(auroc_diff_ci(p, -d, y, c, rng=rng, B=B))
        d0, s = distance_fit_ci(d, y, c, rng=rng, B=B)
        want["distance_model.d0_ci" + suffix] = ci(d0)
        want["distance_model.s_ci" + suffix] = ci(s)
    for key, val in want.items():
        sect, field = key.split(".")
        approx_ci(doc[sect][field], val)


def test_big_world_distance_model(big: Any) -> None:
    w, spec, doc, _ = big
    L = _oracle_L(w, spec)
    fit = fit_distance_model(L.d, L.y)
    assert fit.converged
    dm = doc["distance_model"]
    assert dm["converged"] is True and dm["reason"] is None
    assert (dm["d0"], dm["s"]) == pytest.approx((fit.d0, fit.s), rel=1e-9)
    assert dm["logloss"] == pytest.approx(distance_log_loss(L.d, L.y, fit.d0, fit.s), rel=1e-9)
    assert dm["logloss_placeholder"] == pytest.approx(distance_log_loss(L.d, L.y, 50, 20), rel=1e-9)
    assert dm["auroc"] == pytest.approx(auroc(-L.d, L.y), rel=1e-12)
    assert dm["auroc_diff"] == pytest.approx(auroc(L.p, L.y) - auroc(-L.d, L.y), rel=1e-9)


def test_big_world_terminatornet_against_part_a(big: Any) -> None:
    w, spec, doc, out = big
    rep = w.replicons[CONTIG]
    windows = downstream_windows(w.genes, rep)
    truth = [Interval(CONTIG, t.start, t.end, t.strand) for t in spec.truth]  # type: ignore[arg-type]
    elig = [c for c in spec.calls if window_of(iv(c.start, c.end, c.strand), windows, rep)]
    y = [int(terminator_match(iv(c.start, c.end, c.strand), truth, rep) is not None) for c in elig]
    conf = [c.conf for c in elig]
    cut = choose_cutoff(threshold_curve(conf, y))
    tn = doc["terminatornet"]
    assert (tn["n_calls"], tn["n_eligible"], tn["n_truth"]) == (
        len(spec.calls),
        len(elig),
        len(truth),
    )
    assert (tn["threshold"], tn["rule"], tn["tp"], tn["fp"]) == (
        cut.threshold,
        cut.rule,
        cut.tp,
        cut.fp,
    )
    approx_ci(tn["precision_ci"], list(wilson_ci(cut.tp, cut.tp + cut.fp)))
    assert tn["auroc_diagnostic"] == pytest.approx(auroc(conf, y), rel=1e-12)
    # truth-level recall, oracle written here independently of the driver's helper
    el_truth = [t for t in truth if window_of(t, windows, rep)]
    called = [iv(c.start, c.end, c.strand) for c in elig if c.conf >= cut.threshold]

    def hit(t: Interval) -> bool:
        return any(c.strand == t.strand and interval_gap(c, t, rep) <= 50 for c in called)

    k = sum(hit(t) for t in el_truth)
    assert tn["n_truth_eligible"] == len(el_truth)
    assert tn["recall"] == pytest.approx(k / len(el_truth))
    approx_ci(tn["recall_ci"], list(wilson_ci(k, len(el_truth))))
    assert tn["recall_all"] == pytest.approx(sum(hit(t) for t in truth) / len(truth))
    calls = tsv(out / "terminatornet_calls.tsv")
    assert len(calls) == len(elig)
    assert sum(calls["matched_source_id"] != "") == sum(y)
    starts = [(int(a), int(b)) for a, b in zip(calls["start"], calls["end"], strict=True)]
    assert starts == sorted(starts)


def test_big_world_curve_tsvs(big: Any) -> None:
    w, spec, doc, out = big
    L = _oracle_L(w, spec)
    want = threshold_curve(L.p, L.y)
    curve = tsv(out / "uniop_curve.tsv")
    assert list(curve.columns) == ["threshold", "tp", "fp", "precision", "recall", "f1"]
    th = [float(x) for x in curve["threshold"]]
    assert th == list(want["threshold"])  # repr round-trip, descending
    assert [int(x) for x in curve["tp"]] == list(want["tp"])
    assert [int(x) for x in curve["fp"]] == list(want["fp"])
    P = int(L.y.sum())
    f1 = [2 * a / (a + b + P) for a, b in zip(want["tp"], want["fp"], strict=True)]
    assert [float(x) for x in curve["f1"]] == pytest.approx(f1)
    assert doc["uniop"]["threshold"] in th


# ================================================================== output schema


TOP = {
    "regulondb_version",
    "release_date",
    "truth_sha256",
    "date",
    "seed",
    "B",
    "cluster",
    "n_directons",
    "uniop",
    "distance_model",
    "terminatornet",
}
UNIOP = {
    "pin",
    "pred_sha256",
    "gff_sha256",
    "n_labelled",
    "n_positive",
    "n_negative",
    "n_without_p",
    "auroc",
    "auroc_ci",
    "auroc_ci_iid",
    "threshold",
    "literal_threshold",
    "rule",
    "precision",
    "precision_ci",
    "precision_ci_iid",
    "recall",
    "recall_ci",
    "recall_ci_iid",
    "f1",
    "tp",
    "fp",
    "precision_cv",
    "threshold_boot_q",
}
DISTANCE = {
    "d0",
    "d0_ci",
    "d0_ci_iid",
    "s",
    "s_ci",
    "s_ci_iid",
    "converged",
    "reason",
    "auroc",
    "auroc_ci",
    "auroc_ci_iid",
    "auroc_diff",
    "auroc_diff_ci",
    "auroc_diff_ci_iid",
    "logloss",
    "logloss_placeholder",
}
TERMINATORNET = {
    "commit",
    "terms_sha256",
    "n_calls",
    "n_eligible",
    "n_truth",
    "n_truth_eligible",
    "threshold",
    "literal_threshold",
    "rule",
    "calibrated",
    "precision",
    "precision_ci",
    "precision_is_lower_bound",
    "coincidence_rate",
    "precision_lower_bound",
    "recall",
    "recall_ci",
    "recall_all",
    "recall_all_ci",
    "recall_by_type",
    "auroc_diagnostic",
    "tp",
    "fp",
}


def test_schema(big: Any) -> None:
    w, _spec, doc, out = big
    on_disk = cutoffs(out)  # strict: NaN / Infinity would fail here
    assert on_disk == doc
    assert set(doc) >= TOP
    assert set(doc["uniop"]) >= UNIOP
    assert set(doc["distance_model"]) >= DISTANCE
    assert set(doc["terminatornet"]) >= TERMINATORNET
    assert doc["regulondb_version"] == "14.5.0"
    assert doc["release_date"] == "2026-01-28"
    summary = json.loads((w.regulondb_dir / "truth" / "truth_summary.json").read_text())
    assert doc["truth_sha256"] == summary["sha256"]
    assert isinstance(doc["date"], str) and doc["date"]
    assert (doc["seed"], doc["B"], doc["cluster"]) == (20261014, FAST["B"], "directon")
    assert doc["terminatornet"]["precision_is_lower_bound"] is True
    assert doc["uniop"]["rule"] in ("precision_0.90", "max_f1")
    for sect in ("uniop", "distance_model", "terminatornet"):
        for k, v in doc[sect].items():
            if k.endswith(("_ci", "_ci_iid")):
                assert isinstance(v, list) and len(v) == 2, k
                assert all(x is None or isinstance(x, float | int) for x in v), k
    q = doc["uniop"]["threshold_boot_q"]
    assert isinstance(q, list) and len(q) == 3
    if None not in q:
        assert q[0] <= q[1] <= q[2]
    cv = doc["uniop"]["precision_cv"]
    assert cv is not None and 0.0 <= cv <= 1.0  # 81 directons >= 10 folds


def test_defaults_are_the_spec_values() -> None:
    import inspect

    params = inspect.signature(cc.run).parameters
    assert params["seed"].default == 20261014
    assert params["B"].default == 2000
    assert params["R"].default == 200
    assert params["min_labelled"].default == 500
    assert params["contig"].default == "U00096.3"


# ================================================================== determinism


_DATE = re.compile(r'"date": "[^"]*"')


def test_same_seed_gives_byte_identical_outputs_apart_from_date(tmp_path: Path) -> None:
    w, _ = big_world(tmp_path / "w")
    a, b = tmp_path / "a", tmp_path / "b"
    run(w, a)
    run(w, b)
    ja, jb = ((d / "cutoffs.json").read_text() for d in (a, b))
    assert _DATE.search(ja)
    assert _DATE.sub("", ja) == _DATE.sub("", jb)
    for name in ("uniop_curve.tsv", "terminatornet_curve.tsv", "terminatornet_calls.tsv"):
        assert (a / name).read_bytes() == (b / name).read_bytes(), name


def test_seed_is_used(tmp_path: Path) -> None:
    w, _ = big_world(tmp_path / "w")
    x = run(w, tmp_path / "a", seed=1)
    y = run(w, tmp_path / "b", seed=2)
    assert (x["seed"], y["seed"]) == (1, 2)
    assert x["uniop"]["auroc"] == y["uniop"]["auroc"]
    assert x["uniop"]["auroc_ci"] != y["uniop"]["auroc_ci"]


def test_default_seed_recorded(tmp_path: Path) -> None:
    w = example6(tmp_path / "w")
    assert run(w, tmp_path / "out", min_labelled=2)["seed"] == 20261014


# ================================================================== CLI


@pytest.fixture
def cli_world(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> World:
    w, _ = big_world(tmp_path / "w")
    monkeypatch.setattr(cc, "K12_PINS", cc.Pins(**w.pins))
    return w


def test_cli_success(cli_world: World, tmp_path: Path) -> None:
    out = tmp_path / "out"
    code = cc.main([*cli_world.cli_args(out), "--B", "20", "--seed", "5"])
    assert code == 0
    for name in (
        "cutoffs.json",
        "uniop_curve.tsv",
        "terminatornet_curve.tsv",
        "terminatornet_calls.tsv",
    ):
        assert (out / name).is_file(), name
    doc = cutoffs(out)
    assert (doc["B"], doc["seed"]) == (20, 5)


def test_cli_identity_mismatch_exits_1(
    cli_world: World,
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    capsys: pytest.CaptureFixture[str],
) -> None:
    monkeypatch.setattr(cc, "K12_PINS", cc.Pins(**{**cli_world.pins, "genome_sha256": "0" * 64}))
    out = tmp_path / "out"
    assert cc.main([*cli_world.cli_args(out), "--B", "20"]) == 1
    assert capsys.readouterr().err.strip()
    assert not (out / "cutoffs.json").exists()


def test_cli_uses_the_real_pins_by_default(tmp_path: Path) -> None:
    """Without the monkeypatch the toy recording is not the K-12 recording: exit 1."""
    w = example6(tmp_path / "w")
    out = tmp_path / "out"
    assert cc.main([*w.cli_args(out), "--B", "20"]) == 1
    assert not (out / "cutoffs.json").exists()


def test_cli_escalation_exits_1(
    cli_world: World, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    (cli_world.uniop_dir / "uniop.pred").write_text(
        "operon prediction scores\nGene A\tGene B\tPrediction\n"
    )
    refresh_pins(cli_world)
    monkeypatch.setattr(cc, "K12_PINS", cc.Pins(**cli_world.pins))
    out = tmp_path / "out"
    assert cc.main([*cli_world.cli_args(out), "--B", "20"]) == 1
    assert not (out / "cutoffs.json").exists()


def test_cli_missing_required_flag_is_a_usage_error(tmp_path: Path) -> None:
    with pytest.raises(SystemExit) as info:
        cc.main(["--out-dir", str(tmp_path / "out")])
    assert info.value.code == 2


def test_cli_decompresses_terms_gz(cli_world: World, tmp_path: Path) -> None:
    """terms.tsv.gz is gzip; the parser reads plain text, so the driver decompresses first."""
    assert cli_world.terminatornet.read_bytes()[:2] == b"\x1f\x8b"
    out = tmp_path / "out"
    assert cc.main([*cli_world.cli_args(out), "--B", "20"]) == 0
    assert cutoffs(out)["terminatornet"]["n_calls"] > 0
