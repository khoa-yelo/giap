"""Module 2 outputs (issue M3-16; GIAP_SPEC.md §7.2.9, §9.2, §9.3, §6.5, §6.6).

``write_mapping_outputs`` writes ``peptide_hits.tsv``, ``protein_hits.tsv``, ``mappings.probed``,
``mappings.bed``, ``mappings.gff3``, ``unmapped.tsv`` and ``mapping_fdr.tsv`` into ``outdir``.

Input contract assumed by these tests (resolved from the spec and the merged producers):

* ``peptide_hits``: the ``PEPTIDE_HITS_COLUMNS`` except ``crosses_origin`` (computed from the
  replicon, §6.5 "never stored"), plus an optional ``blocks`` column (M3-08 ``"s1-e1;s2-e2"``,
  ``None`` for one-block hits); ``start``/``end`` 0-based half-open, unwrapped; ``diff_types`` a
  tuple of strings or ``None``; extra columns are ignored.
* ``protein_hits``: M3-13 ``PROTEIN_HITS_COLUMNS`` rows (0-based half-open).
* ``fdr_rows``: M3-10 ``placement`` rows (``section = placement``) concatenated with de novo
  precision rows (``section = denovo_precision``); columns a row does not use are missing/NaN.
* Section sort is lexicographic (``denovo_precision`` before ``placement``); within a section by
  decoy_type, n_diff, then length-class order (numeric lower bound, ``all`` last).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.mapping_outputs_oracle import PROBED_UNIQUENESS, bed12_core, bed_core_of_line

from giap.io.gff3 import read_gff3
from giap.io.tsv import read_tsv
from giap.mapping.contaminants import UnmappedRow
from giap.mapping.outputs import (
    MAPPING_FDR_COLUMNS,
    PEPTIDE_HITS_COLUMNS,
    PROTEIN_HITS_COLUMNS,
    UNMAPPED_COLUMNS,
    write_mapping_outputs,
)
from giap.model import Interval, Replicon

GENOME_SHA = "ab" * 32
CHR1 = Replicon("chr1", 10_000, "circular", "chromosome")
FS1 = Replicon("fs1", 42, "linear")
REPS: dict[str, Replicon] = {"chr1": CHR1, "fs1": FS1}

OUTPUT_FILES = (
    "peptide_hits.tsv",
    "protein_hits.tsv",
    "mappings.probed",
    "mappings.bed",
    "mappings.gff3",
    "unmapped.tsv",
    "mapping_fdr.tsv",
)
PEP_INPUT_COLUMNS = (*tuple(c for c in PEPTIDE_HITS_COLUMNS if c != "crosses_origin"), "blocks")
PLACEMENT_COLS = (
    "section",
    "decoy_type",
    "n_diff",
    "length_class",
    "targets",
    "decoys",
    "pi_wrong",
    "fdr_stratum",
    "n_accepted",
)
DENOVO_COLS = ("section", "score_bin", "precision", "ci_low", "ci_high", "c_value", "flags")


# ------------------------------------------------------------------------------------ builders


def pep(**kw: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "peptide": "PEPTIDEKLLR",
        "source": "database",
        "best_q_input": 0.001,
        "pep_input": 0.01,
        "score_input": 55.5,
        "contig": "chr1",
        "start": 130,
        "end": 160,
        "strand": "+",
        "frame": 2,
        "route": "exact",
        "n_diff": 0,
        "diff_types": None,
        "placement_q": 0.002,
        "length_class": "10-11",
        "locus_class": "cds_inframe",
        "uniqueness": "unique",
        "gene": "GIAP_00010",
        "novel_candidate_id": None,
        "flags": "",
        "blocks": None,
    }
    row.update(kw)
    row.setdefault("peptide_len", len(row["peptide"]))
    return row


def prot(**kw: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "protein_id": "GIAP_00010",
        "protein_len_aa": 120,
        "route": "miniprot",
        "contig": "chr1",
        "start": 1000,
        "end": 1363,
        "strand": "+",
        "rank": 1,
        "identity": 1.0,
        "query_coverage": 1.0,
        "reciprocal_overlap": 1.0,
        "class": "identical",
        "gene": "GIAP_00010",
        "flags": "",
    }
    row.update(kw)
    return row


def pep_frame(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(list(rows), columns=list(PEP_INPUT_COLUMNS))
    for c in ("peptide_len", "start", "end", "frame", "n_diff"):
        df[c] = df[c].astype("int64")
    return df


def prot_frame(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(list(rows), columns=list(PROTEIN_HITS_COLUMNS))
    for c in ("protein_len_aa", "start", "end", "rank"):
        df[c] = df[c].astype("int64")
    return df


def fdr_frame(
    placement: Sequence[Mapping[str, Any]] = (), denovo: Sequence[Mapping[str, Any]] = ()
) -> pd.DataFrame:
    parts = []
    if placement:
        parts.append(pd.DataFrame([{"section": "placement", **r} for r in placement]))
    if denovo:
        parts.append(pd.DataFrame([{"section": "denovo_precision", **r} for r in denovo]))
    if not parts:
        return pd.DataFrame(columns=list(PLACEMENT_COLS))
    return pd.concat(parts, ignore_index=True)


def placement_row(**kw: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "decoy_type": "reverse",
        "n_diff": 0,
        "length_class": "all",
        "targets": 100,
        "decoys": 2,
        "pi_wrong": 0.5,
        "fdr_stratum": 0.015,
        "n_accepted": 95,
    }
    row.update(kw)
    return row


def denovo_row(**kw: Any) -> dict[str, Any]:
    row: dict[str, Any] = {
        "score_bin": 3,
        "precision": 0.8,
        "ci_low": 0.7,
        "ci_high": 0.9,
        "c_value": 0.99,
        "flags": "",
    }
    row.update(kw)
    return row


def write(
    outdir: Path,
    peps: Sequence[Mapping[str, Any]] = (),
    prots: Sequence[Mapping[str, Any]] = (),
    unmapped: Sequence[UnmappedRow] = (),
    fdr: pd.DataFrame | None = None,
    replicons: Mapping[str, Replicon] = REPS,
) -> Path:
    outdir.mkdir(parents=True, exist_ok=True)
    write_mapping_outputs(
        outdir,
        pep_frame(peps),
        prot_frame(prots),
        list(unmapped),
        fdr_frame() if fdr is None else fdr,
        replicons,
        GENOME_SHA,
    )
    return outdir


# ------------------------------------------------------------------------------------- readers


def raw_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    text = path.read_text(encoding="utf-8")
    assert text.endswith("\n")
    lines = text.split("\n")[:-1]
    header = lines[0].split("\t")
    return header, [dict(zip(header, ln.split("\t"), strict=True)) for ln in lines[1:]]


def data_lines(path: Path) -> list[list[str]]:
    text = path.read_text(encoding="utf-8")
    assert text == "" or text.endswith("\n")
    return [
        ln.split("\t")
        for ln in text.splitlines()
        if ln and not ln.startswith("#") and not ln.startswith("track")
    ]


def gff_features(path: Path, replicons: Mapping[str, Replicon] = REPS) -> list[Any]:
    return read_gff3(path, replicons)


def attr(f: Any, key: str) -> str | None:
    v = f.attributes.get(key)
    if v is None or isinstance(v, str):
        return v
    return ",".join(v)


# ----------------------------------------------------------------------------- column tuples


def test_peptide_hits_columns_pinned() -> None:
    assert PEPTIDE_HITS_COLUMNS == (
        "peptide",
        "peptide_len",
        "source",
        "best_q_input",
        "pep_input",
        "score_input",
        "contig",
        "start",
        "end",
        "strand",
        "frame",
        "crosses_origin",
        "route",
        "n_diff",
        "diff_types",
        "placement_q",
        "length_class",
        "locus_class",
        "uniqueness",
        "gene",
        "novel_candidate_id",
        "flags",
    )


def test_protein_hits_columns_pinned() -> None:
    assert PROTEIN_HITS_COLUMNS == (
        "protein_id",
        "protein_len_aa",
        "route",
        "contig",
        "start",
        "end",
        "strand",
        "rank",
        "identity",
        "query_coverage",
        "reciprocal_overlap",
        "class",
        "gene",
        "flags",
    )


def test_unmapped_columns_pinned() -> None:
    assert UNMAPPED_COLUMNS == (
        "sequence_id",
        "sequence",
        "kind",
        "label",
        "reason",
        "contaminant_match",
    )


def test_mapping_fdr_columns_pinned() -> None:
    assert MAPPING_FDR_COLUMNS == (
        "length_class",
        "n_diff",
        "targets",
        "decoys",
        "pi_wrong",
        "fdr_stratum",
        "n_accepted",
        "decoy_type",
        "section",
        "score_bin",
        "precision",
        "ci_low",
        "ci_high",
        "c_value",
        "flags",
    )


def test_writes_exactly_the_seven_files(tmp_path: Path) -> None:
    out = write(tmp_path / "mapping", [pep()], [prot()])
    assert sorted(p.name for p in out.iterdir()) == sorted(OUTPUT_FILES)


def test_tsv_headers_in_order(tmp_path: Path) -> None:
    out = write(tmp_path, [pep()], [prot()], fdr=fdr_frame([placement_row()]))
    for name, cols in (
        ("peptide_hits.tsv", PEPTIDE_HITS_COLUMNS),
        ("protein_hits.tsv", PROTEIN_HITS_COLUMNS),
        ("unmapped.tsv", UNMAPPED_COLUMNS),
        ("mapping_fdr.tsv", MAPPING_FDR_COLUMNS),
    ):
        header, _ = raw_tsv(out / name)
        assert tuple(header) == cols, name


def test_extra_input_columns_not_written(tmp_path: Path) -> None:
    df = pep_frame([pep()])
    df["query"] = "PEPTLDEKLLR"
    df["e_target"] = 0.5
    write_mapping_outputs(tmp_path, df, prot_frame([]), [], fdr_frame(), REPS, GENOME_SHA)
    header, _ = raw_tsv(tmp_path / "peptide_hits.tsv")
    assert tuple(header) == PEPTIDE_HITS_COLUMNS


# ---------------------------------------------------------------------------- Example 1 (TSV)


def test_example_01_exact_hit_tsv_row(tmp_path: Path) -> None:
    out = write(tmp_path, [pep()])
    _, rows = raw_tsv(out / "peptide_hits.tsv")
    (r,) = rows
    assert r["peptide"] == "PEPTIDEKLLR"
    assert r["contig"] == "chr1"
    assert (r["start"], r["end"], r["strand"]) == ("131", "160", "+")
    assert r["frame"] == "2"
    assert r["crosses_origin"] == "false"
    assert r["route"] == "exact"
    assert r["n_diff"] == "0"
    assert r["diff_types"] == "NA"
    assert r["placement_q"] == "0.002"


def test_example_01_scalar_dtypes_written(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(score_input=55.5, best_q_input=0.001, pep_input=0.01)])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert r["peptide_len"] == "11"
    assert r["source"] == "database"
    assert (r["best_q_input"], r["pep_input"], r["score_input"]) == ("0.001", "0.01", "55.5")
    assert r["length_class"] == "10-11"
    assert r["locus_class"] == "cds_inframe"
    assert r["uniqueness"] == "unique"
    assert r["gene"] == "GIAP_00010"
    assert r["novel_candidate_id"] == "NA"


def test_example_01_bed_line(tmp_path: Path) -> None:
    out = write(tmp_path, [pep()])
    (f,) = data_lines(out / "mappings.bed")
    assert len(f) == 12
    assert bed_core_of_line(f) == bed12_core("chr1", 130, 160, "+", "PEPTIDEKLLR", 1000, 10_000)[0]


def test_minus_strand_hit(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(start=301, end=331, strand="-", frame=-3)])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert (r["start"], r["end"], r["strand"], r["frame"]) == ("302", "331", "-", "-3")
    (b,) = data_lines(out / "mappings.bed")
    assert (b[1], b[2], b[5]) == ("301", "331", "-")
    (p,) = data_lines(out / "mappings.probed")
    assert p[5] == "-"
    (g,) = gff_features(out / "mappings.gff3")
    assert g.interval == Interval("chr1", 301, 331, "-")


# ------------------------------------------------------------------------- Example 2 (origin)

ORIGIN_PEP = "PEPTIDEKLLR"


def test_example_02_origin_crossing_tsv(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(start=9990, end=10020, frame=1)])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert (r["start"], r["end"], r["crosses_origin"]) == ("9991", "10020", "true")


def test_example_02_origin_crossing_bed_split(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(start=9990, end=10020, frame=1)])
    lines = data_lines(out / "mappings.bed")
    assert [tuple(f[0:4]) for f in lines] == [
        ("chr1", "0", "20", f"{ORIGIN_PEP}_part2"),  # sorted by start: part 2 first
        ("chr1", "9990", "10000", f"{ORIGIN_PEP}_part1"),
    ]
    assert all(len(f) == 12 for f in lines)
    assert [f[9] for f in lines] == ["1", "1"]


def test_example_02_origin_crossing_probed_split(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(start=9990, end=10020, frame=1)])
    lines = data_lines(out / "mappings.probed")
    assert sorted((f[1], f[2], f[3]) for f in lines) == [
        ("0", "20", f"{ORIGIN_PEP}_part2"),
        ("9990", "10000", f"{ORIGIN_PEP}_part1"),
    ]
    assert all(len(f) == 25 for f in lines)
    assert all(f[13] == ORIGIN_PEP for f in lines)


def test_example_02_origin_crossing_gff3_keeps_end_beyond_length(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(start=9990, end=10020, frame=1)])
    feats = [f for f in gff_features(out / "mappings.gff3") if f.type != "region"]
    (f,) = feats
    assert f.interval == Interval("chr1", 9990, 10020, "+")
    raw = [
        ln.split("\t")
        for ln in (out / "mappings.gff3").read_text().splitlines()
        if ln and not ln.startswith("#") and "\tnucleotide_to_protein_match\t" in ln
    ]
    assert [(c[3], c[4]) for c in raw] == [("9991", "10020")]


def test_origin_crossing_protein_hit_split_in_bed(tmp_path: Path) -> None:
    out = write(tmp_path, prots=[prot(start=9900, end=10050)])
    _, (r,) = raw_tsv(out / "protein_hits.tsv")
    assert (r["start"], r["end"]) == ("9901", "10050")
    names = sorted(f[3] for f in data_lines(out / "mappings.bed"))
    assert names == ["GIAP_00010_part1", "GIAP_00010_part2"]


# --------------------------------------------------------------------- Example 3 (tolerant)


def test_example_03_tolerant_diff_types(tmp_path: Path) -> None:
    row = pep(route="tolerant", n_diff=2, diff_types=("substitution", "swap"), placement_q=0.004)
    out = write(tmp_path, [row])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert r["route"] == "tolerant"
    assert r["n_diff"] == "2"
    assert r["diff_types"] == "substitution;swap"


def test_example_03_tolerant_bed_score_500(tmp_path: Path) -> None:
    row = pep(route="tolerant", n_diff=2, diff_types=("substitution", "swap"))
    out = write(tmp_path, [row])
    (b,) = data_lines(out / "mappings.bed")
    assert b[4] == "500"
    (p,) = data_lines(out / "mappings.probed")
    assert p[4] == "500"


def test_example_03_tolerant_gff3_attributes(tmp_path: Path) -> None:
    row = pep(route="tolerant", n_diff=2, diff_types=("substitution", "swap"), placement_q=0.004)
    out = write(tmp_path, [row])
    (f,) = gff_features(out / "mappings.gff3")
    assert attr(f, "route") == "tolerant"
    assert attr(f, "n_diff") == "2"
    dt = f.attributes["diff_types"]
    items = [dt] if isinstance(dt, str) else list(dt)
    assert [x for s in items for x in s.replace(";", ",").split(",")] == ["substitution", "swap"]


# -------------------------------------------------------------------- Example 4 (frameshift)

FS_ROW = {
    "peptide": "PEPTIDEKLMNA",
    "contig": "fs1",
    "start": 2,
    "end": 39,
    "strand": "+",
    "frame": 3,
    "route": "tolerant",
    "n_diff": 1,
    "diff_types": ("frameshift",),
    "blocks": "2-20;21-39",
}


def test_example_04_frameshift_bed12_blocks(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(**FS_ROW)])
    (f,) = data_lines(out / "mappings.bed")
    assert len(f) == 12
    assert (f[0], f[1], f[2], f[3]) == ("fs1", "2", "39", "PEPTIDEKLMNA")
    assert f[9] == "2"
    assert f[10].rstrip(",") == "18,18"
    assert f[11].rstrip(",") == "0,19"
    assert f[4] == "500"


def test_example_04_frameshift_tsv_row(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(**FS_ROW)])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert (r["start"], r["end"], r["frame"]) == ("3", "39", "3")
    assert r["diff_types"] == "frameshift"
    assert "blocks" not in r


def test_example_04_frameshift_probed_blocks(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(**FS_ROW)])
    (f,) = data_lines(out / "mappings.probed")
    assert len(f) == 25
    assert (f[9], f[10].rstrip(","), f[11].rstrip(",")) == ("2", "18,18", "0,19")


def test_example_04_frameshift_gff3_match_parts(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(**FS_ROW)])
    feats = gff_features(out / "mappings.gff3")
    (parent,) = [f for f in feats if f.type == "nucleotide_to_protein_match"]
    parts = [f for f in feats if f.type == "match_part"]
    assert parent.interval == Interval("fs1", 2, 39, "+")
    assert len(parts) == 2
    pid = attr(parent, "ID")
    assert pid
    assert all(attr(p, "Parent") == pid for p in parts)
    assert sorted((p.interval.start, p.interval.end) for p in parts) == [(2, 20), (21, 39)]
    assert all(p.source == "GIAP" for p in parts)


def test_frameshift_delta_minus_one_first_block_shortened(tmp_path: Path) -> None:
    # δ = -1: second block starts 1 nt before the first ends (M3-08 example 6).
    row = pep(**{**FS_ROW, "end": 37, "blocks": "2-20;19-37"})
    out = write(tmp_path, [row])
    (f,) = data_lines(out / "mappings.bed")
    assert (f[1], f[2], f[9]) == ("2", "37", "2")
    assert f[10].rstrip(",") == "17,18"
    assert f[11].rstrip(",") == "0,17"
    assert (
        bed_core_of_line(f)
        == bed12_core("fs1", 2, 37, "+", "PEPTIDEKLMNA", 500, 42, "2-20;19-37")[0]
    )


def test_frameshift_minus_strand_blocks(tmp_path: Path) -> None:
    # M3-08 example 5: the reverse complement covers [3, 40) with blocks 3-21;22-40.
    row = pep(**{**FS_ROW, "start": 3, "end": 40, "strand": "-", "frame": -3,
                 "blocks": "3-21;22-40"})  # fmt: skip
    out = write(tmp_path, [row])
    (f,) = data_lines(out / "mappings.bed")
    assert (f[1], f[2], f[5], f[9]) == ("3", "40", "-", "2")
    assert f[10].rstrip(",") == "18,18"
    assert f[11].rstrip(",") == "0,19"


# ------------------------------------------------------------------------ Example 5 (%.6g)


def test_example_05_placement_q_six_significant_digits(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(placement_q=0.000123456789)])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert r["placement_q"] == "0.000123457"
    (p,) = data_lines(out / "mappings.probed")
    assert p[17] == "0.000123457"  # fdr


# ------------------------------------------------------------------- Examples 6-7 (fdr rows)


def test_example_06_per_m_total_row(tmp_path: Path) -> None:
    fdr = fdr_frame([placement_row(length_class="all", n_diff=1, targets=40, decoys=1)])
    out = write(tmp_path, fdr=fdr)
    _, (r,) = raw_tsv(out / "mapping_fdr.tsv")
    assert r["length_class"] == "all"
    assert r["section"] == "placement"
    assert (r["n_diff"], r["targets"], r["decoys"], r["n_accepted"]) == ("1", "40", "1", "95")
    assert (r["pi_wrong"], r["fdr_stratum"], r["decoy_type"]) == ("0.5", "0.015", "reverse")
    for c in ("score_bin", "precision", "ci_low", "ci_high", "c_value"):
        assert r[c] == "NA", c


def test_example_07_denovo_bin_row(tmp_path: Path) -> None:
    out = write(tmp_path, fdr=fdr_frame(denovo=[denovo_row(score_bin=3)]))
    _, (r,) = raw_tsv(out / "mapping_fdr.tsv")
    assert r["section"] == "denovo_precision"
    assert r["score_bin"] == "3"
    assert (r["precision"], r["ci_low"], r["ci_high"], r["c_value"]) == (
        "0.8",
        "0.7",
        "0.9",
        "0.99",
    )
    for c in ("length_class", "n_diff", "targets", "decoys", "pi_wrong", "fdr_stratum",
              "n_accepted"):  # fmt: skip
        assert r[c] == "NA", c


def test_denovo_unreported_bin_flags(tmp_path: Path) -> None:
    row = denovo_row(score_bin=0, precision=None, ci_low=None, ci_high=None, flags="unreported")
    out = write(tmp_path, fdr=fdr_frame(denovo=[row]))
    _, (r,) = raw_tsv(out / "mapping_fdr.tsv")
    assert (r["precision"], r["ci_low"], r["ci_high"]) == ("NA", "NA", "NA")
    assert r["flags"] == "unreported"


def test_mapping_fdr_sorted(tmp_path: Path) -> None:
    labels = ["7-9", "10-11", "12-30", "all"]
    placement = [
        placement_row(decoy_type=dt, n_diff=m, length_class=lc)
        for dt in ("shuffle", "reverse")
        for m in (2, 0, 1)
        for lc in ("all", "12-30", "7-9", "10-11")
    ]
    denovo = [denovo_row(score_bin=b) for b in (0, 1, 2)]
    out = write(tmp_path, fdr=fdr_frame(placement, denovo))
    _, rows = raw_tsv(out / "mapping_fdr.tsv")
    assert [r["section"] for r in rows] == ["denovo_precision"] * 3 + ["placement"] * 24
    got = [(r["decoy_type"], r["n_diff"], r["length_class"]) for r in rows[3:]]
    want = [(dt, str(m), lc) for dt in ("reverse", "shuffle") for m in (0, 1, 2) for lc in labels]
    assert got == want
    assert sorted(r["score_bin"] for r in rows[:3]) == ["0", "1", "2"]


# --------------------------------------------------------------------- Example 8 (empty)


def test_example_08_empty_inputs(tmp_path: Path) -> None:
    out = write(tmp_path)
    for name, cols in (
        ("peptide_hits.tsv", PEPTIDE_HITS_COLUMNS),
        ("protein_hits.tsv", PROTEIN_HITS_COLUMNS),
        ("unmapped.tsv", UNMAPPED_COLUMNS),
        ("mapping_fdr.tsv", MAPPING_FDR_COLUMNS),
    ):
        assert (out / name).read_text(encoding="utf-8") == "\t".join(cols) + "\n", name
    assert (out / "mappings.probed").read_text(encoding="utf-8") == "# proBed-version 1.0\n"
    bed = (out / "mappings.bed").read_text(encoding="utf-8").splitlines()
    assert len(bed) <= 1
    assert all(ln.startswith(("#", "track")) for ln in bed)
    assert [f for f in gff_features(out / "mappings.gff3") if f.type != "region"] == []
    assert (out / "mappings.gff3").read_text(encoding="utf-8").startswith("##gff-version 3\n")


# ------------------------------------------------------------------------------- proBed


def test_probed_fields(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(score_input=55.5, placement_q=0.002)])
    lines = (out / "mappings.probed").read_text(encoding="utf-8").splitlines()
    assert lines[0] == "# proBed-version 1.0"
    (f,) = data_lines(out / "mappings.probed")
    assert len(f) == 25
    assert f[0:6] == ["chr1", "130", "160", "PEPTIDEKLLR", "1000", "+"]
    assert f[12] == "GIAP_00010"  # proteinAccession = gene
    assert f[13] == "PEPTIDEKLLR"  # peptideSequence
    assert f[14] == "unique"
    assert f[15] == f"sha256:{GENOME_SHA}"  # genomeReferenceVersion (Q2)
    assert f[16] == "55.5"  # psmScore = score_input
    assert f[17] == "0.002"  # fdr = placement_q (lead decision Q2)


@pytest.mark.parametrize(("uniq", "written"), sorted(PROBED_UNIQUENESS.items()))
def test_probed_uniqueness_vocabulary(tmp_path: Path, uniq: str, written: str) -> None:
    out = write(tmp_path, [pep(uniqueness=uniq, gene=None if uniq == "shared" else "G1")])
    (f,) = data_lines(out / "mappings.probed")
    assert f[14] == written


def test_probed_has_peptides_only(tmp_path: Path) -> None:
    out = write(tmp_path, [pep()], [prot(protein_id="PROT_A", start=2000, end=2300)])
    lines = data_lines(out / "mappings.probed")
    assert [f[3] for f in lines] == ["PEPTIDEKLLR"]
    assert all(len(f) == 25 for f in lines)


def test_probed_every_line_25_fields(tmp_path: Path) -> None:
    peps = [
        pep(),
        pep(peptide="AAAAKLLR", start=9990, end=10014, frame=1),
        pep(**FS_ROW),
        pep(peptide="MSTVKR", start=500, end=518, strand="-", frame=-1, gene=None),
    ]
    out = write(tmp_path, peps)
    lines = data_lines(out / "mappings.probed")
    assert len(lines) == 5  # one origin-crossing peptide split in two
    assert all(len(f) == 25 for f in lines)


# --------------------------------------------------------------------------------- BED


def test_bed_peptides_and_proteins(tmp_path: Path) -> None:
    prots = [
        prot(protein_id="P_ID", start=2000, end=2300, **{"class": "identical"}),
        prot(protein_id="P_NEAR", start=3000, end=3300, identity=0.95, **{"class": "near"}),
        prot(protein_id="P_DIV", start=4000, end=4300, identity=0.6, **{"class": "divergent"}),
    ]
    out = write(tmp_path, [pep()], prots)
    lines = data_lines(out / "mappings.bed")
    assert all(len(f) == 12 for f in lines)
    by_name = {f[3]: f for f in lines}
    assert set(by_name) == {"PEPTIDEKLLR", "P_ID", "P_NEAR", "P_DIV"}
    assert by_name["PEPTIDEKLLR"][4] == "1000"
    assert by_name["P_ID"][4] == "1000"
    assert by_name["P_NEAR"][4] == "500"
    assert by_name["P_DIV"][4] == "500"
    assert (
        bed_core_of_line(by_name["P_ID"])
        == bed12_core("chr1", 2000, 2300, "+", "P_ID", 1000, 10_000)[0]
    )


def test_bed_sorted_by_position(tmp_path: Path) -> None:
    out = write(
        tmp_path,
        [pep(peptide="ZZZK", start=5000, end=5012), pep(peptide="AAAK", start=100, end=112)],
        [prot(protein_id="P1", start=3000, end=3300)],
    )
    starts = [int(f[1]) for f in data_lines(out / "mappings.bed")]
    assert starts == sorted(starts)


# --------------------------------------------------------------------------------- GFF3


def test_gff3_peptide_and_protein_features(tmp_path: Path) -> None:
    out = write(tmp_path, [pep()], [prot(protein_id="PROT_A", start=2000, end=2300, rank=2)])
    feats = [f for f in gff_features(out / "mappings.gff3") if f.type != "region"]
    assert len(feats) == 2
    assert all(f.source == "GIAP" and f.type == "nucleotide_to_protein_match" for f in feats)
    pf, rf = sorted(feats, key=lambda f: f.interval.start)
    assert pf.interval == Interval("chr1", 130, 160, "+")
    assert attr(pf, "Name") == "PEPTIDEKLLR"
    assert attr(pf, "route") == "exact"
    assert attr(pf, "n_diff") == "0"
    assert attr(pf, "placement_q") == "0.002"
    assert attr(pf, "locus_class") == "cds_inframe"
    assert attr(pf, "uniqueness") == "unique"
    assert attr(pf, "gene") == "GIAP_00010"
    assert rf.interval == Interval("chr1", 2000, 2300, "+")
    assert attr(rf, "Name") == "PROT_A"
    assert attr(rf, "identity") == "1"
    assert attr(rf, "query_coverage") == "1"
    assert attr(rf, "class") == "identical"
    assert attr(rf, "rank") == "2"


def test_gff3_ids_unique_for_multi_locus_peptide(tmp_path: Path) -> None:
    peps = [
        pep(uniqueness="multi_identical", gene="G1"),
        pep(uniqueness="multi_identical", gene="G2", start=4000, end=4030, frame=2),
        pep(uniqueness="multi_identical", gene="G3", start=4000, end=4030, strand="-", frame=-1),
    ]
    out = write(tmp_path, peps)
    feats = [f for f in gff_features(out / "mappings.gff3") if f.type != "region"]
    ids = [attr(f, "ID") for f in feats]
    assert len(ids) == 3 and all(ids) and len(set(ids)) == 3
    assert {attr(f, "Name") for f in feats} == {"PEPTIDEKLLR"}


def test_gff3_has_sequence_region_pragmas(tmp_path: Path) -> None:
    out = write(tmp_path, [pep()])
    text = (out / "mappings.gff3").read_text(encoding="utf-8")
    assert "##sequence-region chr1 1 10000" in text
    assert "##sequence-region fs1 1 42" in text


# ------------------------------------------------------------------------------ edge cases


def test_gene_na_written_per_format(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(gene=None, locus_class="intergenic")])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert r["gene"] == "NA"
    (p,) = data_lines(out / "mappings.probed")
    assert p[12] == "."
    (f,) = gff_features(out / "mappings.gff3")
    assert "gene" not in f.attributes


def test_isoleucine_kept_as_written(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(peptide="PEPTIDEKIIR")])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert r["peptide"] == "PEPTIDEKIIR"
    (b,) = data_lines(out / "mappings.bed")
    assert b[3] == "PEPTIDEKIIR"
    (p,) = data_lines(out / "mappings.probed")
    assert (p[3], p[13]) == ("PEPTIDEKIIR", "PEPTIDEKIIR")
    (f,) = gff_features(out / "mappings.gff3")
    assert attr(f, "Name") == "PEPTIDEKIIR"


def test_flags_written(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(flags="low_complexity;same_db")], [prot(flags="fs;st")])
    _, (r,) = raw_tsv(out / "peptide_hits.tsv")
    assert r["flags"] == "low_complexity;same_db"
    _, (q,) = raw_tsv(out / "protein_hits.tsv")
    assert q["flags"] == "fs;st"


def test_empty_flags_read_as_missing(tmp_path: Path) -> None:
    out = write(tmp_path, [pep(flags="")])
    header, rows = read_tsv(out / "peptide_hits.tsv")
    assert rows[0][header.index("flags")] is None


def test_protein_hits_tsv_row(tmp_path: Path) -> None:
    row = prot(route="tblastn", identity=0.97, query_coverage=0.875, reciprocal_overlap=None,
               gene=None, rank=2, strand="-", **{"class": "near"})  # fmt: skip
    out = write(tmp_path, prots=[row])
    _, (r,) = raw_tsv(out / "protein_hits.tsv")
    assert (r["protein_id"], r["protein_len_aa"], r["route"]) == ("GIAP_00010", "120", "tblastn")
    assert (r["contig"], r["start"], r["end"], r["strand"], r["rank"]) == (
        "chr1", "1001", "1363", "-", "2",
    )  # fmt: skip
    assert (r["identity"], r["query_coverage"], r["reciprocal_overlap"]) == ("0.97", "0.875", "NA")
    assert (r["class"], r["gene"]) == ("near", "NA")


def test_unmapped_rows_and_sort(tmp_path: Path) -> None:
    rows = [
        UnmappedRow("pep_b", "KLLR", "peptide", "unmapped", "no_hit", None),
        UnmappedRow("PROT_X", "MSTV", "protein", "contaminant", "contaminant_id", "CON_P1"),
        UnmappedRow("pep_a", "AAAK", "peptide", "contaminant", "no_hit", "CON_TRYP"),
        UnmappedRow("PROT_A", "MKKL", "protein", "unmapped", "below_identity", None),
    ]
    out = write(tmp_path, unmapped=rows)
    _, got = raw_tsv(out / "unmapped.tsv")
    assert [(r["kind"], r["sequence_id"]) for r in got] == [
        ("peptide", "pep_a"),
        ("peptide", "pep_b"),
        ("protein", "PROT_A"),
        ("protein", "PROT_X"),
    ]
    assert got[1] == {
        "sequence_id": "pep_b",
        "sequence": "KLLR",
        "kind": "peptide",
        "label": "unmapped",
        "reason": "no_hit",
        "contaminant_match": "NA",
    }
    assert got[0]["contaminant_match"] == "CON_TRYP"


def test_peptide_hits_sorted(tmp_path: Path) -> None:
    reps = {"zeta": Replicon("zeta", 5000, "circular"), "alpha": Replicon("alpha", 5000, "linear")}
    rows = [
        pep(peptide="BBBK", contig="alpha", start=10, end=22),
        pep(peptide="BBBK", contig="zeta", start=100, end=112, strand="-", frame=-1),
        pep(peptide="CCCK", contig="zeta", start=100, end=112),
        pep(peptide="AAAK", contig="zeta", start=100, end=112),
        pep(peptide="DDDDK", contig="zeta", start=100, end=115),
        pep(peptide="EEEK", contig="zeta", start=40, end=52),
    ]
    out = write(tmp_path, rows, replicons=reps)
    _, got = raw_tsv(out / "peptide_hits.tsv")
    assert [(r["contig"], r["start"], r["end"], r["strand"], r["peptide"]) for r in got] == [
        ("zeta", "41", "52", "+", "EEEK"),
        ("zeta", "101", "112", "+", "AAAK"),
        ("zeta", "101", "112", "+", "CCCK"),
        ("zeta", "101", "112", "-", "BBBK"),
        ("zeta", "101", "115", "+", "DDDDK"),
        ("alpha", "11", "22", "+", "BBBK"),
    ]


def test_protein_hits_sorted(tmp_path: Path) -> None:
    rows = [
        prot(protein_id="P2", start=500, end=800),
        prot(protein_id="P1", start=500, end=800),
        prot(protein_id="P3", start=50, end=350),
        prot(protein_id="P0", contig="fs1", start=0, end=36),
    ]
    out = write(tmp_path, prots=rows)
    _, got = raw_tsv(out / "protein_hits.tsv")
    assert [r["protein_id"] for r in got] == ["P3", "P1", "P2", "P0"]  # chr1 before fs1


def test_byte_identical_reruns_and_no_paths(tmp_path: Path) -> None:
    peps = [pep(), pep(peptide="AAAAKLLR", start=9990, end=10014, frame=1), pep(**FS_ROW)]
    prots = [prot(), prot(protein_id="P2", start=9950, end=10100, **{"class": "near"})]
    unm = [UnmappedRow("pep_x", "KLLR", "peptide", "unmapped", "too_short", None)]
    fdr = fdr_frame([placement_row()], [denovo_row()])
    a = write(tmp_path / "run_a" / "mapping", peps, prots, unm, fdr)
    b = write(tmp_path / "run_b" / "mapping", list(reversed(peps)), list(reversed(prots)), unm, fdr)
    for name in OUTPUT_FILES:
        ta, tb = (a / name).read_bytes(), (b / name).read_bytes()
        assert ta == tb, name
        assert str(tmp_path).encode() not in ta, name
        assert b"\r" not in ta, name
