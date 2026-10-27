"""Hypothesis properties of the Module 2 output writers (issue M3-16; GIPPA_SPEC.md §6.5, §9.2).

Checked against the independent oracle in tests/helpers/mapping_outputs_oracle.py:

* TSV coordinates are the 0-based input shifted to 1-based closed, ``end > L`` kept, and
  ``crosses_origin`` equals ``end > L``;
* BED12 core fields equal the oracle's (origin-crossing split into ``_part1``/``_part2``), sorted;
* proBed has one 25-field line per peptide part, peptides only;
* GFF3 read back with M1-05's reader returns the input intervals;
* rows are sorted by (contig FASTA order, start, end, strand, ID) and the files do not depend on
  the input row order (byte-identical).
"""

from __future__ import annotations

import tempfile
from collections import Counter
from pathlib import Path
from typing import Any

import pandas as pd
from helpers.mapping_outputs_oracle import (
    STRAND_RANK,
    bed12_core,
    bed_core_of_line,
    bed_sort_key,
    fmt6g,
    hit_sort_key,
    tsv_coords,
)
from hypothesis import given, settings
from hypothesis import strategies as st

from gippa.io.gff3 import read_gff3
from gippa.mapping.outputs import PEPTIDE_HITS_COLUMNS, PROTEIN_HITS_COLUMNS, write_mapping_outputs
from gippa.model import Replicon

PROPS = settings(max_examples=60, deadline=None)
SHA = "0f" * 32
FILES = (
    "peptide_hits.tsv",
    "protein_hits.tsv",
    "mappings.probed",
    "mappings.bed",
    "mappings.gff3",
    "unmapped.tsv",
    "mapping_fdr.tsv",
)
AA = "ACDEFGHIKLMNPQRSTVWY"


@st.composite
def replicons(draw: st.DrawFn) -> dict[str, Replicon]:
    n = draw(st.integers(1, 3))
    names = draw(st.lists(st.from_regex(r"[a-z][a-z0-9]{0,5}", fullmatch=True), min_size=n,
                          max_size=n, unique=True))  # fmt: skip
    return {
        name: Replicon(
            name, draw(st.integers(60, 600)), draw(st.sampled_from(["circular", "linear"]))
        )
        for name in names
    }


@st.composite
def interval(draw: st.DrawFn, reps: dict[str, Replicon]) -> tuple[str, int, int, str]:
    rep = reps[draw(st.sampled_from(sorted(reps)))]
    start = draw(st.integers(0, rep.length - 1))
    max_end = start + min(rep.length, 90) if rep.topology == "circular" else rep.length
    end = draw(st.integers(start + 1, max(start + 1, min(max_end, start + 90))))
    if rep.topology == "linear":
        end = min(end, rep.length)
    return rep.id, start, end, draw(st.sampled_from(["+", "-"]))


@st.composite
def dataset(draw: st.DrawFn) -> dict[str, Any]:
    reps = draw(replicons())
    n_pep = draw(st.integers(0, 8))
    peps: list[dict[str, Any]] = []
    seen: set[tuple[str, str, int, int, str]] = set()
    for _ in range(n_pep):
        contig, s, e, strand = draw(interval(reps))
        peptide = draw(st.text(alphabet=AA, min_size=4, max_size=14))
        key = (peptide, contig, s, e, strand)
        if key in seen:
            continue
        seen.add(key)
        exact = draw(st.booleans())
        frame = (s % 3) + 1 if strand == "+" else -(((reps[contig].length - e) % 3) + 1)
        peps.append(
            {
                "peptide": peptide,
                "peptide_len": len(peptide),
                "source": draw(st.sampled_from(["database", "denovo"])),
                "best_q_input": draw(st.none() | st.floats(0, 0.01)),
                "pep_input": draw(st.none() | st.floats(0, 1)),
                "score_input": draw(st.none() | st.floats(0, 300)),
                "contig": contig,
                "start": s,
                "end": e,
                "strand": strand,
                "frame": frame,
                "route": "exact" if exact else "tolerant",
                "n_diff": 0 if exact else 1,
                "diff_types": None if exact else ("substitution",),
                "placement_q": draw(st.floats(0, 0.01)),
                "length_class": "all",
                "locus_class": "cds_inframe",
                "uniqueness": draw(st.sampled_from(["unique", "multi_identical", "shared"])),
                "gene": draw(st.none() | st.sampled_from(["G_00001", "G_00002"])),
                "novel_candidate_id": None,
                "flags": "",
                "blocks": None,
            }
        )
    prots: list[dict[str, Any]] = []
    seen_p: set[tuple[str, str, int, int, str]] = set()
    for _ in range(draw(st.integers(0, 5))):
        contig, s, e, strand = draw(interval(reps))
        pid = draw(st.from_regex(r"P[0-9]{1,3}", fullmatch=True))
        if (pid, contig, s, e, strand) in seen_p:
            continue
        seen_p.add((pid, contig, s, e, strand))
        prots.append(
            {
                "protein_id": pid,
                "protein_len_aa": 50,
                "route": draw(st.sampled_from(["miniprot", "tblastn"])),
                "contig": contig,
                "start": s,
                "end": e,
                "strand": strand,
                "rank": draw(st.integers(1, 3)),
                "identity": 0.9,
                "query_coverage": 0.9,
                "reciprocal_overlap": 0.5,
                "class": draw(st.sampled_from(["identical", "near", "divergent"])),
                "gene": None,
                "flags": "",
            }
        )
    return {"reps": reps, "peps": peps, "prots": prots}


def _pep_df(rows: list[dict[str, Any]]) -> pd.DataFrame:
    cols = [c for c in PEPTIDE_HITS_COLUMNS if c != "crosses_origin"] + ["blocks"]
    df = pd.DataFrame(rows, columns=cols)
    for c in ("peptide_len", "start", "end", "frame", "n_diff"):
        df[c] = df[c].astype("int64")
    return df


def _prot_df(rows: list[dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(rows, columns=list(PROTEIN_HITS_COLUMNS))
    for c in ("protein_len_aa", "start", "end", "rank"):
        df[c] = df[c].astype("int64")
    return df


def _write(
    out: Path, ds: dict[str, Any], peps: list[dict[str, Any]], prots: list[dict[str, Any]]
) -> None:
    out.mkdir(parents=True, exist_ok=True)
    write_mapping_outputs(
        out, _pep_df(peps), _prot_df(prots), [], pd.DataFrame(columns=["section"]), ds["reps"], SHA
    )


def _tsv(path: Path) -> list[dict[str, str]]:
    lines = path.read_text(encoding="utf-8").split("\n")[:-1]
    header = lines[0].split("\t")
    return [dict(zip(header, ln.split("\t"), strict=True)) for ln in lines[1:]]


def _data(path: Path) -> list[list[str]]:
    return [
        ln.split("\t")
        for ln in path.read_text(encoding="utf-8").splitlines()
        if ln and not ln.startswith(("#", "track"))
    ]


@PROPS
@given(dataset())
def test_tsv_coordinates_and_sort(ds: dict[str, Any]) -> None:
    order = {cid: i for i, cid in enumerate(ds["reps"])}
    with tempfile.TemporaryDirectory() as d:
        out = Path(d)
        _write(out, ds, ds["peps"], ds["prots"])
        got = _tsv(out / "peptide_hits.tsv")
        want = sorted(ds["peps"], key=lambda r: hit_sort_key(r, order, "peptide"))
        assert len(got) == len(want)
        for g, w in zip(got, want, strict=True):
            L = ds["reps"][w["contig"]].length
            assert (g["peptide"], g["contig"], g["strand"]) == (
                w["peptide"],
                w["contig"],
                w["strand"],
            )
            assert (g["start"], g["end"], g["crosses_origin"]) == tsv_coords(
                w["start"], w["end"], L
            )
            assert g["placement_q"] == fmt6g(w["placement_q"])
        gotp = _tsv(out / "protein_hits.tsv")
        wantp = sorted(ds["prots"], key=lambda r: hit_sort_key(r, order, "protein_id"))
        assert [(g["protein_id"], g["start"], g["end"]) for g in gotp] == [
            (w["protein_id"], str(w["start"] + 1), str(w["end"])) for w in wantp
        ]


@PROPS
@given(dataset())
def test_bed_matches_oracle(ds: dict[str, Any]) -> None:
    reps = ds["reps"]
    order = {cid: i for i, cid in enumerate(reps)}
    want = []
    for r in ds["peps"]:
        score = 1000 if r["route"] == "exact" else 500
        want += bed12_core(r["contig"], r["start"], r["end"], r["strand"], r["peptide"], score,
                           reps[r["contig"]].length)  # fmt: skip
    for r in ds["prots"]:
        score = 1000 if r["class"] == "identical" else 500
        want += bed12_core(r["contig"], r["start"], r["end"], r["strand"], r["protein_id"], score,
                           reps[r["contig"]].length)  # fmt: skip
    with tempfile.TemporaryDirectory() as d:
        out = Path(d)
        _write(out, ds, ds["peps"], ds["prots"])
        lines = _data(out / "mappings.bed")
        assert all(len(f) == 12 for f in lines)
        got = [bed_core_of_line(f) for f in lines]
        assert Counter(got) == Counter(want)
        keys = [bed_sort_key(c, order) for c in got]
        assert keys == sorted(keys)


@PROPS
@given(dataset())
def test_probed_lines(ds: dict[str, Any]) -> None:
    reps = ds["reps"]
    n_parts = sum(1 + (r["end"] > reps[r["contig"]].length) for r in ds["peps"])
    with tempfile.TemporaryDirectory() as d:
        out = Path(d)
        _write(out, ds, ds["peps"], ds["prots"])
        text = (out / "mappings.probed").read_text(encoding="utf-8")
        assert text.splitlines()[0] == "# proBed-version 1.0"
        lines = _data(out / "mappings.probed")
        assert len(lines) == n_parts
        assert all(len(f) == 25 for f in lines)
        assert {f[13] for f in lines} == {r["peptide"] for r in ds["peps"]}
        assert all(f[15] == f"sha256:{SHA}" for f in lines)
        assert all(f[12] != "NA" for f in lines)


@PROPS
@given(dataset())
def test_gff3_round_trip(ds: dict[str, Any]) -> None:
    reps = ds["reps"]
    with tempfile.TemporaryDirectory() as d:
        out = Path(d)
        _write(out, ds, ds["peps"], ds["prots"])
        feats = [f for f in read_gff3(out / "mappings.gff3", reps)
                 if f.type == "nucleotide_to_protein_match"]  # fmt: skip
        got = Counter(
            (str(f.attributes["Name"]), f.interval.contig, f.interval.start, f.interval.end,
             f.interval.strand)
            for f in feats
        )  # fmt: skip
        want = Counter(
            (r["peptide"], r["contig"], r["start"], r["end"], r["strand"]) for r in ds["peps"]
        ) + Counter(
            (r["protein_id"], r["contig"], r["start"], r["end"], r["strand"]) for r in ds["prots"]
        )
        assert got == want
        ids = [f.attributes["ID"] for f in feats]
        assert len(ids) == len(set(ids))
        assert all(f.source == "GIPPA" for f in feats)


@PROPS
@given(dataset(), st.randoms(use_true_random=False))
def test_input_order_does_not_matter(ds: dict[str, Any], rnd: Any) -> None:
    peps2, prots2 = list(ds["peps"]), list(ds["prots"])
    rnd.shuffle(peps2)
    rnd.shuffle(prots2)
    with tempfile.TemporaryDirectory() as d:
        a, b = Path(d) / "a", Path(d) / "b"
        _write(a, ds, ds["peps"], ds["prots"])
        _write(b, ds, peps2, prots2)
        for name in FILES:
            assert (a / name).read_bytes() == (b / name).read_bytes(), name


def test_strand_rank_oracle_matches_spec() -> None:
    assert sorted(STRAND_RANK, key=STRAND_RANK.__getitem__) == ["+", "-", "."]
