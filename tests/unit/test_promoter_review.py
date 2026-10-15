"""M2-13 lead decisions after code review (PR #76), items 5-7.

5. ``promoter_windows`` runs end to end through ``giap.steps.registry.run_step``: it reads
   ``replicons`` (circular wraps) and writes its chunks into the ``windows/`` directory
   (``windows/windows_000.fa`` ...), never next to it.
6. One test per surviving mutant of the review (local/scratch/review-m2-13/mutants.py): BED merge
   across the origin (and its max score, and touching across the origin stays apart); the
   ``gap % n`` wrap; circular normalisation in ``promoter_between``; the ``max(a, 0)`` clamp; the
   ``lo1`` wrap (a window ending at end3p == n); the pseudogene filter in the step; the header
   check with a TAB-separated foreign header.
7. ``gene_promoter_features(..., replicons)`` requires ``replicons`` (no silent linear default).
"""

from __future__ import annotations

import importlib
import random
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest
from helpers.promoter_oracle import api

from giap.annotate.operons import Tu
from giap.config import GiapConfig
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

N = 100_000
CHR1 = Replicon("chr1", N, "circular", "chromosome")
REPS = {"chr1": CHR1}


def _gene(tag: str, start: int, end: int, strand: str, contig: str = "chr1",
          ftype: str = "cds") -> Gene:  # fmt: skip
    return Gene(
        locus_tag=tag,
        feature_type=ftype,  # type: ignore[arg-type]
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M",
        nt=None,
        start_codon="ATG",
    )


def _scores(rows: list[tuple[str, int, float]], contig: str = "chr1") -> pd.DataFrame:
    return pd.DataFrame(
        [(contig, s, p, x) for s, p, x in rows], columns=["contig", "strand", "end3p", "score"]
    )


def _tu(tu_id: str, genes: list[Gene], start: int, end: int) -> Tu:
    strand = genes[0].interval.strand
    iv = Interval("chr1", start, end, strand)  # type: ignore[arg-type]
    return Tu(tu_id, "chr1", iv, strand, tuple(g.locus_tag for g in genes), iv)  # type: ignore[arg-type]


# ------------------------------------------------------------------ decision 7: replicons required


def test_gene_promoter_features_without_replicons_is_a_type_error() -> None:
    g = _gene("F", 1000, 1300, "+")
    tus = [_tu("TU_1", [g], 1000, 1300)]
    with pytest.raises(TypeError):
        api("gene_promoter_features")([g], tus, _scores([("+", 900, 0.9)]))


# ------------------------------------------------------------------ decision 6: surviving mutants


def test_own_window_ending_at_end3p_n_is_found_lo1_wrap() -> None:
    """+ gene with S = 200: k = -200 puts the window's 3' end on base n (end3p == n)."""
    g = _gene("A", 200, 500, "+")
    df = api("gene_promoter_features")(
        [g], [_tu("TU_1", [g], 200, 500)], _scores([("+", N, 0.9)]), replicons=REPS
    ).set_index("locus_tag")
    assert df.loc["A", "promoter_score"] == pytest.approx(0.9)
    assert df.loc["A", "n_promoter_windows_tu"] == 1


def test_tu_spanning_the_origin_counts_windows_gap_mod_n() -> None:
    """TU F (99_800, +) then G (50, +) across the origin: G's TU range wraps [99_600 .. 70]."""
    f, g = _gene("F", 99_800, 99_950, "+"), _gene("G", 50, 400, "+")
    rows = [("+", 99_700, 0.9), ("+", 60, 0.8), ("+", 500, 0.99)]  # 500: outside every range
    df = api("gene_promoter_features")(
        [f, g], [_tu("TU_1", [f, g], 99_800, N + 400)], _scores(rows), replicons=REPS
    ).set_index("locus_tag")
    assert df.loc["G", "promoter_score"] == pytest.approx(0.8)
    assert df.loc["G", "promoter_score_tu"] == pytest.approx(0.9)
    assert df.loc["G", "n_promoter_windows_tu"] == 2
    assert df.loc["F", "n_promoter_windows_tu"] == 1


def test_promoter_between_clamps_to_downs_own_windows_max_a_0() -> None:
    """up ends (151) before down's first own window (220): windows in between are not down's."""
    up, down = _gene("U", 100, 150, "+"), _gene("D", 420, 700, "+")  # D own end3p in [220, 440]
    got = api("promoter_between")(up, down, _scores([("+", 200, 0.99), ("+", 300, 0.5)]), CHR1)
    assert got == pytest.approx(0.5)


def test_promoter_between_across_the_origin_is_normalised() -> None:
    """up ends at 99_990, down starts at 100 (own end3p in [99_900 .. 100_000] + [1 .. 120])."""
    up, down = _gene("U", 99_800, 99_990, "+"), _gene("D", 100, 500, "+")
    rows = [("+", 99_950, 0.99), ("+", 99_995, 0.7), ("+", 50, 0.8)]  # 99_950: before up's end
    assert api("promoter_between")(up, down, _scores(rows), CHR1) == pytest.approx(0.8)
    assert api("promoter_between")(up, down, _scores([("+", 99_950, 0.99)]), CHR1) is None


def test_bed_merges_overlap_across_the_origin_with_the_max_score() -> None:
    # end3p 30 -> [99_990, 100_030); end3p 60 -> [20, 60): overlap across the origin
    recs = api("promoters_bed")(_scores([("+", 30, 0.6), ("+", 60, 0.95)]), REPS)
    assert [(r.interval.start, r.interval.end, r.score) for r in recs] == [(99_990, N + 60, 950)]


def test_bed_touching_across_the_origin_stays_apart() -> None:
    # end3p 20 -> [99_980, 100_020); end3p 60 -> [20, 60): touch at the origin, no overlap
    recs = api("promoters_bed")(_scores([("+", 20, 0.6), ("+", 60, 0.95)]), REPS)
    assert sorted((r.interval.start, r.interval.end, r.score) for r in recs) == [
        (20, 60, 950),
        (99_980, N + 20, 600),
    ]


def test_parse_promotech_rejects_a_tab_separated_foreign_header(tmp_path: Path) -> None:
    bad = tmp_path / "sequences_predictions.csv"
    bad.write_text("CHROM\tSEQUENCE\tPRED\nw0\t" + "A" * 40 + "\t0.5\n")
    with pytest.raises(ValueError, match="CHROM"):
        api("parse_promotech")(bad)


# Two further survivors outside the lead's list (mutants 26 and 32); 11 and 29 are equivalent.


def test_parse_promotech_rejects_a_row_without_three_fields(tmp_path: Path) -> None:
    bad = tmp_path / "sequences_predictions.csv"
    bad.write_text("CHROM\tSEQ\tPRED\nw0\t" + "A" * 40 + "\t0.5\textra\n")
    with pytest.raises(ValueError, match="3 tab-separated fields"):
        api("parse_promotech")(bad)


def test_write_chunks_rejects_a_window_longer_than_40(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        api("write_chunks")(tmp_path / "windows", [[("w0", "A" * 41)]])


# ------------------------------------------------------- decision 5: the step through run_step


@pytest.fixture
def reg() -> Iterator[ModuleType]:
    from giap.steps import load_step_modules, registry

    with pytest.MonkeyPatch.context() as mp:
        names = load_step_modules()
        mp.setattr(registry, "_REGISTRY", {})
        for n in names:
            if n != "registry":
                importlib.reload(importlib.import_module(f"giap.steps.{n}"))
        yield registry


def _write_inputs(d: Path) -> tuple[dict[str, Path], dict[str, str], list[Gene]]:
    rng = random.Random(1376)
    seqs = {
        c: "".join(rng.choice("ACGT") for _ in range(n)) for c, n in (("c1", 3000), ("p1", 2000))
    }
    d.mkdir(parents=True)
    (d / "genome.fa").write_text("".join(f">{c}\n{s}\n" for c, s in seqs.items()))
    pq.write_table(
        pa.table({"contig": ["c1", "p1"], "length": [3000, 2000],
                  "topology": ["circular", "linear"], "kind": ["chromosome", "plasmid"]}),
        d / "replicons.parquet",
    )  # fmt: skip
    genes = [
        _gene("A", 100, 400, "+", "c1"),  # k < -60 wraps across c1's origin (circular)
        _gene("P", 1500, 1800, "+", "c1", "pseudogene"),  # never windowed
        _gene("B", 2500, 2800, "-", "c1"),
        _gene("C", 50, 350, "+", "p1"),  # linear: only k >= -10 fits
    ]
    cols = {
        "locus_tag": [g.locus_tag for g in genes],
        "feature_type": [g.feature_type for g in genes],
        "contig": [g.interval.contig for g in genes],
        "start": [g.interval.start for g in genes],
        "end": [g.interval.end for g in genes],
        "strand": [g.interval.strand for g in genes],
        "aa": [g.aa for g in genes],
        "nt": [g.nt for g in genes],
        "start_codon": [g.start_codon for g in genes],
    }
    pq.write_table(pa.table(cols), d / "genes.parquet")
    ins = {"genes": d / "genes.parquet", "replicons": d / "replicons.parquet",
           "genome": d / "genome.fa"}  # fmt: skip
    return ins, seqs, genes


def test_promoter_windows_step_end_to_end_through_the_registry(
    tmp_path: Path, reg: ModuleType
) -> None:
    spec = reg.get("promoter_windows")
    assert "replicons" in spec.inputs, "promoter_windows needs replicons (circular wraps)"
    ins, seqs, genes = _write_inputs(tmp_path / "in")
    out = tmp_path / "out"
    outs = {"windows": out / "windows", "windows_map": out / "windows_map.parquet"}
    reg.run_step("promoter_windows", GiapConfig(seed=1), ins, outs)

    assert (out / "windows").is_dir()
    chunks = sorted(p.name for p in (out / "windows").iterdir())
    assert chunks == ["windows_000.fa"]
    assert not list(out.glob("windows_*.fa")), "chunks must be inside windows/"

    wmap = pd.read_parquet(outs["windows_map"])
    assert list(wmap.columns) == ["header", "contig", "strand", "end3p", "owner", "k"]
    assert "P" not in set(wmap["owner"]), "pseudogenes have no promoter windows"
    reps = {"c1": Replicon("c1", 3000, "circular", "chromosome"),
            "p1": Replicon("p1", 2000, "linear", "plasmid")}  # fmt: skip
    want = [
        (w.contig, w.strand, w.end3p, w.owner, w.k)
        for g in genes
        if g.feature_type != "pseudogene"
        for w in api("gene_windows")(g, reps[g.interval.contig], seqs[g.interval.contig])
    ]
    got = list(wmap[["contig", "strand", "end3p", "owner", "k"]].itertuples(index=False, name=None))
    assert got == want
    assert len(got) == 221 + 221 + 31
    a_end3p = wmap.loc[wmap["owner"] == "A", "end3p"]
    assert (a_end3p > 2900).any(), "circular c1: A's upstream windows wrap to the contig end"

    fasta = (out / "windows" / "windows_000.fa").read_text().split(">")[1:]
    recs = {r.split("\n", 1)[0]: r.split("\n", 1)[1].replace("\n", "") for r in fasta}
    assert set(recs) == set(wmap["header"])
    assert all(len(s) == 40 for s in recs.values())
