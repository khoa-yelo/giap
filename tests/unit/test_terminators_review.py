"""M2-15 lead decisions after code review (PR #75), items 4-7 (docs/specs/M2-15-terminators.md).

4. The TerminatorNet parser never skips data lines by prefix: a contig ID starting with ``name``
   parses normally.
5. BacTermFinder (provisional) output directory: empty -> no calls; exactly one CSV -> parsed;
   several CSVs -> ValueError naming them.
6. ``terminator_between(up, down)`` raises ValueError unless both genes share contig and strand.
7. Tests that kill the reviewer's surviving mutants: the ``-L`` shift in ``_max_in`` at the origin;
   the downstream window slack with no next gene; BED score ``round`` (not ``int``); NaN
   ``share_matched`` for a contig without primary calls; cross-check contig order = FASTA order;
   the ``bactermfinder`` gate; ``cross_check_transtermhp: false``; the layer sorted.

Step tests build their own small inputs (replicons.parquet in a non-alphabetical FASTA order,
TerminatorNet TSV, TransTermHP ``.tt``) and resolve step keys by name, as test_terminator_steps.py.
"""

from __future__ import annotations

import importlib
import math
from collections.abc import Iterator
from pathlib import Path
from types import ModuleType
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from giap.annotate.terminators import (
    Terminator,
    crosscheck,
    parse_terminatornet,
    terminator_between,
    terminator_conf_downstream,
    terminators_bed,
)
from giap.config import GiapConfig, resolve_config
from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

CHR1 = Replicon("chr1", 100_000, "circular", "chromosome")
CIRC = Replicon("c", 10_000, "circular", "chromosome")


def gene(tag: str, start: int, end: int, strand: str, contig: str = "chr1") -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M" * max(1, (end - start) // 3 - 1),
        nt=None,
        start_codon="ATG",
    )


def term(
    start: int,
    end: int,
    strand: str,
    conf: float,
    source: str = "TerminatorNet",
    contig: str = "chr1",
) -> Terminator:
    return Terminator(Interval(contig, start, end, strand), conf, source)  # type: ignore[arg-type]


def tn_line(
    contig: str, i: int, start1: int, stop1: int, strand: str, pred: float, gh: float = -10.5
) -> str:
    return "\t".join(
        [contig, f"TERM_{i}", "AAAAAAAA", "GCCGCCAAAAGGCGGC", "TTTTTTTT", "((((((....))))))",
         repr(gh), repr(pred), str(start1), str(stop1), strand]
    )  # fmt: skip


# ------------------------------------------------------------------- item 4: no prefix skip


@pytest.mark.parametrize("contig", ["name_contig", "names", "name1"])
def test_contig_ids_starting_with_name_parse(tmp_path: Path, contig: str) -> None:
    reps = {contig: Replicon(contig, 5_000, "circular", "plasmid")}
    p = tmp_path / "terminatornet.tsv"
    p.write_text(tn_line(contig, 1, 1201, 1240, "-", 0.87, -12.4) + "\n")
    assert parse_terminatornet(p, reps) == [
        Terminator(Interval(contig, 1200, 1240, "-"), 0.87, "TerminatorNet", -12.4)
    ]


def test_name_contig_among_other_contigs_keeps_every_call(tmp_path: Path) -> None:
    reps = {
        "chr1": CHR1,
        "name_contig": Replicon("name_contig", 5_000, "circular", "plasmid"),
    }
    p = tmp_path / "terminatornet.tsv"
    p.write_text(
        tn_line("chr1", 1, 1888, 1907, "+", 0.88)
        + "\n"
        + tn_line("name_contig", 1, 101, 140, "+", 0.6)
        + "\n"
        + "name_contig\tNO TERMINATORS FOUND\n"
        + tn_line("name_contig", 2, 301, 340, "-", 0.7)
        + "\n"
    )
    got = parse_terminatornet(p, reps)
    assert [(t.interval.contig, t.interval.start, t.interval.strand) for t in got] == [
        ("chr1", 1887, "+"),
        ("name_contig", 100, "+"),
        ("name_contig", 300, "-"),
    ]


# ------------------------------------------------------------------- item 6: between preconditions


@pytest.mark.parametrize(
    ("up", "down"),
    [
        (gene("U", 100, 400, "+"), gene("D", 420, 700, "-")),
        (gene("U", 420, 700, "-"), gene("D", 100, 400, "+")),
        (gene("U", 100, 400, "+"), gene("D", 420, 700, "+", contig="chr2")),
    ],
)
def test_between_requires_same_contig_and_strand(up: Gene, down: Gene) -> None:
    ts = [term(360, 400, "+", 0.8), term(360, 400, "-", 0.8)]
    with pytest.raises(ValueError):
        terminator_between(up, down, ts, CHR1)


# ------------------------------------------------------------------- item 7: pure functions


def test_between_window_below_zero_matches_a_call_stored_before_the_origin() -> None:
    """``+`` up gene ending at 30 on a circular contig: window [-20, 100) wraps to [9980, 10000)."""
    up, down = gene("U", 0, 30, "+", "c"), gene("D", 100, 400, "+", "c")
    t = term(9_985, 9_995, "+", 0.7, contig="c")
    assert terminator_between(up, down, [t], CIRC) == 0.7


def test_downstream_minus_gene_near_origin_wraps_below_zero() -> None:
    """``-`` gene at [100, 400) without next gene: window [-200, 150) includes [9800, 10000)."""
    g = gene("M", 100, 400, "-", "c")
    assert (
        terminator_conf_downstream(g, None, [term(9_900, 9_950, "-", 0.4, contig="c")], CIRC) == 0.4
    )
    assert (
        terminator_conf_downstream(g, None, [term(9_799, 9_850, "-", 0.4, contig="c")], CIRC)
        is None
    )


def test_downstream_without_next_gene_includes_the_slack_before_the_stop() -> None:
    """Q3 window is [stop - 50, stop + 300) in transcription direction, as with a next gene."""
    g = gene("G", 100, 400, "+")
    assert terminator_conf_downstream(g, None, [term(360, 400, "+", 0.4)], CHR1) == 0.4
    assert terminator_conf_downstream(g, None, [term(350, 390, "+", 0.4)], CHR1) == 0.4
    assert terminator_conf_downstream(g, None, [term(349, 390, "+", 0.4)], CHR1) is None
    m = gene("M", 1_000, 1_300, "-")  # window [700, 1050)
    assert terminator_conf_downstream(m, None, [term(1_000, 1_040, "-", 0.45)], CHR1) == 0.45
    assert terminator_conf_downstream(m, None, [term(1_010, 1_050, "-", 0.45)], CHR1) == 0.45
    assert terminator_conf_downstream(m, None, [term(1_010, 1_051, "-", 0.45)], CHR1) is None


@pytest.mark.parametrize(("conf", "score"), [(0.8889, 889), (0.4996, 500), (0.9999, 1000)])
def test_bed_score_is_rounded_not_truncated(conf: float, score: int) -> None:
    (rec,) = terminators_bed([term(1000, 1040, "+", conf)], {"chr1": CHR1})
    assert rec.score == score


def test_crosscheck_share_is_nan_for_a_contig_without_primary_calls() -> None:
    prim = [term(1000, 1040, "+", 0.9)]
    oth = [term(1000, 1040, "+", 0.9), term(500, 540, "-", 0.9, "TransTermHP", contig="chr2")]
    rows = {r.contig: r for r in crosscheck(prim, oth)}
    r2 = rows["chr2"]
    assert (r2.n_terminatornet, r2.n_transtermhp, r2.n_matched) == (0, 1, 0)
    assert math.isnan(r2.share_matched)
    assert rows["chr1"].share_matched == 1.0
    assert rows["all"].share_matched == 1.0


# ------------------------------------------------------------------- item 5 and 7: terminator_layer

# FASTA order is deliberately not alphabetical.
REPS = [("zeta_chr", 50_000, "circular", "chromosome"), ("alpha_pl", 8_000, "circular", "plasmid")]
XC_HEADER = ["contig", "n_terminatornet", "n_transtermhp", "n_matched", "share_matched"]
BTF_HEADER = (  # the real table (fix-bactermfinder-parser); SampleName coordinates are 0-based
    "SampleName,probability_binary,probability_ENAC,probability_PS2,probability_NCP,"
    "probability_mean\n"
)
BTF_ZETA = "1_zeta_chr_4000_4040_+,0.9,0.9,0.9,0.9,0.95\n"


def _fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"giap.steps.{n}"))
    return reg


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield _fresh_steps(mp)


@pytest.fixture
def files(tmp_path: Path) -> dict[str, Path]:
    d = tmp_path / "in"
    d.mkdir()
    reps = d / "replicons.parquet"
    pq.write_table(
        pa.table(
            {
                "contig": [r[0] for r in REPS],
                "length": [r[1] for r in REPS],
                "topology": [r[2] for r in REPS],
                "kind": [r[3] for r in REPS],
            }
        ),
        reps,
    )
    tn = d / "terminatornet.tsv"  # unsorted: plasmid first, then descending starts, mixed strands
    tn.write_text(
        "\n".join(
            [
                tn_line("alpha_pl", 1, 5001, 5040, "+", 0.6),
                tn_line("alpha_pl", 2, 101, 140, "-", 0.7),
                tn_line("zeta_chr", 1, 30001, 30040, "+", 0.8),
                tn_line("zeta_chr", 2, 2001, 2040, "+", 0.9),
                tn_line("zeta_chr", 3, 10001, 10040, "-", 0.55),
            ]
        )
        + "\n"
    )
    tt = d / "transtermhp.tt"
    tt.write_text(
        "TransTermHP v2.09\n\n"
        "SEQUENCE zeta_chr  (length 50000)\n\n"
        "  TERM 1         2001 - 2040     + F    93 -11.5 -3.8 | opp_overlap 2001\n"
        "SEQUENCE alpha_pl  (length 8000)\n\n"
        "  TERM 2         5001 - 5040     + F    80 -10.1 -4.0 | opp_overlap 5001\n"
    )
    return {"replicons": reps, "terminatornet": tn, "transtermhp": tt}


def _key_file(key: str, avail: dict[str, Path]) -> Path | None:
    k = key.lower()
    for name in ("replicons", "terminatornet", "bactermfinder"):
        if name in k:
            return avail.get(name)
    if "transterm" in k or k == "tt":
        return avail.get("transtermhp")
    return None


def _run(
    reg: ModuleType, avail: dict[str, Path], out: Path, cfg: GiapConfig
) -> tuple[list[dict[str, Any]], list[list[str]]]:
    spec = reg.get("terminator_layer")
    ins: dict[str, Path] = {}
    for key in spec.inputs:
        p = _key_file(key, avail)
        assert p is not None, f"cannot provide {key!r}"
        ins[key] = p
    for key in spec.optional_inputs:
        p = _key_file(key, avail)
        if p is not None:
            ins[key] = p
    xc_key = next(k for k in spec.outputs if "crosscheck" in k)
    pq_key = next(k for k in spec.outputs if "crosscheck" not in k)
    outs = {xc_key: out / "terminators_crosscheck.tsv", pq_key: out / "terminators.parquet"}
    reg.run_step("terminator_layer", cfg, ins, outs)
    xc = outs[xc_key]
    xc_rows = [ln.split("\t") for ln in xc.read_text().splitlines()] if xc.exists() else []
    return pq.read_table(outs[pq_key]).to_pylist(), xc_rows


def _cfg(*over: str) -> GiapConfig:
    return resolve_config(None, ["seed=1", *over])


def test_crosscheck_rows_follow_fasta_order(
    steps: ModuleType, files: dict[str, Path], tmp_path: Path
) -> None:
    _, xc = _run(steps, files, tmp_path / "o", _cfg())
    assert xc[0] == XC_HEADER
    contigs = [r[0] for r in xc[1:]]
    assert contigs[-1] == "all"
    assert contigs[:-1] == ["zeta_chr", "alpha_pl"]
    assert [int(v) for v in xc[1][1:4]] == [3, 1, 1]
    assert [int(v) for v in xc[2][1:4]] == [2, 1, 1]


def test_cross_check_disabled_writes_header_only(
    steps: ModuleType, files: dict[str, Path], tmp_path: Path
) -> None:
    layer, xc = _run(
        steps, files, tmp_path / "o", _cfg("annotation.terminators.cross_check_transtermhp=false")
    )
    assert xc in ([], [XC_HEADER])
    assert {r["source"] for r in layer} == {"TerminatorNet"}


def test_layer_rows_are_sorted_in_genome_order(
    steps: ModuleType, files: dict[str, Path], tmp_path: Path
) -> None:
    layer, _ = _run(steps, files, tmp_path / "o", _cfg())
    got = [(r["contig"], r["start"], r["end"], r["strand"]) for r in layer]
    assert got == [
        ("zeta_chr", 2000, 2040, "+"),
        ("zeta_chr", 10000, 10040, "-"),
        ("zeta_chr", 30000, 30040, "+"),
        ("alpha_pl", 100, 140, "-"),
        ("alpha_pl", 5000, 5040, "+"),
    ]


def _btf_dir(tmp_path: Path, tables: dict[str, str]) -> Path:
    d = tmp_path / "bactermfinder"
    d.mkdir()
    for name, body in tables.items():
        (d / name).write_text(body)
    return d


def test_bactermfinder_gate_off_ignores_its_output(
    steps: ModuleType, files: dict[str, Path], tmp_path: Path
) -> None:
    d = _btf_dir(tmp_path, {"calls.csv": BTF_HEADER + BTF_ZETA})
    layer, _ = _run(steps, {**files, "bactermfinder": d}, tmp_path / "o", _cfg())
    assert "BacTermFinder" not in {r["source"] for r in layer}


def test_bactermfinder_single_csv_is_parsed(
    steps: ModuleType, files: dict[str, Path], tmp_path: Path
) -> None:
    d = _btf_dir(tmp_path, {"calls.csv": BTF_HEADER + BTF_ZETA})
    cfg = _cfg("annotation.terminators.bactermfinder=true")
    layer, _ = _run(steps, {**files, "bactermfinder": d}, tmp_path / "o", cfg)
    btf = [(r["contig"], r["start"], r["end"], r["strand"], r["confidence"])
           for r in layer if r["source"] == "BacTermFinder"]  # fmt: skip
    assert btf == [("zeta_chr", 4000, 4040, "+", 0.95)]
    assert sum(r["source"] == "TerminatorNet" for r in layer) == 5


def test_bactermfinder_empty_directory_means_no_calls(
    steps: ModuleType, files: dict[str, Path], tmp_path: Path
) -> None:
    d = _btf_dir(tmp_path, {})
    cfg = _cfg("annotation.terminators.bactermfinder=true")
    layer, _ = _run(steps, {**files, "bactermfinder": d}, tmp_path / "o", cfg)
    assert {r["source"] for r in layer} == {"TerminatorNet"}
    assert len(layer) == 5


def test_bactermfinder_several_csvs_raise_naming_them(
    steps: ModuleType, files: dict[str, Path], tmp_path: Path
) -> None:
    from giap.steps.registry import StepError

    d = _btf_dir(
        tmp_path,
        {
            "a_calls.csv": BTF_HEADER + BTF_ZETA,
            "b_calls.csv": BTF_HEADER + "1_alpha_pl_200_240_-,0.6,0.6,0.6,0.6,0.65\n",
        },
    )
    cfg = _cfg("annotation.terminators.bactermfinder=true")
    with pytest.raises((StepError, ValueError)) as ei:
        _run(steps, {**files, "bactermfinder": d}, tmp_path / "o", cfg)
    err = ei.value
    cause = err.cause if isinstance(err, StepError) else err
    assert isinstance(cause, ValueError)
    assert "a_calls.csv" in str(cause) and "b_calls.csv" in str(cause)
