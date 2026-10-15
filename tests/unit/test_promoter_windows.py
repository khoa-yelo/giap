"""Promoter windows, chunks, Promotech command/writer/parser, layer and BED (issue M2-13).

Spec: docs/specs/M2-13-promoter-windows-promotech.md examples 1-6, 9, 10 and edge cases;
GIAP_SPEC.md §7.1.1 (Promotech row and gotcha), §7.1.4, §6.5 (BED, origin split). Lead decision
M2-13 Q1 / index 5: headers ``w<index>`` over unique sequences in first-seen order plus a map frame;
``k = 0`` ends right before the start codon.

What PROMOTECH reads (pinned source ``sequences/process_sequences.py``, see
``tests/data/recorded/promotech/56251ad9/annotate-real-1847900/PROVENANCE.md``): the FASTA goes
through ``Bio.SeqIO.parse(..., "fasta")`` (id = first whitespace token), sequences are upper-cased
and must be exactly 40 nt; the output is ``to_csv(sep='\\t', index=None)`` of CHROM/SEQ/PRED.
"""

from __future__ import annotations

import random
from pathlib import Path

import pandas as pd
import pytest
from helpers.promoter_oracle import api, oracle_windows, revcomp

from giap.model.coords import Interval, Replicon
from giap.model.genes import Gene

L = 100_000
_RNG = random.Random(1307)
GENOME = "".join(_RNG.choice("ACGT") for _ in range(L))
CHR1 = Replicon("chr1", L, "circular", "chromosome")
LIN = Replicon("chr1", L, "linear", "unknown")
REAL = Path(__file__).resolve().parents[1] / "data" / "recorded" / "promotech" / "56251ad9"
REAL_RUN = REAL / "annotate-real-1847900"
MAP_COLUMNS = ["header", "contig", "strand", "end3p", "owner", "k"]


def _gene(tag: str, start: int, end: int, strand: str, contig: str = "chr1") -> Gene:
    return Gene(
        locus_tag=tag,
        feature_type="cds",
        interval=Interval(contig, start, end, strand),  # type: ignore[arg-type]
        aa="M",
        nt=None,
        start_codon="ATG",
    )


def _windows(gene: Gene, rep: Replicon = CHR1, seq: str = GENOME) -> list:
    return api("gene_windows")(gene, rep, seq)


def _by_k(ws: list) -> dict:
    return {w.k: w for w in ws}


# ---------------------------------------------------------------- window geometry (examples 1-4)


def test_example1_plus_gene_window_geometry() -> None:
    ws = _windows(_gene("G1", 1000, 1300, "+"))
    assert len(ws) == 221
    by = _by_k(ws)
    assert (by[0].end3p, by[0].seq) == (1000, GENOME[960:1000])
    assert (by[20].end3p, by[20].seq) == (1020, GENOME[980:1020])
    assert (by[-200].end3p, by[-200].seq) == (800, GENOME[760:800])
    for w in ws:
        assert (w.contig, w.strand, w.owner) == ("chr1", "+", "G1")
        assert len(w.seq) == 40


def test_example2_minus_gene_window_geometry() -> None:
    ws = _windows(_gene("G2", 2000, 2900, "-"))
    assert len(ws) == 221
    by = _by_k(ws)
    assert (by[0].end3p, by[0].seq) == (2901, revcomp(GENOME[2900:2940]))
    assert (by[20].end3p, by[20].seq) == (2881, revcomp(GENOME[2880:2920]))
    assert (by[-200].end3p, by[-200].seq) == (3101, revcomp(GENOME[3100:3140]))
    assert {w.strand for w in ws} == {"-"}


def test_example3_circular_wraps_across_origin() -> None:
    ws = _windows(_gene("G3", 100, 400, "+"))
    assert len(ws) == 221
    by = _by_k(ws)
    assert (by[-200].end3p, by[-200].seq) == (99900, GENOME[99860:99900])
    # k = -80: [-20, 20) straddles the origin
    assert (by[-80].end3p, by[-80].seq) == (20, GENOME[99980:] + GENOME[:20])


def test_minus_gene_near_end_wraps_on_circular() -> None:
    by = _by_k(_windows(_gene("G4", 99600, 99900, "-")))
    assert len(by) == 221
    # k = -200: genome [100100, 100140) -> [100, 140); 3'-most base 0-based 100 -> end3p 101
    assert (by[-200].end3p, by[-200].seq) == (101, revcomp(GENOME[100:140]))


def test_example4_linear_contig_drops_windows_outside_contig() -> None:
    ws = _windows(_gene("G3", 100, 400, "+"), LIN)
    assert len(ws) == 81
    assert sorted(w.k for w in ws) == list(range(-60, 21))
    assert min(w.end3p for w in ws) == 40


def test_linear_minus_gene_near_contig_end_is_mirrored() -> None:
    ws = _windows(_gene("G5", L - 400, L - 100, "-"), LIN)
    assert len(ws) == 81
    assert sorted(w.k for w in ws) == list(range(-60, 21))


def test_windows_are_k_ascending_and_match_oracle() -> None:
    for g in (_gene("A", 1000, 1300, "+"), _gene("B", 2000, 2900, "-"), _gene("C", 30, 90, "-")):
        ws = _windows(g)
        s, e, st = g.interval.start, g.interval.end, g.interval.strand
        assert [(w.k, w.end3p, w.seq) for w in ws] == oracle_windows(s, e, st, GENOME, True)


def test_sorf_windows_extend_into_and_beyond_short_gene() -> None:
    by = _by_k(_windows(_gene("S1", 5000, 5015, "+")))
    assert len(by) == 221
    assert by[20].seq == GENOME[4980:5020]  # k = +20 runs past the 15-nt gene end


def test_window_with_n_is_still_produced() -> None:
    seq = GENOME[:960] + "N" * 40 + GENOME[1000:]
    by = _by_k(_windows(_gene("G1", 1000, 1300, "+"), CHR1, seq))
    assert by[0].seq == "N" * 40


# ------------------------------------------------------------------ dedupe and chunk (5, 6)


def test_example5_shared_windows_dedupe_but_map_keeps_both_owners() -> None:
    # B's windows start at 951..1171, A's at 760..980: 30 genome positions (and 40-mers) shared.
    a, b = _gene("A", 1000, 1100, "+"), _gene("B", 1191, 1500, "+")
    wa, wb = _windows(a), _windows(b)
    chunks, wmap = api("dedupe_and_chunk")(wa + wb)
    records = [r for c in chunks for r in c]
    assert len(records) == 442 - 30
    assert [h for h, _ in records] == [f"w{i}" for i in range(412)]
    assert len({s for _, s in records}) == 412
    assert list(wmap.columns) == MAP_COLUMNS
    assert len(wmap) == 442
    shared = wmap[wmap.duplicated("header", keep=False)]
    assert len(shared) == 60 and set(shared["owner"]) == {"A", "B"}
    seq_of = dict(records)
    for w, row in zip(wa + wb, wmap.itertuples(index=False), strict=True):
        assert (row.contig, row.strand, row.end3p, row.owner, row.k) == (
            w.contig, w.strand, w.end3p, w.owner, w.k,
        )  # fmt: skip
        assert seq_of[row.header] == w.seq


def test_headers_are_first_seen_order() -> None:
    a = _windows(_gene("A", 1000, 1300, "+"))
    chunks, wmap = api("dedupe_and_chunk")(a)
    first = chunks[0][0]
    assert first == ("w0", a[0].seq)
    assert wmap.iloc[0]["header"] == "w0" and wmap.iloc[0]["k"] == -200


def _synthetic(n: int) -> list:
    window = api("Window")
    out = []
    for i in range(n):
        x, s = i, []
        for _ in range(40):
            s.append("ACGT"[x % 4])
            x //= 4
        out.append(window("chr1", "+", i + 1, "".join(s), f"G{i}", 0))
    return out


def test_example6_chunks_of_100000() -> None:
    chunks, wmap = api("dedupe_and_chunk")(_synthetic(250_001))
    assert [len(c) for c in chunks] == [100_000, 100_000, 50_001]
    assert chunks[1][0][0] == "w100000" and chunks[2][-1][0] == "w250000"
    assert len(wmap) == 250_001


def test_write_chunks_names_and_records(tmp_path: Path) -> None:
    chunks, _ = api("dedupe_and_chunk")(_synthetic(7), chunk_size=3)
    assert [len(c) for c in chunks] == [3, 3, 1]
    paths = api("write_chunks")(tmp_path, chunks)
    assert [p.name for p in paths] == ["windows_000.fa", "windows_001.fa", "windows_002.fa"]
    from Bio import SeqIO

    got = [(r.id, str(r.seq)) for p in paths for r in SeqIO.parse(p, "fasta")]
    assert got == [r for c in chunks for r in c]


def test_writer_layout_matches_the_real_promotech_input(tmp_path: Path) -> None:
    """Same layout as the chunk Promotech accepted in the real run: ``>id`` then one 40-nt line."""
    real = (REAL_RUN / "windows_000.fa").read_text().splitlines()
    assert all(ln.startswith(">") for ln in real[0::2])
    assert all(len(ln) == 40 for ln in real[1::2])
    chunks, _ = api("dedupe_and_chunk")(_windows(_gene("A", 1000, 1300, "+")))
    (path,) = api("write_chunks")(tmp_path, chunks)
    lines = path.read_text().splitlines()
    assert len(lines) == 2 * 221
    assert lines[0::2] == [f">w{i}" for i in range(221)]
    assert all(len(ln) == 40 and set(ln) <= set("ACGTN") for ln in lines[1::2])


def test_window_shorter_than_40_is_never_written(tmp_path: Path) -> None:
    with pytest.raises(ValueError):
        api("write_chunks")(tmp_path, [[("w0", "ACGT" * 9 + "ACG")]])


def test_dedupe_and_write_are_deterministic(tmp_path: Path) -> None:
    ws = _windows(_gene("A", 1000, 1300, "+")) + _windows(_gene("B", 2000, 2900, "-"))
    c1, m1 = api("dedupe_and_chunk")(ws)
    c2, m2 = api("dedupe_and_chunk")(list(ws))
    assert c1 == c2
    pd.testing.assert_frame_equal(m1, m2)
    p1 = api("write_chunks")(tmp_path / "a", c1)
    p2 = api("write_chunks")(tmp_path / "b", c2)
    assert [p.read_bytes() for p in p1] == [p.read_bytes() for p in p2]


# ------------------------------------------------------------ Promotech command, parser (10)


def test_promotech_command_targeted_rf_hot() -> None:
    from giap.annotate.promotech_cmd import promotech_command

    home, chunk, out = Path("/env/share/promotech"), Path("/w/windows_000.fa"), Path("/w/p000")
    cmd = promotech_command(chunk, out, home)
    assert cmd[:2] == ["python", str(home / "promotech.py")]
    assert "-s" in cmd and "-g" not in cmd and "-pg" not in cmd
    assert cmd[cmd.index("-m") + 1] == "RF-HOT"
    assert cmd[cmd.index("-f") + 1] == str(chunk)
    assert Path(cmd[cmd.index("-o") + 1]) == out


def _emulate_promotech(fasta: Path, outdir: Path) -> Path:
    """Promotech's own I/O path (process_sequences.py) with a dummy model: SeqIO ids -> CSV."""
    from Bio import SeqIO

    recs = list(SeqIO.parse(fasta, "fasta"))
    df = pd.DataFrame(
        {
            "CHROM": [r.id for r in recs],
            "SEQ": [str(r.seq).upper() for r in recs],
            "PRED": [(i % 1000) / 1000 for i in range(len(recs))],
        }
    )
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "sequences_predictions.csv"
    df.to_csv(path, index=None, sep="\t", columns=None)
    return path


def test_example10_header_round_trip_w0_to_w99999(tmp_path: Path) -> None:
    chunks, _ = api("dedupe_and_chunk")(_synthetic(100_000))
    (fa,) = api("write_chunks")(tmp_path / "in", chunks)
    csv = _emulate_promotech(fa, tmp_path / "out")
    scores = api("parse_promotech")(csv)
    assert list(scores) == [f"w{i}" for i in range(100_000)]
    assert scores["w1234"] == pytest.approx(0.234)


def test_parse_real_promotech_output() -> None:
    from Bio import SeqIO

    scores = api("parse_promotech")(REAL_RUN / "sequences_predictions.csv")
    fasta = [(r.id, str(r.seq)) for r in SeqIO.parse(REAL_RUN / "windows_000.fa", "fasta")]
    assert list(scores) == [h for h, _ in fasta]
    assert len(scores) == 86
    assert scores["KABHIH_0001"] == 0.1773268671071446
    assert scores["KABHIH_0099"] == 0.09589302363189868
    assert all(isinstance(v, float) and 0.0 <= v <= 1.0 for v in scores.values())


def test_parse_promotech_rejects_a_foreign_header(tmp_path: Path) -> None:
    bad = tmp_path / "sequences_predictions.csv"
    bad.write_text("CHROM,SEQ,PRED\nw0," + "A" * 40 + ",0.5\n")
    with pytest.raises(ValueError):
        api("parse_promotech")(bad)


# ------------------------------------------------------------------- layer and BED (example 9)


def test_promoter_layer_unique_positions_with_scores() -> None:
    a, b = _gene("A", 1000, 1100, "+"), _gene("B", 1191, 1500, "+")
    chunks, wmap = api("dedupe_and_chunk")(_windows(a) + _windows(b))
    scores = {h: (int(h[1:]) % 97) / 97 for c in chunks for h, _ in c}
    layer = api("promoter_layer")(wmap, scores)
    assert list(layer.columns) == ["contig", "strand", "end3p", "score"]
    assert not layer.duplicated(["contig", "strand", "end3p"]).any()
    assert len(layer) == 412
    got = {(r.contig, r.strand, r.end3p): r.score for r in layer.itertuples(index=False)}
    for r in wmap.itertuples(index=False):
        assert got[(r.contig, r.strand, r.end3p)] == scores[r.header]
    pd.testing.assert_frame_equal(layer, api("promoter_layer")(wmap, scores))


def _scores(rows: list[tuple[str, str, int, float]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=["contig", "strand", "end3p", "score"])


def _bed(rows: list[tuple[str, str, int, float]], reps: dict | None = None) -> list:
    return api("promoters_bed")(_scores(rows), reps or {"chr1": CHR1})


def test_example9_overlapping_windows_merge_with_max_score() -> None:
    from giap.io.bed import bed_lines

    (rec,) = _bed([("chr1", "+", 1000, 0.6), ("chr1", "+", 1010, 0.9)])
    assert bed_lines(rec, CHR1) == ["chr1\t960\t1010\tpromoter_1\t900\t+"]


def test_bed_drops_windows_below_yes_and_keeps_threshold() -> None:
    recs = _bed([("chr1", "+", 1000, 0.5), ("chr1", "+", 5000, 0.4999)])
    assert [(r.interval.start, r.interval.end, r.score) for r in recs] == [(960, 1000, 500)]


def test_bed_minus_strand_window_interval_and_merge() -> None:
    # '-' end3p p -> genome [p - 1, p + 39); 2901 -> [2900, 2940), 2911 -> [2910, 2950)
    recs = _bed([("chr1", "-", 2901, 0.7), ("chr1", "-", 2911, 0.55), ("chr1", "+", 2940, 0.8)])
    got = {(r.interval.start, r.interval.end, r.interval.strand, r.score) for r in recs}
    assert got == {(2900, 2950, "-", 700), (2900, 2940, "+", 800)}  # strands never merge


def test_bed_chain_of_overlaps_merges_but_touching_does_not() -> None:
    rows = [
        ("chr1", "+", 1000, 0.6), ("chr1", "+", 1030, 0.7),  # [960,1000) + [990,1030) overlap
        ("chr1", "+", 1070, 0.65),  # [1030,1070) only touches [990,1030): separate
    ]  # fmt: skip
    recs = _bed(rows)
    assert [(r.interval.start, r.interval.end, r.score, r.name) for r in recs] == [
        (960, 1030, 700, "promoter_1"),
        (1030, 1070, 650, "promoter_2"),
    ]


def test_bed_numbering_follows_genome_order() -> None:
    reps = {"chr1": CHR1, "p1": Replicon("p1", 5000, "circular", "plasmid")}
    recs = _bed([("p1", "+", 300, 0.9), ("chr1", "-", 501, 0.6), ("chr1", "+", 500, 0.7)], reps)
    order = [(r.interval.contig, r.interval.strand, r.name) for r in recs]
    assert sorted(order, key=lambda t: t[2]) == [
        ("chr1", "+", "promoter_1"), ("chr1", "-", "promoter_2"), ("p1", "+", "promoter_3"),
    ]  # fmt: skip


def test_bed_origin_crossing_record_splits_into_parts() -> None:
    from giap.io.bed import bed_lines

    (rec,) = _bed([("chr1", "+", 20, 0.8)])  # window [-20, 20) -> unwrapped [99980, 100020)
    assert (rec.interval.start, rec.interval.end) == (99980, 100020)
    lines = bed_lines(rec, CHR1)
    assert lines == [
        f"chr1\t99980\t100000\t{rec.name}_part1\t800\t+",
        f"chr1\t0\t20\t{rec.name}_part2\t800\t+",
    ]


# ------------------------------------------------------------------------------------ steps


def test_promoter_steps_are_registered(monkeypatch: pytest.MonkeyPatch) -> None:
    import importlib

    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    names = load_step_modules()
    monkeypatch.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"giap.steps.{n}"))
    for n in ("promoter_windows", "promoter_layer"):
        spec = reg.get(n)
        assert (spec.kind, spec.module) == ("py", "core")
