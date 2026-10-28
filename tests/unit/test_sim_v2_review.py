"""gippasim v2: lead decisions after code review of PR #38 (docs/specs/M2-05-simulator-v2.md, end).

Decision 1: no gene, on either strand, has a second ``AGGAGG`` core with a 5-10 nt spacer before
its start codon besides its true RBS (read here from the genome in the gene's reading orientation).
Decision 2 lives in the M0 / M3-18 hash tests (``truth/features.tsv`` added there).
Decision 3: layer values match the M2-15/M2-18/M2-19 parse steps: region ids ``<type>_<n>``
numbered per type by contig order then start; official tool casing (``TerminatorNet``,
``TransTermHP``, ``geNOMAD``, ``antiSMASH``, ``Bakta``); BGC score ``NA``; every region with an
interval (false positives included) lists the genes it holds by M2-19's rule (>= 50 % of the CDS).
Decision 4: one test per surviving mutant, each checked against the genome sequence.

Conventions pinned here (implementer's docstrings in ``gippasim/signals.py``): a terminator is a
7-bp G/C stem, a 4-nt loop, the reverse-complement stem and an 8-nt U-tract (26 nt) starting
5-15 nt after the TU's last stop codon; ``truth/tus.tsv`` ``start``/``end`` span the TU's genes
(first start codon to last stop codon); the ``-10`` box ``TATAAT`` ends 3-9 nt before the TSS and
the TSS lies 2-10 nt upstream of the first gene's RBS.
"""

from __future__ import annotations

import csv
import dataclasses
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pytest

import gippasim
from gippa.annotate.bakta import parse_bakta_json
from gippasim import simulate

N_SEEDS = 40
RBS = "AGGAGG"
SPACERS = range(5, 11)
STARTS = {"ATG", "GTG", "TTG"}
STOPS = {"TAA", "TAG", "TGA"}
COMP = str.maketrans("ACGTN", "TGCAN")
TERM_SOURCES = {"TerminatorNet", "TransTermHP"}
REGION_SOURCES = {
    "prophage": {"geNOMAD"},
    "bgc": {"antiSMASH"},
    "crispr": {"Bakta"},
    "plasmid": {"Bakta", "geNOMAD", "Bakta+geNOMAD"},
}


def revcomp(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


class Out:
    """One simulator output: genome, reading-orientation access, truth and layer tables."""

    def __init__(self, out: Path) -> None:
        self.out = out
        ann = parse_bakta_json(out / "bakta_like.json")
        self.seqs = {k: v.upper() for k, v in dict(ann.sequences).items()}
        self.order = {c: k for k, c in enumerate(ann.model.replicons)}
        self.genes = list(ann.model.genes)
        self.by_locus = {g.locus_tag: g for g in self.genes}

    def truth(self, name: str) -> list[dict[str, str]]:
        return read_rows(self.out / "truth" / f"{name}.tsv")

    def layer(self, name: str) -> list[dict[str, str]]:
        return read_rows(self.out / "precomputed" / f"{name}.tsv")

    def read(self, contig: str, strand: str) -> str:
        """Contig in reading orientation of ``strand`` (``-``: reverse complement)."""
        s = self.seqs[contig]
        return s if strand == "+" else revcomp(s)

    def rpos(self, contig: str, strand: str, g0: int) -> int:
        """Reading-orientation index of the 0-based genome position ``g0``."""
        n = len(self.seqs[contig])
        return g0 % n if strand == "+" else n - 1 - g0 % n

    def window(self, contig: str, strand: str, r0: int, k: int) -> str:
        """``k`` nt from reading index ``r0`` (may be negative; circular)."""
        r = self.read(contig, strand)
        n = len(r)
        return "".join(r[(r0 + i) % n] for i in range(k))

    def gene_start(self, locus: str) -> tuple[str, str, int]:
        """(contig, strand, reading index of the first base of the start codon)."""
        iv = self.by_locus[locus].interval
        g0 = iv.start if iv.strand == "+" else iv.end - 1
        return iv.contig, iv.strand, self.rpos(iv.contig, iv.strand, g0)


def full_cfg(seed: int) -> Any:
    """SIM100K at full signal strength; odd seeds add an origin-crossing gene, half the seeds a
    high false-positive rate (more false regions for the membership test)."""
    return dataclasses.replace(
        gippasim.SIM100K,
        rbs_strength=1.0,
        promoter_strength=1.0,
        terminator_strength=1.0,
        origin_crossing_gene=bool(seed % 2),
        layer_false_positive_rate=0.6 if seed % 4 < 2 else 0.05,
    )


@pytest.fixture(scope="module")
def runs(tmp_path_factory: pytest.TempPathFactory) -> list[Out]:
    root = tmp_path_factory.mktemp("simv2_review")
    out = []
    for seed in range(1, N_SEEDS + 1):
        simulate(full_cfg(seed), seed, root / f"s{seed}")
        out.append(Out(root / f"s{seed}"))
    simulate(gippasim.SIM100K, gippasim.SIM100K_SEED, root / "sim100k")
    out.append(Out(root / "sim100k"))
    return out


def true_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    assert all(r["true_element"] in ("true", "false") for r in rows)
    return [r for r in rows if r["true_element"] == "true"]


# ---------------------------------------------------------------- decision 1: one SD core


def test_no_second_sd_core_in_spacer_window(runs: list[Out]) -> None:
    bad: list[str] = []
    checked = Counter[str]()
    for run in runs:
        for r in true_rows(run.truth("rbs")):
            contig, strand, p = run.gene_start(r["locus_tag"])
            sp = int(r["spacer_bp"])
            hits = [k for k in SPACERS if run.window(contig, strand, p - k - 6, 6) == RBS]
            checked[strand] += 1
            if hits != [sp]:
                bad.append(
                    f"{run.out.name}:{r['locus_tag']}({strand}) spacer {sp}, cores at {hits}"
                )
    assert checked["+"] > 500 and checked["-"] > 500, checked
    assert not bad, f"{len(bad)} genes: {bad[:10]}"


# ---------------------------------------------------------------- decision 3: layer values


def test_region_ids_numbered_per_type(runs: list[Out]) -> None:
    for run in runs:
        for rows in (run.layer("regions"), run.truth("regions")):
            by_type: dict[str, list[dict[str, str]]] = defaultdict(list)
            for r in rows:
                by_type[r["region_type"]].append(r)
            for kind, rs in by_type.items():
                rs.sort(key=lambda r: (run.order[r["contig"]], int(r["start"]), int(r["end"])))
                ids = [r["region_id"] for r in rs]
                assert ids == [f"{kind}_{k}" for k in range(1, len(rs) + 1)], (run.out.name, ids)


def test_region_tool_names_and_bgc_score(runs: list[Out]) -> None:
    seen = set()
    for run in runs:
        for r in run.layer("regions"):
            kind = r["region_type"]
            seen.add(kind)
            assert r["source_tool"] in REGION_SOURCES[kind], (run.out.name, r)
            if kind == "bgc":
                assert r["score"] == "NA", (run.out.name, r)
    assert {"prophage", "bgc", "crispr", "plasmid"} <= seen


def test_terminator_source_casing(runs: list[Out]) -> None:
    for run in runs:
        srcs = {r["source"] for r in run.layer("terminators")}
        assert srcs and srcs <= TERM_SOURCES, (run.out.name, srcs)


def overlap_frac(run: Out, gene: Any, contig: str, r0: int, r1: int) -> float:
    """Share of ``gene``'s CDS inside the region ``[r0, r1)`` (0-based; origin-aware)."""
    iv = gene.interval
    if iv.contig != contig:
        return 0.0
    n = len(run.seqs[contig])
    ov = sum(max(0, min(iv.end, r1 + k) - max(iv.start, r0 + k)) for k in (-n, 0, n))
    return ov / (iv.end - iv.start)


def test_region_genes_follow_half_overlap_rule(runs: list[Out]) -> None:
    """Region membership (M2-19): a gene is listed iff >= 50 % of its CDS lies inside; false
    positives list their genes like true ones."""
    partial = Counter[str]()
    false_with_genes = 0
    for run in runs:
        for r in run.truth("regions"):
            if r["region_type"] == "crispr" or r["start"] == "NA":
                continue
            r0, r1 = int(r["start"]) - 1, int(r["end"])
            want = set()
            for g in run.genes:
                f = overlap_frac(run, g, r["contig"], r0, r1)
                if f >= 0.5:
                    want.add(g.locus_tag)
                if 0 < f < 1:
                    partial["in" if f >= 0.5 else "out"] += 1
            got = set() if r["genes"] in ("", "NA") else set(r["genes"].split(";"))
            assert got == want, (run.out.name, r["region_id"], sorted(got ^ want))
            false_with_genes += r["true_element"] == "false" and bool(got)
    assert false_with_genes > 0
    assert partial["in"] > 0 and partial["out"] > 0, partial  # the boundary is exercised


# ---------------------------------------------------------------- decision 4: mutants


def terminator_seq(run: Out, r: dict[str, str]) -> str:
    s = run.seqs[r["contig"]][int(r["start"]) - 1 : int(r["end"])]
    return s if r["strand"] == "+" else revcomp(s)


def test_true_terminators_are_hairpin_plus_u_tract(runs: list[Out]) -> None:
    n = 0
    for run in runs:
        for r in true_rows(run.truth("terminators")):
            t = terminator_seq(run, r)
            assert len(t) == 26, (run.out.name, r, t)
            assert set(t[:7]) <= set("GC"), (run.out.name, r, t)
            assert t[11:18] == revcomp(t[:7]), (run.out.name, r, t)
            n += 1
    assert n > 500


def test_true_terminators_end_in_8nt_u_tract(runs: list[Out]) -> None:
    for run in runs:
        for r in true_rows(run.truth("terminators")):
            t = terminator_seq(run, r)
            assert len(t) == 26 and t[18:] == "TTTTTTTT", (run.out.name, r, t)


def tu_rows(run: Out) -> list[dict[str, str]]:
    return run.truth("tus")


def test_terminator_offset_after_last_stop(runs: list[Out]) -> None:
    """The terminator starts 5-15 nt after the TU's last stop codon (all values occur)."""
    offs = Counter[int]()
    for run in runs:
        terms: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
        for r in true_rows(run.truth("terminators")):
            terms[(r["contig"], r["strand"])].append(r)
        for tu in tu_rows(run):
            if tu["terminator"] != "true":
                continue
            key = (tu["contig"], tu["strand"])
            if tu["strand"] == "+":
                cand = [int(t["start"]) - int(tu["end"]) - 1 for t in terms.get(key, [])]
            else:
                cand = [int(tu["start"]) - int(t["end"]) - 1 for t in terms.get(key, [])]
            near = [d for d in cand if 0 <= d <= 40]
            assert len(near) == 1, (run.out.name, tu["tu_id"], sorted(cand, key=abs)[:3])
            offs[near[0]] += 1
    assert set(offs) == set(range(5, 16)), sorted(offs.items())


def test_tus_span_first_start_to_last_stop(runs: list[Out]) -> None:
    n = 0
    for run in runs:
        for tu in tu_rows(run):
            genes = [run.by_locus[x].interval for x in tu["genes"].split(";")]
            assert int(tu["start"]) - 1 == min(g.start for g in genes), (run.out.name, tu)
            assert int(tu["end"]) == max(g.end for g in genes), (run.out.name, tu)
            s = run.seqs[tu["contig"]]
            a, b = int(tu["start"]) - 1, int(tu["end"])
            five, three = (s[a : a + 3], s[b - 3 : b])
            if tu["strand"] == "-":
                five, three = revcomp(s[b - 3 : b]), revcomp(s[a : a + 3])
            assert five in STARTS, (run.out.name, tu, five)
            assert three in STOPS, (run.out.name, tu, three)
            n += 1
    assert n > 500


def test_tss_downstream_of_promoter_and_upstream_of_rbs(runs: list[Out]) -> None:
    """``tss`` (1-based, TU strand): ``TTGACA``+17+``TATAAT`` ends 3-9 nt before it; the first
    gene's RBS begins 2-10 nt after it (each value occurs)."""
    gaps10, utrs = Counter[int](), Counter[int]()
    for run in runs:
        rbs = {r["locus_tag"]: int(r["spacer_bp"]) for r in true_rows(run.truth("rbs"))}
        for tu in tu_rows(run):
            if tu["promoter"] != "true":
                continue
            c, st = tu["contig"], tu["strand"]
            t = run.rpos(c, st, int(tu["tss"]) - 1)
            d = [k for k in range(0, 16) if run.window(c, st, t - k - 29, 29)[:6] == "TTGACA"
                 and run.window(c, st, t - k - 6, 6) == "TATAAT"]  # fmt: skip
            assert d, (run.out.name, tu["tu_id"], run.window(c, st, t - 45, 50))
            gaps10[min(d, key=lambda k: abs(k - 6))] += 1
            first = tu["genes"].split(";")[0]
            if first in rbs:
                _, _, p = run.gene_start(first)
                utrs[(p - rbs[first] - 6 - t) % len(run.seqs[c])] += 1
    assert set(gaps10) == set(range(3, 10)), sorted(gaps10.items())
    assert set(utrs) == set(range(2, 11)), sorted(utrs.items())


def test_true_promoter_layer_rows_are_the_tss(runs: list[Out]) -> None:
    for run in runs:
        want = sorted(
            (t["contig"], t["strand"], t["tss"]) for t in tu_rows(run) if t["promoter"] == "true"
        )
        got = sorted(
            (r["contig"], r["strand"], r["end3p"]) for r in true_rows(run.truth("promoters"))
        )
        assert got == want, run.out.name
