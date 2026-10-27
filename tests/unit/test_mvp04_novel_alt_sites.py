"""gippasim real-genome truth: sample SNVs inside planted novel ORFs (MVP-04 gate 11, bug 3).

``truth/alt_sites.tsv`` covers annotated genes only, yet the sample strain's SNVs also fall in
the planted ``intergenic`` / ``cds_outframe`` / ``antisense`` ORFs, whose peptides then carry
real, gene-less ALTs (5 ``high`` rows on the K-12 run). Design (fewest consumers changed:
``alt_sites.tsv`` keeps its schema and gene-only rows for the viz, compare and pilot readers):
a separate table ``truth/novel_alt_sites.tsv``, columns ``contig, codon_start, strand, orf_id,
ref_aa, alt_aa, n_true_peptides``, read by the gate 11 scorer next to ``alt_sites.tsv``.

Pinned for the three found classes (the oracle below, independent of gippasim):

* one row per codon of a planted ORF, in the ORF's frame, holding >= 1 sample SNV, except its
  start and stop codons (intact on both sides by the novel-ORF rule; the initiator reads M);
* ``codon_start`` = 1-based lowest base, ``ref_aa`` / ``alt_aa`` the code-11 residues of the
  genome / sample codon on the ORF strand; synonymous rows and rows with a stop on either side
  are left out (as in ``alt_sites.tsv``);
* ``n_true_peptides`` = distinct I/L-folded correct peptides of ``truth/peptides.tsv``
  (``locus_tag`` = ``orf_id``) whose 1-based closed interval holds the codon.

Rows of ``nterm_extension`` ORFs (if any) are not pinned beyond naming a truth ORF.
"""

from __future__ import annotations

import importlib
from dataclasses import replace
from pathlib import Path
from typing import Any

import pytest
from helpers import mvp04_data as M

from gippasim.expression import IslandSpec

SEED = 20261023
N_NOVEL = 60
COLUMNS = ["contig", "codon_start", "strand", "orf_id", "ref_aa", "alt_aa", "n_true_peptides"]
FOUND = ("intergenic", "cds_outframe", "antisense")


def _cfg() -> Any:
    base = importlib.import_module("gippasim.scenarios").k12_evidence()
    # 10x the K-12 SNV rate so the 200-kb excerpt's ~180 novel ORFs hold enough SNVs
    strain = replace(base.strain, snv_rate_per_kb=5.0)
    expr = replace(base.expression, islands=(IslandSpec(M.CONTIG, 2, 8, 3.0),))  # excerpt TUs
    return replace(base, strain=strain, expression=expr, n_novel_per_class=N_NOVEL)


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory) -> tuple[M.Excerpt, Path, Any]:
    ex = M.make_excerpt(tmp_path_factory.mktemp("k12_excerpt_nas"))
    out = tmp_path_factory.mktemp("mvp04_nas")
    rg = importlib.import_module("gippasim.real_genome")
    with M.capture_warnings():
        res = rg.simulate_on_genome(_cfg(), SEED, out, genome_fasta=ex.fasta,
                                    bakta_json=ex.bakta_json, tu_set=None,
                                    check_hashes=False)  # fmt: skip
    return ex, out, res


def _oracle(ref: str, sample: str, orfs: list[dict[str, str]],
            peps: list[dict[str, str]]) -> set[tuple[Any, ...]]:  # fmt: skip
    by: dict[str, list[dict[str, str]]] = {}
    for p in peps:
        if p["origin"] != "decoy" and p["is_wrong"] == "false":
            by.setdefault(p["locus_tag"], []).append(p)
    out = set()
    for o in orfs:
        if o["class"] not in FOUND:
            continue
        s1, e1, strand = int(o["start"]), int(o["end"]), o["strand"]
        n_cod = (e1 - s1 + 1) // 3
        for k in range(1, n_cod - 1):  # interior codons only
            c1 = s1 + 3 * k if strand == "+" else e1 - 3 * (k + 1) + 1
            g, s = M.segment(ref, c1, c1 + 2, strand), M.segment(sample, c1, c1 + 2, strand)
            if g == s:
                continue
            ra, aa = M.CODE11[g], M.CODE11[s]
            if ra == aa or "*" in (ra, aa):
                continue
            n = len({M.il(p["peptide"]) for p in by.get(o["orf_id"], ())
                     if p["contig"] == o["contig"] and int(p["start"]) <= c1
                     and c1 + 2 <= int(p["end"])})  # fmt: skip
            out.add((o["contig"], c1, strand, o["orf_id"], ra, aa, n))
    return out


def _rows(path: Path) -> list[dict[str, str]]:
    return M.read_tsv(path)


def test_novel_alt_sites_table_is_written(run: tuple[M.Excerpt, Path, Any]) -> None:
    _, out, res = run
    path = out / "truth" / "novel_alt_sites.tsv"
    assert path.is_file()
    assert M.header_of(path) == COLUMNS
    assert res.paths["novel_alt_sites"] == path


def test_planted_novel_orf_snvs_appear_in_truth(run: tuple[M.Excerpt, Path, Any]) -> None:
    ex, out, _ = run
    sample = M.read_fasta(out / "truth" / "sample_genome.fna")[M.CONTIG]
    orfs = _rows(out / "truth" / "novel_orfs.tsv")
    want = _oracle(ex.seq, sample, orfs, _rows(out / "truth" / "peptides.tsv"))
    assert len(want) >= 10, f"oracle too small to test ({len(want)} sites)"
    assert any(w[6] >= 1 for w in want), "no novel-ORF site has a covering true peptide"
    cls = {o["orf_id"]: o["class"] for o in orfs}
    got = {(r["contig"], int(r["codon_start"]), r["strand"], r["orf_id"], r["ref_aa"],
            r["alt_aa"], int(r["n_true_peptides"]))
           for r in _rows(out / "truth" / "novel_alt_sites.tsv")
           if cls.get(r["orf_id"]) in FOUND}  # fmt: skip
    assert got == want, (sorted(want - got)[:5], sorted(got - want)[:5])


def test_every_row_names_a_truth_orf_once_per_codon(run: tuple[M.Excerpt, Path, Any]) -> None:
    _, out, _ = run
    ids = {o["orf_id"] for o in _rows(out / "truth" / "novel_orfs.tsv")}
    rows = _rows(out / "truth" / "novel_alt_sites.tsv")
    assert rows and all(r["orf_id"] in ids for r in rows)
    keys = [(r["contig"], r["codon_start"], r["strand"], r["orf_id"]) for r in rows]
    assert len(keys) == len(set(keys))


def test_alt_sites_keeps_gene_rows_only(run: tuple[M.Excerpt, Path, Any]) -> None:
    """The design keeps ``alt_sites.tsv`` as it was: every row names a gene."""
    _, out, _ = run
    rows = _rows(out / "truth" / "alt_sites.tsv")
    assert rows and all(r["gene"] not in ("", "NA") for r in rows)
    assert M.header_of(out / "truth" / "alt_sites.tsv") == [
        "contig", "codon_start", "strand", "gene", "ref_aa", "alt_aa", "n_true_peptides",
        "frameshift_upstream"]  # fmt: skip


# ------------------------------------------------ novel_alt_site_rows on hand-made rows

# + strand ORF at 1..15: ATG GCC GCC GCC TAA; the sample changes codon 2 (bases 4..6) GCC -> GTC
_REF = {"c1": "ATGGCCGCCGCCTAA" + "A" * 15}
_SAMPLE = {"c1": "ATGGTCGCCGCCTAA" + "A" * 15}


def _nrow(cls: str) -> list[Any]:
    return ["NOVEL_0001", cls, "c1", 1, 15, "+", "MAAA", 13, 1, 0, "false"]


def _nas(cls: str, peps: dict[str, list[tuple[str, int, int, str]]] | None = None) -> list[Any]:
    rg = importlib.import_module("gippasim.real_genome")
    rows: list[Any] = rg.novel_alt_site_rows([_nrow(cls)], _REF, _SAMPLE, peps or {})
    return rows


@pytest.mark.parametrize("cls", FOUND)
def test_found_class_orf_snv_gives_one_row(cls: str) -> None:
    assert _nas(cls) == [["c1", 4, "+", "NOVEL_0001", "A", "V", 0]]


def test_nterm_extension_orfs_are_left_out() -> None:
    """``nterm_extension`` ORFs are annotated genes, covered by ``alt_sites.tsv``."""
    assert _nas("nterm_extension") == []


def test_n_true_peptides_counts_whole_codon_cover_1_based_closed() -> None:
    """Codon 4..6: peptides 4..6 and 1..6 hold it; 5..9 and 1..5 do not; an I/L twin and another
    contig do not add."""
    peps = {"NOVEL_0001": [("c1", 4, 6, "PEPA"), ("c1", 1, 6, "PEPB"), ("c1", 5, 9, "PEPC"),
                           ("c1", 1, 5, "PEPD"), ("c2", 1, 15, "PEPE")]}  # fmt: skip
    assert _nas("intergenic", peps)[0][6] == 2
