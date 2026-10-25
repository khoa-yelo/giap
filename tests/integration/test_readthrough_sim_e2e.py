"""Planted stop-codon read-through end to end (M7-16 rule 14 and the 16b acceptance item "On
sim100k with ``n_readthrough = 10``, every planted gene with >= 2 distinct spanning true peptides
gives a ``called`` ``cterm_extension`` with ``readthrough_candidate``. Its patch row round-trips
through ``makedb --patch`` with ``X`` at the old stop. Each planted stop with a placed spanning
peptide is a ``stop_codon`` ALT site with the planted residue"; test plan 16b "an e2e local run on
sim100k"; PR 16b-3 per the lead decision in the #335 review round).

Build: giapsim ``SIM100K`` with M3-19 proteomics, M3-20 de novo reads and every precomputed layer
(as ``test_m7_04_variants_run.py``), ``n_readthrough = 10``; ``giap run --executor local
--modules core`` on the **de novo** table (``genome_relation`` default), then ``giap makedb
--patch`` on the run's ``annotation_patch.gff3``.

**Why de novo input (test-author measurement, reported to the lead; main 1549a930, seed
20261014).** Read-through rows are cost-1 tolerant rows (pi = 1 at m >= 1), so the smallest
stratum-1 q is (D + 1) / T_1. With the database table at 100 kb T_1 = 7 (q >= 0.14 > alpha 0.01)
and no read-through row could ever be accepted: the rule-20 floor problem. With the de novo table
the de novo errors fill stratum 1: T_1 = 1361, all accepted, q_min = 7.3e-4, so an error-free
read-through read is accepted at sim100k and the spec's genome size stands (no 1 Mb needed). A
spanning read that carries a de novo error as well is cost 2 (T_2 = 54, q_min 0.019) and is not
accepted, so the "true peptides" of the acceptance item are read here as the spanning true
peptides (>= 10 aa, shorter ones are never placed) that have an error-free de novo read
(``truth/denovo.tsv``). The run takes about 35 s, so the test is not ``slow``.

**Floors** (the test can never pass vacuously): the oracle expects about 3 of 10 planted genes
with >= 2 detected spanning peptides and about 6.5 with >= 1 (detection logit of M3-19 on the 59
eligible genes of this build). The test requires >= 1 called gene and >= 3 stop_codon ALT
sites; a shortfall is reported and answered by raising ``n_readthrough`` (59 eligible genes),
never by dropping an assertion or changing the seed.

Every placed row that consumes an annotated terminal stop (rule 15 geometry) must be at a planted
site: a failure is a finding (a false read-through placement), never a reason to change the seed.
"""

from __future__ import annotations

import dataclasses
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from click.testing import CliRunner
from helpers.readthrough_sim_oracle import (
    annotated_proteins,
    eligible,
    extension,
    read_genome,
    readthrough_table,
    spanning_peptides,
    tsv,
)

N_READTHROUGH = 10
RUN_SEED = 7
FLOOR_CALLED = 1
FLOOR_ALT_SITES = 3
VCF_ALT = {("TGA", "W"): "TGG", ("TAA", "Q"): "CAA", ("TAG", "Q"): "CAG"}


def _il(s: str) -> str:
    return s.replace("I", "L")


def _flags(s: str) -> set[str]:
    return set() if s in ("", "NA") else set(s.split(";"))


def _cli(args: list[str], db: Path) -> Any:
    from giap.cli import cli

    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIAP_DB_DIR", str(db))
        mp.setenv("NXF_ANSI_LOG", "false")
        res = CliRunner().invoke(cli, args)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    assert res.exit_code == 0, res.output[-4000:]
    return res


def _fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    name = ""
    for ln in path.read_text().splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            out[name] = []
        elif ln:
            out[name].append(ln.strip())
    return {k: "".join(v) for k, v in out.items()}


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    from helpers.mapping_steps import write_test_db

    import giapsim
    from giap.mapping.sixframe import load_stop_cells
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig

    tmp = tmp_path_factory.mktemp("rt_e2e")
    cfg = dataclasses.replace(
        giapsim.SIM100K,
        proteomics=ProteomicsConfig(),
        denovo=DenovoConfig(),
        write_precomputed=True,
        n_readthrough=N_READTHROUGH,  # type: ignore[call-arg]
    )
    sim = giapsim.simulate(cfg, giapsim.SIM100K_SEED, tmp / "sim")
    zero = giapsim.simulate(
        dataclasses.replace(cfg, n_readthrough=0),  # type: ignore[call-arg]
        giapsim.SIM100K_SEED,
        tmp / "sim0",
    )
    db = write_test_db(tmp / "test_db")
    out = tmp / "OUT"
    _cli(
        ["run", "--genome", str(sim.genome_fasta),
         "--peptides", str(Path(sim.outdir) / "peptides_denovo.tsv"),
         "--annotation", str(sim.bakta_json), "--config", str(sim.config_yaml),
         "--outdir", str(out), "--modules", "core", "--seed", str(RUN_SEED),
         "--executor", "local"],
        db,
    )  # fmt: skip
    patch = out / "tables" / "annotation_patch.gff3"
    mdb = tmp / "DB"
    _cli(
        ["makedb", "--annotation", str(sim.bakta_json), "--genome", str(sim.genome_fasta),
         "--outdir", str(mdb), "--patch", str(patch), "--no-contaminants"],
        db,
    )  # fmt: skip
    simdir = Path(sim.outdir)
    return {
        "sim": simdir,
        "genome": read_genome(simdir),
        "rt": readthrough_table(simdir),
        "eligible": eligible(Path(zero.outdir)),
        "hits": tsv(out / "mapping" / "peptide_hits.tsv"),
        "variants": tsv(out / "tables" / "variants.tsv"),
        "novel": tsv(out / "tables" / "novel_orfs.tsv"),
        "patch": patch,
        "strain_db": _fasta(mdb / "strain_db.faa"),
        "sites": list(load_stop_cells(out / "work" / "sixframe").sites),
    }


def _interval(row: Any, n: int) -> tuple[int, int]:
    s, e = int(row["start"]) - 1, int(row["end"])
    return s, (e + n if e <= s else e)


def _cell(row: Any, codon_start: int, n: int, circular: bool) -> int | None:
    """Residue index in the row's peptide of the codon at ``codon_start``; None if the row does
    not hold it in frame (rule 15; frameshift rows never consume)."""
    if "frameshift" in _flags(row["diff_types"]):
        return None
    s, e = _interval(row, n)
    for x in (codon_start, codon_start + n) if circular else (codon_start,):
        if s <= x and x + 3 <= e and (x - s) % 3 == 0:
            return (x - s) // 3 if row["strand"] == "+" else (e - x - 3) // 3
    return None


def _consuming(w: dict[str, Any], contig: str, strand: str, codon_start: int) -> pd.DataFrame:
    g, hits = w["genome"], w["hits"]
    rows = hits[(hits["contig"] == contig) & (hits["strand"] == strand)]
    keep = [
        i
        for i, r in rows.iterrows()
        if _cell(r, codon_start, g.n(contig), g.circular[contig]) is not None
    ]
    return rows.loc[keep]


def _three_prime_end(x: Any) -> int:
    """1-based last base of the new stop codon in transcription direction (novel_orfs ``end``)."""
    return x.ext_stop + 3 if x.strand == "+" else x.ext_stop + 1


def _called_genes(w: dict[str, Any]) -> set[str]:
    nv = w["novel"]
    sel = (nv["class"] == "cterm_extension") & (nv["called"].str.lower() == "true")
    return set(nv.loc[sel, "overlapping_gene"])


# ---------------------------------------------------------------- the build


def test_the_build_plants_ten_eligible_genes(world: dict[str, Any]) -> None:
    rt = world["rt"]
    assert len(rt) == N_READTHROUGH
    assert set(rt["locus_tag"]) <= set(world["eligible"])


# ---------------------------------------------------------------- placement


def test_every_consuming_row_is_at_a_planted_site(world: dict[str, Any]) -> None:
    planted = {(x.contig, x.strand, x.codon_start) for x in world["rt"].itertuples(index=False)}
    unexplained = []
    for site in world["sites"]:
        if site.kind != "terminal":
            continue
        rows = _consuming(world, site.contig, site.strand, site.codon_start)
        if len(rows) and (site.contig, site.strand, site.codon_start) not in planted:
            unexplained += [
                (p, site.contig, site.strand, site.codon_start) for p in rows["peptide"]
            ]
    assert not unexplained, f"false read-through placements (a finding): {unexplained}"


def test_rows_over_a_planted_stop_are_flagged_at_their_gene(world: dict[str, Any]) -> None:
    n_rows = 0
    for x in world["rt"].itertuples(index=False):
        for _, r in _consuming(world, x.contig, x.strand, x.codon_start).iterrows():
            assert "readthrough" in _flags(r["flags"]), (x.locus_tag, r["peptide"])
            assert r["gene"] == x.locus_tag, (x.locus_tag, r["peptide"], r["gene"])
            assert r["locus_class"] in ("cterm_extension", "cds_inframe"), r["locus_class"]
            assert r["route"] == "tolerant" and int(r["n_diff"]) >= 1
            n_rows += 1
    assert n_rows >= FLOOR_ALT_SITES


# ---------------------------------------------------------------- variants


def test_placed_spanning_peptides_are_stop_codon_alt_sites(world: dict[str, Any]) -> None:
    g, v = world["genome"], world["variants"]
    sites = 0
    for x in world["rt"].itertuples(index=False):
        rows = _consuming(world, x.contig, x.strand, x.codon_start)
        n, circ = g.n(x.contig), g.circular[x.contig]
        read = [
            r for _, r in rows.iterrows()
            if r["peptide"][_cell(r, x.codon_start, n, circ)] == x.residue  # type: ignore[index]
        ]  # fmt: skip
        if not read:
            continue
        at = v[
            (v["contig"] == x.contig)
            & (v["strand"] == x.strand)
            & (v["codon_start"].astype(int) == x.codon_start + 1)
            & (v["ref_aa"] == "*")
            & (v["alt_aa"] == x.residue)
        ]
        assert len(at) == 1, (x.locus_tag, v[v["gene"] == x.locus_tag].to_dict("records"))
        (r,) = at.to_dict("records")
        assert r["call"] in ("ALT", "mixed"), r
        assert r["gene"] == x.locus_tag and r["ref_codon"] == x.codon, r
        assert {"snv", "stop_codon"} <= _flags(r["flags"]), r
        assert "multi_nt" not in _flags(r["flags"]), r
        assert r["alt_codon_vcf"] == VCF_ALT[(x.codon, x.residue)], r
        sites += 1
    assert sites >= FLOOR_ALT_SITES, (
        f"{sites} planted stops with a placed spanning read; raise n_readthrough and report it"
    )


# ---------------------------------------------------------------- novel ORFs


def test_genes_with_two_cterm_peptides_are_called(world: dict[str, Any]) -> None:
    hits, nv = world["hits"], world["novel"]
    called = 0
    for x in world["rt"].itertuples(index=False):
        mine = hits[(hits["locus_class"] == "cterm_extension") & (hits["gene"] == x.locus_tag)]
        if mine["peptide"].map(_il).nunique() < 2:
            continue
        rows = nv[(nv["class"] == "cterm_extension") & (nv["overlapping_gene"] == x.locus_tag)]
        assert len(rows) == 1, (x.locus_tag, sorted(mine["peptide"]))
        (r,) = rows.to_dict("records")
        assert r["called"].lower() == "true", r
        assert "readthrough_candidate" in _flags(r["flags"]), r
        assert (r["contig"], r["strand"]) == (x.contig, x.strand)
        assert int(r["end"]) == _three_prime_end(x), (x.locus_tag, r["end"], x.ext_stop)
        called += 1
    assert called >= FLOOR_CALLED, (
        f"{called} planted genes with >= 2 cterm_extension peptides; raise n_readthrough and "
        "report it, never drop this assertion"
    )


def test_genes_with_two_readable_spanning_true_peptides_are_called(world: dict[str, Any]) -> None:
    """The acceptance item on the truth side: >= 2 distinct spanning true peptides of >= 10 aa
    with an error-free de novo read give a called cterm_extension."""
    sim = world["sim"]
    dn = tsv(sim / "truth" / "denovo.tsv")
    clean = set(dn.loc[dn["peptide"] == dn["true_peptide"], "true_peptide"])
    want = set()
    for x in world["rt"].itertuples(index=False):
        peps = spanning_peptides(sim, x.locus_tag, x.codon_start, x.contig)
        if len({_il(p) for p in peps if len(p) >= 10 and p in clean}) >= 2:
            want.add(x.locus_tag)
    assert want, "no planted gene has two readable spanning peptides; raise n_readthrough"
    assert want <= _called_genes(world), sorted(want - _called_genes(world))


def test_every_called_cterm_extension_is_a_planted_gene(world: dict[str, Any]) -> None:
    assert _called_genes(world) <= set(world["rt"]["locus_tag"])


# ---------------------------------------------------------------- patch round trip


def test_patch_rows_round_trip_through_makedb_with_x(world: dict[str, Any]) -> None:
    called = _called_genes(world)
    rows = [
        ln.split("\t")
        for ln in world["patch"].read_text().splitlines()
        if ln and not ln.startswith("#") and ln.split("\t")[2] == "CDS"
    ]
    cterm = {
        dict(kv.split("=", 1) for kv in r[8].split(";"))["locus_tag"]: r
        for r in rows
        if "Note=cterm_extension" in r[8].split(";")
    }
    assert set(cterm) == called
    rt = world["rt"].set_index("locus_tag")
    ann = annotated_proteins(world["sim"])
    g = world["genome"]
    for tag in sorted(called):
        x = rt.loc[tag]
        e = world["eligible"][tag]
        want = "M" + ann[tag][1:] + "X" + extension(g, e)
        assert world["strain_db"][f"{tag}_patch"] == want, tag
        assert world["strain_db"][tag] == ann[tag], "the Bakta protein is still written"
        r = cterm[tag]
        assert r[6] == x["strand"]
        three = int(r[4]) if x["strand"] == "+" else int(r[3])
        n = g.n(x["contig"])
        assert (three - 1) % n + 1 == (int(x["ext_stop"]) + 3 if x["strand"] == "+" else
                                       int(x["ext_stop"]) + 1), (tag, r)  # fmt: skip
