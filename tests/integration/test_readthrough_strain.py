"""Strain stop-loss end to end (M7-16b D5, rule 20; spec test plan "e2e, strain stop-loss";
lead decision Q20: a floor of at least 3 of 5 planted rows; statistician review of the 16b
amendment, items 1-2).

gippasim, ``related``, ``StrainConfig(snv_rate_per_kb=5, n_stop_loss=5)``, ``wrong_source =
"other_strain"`` (``pi_wrong`` 0.05), through preflight and the registered steps (``drive``, as
``run_chain`` does) to ``loci`` and ``variants``. For every ``truth/stop_loss.tsv`` row whose
``cterm_peptide`` is accepted:

- its ``peptide_hits.tsv`` row ends on the site's codon in the gene direction (``end =
  codon_start + 3`` on ``+``, ``start = codon_start`` on ``-``; 0-based, unwrapped on circular
  replicons) and is flagged ``readthrough`` (with ``strain_difference``, ``related``);
- ``variants.tsv`` has a ``stop_codon`` ALT site at ``codon_start`` with ``alt_aa =
  reference_residue`` and the ``strain_difference`` flag.

Every placed row that geometrically consumes an annotated-stop site (rule 15) must be at a site
of ``truth/stop_loss.tsv`` (planted or not) or ``truth/readthrough.tsv`` (when the build has
one); a failure is a finding (a false read-through placement), never a reason to change the
seed. The rows flagged ``readthrough`` are exactly the rows that consume a terminal site.

**Genome length (test-author measurement, reported to the lead).** The spec names sim100k. At
100 kb the cost-1 stratum holds only 27-45 target peptides whatever ``pi_wrong`` (0.05-0.3) or
the SNV rate (5-20 /kb) (seeds 11 and 12 on c523c1c0 + the 16b spec commits): with pi = 1 at
m >= 1 the smallest q is (D + 1) / T_1 >= 0.022 > alpha 0.01, so no cost-1 row, and so no planted
stop-loss peptide, is ever accepted, and raising ``pi_wrong`` or ``k`` (the spec's two remedies)
cannot reach T_1 ~ 100. At 1 Mb (seed 11, 30 s for the chain) T_1 = 400, D = 0, q = 0.0025 and
every cost-1 row is accepted, so the test runs at 1 Mb. The realised floor count needs the
implementation and is reported by the implementer's run.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.mapping_steps import drive, fresh_steps, preflight, roots_of

LENGTH = 1_000_000
SEED = 11
K = 5
FLOOR = 3


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _flags(s: str) -> set[str]:
    return set() if s in ("", "NA") else set(s.split(";"))


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    from gippa.config import resolve_config
    from gippa.mapping.sixframe import load_stop_cells
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SimConfig, simulate
    from gippasim.strain import StrainConfig

    tmp = tmp_path_factory.mktemp("rt_strain")
    res = simulate(
        SimConfig(
            genome_length=LENGTH,
            strain=StrainConfig(snv_rate_per_kb=5, n_stop_loss=K),
            proteomics=ProteomicsConfig(
                search_db="other", pi_wrong=0.05, wrong_source="other_strain"
            ),
            write_precomputed=True,
        ),
        SEED,
        tmp / "sim",
    )
    run = preflight(res, tmp / "pf")
    cfg = resolve_config(
        run.config_yaml, [f"seed={SEED}", "genome_relation=related", "sample_name=S1"]
    )
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, tmp / "work", roots_of(run), targets=("loci", "variants"))
    out = Path(res.outdir)
    sl = _tsv(out / "truth" / "stop_loss.tsv")
    sl["codon_start"] = sl["codon_start"].astype(int)
    sl["planted"] = sl["planted"].str.lower() == "true"
    rt_path = out / "truth" / "readthrough.tsv"
    from gippa.io.fasta import read_genome_fasta

    lengths = {r.id: len(r.sequence) for r in read_genome_fasta(Path(res.genome_fasta))}
    return {
        "res": res,
        "avail": avail,
        "stop_loss": sl,
        "readthrough": _tsv(rt_path) if rt_path.exists() else None,
        "hits": _tsv(avail["mapping"] / "peptide_hits.tsv"),
        "variants": _tsv(avail["variants"]),
        "placement": pd.read_parquet(avail["placement"]).drop_duplicates("peptide"),
        "sites": list(load_stop_cells(avail["sixframe"]).sites),
        "lengths": lengths,
    }


def accepted(w: dict[str, Any]) -> set[str]:
    p = w["placement"]
    return set(p.loc[p["accepted"].astype(bool), "peptide"])


def ends_on(row: Any, codon_start: int, n: int) -> bool:
    s0, e = int(row["start"]) - 1, int(row["end"])
    if row["strand"] == "+":
        return e % n == (codon_start + 3) % n
    return s0 % n == codon_start


def consumed(row: Any, sites: list[Any], n: int) -> list[Any]:
    """Rule 15 geometry on a 1-based TSV row (circular ends unwrapped)."""
    if "frameshift" in str(row["diff_types"]).split(";"):  # route is exact|tolerant only
        return []
    s, e = int(row["start"]) - 1, int(row["end"])
    if e <= s:
        e += n
    out = []
    for site in sites:
        if site.contig != row["contig"] or site.strand != row["strand"]:
            continue
        if any(
            s <= x and x + 3 <= e and (x - s) % 3 == 0
            for x in (site.codon_start, site.codon_start + n)
        ):
            out.append(site)
    return out


# ---------------------------------------------------------------- the floor


def test_at_least_three_of_five_planted_rows_are_accepted(world: dict[str, Any]) -> None:
    planted = world["stop_loss"][world["stop_loss"]["planted"]]
    assert len(planted) == K
    got = [p for p in planted["cterm_peptide"] if p in accepted(world)]
    assert len(got) >= FLOOR, (
        f"{len(got)} of {K} planted stop-loss peptides accepted: {sorted(got)}; raise pi_wrong or"
        " k and report it, never drop this assertion (spec test plan)"
    )


# ---------------------------------------------------------------- per accepted stop-loss row


def test_accepted_stop_loss_rows_are_flagged_and_end_on_the_codon(world: dict[str, Any]) -> None:
    hits, acc = world["hits"], accepted(world)
    checked = 0
    for x in world["stop_loss"].itertuples(index=False):
        if not x.cterm_peptide or x.cterm_peptide not in acc:
            continue
        rows = hits[hits["peptide"] == x.cterm_peptide]
        n = world["lengths"][x.contig]
        on = [
            r
            for _, r in rows.iterrows()
            if r["contig"] == x.contig and r["strand"] == x.strand and ends_on(r, x.codon_start, n)
        ]
        assert on, (
            x.locus_tag,
            x.cterm_peptide,
            rows[["start", "end", "strand"]].to_dict("records"),
        )
        for r in on:
            assert {"readthrough", "strain_difference"} <= _flags(r["flags"]), (x.locus_tag, r)
        checked += 1
    assert checked >= FLOOR


def test_accepted_stop_loss_rows_are_stop_codon_alt_sites(world: dict[str, Any]) -> None:
    v, acc = world["variants"], accepted(world)
    checked = 0
    for x in world["stop_loss"].itertuples(index=False):
        if not x.cterm_peptide or x.cterm_peptide not in acc:
            continue
        site = v[
            (v["contig"] == x.contig)
            & (v["strand"] == x.strand)
            & (v["codon_start"].astype(int) == x.codon_start + 1)
        ]
        alts = site[(site["call"] == "ALT") & (site["ref_aa"] == "*")]
        assert len(alts) >= 1, (x.locus_tag, site.to_dict("records"))
        assert x.reference_residue in set(alts["alt_aa"]), (x.locus_tag, alts.to_dict("records"))
        for r in alts[alts["alt_aa"] == x.reference_residue].itertuples(index=False):
            assert {"stop_codon", "strain_difference"} <= _flags(r.flags), r
            assert r.gene == x.locus_tag
        checked += 1
    assert checked >= FLOOR


# ---------------------------------------------------------------- no unexplained read-through


def test_every_consuming_row_is_at_a_true_stop_loss_or_readthrough(world: dict[str, Any]) -> None:
    truth = {(x.contig, x.strand, x.codon_start) for x in world["stop_loss"].itertuples()}
    rt = world["readthrough"]
    if rt is not None:
        truth |= {(r["contig"], r["strand"], int(r["codon_start"])) for _, r in rt.iterrows()}
    unexplained = []
    for _, r in world["hits"].iterrows():
        for site in consumed(r, world["sites"], world["lengths"][r["contig"]]):
            if (site.contig, site.strand, site.codon_start) not in truth:
                unexplained.append((r["peptide"], site.contig, site.strand, site.codon_start))
    assert not unexplained, f"false read-through placements (a finding): {unexplained}"


def test_readthrough_flag_equals_terminal_site_consumption(world: dict[str, Any]) -> None:
    terminal = [s for s in world["sites"] if s.kind == "terminal"]
    flagged, geometric = set(), set()
    for i, r in world["hits"].iterrows():
        if "readthrough" in _flags(r["flags"]):
            flagged.add(i)
        if consumed(r, terminal, world["lengths"][r["contig"]]):
            geometric.add(i)
    assert flagged == geometric
    assert len(flagged) >= FLOOR
