"""Recoded Sec codons end to end on giapsim (M7-16c; spec rules 21-23, example 32; test plan C12).

giapsim sim100k, same strain, ``n_selenoproteins = 5``, ``search_db = same``, through preflight
and the registered steps (``drive``) to ``loci`` and ``variants``
(``sixframe_index`` -> ``map_exact`` -> ``map_tolerant`` -> ``placement`` -> ``loci`` ->
``variants``):

- the six-frame index holds exactly the 5 truth codons as ``recoded`` TGA sites;
- every true peptide (>= 10 aa) of the input that spans a truth codon has a cost-1 target row
  over its own interval in ``hits_tolerant`` (the search side, independent of the placement FDR);
- every truth site whose spanning true peptide is accepted has a ``sec_site`` row in
  ``peptide_hits.tsv`` (never ``readthrough``) and an ALT ``U`` site at its codon in
  ``variants.tsv`` (``stop_codon;vcf_unrepresentable``);
- no accepted row consumes a recoded cell absent from ``truth/selenoproteins.tsv``;
- the site set ``variants`` re-aligns with equals the six-frame index's sites.

Note for the lead (test-author reading, not measured on an implementation): at 100 kb the cost-1
stratum under ``same`` is small (16b measured 27-45 target peptides at 100 kb under ``other``,
``test_readthrough_strain.py``), so the accepted-row checks may hold vacuously; the number of
accepted Sec peptides is printed and the search-side check above is the non-vacuous one.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.mapping_steps import drive, fresh_steps, preflight, roots_of

LENGTH = 100_000
SEED = 1
K = 5


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _flags(s: str) -> set[str]:
    return set() if s in ("", "NA") else set(s.split(";"))


@pytest.fixture(scope="module")
def world(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    import giap.variants.align as align_mod
    from giap.config import resolve_config
    from giap.io.fasta import read_genome_fasta
    from giap.mapping.sixframe import load_stop_cells
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SimConfig, simulate

    tmp = tmp_path_factory.mktemp("rec_sim")
    res = simulate(
        SimConfig(
            genome_length=LENGTH,
            n_selenoproteins=K,
            proteomics=ProteomicsConfig(search_db="same"),
            write_precomputed=True,
        ),
        SEED,
        tmp / "sim",
    )
    run = preflight(res, tmp / "pf")
    cfg = resolve_config(
        run.config_yaml, [f"seed={SEED}", "genome_relation=same", "sample_name=S1"]
    )
    seen: list[frozenset[tuple[str, str, int]]] = []
    real = align_mod.align_hits

    def spy(*a: Any, **k: Any) -> Any:
        seen.append(frozenset(k.get("stop_sites", frozenset())))
        return real(*a, **k)

    with pytest.MonkeyPatch.context() as mp:
        mp.setattr(align_mod, "align_hits", spy)
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, tmp / "work", roots_of(run), targets=("loci", "variants"))
    out = Path(res.outdir)
    sec = _tsv(out / "truth" / "selenoproteins.tsv")
    sec["codon_start"] = sec["codon_start"].astype(int)
    lengths = {r.id: len(r.sequence) for r in read_genome_fasta(Path(res.genome_fasta))}
    return {
        "res": res,
        "run": run,
        "avail": avail,
        "sec": sec,
        "truth_peptides": _tsv(out / "truth" / "peptides.tsv"),
        "inputs": set(pd.read_parquet(run.peptides)["peptide"]),
        "hits_tolerant": pd.read_parquet(avail["hits_tolerant"]),
        "hits": _tsv(avail["mapping"] / "peptide_hits.tsv"),
        "variants": _tsv(avail["variants"]),
        "placement": pd.read_parquet(avail["placement"]).drop_duplicates("peptide"),
        "sites": list(load_stop_cells(avail["sixframe"]).sites),
        "lengths": lengths,
        "variant_sets": seen,
    }


def accepted(w: dict[str, Any]) -> set[str]:
    p = w["placement"]
    return set(p.loc[p["accepted"].astype(bool), "peptide"])


def spanning(w: dict[str, Any], x: Any) -> list[tuple[str, int, int]]:
    """(peptide, 0-based start, end) of the true gene peptides (>= 10 aa) over a truth codon."""
    t = w["truth_peptides"]
    t = t[(t["is_wrong"] == "false") & (t["origin"] == "gene") & (t["locus_tag"] == x.locus_tag)]
    out = []
    for r in t.itertuples(index=False):
        s0, e = int(r.start) - 1, int(r.end)
        if len(r.peptide) >= 10 and s0 <= x.codon_start and x.codon_start + 3 <= e:
            out.append((r.peptide, s0, e))
    return out


def consumed(row: Any, sites: list[Any], n: int) -> list[Any]:
    """Rule 15 geometry on a 1-based TSV row (circular ends unwrapped)."""
    if "frameshift" in str(row["diff_types"]).split(";"):
        return []
    s, e = int(row["start"]) - 1, int(row["end"])
    if e <= s:
        e += n
    return [
        site
        for site in sites
        if site.contig == row["contig"]
        and site.strand == row["strand"]
        and any(
            s <= x and x + 3 <= e and (x - s) % 3 == 0
            for x in (site.codon_start, site.codon_start + n)
        )
    ]


# ---------------------------------------------------------------- the index


def test_the_index_holds_the_truth_codons_as_recoded_sites(world: dict[str, Any]) -> None:
    rec = {(s.contig, s.strand, s.codon_start): s for s in world["sites"] if s.kind == "recoded"}
    want = {(x.contig, x.strand, x.codon_start) for x in world["sec"].itertuples(index=False)}
    assert len(want) == K and set(rec) == want
    assert {s.codon for s in rec.values()} == {"TGA"}
    for x in world["sec"].itertuples(index=False):
        assert x.locus_tag in rec[(x.contig, x.strand, x.codon_start)].locus_tags


def test_variants_realigns_with_the_index_site_set(world: dict[str, Any]) -> None:
    assert world["variant_sets"], "align_hits was not called by the variants step"
    keys = frozenset((s.contig, s.strand, s.codon_start) for s in world["sites"])
    assert all(v == keys for v in world["variant_sets"])


# ---------------------------------------------------------------- the search side


def test_spanning_true_peptides_get_a_cost_1_target_row(world: dict[str, Any]) -> None:
    ht = world["hits_tolerant"]
    ht = ht[ht["search"] == "target"]
    checked = 0
    for x in world["sec"].itertuples(index=False):
        for pep, s0, e in spanning(world, x):
            if pep not in world["inputs"]:
                continue
            rows = ht[
                (ht["peptide"] == pep) & (ht["contig"] == x.contig) & (ht["strand"] == x.strand)
            ]
            n = world["lengths"][x.contig]
            here = rows[
                (rows["start"].astype(int) % n == s0 % n)
                & (rows["end"].astype(int) - rows["start"].astype(int) == e - s0)
            ]
            assert len(here) == 1, (x.locus_tag, pep, rows.to_dict("records"))
            assert int(here.iloc[0]["n_diff"]) == 1 and here.iloc[0]["diff_types"] == "substitution"
            checked += 1
    assert checked >= 1, "no input peptide spans a planted Sec codon"


# ---------------------------------------------------------------- accepted rows


def test_accepted_spanning_peptides_have_sec_site_rows_and_u_alt_sites(
    world: dict[str, Any],
) -> None:
    acc = accepted(world)
    hits, v = world["hits"], world["variants"]
    n_acc = 0
    for x in world["sec"].itertuples(index=False):
        peps = [p for p, _s, _e in spanning(world, x) if p in acc]
        if not peps:
            continue
        n_acc += len(peps)
        for pep in peps:
            rows = hits[(hits["peptide"] == pep) & (hits["contig"] == x.contig)]
            n = world["lengths"][x.contig]
            over = [
                r
                for _, r in rows.iterrows()
                if any(s.codon_start == x.codon_start for s in consumed(r, world["sites"], n))
            ]
            assert over, (x.locus_tag, pep, rows.to_dict("records"))
            for r in over:
                assert "sec_site" in _flags(r["flags"]) and "readthrough" not in _flags(r["flags"])
        site = v[
            (v["contig"] == x.contig)
            & (v["strand"] == x.strand)
            & (v["codon_start"].astype(int) == x.codon_start + 1)
        ]
        u = site[(site["call"] == "ALT") & (site["alt_aa"] == "U")]
        assert len(u) >= 1, (x.locus_tag, site.to_dict("records"))
        for r in u.itertuples(index=False):
            assert {"stop_codon", "vcf_unrepresentable"} <= _flags(r.flags)
            assert r.ref_aa == "*" and r.gene == x.locus_tag
    print(f"accepted Sec-spanning peptides: {n_acc}")


def test_no_accepted_row_consumes_a_recoded_cell_outside_the_truth(world: dict[str, Any]) -> None:
    truth = {(x.contig, x.strand, x.codon_start) for x in world["sec"].itertuples(index=False)}
    rec = [s for s in world["sites"] if s.kind == "recoded"]
    for _, r in world["hits"].iterrows():
        n = world["lengths"][r["contig"]]
        for s in consumed(r, rec, n):
            assert (s.contig, s.strand, s.codon_start) in truth, (r.to_dict(), s)


def test_sec_site_rows_are_exactly_the_rows_consuming_a_recoded_tga(world: dict[str, Any]) -> None:
    rec = [s for s in world["sites"] if s.kind == "recoded"]
    for _, r in world["hits"].iterrows():
        n = world["lengths"][r["contig"]]
        has = any(s.codon == "TGA" for s in consumed(r, rec, n))
        assert ("sec_site" in _flags(r["flags"])) == has, r.to_dict()
        assert "pyl_site" not in _flags(r["flags"])
