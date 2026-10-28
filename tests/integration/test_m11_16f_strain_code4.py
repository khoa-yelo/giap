"""AC 12, second half (M11-16f; moved here from the 16e review): "a code-4 sim strain run has no
``stop_codon`` flag on TGA -> TGG sites" (docs/specs/M11-16-translation-tables-4-25.md example 12,
AC 12).

gippasim at ``genetic_code = 4``, ``related``, ``StrainConfig(snv_rate_per_kb=5)``, other-strain
wrong peptides (the ``test_readthrough_strain.py`` world without planted stop-losses), 1 Mb, through
preflight and the registered steps to ``variants`` with ``annotation.translation_table = 4``.

- The strain SNVs that turn a sample TGA into a reference TGG (or back) exist (non-vacuous): at
  code 4 both read W, so they are synonymous.
- No ``variants.tsv`` row exists at such a codon, or every row there has ``ref_aa = W``, an
  ``alt_aa`` other than W and no ``stop_codon`` flag (statistician review of #481, N4; the row
  count is logged). No row at any TGA codon has ``ref_aa = *`` or the ``stop_codon`` flag: TGA is
  read as W (>= 5 rows at TGA codons exist and say W).
- Control: the same simulated world and preflight driven at ``annotation.translation_table = 11``
  gives at least one ``stop_codon`` row at a TGA codon, so the flag the checks above look for is
  one this world can produce.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.mapping_steps import drive, fresh_steps, preflight, roots_of

LENGTH = 1_000_000
SEED = 11
TABLE = 4
_RC = str.maketrans("ACGT", "TGCA")


def _tsv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def _flags(s: str) -> set[str]:
    return set() if s in ("", "NA") else set(s.split(";"))


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Any]:
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SimConfig, simulate
    from gippasim.strain import StrainConfig

    tmp = tmp_path_factory.mktemp("strain_code4")
    res = simulate(
        SimConfig(
            genome_length=LENGTH,
            genetic_code=TABLE,
            strain=StrainConfig(snv_rate_per_kb=5),
            proteomics=ProteomicsConfig(
                search_db="other", pi_wrong=0.05, wrong_source="other_strain"
            ),
            write_precomputed=True,
        ),
        SEED,
        tmp / "sim",
    )
    return {"res": res, "run": preflight(res, tmp / "pf"), "tmp": tmp}


def _variants(sim: dict[str, Any], table: int) -> pd.DataFrame:
    """``variants.tsv`` of the simulated world driven at ``annotation.translation_table``."""
    from gippa.config import resolve_config

    cfg = resolve_config(
        sim["run"].config_yaml,
        [
            f"seed={SEED}",
            "genome_relation=related",
            "sample_name=S1",
            f"annotation.translation_table={table}",
        ],
    )
    assert cfg.annotation.translation_table == table
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, sim["tmp"] / f"work{table}", roots_of(sim["run"]),
                      targets=("loci", "variants"))  # fmt: skip
    return _tsv(avail["variants"])


@pytest.fixture(scope="module")
def world(sim: dict[str, Any]) -> dict[str, Any]:
    from gippa.io.fasta import read_genome_fasta
    from gippa.model import Replicon
    from gippa.model.genes import gene_model_from_bakta_json

    res = sim["res"]
    out = Path(res.outdir)
    sample = {r.id: r.sequence for r in read_genome_fasta(Path(res.genome_fasta))}
    reference = {r.id: r.sequence for r in read_genome_fasta(out / "reference.fasta")}
    reps = {c: Replicon(c, len(s), "circular") for c, s in sample.items()}
    model = gene_model_from_bakta_json(res.bakta_json, reps)
    return {
        "variants": _variants(sim, TABLE),
        "diffs": _tsv(out / "truth" / "strain_differences.tsv"),
        "sample": sample,
        "reference": reference,
        "genes": {g.locus_tag: g for g in model.genes},
    }


def _codon(seq: str, start: int, strand: str) -> str:
    nt = seq[start : start + 3]
    return nt if strand == "+" else nt.translate(_RC)[::-1]


def tga_tgg_sites(w: dict[str, Any]) -> set[tuple[str, int, str]]:
    """(contig, 0-based codon start, strand) of the strain SNVs TGA <-> TGG in a CDS frame."""
    out = set()
    for d in w["diffs"].itertuples(index=False):
        if d.type != "snv" or d.gene not in w["genes"]:
            continue
        g = w["genes"][d.gene]
        iv, pos = g.interval, int(d.pos_sample) - 1  # 1-based; ``ref`` is the reference strain's
        assert pos == int(d.pos_reference) - 1, d  # SNV-only builds: coordinates agree
        assert (w["reference"][d.contig][pos], w["sample"][d.contig][pos]) == (d.ref, d.alt), d
        if not iv.start <= pos < iv.end:
            continue
        cs = (
            pos - (pos - iv.start) % 3
            if iv.strand == "+"
            else iv.end - 3 * ((iv.end - 1 - pos) // 3) - 3
        )
        pair = {_codon(w["sample"][d.contig], cs, iv.strand),
                _codon(w["reference"][d.contig], cs, iv.strand)}  # fmt: skip
        if pair == {"TGA", "TGG"}:
            out.add((d.contig, cs, iv.strand))
    return out


def test_tga_tgg_strain_sites_exist(world: dict[str, Any]) -> None:
    assert len(tga_tgg_sites(world)) >= 5


def test_no_stop_codon_flag_on_tga_tgg_sites(world: dict[str, Any]) -> None:
    v = world["variants"]
    sites = {(c, s + 1, st) for c, s, st in tga_tgg_sites(world)}  # variants.tsv is 1-based
    at = v[[(r.contig, int(r.codon_start), r.strand) in sites for r in v.itertuples()]]
    print(f"variant rows at the {len(sites)} TGA<->TGG strain codons: len(at) = {len(at)}")
    bad = [r for r in at.itertuples()
           if "stop_codon" in _flags(r.flags) or r.ref_aa != "W" or r.alt_aa == "W"]  # fmt: skip
    assert not bad, (len(at), bad[:5])


def test_tga_is_read_as_w_in_every_variant_row(world: dict[str, Any]) -> None:
    v = world["variants"]
    tga = v[v["ref_codon"] == "TGA"]
    assert len(tga) >= 5, f"{len(tga)} variant rows at TGA codons (< 5): the check is weak"
    assert set(tga["ref_aa"]) == {"W"}, tga["ref_aa"].value_counts().to_dict()
    flagged = [r for r in tga.itertuples() if "stop_codon" in _flags(r.flags)]
    assert not flagged, flagged[:5]


def test_control_table_11_flags_stop_codon_at_tga(sim: dict[str, Any]) -> None:
    """The same world at table 11: TGA is a stop there, so the flag can occur (positive control)."""
    v = _variants(sim, 11)
    tga = v[v["ref_codon"] == "TGA"]
    flagged = [r for r in tga.itertuples() if "stop_codon" in _flags(r.flags)]
    print(f"table 11 control: {len(tga)} rows at TGA codons, {len(flagged)} with stop_codon")
    assert flagged, f"no stop_codon row at a TGA codon at table 11 ({len(tga)} TGA rows)"
