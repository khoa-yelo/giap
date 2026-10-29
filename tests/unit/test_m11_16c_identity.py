"""M11-16c: gippasim at genetic code 11 is byte-identical to the simulator before M11-16.

docs/specs/M11-16-translation-tables-4-25.md item 8 ("With 11, every output is byte-identical to
today for every existing scenario and seed. This is the first test of PR 16c") and AC 14
(``SimConfig()`` output, all existing scenarios, seeds 1-3, byte-identical to main).

The golden digests in ``tests/data/m11_16c/table11_digests.json`` were produced from the
simulator at main 1f7203d3 (branch point of ``spec/translation-tables``, before any 16c change)
with the ``_digests`` / ``_layout_digests`` / ``_fixture_digest`` functions of this module. They
pin the bytes of every file ``gippasim.simulate`` writes (the output directory path, which
``config.yaml`` embeds, is replaced by ``<OUT>``). **They must not be regenerated in 16c**: a
mismatch is a 16c regression, not a stale golden file.

Every ``simulate`` case passes ``genetic_code=11`` explicitly (``dataclasses.replace``), so it
also checks that the explicit default is the implicit one. The layout pins run twice: with the
implicit default and with ``LayoutConfig(genetic_code=11)`` (spec item 8 lists ``layout`` among
the modules that take the table from their config; lead decision on #474). The M11-14w fixture
pins (codes 4 and 11) are pure regression guards and pass before 16c.
"""

from __future__ import annotations

import hashlib
import json
import sys
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np
import pytest

GOLDEN = Path(__file__).resolve().parents[1] / "data" / "m11_16c" / "table11_digests.json"
SEEDS = (1, 2, 3)
CASES = (
    "m0_default",
    "m3_18_features",
    "m7_16_sec_readthrough",
    "m7_16_stop_loss",
    "m2_05_sim100k",
    "m7_12_variant_f1",
    "m7_12_novel_orf_power",
    "m7_12_same_strain_insilico",
    "m7_12_isobaric_attack",
)


def _cases() -> dict[str, Callable[[], Any]]:
    """Every existing gippasim scenario family (M0, M3-18, M7-16, M2-05, M7-12)."""
    from gippasim import scenarios
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SIM100K, SimConfig
    from gippasim.strain import StrainConfig

    return {
        "m0_default": lambda: SimConfig(),
        "m3_18_features": lambda: SimConfig(
            genome_length=60_000,
            n_plasmids=1,
            plasmid_length=(6_000, 6_000),
            linear_contig_length=5_000,
            origin_crossing_gene=True,
            n_rrna_operons=1,
            n_trna=3,
            n_pseudogenes=2,
            n_sorfs=2,
            n_nterm_extensions=2,
            n_novel_orfs=(("intergenic", 2), ("cds_outframe", 2), ("antisense", 2)),
            strain=StrainConfig(),
        ),
        "m7_16_sec_readthrough": lambda: SimConfig(
            n_selenoproteins=2, n_readthrough=2, proteomics=ProteomicsConfig()
        ),
        "m7_16_stop_loss": lambda: SimConfig(
            strain=StrainConfig(snv_rate_per_kb=2.0, n_stop_loss=2), proteomics=ProteomicsConfig()
        ),
        "m2_05_sim100k": lambda: SIM100K,
        "m7_12_variant_f1": lambda: scenarios.variant_f1(0),
        "m7_12_novel_orf_power": lambda: scenarios.novel_orf_power(0),
        "m7_12_same_strain_insilico": lambda: scenarios.same_strain_insilico(0),
        "m7_12_isobaric_attack": lambda: scenarios.isobaric_attack(0),
    }


def _digests(outdir: Path) -> dict[str, str]:
    """sha256 of every file under ``outdir`` (relative POSIX path -> hex), path-normalised."""
    root = str(outdir).encode()
    out: dict[str, str] = {}
    for p in sorted(q for q in outdir.rglob("*") if q.is_file()):
        data = p.read_bytes().replace(root, b"<OUT>")
        out[p.relative_to(outdir).as_posix()] = hashlib.sha256(data).hexdigest()
    return out


def _layout_digests(layout: Any) -> dict[str, str]:
    """sha256 of each frame of a ``gippasim.layout.Layout`` (CSV text)."""
    names = ("replicons", "genes", "tus", "operon_pairs", "regions", "region_genes", "truth")
    names += ("xrefs", "features")
    return {n: hashlib.sha256(getattr(layout, n).to_csv().encode()).hexdigest() for n in names}


def _fixture_digest(seed: int, gc: float, code: int) -> str:
    from helpers.genetic_code_fixture import synthetic_genome

    return hashlib.sha256(synthetic_genome(seed, gc, code).encode()).hexdigest()


def _golden() -> dict[str, Any]:
    return json.loads(GOLDEN.read_text())


def test_default_genetic_code_is_11() -> None:
    from gippasim.simulate import SimConfig

    assert SimConfig().genetic_code == 11


@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("case", CASES)
def test_simulate_table11_byte_identical(case: str, seed: int, tmp_path: Path) -> None:
    from gippasim.simulate import simulate

    cfg = replace(_cases()[case](), genetic_code=11)
    simulate(cfg, seed, tmp_path)
    got = _digests(tmp_path)
    want = _golden()["simulate"][case][str(seed)]
    assert sorted(got) == sorted(want), "file set changed at genetic_code=11"
    changed = sorted(k for k in want if got[k] != want[k])
    assert not changed, f"{case} seed {seed}: bytes changed at genetic_code=11: {changed}"


def test_golden_covers_every_case() -> None:
    assert sorted(_golden()["simulate"]) == sorted(_cases()) == sorted(CASES)


# operon_pairs.uniop_p is an rng.beta draw: its last bits follow the platform's libm, so the
# golden CSV bytes (made on Linux) differ on macOS. The pin guards code changes; Linux CI runs it.
@pytest.mark.skipif(sys.platform != "linux", reason="golden layout bytes are Linux libm (rng.beta)")
@pytest.mark.parametrize("explicit", [False, True], ids=["implicit", "genetic_code_11"])
@pytest.mark.parametrize("seed", SEEDS)
@pytest.mark.parametrize("name", ["default", "draft_origin"])
def test_layout_byte_identical(name: str, seed: int, explicit: bool) -> None:
    """``simulate_layout`` at code 11 (implicit, and ``genetic_code=11``) is unchanged."""
    from gippasim.layout import LayoutConfig, simulate_layout

    cfg = {
        "default": LayoutConfig(n_genes_chromosome=300),
        "draft_origin": LayoutConfig(
            n_genes_chromosome=300,
            plasmids=(20,),
            linear_contigs=(10, 10),
            origin_crossing=True,
            oric_frac=0.8,
        ),
    }[name]
    if explicit:
        cfg = replace(cfg, genetic_code=11)
    got = _layout_digests(simulate_layout(cfg, np.random.default_rng(seed)))
    assert got == _golden()["layout"][name][str(seed)]


@pytest.mark.parametrize(
    ("seed", "gc", "code"),
    [(s, gc, c) for s in (1, 7) for gc in (0.28, 0.5) for c in (4, 11)],
)
def test_m11_14w_fixture_code4_and_11_unchanged(seed: int, gc: float, code: int) -> None:
    """Item 8: "The code-4 path keeps its seeds and values" (and code 11 with it)."""
    key = f"{seed}:{gc}:{code}"
    assert _fixture_digest(seed, gc, code) == _golden()["m11_14w_fixture"][key]
