"""The giapsim ``SIM5M`` scenario (docs/specs/M10-12-benchmark-standard.md, criterion 8).

``SIM5M`` is implementer tooling (AGENTS.md: the implementer writes giapsim, the test-author tests
it), so only its contract is tested: the name and seed are exported next to ``SIM100K``, the
options are ``SIM100K``'s with the spec's overrides, the real tools must run
(``write_precomputed=False``), and two ``simulate(SIM5M, SIM5M_SEED)`` calls are byte-identical
(the summary's ``inputs_differ`` check relies on it). The determinism test simulates a 5-Mb
genome twice and is marked ``slow``.
"""

from __future__ import annotations

import dataclasses
import importlib
from pathlib import Path

import pytest

SEED = 20261014


def _sim():  # type: ignore[no-untyped-def]
    mod = importlib.import_module("giapsim.simulate")
    assert hasattr(mod, "SIM5M"), "giapsim.simulate.SIM5M does not exist (M10-12)"
    return mod


def test_sim5m_exported_next_to_sim100k() -> None:
    mod = _sim()
    import giapsim

    for name in ("SIM5M", "SIM5M_SEED"):
        assert name in mod.__all__, f"giapsim.simulate.__all__ lacks {name}"
        assert name in giapsim.__all__, f"giapsim.__all__ lacks {name}"
    assert giapsim.SIM5M is mod.SIM5M and giapsim.SIM5M_SEED == mod.SIM5M_SEED


def test_sim5m_seed() -> None:
    assert _sim().SIM5M_SEED == SEED


def test_sim5m_is_sim100k_with_the_spec_overrides() -> None:
    mod = _sim()
    from giapsim.proteomics import ProteomicsConfig

    expected = dataclasses.replace(
        mod.SIM100K,
        genome_length=5_000_000,
        n_plasmids=1,
        plasmid_length=(80_000, 80_000),
        n_rrna_operons=7,
        n_trna=86,
        prophage_bp=60_000,
        bgc_bp=40_000,
        n_pseudogenes=40,
        n_sorfs=30,
        write_precomputed=False,
        proteomics=ProteomicsConfig(),
    )
    assert expected == mod.SIM5M, {
        f.name: (getattr(mod.SIM5M, f.name), getattr(expected, f.name))
        for f in dataclasses.fields(expected)
        if getattr(mod.SIM5M, f.name) != getattr(expected, f.name)
    }
    assert mod.SIM5M.write_precomputed is False, "the real tools must run"
    assert mod.SIM5M.genome_length == 5_000_000


def _records(fasta: Path) -> list[tuple[str, int]]:
    out: list[tuple[str, int]] = []
    for line in fasta.read_text().splitlines():
        if line.startswith(">"):
            out.append((line, 0))
        elif out:
            out[-1] = (out[-1][0], out[-1][1] + len(line.strip()))
    return out


@pytest.mark.slow
def test_sim5m_is_deterministic(tmp_path: Path) -> None:
    mod = _sim()
    a = mod.simulate(mod.SIM5M, mod.SIM5M_SEED, tmp_path / "a")
    b = mod.simulate(mod.SIM5M, mod.SIM5M_SEED, tmp_path / "b")
    for name in ("genome.fasta", "peptides.tsv", "bakta_like.json"):
        pa, pb = tmp_path / "a" / name, tmp_path / "b" / name
        assert pa.is_file(), name
        assert pa.read_bytes() == pb.read_bytes(), f"{name} differs between two calls"
    assert a.peptides is not None and b.peptides is not None
    recs = _records(tmp_path / "a" / "genome.fasta")
    assert len(recs) == 2, "one chromosome and one plasmid"
    assert recs[0][1] == 5_000_000, recs[0]
    assert recs[1][1] == 80_000, recs[1]
    assert all("topology=circular" in h for h, _ in recs), "both replicons circular"
    assert not (tmp_path / "a" / "config.yaml").exists(), "no precomputed layers"
