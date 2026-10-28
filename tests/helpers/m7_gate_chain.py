"""One M7-13a gate replicate: a gippasim scenario through preflight, the mapping chain and the
registered ``variants`` step (docs/specs/M7-13-calibration-gate.md rule 1; test plumbing, no
statistics). Used by tests/calibration/test_variants_gate.py and the attack file.

``run_chain`` simulates ``sim_cfg`` at ``seed``; with a strain the sample's peptides are mapped to
the reference strain (``reference.fasta``, ``reference_bakta_like.json``), as M7-10's
``helpers.class_error_chain.build_chain`` does. ``peptides`` / ``reference`` name gippasim output
files (``None``: no reference peptides); ``rewrite`` may replace either input file before preflight
(the A6c systematic swaps, the reference-proteome arm), keyed by the file name. A test
``$GIPPA_DB_DIR`` is used (never the real one).
"""

from __future__ import annotations

import dataclasses
import shutil
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from gippa.config import resolve_config
from helpers.class_error_chain import STEP_KEYS, library_preflight, oracle_sites, run_step
from helpers.mapping_steps import drive, fresh_steps

Rewrite = Callable[[Path, Any], Path]  # (gippasim file, SimResult) -> file to use instead


@dataclasses.dataclass
class GateRun:
    res: Any  # gippasim SimResult (genome/bakta replaced by the reference strain's with a strain)
    cfg: Any
    avail: dict[str, Path]
    variants: pd.DataFrame  # variants.tsv as written (strings; codon_start 1-based)
    variants_path: Path
    peptides_tsv: Path
    reference_tsv: Path | None

    @property
    def outdir(self) -> Path:
        return Path(self.res.outdir)

    def site_peptides(self) -> dict[tuple[str, int, str, str, str], frozenset[str]]:
        """``(contig, codon_start 1-based, strand, alt_type, alt_aa) -> ALT peptides`` of every
        ``sites`` row (the M7-04 chain rebuilt from merged pieces)."""
        st, peps = oracle_sites(self.avail, self.cfg)
        return {
            (str(r["contig"]), int(r["codon_start"]) + 1, str(r["strand"]), str(r["alt_type"]),
             str(r["alt_aa"])): frozenset(map(str, p))
            for r, p in zip(st.to_dict("records"), peps, strict=True)
        }  # fmt: skip


def run_chain(
    sim_cfg: Any,
    seed: int,
    root: Path,
    *,
    peptides: str = "peptides_denovo.tsv",
    reference: str | None = "reference_peptides.tsv",
    overrides: Sequence[str] = (),
    rewrite: Mapping[str, Rewrite] | None = None,
) -> GateRun:
    from gippasim.simulate import simulate

    res = simulate(sim_cfg, seed, root / "sim")
    if res.reference_fasta is not None and res.reference_json is not None:
        res = dataclasses.replace(res, genome_fasta=res.reference_fasta,
                                  bakta_json=res.reference_json)  # fmt: skip
    rewrite = dict(rewrite or {})
    pep = Path(res.outdir) / peptides
    if peptides in rewrite:
        pep = rewrite[peptides](pep, res)
    ref = None if reference is None else Path(res.outdir) / reference
    if reference is not None and reference in rewrite:
        assert ref is not None
        ref = rewrite[reference](ref, res)
    _, work = library_preflight(res.genome_fasta, pep, res.bakta_json, ref, root / "pf",
                                res.config_yaml)  # fmt: skip
    cfg = resolve_config(res.config_yaml, [f"seed={seed}", "sample_name=S1", *overrides])
    roots: dict[str, Path] = {"genome": work / "genome.norm.fa",
                              "peptides": work / "peptides.norm.parquet",
                              "bakta_json": res.bakta_json}  # fmt: skip
    if ref is not None:
        roots["reference_peptides"] = ref
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, root / "work", roots)
    inputs = {k: avail[k] for k in STEP_KEYS}
    if ref is not None:
        inputs["reference_peptides"] = ref
    spectra = work / "peptide_spectra.parquet"
    if spectra.exists():
        inputs["peptide_spectra"] = spectra
    paths, _ = run_step(inputs, root / "out", cfg)
    vt = pd.read_csv(paths["variants"], sep="\t", dtype=str, keep_default_na=False)
    return GateRun(res, cfg, avail, vt, paths["variants"], pep, ref)


def cleanup(root: Path) -> None:
    shutil.rmtree(root, ignore_errors=True)
