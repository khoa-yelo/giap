"""gippasim chains and a brute-force oracle for the M7-10 class error rates
(docs/specs/M7-10-class-error-rates.md; used by tests/unit/steps/test_variants_step_class_error.py
and tests/calibration/test_variant_class_error.py).

``build_chain`` simulates a genome with M3-19 proteomics and M3-20 de novo reads (optionally an
M3-18 strain: the sample's peptides on the reference genome, as tests/unit/steps/
test_variants_step.py does), runs preflight **through the library** with the de novo reads as
``--peptides`` and gippasim's ``reference_peptides.tsv`` as ``reference_peptides`` (lead
decision 1: preflight writes ``peptide_spectra.parquet``), then the mapping chain through
``helpers.mapping_steps.drive`` with a test ``$GIPPA_DB_DIR``.

``oracle_sites`` rebuilds the M7-04 ``run_variants`` chain from merged public pieces (M7-01/02/05)
to get each site's supporting peptides and class; ``oracle_rates`` counts rule 1-3 by brute force
from the *raw* gippasim tables (never from GIPPA's C17), so it is independent of the code under
test.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from gippa.config import resolve_config
from helpers.mapping_steps import drive, fresh_steps, read_tsv_df, write_test_db

Class = tuple[str, bool]
STEP_KEYS = ("mapping", "genome", "genes", "replicons", "peptides", "hits_tolerant")
REF_MAX_Q = 0.01


@dataclasses.dataclass
class ClassErrorChain:
    res: Any  # gippasim SimResult (genome/bakta replaced by the reference strain's when strain)
    cfg: Any
    avail: dict[str, Path]
    preflight_work: Path
    preflight_report: Any
    denovo_tsv: Path
    reference_tsv: Path

    @property
    def peptide_spectra(self) -> Path:
        return self.preflight_work / "peptide_spectra.parquet"

    def step_inputs(self, *, with_spectra: bool = True) -> dict[str, Path]:
        out = {k: self.avail[k] for k in STEP_KEYS}
        out["reference_peptides"] = self.reference_tsv
        if with_spectra:
            out["peptide_spectra"] = self.peptide_spectra
        return out


def library_preflight(
    genome: Path, peptides: Path, annotation: Path, reference: Path | None, out: Path,
    config_yaml: Path | None,
) -> tuple[Any, Path]:  # fmt: skip
    """``run_preflight`` + ``write_outputs`` with ``reference_peptides`` (mode ``preflight``)."""
    from gippa.exitcodes import ExitCode
    from gippa.preflight.run import PreflightInputs, exit_code, run_preflight, write_outputs

    cfg = resolve_config(config_yaml) if config_yaml is not None else resolve_config()
    inp = PreflightInputs(genome, None, peptides, None, annotation, reference_peptides=reference)
    db = write_test_db(out.with_name(out.name + "_db"))
    with pytest.MonkeyPatch.context() as mp:
        mp.setenv("GIPPA_DB_DIR", str(db))
        report, norm = run_preflight(inp, cfg, None, "preflight", out)
    errs = [(c.rule_id, c.message) for c in report.checks if c.level == "error"]
    assert exit_code(report, False) == ExitCode.SUCCESS, errs
    write_outputs(out, report, norm)
    return report, out / "work" / "preflight"


def build_chain(
    root: Path,
    seed: int,
    *,
    length: int = 50_000,
    pi_wrong: float = 0.05,
    strain: Any = None,
    denovo: Any = None,
    reference: Path | None = None,
) -> ClassErrorChain:
    """The chain above; ``reference`` replaces gippasim's ``reference_peptides.tsv``."""
    from gippasim.denovo import DenovoConfig
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SimConfig, simulate

    cfg_sim = SimConfig(genome_length=length, proteomics=ProteomicsConfig(pi_wrong=pi_wrong),
                        strain=strain, denovo=denovo or DenovoConfig(),
                        write_precomputed=True)  # fmt: skip
    res = simulate(cfg_sim, seed, root / "sim")
    if strain is not None:
        assert res.reference_fasta is not None and res.reference_json is not None
        res = dataclasses.replace(res, genome_fasta=res.reference_fasta,
                                  bakta_json=res.reference_json)  # fmt: skip
    dn = res.outdir / "peptides_denovo.tsv"
    ref = reference or res.outdir / "reference_peptides.tsv"
    report, work = library_preflight(res.genome_fasta, dn, res.bakta_json, ref, root / "pf",
                                     res.config_yaml)  # fmt: skip
    cfg = resolve_config(res.config_yaml, [f"seed={seed}", "sample_name=S1"])
    roots = {"genome": work / "genome.norm.fa", "peptides": work / "peptides.norm.parquet",
             "bakta_json": res.bakta_json, "reference_peptides": ref}  # fmt: skip
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, root / "work", roots)
    return ClassErrorChain(res, cfg, avail, work, report, dn, ref)


# ---------------------------------------------------------------------------- oracle


def oracle_sites(avail: Mapping[str, Path], cfg: Any) -> tuple[pd.DataFrame, list[frozenset]]:
    """M7-04 ``run_variants``'s ``sites`` (with M7-05 flags) and each row's ALT peptides."""
    from gippa.io.fasta import read_genome_fasta
    from gippa.io.work_tables import load_tables
    from gippa.mapping.initiator import annotated_starts, start_keys
    from gippa.steps.mapping import read_gene_model
    from gippa.steps.variants import support_of
    from gippa.variants.align import align_hits
    from gippa.variants.flags import flags_for, left_block_site, start_codons
    from gippa.variants.pileup import eligible, pileup, sites

    mapping = Path(avail["mapping"])
    frames, _ = load_tables({"peptide_hits": mapping / "peptide_hits.tsv"}, ["peptide_hits"])
    ph = frames["peptide_hits"]
    ph["diff_types"] = [";".join(v) if isinstance(v, list) else v for v in ph["diff_types"]]
    ht = pd.read_parquet(avail["hits_tolerant"])
    genome = {r.id: r.sequence for r in read_genome_fasta(Path(avail["genome"]))}
    gm = read_gene_model(Path(avail["genes"]), Path(avail["replicons"]))
    peptides = pd.read_parquet(avail["peptides"])
    support = support_of(ph, peptides, {g.locus_tag for g in gm.genes})
    ann = cfg.annotation  # M11-20: initiator start-site keys, as run_variants_qc
    init_sites, _ = annotated_starts(gm, genome, ann.start_codons, table=ann.translation_table)
    inits = start_keys(init_sites)
    aligns = align_hits(ph, ht, genome, gm.replicons, starts=inits)
    pile = pileup(eligible(aligns, support, cfg), gm, genome, cfg, support=support)
    st = sites(pile, cfg)
    peps = list(st.attrs["alt_peptides"])
    starts = start_codons(gm)
    tol = cfg.variants.mass_equivalent_tol_da
    st["flags"] = [  # ADR-0049: the left 2->1 block rule and swaps are flagged mass_equivalent
        tuple(flags_for(r, starts, cfg, {(support[p].source, support[p].search_db) for p in ps},
                        left_block=left_block_site(
                            r, [a for a in aligns if a.peptide in ps], genome, gm.replicons, tol)))
        for r, ps in zip(st.to_dict("records"), peps, strict=True)
    ]  # fmt: skip
    return st, peps


def site_class(row: Mapping[str, Any]) -> Class:
    flags = row.get("flags") or ()
    return (str(row["alt_type"]), "mass_equivalent" in tuple(flags))


def peptide_classes(st: pd.DataFrame, peps: Iterable[frozenset]) -> dict[str, set[Class]]:
    out: dict[str, set[Class]] = {}
    for row, ps in zip(st.to_dict("records"), peps, strict=True):
        for p in ps:
            out.setdefault(str(p), set()).add(site_class(row))
    return out


def _fold(p: str) -> str:
    return p.replace("I", "L")


def oracle_rates(
    classes: Mapping[str, set[Class]], denovo_tsv: Path, reference_tsv: Path,
    q_max: float = REF_MAX_Q,
) -> dict[Class, tuple[int, int]]:  # fmt: skip
    """Rules 1-3 by brute force on the raw tables: ``{class: (n, errors)}``."""
    dn = read_tsv_df(denovo_tsv)
    ref = read_tsv_df(reference_tsv)
    best: dict[tuple[str, str], str] = {}
    if "spectrum_id" in ref.columns:
        for r in ref.itertuples(index=False):
            if pd.isna(r.spectrum_id) or str(r.is_decoy).lower() == "true":
                continue
            if pd.isna(r.q_value) or float(r.q_value) > q_max:
                continue
            best[(str(r.run), str(r.spectrum_id))] = str(r.peptide)
    out: dict[Class, list[int]] = {}
    for r in dn.itertuples(index=False):
        if r.source != "denovo" or pd.isna(r.spectrum_id) or str(r.peptide) not in classes:
            continue
        hit = best.get((str(r.run), str(r.spectrum_id)))
        if hit is None:
            continue
        err = _fold(str(r.peptide)) != _fold(hit)
        for c in classes[str(r.peptide)]:
            v = out.setdefault(c, [0, 0])
            v[0] += 1
            v[1] += int(err)
    return {k: (v[0], v[1]) for k, v in out.items()}


def truth_rates(
    classes: Mapping[str, set[Class]], denovo_truth: Path, *, paired: set | None = None
) -> dict[Class, tuple[int, int]]:
    """True per-class error counts from gippasim's ``truth/denovo.tsv``: every ALT-supporting de
    novo read, an error when it differs from its spectrum's true peptide after I/L folding.
    ``paired`` restricts to those ``(run, spectrum_id)``."""
    tr = read_tsv_df(denovo_truth)
    out: dict[Class, list[int]] = {}
    for r in tr.itertuples(index=False):
        if str(r.peptide) not in classes:
            continue
        if paired is not None and (str(r.run), str(r.spectrum_id)) not in paired:
            continue
        err = _fold(str(r.peptide)) != _fold(str(r.true_peptide))
        for c in classes[str(r.peptide)]:
            v = out.setdefault(c, [0, 0])
            v[0] += 1
            v[1] += int(err)
    return {k: (v[0], v[1]) for k, v in out.items()}


def confident_pairs(reference_tsv: Path, q_max: float = REF_MAX_Q) -> set[tuple[str, str]]:
    ref = read_tsv_df(reference_tsv)
    keep = (ref["is_decoy"].str.lower() != "true") & (ref["q_value"].astype(float) <= q_max)
    return {(str(a), str(b)) for a, b in zip(ref.loc[keep, "run"], ref.loc[keep, "spectrum_id"],
                                             strict=True)}  # fmt: skip


def run_step(inputs: Mapping[str, Path], out: Path, cfg: Any) -> tuple[dict[str, Path], list]:
    """The registered ``variants`` step through ``run_step``; returns (outputs, log records)."""
    import logging

    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    paths = {"pileup": out / "tables/pileup.tsv.gz", "variants": out / "tables/variants.tsv",
             "variants_vcf": out / "tables/variants.vcf",
             "variants_qc": out / "work/variants_qc.parquet"}  # M7-14  # fmt: skip
    records: list[logging.LogRecord] = []

    class _H(logging.Handler):
        def emit(self, record: logging.LogRecord) -> None:
            records.append(record)

    h = _H(level=logging.DEBUG)
    log = logging.getLogger("gippa.steps.variants")
    old = log.level
    log.addHandler(h)
    log.setLevel(logging.DEBUG)
    try:
        reg.run_step("variants", cfg, dict(inputs), paths)
    finally:
        log.removeHandler(h)
        log.setLevel(old)
    return paths, records


def read_variants(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False)


def table_classes(variants: pd.DataFrame) -> list[Class]:
    return [(str(t), "mass_equivalent" in str(f).split(";"))
            for t, f in zip(variants["alt_type"], variants["flags"], strict=True)]  # fmt: skip


def measured_rates(variants: pd.DataFrame) -> dict[Class, tuple[float, int]]:
    """``{class: (rate, n)}`` of the estimated rows of a written ``variants.tsv``; asserts every
    row of one class carries the same value (rule 4: pooled per class)."""
    out: dict[Class, tuple[float, int]] = {}
    for c, r, n in zip(table_classes(variants), variants["class_error_rate"],
                       variants["class_error_n"], strict=True):  # fmt: skip
        if r == "unestimated":
            assert n in ("NA", ""), (c, r, n)
            continue
        v = (float(r), int(n))
        assert out.setdefault(c, v) == v, f"class {c}: {out[c]} and {v} on two rows"
    return out
