"""M7-15 integration hardening on the 50-kb M7 fixture (docs/specs/M7-15-integration-hardening.md
items 2, 5, 6 and 7; examples 4; GIAP_SPEC.md §6.3, §9.5, §11.2, §11.6 TE7).

Inputs: tests/helpers/m7_chain.py (the sim50k_m7 golden fixture, any seed), a test
``$GIAP_DB_DIR``, ``giap run --executor local`` unless stated.

* **Item 2 / example 4 -- in-pipeline contract checks.** The ``variants`` and ``novel_orfs`` steps
  validate what they wrote before returning: the TSVs against their Frictionless resources, the VCF
  and the GFF3s with the format readers. The test hook wraps the step module's writer
  (``m7_chain.after``) so the file is damaged after a correct write; the step must then fail with
  a message naming the output (in-process ``run_step``), and ``giap run`` must exit 4 naming the
  step itself (``step 'variants' failed``, not a downstream reader such as ``plots``, which is
  what happens today) and ``logs/<step>.log``. Damage: an integer column set to ``x`` in
  ``pileup.tsv.gz`` (gzip kept), ``variants.tsv`` and ``novel_orfs.tsv``; a record with a
  non-numeric POS and too few columns appended to ``variants.vcf``; a CDS with start > end
  appended to ``novel_orfs.gff3`` and ``annotation_patch.gff3``.
* **Item 5 -- VCF FILTER** (M7-03 decision 6): every record ``PASS``; no generic
  ``##FILTER=<...Description="GIAP filter ...">`` line; any ``##FILTER`` line declares ``PASS``.
* **Item 6 -- ``novel_candidate_id``** (M7-06 lead decision 1, question 6: the proposed default,
  back-filled in M7-15). The published ``mapping/peptide_hits.tsv`` carries, for every peptide an
  ORF of ``novel_orfs.tsv`` lists, that ORF's ``orf_id`` on a row inside the ORF's span (I = L
  folded); every non-NA value is an ``orf_id``; ``cds_inframe`` / ``rna_overlap`` rows stay NA;
  the file still validates.
* **Item 7 -- TE7**: two local runs with the same seed give the seven files byte-identical; the
  Nextflow executor (``--profile farmshare_local``, every layer precomputed, no tool runs) gives
  the same seven files and the same ``mapping/peptide_hits.tsv``.
"""

from __future__ import annotations

import re
import shutil
import sys
from collections.abc import Callable
from pathlib import Path
from typing import Any

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "helpers"))
import m7_chain as m7

SEED = 11


@pytest.fixture(scope="module")
def sim(tmp_path_factory: pytest.TempPathFactory) -> Any:
    return m7.simulate(tmp_path_factory.mktemp("m7_15_sim") / "sim", SEED)


@pytest.fixture(scope="module")
def root(tmp_path_factory: pytest.TempPathFactory) -> Path:
    return tmp_path_factory.mktemp("m7_15")


@pytest.fixture(scope="module")
def local(sim: Any, root: Path) -> Path:
    with m7.db_env(root):
        return m7.run_local(sim, root / "A" / "OUT", SEED)


# ----------------------------------------------------------------------------- item 2 (in-process)


def _run_step(name: str, out: Path, dest: Path) -> dict[str, Path]:
    from helpers import novel_steps_fixture as nf

    cfg = m7.resolved_config(out)
    inputs = nf.step_inputs(name, out / "work")
    dest.mkdir(parents=True)
    return nf.run(name, cfg, inputs, dest)


def _tsv_damage(column: str) -> Callable[[Path], None]:
    return lambda p: m7.corrupt_tsv(p, column)


DAMAGE: dict[str, tuple[str, str, str, Callable[[Path], None]]] = {
    # case: (step, module attribute, output key, damage)
    "pileup": ("variants", "write_pileup", "pileup", _tsv_damage("codon_start")),
    "variants": ("variants", "write_variants", "variants", _tsv_damage("codon_start")),
    "variants_vcf": ("variants", "write_vcf", "variants_vcf",
                     lambda p: m7.append_line(p, m7.BAD_VCF_RECORD)),
    "novel_orfs": ("novel_orfs", "write_novel_orfs", "novel_orfs", _tsv_damage("end")),
    "novel_orfs_gff3": ("novel_orfs", "write_novel_gff3", "novel_orfs_gff3",
                        lambda p: m7.append_line(p, m7.BAD_GFF3_LINE)),
    "annotation_patch": ("novel_orfs", "write_annotation_patch", "annotation_patch",
                         lambda p: m7.append_line(p, m7.BAD_GFF3_LINE)),
}  # fmt: skip
_MODULES = {"variants": "giap.steps.variants", "novel_orfs": "giap.steps.novel"}


def test_clean_steps_pass_their_own_checks(local: Path, tmp_path: Path) -> None:
    """Without damage both steps re-run in-process on the run's work files and succeed."""
    for step in ("variants", "novel_orfs"):
        outs = _run_step(step, local, tmp_path / step)
        assert all(p.exists() for p in outs.values()), outs


@pytest.mark.parametrize("case", sorted(DAMAGE))
def test_example_4_damaged_output_fails_its_step(
    case: str, local: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import importlib

    step, attr, key, damage = DAMAGE[case]
    mod = importlib.import_module(_MODULES[step])
    monkeypatch.setattr(mod, attr, m7.after(getattr(mod, attr), damage))
    with pytest.raises(Exception) as exc:  # the step's own check: any error type
        _run_step(step, local, tmp_path / step)
    assert key in str(exc.value), f"the failure does not name {key}: {exc.value!r}"


def test_example_4_damaged_variants_tsv_makes_giap_run_exit_4(
    sim: Any, root: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import giap.steps.variants as sv
    from giap.exitcodes import ExitCode

    monkeypatch.setattr(sv, "write_variants",
                        m7.after(sv.write_variants, _tsv_damage("codon_start")))  # fmt: skip
    out = tmp_path / "OUT"
    with m7.db_env(tmp_path):
        code, log = m7.invoke(m7.argv(sim, out, SEED, "--executor", "local"))
    assert code == ExitCode.STEP_FAILED, log[-3000:]
    assert re.search(r"step '?variants'? failed", log), log[-3000:]
    assert "logs/variants.log" in log
    assert (out / "logs" / "variants.log").is_file()


def test_example_4_damaged_annotation_patch_makes_giap_run_exit_4(
    sim: Any, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    import giap.steps.novel as sn
    from giap.exitcodes import ExitCode

    monkeypatch.setattr(sn, "write_annotation_patch",
                        m7.after(sn.write_annotation_patch,
                                 lambda p: m7.append_line(p, m7.BAD_GFF3_LINE)))  # fmt: skip
    out = tmp_path / "OUT"
    with m7.db_env(tmp_path):
        code, log = m7.invoke(m7.argv(sim, out, SEED, "--executor", "local"))
    assert code == ExitCode.STEP_FAILED, log[-3000:]
    assert re.search(r"step '?novel_orfs'? failed", log), log[-3000:]
    assert "logs/novel_orfs.log" in log
    assert (out / "logs" / "novel_orfs.log").is_file()


# ----------------------------------------------------------------------------- §9.5 on the run


def test_published_tsvs_validate(local: Path) -> None:
    m7.validate_tsv(local / "tables" / "pileup.tsv.gz", "pileup")
    m7.validate_tsv(local / "tables" / "variants.tsv", "variants")
    m7.validate_tsv(local / "tables" / "novel_orfs.tsv", "novel_orfs")


def test_published_vcf_and_gff3_read_back(local: Path) -> None:
    from giap.io.gff3 import read_gff3
    from giap.steps.common import read_replicons

    reps = read_replicons(local / "work" / "replicons")
    for name in ("novel_orfs.gff3", "annotation_patch.gff3"):
        assert read_gff3(local / "tables" / name, reps), name
    pysam = pytest.importorskip("pysam")
    with pysam.VariantFile(str(local / "tables" / "variants.vcf")) as vf:
        assert list(vf.header.samples) == ["peptides_denovo"]
        assert list(vf), "no VCF record on the strain fixture"


# ----------------------------------------------------------------------------- item 5


def test_item_5_vcf_filter_is_pass_only(local: Path) -> None:
    lines = (local / "tables" / "variants.vcf").read_text(encoding="utf-8").splitlines()
    recs = [ln.split("\t") for ln in lines if ln and not ln.startswith("#")]
    assert recs
    assert {r[6] for r in recs} == {"PASS"}
    filters = [ln for ln in lines if ln.startswith("##FILTER=")]
    assert not any('Description="GIAP filter ' in ln for ln in filters), filters
    assert all(ln.startswith("##FILTER=<ID=PASS,") for ln in filters), filters


# ----------------------------------------------------------------------------- item 6


def _il(s: str) -> str:
    return s.replace("I", "L")


def test_item_6_novel_candidate_id_back_filled(local: Path) -> None:
    hits = m7._rows(local / "mapping" / "peptide_hits.tsv")
    orfs = m7._rows(local / "tables" / "novel_orfs.tsv")
    assert orfs, "no novel ORF on the fixture"
    ids = {o["orf_id"] for o in orfs}
    filled = [h for h in hits if h["novel_candidate_id"] not in ("NA", "")]
    assert filled, "novel_candidate_id is NA on every peptide_hits.tsv row"
    stray = sorted({h["novel_candidate_id"] for h in filled} - ids)
    assert not stray, f"novel_candidate_id values that are no orf_id: {stray[:5]}"
    annotated = [h for h in filled if h["locus_class"] in ("cds_inframe", "rna_overlap")]
    assert not annotated, annotated[:3]
    for o in orfs:
        ends = [int(o["end"])] + ([int(o["start_longest"])] if o["start_longest"] != "NA" else [])
        lo, hi = min(ends), max(ends)
        mine = [h for h in filled if h["novel_candidate_id"] == o["orf_id"]]
        for h in mine:
            assert (h["contig"], h["strand"]) == (o["contig"], o["strand"]), (o["orf_id"], h)
            assert lo <= int(h["start"]) and int(h["end"]) <= hi, (o["orf_id"], h)
        have = {_il(h["peptide"]) for h in mine}
        want = {_il(p) for p in o["peptides"].split(";") if p}
        assert want <= have, (o["orf_id"], sorted(want - have)[:3])


def test_item_6_peptide_hits_still_validates(local: Path) -> None:
    m7.validate_tsv(local / "mapping" / "peptide_hits.tsv", "peptide_hits")


# ----------------------------------------------------------------------------- item 7 (TE7)


@pytest.mark.slow
def test_item_7_two_local_runs_byte_identical(sim: Any, local: Path, root: Path) -> None:
    with m7.db_env(root / "b"):
        again = m7.run_local(sim, root / "B" / "OUT", SEED)
    names = (*m7.SEVEN, "mapping/peptide_hits.tsv")
    a, b = m7.files(local, names), m7.files(again, names)
    assert [k for k in a if a[k] != b[k]] == []


@pytest.mark.slow
@pytest.mark.needs_nextflow
def test_item_7_nextflow_gives_the_same_files(
    sim: Any, local: Path, root: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    if shutil.which("nextflow") is None:
        pytest.skip("nextflow not on PATH")
    monkeypatch.setenv("NXF_ANSI_LOG", "false")
    out = root / "NF" / "OUT"
    with m7.db_env(root / "nf"):
        code, log = m7.invoke(m7.argv(sim, out, SEED, "--executor", "nextflow",
                                      "--profile", "farmshare_local"))  # fmt: skip
    assert code == 0, log[-4000:]
    names = (*m7.SEVEN, "mapping/peptide_hits.tsv")
    a, b = m7.files(local, names), m7.files(out, names)
    assert [k for k in a if a[k] != b[k]] == []
