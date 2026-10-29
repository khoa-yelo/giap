"""Genetic codes 4 and 25: the M11-16e review leftovers N4 and N5 (M11-16h;
docs/specs/M11-16-translation-tables-4-25.md 6.1, example 10, AC 5, AC 11, Q-L3).

Readings pinned by the test author where the spec is silent:

* N3 (codon-usage internal-stop eligibility) is not here: the rule would change table-11
  output (AC 1), e.g. selenoprotein CDSs with an internal TGA; parked for Khoa (lead,
  2026-10-28). The test author's N3 tests are in ce93f7f1.
* N4 (Q-L3, "Proposed: add it only when != 11"): ``write_vcf(..., translation_table=<t>)``
  (keyword, default 11) writes ``##gippa_translation_table=<t>`` at 4 and 25, right after
  ``##reference=`` and before the first ``##contig``; at 11 there is no such line and the file is
  byte-identical to today's. The variants step passes ``cfg.annotation.translation_table``.
* N5 (AC 5, one table source): ``novel.candidates._stops`` and
  ``genetable.covariates._translation_maps`` agree with ``gippa.model.genetic_code`` at 4, 11 and
  25 and, once they read it, raise ``ValueError`` for the tables it refuses (1, 2, 12, 26).
  No ``src/gippa`` module imports ``Bio.Data.CodonTable`` (no such guard existed in scripts/ci or
  tests/unit).

No statistical claim is made here, so there is no calibration test (spec 8.4).
"""

from __future__ import annotations

import ast
from pathlib import Path
from typing import Any

import pytest

from gippa.config import GippaConfig
from gippa.io.vcf import VcfHeaderField, VcfRecord, write_vcf
from gippa.model import Replicon
from gippa.model.genetic_code import codon_table, stop_codons

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src" / "gippa"
TABLES = (4, 11, 25)
UNSUPPORTED = (1, 2, 12, 26)


def cfg_at(table: int, **kw: Any) -> GippaConfig:
    return GippaConfig.model_validate({"seed": 1, "annotation": {"translation_table": table}, **kw})


# ======================================================================== N4: VCF header

C = Replicon("chr", 2000, "circular")
SHA = "0123456789abcdef" * 4
INFO = [VcfHeaderField("AA_REF", "1", "String", "Reference residue")]
FORMAT = [
    VcfHeaderField("GT", "1", "String", "Genotype"),
    VcfHeaderField("AD", "R", "Integer", "Distinct peptides for REF and each ALT"),
    VcfHeaderField("DP", "1", "Integer", "Total distinct peptides"),
]
REC = VcfRecord(contig="chr", pos0=999, ref="TGA", alts=("TGG",), info={"AA_REF": "W"},
                gt="1", ad=(0, 2), dp=2)  # fmt: skip


def _vcf(path: Path, **kw: Any) -> list[str]:
    write_vcf(path, [REC], [C], gippa_version="0.1.0", genome_sha256=SHA, sample="S1",
              info_fields=INFO, format_fields=FORMAT, **kw)  # fmt: skip
    return path.read_text().splitlines()


@pytest.mark.parametrize("table", (4, 25))
def test_vcf_header_records_a_non_default_table_after_reference(table: int, tmp_path: Path) -> None:
    """Q-L3: ``##gippa_translation_table=<t>`` right after ``##reference=``, once."""
    lines = _vcf(tmp_path / "v.vcf", translation_table=table)
    meta = [ln for ln in lines if ln.startswith("##")]
    assert meta[:4] == [
        "##fileformat=VCFv4.3",
        "##source=gippa-0.1.0",
        f"##reference=sha256:{SHA}",
        f"##gippa_translation_table={table}",
    ]
    assert meta[4] == "##contig=<ID=chr,length=2000>"
    assert sum(ln.startswith("##gippa_translation_table") for ln in lines) == 1


def test_vcf_header_at_11_has_no_table_line_and_is_unchanged(tmp_path: Path) -> None:
    """Q-L3 "only when != 11": the code-11 VCF keeps today's bytes (explicit 11 = default)."""
    explicit = _vcf(tmp_path / "a.vcf", translation_table=11)
    default = _vcf(tmp_path / "b.vcf")
    assert explicit == default
    assert not any("translation_table" in ln for ln in explicit)
    assert explicit[:4] == [
        "##fileformat=VCFv4.3",
        "##source=gippa-0.1.0",
        f"##reference=sha256:{SHA}",
        "##contig=<ID=chr,length=2000>",
    ]


@pytest.mark.parametrize("table", (4, 25))
def test_vcf_with_table_line_passes_the_output_check(table: int, tmp_path: Path) -> None:
    from gippa.io.output_checks import check_vcf

    p = tmp_path / "v.vcf"
    _vcf(p, translation_table=table)
    check_vcf("variants_vcf", p)


@pytest.fixture(scope="module")
def variants_chain(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """A 30-kb gippasim genome mapped at 11 (``search_db: same``); the variants step reruns on
    it at each table, which only changes the VCF header and the TGA reading."""
    from helpers.mapping_steps import drive, fresh_steps, preflight, roots_of

    from gippa.config import resolve_config
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.simulate import SimConfig, simulate

    root = Path(tmp_path_factory.mktemp("v16h"))
    res = simulate(SimConfig(genome_length=30_000, proteomics=ProteomicsConfig(),
                             write_precomputed=False), 11, root / "sim")  # fmt: skip
    run = preflight(res, root / "pf")
    cfg = resolve_config(run.config_yaml, ["seed=11"])
    with pytest.MonkeyPatch.context() as mp:
        reg = fresh_steps(mp)
        avail = drive(reg, cfg, root / "work", roots_of(run))
    keys = ("mapping", "genome", "genes", "replicons", "peptides", "hits_tolerant")
    return {k: avail[k] for k in keys}


@pytest.mark.parametrize("table", TABLES)
def test_variants_step_writes_the_run_table_into_the_vcf(
    table: int, variants_chain: dict[str, Path], tmp_path: Path
) -> None:
    from gippa.steps import load_step_modules
    from gippa.steps import registry as reg

    load_step_modules()
    outs = {
        "pileup": tmp_path / "pileup.tsv.gz",
        "variants": tmp_path / "variants.tsv",
        "variants_vcf": tmp_path / "variants.vcf",
        "variants_qc": tmp_path / "variants_qc.parquet",
    }
    reg.run_step("variants", cfg_at(table, sample_name="S1"), dict(variants_chain), outs)
    meta = [ln for ln in outs["variants_vcf"].read_text().splitlines() if ln.startswith("##")]
    line = f"##gippa_translation_table={table}"
    if table == 11:
        assert not any(ln.startswith("##gippa_translation_table") for ln in meta)
    else:
        assert meta[3] == line and meta[2].startswith("##reference=sha256:")


# ======================================================================== N5: one table source


def _model_maps(table: int) -> tuple[dict[str, str], frozenset[str]]:
    ct = codon_table(table)
    return {c: a for c, a in ct.items() if a != "*"}, stop_codons(table)


@pytest.mark.parametrize("table", TABLES)
def test_novel_candidate_stops_equal_the_model(table: int) -> None:
    from gippa.novel.candidates import _stops

    assert _stops(table) == stop_codons(table)


@pytest.mark.parametrize("table", TABLES)
def test_covariate_translation_maps_equal_the_model(table: int) -> None:
    from gippa.genetable.covariates import _translation_maps

    forward, stops = _translation_maps(table)
    want_forward, want_stops = _model_maps(table)
    assert dict(forward) == want_forward
    assert frozenset(stops) == want_stops


@pytest.mark.parametrize("table", UNSUPPORTED)
def test_novel_candidate_stops_refuse_tables_the_model_refuses(table: int) -> None:
    from gippa.novel.candidates import _stops

    with pytest.raises(ValueError, match="translation table"):
        _stops(table)


@pytest.mark.parametrize("table", UNSUPPORTED)
def test_covariate_translation_maps_refuse_tables_the_model_refuses(table: int) -> None:
    from gippa.genetable.covariates import _translation_maps

    with pytest.raises(ValueError, match="translation table"):
        _translation_maps(table)


@pytest.mark.parametrize("table", TABLES)
def test_protein_sequence_reads_tga_with_the_model_table(table: int) -> None:
    """Pin (passes today): covariates translate TGA as the model does (W at 4, G at 25)."""
    from gippa.genetable.covariates import protein_sequence

    got = protein_sequence(None, "ATGAAATGAGCCTAA", "cds", table)
    assert got == {11: "MK", 4: "MKWA", 25: "MKGA"}[table]


_ALLOWED = {"gippa/model/genetic_code.py"}


def _codontable_imports(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    hits: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            hits += [a.name for a in node.names if a.name.startswith("Bio.Data.CodonTable")]
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            if node.module.startswith("Bio.Data.CodonTable"):
                hits.append(node.module)
            elif node.module in ("Bio.Data", "Bio") and any(
                a.name in ("CodonTable", "Data") for a in node.names
            ):
                hits.append(f"{node.module}.{','.join(a.name for a in node.names)}")
    return hits


def test_guard_finds_each_import_form(tmp_path: Path) -> None:
    """The guard's own oracle: every way to reach ``Bio.Data.CodonTable`` is caught."""
    forms = [
        "from Bio.Data import CodonTable",
        "import Bio.Data.CodonTable",
        "import Bio.Data.CodonTable as ct",
        "from Bio.Data.CodonTable import unambiguous_dna_by_id",
        "from Bio import Data",
        "def f():\n    from Bio.Data import CodonTable\n",
    ]
    for i, src in enumerate(forms):
        p = tmp_path / f"m{i}.py"
        p.write_text(src + "\n")
        assert _codontable_imports(p), src
    ok = tmp_path / "ok.py"
    ok.write_text("from Bio import SeqIO\nfrom Bio.SeqUtils.ProtParam import ProteinAnalysis\n")
    assert _codontable_imports(ok) == []


def test_no_module_but_genetic_code_imports_biopython_codon_tables() -> None:
    """N5: ``gippa.model.genetic_code`` is the one source of the genetic code in ``src/gippa``."""
    bad = {
        rel: hits
        for path in sorted(SRC.rglob("*.py"))
        if (rel := str(path.relative_to(SRC.parent))) not in _ALLOWED
        and (hits := _codontable_imports(path))
    }
    assert bad == {}
