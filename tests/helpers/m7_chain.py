"""Driver and checks shared by the M7-15 tests (docs/specs/M7-15-integration-hardening.md).

**The 50-kb M7 fixture (lead decision 1, 2026-10-22).** A second golden directory,
``tests/data/golden/sim50k_m7/``; the sim50k goldens (``build_m4_fixture``) are not touched. The
fixture is a giapsim nucleotide genome of 50 kb with the M7-12 settings at that size:

* simulator-v2 signals at strength 0.8 and every precomputed layer (``write_precomputed``), so
  ``giap run`` runs no tool under either executor (the M7-08 integration pattern);
* a strain difference: ``StrainConfig(snv_rate_per_kb=1.0, indel_rate_per_kb=0.0)``. GIAP is given
  the related **reference** genome (``reference.fasta`` / ``reference_bakta_like.json``, M7-12 rule
  1) and the sample's peptides, with ``--genome-relation related``. Indels are left out (M7-12's
  ``variant_f1`` uses 0.2 / kb): the precomputed layers are written on the sample's coordinates,
  and SNVs alone keep them valid on the reference (sample and reference coordinates are equal);
* de novo peptides (``DenovoConfig()``, every error type) as ``--peptides`` and their paired
  ``reference_peptides.tsv`` as ``--reference-peptides`` (M7-10);
* planted novel ORFs of the three classes and N-terminal extensions; with the M7-12 knobs
  (``novel_orf_len_aa = (30, 80)``, ``plant_novel_rbs``) when giapsim has them. The golden test
  requires them (``has_m7_12``); the other M7-15 tests run with or without them.
* **class-FDR power (M7-11b).** The class q is computed per class, and its floor is c / T in the
  de novo regime (T = target peptides of the class), so a class with fewer than about c / alpha
  target peptides cannot be called (``novel_orfs.class_q_floor_above_alpha``). With 3 ORFs per
  class every class stayed under that (10-43 target peptides at GOLDEN_SEED) and no planted ORF
  was called. The fixture therefore plants 28 ``antisense`` ORFs (133 target peptides, floor
  ~ 0.0065 at GOLDEN_SEED); the other classes keep 3 and still emit the floor warning. alpha
  (``stats.alpha.novel_orfs``) and the call rule are unchanged.

The run config is the simulator's ``config.yaml`` with ``annotation.precomputed.bakta_json``
pointed at the reference JSON (``run_config``); locus tags are the same on both genomes.

**Golden comparison.** Byte for byte (gzip compared decompressed, as M5-11), except:

* the VCF ``##source=giap-<version>`` line is compared with the version masked (``mask_version``),
  so a version bump alone does not regenerate the golden;
* BLAS-dependent numeric columns follow the M5-11 decision 44 pattern (``math.isclose``, rel and
  abs 1e-3) through ``BLAS_NUMERIC``. It is empty: no column of the five files comes from a
  BLAS-backed fit (variants: counts and codes; novel ORFs: precomputed layers, quantiles of
  observed TIRs, counts; M7-10 / M7-11 rates are count-based). A column joins it only with a
  recorded cross-platform difference, as decision 44 did.

**In-pipeline contract hooks (M7-15 item 2, example 4).** The steps write through module-level
names, which the tests wrap to corrupt the written file after the real writer ran:
``giap.steps.variants.{write_pileup, write_variants, write_vcf}`` and
``giap.steps.novel.{write_novel_orfs, write_novel_gff3, write_annotation_patch}``. The step's
own check must then fail the step with a message naming the output.
"""

from __future__ import annotations

import gzip
import json
import math
import os
import sys
from collections.abc import Callable, Iterator, Mapping, Sequence
from contextlib import contextmanager
from dataclasses import fields
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import m5_chain as m5

REPO = Path(__file__).resolve().parents[2]
GOLDEN_DIR = REPO / "tests" / "data" / "golden" / "sim50k_m7"
SCHEMA_DIR = REPO / "src" / "giap" / "schemas" / "tables"
GENOME_LENGTH = 50_000
#: 28 antisense ORFs give one class enough target peptides to pass the M7-11b class q floor
#: (module doc, "class-FDR power"); the other classes keep 3
NOVEL_COUNTS = (("intergenic", 3), ("cds_outframe", 3), ("antisense", 28))
N_NTERM = 3
#: the M7-12 SimConfig knobs the golden needs (M7-12 "novel_orf_power" settings)
M7_12_KNOBS: Mapping[str, Any] = {"novel_orf_len_aa": (30, 80), "plant_novel_rbs": True}

#: the seven files Modules 5 and 6 publish (M7 index "End state"), relative to OUT
SEVEN: tuple[str, ...] = (
    "tables/pileup.tsv.gz",
    "tables/variants.tsv",
    "tables/variants.vcf",
    "tables/novel_orfs.tsv",
    "tables/novel_orfs.gff3",
    "tables/novel_orfs.faa",
    "tables/annotation_patch.gff3",
)
#: the golden files (spec item 1): logical name (relative to OUT) -> committed file name
GOLDEN_NAMES: Mapping[str, str] = {
    "tables/pileup.tsv.gz": "pileup.tsv.gz",
    "tables/variants.tsv": "variants.tsv",
    "tables/variants.vcf": "variants.vcf",
    "tables/novel_orfs.tsv": "novel_orfs.tsv",
    "tables/annotation_patch.gff3": "annotation_patch.gff3",
}
#: BLAS-dependent numeric columns per golden TSV (M5-11 decision 44 pattern); see module doc
BLAS_NUMERIC: Mapping[str, tuple[str, ...]] = {}
REL_TOL, ABS_TOL = 1e-3, 1e-3
MAX_GOLDEN_BYTES = m5.MAX_GOLDEN_BYTES


# --------------------------------------------------------------------------- fixture


def has_m7_12() -> bool:
    """giapsim has the M7-12 knobs (``novel_orf_len_aa``, ``plant_novel_rbs``)."""
    from giapsim.simulate import SimConfig

    names = {f.name for f in fields(SimConfig)}
    return all(k in names for k in M7_12_KNOBS)


def sim_config() -> Any:
    from giapsim.denovo import DenovoConfig
    from giapsim.proteomics import ProteomicsConfig
    from giapsim.simulate import SimConfig
    from giapsim.strain import StrainConfig

    kw: dict[str, Any] = {
        "genome_length": GENOME_LENGTH,
        "rbs_strength": 0.8,
        "promoter_strength": 0.8,
        "terminator_strength": 0.8,
        "write_precomputed": True,
        "strain": StrainConfig(snv_rate_per_kb=1.0, indel_rate_per_kb=0.0),
        "proteomics": ProteomicsConfig(),
        "denovo": DenovoConfig(),
        "n_novel_orfs": NOVEL_COUNTS,
        "n_nterm_extensions": N_NTERM,
    }
    if has_m7_12():
        kw.update(M7_12_KNOBS)
    return SimConfig(**kw)


def simulate(out: Path, seed: int) -> Any:
    import giapsim

    sim = giapsim.simulate(sim_config(), seed, out)
    assert sim.reference_fasta is not None and sim.reference_json is not None
    assert sim.config_yaml is not None
    return sim


def denovo(sim: Any) -> Path:
    p = Path(sim.outdir) / "peptides_denovo.tsv"
    assert p.is_file(), "giapsim writes peptides_denovo.tsv with SimConfig.denovo"
    return p


def reference_peptides(sim: Any) -> Path:
    p = Path(sim.outdir) / "reference_peptides.tsv"
    assert p.is_file(), "giapsim writes reference_peptides.tsv with SimConfig.denovo"
    return p


def run_config(sim: Any, dest: Path) -> Path:
    """The simulator's ``config.yaml`` with ``bakta_json`` pointed at the reference JSON."""
    src = Path(sim.config_yaml).read_text(encoding="utf-8")
    text = src.replace(str(Path(sim.bakta_json).resolve()), str(Path(sim.reference_json).resolve()))
    assert text != src, "config.yaml does not name the sample's Bakta JSON"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(text, encoding="utf-8")
    return dest


def argv(sim: Any, out: Path, seed: int, *extra: str, config: Path | None = None) -> list[str]:
    """``giap run`` on the reference genome with the sample's de novo peptides."""
    cfg = config if config is not None else run_config(sim, out.parent / "run_config.yaml")
    args = ["run", "--genome", str(sim.reference_fasta), "--peptides", str(denovo(sim)),
            "--reference-peptides", str(reference_peptides(sim)),
            "--annotation", str(sim.reference_json), "--config", str(cfg),
            "--genome-relation", "related", "--outdir", str(out), "--modules", "core",
            "--seed", str(seed)]  # fmt: skip
    return [*args, *extra]


def invoke(args: Sequence[str]) -> tuple[int, str]:
    """``giap run`` through the click runner: (exit code, stdout + stderr)."""
    from click.testing import CliRunner

    from giap.cli import cli

    res = CliRunner().invoke(cli, list(args))
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    err = res.stderr if res.stderr_bytes is not None else ""
    return res.exit_code, res.output + err


@contextmanager
def db_env(root: Path) -> Iterator[Path]:
    """A test ``$GIAP_DB_DIR`` (never the real one)."""
    from helpers.mapping_steps import write_test_db

    old = os.environ.get("GIAP_DB_DIR")
    db = write_test_db(root / "test_db")
    os.environ["GIAP_DB_DIR"] = str(db)
    try:
        yield db
    finally:
        if old is None:
            os.environ.pop("GIAP_DB_DIR", None)
        else:
            os.environ["GIAP_DB_DIR"] = old


def run_local(sim: Any, out: Path, seed: int, *extra: str) -> Path:
    code, log = invoke(argv(sim, out, seed, "--executor", "local", *extra))
    assert code == 0, log[-4000:]
    return out


def files(out: Path, names: Sequence[str] = SEVEN) -> dict[str, bytes]:
    got: dict[str, bytes] = {}
    for rel in names:
        p = out / rel
        assert p.is_file(), f"{rel} not published"
        got[rel] = p.read_bytes()
    return got


def resolved_config(out: Path) -> Any:
    """The run's resolved config (``run_manifest.json`` > ``config``)."""
    from giap.config import GiapConfig

    data = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    return GiapConfig.model_validate(data["config"])


# --------------------------------------------------------------------------- golden


def mask_version(data: bytes) -> bytes:
    """``##source=giap-<version>`` -> ``##source=giap-VERSION`` (VCF header only)."""
    out = []
    for ln in data.split(b"\n"):
        if ln.startswith(b"##source=giap-"):
            ln = b"##source=giap-VERSION"
        out.append(ln)
    return b"\n".join(out)


def produced(out: Path) -> dict[str, bytes]:
    """The run's bytes per golden logical name (gzip decompressed, VCF version masked)."""
    got: dict[str, bytes] = {}
    for rel in GOLDEN_NAMES:
        raw = (out / rel).read_bytes()
        data = gzip.decompress(raw) if rel.endswith(".gz") else raw
        got[rel] = mask_version(data) if rel.endswith(".vcf") else data
    return got


def _close(a: str, b: str) -> bool:
    if a == b:
        return True
    try:
        return math.isclose(float(a), float(b), rel_tol=REL_TOL, abs_tol=ABS_TOL)
    except ValueError:
        return False


def tsv_equal(want: bytes, got: bytes, loose: Sequence[str] = ()) -> bool:
    """Byte equality, except cells of the ``loose`` columns compared with ``math.isclose``."""
    if not loose:
        return want == got
    w, g = want.decode().split("\n"), got.decode().split("\n")
    if len(w) != len(g) or not w or w[0] != g[0]:
        return False
    header = w[0].split("\t")
    idx = {header.index(c) for c in loose if c in header}
    for lw, lg in zip(w[1:], g[1:], strict=True):
        cw, cg = lw.split("\t"), lg.split("\t")
        if len(cw) != len(cg):
            return False
        for i, (a, b) in enumerate(zip(cw, cg, strict=True)):
            if not (_close(a, b) if i in idx else a == b):
                return False
    return True


def golden_equal(rel: str, want: bytes, got: bytes) -> bool:
    return tsv_equal(want, got, BLAS_NUMERIC.get(GOLDEN_NAMES[rel], ()))


read_golden = m5.read_golden
write_golden = m5.write_golden
updating = m5.updating
first_diff = m5.first_diff


# --------------------------------------------------------------------------- non-vacuity


def _rows(path: Path) -> list[dict[str, str]]:
    import csv

    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def strain_snvs(sim: Any) -> list[tuple[str, int]]:
    """``(contig, 1-based reference position)`` of every planted SNV."""
    rows = _rows(Path(sim.strain_truth))
    return [(r["contig"], int(r["pos_reference"])) for r in rows if r["type"].lower() == "snv"]


def golden_problems(sim: Any, out: Path) -> list[str]:
    """Why the run is not a useful golden (empty: fine). Checked before freezing and on every
    comparison: a site at a planted SNV, a called novel ORF on a planted one, patch features."""
    bad: list[str] = []
    var = _rows(out / "tables" / "variants.tsv")
    snvs = strain_snvs(sim)
    if not snvs:
        bad.append("the fixture plants no SNV")
    at_snv = [r for r in var
              if any(c == r["contig"] and int(r["codon_start"]) <= p < int(r["codon_start"]) + 6
                     for c, p in snvs)]  # fmt: skip
    if not at_snv:
        bad.append("no variants.tsv row at a planted strain SNV")
    if not any("strain_difference" in r["flags"].split(";") for r in var):
        bad.append("no variants.tsv row flagged strain_difference (--genome-relation related)")
    novel = _rows(out / "tables" / "novel_orfs.tsv")
    called = [r for r in novel if r["called"] == "true"]
    if not called:
        bad.append("no called novel ORF")
    planted = _rows(Path(sim.novel_orfs_truth))
    hit = [p for p in planted for r in called
           if r["contig"] == p["contig"] and r["strand"] == p["strand"]
           and _overlaps(r, p)]  # fmt: skip
    if planted and not hit:
        bad.append("no called novel ORF overlaps a planted one")
    patch = [ln for ln in (out / "tables" / "annotation_patch.gff3").read_text().splitlines()
             if ln and not ln.startswith("#") and ln.split("\t")[2] == "CDS"]  # fmt: skip
    if not patch:
        bad.append("annotation_patch.gff3 has no CDS")
    if not any(e.get("floor_above_alpha") is False for e in class_fdr_events(out)):
        bad.append("no novel-ORF class has a class q floor <= alpha (M7-11b rule 8)")
    pile = gzip.decompress((out / "tables" / "pileup.tsv.gz").read_bytes()).splitlines()
    if len(pile) < 2:
        bad.append("pileup.tsv.gz has no row")
    return bad


def class_fdr_events(out: Path) -> list[dict[str, Any]]:
    """The ``data`` of every ``novel_orfs.class_fdr`` event in ``logs/events.jsonl``."""
    path = out / "logs" / "events.jsonl"
    if not path.is_file():
        return []
    rows = [json.loads(ln) for ln in path.read_text(encoding="utf-8").splitlines() if ln.strip()]
    return [r.get("data") or {} for r in rows if r.get("event") == "novel_orfs.class_fdr"]


def _overlaps(call: Mapping[str, str], truth: Mapping[str, str]) -> bool:
    ends = [int(call["end"])] + [int(call[k]) for k in ("start_longest",) if call[k] != "NA"]
    lo, hi = min(ends), max(ends)
    return int(truth["start"]) <= hi and lo <= int(truth["end"])


# --------------------------------------------------------------------------- corruption hooks


def corrupt_tsv(path: Path, column: str, value: str = "x") -> None:
    """Set ``column`` of the first data row to ``value`` (gzip kept when the file is gzip)."""
    raw = path.read_bytes()
    gz = raw[:2] == b"\x1f\x8b"
    lines = (gzip.decompress(raw) if gz else raw).decode("utf-8").split("\n")
    header = lines[0].split("\t")
    assert len(lines) > 2, f"{path}: no data row to corrupt"
    cells = lines[1].split("\t")
    cells[header.index(column)] = value
    lines[1] = "\t".join(cells)
    data = "\n".join(lines).encode("utf-8")
    path.write_bytes(gzip.compress(data, mtime=0) if gz else data)


BAD_VCF_RECORD = "chr1\tnot_a_position\t.\tA\n"
BAD_GFF3_LINE = "chr1\tGIAP\tCDS\t200\t100\t.\t+\t0\tID=corrupt\n"


def append_line(path: Path, line: str) -> None:
    with path.open("a", encoding="utf-8") as fh:
        fh.write(line)


def after(writer: Callable[..., Any], damage: Callable[[Path], None]) -> Callable[..., Any]:
    """``writer`` followed by ``damage(path)`` on the path it wrote (its first argument)."""

    def wrapped(path: Any, *a: Any, **kw: Any) -> Any:
        res = writer(path, *a, **kw)
        damage(Path(path))
        return res

    return wrapped


def validate_tsv(path: Path, stem: str) -> None:
    """Frictionless validation of a written TSV against its schema (spec §9.5)."""
    import pytest

    frictionless = pytest.importorskip("frictionless")
    d = json.loads((SCHEMA_DIR / f"{stem}.schema.json").read_text(encoding="utf-8"))
    d = {k: v for k, v in d.items() if k != "coordinates"}
    res = frictionless.Resource(path=path.name, basepath=str(path.parent), format="tsv",
                                schema=frictionless.Schema.from_descriptor(d))  # fmt: skip
    report = res.validate()
    assert report.valid, (stem, report.flatten(["rowNumber", "fieldName", "type", "note"]))
