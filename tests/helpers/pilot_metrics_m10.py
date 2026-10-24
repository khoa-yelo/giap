"""Builders for the M10-05 pilot-metrics tests (docs/specs/M10-05-pilot-metrics.md).

Everything here is synthetic: a giapsim SIM100K genome and Bakta JSON (deterministic, built in
``tmp``) stand in for Nissle, and every run directory, peptide table, ``convert.json`` and
``summary.json`` is written by hand. No network, nothing from ``$GIAP_ROOT``, nothing derived
from a DIA-NN report.

``build_pilot(root, world, ...)`` writes the inputs of ``metrics.py compute``:

* ``out_a1``, ``out_a2``: ``tables/genes.tsv``, ``tables/peptide_hits.tsv``,
  ``tables/unmapped.tsv``, ``tables/detection_summary.tsv`` (every table under its schema
  header), ``run_manifest.json`` (built by ``giap.manifest.build_manifest`` and checked against
  ``run_manifest.schema.json``), ``logs/trace.tsv``, ``preflight/preflight.json``;
* ``cmp``: ``compare_detection.tsv``, ``compare_abundance.tsv``, ``compare_stability.tsv``;
* ``pep_a1``, ``pep_a2`` (gzip TSV in the generic peptide-table layout), ``conv_a*``,
  ``sum_a*``, ``k12`` (a small FASTA), ``bakta`` and ``genome``.

The A2 scenario: 20 strain genes each with one unique tryptic peptide (protein_ids = its locus
tag), one decoy row, one ``Cont_``-only peptide. ``te1_misses`` of the 20 peptides have no
placement (``unmapped.tsv`` reason ``no_hit``); the others have one exact ``cds_inframe`` row at
their own gene. A1: the same peptides with UniProt-like accessions and ``search_db = other``; two
have tolerant rows only, one is unmapped, the rest exact, and ``detection_summary.tsv`` carries
``genome_inconsistent_rate`` = 3 / 20.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib
import json
import os
import subprocess
import sys
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[2]
METRICS_PY = REPO / "validation" / "pilot" / "metrics.py"
TEMPLATE = REPO / "validation" / "pilot" / "report_template.md"

N_PEPTIDES = 20
CONT_PEPTIDE = "WWWGGGWWK"
DECOY_PEPTIDE = "KEDITPEPWW"
HOST = "sh04-01n01"


def load_metrics() -> Any:
    """``validation.pilot.metrics`` (a namespace package under the checkout, as te_proxies)."""
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return importlib.import_module("validation.pilot.metrics")


def load_te_proxies() -> Any:
    if str(REPO) not in sys.path:
        sys.path.insert(0, str(REPO))
    return importlib.import_module("validation.proxies.te_proxies")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_cli(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess[str]:
    """``python validation/pilot/metrics.py ARGS`` with ``PYTHONPATH=$REPO/src:$REPO``.

    Fails the calling test when the script is absent: Python's own "can't open file" exit code
    is 2, the code the precondition tests expect."""
    import pytest

    if not METRICS_PY.is_file():
        pytest.fail(f"{METRICS_PY.relative_to(REPO)} does not exist")
    env = dict(os.environ)
    env["PYTHONPATH"] = os.pathsep.join([str(REPO / "src"), str(REPO), str(REPO / "tests" / "sim")])
    return subprocess.run(
        [sys.executable, str(METRICS_PY), *args],
        capture_output=True,
        text=True,
        env=env,
        cwd=str(cwd) if cwd else None,
        check=False,
    )


# ------------------------------------------------------------------ the synthetic genome


@dataclass(frozen=True)
class World:
    genome: Path
    bakta: Path
    genes: list[Any]  # giap.model.genes.Gene (coding genes, read_coding_genes)
    circular: dict[str, int | None]
    regions: dict[str, set[str]]  # region id -> locus tags (giapsim truth)


def sim_world(root: Path) -> World:
    """giapsim SIM100K in ``root`` (no proteomics needed) and its coding genes."""
    sys.path.insert(0, str(REPO / "tests" / "sim"))
    import giapsim

    sim = giapsim.simulate(giapsim.SIM100K, giapsim.SIM100K_SEED, root)
    te = load_te_proxies()
    genes, circ = te.read_coding_genes(Path(sim.bakta_json), Path(sim.genome_fasta))
    regions: dict[str, set[str]] = {}
    with Path(sim.truth["regions"]).open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["genes"] not in ("", "NA"):
                regions[r["region_id"]] = set(r["genes"].split(";"))
    return World(Path(sim.genome_fasta), Path(sim.bakta_json), list(genes), dict(circ), regions)


def unique_peptides(genes: Sequence[Any], n: int) -> list[tuple[str, str]]:
    """``n`` (locus tag, peptide): the first digest product unique to each cds/sorf gene."""
    from giap.mapping.digest import digest

    owners: dict[str, set[str]] = {}
    for g in genes:
        for _, _, p in digest(g.aa):
            owners.setdefault(p.replace("I", "L"), set()).add(g.locus_tag)
    out: list[tuple[str, str]] = []
    for g in sorted(genes, key=lambda x: x.locus_tag):
        if g.feature_type not in ("cds", "sorf"):
            continue
        for _, _, p in digest(g.aa):
            if owners[p.replace("I", "L")] == {g.locus_tag} and not p.startswith("M"):
                out.append((g.locus_tag, p))
                break
        if len(out) == n:
            return out
    raise AssertionError("not enough unique peptides in the synthetic genome")


# ------------------------------------------------------------------ table writers


def _schema_names(stem: str) -> list[str]:
    path = REPO / "src" / "giap" / "schemas" / "tables" / f"{stem}.schema.json"
    return [f["name"] for f in json.loads(path.read_text(encoding="utf-8"))["fields"]]


def write_table(path: Path, stem: str, rows: Sequence[Mapping[str, object]]) -> Path:
    from giap.io.tsv import write_tsv

    cols = _schema_names(stem)
    path.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(path, cols, [[r.get(c) for c in cols] for r in rows])
    return path


def gene_row(g: Any, *, status: str, regions: Mapping[str, set[str]]) -> dict[str, object]:
    """A ``genes.tsv`` row (1-based ``start`` on disk) for a coding gene."""
    reg = sorted(rid for rid, tags in regions.items() if g.locus_tag in tags)
    iv = g.interval
    return {
        "locus_tag": g.locus_tag, "feature_type": g.feature_type, "hypothetical": False,
        "contig": iv.contig, "start": iv.start + 1, "end": iv.end, "strand": iv.strand,
        "length_aa": len(g.aa), "partial": False, "crosses_origin": False,
        "replicon_type": "plasmid" if iv.contig.startswith("plasmid") else "chromosome",
        "tu_id": f"TU_{g.locus_tag}", "region_ids": ";".join(reg) if reg else None,
        "detected": status == "detected", "detection_status": status,
        "abundance_log2": 10.0 if status == "detected" else None,
        "abundance_imputed": False if status == "detected" else None,
        "abundance_type": "maxlfq", "abundance_scale": "relative",
        "annotation_source": "bakta_full",
    }  # fmt: skip


def hit_row(pep: str, g: Any, *, route: str = "exact", cls: str = "cds_inframe") -> dict[str, Any]:
    iv = g.interval
    return {
        "peptide": pep, "peptide_len": len(pep), "source": "database", "best_q_input": 0.001,
        "contig": iv.contig, "start": iv.start + 1, "end": iv.start + 3 * len(pep),
        "strand": iv.strand, "frame": 0, "crosses_origin": False, "route": route,
        "n_diff": 0 if route == "exact" else 1, "placement_q": 0.0, "length_class": "7-40",
        "locus_class": cls, "uniqueness": "unique", "gene": g.locus_tag,
    }  # fmt: skip


def unmapped_row(pep: str, *, label: str = "unmapped", reason: str = "no_hit") -> dict[str, Any]:
    return {"sequence_id": pep, "sequence": pep, "kind": "peptide", "label": label,
            "reason": reason}  # fmt: skip


PEPTIDE_COLS = ("peptide", "modified_peptide", "protein_ids", "q_value", "pep", "is_decoy",
                "intensity", "run", "charge", "source", "search_db")  # fmt: skip


def write_peptides(path: Path, rows: Sequence[Mapping[str, object]]) -> Path:
    """The generic peptide table (M10-03 layout), gzip with a fixed mtime."""
    from giap.io.tsv import format_value

    lines = ["\t".join(PEPTIDE_COLS)]
    lines += ["\t".join(format_value(r.get(c)) for c in PEPTIDE_COLS) for r in rows]
    data = ("\n".join(lines) + "\n").encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fh:
        fh.write(data)
    return path


def pep_row(pep: str, ids: str, *, q: float | None = 0.001, decoy: bool = False,
            run: str = "WT_1", search_db: str = "same") -> dict[str, object]:  # fmt: skip
    return {"peptide": pep, "modified_peptide": pep, "protein_ids": ids, "q_value": q,
            "pep": 0.001, "is_decoy": decoy, "intensity": 1000.0, "run": run, "charge": 2,
            "source": "database", "search_db": search_db}  # fmt: skip


GIB = 1024**3


def write_manifest(
    out: Path,
    *,
    genome_sha: str,
    peptides_sha: str,
    exit_code: int = 0,
    steps: Sequence[Mapping[str, Any]] | None = None,
) -> Path:
    """A schema-valid ``run_manifest.json`` (``giap.manifest.build_manifest``)."""
    from helpers.manifest_data import build, preflight_doc, validate

    pf = preflight_doc(inputs=[
        {"role": "genome", "file": "genome.fasta", "sha256": genome_sha, "size": 1000},
        {"role": "peptides", "file": "peptides.tsv.gz", "sha256": peptides_sha, "size": 100},
    ])  # fmt: skip
    m = build(preflight=pf, exit_code=exit_code, started_at="2026-10-22T10:00:00+00:00",
              finished_at="2026-10-22T10:42:30+00:00")  # fmt: skip
    m["steps"] = [dict(s) for s in (steps if steps is not None else default_steps())]
    validate(m)
    path = out / "run_manifest.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(m, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return path


def default_steps() -> list[dict[str, Any]]:
    """Example 11's steps."""
    return [
        {"name": "mapping", "status": "completed", "n_tasks": 1, "n_retries": 0,
         "wall_time_s": 120.0, "cpu_time_s": 100.0, "peak_rss_bytes": 2 * GIB},
        {"name": "bakta", "status": "cached", "n_tasks": 1, "n_retries": 0},
        {"name": "gene_table", "status": "completed", "n_tasks": 1, "n_retries": 0,
         "wall_time_s": 30.0, "cpu_time_s": 25.0},
    ]  # fmt: skip


def write_trace(out: Path, hosts: Sequence[str] = (HOST,)) -> Path:
    from helpers.manifest_data import trace_row

    rows = [trace_row(f"GIAP:STEP_{i}", hostname=h) for i, h in enumerate(hosts)]
    path = out / "logs" / "trace.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return path


@dataclass(frozen=True)
class Pilot:
    root: Path
    out_a1: Path
    out_a2: Path
    cmp: Path
    pep_a1: Path
    pep_a2: Path
    conv_a1: Path
    conv_a2: Path
    sum_a1: Path
    sum_a2: Path
    bakta: Path
    genome: Path
    k12: Path
    peptides: tuple[tuple[str, str], ...]  # (locus tag, peptide) of the A2 scenario

    def compute_args(self, out: Path, *extra: str, summaries: bool = True) -> list[str]:
        args = [
            "compute", "--a1", str(self.out_a1), "--a2", str(self.out_a2), "--cmp", str(self.cmp),
            "--peptides-a1", str(self.pep_a1), "--peptides-a2", str(self.pep_a2),
            "--convert-a1", str(self.conv_a1), "--convert-a2", str(self.conv_a2),
            "--bakta-json", str(self.bakta), "--genome", str(self.genome),
            "--k12-fasta", str(self.k12), "--out", str(out),
        ]  # fmt: skip
        if summaries:
            args += ["--summary-a1", str(self.sum_a1), "--summary-a2", str(self.sum_a2)]
        return [*args, *extra]


def write_json(path: Path, doc: Mapping[str, Any]) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    return path


def build_pilot(
    root: Path,
    world: World,
    *,
    te1_misses: int = 0,
    q_crossings: int = 0,
    raw_vs_written: int = 0,
    conv_offset: int = 0,
    exit_a1: int = 0,
    genome_sha_a1: str | None = None,
    stability_rows: Sequence[Mapping[str, object]] | None = None,
) -> Pilot:
    """The synthetic pilot (module docstring). ``raw_vs_written``: ``summary.json``'s
    ``distinct_stripped_q01`` is ``peptides.all - raw_vs_written``; ``conv_offset`` is added to
    ``convert.json``'s ``distinct_peptides_target_q01``."""
    by_tag = {g.locus_tag: g for g in world.genes}
    peps = unique_peptides(world.genes, N_PEPTIDES)
    genome_sha = sha256(world.genome)

    # peptide tables
    a2_rows = [pep_row(p, t) for t, p in peps]
    a2_rows += [pep_row(DECOY_PEPTIDE, "rev_" + peps[0][0], q=0.005, decoy=True),
                pep_row(CONT_PEPTIDE, "Cont_P00761", q=0.002)]  # fmt: skip
    a1_rows = [pep_row(p, f"P{i:05d}", search_db="other") for i, (_, p) in enumerate(peps)]
    a1_rows += [pep_row(DECOY_PEPTIDE, "rev_P00000", q=0.005, decoy=True, search_db="other"),
                pep_row(CONT_PEPTIDE, "Cont_P00761", q=0.002, search_db="other")]  # fmt: skip
    pep_a2 = write_peptides(root / "inputs" / "a2.tsv.gz", a2_rows)
    pep_a1 = write_peptides(root / "inputs" / "a1.tsv.gz", a1_rows)
    n_all = N_PEPTIDES + 1  # the strain peptides and the contaminant (targets, q <= 0.01)
    convs, sums = [], []
    for arm in ("a1", "a2"):
        counts = {"distinct_peptides_target_q01": n_all + conv_offset, "q_crossings": q_crossings}
        convs.append(write_json(root / "inputs" / f"{arm}.convert.json", {"counts": counts}))
        sums.append(write_json(root / "diann" / arm / "summary.json",
                               {"distinct_stripped_q01": n_all - raw_vs_written}))  # fmt: skip

    # run directories
    missed = {p for _, p in peps[:te1_misses]}
    detected_a2 = {t for t, p in peps if p not in missed}
    a1_tolerant = {p for _, p in peps[:2]}
    a1_unmapped = {peps[2][1]}
    detected_a1 = {t for t, p in peps if p not in a1_tolerant | a1_unmapped}
    out_a2, out_a1 = root / "runs" / "out_a2", root / "runs" / "out_a1"
    for out, detected in ((out_a2, detected_a2), (out_a1, detected_a1)):
        rows = [gene_row(g, status="detected" if g.locus_tag in detected else "not_detected",
                         regions=world.regions) for g in world.genes]  # fmt: skip
        write_table(out / "tables" / "genes.tsv", "genes", rows)

    hits2 = [hit_row(p, by_tag[t]) for t, p in peps if p not in missed]
    unm2 = [unmapped_row(p) for p in sorted(missed)]
    unm2.append(unmapped_row(CONT_PEPTIDE, label="contaminant", reason="contaminant_id"))
    write_table(out_a2 / "tables" / "peptide_hits.tsv", "peptide_hits", hits2)
    write_table(out_a2 / "tables" / "unmapped.tsv", "unmapped", unm2)
    hits1 = [hit_row(p, by_tag[t], route="tolerant" if p in a1_tolerant else "exact")
             for t, p in peps if p not in a1_unmapped]  # fmt: skip
    unm1 = [unmapped_row(p) for p in sorted(a1_unmapped)]
    unm1.append(unmapped_row(CONT_PEPTIDE, label="contaminant", reason="contaminant_id"))
    write_table(out_a1 / "tables" / "peptide_hits.tsv", "peptide_hits", hits1)
    write_table(out_a1 / "tables" / "unmapped.tsv", "unmapped", unm1)
    qc = {"section": "qc", "group_type": "run", "group": "all",
          "metric": "genome_inconsistent_rate"}  # fmt: skip
    write_table(out_a1 / "tables" / "detection_summary.tsv", "detection_summary",
                [{**qc, "value": 3 / N_PEPTIDES, "n": N_PEPTIDES}])  # fmt: skip
    write_table(out_a2 / "tables" / "detection_summary.tsv", "detection_summary",
                [{**qc, "value": None, "n": None}])  # fmt: skip

    from helpers.compare_data import write_preflight

    for out, pep, code, gsha in (
        (out_a2, pep_a2, 0, genome_sha),
        (out_a1, pep_a1, exit_a1, genome_sha_a1 or genome_sha),
    ):
        write_manifest(out, genome_sha=gsha, peptides_sha=sha256(pep), exit_code=code)
        write_trace(out)
        write_preflight(out, gsha)

    # compare directory
    cmp = root / "runs" / "cmp_a1_a2"
    det_rows = []
    for g in world.genes:
        for label, detected in (("A1", detected_a1), ("A2", detected_a2)):
            d = g.locus_tag in detected
            det_rows.append({"section": "matrix", "locus_tag": g.locus_tag, "run": label,
                             "detected": d,
                             "detection_status": "detected" if d else "not_detected"})  # fmt: skip
    write_table(cmp / "compare_detection.tsv", "compare_detection", det_rows)
    write_table(cmp / "compare_abundance.tsv", "compare_abundance", [ABUNDANCE_ROW])
    write_table(cmp / "compare_stability.tsv", "compare_stability",
                stability_rows if stability_rows is not None else STABILITY_ROWS)  # fmt: skip

    # K-12 side: the proteins of the first five strain genes are "shared"
    k12 = root / "k12" / "UP000000625.fasta"
    k12.parent.mkdir(parents=True, exist_ok=True)
    k12.write_text("".join(f">sp|P9{i:04d}|K12_{i}\n{by_tag[t].aa}\n"
                           for i, (t, _) in enumerate(peps[:5])), encoding="utf-8")  # fmt: skip
    return Pilot(root, out_a1, out_a2, cmp, pep_a1, pep_a2, convs[0], convs[1], sums[0], sums[1],
                 world.bakta, world.genome, k12, tuple(peps))  # fmt: skip


ABUNDANCE_ROW: dict[str, object] = {
    "run_a": "A1", "run_b": "A2", "n_shared": 17, "n_units": 17, "spearman_rho": 0.83,
    "ci_low": 0.76, "ci_high": 0.88, "ba_bias_log2": 0.1, "ba_loa_low": -0.9,
    "ba_loa_high": 1.1, "method": "spearman_tu_bootstrap", "n_boot": 1000, "flags": None,
}  # fmt: skip

#: Example 10's stability rows (islands J 0.5; operon flags both empty; function module off in A1)
STABILITY_ROWS: list[dict[str, object]] = [
    {"family": "islands", "run_a": "A1", "run_b": "A2", "n_a": 10, "n_b": 8, "n_shared": 6,
     "jaccard": 0.5, "flags": None},
    {"family": "operon_flags", "run_a": "A1", "run_b": "A2", "n_a": 0, "n_b": 0, "n_shared": 0,
     "jaccard": None, "flags": "both_empty"},
    {"family": "enriched_categories", "run_a": "A1", "run_b": "A2", "n_a": None, "n_b": 3,
     "n_shared": None, "jaccard": None, "flags": "missing_in:A1"},
    {"family": "novel_orfs", "run_a": "A1", "run_b": "A2", "n_a": 2, "n_b": 4, "n_shared": 1,
     "jaccard": 0.2, "flags": None},
    {"family": "alt_sites", "run_a": "A1", "run_b": "A2", "n_a": 0, "n_b": 0, "n_shared": 0,
     "jaccard": None, "flags": "both_empty"},
]  # fmt: skip


def write_a3(root: Path, genome_sha: str, exit_code: int = 0) -> Path:
    """An A3 run directory with only its manifest (the checklist's first bullet reads it)."""
    out = root / "runs" / "out_a3"
    write_manifest(out, genome_sha=genome_sha, peptides_sha="c" * 64, exit_code=exit_code)
    return out


def output_files(out: Path) -> dict[str, bytes]:
    return {p.name: p.read_bytes() for p in sorted(out.glob("*")) if p.is_file()}
