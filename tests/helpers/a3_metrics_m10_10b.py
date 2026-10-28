"""Builders for the M10-10b tests (metric 6 and hypothesis 4; docs/specs/M10-10-a3-run.md).

Everything is synthetic: every ``novel_orfs.tsv``, ``variants.tsv``, ``mapping_fdr.tsv`` and
``compare_stability.tsv`` frame is written under its schema with GIPPA's own TSV writer and read
back with ``gippa.io.run_view.read_schema_tsv`` (the reader ``metrics.py``'s ``_table`` uses, so
the dtypes are the ones the code sees; coordinates stay 1-based as written). The de novo
peptides are made-up strings (``synthetic_peptides``), never taken from a real run; nothing
comes from ``$GIPPA_ROOT``.

``build_a3(pilot, world)`` extends M10-05's synthetic pilot (``helpers.pilot_metrics_m10``) with
the A3 inputs of ``metrics.py compute --a3``: ``out_a3`` (manifest with the ``genome``,
``peptides`` and ``reference_peptides`` input rows; ``tables/novel_orfs.tsv``,
``tables/variants.tsv``; ``mapping/mapping_fdr.tsv``), ``a3.tsv.gz`` (the Cascadia converter's
columns), ``a3.convert.json`` (10a's ``summary`` layout) and ``cmp_a1_a2_a3/tables/*`` (three
labels). It also writes ``novel_orfs.tsv``/``variants.tsv`` into A1's and A2's run directories.
Examples 9 and 10 of the spec are the A2/A3 ORF and ALT sets.
"""

from __future__ import annotations

import gzip
import itertools
import json
import math
import shutil
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
from scipy.stats import norm

from helpers.pilot_metrics_m10 import (
    ABUNDANCE_ROW,
    STABILITY_ROWS,
    Pilot,
    World,
    sha256,
    write_json,
    write_manifest,
    write_table,
)

REPO = Path(__file__).resolve().parents[2]

NOVEL_CLASSES = (
    "nterm_extension",
    "cterm_extension",
    "intergenic",
    "cds_outframe",
    "antisense",
    "pseudogene",
    "sorf",
)
BUCKETS = ("called", "unestimated", "not_called")
ALT_TYPES = ("substitution", "swap", "block", "frameshift")
CONFIDENCE = ("high", "candidate")
VERDICTS = ("not_supported", "consistent_with_lower_bound", "not_evaluable")


# ------------------------------------------------------------------ independent oracles


def wilson(k: int, n: int, level: float = 0.95) -> tuple[float, float]:
    """The Wilson score interval written out from its textbook formula (an oracle independent of
    ``gippa.stats.cutoffs.wilson_ci``); exact 0 / 1 ends at k = 0 / n."""
    z = float(norm.ppf(1.0 - (1.0 - level) / 2.0))
    p = k / n
    d = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / d
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / d
    return (0.0 if k == 0 else centre - half, 1.0 if k == n else centre + half)


def z_ratio(b: int) -> float:
    """z_B / z of definition 7: Φ⁻¹(1 - 0.025/B) / Φ⁻¹(0.975)."""
    return float(norm.ppf(1.0 - 0.025 / b) / norm.ppf(0.975))


def upper(precision: float, ci_high: float, b: int) -> float:
    """U_b of definition 7."""
    return min(1.0, precision + (ci_high - precision) * z_ratio(b))


# ------------------------------------------------------------------ schema-typed frames

_SERIAL = itertools.count()


def frame(stem: str, rows: Sequence[Mapping[str, object]], tmp: Path) -> pd.DataFrame:
    """``rows`` written as ``<stem>.tsv`` under its schema and read back as metrics.py reads it."""
    from gippa.io.run_view import read_schema_tsv

    path = write_table(tmp / f"{stem}.{next(_SERIAL)}.tsv", stem, rows)
    schema = json.loads(
        (REPO / "src" / "gippa" / "schemas" / "tables" / f"{stem}.schema.json").read_text("utf-8")
    )
    return read_schema_tsv(path, stem, schema)


def orf(
    contig: str,
    strand: str,
    end: int,
    cls: str = "intergenic",
    *,
    called: bool = True,
    class_fdr: str = "0.01",
    lowest: int | None = None,
    n_peptides: int = 2,
    **kw: object,
) -> dict[str, object]:
    """A ``novel_orfs.tsv`` row (1-based ``end``; ``orf_id`` from ``lowest``, default a fixed
    offset below ``end``)."""
    low = lowest if lowest is not None else max(1, end - 300)
    row: dict[str, object] = {
        "orf_id": f"NOVEL_{contig}_{low}_{'p' if strand == '+' else 'm'}",
        "class": cls,
        "contig": contig,
        "strand": strand,
        "frame": 0,
        "start_longest": low,
        "end": end,
        "n_peptides": n_peptides,
        "class_fdr": class_fdr,
        "called": called,
    }
    row.update(kw)
    return row


def variant(
    contig: str,
    codon_start: int,
    strand: str,
    alt_aa: str,
    *,
    i: int = 1,
    alt_type: str = "substitution",
    call: str = "ALT",
    confidence: str = "high",
    class_error_rate: str = "unestimated",
    source: str = "database",
) -> dict[str, object]:
    """A ``variants.tsv`` row (1-based ``codon_start`` as written)."""
    return {
        "variant_id": f"VAR_{i}",
        "contig": contig,
        "codon_start": codon_start,
        "strand": strand,
        "ref_codon": "GAT",
        "ref_aa": "D" if alt_aa != "D" else "N",
        "alt_aa": alt_aa,
        "alt_type": alt_type,
        "call": call,
        "confidence": confidence,
        "n_ref": 0,
        "n_alt": 2,
        "class_error_rate": class_error_rate,
        "source": source,
    }


def stab(
    family: str,
    a: str,
    b: str,
    n_a: int | None,
    n_b: int | None,
    n_shared: int | None,
    jaccard: float | None,
    flags: str | None = None,
) -> dict[str, object]:
    return {
        "family": family,
        "run_a": a,
        "run_b": b,
        "n_a": n_a,
        "n_b": n_b,
        "n_shared": n_shared,
        "jaccard": jaccard,
        "flags": flags,
    }


def fdr_row(
    score_bin: int,
    precision: float | None,
    ci_low: float | None,
    ci_high: float | None,
    *,
    c_value: float = 0.9,
    flags: str | None = None,
) -> dict[str, object]:
    return {
        "section": "denovo_precision",
        "decoy_type": "shuffle",
        "score_bin": score_bin,
        "precision": precision,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "c_value": c_value,
        "flags": flags,
    }


PLACEMENT_ROW: dict[str, object] = {
    "section": "placement",
    "length_class": "7-40",
    "n_diff": 0,
    "targets": 40,
    "decoys": 1,
    "pi_wrong": 1.0,
    "fdr_stratum": 0.025,
    "n_accepted": 40,
    "decoy_type": "reverse",
}


# ------------------------------------------------------------------ examples 9 and 10

C1, C2 = "CP171242.1", "CP171243.1"


def example9_novel() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Example 9: A2 {(C1, +, 1203), (C1, -, 88001)}; A3 {(C1, +, 1203), (C2, +, 3011),
    (C1, -, 5002)}; the shared ORF is ``intergenic`` in A2 and ``nterm_extension`` in A3, with
    ``orf_id`` ``NOVEL_CP171242.1_1001_p`` vs ``NOVEL_CP171242.1_1050_p``."""
    a2 = [
        orf(C1, "+", 1203, "intergenic", lowest=1001),
        orf(C1, "-", 88001, "antisense", lowest=87700),
    ]
    a3 = [
        orf(C1, "+", 1203, "nterm_extension", lowest=1050),
        orf(C2, "+", 3011, "intergenic", lowest=2700),
        orf(C1, "-", 5002, "antisense", lowest=4700),
    ]
    return a2, a3


def example10_variants() -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    """Example 10: A2 {(C1, 4501, +, N)}, A3 {(C1, 4501, +, N), (C1, 7702, -, D)}."""
    a2 = [variant(C1, 4501, "+", "N", i=1)]
    a3 = [
        variant(C1, 4501, "+", "N", i=1, source="denovo"),
        variant(C1, 7702, "-", "D", i=2, source="denovo"),
    ]
    return a2, a3


def example_stability(
    novel_j: float | None = 0.0, alt: Sequence[Any] | None = None
) -> list[dict[str, object]]:
    """``compare_stability`` rows of a three-label compare: (A2, A3) ``novel_orfs`` keyed on
    ``orf_id`` (2, 3, 0 shared: Jaccard 0, example 9) and ``alt_sites`` (1, 2, 1: 0.5, example
    10) unless overridden."""
    rows = list(STABILITY_ROWS)
    for a, b in (("A1", "A3"), ("A2", "A3")):
        for fam in ("islands", "operon_flags", "enriched_categories"):
            rows.append(stab(fam, a, b, 0, 0, 0, None, "both_empty"))
    rows.append(stab("novel_orfs", "A1", "A3", 0, 3, 0, 0.0))
    rows.append(stab("alt_sites", "A1", "A3", 0, 2, 0, 0.0))
    rows.append(stab("novel_orfs", "A2", "A3", 2, 3, 0, novel_j))
    a = list(alt) if alt is not None else [1, 2, 1, 0.5, None]
    rows.append(stab("alt_sites", "A2", "A3", *a))
    return rows


# ------------------------------------------------------------------ peptides


AA = "ADEFGHNQSTVWY"  # no K/R (no internal tryptic site), no I/L, no C/M


def synthetic_peptides(n: int, length: int = 9) -> list[str]:
    """``n`` distinct made-up tryptic-looking peptides of ``length`` residues (ending in K, no
    I/L): equal length, so none is a substring of another."""
    out = []
    for combo in itertools.product(AA, repeat=length - 1):
        out.append("".join(combo) + "K")
        if len(out) == n:
            return out
    raise AssertionError("not enough synthetic peptides")


def queries_frame(peptides: Sequence[str], scores: Sequence[float]) -> pd.DataFrame:
    """The ``denovo_queries`` frame layout ``hypothesis4`` reads: ``peptide``, ``score``
    (one row per peptide), plus ``source = denovo``."""
    return pd.DataFrame(
        {
            "peptide": list(peptides),
            "score": [float(s) for s in scores],
            "source": ["denovo"] * len(peptides),
        }
    )


def equal_bins(n_bins: int, per_bin: int) -> tuple[list[str], list[float], list[list[str]]]:
    """``n_bins x per_bin`` distinct peptides with distinct scores in (0.8, 1]; ``score_bins``
    puts exactly ``per_bin`` in each decile bin (checked). Returns peptides, scores and the
    peptides of each bin."""
    import numpy as np

    from gippa.mapping.denovo_precision import score_bins

    n = n_bins * per_bin
    peps = synthetic_peptides(n)
    scores = [0.8 + 0.2 * (i + 1) / n for i in range(n)]
    bins = score_bins(np.asarray(scores), n_bins)
    by_bin = [[p for p, b in zip(peps, bins, strict=True) if b == j] for j in range(n_bins)]
    assert [len(x) for x in by_bin] == [per_bin] * n_bins, [len(x) for x in by_bin]
    return peps, scores, by_bin


PEPTIDE_COLS_A3 = ("peptide", "modified_peptide", "score", "run", "charge", "spectrum_id", "source")


def write_denovo_table(path: Path, rows: Sequence[Mapping[str, object]]) -> Path:
    """A de novo peptide table in ``gippa convert --from cascadia``'s columns (gzip, mtime 0)."""
    from gippa.io.tsv import format_value

    lines = ["\t".join(PEPTIDE_COLS_A3)]
    lines += ["\t".join(format_value(r.get(c)) for c in PEPTIDE_COLS_A3) for r in rows]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as fh:
        fh.write(("\n".join(lines) + "\n").encode("utf-8"))
    return path


def dn_row(
    pep: str, score: float, *, run: str = "WT_1", charge: int = 2, scan: int = 1
) -> dict[str, object]:
    return {
        "peptide": pep,
        "modified_peptide": pep,
        "score": score,
        "run": run,
        "charge": charge,
        "spectrum_id": f"{run}:scan={scan}",
        "source": "denovo",
    }


# ------------------------------------------------------------------ the CLI fixture


@dataclass(frozen=True)
class A3:
    out_a3: Path
    pep_a3: Path
    conv_a3: Path
    cmp3: Path
    strain_db: Path
    denovo: tuple[str, ...]  # the A3 query peptides (synthetic and A2-confirmed)
    n_confirmed: int

    def args(self, *, skip: Sequence[str] = ()) -> list[str]:
        """``--a3 … --peptides-a3 … --convert-a3 … --cmp3 …`` minus the options in ``skip``."""
        opts = {
            "--a3": self.out_a3,
            "--peptides-a3": self.pep_a3,
            "--convert-a3": self.conv_a3,
            "--cmp3": self.cmp3,
        }
        return [x for k, v in opts.items() if k not in skip for x in (k, str(v))]


N_SYNTH = 20
N_CONFIRMED = 10


def _file_entry(path: Path, rows: int) -> dict[str, object]:
    return {"path": str(path), "sha256": sha256(path), "bytes": path.stat().st_size, "rows": rows}


def write_arm_tables(
    out: Path, novel: Sequence[Mapping[str, object]], variants: Sequence[Mapping[str, object]]
) -> None:
    write_table(out / "tables" / "novel_orfs.tsv", "novel_orfs", novel)
    write_table(out / "tables" / "variants.tsv", "variants", variants)


def build_a3(
    pilot: Pilot,
    world: World,
    *,
    exit_code: int = 0,
    genome_sha: str | None = None,
    peptides_sha: str | None = None,
    reference_sha: str | None = None,
    conv_failures: Sequence[str] = (),
    conv_peptides_sha: str | None = None,
    stability: Sequence[Mapping[str, object]] | None = None,
    fdr_extra_bin: bool = False,
    drop_cmp_table: str | None = None,
) -> A3:
    """The A3 inputs (module docstring). The overrides break one precondition or cross-check
    each: ``reference_sha`` (default: sha256 of A2's peptides table), ``conv_peptides_sha``
    (default: sha256 of ``a3.tsv.gz``), ``stability`` (``cmp3``'s ``compare_stability`` rows),
    ``fdr_extra_bin`` (a ``denovo_precision`` row with a bin label ``score_bins`` never gives),
    ``drop_cmp_table`` (a ``cmp3`` table left out)."""
    from helpers.compare_data import write_preflight

    root = pilot.root
    gsha = sha256(world.genome)
    # A1/A2 tables (examples 9 and 10 for A2; one row each for A1)
    n2, n3 = example9_novel()
    v2, v3 = example10_variants()
    write_arm_tables(pilot.out_a2, n2, v2)
    write_arm_tables(
        pilot.out_a1, [orf(C1, "+", 1203, "intergenic", lowest=1001)], [variant(C1, 4501, "+", "N")]
    )

    # de novo table: 20 synthetic peptides plus 10 of A2's peptides, I written as L (Cascadia)
    a2_peps = [p for _, p in pilot.peptides][:N_CONFIRMED]
    denovo = synthetic_peptides(N_SYNTH) + [p.replace("I", "L") for p in a2_peps]
    rows = [dn_row(p, 0.8 + 0.19 * (i + 1) / len(denovo), scan=i + 1) for i, p in enumerate(denovo)]
    rows.append(dn_row(denovo[0], 0.801, run="WT_2", scan=999))  # a second, lower-scored row
    pep_a3 = write_denovo_table(root / "inputs" / "a3.tsv.gz", rows)
    log = root / "inputs" / "a3.convert.log"
    log.write_text("", encoding="utf-8")
    conv = {
        "failures": list(conv_failures),
        "outputs": {"peptides": _file_entry(pep_a3, len(rows)), "log": _file_entry(log, 0)},
        "counts": {"rows": len(rows), "distinct_peptides": len(denovo)},
    }
    if conv_peptides_sha is not None:
        conv["outputs"]["peptides"]["sha256"] = conv_peptides_sha  # type: ignore[index]
    conv_a3 = write_json(root / "inputs" / "a3.convert.json", conv)

    # out_a3
    out = root / "runs" / "out_a3"
    write_manifest(
        out,
        genome_sha=genome_sha or gsha,
        peptides_sha=peptides_sha or sha256(pep_a3),
        exit_code=exit_code,
    )
    mp = out / "run_manifest.json"
    doc = json.loads(mp.read_text(encoding="utf-8"))
    doc["inputs"].append(
        {
            "role": "reference_peptides",
            "file": "a2.tsv.gz",
            "path": str(pilot.pep_a2),
            "sha256": reference_sha or sha256(pilot.pep_a2),
            "size": 2048,
        }
    )
    mp.write_text(json.dumps(doc, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    write_preflight(out, genome_sha or gsha)
    write_arm_tables(out, n3, v3)
    fdr = [PLACEMENT_ROW] + [fdr_row(b, 0.5, 0.3, 0.7) for b in range(10)]
    if fdr_extra_bin:
        fdr.append(fdr_row(10, 0.5, 0.3, 0.7))
    write_table(out / "mapping" / "mapping_fdr.tsv", "mapping_fdr", fdr)

    # cmp_a1_a2_a3
    cmp3 = root / "runs" / "cmp_a1_a2_a3" / "tables"
    cmp3.mkdir(parents=True, exist_ok=True)
    shutil.copy(pilot.cmp / "tables" / "compare_detection.tsv", cmp3 / "compare_detection.tsv")
    ab = [dict(ABUNDANCE_ROW)]
    for a, b in (("A1", "A3"), ("A2", "A3")):
        ab.append(
            {"run_a": a, "run_b": b, "method": "spearman_tu_bootstrap", "flags": "missing_in:A3"}
        )
    write_table(cmp3 / "compare_abundance.tsv", "compare_abundance", ab)
    write_table(
        cmp3 / "compare_stability.tsv",
        "compare_stability",
        list(stability) if stability is not None else example_stability(),
    )
    if drop_cmp_table is not None:
        (cmp3 / f"{drop_cmp_table}.tsv").unlink()

    # A2's search FASTA stand-in for --strain-db (the proteins of the A2 peptides' genes)
    by_tag = {g.locus_tag: g for g in world.genes}
    strain_db = root / "inputs" / "strain_db.faa"
    strain_db.write_text(
        "".join(f">{t}\n{by_tag[t].aa}\n" for t, _ in pilot.peptides), encoding="utf-8"
    )
    return A3(out, pep_a3, conv_a3, cmp3.parent, strain_db, tuple(denovo), N_CONFIRMED)
