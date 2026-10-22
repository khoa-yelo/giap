"""Row builders for the M9-05 mapping plots (V5 ``mapping_sankey``, V6 ``mapping_fdr``).

Rows are written with ``helpers.viz_out_tree.write_schema_table`` under the M3 schemas
(``peptide_hits``, ``protein_hits``, ``unmapped``, ``mapping_fdr``); coordinates are the 1-based
closed ones of the published TSVs. Nothing here touches ``$GIAP_DB_DIR`` or the network.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

from helpers.viz_out_tree import make_out, write_schema_table

#: ``mapping.length_classes`` defaults as ``mapping_fdr`` labels them
LENGTH_LABELS = ("7", "8", "9", "10-11", "12-14", "15-19", "20-40")
STRATA = (0, 1, 2)


def peptide_hit(peptide: str, start: int = 101, route: str = "exact", **kw: object) -> dict:
    row: dict[str, object] = {
        "peptide": peptide,
        "peptide_len": len(peptide),
        "source": "database",
        "contig": "chr1",
        "start": start,
        "end": start + 3 * len(peptide) - 1,
        "strand": "+",
        "frame": 1,
        "crosses_origin": False,
        "route": route,
        "n_diff": 0 if route == "exact" else 1,
        "uniqueness": "unique",
    }
    row.update(kw)
    return row


def protein_hit(protein_id: str, route: str, rank: int, start: int = 1001, **kw: object) -> dict:
    row: dict[str, object] = {
        "protein_id": protein_id,
        "protein_len_aa": 100,
        "route": route,
        "contig": "chr1",
        "start": start,
        "end": start + 299,
        "strand": "+",
        "rank": rank,
        "identity": 0.98,
        "query_coverage": 1.0,
        "class": "near",
    }
    row.update(kw)
    return row


def unmapped(
    seq_id: str, reason: str, kind: str = "peptide", label: str = "unmapped", **kw: object
) -> dict:
    row: dict[str, object] = {
        "sequence_id": seq_id,
        "sequence": seq_id,
        "kind": kind,
        "label": label,
        "reason": reason,
    }
    row.update(kw)
    return row


def contaminant(seq_id: str, kind: str = "peptide") -> dict:
    return unmapped(seq_id, "contaminant_id", kind, "contaminant", contaminant_match="CON__P1")


def fdr_row(
    length_class: str,
    n_diff: int,
    targets: int,
    decoys: int,
    fdr: float | None,
    accepted: int,
    decoy_type: str = "reverse",
) -> dict:
    return {
        "length_class": length_class,
        "n_diff": n_diff,
        "targets": targets,
        "decoys": decoys,
        "pi_wrong": 1.0,
        "fdr_stratum": fdr,
        "n_accepted": accepted,
        "decoy_type": decoy_type,
        "section": "placement",
    }


def denovo_row(
    score_bin: int,
    precision: float | None,
    ci_low: float | None,
    ci_high: float | None,
    flags: str | None = None,
) -> dict:
    return {
        "section": "denovo_precision",
        "decoy_type": "shuffle",
        "score_bin": score_bin,
        "precision": precision,
        "ci_low": ci_low,
        "ci_high": ci_high,
        "c_value": 0.99,
        "flags": flags,
    }


def full_placement(
    totals: Mapping[int, tuple[int, int, float | None, int]],
    per_class: Mapping[tuple[str, int], tuple[int, int]] | None = None,
    decoy_type: str = "reverse",
) -> list[dict]:
    """Every length-class row plus the ``all`` row per m, as ``placement_tables`` writes them.

    ``totals[m] = (targets, decoys, fdr_stratum, n_accepted)`` for the ``all`` rows; per-class
    ``(targets, decoys)`` default to 0 (fdr NA, accepted 0)."""
    per_class = per_class or {}
    rows = []
    for m in STRATA:
        for label in LENGTH_LABELS:
            t, d = per_class.get((label, m), (0, 0))
            rows.append(fdr_row(label, m, t, d, (d / t) if t else None, 0, decoy_type))
        t, d, fdr, acc = totals.get(m, (0, 0, None, 0))
        rows.append(fdr_row("all", m, t, d, fdr, acc, decoy_type))
    return rows


#: example 4 (M9-05): all/0, all/1, all/2
EXAMPLE4_TOTALS: dict[int, tuple[int, int, float | None, int]] = {
    0: (900, 3, 0.0033, 880),
    1: (120, 2, 0.025, 40),
    2: (0, 0, None, 0),
}


def mapping_out(
    tmp: Path,
    *,
    hits: Sequence[Mapping[str, object]] | None = None,
    proteins: Sequence[Mapping[str, object]] | None = None,
    unmapped_rows: Sequence[Mapping[str, object]] | None = None,
    fdr: Sequence[Mapping[str, object]] | None = None,
) -> Path:
    """An ``OUT/`` with the given ``mapping/`` tables; ``None`` leaves the file out, ``[]``
    writes it header-only."""
    out = make_out(tmp)
    tables = (
        ("peptide_hits", hits),
        ("protein_hits", proteins),
        ("unmapped", unmapped_rows),
        ("mapping_fdr", fdr),
    )
    for name, rows in tables:
        if rows is not None:
            write_schema_table(out / "mapping" / f"{name}.tsv", name, rows)
    return out
