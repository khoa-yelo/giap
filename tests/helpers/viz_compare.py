# ruff: noqa: RUF002, RUF003  (ρ is the spec's own character)
"""Builders for the M9-23 V16 ``compare`` tests (docs/specs/M9-23-v16-compare-plot.md).

``inputs(genes)`` turns per-run ``genes.tsv`` frames (``helpers.compare_data.genes_frame``) into
the two ``giap compare`` tables V16 reads, through the real M9-22 functions:
``compare_detection`` (matrix + intersection sections) and ``compare_abundance`` (one row per
unordered pair, TUs of the pair's first run, lead decision 7 of M9-22). Tests that need exact
ρ / CI values overwrite cells of the abundance frame (``set_rho``).

``detect_runs(sets, n_genes)`` builds runs from detected sets ({label: set of gene indices});
every detected gene is measured with a log2 that differs per run, so every pair has a defined
Spearman ρ when it shares at least three genes.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from helpers.compare_data import gene, genes_frame


def detect_runs(
    sets: Mapping[str, set[int]], n_genes: int, *, prefix: str = "g"
) -> dict[str, pd.DataFrame]:
    """``{label: genes frame}``; gene ``i`` is ``f"{prefix}{i + 1}"`` at ``start = 1000 * i``."""
    out: dict[str, pd.DataFrame] = {}
    for r, (label, hit) in enumerate(sets.items()):
        rows = [
            gene(f"{prefix}{i + 1}", start=1000 * i, detected=i in hit,
                 log2=(float(10 + np.sin(i * 1.3) * 3 + 0.25 * r * (i % 4)) if i in hit else None),
                 tu=f"TU_{i // 2}")
            for i in range(n_genes)
        ]  # fmt: skip
        out[label] = genes_frame(rows)
    return out


def inputs(genes: Mapping[str, pd.DataFrame], B: int = 20) -> tuple[pd.DataFrame, pd.DataFrame]:
    """``(compare_detection, compare_abundance)`` as ``giap compare`` builds them."""
    from giap.compare.abundance import abundance_pairs
    from giap.compare.detection import detection_matrix, intersections

    labels = list(genes)
    matrix = detection_matrix(genes)
    detection = pd.concat([matrix, intersections(matrix, labels)], ignore_index=True)
    rng = np.random.default_rng(7)
    pairs = []
    for i, la in enumerate(labels):
        tus = pd.Series(genes[la]["tu_id"].to_numpy(), index=genes[la]["locus_tag"].to_numpy())
        for lb in labels[i + 1 :]:
            pairs.append(abundance_pairs({la: genes[la], lb: genes[lb]}, tus, B, rng))
    return detection, pd.concat(pairs, ignore_index=True)


def set_rho(
    abundance: pd.DataFrame,
    run_a: str,
    run_b: str,
    rho: float | None,
    lo: float | None,
    hi: float | None,
    *,
    flags: str | None = None,
) -> pd.DataFrame:
    """A copy with the (run_a, run_b) row's ρ, CI and flags replaced."""
    out = abundance.copy()
    m = (out["run_a"] == run_a) & (out["run_b"] == run_b)
    assert int(m.sum()) == 1, (run_a, run_b)
    for col, v in (("spearman_rho", rho), ("ci_low", lo), ("ci_high", hi)):
        out[col] = out[col].astype(float)
        out.loc[m, col] = np.nan if v is None else float(v)
    out["flags"] = out["flags"].astype(object)
    out.loc[m, "flags"] = flags
    return out


def intersection_rows(detection: pd.DataFrame) -> pd.DataFrame:
    return detection[detection["section"] == "intersection"].reset_index(drop=True)


def example_1() -> dict[str, pd.DataFrame]:
    """M9-22 example 1: A1 {g1, g2, g3}, A2 {g2, g3, g4}, A3 {g3} on g1-g5."""
    return detect_runs({"A1": {0, 1, 2}, "A2": {1, 2, 3}, "A3": {2}}, 5)


def many_intersections(n_inter: int, n_runs: int) -> dict[str, pd.DataFrame]:
    """``n_inter`` distinct non-empty exact intersections over ``n_runs`` runs; intersection ``j``
    (the ``j``-th non-empty subset in bit order) holds ``j + 1`` genes, so the counts are
    distinct and the 30 largest are unambiguous."""
    labels = [f"R{i + 1:02d}" for i in range(n_runs)]
    masks = list(range(1, 2**n_runs))[:n_inter]
    assert len(masks) == n_inter
    sets: dict[str, set[int]] = {lab: set() for lab in labels}
    g = 0
    for j, mask in enumerate(masks):
        for _ in range(j + 1):
            for b, lab in enumerate(labels):
                if mask >> b & 1:
                    sets[lab].add(g)
            g += 1
    return detect_runs(sets, g)


def labels_of(mask: int, labels: Sequence[str]) -> str:
    return "&".join(lab for b, lab in enumerate(labels) if mask >> b & 1)
