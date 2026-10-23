"""Calibration sanity of the M7-12 generators (issue M7-12 acceptance; GIAP_SPEC.md §7.6.6, §8.4).

The M7-13 gate G3 must test GIAP, not the simulator, so the planted novel ORFs must be callable:

* **RBS support** (spec acceptance): over the planted starts (``rbs_planted = true`` rows of
  ``truth/novel_orfs.tsv``: intergenic starts and N-terminal true starts), the ``tir_ostir`` of
  ``precomputed/rbs_novel.tsv`` at ``(contig, strand, true_start)`` is ≥ the 10 % quantile
  (``numpy.quantile``, linear; M7-09 rule 2, ``novel_orfs.rbs_min_quantile`` default 0.10) of the
  first-round ``tir_ostir`` of detected genes in ≥ 90 % of planted starts, pooled over replicates.
  The detected-gene distribution is built as GIAP builds it: ``precomputed/rbs.tsv`` with the
  last row per ``locus_tag`` (``giap.annotate.layers._by_tag``), restricted to annotated non-pseudo
  CDS/sORFs that are not origin-crossing and have ≥ 1 correct peptide in ``truth/peptides.tsv``
  (the simulator's ground truth stands in for GIAP's ``detected``). Each replicate needs ≥ 20 such
  genes (M7-09 decision 1), else the threshold is unestimated and the test fails.
* **Peptide support** (example 2 "``n_true_peptides ≥ 2`` for most"): pooled over replicates and
  the four classes, more than half of the planted ORFs have ≥ 2 true peptides.

Replicates use seeds ``seeds.calibration_base + r`` (tolerances.yaml). Fast: 4 replicates at
100 kb (each < 2 s, rule 4). Full (``calibration_full``): ``replicates_fast.novel_orfs`` (50)
replicates at 100 kb plus 5 at 1 Mb. The 90 % share is the spec's stated claim, checked as a point
estimate (no SE allowance): the generator is supposed to clear it by design. Suite ``novel_orfs``.
"""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Any

import numpy as np
import pytest
from helpers import m7_sim_oracle as O

import giapsim

pytestmark = pytest.mark.calsuite("novel_orfs")

RBS_Q = 0.10
MIN_SHARE = 0.90
MIN_GENES = 20
N_FAST = 4


def _sim(d: Path, seed: int, size: str) -> Path:
    import giapsim.scenarios as sc

    cfg = dataclasses.replace(sc.novel_orf_power(seed, size), planted_tir=(9.0, 0.5))
    giapsim.simulate(cfg, seed, d)  # the decision 8 sanity arm (decision 14; lead edit)
    return d


def _detected_tirs(out: Path) -> np.ndarray:
    data = json.loads((out / "bakta_like.json").read_text())
    lengths = {s["id"]: int(s["length"]) for s in data["sequences"]}
    ok = set()
    for f in data["features"]:
        if f["type"] not in ("cds", "sorf") or f.get("pseudo"):
            continue
        contig = f.get("contig", f.get("sequence"))
        if int(f["stop"]) < int(f["start"]) or int(f["stop"]) > lengths[contig]:
            continue  # origin-crossing: partial in the gene table
        ok.add(str(f["locus"]))
    detected = {r["locus_tag"] for r in O.correct_peptide_rows(out)} & ok
    tir: dict[str, float] = {}
    for r in O.rows(out / "precomputed" / "rbs.tsv"):
        if r["tir_ostir"] not in ("", "NA"):
            tir[r["locus_tag"]] = float(r["tir_ostir"])  # last row per locus wins
    return np.array([tir[t] for t in sorted(detected) if t in tir], dtype=float)


def _planted_tirs(out: Path) -> list[float]:
    layer = {(r["contig"], r["strand"], int(r["start"])): float(r["tir_ostir"])
             for r in O.rows(out / "precomputed" / "rbs_novel.tsv")}  # fmt: skip
    vals = []
    for r in O.rows(out / "truth" / "novel_orfs.tsv"):
        if O.truthy(r["rbs_planted"]):
            vals.append(layer[(r["contig"], r["strand"], int(r["true_start"]))])
    return vals


def _replicate(out: Path) -> tuple[int, int, int, int]:
    """(planted starts passing the threshold, planted starts, ORFs with >= 2 peptides, ORFs)."""
    det = _detected_tirs(out)
    assert len(det) >= MIN_GENES, f"{len(det)} detected genes with a TIR in {out}"
    thr = float(np.quantile(det, RBS_Q, method="linear"))
    planted = _planted_tirs(out)
    assert planted, "no planted RBS start"
    rows = O.rows(out / "truth" / "novel_orfs.tsv")
    two = sum(int(r["n_true_peptides"]) >= 2 for r in rows)
    return sum(t >= thr for t in planted), len(planted), two, len(rows)


def _check(results: list[tuple[int, int, int, int]]) -> None:
    passed, planted, two, n = (sum(x[i] for x in results) for i in range(4))
    assert passed / planted >= MIN_SHARE, (passed, planted)
    assert two / n > 0.5, (two, n)


@pytest.fixture(scope="module")
def base(tolerances: dict[str, Any]) -> int:
    return int(tolerances["seeds"]["calibration_base"])  # type: ignore[index]


@pytest.mark.calibration
def test_planted_rbs_pass_detected_gene_quantile_fast(
    tmp_path_factory: pytest.TempPathFactory, base: int
) -> None:
    res = [_replicate(_sim(tmp_path_factory.mktemp(f"m7rbs{r}"), base + r, "100k"))
           for r in range(N_FAST)]  # fmt: skip
    _check(res)


@pytest.mark.calibration_full
def test_planted_rbs_pass_detected_gene_quantile_full(
    tmp_path_factory: pytest.TempPathFactory, base: int, tolerances: dict[str, Any]
) -> None:
    n = int(tolerances["replicates_fast"]["novel_orfs"])  # type: ignore[index]
    res = [_replicate(_sim(tmp_path_factory.mktemp(f"m7rbsF{r}"), base + 1000 + r, "100k"))
           for r in range(n)]  # fmt: skip
    res += [_replicate(_sim(tmp_path_factory.mktemp(f"m7rbsM{r}"), base + 2000 + r, "1M"))
            for r in range(5)]  # fmt: skip
    _check(res)
