"""Read-through rows of the ``related_other_strain`` suite, split by origin (M7-16b D5, rule 20;
spec test plan "Calibration observer (D5)"; statistician review of the 16b amendment, item 3;
methods note ``placement-fdr.md`` 16a items 6 and 8).

A diagnostic observer, not a gate: no tolerance key, bound, alpha or R is added or changed (R3).
Suite ``placement``, full tier only (R = ``replicates_full.placement``, 1 Mb, Slurm). It replays
the ``placement.sim_related`` replicates of ``test_placement_fdr_sim.py`` (same stream, same
seeds, same ``run_chain`` arguments) and, per replicate:

- counts the rows the methods note counted (0-5 per replicate at 1 Mb): best-cost target rows
  whose window holds a target cell (``holds_cell``), among peptides with ``m_target >= 1`` and
  ``placement_q <= max(levels)``;
- splits them by the suite's own truth (``database_truth``: correct peptides at their gene,
  other-strain wrong peptides at their lifted ortholog) into rows placed at their peptide's
  origin locus and the rest;
- **asserts** that every origin row's site (the annotated stop codon in its window, rule 15
  geometry) reads a stop in the sample and a different, sense codon in the reference at the same
  coordinates: a strain stop-loss. The gippasim strain model is SNV-only, so the coordinates
  agree, and this is implied by the model (rule 20), so it is asserted without a tolerance;
- **reports** the off-origin rows (chance read-through placements, already false placements in
  the FDP; expected about alpha x T_acc(m >= 1) x rho) to
  ``<log_dir>/placement_m3_21_diagnostics.jsonl`` and does not assert on them.

A non-vacuity guard: over all replicates at least one row is counted (the 16a full-tier run
counted 200 / 221 / 224 rows at q 0.01 / 0.05 / 0.10). Replicates 0-6 of the stream count
none (the first rows are at replicate 7; local/stat/m7-16-observer-reconcile.md, job 1891223),
so a reduced-R smoke run of this file must go past replicate 7 (0-9 counts 4 rows); R itself
stays ``replicates_full.placement``.
"""

from __future__ import annotations

import shutil
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pandas as pd
import pytest
from _placement_suite import write_diagnostics
from helpers.placement_sims import hits_of, run_chain
from helpers.readthrough_counts import best_cost_rows, holds_cell

pytestmark = [
    pytest.mark.calibration,
    pytest.mark.calibration_full,
    pytest.mark.calsuite("placement"),
    pytest.mark.calibration_unshardable,
]

SUITE = "placement"
STREAM = "placement.sim_related"  # test_placement_fdr_sim.py::_related
LENGTH = 1_000_000
STOPS = frozenset({"TAA", "TAG", "TGA"})
_COMP = str.maketrans("ACGT", "TGCA")


def _fasta(path: Path) -> dict[str, str]:
    out: dict[str, list[str]] = {}
    name = ""
    for ln in Path(path).read_text().splitlines():
        if ln.startswith(">"):
            name = ln[1:].split()[0]
            out[name] = []
        elif ln:
            out[name].append(ln.strip().upper())
    return {k: "".join(v) for k, v in out.items()}


def _codon(seq: str, a: int, strand: str) -> str:
    c = "".join(seq[(a + i) % len(seq)] for i in range(3))
    return c if strand == "+" else c.translate(_COMP)[::-1]


def _sites_in(row: Any, sites: list[Any], n: int, circular: bool) -> list[Any]:
    """Annotated-stop sites in the row's window (rule 15 geometry; ``start`` may be unwrapped)."""
    out = []
    s, e = int(row.start), int(row.end)
    for site in sites:
        if site.contig != row.contig or site.strand != row.strand:
            continue
        xs = (site.codon_start, site.codon_start + n) if circular else (site.codon_start,)
        if any(s <= x and x + 3 <= e and (x - s) % 3 == 0 for x in xs):
            out.append(site)
    return out


def observe(ch: Any, q_max: float) -> dict[str, Any]:
    from gippa.mapping.sixframe import load_stop_cells

    cells = load_stop_cells(ch.avail["sixframe"], table=11)
    hits = best_cost_rows(hits_of(ch.avail))
    hits = hits[hits["search"] == "target"]
    hits = hits[holds_cell(hits, cells)]
    per = pd.read_parquet(ch.avail["placement"]).drop_duplicates("peptide")
    keep = per[(per["m_target"].astype(int) >= 1) & (per["placement_q"] <= q_max)]
    hits = hits[hits["peptide"].isin(set(keep["peptide"]))]
    sample = _fasta(ch.res.genome_fasta)
    reference = _fasta(Path(ch.res.outdir) / "reference.fasta")
    circ = ch.truth.circular
    origin, off, bad = [], [], []
    for r in hits.itertuples(index=False):
        iv = (str(r.contig), int(r.start), int(r.end), str(r.strand))
        rec = {
            "peptide": str(r.peptide),
            "contig": iv[0],
            "start": iv[1],
            "end": iv[2],
            "strand": iv[3],
        }
        if not ch.truth.overlaps(iv, str(r.peptide)):
            off.append(rec)
            continue
        origin.append(rec)
        n = len(sample[iv[0]])
        found = _sites_in(r, list(cells.sites), n, circ.get(iv[0], True))
        assert found, f"a row holding a cell has no site in its window: {rec}"
        for site in found:
            sc = _codon(sample[site.contig], site.codon_start, site.strand)
            rc_ = _codon(reference[site.contig], site.codon_start, site.strand)
            if not (sc in STOPS and rc_ not in STOPS and sc != rc_):
                bad.append({**rec, "codon_start": site.codon_start, "sample": sc, "reference": rc_})
    return {"origin": origin, "off": off, "bad": bad}


def run_observer(tmp: Path) -> dict[str, Any]:
    from gippasim.proteomics import ProteomicsConfig
    from gippasim.strain import StrainConfig

    levels = C.load_tolerances().levels
    q_max = float(max(levels))
    R = C.replicates(SUITE, True)
    per_rep = []
    bad: list[dict[str, Any]] = []
    off_rows: list[dict[str, Any]] = []
    for i in range(R):
        seed = int(C.replicate_rng(STREAM, i).integers(0, 2**31 - 1))
        d = tmp / f"r{i}"
        ch = run_chain(
            "database",
            seed,
            d,
            length=LENGTH,
            strain=StrainConfig(snv_rate_per_kb=5),
            proteomics=ProteomicsConfig(
                search_db="other", pi_wrong=0.05, wrong_source="other_strain"
            ),
            overrides=["genome_relation=related"],
        )
        got = observe(ch, q_max)
        per_rep.append((len(got["origin"]), len(got["off"])))
        bad += [{**b, "replicate": i} for b in got["bad"]]
        off_rows += [{**b, "replicate": i} for b in got["off"]]
        shutil.rmtree(d, ignore_errors=True)
    o = np.array([a for a, _ in per_rep])
    f = np.array([b for _, b in per_rep])
    rec = {
        "suite": "placement.readthrough_related_observer",
        "R": R,
        "q_max": q_max,
        "origin_total": int(o.sum()),
        "off_origin_total": int(f.sum()),
        "origin_range": [int(o.min()), int(o.max())],
        "off_origin_range": [int(f.min()), int(f.max())],
        "replicates_with_rows": int(((o + f) > 0).sum()),
        "off_origin_rows": off_rows[:50],
        "non_stop_loss_origin_rows": bad[:50],
    }
    write_diagnostics(rec)
    return {"rec": rec, "bad": bad}


def test_related_other_strain_readthrough_rows_at_origin_are_stop_losses(tmp_path: Path) -> None:
    out = run_observer(tmp_path)
    rec = out["rec"]
    assert rec["origin_total"] + rec["off_origin_total"] >= 1, (
        f"vacuous: no accepted read-through row in {rec['R']} replicates"
    )
    assert not out["bad"], (
        f"{len(out['bad'])} origin rows whose site is not a strain stop-loss: {out['bad'][:10]}"
    )
