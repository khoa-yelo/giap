"""M7-10 calibration: the class error rate is unbiased for the true per-class de novo error rate
when the sample carries no variant (docs/specs/M7-10-class-error-rates.md acceptance 2, rule 5,
lead decision 3; GIAP_SPEC.md §7.5.4, §8.1 row "ALT call is real", §8.4; M7-index "class error
rate bias | M7-10 | variant_f1 without variants | |bias| <= 0.05").

R comes from ``tolerances.yaml``: ``replicates_fast.variants`` / ``replicates_full.variants``
(ADR-0042; never hard-coded here). Bound: ``denovo_precision_max_abs_bias`` (0.05), as the
spec says. Seeds: ``replicate_rng("variants.class_error", i)`` (planted arm:
``variants.class_error.planted``).

Replicate: a 50-kb giapsim genome, M3-19 database identifications (pi_wrong 0.05) and their M3-20
de novo reads as ``--peptides``, giapsim's paired ``reference_peptides.tsv`` as
``--reference-peptides``; preflight (writing C17) -> mapping chain -> the registered ``variants``
step with ``peptide_spectra``. The estimate is the written ``variants.tsv`` value per class.

**Estimand.** The true error rate of class c is the share of erroneous reads (giapsim
``truth/denovo.tsv``: read != true peptide after I/L folding) among *all* de novo reads that
support an ALT of class c (paired or not), the sites and classes taken from the merged M7-02/05
chain. Per class, bias = mean over the replicates where c is estimated of (rate_hat - rate_true);
classes estimated in at least half of the replicates are gated.

**Scenario.** M7-12's ``variant_f1`` scenario is not merged yet; this is its no-variant arm built
from the same giapsim pieces (M3-19 + M3-20, no strain layer). When M7-12 lands, the replicate
should switch to ``giapsim.scenarios.variant_f1(seed)`` without its strain layer.

**Planted arm (recorded, not gated; rule 5).** The same with an M3-18 strain layer (sample
peptides on the reference genome): the measured bias is written to the diagnostics only.
giapsim's reference peptides are the sample's own identifications, so this arm does not yet model
a reference search against a database without the variant (M7-12/M7-13 to refine).
"""

from __future__ import annotations

import json
import shutil
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import _criteria as C
import numpy as np
import pytest
import yaml
from helpers.class_error_chain import (
    build_chain,
    confident_pairs,
    measured_rates,
    oracle_sites,
    peptide_classes,
    read_variants,
    run_step,
    truth_rates,
)

pytestmark = [pytest.mark.calibration, pytest.mark.calsuite("variants")]

SUITE = "variants"
STREAM = "variants.class_error"
MIN_GATED_CLASSES = 3  # substitution (both), block (mass_equivalent), swap at 50 kb


def _max_bias() -> float:
    raw = yaml.safe_load(C.TOLERANCES_PATH.read_text())
    return float(raw["denovo_precision_max_abs_bias"])


def _replicate(seed: int, d: Path, *, planted: bool) -> dict[str, Any]:
    strain = None
    if planted:
        from giapsim.strain import StrainConfig

        strain = StrainConfig(snv_rate_per_kb=5.0, indel_rate_per_kb=0.5, indel_len=(1, 1))
    ch = build_chain(d, seed, strain=strain)
    paths, _ = run_step(ch.step_inputs(), d / "out", ch.cfg)
    est = measured_rates(read_variants(paths["variants"]))
    st, peps = oracle_sites(ch.avail, ch.cfg)
    classes = peptide_classes(st, peps)
    truth_file = ch.res.outdir / "truth" / "denovo.tsv"
    true_all = truth_rates(classes, truth_file)
    true_paired = truth_rates(classes, truth_file, paired=confident_pairs(ch.reference_tsv))
    shutil.rmtree(d, ignore_errors=True)
    out: dict[str, Any] = {}
    for c, (n, e) in true_all.items():
        key = f"{c[0]}|{c[1]}"
        tp = true_paired.get(c)
        out[key] = {
            "est": est[c][0] if c in est else None,
            "n_est": est[c][1] if c in est else 0,
            "true": e / n,
            "n_true": n,
            "true_paired": (tp[1] / tp[0]) if tp and tp[0] else None,
        }
    for c in est:
        assert c in true_all, f"class {c} estimated but has no ALT-supporting de novo read"
    return out


def _summarise(reps: list[Mapping[str, Any]]) -> dict[str, dict[str, Any]]:
    keys = sorted({k for r in reps for k in r})
    out: dict[str, dict[str, Any]] = {}
    for k in keys:
        d = [r[k]["est"] - r[k]["true"] for r in reps if k in r and r[k]["est"] is not None]
        dp = [
            r[k]["est"] - r[k]["true_paired"]
            for r in reps
            if k in r and r[k]["est"] is not None and r[k]["true_paired"] is not None
        ]
        out[k] = {
            "R_estimated": len(d),
            "bias": float(np.mean(d)) if d else None,
            "se": float(np.std(d, ddof=1) / np.sqrt(len(d))) if len(d) > 1 else None,
            "bias_vs_paired_truth": float(np.mean(dp)) if dp else None,
            "mean_true": float(np.mean([r[k]["true"] for r in reps if k in r])),
            "mean_n": float(np.mean([r[k]["n_est"] for r in reps if k in r])),
        }
    return out


def _diagnostics(record: dict[str, Any]) -> None:
    d = C.log_dir()
    d.mkdir(parents=True, exist_ok=True)
    with (d / "variants_class_error_diagnostics.jsonl").open("a") as fh:
        fh.write(json.dumps(record, default=float) + "\n")


def _reps(stream: str, R: int, tmp: Path, *, planted: bool) -> list[dict[str, Any]]:
    return [
        _replicate(
            int(C.replicate_rng(stream, i).integers(0, 2**31 - 1)), tmp / f"r{i}", planted=planted
        )
        for i in C.replicate_range(R)
    ]


def _run(full: bool, tmp: Path) -> None:
    R = C.replicates(SUITE, full)
    max_bias = _max_bias()
    reps = _reps(STREAM, R, tmp / "null", planted=False)
    summary = _summarise(reps)
    planted = _summarise(_reps(f"{STREAM}.planted", R, tmp / "planted", planted=True))
    _diagnostics({"suite": STREAM, "full": full, "R_run": len(reps), "max_abs_bias": max_bias,
                  "no_variants": summary, "planted_recorded": planted})  # fmt: skip
    gated = {k: s for k, s in summary.items() if s["R_estimated"] >= max(2, len(reps) // 2)}
    assert len(gated) >= MIN_GATED_CLASSES, f"too few estimated classes: {summary}"
    bad = {k: round(s["bias"], 4) for k, s in gated.items() if abs(s["bias"]) > max_bias}
    assert not bad, f"|bias| > {max_bias} for classes {bad}; summary {summary}"
    for k, s in planted.items():  # recorded only, but it must be computable
        assert s["R_estimated"] == 0 or s["bias"] is not None, k


def test_replicate_counts_come_from_tolerances() -> None:
    raw = yaml.safe_load(C.TOLERANCES_PATH.read_text())
    assert C.replicates(SUITE, False) == int(raw["replicates_fast"]["variants"])
    assert C.replicates(SUITE, True) == int(raw["replicates_full"]["variants"])


def test_class_error_rate_unbiased_without_variants(tmp_path: Path) -> None:
    _run(False, tmp_path)


@pytest.mark.calibration_full
def test_class_error_rate_unbiased_without_variants_full(tmp_path: Path) -> None:
    _run(True, tmp_path)
