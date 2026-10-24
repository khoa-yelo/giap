"""Smoke check of the wired novel-ORF class FDR on a power-scenario run (docs/specs/
M7-11b-class-fdr-wiring.md acceptance criterion 6; GIAP_SPEC.md §7.6.4, §8.1).

Input: one giapsim ``novel_orf_power(seed, "100k")`` run (planted ORFs of the four tagged
classes, makedb six-frame ``rev_GIAPORF_...`` decoys in the database peptide table), through
``giap run`` with the local executor and the database peptides (no de novo input).

Claims (the three bullets of criterion 6):

* the ``intergenic``, ``cds_outframe``, ``antisense`` and ``nterm_extension`` ORFs have no
  ``unestimated`` ``class_fdr`` other than ``no_makedb_pair`` caused by a gene-model miss (zero
  is expected);
* at least one ``novel_orfs.class_fdr`` event has ``estimated: true``;
* every class summary's ``floor_above_alpha`` equals the floor computed from its own counts.

Readings pinned by the test author:

* the ORF's class is the side table's ``class`` (the reconstructed class of ``novel_orfs.tsv``);
  "unestimated" is the ``class_fdr`` cell;
* a "gene-model miss" is an ``nterm_extension`` ORF whose overlapping gene has no
  ``nterm_entries`` record in the run's gene model, or whose record's ``db_id`` is not among the
  searched protein IDs (the search database and the run's gene model disagree about that gene);
* "every summary" is every ``novel_orfs.class_fdr`` event in ``logs/events.jsonl``; the
  manifest's ``novel_orfs.classes`` must hold the same summaries. The floor is
  ``1 / n_target_wins`` (picked, None when no win) or ``c / n_targets`` (de novo, None when no
  target); ``floor_above_alpha`` is ``floor > alpha`` with ``alpha`` from the event.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers import novel_steps_fixture as nf

pytestmark = pytest.mark.slow

TAGGED = ("intergenic", "cds_outframe", "antisense", "nterm_extension")
CLASS_EVENT = "novel_orfs.class_fdr"


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory) -> tuple[Any, Path]:
    import giapsim
    import giapsim.scenarios as sc

    root = tmp_path_factory.mktemp("m7_11b_smoke")
    seed = giapsim.SIM100K_SEED
    sim = giapsim.simulate(sc.novel_orf_power(seed, "100k"), seed, root / "sim")
    out = root / "OUT"
    with nf.db_env(root):
        code, log = nf.run_giap(sim, out)
    assert code == 0, log
    return sim, out


def _events(out: Path) -> list[dict[str, Any]]:
    lines = (out / "logs" / "events.jsonl").read_text().splitlines()
    return [
        e for e in (json.loads(ln) for ln in lines if ln.strip()) if e.get("event") == CLASS_EVENT
    ]


def _db_ids(peptides: pd.DataFrame) -> set[str]:
    out: set[str] = set()
    for ids in peptides["protein_ids"]:
        if isinstance(ids, str):
            out.update(i for i in ids.split(";") if i)
        elif ids is not None:
            out.update(str(i) for i in ids)
    return out


def test_tagged_classes_have_no_unestimated_orf_but_gene_model_misses(
    run: tuple[Any, Path],
) -> None:
    from giap.io.fasta import read_genome_fasta
    from giap.makedb.build import nterm_entries
    from giap.steps.mapping import read_gene_model
    from giap.steps.novel import read_candidates

    _, out = run
    work = out / "work"
    side = pd.read_parquet(out / "novel" / "novel_class_fdr.parquet")
    tagged = side[side["class"].isin(TAGGED)]
    assert set(tagged["class"]) == set(TAGGED), sorted(set(side["class"]))
    assert set(side["regime"]) <= {"picked_decoy", "unestimated"}  # database input only
    orfs, _ = read_candidates(work / "novel_candidates")
    gene_of = {o.orf_id: o.overlapping_gene for o in orfs}
    gm = read_gene_model(work / "genes", work / "replicons", work / "features")
    genome = {r.id: r.sequence for r in read_genome_fasta(work / "genome")}
    nterm = {e.locus_tag: e.db_id for e in nterm_entries(gm, genome)}
    ids = _db_ids(pd.read_parquet(work / "peptides"))

    def gene_model_miss(r: dict[str, Any]) -> bool:
        g = gene_of[r["orf_id"]]
        return (r["class"] == "nterm_extension" and r["reason"] == "no_makedb_pair"
                and (g is None or g not in nterm or nterm[g] not in ids))  # fmt: skip

    bad = [r for r in tagged.to_dict("records")
           if r["class_fdr"] == "unestimated" and not gene_model_miss(r)]  # fmt: skip
    assert not bad, f"{len(bad)} of {len(tagged)} tagged-class ORFs unestimated: " + "; ".join(
        f"{r['orf_id']} {r['class']} {r['regime']} {r['reason']}" for r in bad[:20]
    )
    misses = sum(gene_model_miss(r) for r in tagged.to_dict("records"))
    assert misses == 0, f"{misses} nterm_extension gene-model misses (expected zero; M7-13 bug)"


def test_at_least_one_class_event_is_estimated(run: tuple[Any, Path]) -> None:
    evs = _events(run[1])
    assert evs, "no novel_orfs.class_fdr event"
    assert any(e["data"]["estimated"] is True for e in evs), [e["data"] for e in evs]


def _floor(d: dict[str, Any]) -> float | None:
    if d["regime"] == "picked_decoy":
        w = d["n_target_wins"]
        return 1 / w if w else None
    return d["c"] / d["n_targets"] if d["n_targets"] else None


def test_every_summary_floor_matches_its_own_counts(run: tuple[Any, Path]) -> None:
    _, out = run
    evs = [e["data"] for e in _events(out)]
    assert evs
    for d in evs:
        f = _floor(d)
        want = f is not None and f > d["alpha"]
        assert d["floor_above_alpha"] is want, d
        if f is None:
            assert d["q_floor"] is None or math.isnan(d["q_floor"]), d
        else:
            assert d["q_floor"] == pytest.approx(f, rel=1e-12), d
    man = json.loads((out / "run_manifest.json").read_text())
    block = man["novel_orfs"]
    assert block is not None
    key = lambda d: (d["regime"], d["class"])  # noqa: E731
    assert sorted(map(key, block["classes"])) == sorted(map(key, evs))
    for d in block["classes"]:
        f = _floor(d)
        assert d["floor_above_alpha"] is (f is not None and f > d["alpha"]), d
