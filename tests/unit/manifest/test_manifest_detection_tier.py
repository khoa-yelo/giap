"""The manifest's ``detection_tier`` is the tier the run used (review of #359).

Spec: GIAP_SPEC.md §7.3.2 "The tier is chosen once per run from the input and recorded in the
manifest and in ``genes.tsv`` (``detection_tier``)": one tier per run, the same in both places.
docs/specs/M4-24-decoy-unscored.md rules 2-4: ``choose_tier`` may settle on another tier than the
columns suggest (unscored decoys, S = ∅ → the no-decoy chain, ``mayu`` or ``two_peptide``) and
raises ``decoy_unscored``. docs/specs/MVP-04-ecoli-evidence-run.md acceptance row 2 reads the
detection warnings from the manifest's ``warnings`` (backlog: "gene_table logs ... with an event,
so the manifest carries them").

Lead decision (2026-10-25, reversible default): the manifest's ``detection_tier`` comes from the
detection step's output and equals ``tables/genes.tsv`` ``detection_tier`` and the
``tables/detection_summary.tsv`` ``tier`` rows' ``group``. Preflight's ``evidence_tier`` stays a
provisional value in ``preflight.json`` only. No new manifest field. (MVP-03 D-C names
``preflight.json`` ``evidence_tier`` as the source; that row predates M4-24, under which the two
can differ.)

The bug (quickstart demo, Slurm 1894327): preflight said ``decoy_picked``, the decoys carry no
q_value, PEP or score, detection rule 4(c) chose ``mayu``; genes.tsv and detection_summary.tsv say
``mayu``, the manifest said ``decoy_picked`` and had no ``decoy_unscored`` warning.

Interface used (all on main): ``giap.steps.manifest.run_tail(outdir, base, exit_code, ...)``
writes ``OUT/run_manifest.json`` from the files of ``OUT``; the tables are written with the real
writers (``write_genes_tsv``, ``write_detection_summary``) from the producers' output on the
M4-15 hand fixture (``tests/helpers/gene_table_fixture.py``). The tests pin no new source field
or file name: the run directory holds what a finished run publishes.
"""

from __future__ import annotations

import json
import logging
import math
import sys
from collections.abc import Callable, Iterator
from pathlib import Path
from typing import Any

import pandas as pd
import pytest
from helpers.manifest_data import preflight_doc, spec, validate

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "helpers"))
import gene_table_fixture as gt

UNSCORED = "decoy_unscored"
PREFLIGHT_TIER = "decoy_picked"  # T-PRE-18: decoy rows with protein_ids exist

# ----------------------------------------------------------------------------- fixtures


def _unscored(p: pd.DataFrame) -> None:
    """Decoys without q_value, PEP and score (the demo's shape): S = ∅ → ``mayu``."""
    for c in ("q_value", "pep", "score"):
        p.loc[p["is_decoy"], c] = math.nan


def _unscored_no_target_q(p: pd.DataFrame) -> None:
    """As above, and the targets carry only a score: S = ∅ → ``two_peptide``."""
    _unscored(p)
    for c in ("q_value", "pep"):
        p[c] = math.nan


def _scored(p: pd.DataFrame) -> None:
    """The decoy carries q like the targets (q 0.003 ≤ q_value_max): ``decoy_picked``."""
    p.loc[p["is_decoy"], "q_value"] = 0.003


SHAPES: dict[str, tuple[Callable[[pd.DataFrame], None], str]] = {
    "unscored_decoys": (_unscored, "mayu"),
    "unscored_decoys_no_target_q": (_unscored_no_target_q, "two_peptide"),
    "scored_decoys": (_scored, "decoy_picked"),
}


def _world(shape: str) -> Any:
    edit, tier = SHAPES[shape]
    frames = gt.frames_with()
    edit(frames["peptides"])
    w = gt.build_world(frames=frames)
    got = w.detection.table["detection_tier"].astype(str).unique().tolist()
    assert got == [tier], f"fixture {shape}: detection chose {got}, the test expects [{tier!r}]"
    return w


def _write_tables(out: Path, w: Any) -> None:
    """``OUT/tables/genes.tsv`` and ``OUT/tables/detection_summary.tsv`` as a finished run has."""
    from giap.io.genes_tsv import write_genes_tsv
    from giap.io.summary_tsv import write_detection_summary

    (out / "tables").mkdir(parents=True, exist_ok=True)
    write_genes_tsv(out / "tables" / "genes.tsv", w.assemble())
    write_detection_summary(out / "tables" / "detection_summary.tsv", w.detection.summary_rows)


def _base(out: Path, executor: str, tier: str = PREFLIGHT_TIER) -> dict[str, Any]:
    """Launch-time manifest sources (the ``launch_sources`` keys the builder needs)."""
    return {
        "cfg": gt.cfg(),
        "preflight": preflight_doc(tier=tier),
        "registered": [spec("gene_table")],
        "executor": executor,
        "profile": None,
        "resume": False,
        "outdir": str(out),
        "command_line": ["giap", "run"],
        "started_at": "2026-10-25T10:00:00+00:00",
        "input_paths": {"genome": "/data/genome.fa", "peptides": "/data/peptides.tsv"},
    }


def _manifest(out: Path, executor: str, tier: str = PREFLIGHT_TIER) -> dict[str, Any]:
    from giap.steps.manifest import run_tail

    steps = [{"name": "gene_table", "status": "completed", "wall_time_s": 1.0,
              "cpu_time_s": 1.0}] if executor == "local" else []  # fmt: skip
    run_tail(out, _base(out, executor, tier), 0, local_steps=steps)
    doc = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert isinstance(doc, dict)
    return doc


def _genes_tiers(out: Path) -> set[str]:
    df = pd.read_csv(out / "tables" / "genes.tsv", sep="\t", dtype=str, keep_default_na=False)
    return set(df["detection_tier"])


def _summary_tiers(out: Path) -> set[str]:
    df = pd.read_csv(out / "tables" / "detection_summary.tsv", sep="\t", dtype=str,
                     keep_default_na=False)  # fmt: skip
    return set(df.loc[(df["section"] == "tier") & (df["group_type"] == "tier"), "group"])


@pytest.fixture
def gene_table_events(tmp_path: Path) -> Iterator[Callable[[Path, Any], None]]:
    """Run the real ``gene_table`` step body on a fixture world, its log going to
    ``OUT/logs/events.jsonl`` through the run's ``EventsHandler`` (as ``giap run`` does)."""
    from giap.io.events import EventsHandler
    from giap.steps import gene_table as G

    log = logging.getLogger("giap.steps.gene_table")
    handlers: list[logging.Handler] = []
    old_level = log.level

    def run(out: Path, w: Any) -> None:
        h = EventsHandler(out / "logs")
        handlers.append(h)
        log.addHandler(h)
        log.setLevel(logging.INFO)
        mp = pytest.MonkeyPatch()
        try:
            mp.setattr(G, "load_m4_inputs", lambda paths: (w.inputs, []))
            pre = G.build_genes_pre({}, w.cfg, log)
        finally:
            mp.undo()
        got = pre.table["detection_tier"].astype(str).unique().tolist()
        want = w.detection.table["detection_tier"].astype(str).unique().tolist()
        assert got == want, "the step body chose another tier than the fixture world"

    yield run
    for h in handlers:
        log.removeHandler(h)
    log.setLevel(old_level)


# ----------------------------------------------------------------------------- tier


@pytest.mark.parametrize("executor", ["local", "nextflow"])
@pytest.mark.parametrize("shape", ["unscored_decoys", "unscored_decoys_no_target_q"])
def test_manifest_tier_is_the_detection_tier_not_preflights(
    tmp_path: Path, shape: str, executor: str
) -> None:
    """Preflight says ``decoy_picked``; detection settled on another tier → the manifest says
    the detection tier (§7.3.2: one tier, recorded in the manifest and in genes.tsv)."""
    out = tmp_path / "OUT"
    w = _world(shape)
    _write_tables(out, w)
    m = _manifest(out, executor)
    want = SHAPES[shape][1]
    assert _genes_tiers(out) == {want} and _summary_tiers(out) == {want}
    assert m["detection_tier"] == want, (
        f"manifest detection_tier {m['detection_tier']!r} is preflight's provisional "
        f"{PREFLIGHT_TIER!r}, not the tier genes.tsv and detection_summary.tsv record ({want!r})"
    )
    validate(m)


@pytest.mark.parametrize("executor", ["local", "nextflow"])
@pytest.mark.parametrize("shape", list(SHAPES))
def test_manifest_tier_equals_genes_tsv_and_detection_summary(
    tmp_path: Path, shape: str, executor: str
) -> None:
    """The invariant over every fixture shape: manifest == genes.tsv == detection_summary.tsv."""
    out = tmp_path / "OUT"
    _write_tables(out, _world(shape))
    m = _manifest(out, executor)
    assert {m["detection_tier"]} == _genes_tiers(out) == _summary_tiers(out)


@pytest.mark.parametrize("executor", ["local", "nextflow"])
def test_agreeing_tier_is_kept(tmp_path: Path, executor: str) -> None:
    """Preflight and detection agree on ``decoy_picked``: the manifest says ``decoy_picked``."""
    out = tmp_path / "OUT"
    _write_tables(out, _world("scored_decoys"))
    m = _manifest(out, executor)
    assert m["detection_tier"] == "decoy_picked"
    validate(m)


@pytest.mark.parametrize("pre_tier", ["mayu", "two_peptide", "protein_list"])
def test_preflight_tier_never_overrides_the_detection_tier(tmp_path: Path, pre_tier: str) -> None:
    """Whatever preflight provisionally said, the run's tier is the detection step's."""
    out = tmp_path / "OUT"
    _write_tables(out, _world("scored_decoys"))
    m = _manifest(out, "local", tier=pre_tier)
    assert m["detection_tier"] == "decoy_picked", (
        f"manifest copied preflight's {pre_tier!r}; detection used 'decoy_picked'"
    )


def test_preflight_report_keeps_its_provisional_tier(tmp_path: Path) -> None:
    """The provisional tier is preflight's own record, untouched by the manifest writer."""
    out = tmp_path / "OUT"
    pre = out / "preflight" / "preflight.json"
    pre.parent.mkdir(parents=True)
    pre.write_text(json.dumps(preflight_doc(tier=PREFLIGHT_TIER)), encoding="utf-8")
    before = pre.read_bytes()
    _write_tables(out, _world("unscored_decoys"))
    m = _manifest(out, "local")
    assert m["detection_tier"] == "mayu"
    assert pre.read_bytes() == before
    assert "evidence_tier" not in m, "no new manifest field (lead decision)"


# ----------------------------------------------------------------------------- warnings


@pytest.mark.parametrize("executor", ["local", "nextflow"])
def test_decoy_unscored_reaches_the_manifest_warnings(
    tmp_path: Path, executor: str, gene_table_events: Callable[[Path, Any], None]
) -> None:
    """The demo case end to end: the gene_table step's ``decoy_unscored`` reaches
    ``run_manifest.json`` ``warnings`` (MVP-04 acceptance 2 reads it there), and the tier is
    ``mayu``."""
    out = tmp_path / "OUT"
    w = _world("unscored_decoys")
    assert UNSCORED in w.detection.warnings
    gene_table_events(out, w)
    _write_tables(out, w)
    m = _manifest(out, executor)
    hits = [x for x in m["warnings"] if UNSCORED in x]
    assert hits, f"no {UNSCORED!r} among the manifest warnings {m['warnings']}"
    assert m["warnings"] == sorted(set(m["warnings"])), "warnings stay sorted and unique"
    assert m["detection_tier"] == "mayu"
    validate(m)


def test_scored_decoys_raise_no_decoy_unscored_warning(
    tmp_path: Path, gene_table_events: Callable[[Path, Any], None]
) -> None:
    """No false warning: scored decoys (S = C_T) give no ``decoy_unscored`` in the manifest."""
    out = tmp_path / "OUT"
    w = _world("scored_decoys")
    assert UNSCORED not in w.detection.warnings
    gene_table_events(out, w)
    _write_tables(out, w)
    m = _manifest(out, "local")
    assert not [x for x in m["warnings"] if UNSCORED in x], m["warnings"]
