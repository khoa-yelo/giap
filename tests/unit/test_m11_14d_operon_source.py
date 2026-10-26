"""M11-14d: the ``operon_pairs`` reasons and flags, the manifest ``operon_source*`` and the report.

docs/specs/M11-14-container-install.md item 1 option (c) "in detail (14d)": "``operon_pairs``
takes its reason from the resolved config, not from the missing input" (three reasons), the run
flags ``operon_source_distance_model`` and ``layer_missing:uniop``, "Manifest (schema minor) gains
``operon_source`` ... and ``operon_source_reason``, copied from ``operon_meta.json``", and the
"Report" sentence; examples 8 and 9b; AC 8 (each of the three reasons in its case; the report
fixture with the numbers from ``calibration.yaml``).

Interface pinned here (beyond the spec text; see also ``test_m11_14d_uniop_enabled.py``):

- ``operon_meta.json`` keeps its keys ``operon_source`` and ``reason``; the reasons are exactly
  ``UniOP not installed`` (``auto`` resolved to false), ``UniOP disabled
  (annotation.uniop.enabled=false)`` (an explicit false) and ``UniOP failed (...)`` (UniOP ran
  and produced no pair; today's text). The step reads the resolution from its config: the tests
  feed it the ``resolved_config.yaml`` that ``write_params`` writes for
  ``resolve_uniop(...).cfg``, reloaded with ``resolve_config`` (as ``giap _step`` does).
- Both run flags are carried by the step's §6.7 events (``data.flag`` or ``data.flags``), so the
  manifest's ``flags`` and ``layers_missing`` pick them up unchanged.
- ``giap.manifest.ManifestSources`` gains ``operon_meta: Mapping | None = None``, the parsed
  ``operon_meta.json``. ``build_manifest`` writes the top-level ``operon_source`` and
  ``operon_source_reason`` (``None`` without it). ``run_tail`` reads it like ``annotation_meta``:
  Nextflow from the ``operon_pairs`` task directory of the trace, local from
  ``OUT/work/operon_meta``. The schema requires both keys and allows ``uniop`` |
  ``distance_model`` | ``precomputed`` | null.
- The report reads ``operon_source``/``operon_source_reason`` from ``run_manifest.json`` and the
  numbers through ``giap.data.calibration.load()``, called at render time as a module attribute
  (the tests patch it). The 14s2 keys are flat under ``operon_distance`` (PR #382):
  ``precision_at_pair_threshold`` with ``_ci_low``/``_ci_high``, and the same for
  ``recall_at_pair_threshold``. Numbers are written with three decimals, as in ADR-0061 and the
  14s2 addendum ("precision 0.828 [0.805, 0.847]"). The sentence goes into the Summary and the
  Methods sections.
"""

from __future__ import annotations

import csv
import importlib
import json
import logging
from collections.abc import Iterator, Mapping
from pathlib import Path
from types import ModuleType
from typing import Any

import pytest

import giapsim
from giap.config import GiapConfig, resolve_config

REASON_NOT_INSTALLED = "UniOP not installed"
REASON_DISABLED = "UniOP disabled (annotation.uniop.enabled=false)"
REASON_FAILED_PREFIX = "UniOP failed ("
DISTANCE_FLAG = "operon_source_distance_model"
LAYER_FLAG = "layer_missing:uniop"

PARSE_OUT = {
    "genes": "genes.parquet", "features": "features.parquet", "replicons": "replicons.parquet",
    "rrna16s": "rrna16s.fa", "trna_counts": "trna_counts.tsv", "xrefs": "xrefs.parquet",
    "annotation_meta": "annotation_meta.json",
}  # fmt: skip
UNIOP_IN_OUT = {"uniop_cds": "uniop_cds.gff", "uniop_map": "uniop_map.tsv"}
PAIRS_OUT = {"pairs_raw": "pairs_raw.parquet", "tus": "tus.parquet", "tu_upstream": "tu_upstream",
             "operon_meta": "operon_meta.json"}  # fmt: skip


# ---------------------------------------------------------------- step harness (as M2-12 tests)


def _fresh_steps(mp: pytest.MonkeyPatch) -> ModuleType:
    from giap.steps import load_step_modules
    from giap.steps import registry as reg

    names = load_step_modules()
    mp.setattr(reg, "_REGISTRY", {})
    for n in names:
        if n != "registry":
            importlib.reload(importlib.import_module(f"giap.steps.{n}"))
    return reg


@pytest.fixture
def steps() -> Iterator[ModuleType]:
    with pytest.MonkeyPatch.context() as mp:
        yield _fresh_steps(mp)


@pytest.fixture(scope="module")
def work(tmp_path_factory: pytest.TempPathFactory) -> dict[str, Path]:
    """SIM100K after ``parse_bakta`` and ``uniop_input``: the inputs ``operon_pairs`` needs."""
    out = tmp_path_factory.mktemp("sim100k_14d")
    res = giapsim.simulate(giapsim.SIM100K, giapsim.SIM100K_SEED, out / "sim")
    w = out / "work"
    avail = {k: w / v for k, v in PARSE_OUT.items()}
    with pytest.MonkeyPatch.context() as mp:
        reg = _fresh_steps(mp)
        reg.run_step("parse_bakta", GiapConfig(seed=1),
                     {"bakta_json": res.bakta_json, "genome": res.genome_fasta}, avail)  # fmt: skip
        ui = {k: w / v for k, v in UNIOP_IN_OUT.items()}
        reg.run_step("uniop_input", GiapConfig(seed=1),
                     {"genes": avail["genes"], "replicons": avail["replicons"]}, ui)  # fmt: skip
    return {**avail, **ui, "genome": res.genome_fasta}


def _fake_pred(gff: Path, path: Path, p: float) -> Path:
    """UniOP 3ed14dfc's ``uniop.pred`` for ``gff`` with every same-strand pair at ``p``."""
    rows = []
    for ln in gff.read_text().splitlines():
        if ln and not ln.startswith("#"):
            f = ln.split("\t")
            rows.append((f[0], int(f[3]), f[6]))
    order = sorted(range(len(rows)), key=lambda i: rows[i][1])
    body = []
    for k, i in enumerate(order[:-1], start=1):
        nxt = next((rows[j] for j in order[k:] if rows[j][0] == rows[i][0]), None)
        cell = f"{p:.6f} " if nxt is not None and nxt[2] == rows[i][2] else " " * 9
        body.append(f"{k}\t{k + 1}\t{cell}")
    path.write_text("\n".join(["operon prediction scores", "Gene A\tGene B\tPrediction", *body]))
    return path


class _Recorder(logging.Handler):
    def __init__(self) -> None:
        super().__init__(level=logging.DEBUG)
        self.records: list[logging.LogRecord] = []

    def emit(self, record: logging.LogRecord) -> None:
        self.records.append(record)


def _run_pairs(
    reg: ModuleType, work: dict[str, Path], out: Path, cfg: GiapConfig, pred: Path | None = None
) -> tuple[dict[str, Any], set[str]]:
    """``operon_pairs`` → (``operon_meta.json``, run flags of its events)."""
    spec = reg.get("operon_pairs")
    ins = {k: work[k] for k in spec.inputs}
    if pred is not None:
        ins["uniop_pred"] = pred
    outs = {k: out / v for k, v in PAIRS_OUT.items()}
    rec = _Recorder()
    log = logging.getLogger("giap")
    old = log.level
    log.addHandler(rec)
    log.setLevel(logging.DEBUG)
    try:
        reg.run_step("operon_pairs", cfg, ins, outs)
    finally:
        log.removeHandler(rec)
        log.setLevel(old)
    flags: set[str] = set()
    for r in rec.records:
        data = getattr(r, "data", None)
        if isinstance(data, Mapping):
            if isinstance(data.get("flag"), str):
                flags.add(data["flag"])
            flags |= {f for f in data.get("flags") or [] if isinstance(f, str)}
    meta = json.loads(outs["operon_meta"].read_text(encoding="utf-8"))
    assert isinstance(meta, dict)
    return meta, flags


def _step_cfg(tmp: Path, enabled: bool | str, installed: bool) -> GiapConfig:
    """The config ``giap _step`` sees: ``resolve_uniop`` → ``write_params`` → reload."""
    from giap.preflight.environment import resolve_uniop  # type: ignore[attr-defined]
    from giap.runner.nextflow import write_params

    cfg = resolve_config(overrides=[f"annotation.uniop.enabled={str(enabled).lower()}"])
    res = resolve_uniop(cfg, installed=installed, publication=False)
    assert res.level != "error", res
    _, resolved = write_params(tmp / "OUT", res.cfg, {})
    return resolve_config(resolved)


# ---------------------------------------------------------------- the three reasons (AC 8)


def test_auto_resolved_false_reason_is_not_installed(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path
) -> None:
    cfg = _step_cfg(tmp_path, "auto", installed=False)
    meta, flags = _run_pairs(steps, work, tmp_path / "op", cfg)
    assert meta["operon_source"] == "distance_model"
    assert meta["reason"] == REASON_NOT_INSTALLED
    assert {DISTANCE_FLAG, LAYER_FLAG} <= flags, flags


@pytest.mark.parametrize("installed", [False, True])
def test_explicit_false_reason_is_disabled(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path, installed: bool
) -> None:
    cfg = _step_cfg(tmp_path, False, installed=installed)
    meta, flags = _run_pairs(steps, work, tmp_path / "op", cfg)
    assert meta["operon_source"] == "distance_model"
    assert meta["reason"] == REASON_DISABLED
    assert {DISTANCE_FLAG, LAYER_FLAG} <= flags, flags


def test_explicit_false_from_the_user_config_without_preflight(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path
) -> None:
    """A config that says ``false`` itself (as a user's resolved_config would) is "disabled"."""
    cfg = resolve_config(overrides=["annotation.uniop.enabled=false"])
    meta, _ = _run_pairs(steps, work, tmp_path / "op", cfg)
    assert meta["reason"] == REASON_DISABLED


def test_explicit_false_ignores_a_stray_uniop_pred(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path
) -> None:
    cfg = _step_cfg(tmp_path, False, installed=True)
    pred = _fake_pred(work["uniop_cds"], tmp_path / "uniop.pred", 0.9)
    meta, flags = _run_pairs(steps, work, tmp_path / "op", cfg, pred)
    assert meta["operon_source"] == "distance_model"
    assert meta["reason"] == REASON_DISABLED
    assert {DISTANCE_FLAG, LAYER_FLAG} <= flags


@pytest.mark.parametrize("content", ["", "operon prediction scores\nGene A\tGene B\tPrediction\n"])
def test_enabled_true_uniop_ran_without_pairs_is_failed(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path, content: str
) -> None:
    cfg = _step_cfg(tmp_path, True, installed=True)
    pred = tmp_path / "uniop.pred"
    pred.write_text(content)
    meta, flags = _run_pairs(steps, work, tmp_path / "op", cfg, pred)
    assert meta["operon_source"] == "distance_model"
    assert meta["reason"].startswith(REASON_FAILED_PREFIX), meta
    assert {DISTANCE_FLAG, LAYER_FLAG} <= flags, flags


def test_enabled_true_without_pred_is_failed(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path
) -> None:
    meta, flags = _run_pairs(steps, work, tmp_path / "op", GiapConfig(seed=1))
    assert meta["reason"].startswith(REASON_FAILED_PREFIX), meta
    assert {DISTANCE_FLAG, LAYER_FLAG} <= flags, flags


@pytest.mark.parametrize("enabled", [True, "auto"])
def test_uniop_success_has_no_fallback_flag(
    steps: ModuleType, work: dict[str, Path], tmp_path: Path, enabled: bool | str
) -> None:
    cfg = _step_cfg(tmp_path, enabled, installed=True)
    pred = _fake_pred(work["uniop_cds"], tmp_path / "uniop.pred", 0.9)
    meta, flags = _run_pairs(steps, work, tmp_path / "op", cfg, pred)
    assert meta["operon_source"] == "uniop"
    assert not flags & {DISTANCE_FLAG, LAYER_FLAG}, flags


def test_the_three_reasons_are_distinct() -> None:
    assert len({REASON_NOT_INSTALLED, REASON_DISABLED, REASON_FAILED_PREFIX}) == 3


# ---------------------------------------------------------------- manifest builder and schema


def _build(**over: Any) -> dict[str, Any]:
    from helpers.manifest_data import build

    return build(**over)


@pytest.mark.parametrize(
    ("meta", "source", "reason"),
    [
        ({"operon_source": "distance_model", "reason": REASON_NOT_INSTALLED}, "distance_model",
         REASON_NOT_INSTALLED),
        ({"operon_source": "distance_model", "reason": REASON_DISABLED}, "distance_model",
         REASON_DISABLED),
        ({"operon_source": "uniop", "reason": "UniOP probabilities (pin 3ed14dfc)"}, "uniop",
         "UniOP probabilities (pin 3ed14dfc)"),
        ({"operon_source": "precomputed", "reason": "annotation.precomputed.operon_pairs"},
         "precomputed", "annotation.precomputed.operon_pairs"),
    ],
)  # fmt: skip
def test_manifest_copies_operon_meta(meta: dict[str, str], source: str, reason: str) -> None:
    from helpers.manifest_data import validate

    m = _build(operon_meta=meta)
    assert m["operon_source"] == source
    assert m["operon_source_reason"] == reason
    validate(m)


def test_manifest_without_operon_meta_has_null_source() -> None:
    from helpers.manifest_data import validate

    m = _build()
    assert "operon_source" in m and m["operon_source"] is None
    assert "operon_source_reason" in m and m["operon_source_reason"] is None
    validate(m)


def test_schema_requires_both_keys() -> None:
    from helpers.manifest_data import schema

    req = set(schema()["required"])
    assert {"operon_source", "operon_source_reason"} <= req, sorted(req)


def test_schema_rejects_an_unknown_source() -> None:
    from helpers.manifest_data import schema_errors

    m = _build(operon_meta={"operon_source": "distance_model", "reason": REASON_DISABLED})
    m["operon_source"] = "guess"
    assert schema_errors(m), "operon_source 'guess' validated"


def test_operon_source_is_outside_the_provenance_core_inputs() -> None:
    """The source is a run outcome (not known before launch): the digest of two runs with the
    same inputs and config does not depend on it."""
    a = _build(operon_meta={"operon_source": "uniop", "reason": "x"})
    b = _build(operon_meta={"operon_source": "distance_model", "reason": REASON_NOT_INSTALLED})
    assert a["provenance_digest"] == b["provenance_digest"]


def test_example_8_layers_missing_and_flag_from_the_step_events() -> None:
    from helpers.manifest_data import event

    ev = event("operon_pairs", "operon_pairs.run_flag",
               {"flags": [DISTANCE_FLAG, LAYER_FLAG]}, level="WARNING")  # fmt: skip
    m = _build(events=[ev], operon_meta={"operon_source": "distance_model",
                                         "reason": REASON_NOT_INSTALLED})  # fmt: skip
    assert "uniop" in m["layers_missing"]
    assert DISTANCE_FLAG in m["flags"] and LAYER_FLAG in m["flags"]
    assert m["operon_source"] == "distance_model"
    assert m["operon_source_reason"] == REASON_NOT_INSTALLED


# ---------------------------------------------------------------- run_tail reads operon_meta


def _base(out: Path, executor: str) -> dict[str, Any]:
    from helpers.manifest_data import preflight_doc, spec

    return {
        "cfg": GiapConfig(seed=1), "preflight": preflight_doc(),
        "registered": [spec("operon_pairs")], "executor": executor, "profile": None,
        "resume": False, "outdir": str(out), "command_line": ["giap", "run"],
        "started_at": "2026-10-25T10:00:00+00:00",
    }  # fmt: skip


META = {"operon_source": "distance_model", "reason": REASON_DISABLED}


def test_run_tail_local_reads_work_operon_meta(tmp_path: Path) -> None:
    from giap.steps.manifest import run_tail

    out = tmp_path / "OUT"
    (out / "work").mkdir(parents=True)
    (out / "work" / "operon_meta").write_text(json.dumps(META), encoding="utf-8")
    steps = [{"name": "operon_pairs", "status": "completed", "wall_time_s": 1.0,
              "cpu_time_s": 1.0}]  # fmt: skip
    run_tail(out, _base(out, "local"), 0, local_steps=steps)
    doc = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert doc["operon_source"] == "distance_model"
    assert doc["operon_source_reason"] == REASON_DISABLED


def test_run_tail_nextflow_reads_the_operon_pairs_task_dir(tmp_path: Path) -> None:
    from helpers.manifest_data import trace_row

    from giap.steps.manifest import run_tail

    out = tmp_path / "OUT"
    task = out / "work" / "ab" / "cdef12"
    task.mkdir(parents=True)
    (task / "operon_meta.json").write_text(json.dumps(META), encoding="utf-8")
    row = trace_row("ANNOTATE:OPERON_PAIRS (operon_pairs)", workdir=str(task))
    (out / "logs").mkdir(parents=True)
    with (out / "logs" / "trace.tsv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
        w.writeheader()
        w.writerow(row)
    run_tail(out, _base(out, "nextflow"), 0)
    doc = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert doc["operon_source"] == "distance_model"
    assert doc["operon_source_reason"] == REASON_DISABLED


def test_run_tail_without_operon_meta_is_null(tmp_path: Path) -> None:
    from giap.steps.manifest import run_tail

    out = tmp_path / "OUT"
    (out / "work").mkdir(parents=True)
    run_tail(out, _base(out, "local"), 0)
    doc = json.loads((out / "run_manifest.json").read_text(encoding="utf-8"))
    assert doc["operon_source"] is None and doc["operon_source_reason"] is None


# ---------------------------------------------------------------- report text

#: the 14s2 operating point (PR #382, calibration.yaml operon_distance)
OP_14S2 = {
    "precision_at_pair_threshold": 0.8284518828451883,
    "precision_at_pair_threshold_ci_low": 0.8047682923680792,
    "precision_at_pair_threshold_ci_high": 0.8474587186405749,
    "recall_at_pair_threshold": 0.915606936416185,
    "recall_at_pair_threshold_ci_low": 0.9008000451524601,
    "recall_at_pair_threshold_ci_high": 0.930507500774723,
}
OTHER = {
    "precision_at_pair_threshold": 0.5, "precision_at_pair_threshold_ci_low": 0.4,
    "precision_at_pair_threshold_ci_high": 0.6, "recall_at_pair_threshold": 0.7,
    "recall_at_pair_threshold_ci_low": 0.65, "recall_at_pair_threshold_ci_high": 0.75,
}  # fmt: skip


def _sentence(reason: str, nums: Mapping[str, float]) -> str:
    def f(key: str) -> str:
        return f"{nums[key]:.3f}"

    p, r = "precision_at_pair_threshold", "recall_at_pair_threshold"
    return (
        f"Operons/TUs: intergenic-distance model ({reason}); "
        f"precision {f(p)} [{f(p + '_ci_low')}, {f(p + '_ci_high')}], "
        f"recall {f(r)} [{f(r + '_ci_low')}, {f(r + '_ci_high')}] on RegulonDB K-12 at the pair "
        "threshold; TU-based tests use these units"
    )


@pytest.fixture
def calib(monkeypatch: pytest.MonkeyPatch) -> Any:
    """Patch ``giap.data.calibration.load`` to carry a given operating point."""
    from giap.data import calibration

    real = calibration.load

    def use(nums: Mapping[str, float]) -> None:
        def fake() -> dict[str, Any]:
            data = real()
            data["operon_distance"] = {**data["operon_distance"], **nums}
            return data

        monkeypatch.setattr(calibration, "load", fake)

    return use


def _page(tmp_path: Path, source: str | None, reason: str | None) -> Any:
    from helpers.report_content import manifest, page, write_manifest, write_pairs
    from helpers.report_fixtures import run_out

    out = run_out(tmp_path)
    write_pairs(out)
    write_manifest(out, manifest(operon_source=source, operon_source_reason=reason))
    return page(out)


@pytest.mark.parametrize("reason", [REASON_NOT_INSTALLED, REASON_DISABLED,
                                    "UniOP failed (no uniop.pred)"])  # fmt: skip
@pytest.mark.parametrize("sid", ["summary", "methods"])
def test_report_sentence_with_the_14s2_numbers(
    tmp_path: Path, calib: Any, reason: str, sid: str
) -> None:
    from helpers.report_content import section_text

    calib(OP_14S2)
    text = section_text(_page(tmp_path, "distance_model", reason), sid)
    want = _sentence(reason, OP_14S2)
    assert want in text, f"{want!r} not in the {sid} section: {text[:2000]!r}"
    assert "precision 0.828 [0.805, 0.847], recall 0.916 [0.901, 0.931]" in want


@pytest.mark.parametrize("sid", ["summary", "methods"])
def test_report_numbers_come_from_calibration_yaml(tmp_path: Path, calib: Any, sid: str) -> None:
    from helpers.report_content import section_text

    calib(OTHER)
    text = section_text(_page(tmp_path, "distance_model", REASON_DISABLED), sid)
    assert _sentence(REASON_DISABLED, OTHER) in text, text[:2000]


@pytest.mark.parametrize("source", ["uniop", "precomputed", None])
def test_report_has_no_operon_source_text_otherwise(
    tmp_path: Path, calib: Any, source: str | None
) -> None:
    from helpers.report_content import words

    calib(OP_14S2)
    reason = None if source is None else "x"
    p = _page(tmp_path, source, reason)
    assert "intergenic-distance model" not in words(p.root.text())


def test_example_9b_report_caveat_under_publication(tmp_path: Path, calib: Any) -> None:
    """Example 9b: explicit ``false`` under ``--publication``: the caveat is present."""
    from helpers.report_content import manifest, page, section_text, write_manifest, write_pairs
    from helpers.report_fixtures import run_out

    calib(OP_14S2)
    out = run_out(tmp_path)
    write_pairs(out)
    doc = manifest(operon_source="distance_model", operon_source_reason=REASON_DISABLED,
                   publication=True, flags=[DISTANCE_FLAG, LAYER_FLAG],
                   layers_missing=["uniop"])  # fmt: skip
    write_manifest(out, doc)
    assert _sentence(REASON_DISABLED, OP_14S2) in section_text(page(out), "methods")


# -------------------------------------------------- run tail warning (Khoa 2026-10-26)
# A run-time UniOP failure cannot be seen by preflight: ``giap run`` flags it at the end of the
# run, with one line on stderr naming the reason and the distance model, e.g.
# ``giap run: warning: operons used the distance model: UniOP failed (no uniop.pred)``.

FAKE_NF = """#!/bin/sh
prev=""
for a in "$@"; do
  if [ "$prev" = "-with-trace" ]; then mkdir -p "$(dirname "$a")"; cp "$FAKE_NF_TRACE" "$a"; fi
  prev="$a"
done
exit 0
"""


def _cli_run(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, meta: Mapping[str, str]) -> str:
    import os
    import random
    import stat

    from click.testing import CliRunner
    from helpers.fake_env import fake_giap_root, write_setup_manifest
    from helpers.manifest_data import trace_row

    from giap.cli import cli
    from giap.preflight import environment

    monkeypatch.setattr(environment, "probe_uniop", lambda *a, **k: (0, ""), raising=False)
    monkeypatch.setenv("GIAP_ROOT", str(fake_giap_root(tmp_path / "root")))
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR", raising=False)
    monkeypatch.setenv("GIAP_DB_DIR", str(write_setup_manifest(tmp_path / "db")))
    task = tmp_path / "task" / "ab" / "cdef12"
    task.mkdir(parents=True)
    (task / "operon_meta.json").write_text(json.dumps(dict(meta)), encoding="utf-8")
    row = trace_row("ANNOTATE:OPERON_PAIRS (operon_pairs)", workdir=str(task))
    trace = tmp_path / "trace.tsv"
    with trace.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(row), delimiter="\t")
        w.writeheader()
        w.writerow(row)
    bindir = tmp_path / "bin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_TRACE", str(trace))
    rng = random.Random(3)
    genome = tmp_path / "in" / "genome.fa"
    genome.parent.mkdir(parents=True)
    genome.write_text(
        f">chr1 [topology=circular]\n{''.join(rng.choice('ACGT') for _ in range(6000))}\n"
    )
    prot = tmp_path / "in" / "p.faa"
    prot.write_text(
        ">P1\nM" + "".join(rng.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(150)) + "\n"
    )
    argv = ["run", "--genome", str(genome), "--proteins", str(prot), "--outdir",
            str(tmp_path / "OUT"), "--profile", "farmshare,apptainer",
            "--set", "preflight.min_n50_warn=1",
            "--set", "preflight.min_largest_contig_warn=1"]  # fmt: skip
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    doc = json.loads((tmp_path / "OUT" / "run_manifest.json").read_text(encoding="utf-8"))
    assert doc["operon_source"] == meta["operon_source"], doc.get("operon_source")
    return res.output + (res.stderr if res.stderr_bytes is not None else "")


def test_run_tail_warns_on_a_run_time_uniop_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    reason = "UniOP failed (no uniop.pred)"
    text = _cli_run(tmp_path, monkeypatch, {"operon_source": "distance_model", "reason": reason})
    lines = [ln for ln in text.splitlines() if reason in ln]
    assert lines, text
    assert any("warning" in ln.lower() and "distance model" in ln for ln in lines), lines


def test_run_tail_no_operon_warning_when_uniop_ran(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    meta = {"operon_source": "uniop", "reason": "UniOP probabilities (pin 3ed14dfc)"}
    text = _cli_run(tmp_path, monkeypatch, meta)
    assert "distance model" not in text, text
