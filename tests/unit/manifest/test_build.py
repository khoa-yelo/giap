"""MVP-03a examples 1-6 and the run-free edge cases against the pure manifest builder
(docs/specs/MVP-03-manifest-finalize.md D-C, D-E, "Examples", "Edge cases"; GIAP_SPEC.md §6.6,
§6.7, §9.1, §9.4).

Interface assumed: ``giap.manifest.build_manifest(sources: ManifestSources) -> dict`` (no file
I/O); the ``ManifestSources`` fields are listed in ``tests/helpers/manifest_data.py``.

Builder rules pinned here (spec D-C, lead decisions 7, 8 and 11 and question 1):

* ``flags`` = sorted unique union of every event's ``data.flag`` (a string) and ``data.flags`` (a
  list of strings) and of ``annotation_meta.flags``, minus ``layer_missing:genomad`` when
  ``mobile`` is off and ``layer_missing:antismash`` when ``bgc`` is off (``cog``/``go``/``kegg``
  always kept). ``layers_missing`` = the ``<layer>`` of each kept ``layer_missing:<layer>`` flag.
* ``warnings`` = sorted unique ``preflight:<rule_id>`` for each ``warn`` check and
  ``<step>:<event>`` for each WARNING event (machine IDs, never messages).
* ``steps`` (Nextflow) = the trace grouped by step name (``giap.runner.nextflow``'s rule: process
  name lowercased, the tag for ``GIAP_STEP``), last attempt of each task only: ``wall_time_s``
  = sum of realtime, ``cpu_time_s`` = sum of realtime * %cpu / 100, ``peak_rss_bytes`` = max
  peak_rss (binary units, an int), ``n_tasks`` = distinct task names, ``n_retries`` = earlier
  attempts. A metric ``-`` or empty is left out. ``steps`` (local) = ``local_steps`` with
  ``n_tasks 1``, ``n_retries 0`` and no ``peak_rss_bytes``.
* ``seeds`` = ``{"master": cfg.seed}`` plus ``crc32(step)`` for every executed py step.
* ``skipped`` = registered minus executed, with ``reason`` ``module_off:<m>`` / ``opt_in_off`` /
  ``precomputed:<key>`` / ``not_run`` (that precedence), sorted by name.
"""

from __future__ import annotations

import json
import random
import zlib
from typing import Any

import pytest
from helpers.manifest_data import (
    SHA_LAYER,
    build,
    check,
    event,
    layer_missing,
    preflight_doc,
    spec,
    trace_row,
    validate,
)
from hypothesis import given, settings
from hypothesis import strategies as st

from giap.config import GiapConfig

KEYS_03A = {
    "schema_version", "giap_version", "git_commit", "git_dirty", "command_line", "executor",
    "profile", "resume", "outdir", "started_at", "finished_at", "status", "exit_code",
    "failed_step", "failed_log", "publication", "config", "inputs", "inputs_absent", "tools",
    "databases", "annotation", "detection_tier", "seeds", "steps", "skipped", "modules_off",
    "flags", "layers_missing", "warnings", "volatile", "provenance_digest",
}  # fmt: skip


# ----------------------------------------------------------------------------- shape


#: The tools/databases readers moved to MVP-03c (c06f140); 03c reverts it and drops this mark.


def test_builder_output_has_every_03a_field_and_validates() -> None:
    m = build()
    missing = KEYS_03A - set(m)
    assert not missing, f"manifest lacks {sorted(missing)}"
    validate(m)


def test_builder_output_is_strict_json() -> None:
    m = build()
    text = json.dumps(m, allow_nan=False, sort_keys=True)
    assert json.loads(text) == m


def test_completed_run_fields() -> None:
    m = build()
    assert m["status"] == "completed"
    assert m["exit_code"] == 0
    assert m["failed_step"] is None and m["failed_log"] is None
    assert m["executor"] == "local"
    assert m["schema_version"] == GiapConfig().schema_version
    assert m["config"] == GiapConfig(seed=7).model_dump(mode="json")
    assert m["detection_tier"] == "peptides"


def test_failed_run_fields() -> None:
    m = build(exit_code=4, failed_step="boom", failed_log="/runs/OUT/logs/boom.log")
    assert (m["status"], m["exit_code"]) == ("failed", 4)
    assert (m["failed_step"], m["failed_log"]) == ("boom", "/runs/OUT/logs/boom.log")
    validate(m)


# ----------------------------------------------------------------------------- example 1


def test_example_1_resume_duplicate_event_counted_once() -> None:
    ev = layer_missing("gene_table", "cog")
    m = build(events=[ev, dict(ev)])
    assert m["flags"] == ["layer_missing:cog"]
    assert m["layers_missing"] == ["cog"]
    assert m["warnings"] == ["gene_table:gene_table.layer_missing"]


# ----------------------------------------------------------------------------- example 2


def test_example_2_annotation_meta_flags_and_block() -> None:
    meta = {"db_type": "light", "flags": ["anti_sd_default", "bakta_db_light"],
            "bakta_version": "1.11.4", "db_version": "6.0", "annotation_source": "bakta_light",
            "anti_sd": "CCTCCT", "n_genes": 10, "n_features": 12, "skipped_types": {}}  # fmt: skip
    m = build(annotation_meta=meta)
    assert {"anti_sd_default", "bakta_db_light"} <= set(m["flags"])
    assert m["annotation"]["db_type"] == "light"
    assert m["annotation"]["bakta_version"] == "1.11.4"
    assert m["annotation"]["db_version"] == "6.0"
    assert m["annotation"]["source"] == "bakta_light"
    validate(m)


# ----------------------------------------------------------------------------- example 3


def test_example_3_off_module_explains_its_layer() -> None:
    cfg = GiapConfig(seed=7, modules=["core", "mobile"])
    m = build(cfg=cfg, events=[layer_missing("context_layers", "antismash")])
    assert m["modules_off"] == ["bgc", "function"]
    assert "antismash" not in m["layers_missing"]
    assert "layer_missing:antismash" not in m["flags"]


def test_decision_7_genomad_dropped_when_mobile_off() -> None:
    cfg = GiapConfig(seed=7, modules=["core"])
    m = build(cfg=cfg, events=[layer_missing("context_layers", "genomad")])
    assert "layer_missing:genomad" not in m["flags"]
    assert m["layers_missing"] == []


def test_decision_7_off_module_layer_kept_when_module_on() -> None:
    m = build(events=[layer_missing("context_layers", "antismash")])  # every module on
    assert "layer_missing:antismash" in m["flags"]
    assert m["layers_missing"] == ["antismash"]


@pytest.mark.parametrize("layer", ["cog", "go", "kegg"])
def test_decision_7_cog_go_kegg_always_kept(layer: str) -> None:
    cfg = GiapConfig(seed=7, modules=["core"])
    m = build(cfg=cfg, events=[layer_missing("gene_table", layer)])
    assert f"layer_missing:{layer}" in m["flags"]
    assert layer in m["layers_missing"]


def test_modules_off_sorted_for_core_only() -> None:
    m = build(cfg=GiapConfig(seed=7, modules=["core"]))
    assert m["modules_off"] == ["bgc", "function", "mobile"]


def test_modules_off_empty_with_every_module() -> None:
    assert build()["modules_off"] == []


# ----------------------------------------------------------------------------- example 4


def test_example_4_one_trace_row() -> None:
    row = trace_row("ANNOTATE:PARSE_BAKTA (parse_bakta)", realtime="5.5s", **{"%cpu": "65.8%"},
                    peak_rss="269.9 MB")  # fmt: skip
    m = build(executor="nextflow", trace=[row], local_steps=(),
              registered=[spec("parse_bakta")])  # fmt: skip
    assert len(m["steps"]) == 1
    (s,) = m["steps"]
    assert set(s) == {"name", "status", "n_tasks", "n_retries", "wall_time_s", "cpu_time_s",
                      "peak_rss_bytes"}  # fmt: skip
    assert (s["name"], s["status"], s["n_tasks"], s["n_retries"]) == ("parse_bakta", "completed",
                                                                      1, 0)  # fmt: skip
    assert s["wall_time_s"] == pytest.approx(5.5)
    assert s["cpu_time_s"] == pytest.approx(3.619)
    assert s["peak_rss_bytes"] == 283010662
    assert isinstance(s["peak_rss_bytes"], int)
    validate(m)


def test_trace_units_are_binary_and_durations_compound() -> None:
    row = trace_row("PROMOTECH (chunk_1)", realtime="1m 2.5s", **{"%cpu": "200.0%"},
                    peak_rss="1.5 GB")  # fmt: skip
    m = build(executor="nextflow", trace=[row], local_steps=(),
              registered=[spec("promotech", kind="ext")])  # fmt: skip
    (s,) = m["steps"]
    assert s["wall_time_s"] == pytest.approx(62.5)
    assert s["cpu_time_s"] == pytest.approx(125.0)
    assert s["peak_rss_bytes"] == int(1.5 * 1024**3)


@pytest.mark.parametrize("blank", ["-", ""])
def test_missing_trace_metric_is_left_out_never_zero(blank: str) -> None:
    row = trace_row("PROMOTECH (chunk_1)", realtime="5s", **{"%cpu": blank}, peak_rss=blank)
    m = build(executor="nextflow", trace=[row], local_steps=(),
              registered=[spec("promotech", kind="ext")])  # fmt: skip
    (s,) = m["steps"]
    assert s["wall_time_s"] == pytest.approx(5.0)
    assert "cpu_time_s" not in s
    assert "peak_rss_bytes" not in s
    validate(m)


# ----------------------------------------------------------------------------- example 5


def _example_5_trace() -> list[dict[str, str]]:
    return [
        trace_row("PROMOTECH (chunk_1)", "FAILED", exit="137", attempt="1", realtime="100s",
                  **{"%cpu": "100.0%"}, peak_rss="9 GB"),
        trace_row("PROMOTECH (chunk_2)", "COMPLETED", attempt="1", realtime="20s",
                  **{"%cpu": "50.0%"}, peak_rss="1 GB"),
        trace_row("PROMOTECH (chunk_1)", "COMPLETED", attempt="2", realtime="30s",
                  **{"%cpu": "100.0%"}, peak_rss="2 GB"),
    ]  # fmt: skip


def test_example_5_retry_counts_last_attempt_only() -> None:
    m = build(executor="nextflow", trace=_example_5_trace(), local_steps=(),
              registered=[spec("promotech", kind="ext")])  # fmt: skip
    (s,) = [s for s in m["steps"] if s["name"] == "promotech"]
    assert len(m["steps"]) == 1
    assert (s["n_tasks"], s["n_retries"], s["status"]) == (2, 1, "completed")
    assert s["wall_time_s"] == pytest.approx(50.0)  # 30 s (attempt 2) + 20 s (chunk_2)
    assert s["cpu_time_s"] == pytest.approx(30.0 + 10.0)
    assert s["peak_rss_bytes"] == 2 * 1024**3  # the failed 9 GB attempt is not counted


def test_failed_last_attempt_gives_status_failed() -> None:
    trace = [
        trace_row("PROMOTECH (chunk_1)", "FAILED", exit="137", attempt="1", realtime="1s"),
        trace_row("PROMOTECH (chunk_1)", "FAILED", exit="137", attempt="2", realtime="1s"),
        trace_row("ANNOTATE:GIAP_STEP (parse_bakta)", "COMPLETED", realtime="2s"),
    ]
    m = build(executor="nextflow", trace=trace, local_steps=(), exit_code=4,
              failed_step="promotech", failed_log="/w/ab/cdef12/.command.log",
              registered=[spec("promotech", kind="ext"), spec("parse_bakta")])  # fmt: skip
    by = {s["name"]: s for s in m["steps"]}
    assert by["promotech"]["status"] == "failed"
    assert by["promotech"]["n_retries"] == 1
    assert by["parse_bakta"]["status"] == "completed"
    validate(m)


def test_cached_rows_give_status_cached() -> None:
    trace = [trace_row("ANNOTATE:GIAP_STEP (parse_bakta)", "CACHED", realtime="2s")]
    m = build(executor="nextflow", trace=trace, local_steps=(), resume=True,
              registered=[spec("parse_bakta")])  # fmt: skip
    (s,) = m["steps"]
    assert s["status"] == "cached"
    assert s["wall_time_s"] == pytest.approx(2.0)  # the original metrics
    assert m["resume"] is True


def test_steps_keep_trace_order() -> None:
    trace = [trace_row("X:GIAP_STEP (zeta)", realtime="1s"),
             trace_row("X:GIAP_STEP (alpha)", realtime="1s")]  # fmt: skip
    m = build(executor="nextflow", trace=trace, local_steps=(),
              registered=[spec("zeta"), spec("alpha")])  # fmt: skip
    assert [s["name"] for s in m["steps"]] == ["zeta", "alpha"]


def test_trace_absent_gives_no_steps() -> None:
    m = build(executor="nextflow", trace=None, local_steps=(), exit_code=4)
    assert m["steps"] == []
    assert m["status"] == "failed"
    assert m["failed_step"] is None
    validate(m)


def test_local_steps_in_plan_order_without_rss() -> None:
    m = build()
    assert [s["name"] for s in m["steps"]] == ["detectability", "gene_table"]
    for s in m["steps"]:
        assert (s["n_tasks"], s["n_retries"], s["status"]) == (1, 0, "completed")
        assert "peak_rss_bytes" not in s
    assert m["steps"][0]["wall_time_s"] == pytest.approx(1.5)
    assert m["steps"][0]["cpu_time_s"] == pytest.approx(1.25)


# ----------------------------------------------------------------------------- example 6


def test_example_6_seeds() -> None:
    m = build(cfg=GiapConfig(seed=7))
    assert m["seeds"] == {"master": 7, "detectability": 3987542110, "gene_table": 393846822}


def test_seeds_are_crc32_of_step_names() -> None:
    assert zlib.crc32(b"detectability") == 3987542110
    assert zlib.crc32(b"gene_table") == 393846822


def test_seeds_skip_external_steps() -> None:
    trace = [trace_row("PROMOTECH (chunk_1)", realtime="1s"),
             trace_row("X:GIAP_STEP (gene_table)", realtime="1s")]  # fmt: skip
    m = build(executor="nextflow", trace=trace, local_steps=(), cfg=GiapConfig(seed=3),
              registered=[spec("promotech", kind="ext"), spec("gene_table")])  # fmt: skip
    assert m["seeds"] == {"master": 3, "gene_table": zlib.crc32(b"gene_table")}


# ----------------------------------------------------------------------------- skipped


def test_skipped_reasons_and_order() -> None:
    cfg = GiapConfig.model_validate(
        {"seed": 7, "modules": ["core"],
         "annotation": {"precomputed": {"terminators": "/data/terminators.tsv"}}}
    )  # fmt: skip
    registered = [
        spec("detectability"), spec("gene_table"),
        spec("genomad", kind="ext", module="mobile"),
        spec("bactermfinder", kind="ext", module="opt_in"),
        spec("terminatornet", kind="ext", precomputed="terminators"),
        spec("context_stats"),
    ]  # fmt: skip
    m = build(cfg=cfg, registered=registered)
    assert m["skipped"] == [
        {"name": "bactermfinder", "reason": "opt_in_off"},
        {"name": "context_stats", "reason": "not_run"},
        {"name": "genomad", "reason": "module_off:mobile"},
        {"name": "terminatornet", "reason": "precomputed:terminators"},
    ]
    assert "opt_in" not in m["modules_off"]
    validate(m)


def test_executed_steps_are_not_skipped() -> None:
    m = build()
    assert m["skipped"] == []


# ----------------------------------------------------------------------------- inputs


def test_inputs_reuse_preflight_hashes_and_add_paths() -> None:
    extra = [{"role": "precomputed:terminators", "file": "terminators.tsv",
              "path": "/data/terminators.tsv", "sha256": SHA_LAYER, "size": 42}]  # fmt: skip
    m = build(extra_inputs=extra)
    by = {r["role"]: r for r in m["inputs"]}
    assert by["genome"] == {"role": "genome", "file": "genome.fa", "path": "/data/genome.fa",
                            "sha256": "a" * 64, "size": 1200}  # fmt: skip
    assert by["peptides"]["path"] == "/data/peptides.tsv"
    assert by["precomputed:terminators"]["sha256"] == SHA_LAYER
    assert len(by) == len(m["inputs"]), "a role appears twice"
    validate(m)


def test_inputs_absent_sorted() -> None:
    m = build(inputs_absent=["replicon_depth", "abundance", "proteins"])
    assert m["inputs_absent"] == ["abundance", "proteins", "replicon_depth"]


# ----------------------------------------------------------------------------- flags / warnings


def test_data_flags_list_is_read_without_a_special_case() -> None:
    ev = event("function", "function.run_flags", {"flags": ["tu_recal_skipped",
                                                             "p_detect_clipped"]})  # fmt: skip
    m = build(events=[ev])
    assert {"tu_recal_skipped", "p_detect_clipped"} <= set(m["flags"])
    assert m["warnings"] == []  # INFO events give no warning


def test_run_flag_event_shape() -> None:
    ev = event("placement", "placement.run_flag", {"flag": "c_assumed"}, level="WARNING")
    m = build(events=[ev])
    assert m["flags"] == ["c_assumed"]
    assert m["warnings"] == ["placement:placement.run_flag"]


def test_events_without_flags_add_nothing() -> None:
    ev = event("preflight", "environment.resources", {"free_disk_bytes": 1, "ram_bytes": 2})
    m = build(events=[ev])
    assert m["flags"] == [] and m["warnings"] == []


def test_preflight_warn_checks_are_rule_ids() -> None:
    doc = preflight_doc(checks=[check("G03", "warn"), check("G01", "ok"), check("P02", "warn"),
                                check("P05", "info")])  # fmt: skip
    m = build(preflight=doc)
    assert m["warnings"] == ["preflight:G03", "preflight:P02"]


def test_flags_and_warnings_sorted_unique() -> None:
    evs = [event("b", "b.run_flag", {"flag": "zz_flag"}, level="WARNING"),
           event("a", "a.run_flag", {"flag": "aa_flag"}, level="WARNING"),
           event("b", "b.run_flag", {"flag": "zz_flag"}, level="WARNING")]  # fmt: skip
    m = build(events=evs)
    assert m["flags"] == ["aa_flag", "zz_flag"]
    assert m["warnings"] == ["a:a.run_flag", "b:b.run_flag"]


def _flag_events() -> list[dict[str, Any]]:
    return [
        layer_missing("gene_table", "cog"),
        layer_missing("context_layers", "genomad"),
        event("placement", "placement.run_flag", {"flag": "c_assumed"}, level="WARNING"),
        event("function", "function.run_flags", {"flags": ["null_degenerate"]}),
        event("mapping", "mapping.search_db_resolved", {"search_db": "same", "n_queries": 3}),
        event("preflight", "environment.resources", {"ram_bytes": 2}),
    ]


@settings(max_examples=60, deadline=None)
@given(st.randoms(use_true_random=False))
def test_property_builder_invariant_under_event_permutation(rnd: random.Random) -> None:
    evs = _flag_events()
    want = build(events=evs)
    shuffled = list(evs)
    rnd.shuffle(shuffled)
    assert build(events=shuffled) == want


# ----------------------------------------------------------------------------- databases


def _setup_manifest() -> dict[str, Any]:
    item = {"url": "https://example.invalid/x", "path": "x", "message": ""}
    return {"schema_version": "1.0", "items": [
        {**item, "name": "kegg", "version": "2026-09", "checksum": "sha256:" + "1" * 64,
         "size_bytes": 10, "retrieved": "2026-10-01", "status": "ok"},
        {**item, "name": "bakta_light", "version": "6.1", "checksum": "sha256:" + "2" * 64,
         "size_bytes": 20, "retrieved": "2026-09-01", "status": "ok"},
        {**item, "name": "genomad", "version": "1.9", "checksum": "", "size_bytes": 0,
         "retrieved": "2026-09-01", "status": "failed"},
    ]}  # fmt: skip


def test_databases_from_setup_manifest_ok_items_sorted() -> None:
    m = build(db_manifest=_setup_manifest())
    assert [d["name"] for d in m["databases"]] == ["bakta_light", "kegg"]
    kegg = m["databases"][1]
    assert kegg == {"name": "kegg", "version": "2026-09", "checksum": "sha256:" + "1" * 64,
                    "retrieved": "2026-10-01", "size_bytes": 10}  # fmt: skip
    validate(m)


def test_databases_from_provisional_form() -> None:
    m = build(db_manifest={"cog2020_def": {"version": "2020"}, "bakta_db": {"version": "6.0"}})
    assert [d["name"] for d in m["databases"]] == ["bakta_db", "cog2020_def"]
    for d in m["databases"]:
        assert d["checksum"] is None and d["retrieved"] is None
    assert m["databases"][0]["version"] == "6.0"
    validate(m)


def test_no_db_dir_gives_no_databases() -> None:
    assert build(db_manifest=None)["databases"] == []


# ----------------------------------------------------------------------------- tools


def _report(process: str, name: str, version: str, digest: str | None = None,
            lock: str | None = None) -> dict[str, Any]:  # fmt: skip
    return {"process": process, "name": name, "version": version,
            "container_digest": digest, "conda_lock_hash": lock}  # fmt: skip


def test_tools_one_row_per_name_and_version() -> None:
    reports = [
        _report("ANNOTATE:BAKTA", "bakta", "1.11.4", "sha256:" + "d" * 64),
        _report("NOVEL:BAKTA", "bakta", "1.11.4", "sha256:" + "d" * 64),
        _report("ANNOTATE:PROMOTECH", "promotech", "1.0", None, "e" * 64),
        _report("OTHER:PROMOTECH", "promotech", "1.1", None, "f" * 64),
        _report("ANNOTATE:GIAP_STEP", "giap", "0.1.0"),
    ]
    m = build(executor="nextflow", trace=[], local_steps=(), tool_reports=reports)
    assert [(t["name"], t["version"]) for t in m["tools"]] == [
        ("bakta", "1.11.4"), ("promotech", "1.0"), ("promotech", "1.1")]  # fmt: skip
    bakta = m["tools"][0]
    assert sorted(bakta["processes"]) == ["ANNOTATE:BAKTA", "NOVEL:BAKTA"]
    assert bakta["container_digest"] == "sha256:" + "d" * 64
    assert m["tools"][1]["conda_lock_hash"] == "e" * 64
    assert m["tools"][1]["container_digest"] is None
    validate(m)


def test_giap_itself_is_not_a_tool() -> None:
    m = build(executor="nextflow", trace=[], local_steps=(),
              tool_reports=[_report("X:GIAP_STEP", "giap", "0.1.0")])  # fmt: skip
    assert m["tools"] == []


def test_local_has_no_tools() -> None:
    assert build()["tools"] == []


# ----------------------------------------------------------------------------- git, volatile


def test_non_git_install() -> None:
    m = build(git_commit=None, git_dirty=None)
    assert m["git_commit"] is None and m["git_dirty"] is None
    validate(m)


def test_git_checkout_fields() -> None:
    m = build(git_commit="0" * 40, git_dirty=True)
    assert (m["git_commit"], m["git_dirty"]) == ("0" * 40, True)
    validate(m)


def test_volatile_field_lists_the_constant() -> None:
    from giap.manifest import VOLATILE

    assert build()["volatile"] == list(VOLATILE)


def test_digest_field_matches_the_core() -> None:
    from giap.manifest import provenance_core, provenance_digest

    m = build()
    assert m["provenance_digest"] == provenance_digest(provenance_core(m))


def test_publication_flag_recorded() -> None:
    assert build(publication=True)["publication"] is True
    assert build()["publication"] is False


def test_builder_does_not_mutate_its_sources() -> None:
    from helpers.manifest_data import sources

    from giap.manifest import build_manifest

    evs = _flag_events()
    before = json.dumps(evs, sort_keys=True)
    build_manifest(sources(events=evs))
    assert json.dumps(evs, sort_keys=True) == before
