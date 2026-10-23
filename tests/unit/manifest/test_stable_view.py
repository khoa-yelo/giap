"""MVP-03a ``VOLATILE``, ``stable_view``, the provenance core and ``provenance_digest``
(docs/specs/MVP-03-manifest-finalize.md D-C "Volatile" column, D-D, lead decision 6; GIAP_SPEC.md
§6.4, §6.7).

Interface assumed (test-author names; the spec names ``VOLATILE``, ``stable_view`` and
``provenance_digest``):

* ``giap.manifest.VOLATILE: tuple[str, ...]`` holds RFC 6901 JSON pointers; ``*`` matches every
  array index or object key at that level. Exactly the D-C "Volatile" rows: ``/command_line``,
  ``/executor``, ``/profile``, ``/resume``, ``/outdir``, ``/started_at``, ``/finished_at``,
  ``/failed_log``, ``/steps``, ``/inputs/*/path``, ``/config/annotation/precomputed/*``.
* ``stable_view(manifest, volatile=VOLATILE) -> dict`` returns a new dict without those pointers
  and never mutates its argument. A precomputed config path may either be removed or reduced to
  its file name (D-D "with the same path normalisation"); both readings are accepted.
* ``canonical_json(obj) -> str`` = ``json.dumps(obj, sort_keys=True, indent=2,
  ensure_ascii=False, allow_nan=False) + "\\n"`` (D-C "Canonical JSON").
* ``provenance_core(manifest) -> dict`` = ``schema_version``, ``giap_version``, ``git_commit``,
  ``git_dirty``, ``publication``, ``config`` (each non-null ``annotation.precomputed.*`` value
  replaced by its file name), ``inputs`` (each row without ``path``), ``databases`` and
  ``seeds: {"master": ...}``.
* ``provenance_digest(core) -> str`` = ``sha256(canonical_json(core).encode()).hexdigest()``.
"""

from __future__ import annotations

import copy
import hashlib
import json
import math
from pathlib import PurePosixPath
from typing import Any

import pytest
from helpers.manifest_data import SHA_LAYER, build, trace_row

from giap.config import GiapConfig

WANT_VOLATILE = {
    "/command_line", "/executor", "/profile", "/resume", "/outdir", "/started_at",
    "/finished_at", "/failed_log", "/steps", "/inputs/*/path",
    "/config/annotation/precomputed/*",
}  # fmt: skip
CORE_KEYS = {"schema_version", "giap_version", "git_commit", "git_dirty", "publication", "config",
             "inputs", "databases", "seeds"}  # fmt: skip


def _cfg_with(where: str = "/data/in", **top: Any) -> GiapConfig:
    pre = {"bakta_json": f"{where}/bakta.json", "terminators": f"{where}/term.tsv"}
    return GiapConfig.model_validate({"seed": 7, **top, "annotation": {"precomputed": pre}})


def _cfg() -> GiapConfig:
    return _cfg_with()


def _manifest(**over: Any) -> dict[str, Any]:
    extra = [{"role": "precomputed:bakta_json", "file": "bakta.json",
              "path": "/data/in/bakta.json", "sha256": SHA_LAYER, "size": 9}]  # fmt: skip
    base: dict[str, Any] = {"cfg": _cfg(), "extra_inputs": extra}
    base.update(over)
    return build(**base)


def _delete(doc: Any, parts: list[str]) -> None:
    """Test oracle: delete RFC 6901 pointer ``parts`` (``*`` = every child) from ``doc``."""
    if not parts:
        return
    head, rest = parts[0], parts[1:]
    if isinstance(doc, dict):
        keys = list(doc) if head == "*" else [head] if head in doc else []
        for k in keys:
            if rest:
                _delete(doc[k], rest)
            else:
                del doc[k]
    elif isinstance(doc, list):
        idx = range(len(doc)) if head == "*" else [int(head)] if head.isdigit() else []
        if rest:
            for i in idx:
                _delete(doc[i], rest)
        else:
            for i in sorted(idx, reverse=True):
                del doc[i]


def _oracle_view(m: dict[str, Any], volatile: list[str]) -> dict[str, Any]:
    out = copy.deepcopy(m)
    for ptr in volatile:
        _delete(out, ptr.split("/")[1:])
    return out


def _without_precomputed(view: dict[str, Any]) -> dict[str, Any]:
    v = copy.deepcopy(view)
    v.get("config", {}).get("annotation", {}).pop("precomputed", None)
    return v


# ----------------------------------------------------------------------------- VOLATILE


def test_volatile_is_exactly_the_spec_list() -> None:
    from giap.manifest import VOLATILE

    assert isinstance(VOLATILE, tuple)
    assert len(VOLATILE) == len(set(VOLATILE)), "duplicate pointer"
    assert set(VOLATILE) == WANT_VOLATILE


def test_every_volatile_pointer_names_a_manifest_field() -> None:
    from giap.manifest import VOLATILE

    m = _manifest()
    for ptr in VOLATILE:
        assert ptr.split("/")[1] in m, ptr


# ----------------------------------------------------------------------------- stable_view


def test_stable_view_removes_the_volatile_fields() -> None:
    from giap.manifest import stable_view

    m = _manifest()
    view = stable_view(m)
    for top in ("command_line", "executor", "profile", "resume", "outdir", "started_at",
                "finished_at", "failed_log", "steps"):  # fmt: skip
        assert top not in view, top
    for row in view["inputs"]:
        assert "path" not in row and "sha256" in row and "role" in row


def test_stable_view_removes_exactly_volatile() -> None:
    from giap.manifest import VOLATILE, stable_view

    m = _manifest()
    view = stable_view(m)
    want = _oracle_view(m, [p for p in VOLATILE if not p.startswith("/config/")])
    assert _without_precomputed(view) == _without_precomputed(want)


def test_stable_view_precomputed_paths_removed_or_reduced_to_names() -> None:
    from giap.manifest import stable_view

    view = stable_view(_manifest())
    pre = view["config"]["annotation"].get("precomputed", {})
    for key, value in pre.items():
        assert value is None or value == PurePosixPath(str(value)).name, (key, value)
        assert value is None or "/" not in str(value), (key, value)


def test_stable_view_does_not_mutate() -> None:
    from giap.manifest import stable_view

    m = _manifest()
    before = json.dumps(m, sort_keys=True)
    stable_view(m)
    assert json.dumps(m, sort_keys=True) == before


def test_stable_view_equal_across_volatile_differences() -> None:
    from giap.manifest import stable_view

    a = _manifest()
    b = _manifest(executor="nextflow", profile="farmshare_local", resume=True,
                  outdir="/elsewhere/OUT", command_line=["giap", "run", "-x"],
                  started_at="2027-01-01T00:00:00+00:00",
                  finished_at="2027-01-01T01:00:00+00:00",
                  input_paths={"genome": "/x/genome.fa", "peptides": "/x/peptides.tsv"},
                  cfg=_cfg_with(where="/other"),
                  extra_inputs=[{"role": "precomputed:bakta_json", "file": "bakta.json",
                                 "path": "/other/bakta.json", "sha256": SHA_LAYER, "size": 9}],
                  trace=[trace_row("X:GIAP_STEP (detectability)", realtime="3s"),
                         trace_row("X:GIAP_STEP (gene_table)", realtime="4s")],
                  local_steps=())  # fmt: skip
    assert stable_view(a) == stable_view(b)


def test_stable_view_differs_on_a_stable_field() -> None:
    from giap.manifest import stable_view

    a = _manifest()
    b = _manifest(cfg=_cfg_with(seed=8))
    assert stable_view(a) != stable_view(b)


def test_stable_view_accepts_a_custom_list() -> None:
    from giap.manifest import VOLATILE, stable_view

    m = _manifest()
    fewer = [p for p in VOLATILE if p != "/executor"]
    assert stable_view(m, volatile=fewer)["executor"] == "local"
    assert "executor" not in stable_view(m, volatile=VOLATILE)


# ----------------------------------------------------------------------------- canonical JSON


def test_canonical_json_form() -> None:
    from giap.manifest import canonical_json

    obj = {"b": [1, {"z": "é", "a": None}], "a": True}
    text = canonical_json(obj)
    assert text == json.dumps(obj, sort_keys=True, indent=2, ensure_ascii=False,
                              allow_nan=False) + "\n"  # fmt: skip
    assert text.endswith("}\n") and not text.endswith("\n\n")


@pytest.mark.parametrize("bad", [math.nan, math.inf, -math.inf])
def test_canonical_json_rejects_non_finite(bad: float) -> None:
    from giap.manifest import canonical_json

    with pytest.raises(ValueError):
        canonical_json({"x": bad})


def test_round_trip_of_a_built_manifest() -> None:
    from giap.manifest import canonical_json

    m = _manifest()
    text = canonical_json(m)
    assert json.loads(text) == m
    assert "NaN" not in text and "Infinity" not in text


# ----------------------------------------------------------------------------- core and digest


def test_provenance_core_fields() -> None:
    from giap.manifest import provenance_core

    m = _manifest()
    core = provenance_core(m)
    assert set(core) == CORE_KEYS
    assert core["seeds"] == {"master": 7}
    for row in core["inputs"]:
        assert "path" not in row
    assert {r["role"] for r in core["inputs"]} == {r["role"] for r in m["inputs"]}
    pre = core["config"]["annotation"]["precomputed"]
    assert pre["bakta_json"] == "bakta.json"
    assert pre["terminators"] == "term.tsv"
    assert pre["promoters"] is None
    assert core["databases"] == m["databases"]
    for key in ("schema_version", "giap_version", "git_commit", "git_dirty", "publication"):
        assert core[key] == m[key], key


def test_provenance_core_does_not_mutate() -> None:
    from giap.manifest import provenance_core

    m = _manifest()
    before = json.dumps(m, sort_keys=True)
    provenance_core(m)
    assert json.dumps(m, sort_keys=True) == before


def test_digest_is_sha256_of_canonical_core() -> None:
    from giap.manifest import canonical_json, provenance_core, provenance_digest

    core = provenance_core(_manifest())
    want = hashlib.sha256(canonical_json(core).encode("utf-8")).hexdigest()
    assert provenance_digest(core) == want


def test_digest_stable_under_key_order() -> None:
    from giap.manifest import provenance_core, provenance_digest

    core = provenance_core(_manifest())

    def reverse(x: Any) -> Any:
        if isinstance(x, dict):
            return {k: reverse(x[k]) for k in reversed(list(x))}
        if isinstance(x, list):
            return [reverse(v) for v in x]
        return x

    assert provenance_digest(reverse(core)) == provenance_digest(core)


def test_digest_changes_with_one_input_sha256() -> None:
    from giap.manifest import provenance_core, provenance_digest

    core = provenance_core(_manifest())
    other = copy.deepcopy(core)
    row = next(r for r in other["inputs"] if r["role"] == "genome")
    row["sha256"] = "f" * 64
    assert provenance_digest(other) != provenance_digest(core)


def test_digest_ignores_tools_steps_flags_and_volatile() -> None:
    a = _manifest()
    b = _manifest(executor="nextflow", local_steps=(), outdir="/elsewhere",
                  trace=[{"name": "X:GIAP_STEP (detectability)", "status": "COMPLETED",
                          "attempt": "1", "realtime": "9s", "%cpu": "-", "peak_rss": "-"}],
                  events=[{"ts": "t", "step": "gene_table", "level": "WARNING",
                           "event": "gene_table.layer_missing",
                           "data": {"layer": "cog", "flag": "layer_missing:cog"}}],
                  tool_reports=[{"process": "A:BAKTA", "name": "bakta", "version": "1.11.4",
                                 "container_digest": None, "conda_lock_hash": None}],
                  started_at="2030-01-01T00:00:00+00:00")  # fmt: skip
    assert a["provenance_digest"] == b["provenance_digest"]


def test_digest_changes_with_the_config() -> None:
    a = _manifest()
    b = _manifest(cfg=_cfg_with(threads=3))
    assert a["provenance_digest"] != b["provenance_digest"]


def test_digest_independent_of_precomputed_directory() -> None:
    a = _manifest()
    b = _manifest(cfg=_cfg_with(where="/moved"))
    assert a["provenance_digest"] == b["provenance_digest"]
