"""Missing Bakta GBFF warning (M11-17; docs/specs/M11-17-preflight-gbff-warning.md).

With a precomputed Bakta JSON, Bakta does not run and antiSMASH (bgc module) reads the GBFF found
next to the JSON (``<stem>.gbff``, else ``<stem>.gbff.gz``; workflow/subworkflows/local/
annotate.nf). Without one, antiSMASH is skipped (``layer_missing:antismash``). This rule tells
the user in preflight. It never errors and never changes results.
"""

from __future__ import annotations

import re
from pathlib import Path

from gippa.config import GippaConfig
from gippa.preflight.checks import Check, JsonValue

__all__ = ["RULE_ID", "RULE_TEXT", "bgc_gbff_check", "expected_gbff"]

RULE_ID = "inputs.bakta_gbff"
RULE_TEXT = (
    "bgc module on and a precomputed Bakta JSON <stem>.json[.gz] -> <stem>.gbff or "
    "<stem>.gbff.gz must exist next to it, else warn (antiSMASH is skipped: layer_missing)"
)
MSG_BGC_OFF = "not evaluated: the bgc module is off"
MSG_REGIONS = "not evaluated: annotation.precomputed.regions replaces antiSMASH"

_JSON_SUFFIX = re.compile(r"(?i)\.json(\.gz)?$")


def _is_bakta_json(path: Path) -> bool:
    """Named like a Bakta JSON: ``.json`` or ``.json.gz``, any case (as ``gippa.cli``)."""
    return path.name.lower().endswith((".json", ".json.gz"))


def expected_gbff(bakta_json: Path) -> tuple[Path, Path]:
    """``(<stem>.gbff, <stem>.gbff.gz)`` for ``bakta_json`` with a final ``.json`` or
    ``.json.gz`` (any case) removed; pure, a relative path stays relative."""
    stem = _JSON_SUFFIX.sub("", str(bakta_json))
    return Path(stem + ".gbff"), Path(stem + ".gbff.gz")


def _bakta_json(cfg: GippaConfig, annotation: Path | None) -> Path | None:
    configured = cfg.annotation.precomputed.bakta_json
    if configured is not None:
        return configured
    if annotation is not None and _is_bakta_json(annotation):
        return annotation
    return None


def _message(js: Path, gbff: Path, gbff_gz: Path) -> str:
    return (
        f"no GBFF next to the precomputed Bakta JSON {js}: expected {gbff} (or {gbff_gz}); "
        "antiSMASH will be skipped (layer_missing:antismash) and the bgc module will report no "
        "biosynthetic gene clusters. Put Bakta's GBFF output there, or leave bgc out of "
        "--modules."
    )


def bgc_gbff_check(cfg: GippaConfig, annotation: Path | None) -> Check | None:
    """The ``inputs.bakta_gbff`` check, or None when the run has no Bakta JSON (Bakta runs and
    writes its own GBFF). The configured ``annotation.precomputed.bakta_json`` wins over a
    ``.json[.gz]`` ``annotation``; the level is ``ok``, ``info`` or ``warn``."""
    js = _bakta_json(cfg, annotation)
    if js is None:
        return None
    if "bgc" not in cfg.modules:
        return Check(RULE_ID, "inputs", "info", RULE_TEXT, None, MSG_BGC_OFF)
    if cfg.annotation.precomputed.regions is not None:
        return Check(RULE_ID, "inputs", "info", RULE_TEXT, None, MSG_REGIONS)
    gbff, gbff_gz = expected_gbff(js)
    found = next((p for p in (gbff, gbff_gz) if p.is_file()), None)
    value: dict[str, JsonValue] = {
        "bakta_json": str(js),
        "expected": [str(gbff), str(gbff_gz)],
        "gbff": None if found is None else str(found),
    }
    if found is not None:
        return Check(RULE_ID, "inputs", "ok", RULE_TEXT, value, "")
    return Check(RULE_ID, "inputs", "warn", RULE_TEXT, value, _message(js, gbff, gbff_gz))
