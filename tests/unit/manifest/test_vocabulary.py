"""MVP-03a run-flag vocabulary (docs/specs/MVP-03-manifest-finalize.md D-E, "Needs Khoa" 1, lead
decisions 2026-10-22; GIAP_SPEC.md §9.4 Run row, ADR-0036 §2, ADR-0037).

Interface assumed: ``giap.manifest.RUN_FLAGS`` and ``giap.manifest.RUN_FLAGS_PENDING`` are
disjoint ``frozenset[str]``. An entry ending in ``:`` is a prefix (``layer_missing:`` covers every
``layer_missing:<layer>``). Per the lead decisions the run-scope tokens beyond §9.4 —
``input_prefiltered``, ``cai_reference_missing``, ``no_detection_limit_model`` and the
``uncalibrated:`` prefix — sit in ``RUN_FLAGS_PENDING`` until Khoa answers.

"Every token emitted anywhere in ``src/giap``": the run-scope constants named in D-E (imported
here, so a rename is caught), the ``bakta_db_light`` and ``c_assumed`` literals, the function
module's run flags (ADR-0036 §2, ADR-0037) and every ``layer_missing:<word>`` literal found in
the source tree.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import pytest

SRC = Path(__file__).resolve().parents[3] / "src" / "giap"
FLAG_PATTERN = re.compile(r"^[a-z0-9_]+(:[a-z0-9_.-]+)?$")

#: (module, constant) of every run-scope token constant in src/giap (D-E).
CONSTANTS = [
    ("giap.annotate.bakta_layers", "ANTI_SD_DEFAULT_FLAG"),
    ("giap.annotate.operons", "DISTANCE_FLAG"),
    ("giap.mapping.placement_fdr", "PREFILTERED_FLAG"),
    ("giap.annotate.codon_usage", "CAI_REFERENCE_MISSING"),
    ("giap.operons.discordance", "NO_DL_MODEL"),
    ("giap.operons.flags", "LAYER_MISSING"),
]
LITERALS = ["bakta_db_light", "c_assumed"]
FUNCTION_FLAGS = ["tu_recal_skipped", "too_few_ranked_tus", "no_eligible_genes",
                  "p_detect_clipped", "null_degenerate", "tu_recal_reduced"]  # fmt: skip
SECTION_9_4 = ["layer_missing:", "anti_sd_default", "bakta_db_light", "c_assumed",
               "operon_source_distance_model"]  # fmt: skip
PENDING = ["input_prefiltered", "cai_reference_missing", "no_detection_limit_model",
           "uncalibrated:"]  # fmt: skip


def _vocab() -> tuple[frozenset[str], frozenset[str]]:
    from giap.manifest import RUN_FLAGS, RUN_FLAGS_PENDING

    return frozenset(RUN_FLAGS), frozenset(RUN_FLAGS_PENDING)


def _covered(token: str, vocab: frozenset[str]) -> bool:
    if token in vocab:
        return True
    return any(v.endswith(":") and token.startswith(v) and len(token) > len(v) for v in vocab)


def test_vocabularies_are_disjoint_frozensets() -> None:
    from giap.manifest import RUN_FLAGS, RUN_FLAGS_PENDING

    assert isinstance(RUN_FLAGS, frozenset) and isinstance(RUN_FLAGS_PENDING, frozenset)
    assert not RUN_FLAGS & RUN_FLAGS_PENDING


@pytest.mark.parametrize("token", SECTION_9_4 + FUNCTION_FLAGS)
def test_registered_tokens(token: str) -> None:
    flags, _ = _vocab()
    assert token in flags


@pytest.mark.parametrize("token", PENDING)
def test_pending_tokens_wait_for_khoa(token: str) -> None:
    flags, pending = _vocab()
    assert token in pending
    assert token not in flags


@pytest.mark.parametrize(("module", "name"), CONSTANTS)
def test_every_run_flag_constant_is_known(module: str, name: str) -> None:
    flags, pending = _vocab()
    token = getattr(importlib.import_module(module), name)
    assert isinstance(token, str)
    assert token in flags | pending, f"{module}.{name} = {token!r} is in neither vocabulary"


@pytest.mark.parametrize("token", LITERALS + FUNCTION_FLAGS)
def test_literal_tokens_are_known(token: str) -> None:
    flags, pending = _vocab()
    assert _covered(token, flags | pending)


def test_every_layer_missing_literal_in_src_is_covered() -> None:
    flags, pending = _vocab()
    found = set()
    for path in SRC.rglob("*.py"):
        found |= set(re.findall(r"layer_missing:[a-z0-9_]+", path.read_text(encoding="utf-8")))
    assert found, "no layer_missing:<layer> literal found under src/giap"
    for token in sorted(found):
        assert _covered(token, flags | pending), token


def test_tokens_match_the_schema_pattern() -> None:
    flags, pending = _vocab()
    for v in sorted(flags | pending):
        sample = v + "x" if v.endswith(":") else v
        assert FLAG_PATTERN.fullmatch(sample), v


def test_schema_flag_pattern_is_the_spec_one() -> None:
    from helpers.manifest_data import schema

    items = schema()["properties"]["flags"]["items"]
    assert items.get("pattern") == FLAG_PATTERN.pattern
