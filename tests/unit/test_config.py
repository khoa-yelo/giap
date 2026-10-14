"""Configuration model, defaults sync and resolution (GIAP_SPEC.md §5.2)."""

from __future__ import annotations

import math
from pathlib import Path

import pytest
import yaml

DEFAULTS_YAML = Path(__file__).resolve().parents[2] / "src" / "giap" / "config" / "defaults.yaml"


def _raw_defaults() -> dict[str, object]:
    with DEFAULTS_YAML.open() as fh:
        data: dict[str, object] = yaml.safe_load(fh)
    return data


def test_load_defaults_returns_parsed_yaml() -> None:
    from giap.config import load_defaults

    assert load_defaults() == _raw_defaults()


def test_model_defaults_equal_defaults_yaml() -> None:
    from giap.config import GiapConfig, load_defaults

    assert GiapConfig().model_dump(mode="python") == load_defaults()


def test_distance_bins_keep_infinities() -> None:
    from giap.config import GiapConfig

    bins = GiapConfig().model_dump(mode="python")["operons"]["distance_bins_bp"]
    assert bins[0] == -math.inf
    assert bins[-1] == math.inf


def test_tolerant_max_diff_keys_are_length_ranges() -> None:
    from giap.config import GiapConfig

    dump = GiapConfig().model_dump(mode="python")
    assert dump["mapping"]["tolerant"]["max_diff"] == {"10-11": 1, "12-40": 2}


def test_unknown_top_level_key_rejected() -> None:
    from giap.config import GiapConfig

    with pytest.raises(ValueError):
        GiapConfig.model_validate({"no_such_key": 1})


def test_unknown_nested_key_rejected() -> None:
    from giap.config import GiapConfig

    with pytest.raises(ValueError):
        GiapConfig.model_validate({"stats": {"permutationz": 10}})


def test_unknown_deeply_nested_key_rejected() -> None:
    from giap.config import GiapConfig

    with pytest.raises(ValueError):
        GiapConfig.model_validate({"stats": {"alpha": {"nonsense": 0.1}}})


def test_resolve_without_arguments_gives_defaults() -> None:
    from giap.config import load_defaults, resolve_config

    assert resolve_config().model_dump(mode="python") == load_defaults()


def test_resolve_config_file_deep_merges(tmp_path: Path) -> None:
    from giap.config import resolve_config

    cfg = tmp_path / "c.yaml"
    cfg.write_text("seed: 7\nstats:\n  alpha:\n    context: 0.01\n")
    c = resolve_config(config_file=cfg)
    assert c.seed == 7
    assert c.stats.alpha.context == 0.01
    # sibling keys keep their defaults
    assert c.stats.alpha.enrichment == 0.05
    assert c.stats.permutations == 1000


def test_resolve_override_integer() -> None:
    from giap.config import resolve_config

    assert resolve_config(overrides=["stats.permutations=50"]).stats.permutations == 50


def test_resolve_override_yaml_list() -> None:
    from giap.config import resolve_config

    assert list(resolve_config(overrides=["modules=[core, bgc]"]).modules) == ["core", "bgc"]


def test_override_wins_over_config_file(tmp_path: Path) -> None:
    from giap.config import resolve_config

    cfg = tmp_path / "c.yaml"
    cfg.write_text("seed: 7\n")
    assert resolve_config(config_file=cfg, overrides=["seed=9"]).seed == 9


def test_resolve_override_unknown_key_raises() -> None:
    from giap.config import resolve_config

    with pytest.raises(ValueError, match="nonexistent"):
        resolve_config(overrides=["stats.nonexistent=1"])


def test_resolve_config_file_unknown_key_raises(tmp_path: Path) -> None:
    from giap.config import resolve_config

    cfg = tmp_path / "c.yaml"
    cfg.write_text("bogus_section: {a: 1}\n")
    with pytest.raises(ValueError):
        resolve_config(config_file=cfg)


def test_resolve_override_without_equals_raises() -> None:
    from giap.config import resolve_config

    with pytest.raises(ValueError):
        resolve_config(overrides=["seed"])


def test_publication_sets_permutations_and_bootstrap() -> None:
    from giap.config import resolve_config

    c = resolve_config(publication=True)
    assert c.stats.permutations == 10000
    assert c.stats.bootstrap == 10000


def test_genome_relation_rejects_other_values() -> None:
    from giap.config import resolve_config

    with pytest.raises(ValueError):
        resolve_config(overrides=["genome_relation=distant"])


def test_genome_relation_accepts_same() -> None:
    from giap.config import resolve_config

    assert resolve_config(overrides=["genome_relation=same"]).genome_relation == "same"


def test_abundance_missing_rejects_other_values() -> None:
    from giap.config import resolve_config

    with pytest.raises(ValueError):
        resolve_config(overrides=["abundance.missing=knn"])


@pytest.mark.parametrize("value", ["median", "exclude", "min", "zero", "mean"])
def test_abundance_missing_accepts_allowed_values(value: str) -> None:
    from giap.config import resolve_config

    assert resolve_config(overrides=[f"abundance.missing={value}"]).abundance.missing == value


def test_core_module_added_when_missing() -> None:
    from giap.config import resolve_config

    modules = list(resolve_config(overrides=["modules=[bgc]"]).modules)
    assert "core" in modules
    assert "bgc" in modules


def test_json_schema_lists_seed() -> None:
    from giap.config import json_schema

    schema = json_schema()
    assert isinstance(schema, dict)
    assert "seed" in schema["properties"]
