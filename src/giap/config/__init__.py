"""Configuration: pydantic model, bundled defaults and resolution (GIAP_SPEC.md §5.2)."""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from importlib.resources import files
from pathlib import Path
from typing import Any

import yaml

from giap.config.model import GiapConfig

__all__ = ["GiapConfig", "json_schema", "load_defaults", "resolve_config", "resolve_sample_name"]

PUBLICATION_RESAMPLES = 10000


def _load_yaml_mapping(text: str, source: str) -> dict[str, Any]:
    try:
        data = yaml.safe_load(text)
    except yaml.YAMLError as exc:
        raise ValueError(f"{source}: malformed YAML: {exc}") from exc
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ValueError(f"{source}: top level must be a mapping")
    return data


def load_defaults() -> dict[str, Any]:
    """Parse the bundled ``defaults.yaml``."""
    text = files("giap.config").joinpath("defaults.yaml").read_text(encoding="utf-8")
    return _load_yaml_mapping(text, "defaults.yaml")


def _deep_merge(base: dict[str, Any], update: Mapping[str, Any]) -> dict[str, Any]:
    out = dict(base)
    for key, value in update.items():
        if isinstance(value, Mapping) and isinstance(out.get(key), dict):
            out[key] = _deep_merge(out[key], value)
        else:
            out[key] = value
    return out


def _set_dotted(data: dict[str, Any], override: str) -> None:
    key, sep, raw = override.partition("=")
    if not sep or not key.strip():
        raise ValueError(f"--set expects key=value, got {override!r}")
    *parents, leaf = key.strip().split(".")
    node = data
    for part in parents:
        child = node.setdefault(part, {})
        if not isinstance(child, dict):
            raise ValueError(f"--set {key}: {part!r} is not a section")
        node = child
    try:
        node[leaf] = yaml.safe_load(raw)
    except yaml.YAMLError as exc:
        raise ValueError(f"--set {key}: malformed YAML value {raw!r}: {exc}") from exc


def resolve_config(
    config_file: Path | str | None = None,
    overrides: Iterable[str] = (),
    publication: bool = False,
) -> GiapConfig:
    """Defaults ← config file (deep merge) ← ``--publication`` ← ``--set key=value`` overrides.

    Raises ``ValueError`` (pydantic ``ValidationError`` included) naming any unknown key.
    """
    data = load_defaults()
    if config_file is not None:
        path = Path(config_file)
        data = _deep_merge(data, _load_yaml_mapping(path.read_text(encoding="utf-8"), str(path)))
    if publication:
        data = _deep_merge(
            data,
            {"stats": {"permutations": PUBLICATION_RESAMPLES, "bootstrap": PUBLICATION_RESAMPLES}},
        )
    for override in overrides:
        _set_dotted(data, override)
    return GiapConfig.model_validate(data)


def json_schema() -> dict[str, Any]:
    """JSON Schema of :class:`GiapConfig` (written to ``giap/schemas/config.schema.json``)."""
    return GiapConfig.model_json_schema()


_VCF_UNSAFE = re.compile(r"[^A-Za-z0-9._-]")


def resolve_sample_name(cfg: GiapConfig, peptides: Path | None, proteins: Path | None) -> str:
    """``cfg.sample_name`` with ``auto`` resolved (§5.2, M7-04 rule 4): the peptide file's name
    without ``.gz`` and its last suffix, else the protein file's, else ``sample``; then every
    character outside ``[A-Za-z0-9._-]`` becomes ``_`` (VCF-safe, as ADR-0005 maps contig IDs).
    Pure."""
    name = cfg.sample_name
    if name == "auto":
        src = peptides if peptides is not None else proteins
        base = "" if src is None else src.name.removesuffix(".gz")
        name = Path(base).stem or "sample"
    return _VCF_UNSAFE.sub("_", name)
