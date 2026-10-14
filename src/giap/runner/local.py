"""Developer-only local executor (GIAP_SPEC.md §6.4; M2-02).

``giap run --executor local`` runs the Python steps in topological order in one process. External
steps never run here: each must be replaced by its ``annotation.precomputed`` layer, and a missing
layer is an error naming the key, raised before any step runs. Every step key is stored at
``workdir / key``; a consumer of an output of a replaced step receives the precomputed file itself.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

from giap.config import GiapConfig
from giap.steps import registry as reg


class LayerMissing(Exception):
    """An external step has no ``annotation.precomputed`` layer under the local executor."""

    def __init__(self, key: str) -> None:
        super().__init__(f"{key} is required with --executor local")
        self.key = key


def _layer(cfg: GiapConfig, key: str) -> Path | None:
    path: Path | None = getattr(cfg.annotation.precomputed, key)
    return path


def _enabled(cfg: GiapConfig, spec: reg.StepSpec, targets: Sequence[str]) -> bool:
    if spec.module == "opt_in":
        return spec.name in targets
    return spec.module in cfg.modules


def _walk(cfg: GiapConfig, targets: Sequence[str]) -> tuple[list[str], dict[str, Path]]:
    """Planned py steps and the step keys supplied by precomputed layers."""
    specs = {s.name: s for s in reg.all_steps()}
    producers: dict[str, list[str]] = {}
    for s in specs.values():
        for key in s.outputs:
            producers.setdefault(key, []).append(s.name)
    planned: set[str] = set()
    supplied: dict[str, Path] = {}
    seen: set[str] = set()
    stack = sorted(set(targets), reverse=True)
    while stack:
        name = stack.pop()
        if name in seen:
            continue
        seen.add(name)
        spec = reg.get(name)
        if not _enabled(cfg, spec, targets):
            continue
        layer = _layer(cfg, spec.precomputed) if spec.precomputed else None
        if layer is not None:
            for key in spec.outputs:
                supplied[key] = layer
            continue
        if spec.kind == "ext":
            if spec.precomputed is None:
                raise reg.WiringError(f"external step '{name}' cannot run with --executor local")
            raise LayerMissing(f"annotation.precomputed.{spec.precomputed}")
        planned.add(name)
        for key in sorted({*spec.inputs, *spec.optional_inputs}, reverse=True):
            stack.extend(sorted(producers.get(key, ()), reverse=True))
    return reg.topo_order(planned), supplied


def plan_local(cfg: GiapConfig, targets: Sequence[str]) -> list[str]:
    """Python steps needed for ``targets``, in run order; raises ``LayerMissing``."""
    return _walk(cfg, targets)[0]


def run_local(cfg: GiapConfig, workdir: Path, targets: Sequence[str]) -> None:
    """Plan (validating every layer first), then run each Python step under ``workdir``."""
    order, supplied = _walk(cfg, targets)
    workdir.mkdir(parents=True, exist_ok=True)
    for name in order:
        spec = reg.get(name)
        inputs: dict[str, Path] = {}
        for key in (*spec.inputs, *spec.optional_inputs):
            path = workdir / key
            if not path.exists() and key in supplied:
                path = supplied[key]
            if key in spec.inputs or path.exists():
                inputs[key] = path
        outputs = {key: workdir / key for key in spec.outputs}
        reg.run_step(name, cfg, inputs, outputs)
