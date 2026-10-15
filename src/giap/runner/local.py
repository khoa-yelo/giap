"""Developer-only local executor (GIAP_SPEC.md §6.4; M2-02).

``giap run --executor local`` runs the Python steps in topological order in one process. External
steps never run here: each must be replaced by its ``annotation.precomputed`` layer, and a missing
layer is an error naming the key, raised before any step runs. Every step key is stored at
``workdir / key``; a consumer of an output of a replaced step receives the precomputed file itself,
even when an older ``workdir / key`` exists. Opt-in steps run only when their configuration enables
them (e.g. ``annotation.terminators.bactermfinder``).
"""

from __future__ import annotations

import shutil
from collections.abc import Callable, Collection, Mapping, Sequence
from pathlib import Path

from giap.config import GiapConfig
from giap.steps import registry as reg


class LayerMissing(Exception):
    """An external step has no ``annotation.precomputed`` layer under the local executor."""

    def __init__(self, key: str) -> None:
        super().__init__(f"{key} is required with --executor local")
        self.key = key


class ToolUnavailable(reg.WiringError):
    """An external step without a precomputed layer is needed (``giap run`` exit code 4)."""

    def __init__(self, name: str) -> None:
        super().__init__(
            f"step '{name}' runs the external tool {name}, which --executor local cannot run "
            "(GIAP_SPEC.md §6.4); use --executor nextflow"
        )
        self.name = name


def _layer(cfg: GiapConfig, key: str) -> Path | None:
    path: Path | None = getattr(cfg.annotation.precomputed, key)
    return path


# Opt-in steps follow configuration, not the target list; an unknown opt-in step stays off.
_OPT_IN: dict[str, Callable[[GiapConfig], bool]] = {
    "bactermfinder": lambda cfg: cfg.annotation.terminators.bactermfinder,
}


def _enabled(cfg: GiapConfig, spec: reg.StepSpec) -> bool:
    if spec.module == "opt_in":
        flag = _OPT_IN.get(spec.name)
        return flag is not None and bool(flag(cfg))
    return spec.module in cfg.modules


def _feasible(
    cfg: GiapConfig,
    specs: Mapping[str, reg.StepSpec],
    producers: Mapping[str, list[str]],
    roots: Collection[str],
) -> Callable[[str], bool]:
    """Whether a step can run: enabled, replaced by a layer, or every required input is a root
    key or has a feasible producer (e.g. the protein route needs ``proteins_long``)."""
    memo: dict[str, bool] = {}

    def ok(name: str, stack: frozenset[str] = frozenset()) -> bool:
        if name in memo:
            return memo[name]
        if name in stack:
            return False
        spec = specs[name]
        if not _enabled(cfg, spec):
            res = False
        elif spec.precomputed and _layer(cfg, spec.precomputed) is not None:
            res = True
        else:
            res = all(
                key in roots or any(ok(p, stack | {name}) for p in producers.get(key, ()))
                for key in spec.inputs
            )
        memo[name] = res
        return res

    return ok


def _walk(
    cfg: GiapConfig, targets: Sequence[str], roots: Collection[str] | None = None
) -> tuple[list[str], dict[str, Path]]:
    """Planned py steps and the step keys supplied by precomputed layers.

    With ``roots`` (the staged root keys), steps whose required inputs can be neither staged nor
    produced are left out instead of failing later.
    """
    specs = {s.name: s for s in reg.all_steps()}
    producers: dict[str, list[str]] = {}
    for s in specs.values():
        for key in s.outputs:
            producers.setdefault(key, []).append(s.name)
    feasible = None if roots is None else _feasible(cfg, specs, producers, roots)
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
        if not _enabled(cfg, spec) or (feasible is not None and not feasible(name)):
            continue
        layer = _layer(cfg, spec.precomputed) if spec.precomputed else None
        if layer is not None:
            for key in spec.outputs:
                supplied[key] = layer
            continue
        if spec.kind == "ext":
            if spec.precomputed is None:
                raise ToolUnavailable(name)
            raise LayerMissing(f"annotation.precomputed.{spec.precomputed}")
        planned.add(name)
        for key in sorted({*spec.inputs, *spec.optional_inputs}, reverse=True):
            stack.extend(sorted(producers.get(key, ()), reverse=True))
    return reg.topo_order(planned), supplied


def plan_local(
    cfg: GiapConfig, targets: Sequence[str], roots: Collection[str] | None = None
) -> list[str]:
    """Python steps needed for ``targets``, in run order; raises ``LayerMissing`` or
    ``ToolUnavailable``."""
    return _walk(cfg, targets, roots)[0]


def stage_roots(workdir: Path, roots: Mapping[str, Path]) -> None:
    """Copy root input files (normalised by preflight) to ``workdir / key``."""
    workdir.mkdir(parents=True, exist_ok=True)
    for key, src in sorted(roots.items()):
        shutil.copyfile(src, workdir / key)


def run_local(
    cfg: GiapConfig,
    workdir: Path,
    targets: Sequence[str],
    roots: Mapping[str, Path] | None = None,
) -> None:
    """Plan (validating every layer first), stage ``roots``, then run each Python step under
    ``workdir``. Without ``roots`` every planned step must find its inputs in ``workdir``."""
    order, supplied = _walk(cfg, targets, None if roots is None else set(roots))
    if roots is not None:
        stage_roots(workdir, roots)
    workdir.mkdir(parents=True, exist_ok=True)
    for name in order:
        spec = reg.get(name)
        inputs: dict[str, reg.InputPath] = {}
        for key in (*spec.inputs, *spec.optional_inputs):
            # a configured layer always wins over a stale ``workdir / key`` from an earlier run
            path = supplied.get(key, workdir / key)
            if key in spec.inputs or path.exists():
                inputs[key] = path
        outputs = {key: workdir / key for key in spec.outputs}
        reg.run_step(name, cfg, inputs, outputs)
