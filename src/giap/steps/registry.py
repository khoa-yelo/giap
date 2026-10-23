"""Step registry, topological order and ``run_step`` (GIAP_SPEC.md §6.2, §6.6; M2-02).

Every pipeline step declares its file inputs and outputs by *key* (the work-directory file name,
see the M2-index contract table). Python steps (``kind="py"``) carry a function and run through
``giap _step``; external steps (``kind="ext"``) are Nextflow processes and are only declared here
so the graph is complete. A step with a ``precomputed`` key is replaced by that
``annotation.precomputed`` layer when the layer is given.
"""

from __future__ import annotations

import hashlib
import heapq
import json
import logging
import shutil
import traceback
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np

from giap.config import GiapConfig
from giap.rng import for_step

StepKind = Literal["py", "ext"]
StepModule = Literal["core", "mobile", "bgc", "function", "opt_in"]
InputPath = Path | list[Path]  # a repeated ``--in key=...`` gives a list (e.g. Promotech chunks)


@dataclass(frozen=True, slots=True)
class StepContext:
    """What a Python step sees: config, its declared files, its own RNG stream and a logger."""

    cfg: GiapConfig
    inputs: Mapping[str, InputPath]
    outputs: Mapping[str, Path]
    rng: np.random.Generator
    log: logging.Logger


StepFn = Callable[[StepContext], None]


@dataclass(frozen=True, slots=True)
class StepSpec:
    name: str
    kind: StepKind
    inputs: tuple[str, ...]
    outputs: tuple[str, ...]
    optional_inputs: tuple[str, ...]
    module: StepModule
    precomputed: str | None  # annotation.precomputed key that replaces this step
    func: StepFn | None  # None for ext
    #: subset of ``outputs`` a caller may leave unwired; the step then does not write it
    optional_outputs: tuple[str, ...] = ()


class StepError(Exception):
    """A step function raised (``giap _step`` exit code 4)."""

    def __init__(self, name: str, cause: BaseException) -> None:
        super().__init__(f"step '{name}' failed: {type(cause).__name__}: {cause}")
        self.name = name
        self.cause = cause
        self.traceback = "".join(traceback.format_exception(cause))

    def write_log(self, log_dir: Path) -> Path:
        """Write ``log_dir/<step>.log`` with the traceback; returns its path."""
        log_dir.mkdir(parents=True, exist_ok=True)
        path = log_dir / f"{self.name}.log"
        path.write_text(f"{self}\n\n{self.traceback}")
        return path


class WiringError(Exception):
    """Unknown step or key, missing input or missing output (``giap _step`` exit code 5)."""


_REGISTRY: dict[str, StepSpec] = {}


def _register(spec: StepSpec) -> None:
    if spec.name in _REGISTRY:
        raise ValueError(f"step '{spec.name}' already registered")
    for other in _REGISTRY.values():
        dup = sorted(set(spec.outputs) & set(other.outputs))
        if dup:
            raise ValueError(
                f"step '{spec.name}': output key {', '.join(dup)} already declared by "
                f"step '{other.name}'"
            )
    _REGISTRY[spec.name] = spec


def step(
    name: str,
    *,
    inputs: Sequence[str],
    outputs: Sequence[str],
    optional_inputs: Sequence[str] = (),
    optional_outputs: Sequence[str] = (),
    module: StepModule = "core",
    precomputed: str | None = None,
) -> Callable[[StepFn], StepFn]:
    """Decorator registering a Python step; returns the function unchanged.

    ``optional_outputs`` must be a subset of ``outputs``: a caller may leave such a key
    unwired (the step then skips that artifact), which lets one step be wired in stages while
    the keys stay reserved in the graph.
    """
    _check_optional(name, outputs, optional_outputs)

    def deco(func: StepFn) -> StepFn:
        _register(
            StepSpec(
                name=name,
                kind="py",
                inputs=tuple(inputs),
                outputs=tuple(outputs),
                optional_inputs=tuple(optional_inputs),
                module=module,
                precomputed=precomputed,
                func=func,
                optional_outputs=tuple(optional_outputs),
            )
        )
        return func

    return deco


def external(
    name: str,
    *,
    inputs: Sequence[str],
    outputs: Sequence[str],
    optional_inputs: Sequence[str] = (),
    optional_outputs: Sequence[str] = (),
    module: StepModule = "core",
    precomputed: str | None = None,
) -> None:
    """Declare an external (Nextflow process) step."""
    _check_optional(name, outputs, optional_outputs)
    _register(
        StepSpec(
            name=name,
            kind="ext",
            inputs=tuple(inputs),
            outputs=tuple(outputs),
            optional_inputs=tuple(optional_inputs),
            module=module,
            precomputed=precomputed,
            func=None,
            optional_outputs=tuple(optional_outputs),
        )
    )


def _check_optional(name: str, outputs: Sequence[str], optional: Sequence[str]) -> None:
    extra = sorted(set(optional) - set(outputs))
    if extra:
        raise ValueError(f"step '{name}': optional_outputs not in outputs: {', '.join(extra)}")


def get(name: str) -> StepSpec:
    """The registered step ``name``; ``KeyError`` if unknown."""
    try:
        return _REGISTRY[name]
    except KeyError:
        raise KeyError(f"unknown step '{name}'") from None


def all_steps() -> list[StepSpec]:
    """Every registered step, sorted by name."""
    return [_REGISTRY[n] for n in sorted(_REGISTRY)]


#: top-level config sections read only by the named steps (M7-15 lead decision 7). Any other
#: section, and every other top-level value except ``_NOT_HASHED``, is in every step's slice,
#: so a change re-runs too much rather than too little.
SECTION_READERS: Mapping[str, frozenset[str]] = {
    "variants": frozenset({"variants"}),
    "novel_orfs": frozenset({"novel_orfs"}),
}
_NOT_HASHED = frozenset({"threads"})  # resources only; outputs never depend on it (§6.6)


def config_slice(dumped: Mapping[str, object], name: str) -> dict[str, object]:
    """The part of a dumped config (``model_dump(mode="json")``) that step ``name`` reads. Pure."""
    mine = frozenset({name})
    return {
        k: v
        for k, v in dumped.items()
        if k not in _NOT_HASHED and name in SECTION_READERS.get(k, mine)
    }


def config_hashes(dumped: Mapping[str, object], names: Iterable[str]) -> dict[str, str]:
    """Per step, the SHA-256 of its ``config_slice`` (canonical JSON), for the task hash. Pure."""
    return {
        n: hashlib.sha256(
            json.dumps(config_slice(dumped, n), sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        for n in names
    }


def topo_order(names: Iterable[str], members: Mapping[str, str] | None = None) -> list[str]:
    """Order ``names`` so that a step producing a key runs before every step consuming it.

    Edges come from output keys to (required or optional) input keys among ``names`` only (a key
    in ``members``, a file inside another output key, is produced by that key's producer); ties
    are broken by name (Kahn's algorithm with a min-heap), so the result does not depend on
    registration or argument order. A cycle raises ``ValueError`` naming the steps involved.
    """
    specs = [get(n) for n in sorted(set(names))]
    producer: dict[str, str] = {}
    for s in specs:
        for key in s.outputs:
            producer.setdefault(key, s.name)
    for key, parent in (members or {}).items():
        if parent in producer:
            producer.setdefault(key, producer[parent])
    succ: dict[str, set[str]] = {s.name: set() for s in specs}
    indeg: dict[str, int] = dict.fromkeys(succ, 0)
    for s in specs:
        for key in (*s.inputs, *s.optional_inputs):
            p = producer.get(key)
            if p is not None and p != s.name and s.name not in succ[p]:
                succ[p].add(s.name)
                indeg[s.name] += 1
    heap = [n for n, d in indeg.items() if d == 0]
    heapq.heapify(heap)
    order: list[str] = []
    while heap:
        n = heapq.heappop(heap)
        order.append(n)
        for m in succ[n]:
            indeg[m] -= 1
            if indeg[m] == 0:
                heapq.heappush(heap, m)
    if len(order) < len(succ):
        stuck = sorted(n for n, d in indeg.items() if d > 0)
        raise ValueError(f"step graph has a cycle among: {', '.join(stuck)}")
    return order


def _check_keys(
    spec: StepSpec, inputs: Mapping[str, InputPath], outputs: Mapping[str, Path]
) -> None:
    allowed_in = set(spec.inputs) | set(spec.optional_inputs)
    for key in inputs:
        if key not in allowed_in:
            raise WiringError(f"step '{spec.name}': undeclared input key '{key}'")
    for key in outputs:
        if key not in spec.outputs:
            raise WiringError(f"step '{spec.name}': undeclared output key '{key}'")
    for key in spec.inputs:
        if key not in inputs:
            raise WiringError(f"step '{spec.name}': required input '{key}' not given")
    for key, given in inputs.items():
        for path in given if isinstance(given, list) else [given]:
            if not path.exists():
                raise WiringError(f"step '{spec.name}': input '{key}' does not exist: {path}")
    for key in spec.outputs:
        if key not in outputs and key not in spec.optional_outputs:
            raise WiringError(f"step '{spec.name}': output '{key}' not given")


def _remove(path: Path) -> None:
    if path.is_dir() and not path.is_symlink():
        shutil.rmtree(path)
    else:
        path.unlink(missing_ok=True)


def _tmp_path(name: str, path: Path) -> Path:
    """Temporary sibling of ``path`` (same file system, same suffix) for an atomic rename."""
    return path.with_name(f".tmp.{name}.{path.name}")


def run_step(
    name: str, cfg: GiapConfig, inputs: Mapping[str, InputPath], outputs: Mapping[str, Path]
) -> None:
    """Run the Python step ``name`` on the given files (§6.2).

    Declared outputs are removed first, so a stale file never passes the output check; the step
    writes temporary paths that are renamed into place only when it succeeds and wrote them all.
    Wiring problems (unknown step, undeclared key, missing input, missing output after the run,
    an external step) raise ``WiringError``; an exception from the step raises ``StepError``.
    """
    spec = _REGISTRY.get(name)
    if spec is None:
        raise WiringError(f"unknown step '{name}'")
    if spec.func is None:
        raise WiringError(f"step '{name}' is external and has no Python function")
    _check_keys(spec, inputs, outputs)
    tmp = {key: _tmp_path(name, path) for key, path in outputs.items()}
    for key, path in outputs.items():
        path.parent.mkdir(parents=True, exist_ok=True)
        _remove(path)
        _remove(tmp[key])
    ctx = StepContext(
        cfg=cfg,
        inputs=dict(inputs),
        outputs=dict(tmp),
        rng=for_step(name, cfg.seed),
        log=logging.getLogger(f"giap.steps.{name}"),
    )
    try:
        try:
            spec.func(ctx)
        except Exception as exc:
            raise StepError(name, exc) from exc
        for key, path in tmp.items():
            if not path.exists():
                raise WiringError(
                    f"step '{name}': declared output '{key}' was not written: {outputs[key]}"
                )
    except BaseException:
        for path in tmp.values():
            _remove(path)
        raise
    for key, path in tmp.items():
        path.replace(outputs[key])
