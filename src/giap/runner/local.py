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


# Root keys a run must stage when a step needs them; other roots (``proteins_long``,
# ``reference_peptides``, ``contaminants``, ``peptides``) are optional and prune silently.
REQUIRED_ROOTS: tuple[str, ...] = ("bakta_json",)


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


# An ``annotation.precomputed`` layer arrives as its M2-02 TSV (giap.annotate.precomputed), not as
# the M4-06 contract Parquet: layer key -> the contract keys it then supplies as that TSV. M4 steps
# ignoring such a key still run (``gene_table`` does), but a step *requiring* one cannot, so it is
# left out of the plan (MVP-01). Under Nextflow the same steps are skipped because ANNOTATE has no
# C7/C8 Parquet to emit then.
_LAYER_TSV_KEYS: Mapping[str, tuple[str, ...]] = {"regions": ("regions", "region_genes")}


def _tsv_only(cfg: GiapConfig) -> frozenset[str]:
    """Contract keys that only a precomputed layer's TSV would supply in this run."""
    return frozenset(
        key
        for layer, keys in _LAYER_TSV_KEYS.items()
        if _layer(cfg, layer) is not None
        for key in keys
    )


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
    tsv_only = _tsv_only(cfg)
    pruned: list[str] = []
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
        if not _enabled(cfg, spec):
            continue
        if not tsv_only.isdisjoint(spec.inputs):  # needs a contract table only a layer supplies
            continue
        if feasible is not None and not feasible(name):
            pruned.append(name)
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
    if roots is not None:
        _check_required_roots(cfg, specs, producers, roots, pruned)
    return reg.topo_order(planned), supplied


def _check_required_roots(
    cfg: GiapConfig,
    specs: Mapping[str, reg.StepSpec],
    producers: Mapping[str, list[str]],
    roots: Collection[str],
    pruned: Sequence[str],
) -> None:
    """Raise ``LayerMissing`` when a pruned step would be feasible with one more REQUIRED root
    (M3-17 decision 13); steps missing an optional root stay pruned."""
    for key in REQUIRED_ROOTS:
        if key in roots:
            continue
        ok = _feasible(cfg, specs, producers, {*roots, key})
        for name in sorted(pruned):
            if ok(name):  # exit 1 naming the key
                raise LayerMissing(
                    f"annotation.precomputed.{key} (or --annotation <Bakta JSON>; step {name})"
                )


def plan_local(
    cfg: GiapConfig, targets: Sequence[str], roots: Collection[str] | None = None
) -> list[str]:
    """Python steps needed for ``targets``, in run order; raises ``LayerMissing`` or
    ``ToolUnavailable``."""
    return _walk(cfg, targets, roots)[0]


def stage_roots(workdir: Path, roots: Mapping[str, Path]) -> None:
    """Copy root input files (normalised by preflight) to ``workdir / key``. Only ``preflight``,
    whose ``workdir / preflight`` is preflight's own work directory, goes inside that directory
    (M4-15 decisions 19, 21); any other key onto a directory raises ``IsADirectoryError``."""
    workdir.mkdir(parents=True, exist_ok=True)
    for key, src in sorted(roots.items()):
        dst = workdir / key
        if dst.is_dir() and key == "preflight":
            dst = dst / src.name
        elif dst.is_dir():
            raise IsADirectoryError(f"cannot stage root '{key}': {dst} is a directory")
        shutil.copyfile(src, dst)


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
    # with staged roots an optional input is a root of this run or an output of a planned step,
    # never a stale ``workdir / key`` left by an earlier run (M3-17 review decision 8)
    fresh = None if roots is None else {*roots, *(k for n in order for k in reg.get(n).outputs)}
    for name in order:
        spec = reg.get(name)
        inputs: dict[str, reg.InputPath] = {}
        for key in (*spec.inputs, *spec.optional_inputs):
            # a configured layer always wins over a stale ``workdir / key`` from an earlier run
            path = supplied.get(key, workdir / key)
            stale = fresh is not None and key not in supplied and key not in fresh
            if key in spec.inputs or (path.exists() and not stale):
                inputs[key] = path
        outputs = {key: workdir / key for key in spec.outputs}
        reg.run_step(name, cfg, inputs, outputs)


# §9.1 published directories: work key -> OUT subdirectory
PUBLISHED: dict[str, str] = {"annotation": "annotation", "mapping": "mapping"}


def publish(workdir: Path, outdir: Path) -> list[Path]:
    """Copy the published step outputs of ``workdir`` to ``outdir`` (§9.1), like Nextflow."""
    done: list[Path] = []
    for key, sub in sorted(PUBLISHED.items()):
        src, dst = workdir / key, outdir / sub
        if src.is_dir():
            shutil.copytree(src, dst, dirs_exist_ok=True)
            done.append(dst)
    return done


def contaminant_library(db_dir: Path | None) -> Path | None:
    """``contaminants/contaminants.fasta`` when ``db_dir``'s manifest lists ``contaminants``
    (M3-15), else ``None``."""
    from giap.preflight.environment import read_manifest

    if db_dir is None:
        return None
    found = read_manifest(str(db_dir)) or {}
    lib = db_dir / "contaminants" / "contaminants.fasta"
    return lib if "contaminants" in found and lib.is_file() else None
