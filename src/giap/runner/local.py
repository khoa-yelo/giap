"""Developer-only local executor (GIAP_SPEC.md §6.4; M2-02).

``giap run --executor local`` runs the Python steps in topological order in one process. External
steps never run here: each must be replaced by its ``annotation.precomputed`` layer, and a missing
layer is an error naming the key, raised before any step runs. Every step key is stored at
``workdir / key``; a consumer of an output of a replaced step receives the precomputed file itself,
even when an older ``workdir / key`` exists, except a ``LAYER_CONVERTERS`` step, which still
runs on the layer to write the contract Parquet. A ``MEMBER_KEYS`` key is a file inside another
step output (M3's tables inside ``mapping/``). Opt-in steps run only when their configuration
enables them (e.g. ``annotation.terminators.bactermfinder``).
"""

from __future__ import annotations

import json
import shutil
import time
from collections.abc import Callable, Collection, Mapping, Sequence
from pathlib import Path
from typing import Any

from giap.config import GiapConfig
from giap.io.work_tables import FILE_NAMES
from giap.steps import registry as reg


class LayerMissing(Exception):
    """An external step has no ``annotation.precomputed`` layer under the local executor."""

    def __init__(self, key: str) -> None:
        super().__init__(f"{key} is required with --executor local")
        self.key = key


#: Run-tail steps: registered, never planned; ``giap run`` calls them after the graph (ADR-0040).
TAIL_STEPS: tuple[str, ...] = ("manifest",)

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


# Steps that still run when their ``annotation.precomputed`` layer is set: they turn the layer's
# M2-02 TSV into the M4-06 contract Parquet (``region_layer`` -> C7/C8, M4-21 lead decision 2),
# so ``context_stats`` and ``sensitivity`` are planned on precomputed runs too. Only the external
# tools they would otherwise read are replaced.
LAYER_CONVERTERS: frozenset[str] = frozenset({"region_layer"})

# Keys that are a file inside another step output (key -> (output key, file name)): M3's final
# tables live in the ``loci`` output directory ``mapping/`` (M4-21 lead decision 4), the way
# workflow/subworkflows/local/context.nf picks them out of the LOCI outputs.
MEMBER_KEYS: Mapping[str, tuple[str, str]] = {
    k: ("mapping", FILE_NAMES[k]) for k in ("peptide_hits", "protein_hits", "unmapped")
}
_PARENT: Mapping[str, str] = {k: parent for k, (parent, _) in MEMBER_KEYS.items()}


def _replacing_layer(cfg: GiapConfig, spec: reg.StepSpec) -> Path | None:
    """The layer that replaces ``spec`` in this run (``None``: the step runs)."""
    if not spec.precomputed or spec.name in LAYER_CONVERTERS:
        return None
    return _layer(cfg, spec.precomputed)


def _producers(specs: Mapping[str, reg.StepSpec]) -> dict[str, list[str]]:
    """Key -> producing steps, a member key produced by whoever produces its output key."""
    producers: dict[str, list[str]] = {}
    for s in specs.values():
        for key in s.outputs:
            producers.setdefault(key, []).append(s.name)
    for key, parent in _PARENT.items():
        if parent in producers:
            producers.setdefault(key, list(producers[parent]))
    return producers


def _key_path(workdir: Path, key: str) -> Path:
    """Where ``key`` lives under ``workdir`` (a member key inside its output directory)."""
    if key in MEMBER_KEYS:
        parent, name = MEMBER_KEYS[key]
        return workdir / parent / name
    return workdir / key


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
        elif _replacing_layer(cfg, spec) is not None:
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
    producers = _producers(specs)
    feasible = None if roots is None else _feasible(cfg, specs, producers, roots)
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
        if name in TAIL_STEPS or not _enabled(cfg, spec):
            continue
        if feasible is not None and not feasible(name):
            pruned.append(name)
            continue
        layer = _replacing_layer(cfg, spec)
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
    return reg.topo_order(planned, _PARENT), supplied


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
    timings: list[dict[str, Any]] | None = None,
    live: Mapping[str, Path] | None = None,
) -> None:
    """Plan (validating every layer first), stage ``roots``, then run each Python step under
    ``workdir``. Without ``roots`` every planned step must find its inputs in ``workdir``.
    ``timings`` receives ``{name, status, wall_time_s, cpu_time_s}`` per step run (a failed step
    included, as ``failed``) for the manifest. ``live`` keys are passed as these paths, never
    copied (the run's growing events log)."""
    order, supplied = _walk(cfg, targets, None if roots is None else set(roots))
    supplied = {**supplied, **(live or {})}
    if roots is not None:
        stage_roots(workdir, roots)
    workdir.mkdir(parents=True, exist_ok=True)
    # with staged roots an optional input is a root of this run or an output of a planned step,
    # never a stale ``workdir / key`` left by an earlier run (M3-17 review decision 8)
    fresh = None if roots is None else {*roots, *(k for n in order for k in reg.get(n).outputs)}
    if fresh is not None:
        fresh |= {k for k, parent in _PARENT.items() if parent in fresh}
    for name in order:
        spec = reg.get(name)
        inputs: dict[str, reg.InputPath] = {}
        for key in (*spec.inputs, *spec.optional_inputs):
            # a configured layer always wins over a stale ``workdir / key`` from an earlier run
            path = supplied.get(key, _key_path(workdir, key))
            stale = fresh is not None and key not in supplied and key not in fresh
            if key in spec.inputs or (path.exists() and not stale):
                inputs[key] = path
        if name == "plots":  # M9-17: the analysis steps' support directories, as under Nextflow
            staged = [workdir / k for k in SUPPORT_OUTPUTS if (workdir / k / "support").is_dir()
                      and (fresh is None or k in fresh)]  # fmt: skip
            if staged:
                inputs["plots_support"] = staged
        outputs = {key: workdir / key for key in spec.outputs}
        wall, cpu, status = time.perf_counter(), time.process_time(), "failed"
        try:
            reg.run_step(name, cfg, inputs, outputs)
            status = "completed"
        finally:
            if timings is not None:
                timings.append({"name": name, "status": status,
                                "wall_time_s": time.perf_counter() - wall,
                                "cpu_time_s": time.process_time() - cpu})  # fmt: skip


# §9.1 published directories: work key -> OUT subdirectory
PUBLISHED: dict[str, str] = {
    "annotation": "annotation",
    "mapping": "mapping",
    "plots": "plots",
    "report": "report",
}
#: §9.1 published files: work key -> OUT path, the ``key=path`` names the Nextflow subworkflows
#: give GIAP_STEP (pinned against workflow/ by tests/unit/test_local_publish_names.py)
PUBLISHED_FILES: dict[str, str] = {
    "autocorrelation": "tables/autocorrelation.tsv",
    "category_summary": "tables/category_summary.tsv",
    "detectability_summary": "genetable/detection_summary_parts/detectability.parquet",
    "detection_summary": "tables/detection_summary.tsv",
    "determinants_attribution": "tables/determinants_attribution.tsv",
    "determinants_features": "tables/determinants_features.tsv",
    "determinants_model": "tables/determinants_model.tsv",
    "enrichment_detection": "tables/enrichment_detection.tsv",
    "enrichment_gsea": "tables/enrichment_gsea.tsv",
    "gene_context": "context/gene_context.parquet",
    "gene_table": "genetable/genes.parquet",
    "genes_pre": "genetable/genes_pre.parquet",
    "genes_tsv": "tables/genes.tsv",
    "island_diagnostics": "context/island_diagnostics.parquet",
    "islands": "tables/islands.tsv",
    "novel_orfs": "tables/novel_orfs.tsv",
    "novel_orfs_faa": "tables/novel_orfs.faa",
    "novel_orfs_gff3": "tables/novel_orfs.gff3",
    "annotation_patch": "tables/annotation_patch.gff3",
    "operon_flags": "tables/operon_flags.tsv",
    "operons": "tables/operons.tsv",
    "pairs": "tables/pairs.tsv",
    "regions_stats": "tables/regions_stats.tsv",
    "pileup": "tables/pileup.tsv.gz",
    "sensitivity": "tables/sensitivity.tsv",
    "tier_summary": "genetable/detection_summary_parts/tier.parquet",
    "tu_folds": "genetable/tu_folds.parquet",
    "variants": "tables/variants.tsv",
    "variants_vcf": "tables/variants.vcf",
    "variants_qc": "work/variants_qc.parquet",  # an intermediate: never published (M7-14)
}


#: M9-17: the analysis steps' optional support directories; ``publish`` merges their
#: ``support/`` subdirectories (``giap.io.plot_support.write_support``) into ``OUT/plots/support``
SUPPORT_OUTPUTS: tuple[str, ...] = (
    "context_stats_support", "determinants_support", "function_support",
)  # fmt: skip


#: outputs named ``work/<file>`` are intermediates: GIAP_STEP and :func:`publish` skip them
UNPUBLISHED_PREFIX = "work/"


def publish(workdir: Path, outdir: Path) -> list[Path]:
    """Copy the published step outputs of ``workdir`` to ``outdir`` (§9.1), like Nextflow."""
    done: list[Path] = []
    for key, sub in sorted(PUBLISHED.items()):
        src, dst = workdir / key, outdir / sub
        if src.is_dir():
            shutil.copytree(src, dst, dirs_exist_ok=True)
            done.append(dst)
    for key, rel in sorted(PUBLISHED_FILES.items()):
        if rel.startswith(UNPUBLISHED_PREFIX):
            continue
        src, dst = workdir / key, outdir / rel
        if src.is_file():
            dst.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(src, dst)
            done.append(dst)
    support = outdir / "plots" / "support"
    for key in SUPPORT_OUTPUTS:
        if (workdir / key / "support").is_dir():
            shutil.copytree(workdir / key / "support", support, dirs_exist_ok=True)
            done += [] if support in done else [support]
    return done


def write_datapackage(outdir: Path) -> Path | None:
    """Copy ``schemas/datapackage.json`` and the table schemas it names to ``OUT/schemas/``
    when a resource it lists was published (§9.1, §9.5; M7-03 decision 3); resource paths are
    ``tables/...`` under ``OUT``. Returns the written path, or ``None``."""
    from importlib.resources import files

    src = files("giap.schemas")
    package = json.loads(src.joinpath("datapackage.json").read_text(encoding="utf-8"))
    resources = package["resources"]
    if not any((outdir / r["path"]).is_file() for r in resources):
        return None
    dst = outdir / "schemas"
    for rel in sorted({str(r["schema"]) for r in resources if isinstance(r.get("schema"), str)}):
        (dst / rel).parent.mkdir(parents=True, exist_ok=True)
        (dst / rel).write_bytes(src.joinpath(rel).read_bytes())
    (dst / "datapackage.json").write_bytes(src.joinpath("datapackage.json").read_bytes())
    return dst / "datapackage.json"


def contaminant_library(db_dir: Path | None) -> Path | None:
    """``contaminants/contaminants.fasta`` when ``db_dir``'s manifest lists ``contaminants``
    (M3-15), else ``None``."""
    from giap.preflight.environment import read_manifest

    if db_dir is None:
        return None
    found = read_manifest(str(db_dir)) or {}
    lib = db_dir / "contaminants" / "contaminants.fasta"
    return lib if "contaminants" in found and lib.is_file() else None
