"""Static checks on the external-tool Nextflow modules (M2-07 example 7 and acceptance criteria).

GIAP_SPEC.md §6.3 "Process conventions" (tag, label, conda + container, versions.yml, stub:,
explicit memory for memory-heavy tools), §7.1.1 (one process per tool), R5 (licence in meta.yml),
R14 / M2-06 lead decision 3 (farmshare: per-task cpus <= 16, memory <= 64 GB, cpus >=
ceil(memory / 4 GB), queueSize x max cpus <= 128). No Nextflow needed: reads committed files only.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path
from typing import Any

import pytest
import yaml

from giap.annotate.tools import TOOLS, read_containers_config

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
MODULES = WF / "modules" / "local"
SUBWF = WF / "subworkflows" / "local" / "annotate.nf"
CONTAINERS = WF / "conf" / "containers.config"
BASE = WF / "conf" / "base.config"
FARMSHARE = WF / "conf" / "farmshare.config"
STUB_BAKTA = WF / "assets" / "stub" / "bakta.json"

# M2-07 "Public interface" table: tool -> label.
LABELS = {
    "bakta": "process_high",
    "uniop": "process_low",
    "promotech": "process_medium",
    "promoter_calculator": "process_single",
    "terminatornet": "process_medium",
    "transtermhp": "process_single",
    "bactermfinder": "process_medium",
    "ostir": "process_medium",
    "genomad": "process_high",
    "antismash": "process_high",
}
TOOL_NAMES = sorted(LABELS)
# "Memory-heavy processes set memory explicitly (Bakta 48.GB, geNOMAD 48.GB, antiSMASH 32.GB)".
EXPLICIT_MEMORY_GB = {"bakta": 48, "genomad": 48, "antismash": 32}
# M2-01 lead decision 1: no published image yet, conda env only (listed in containers.config).
CONDA_ONLY = {
    "uniop",
    "promotech",
    "promoter_calculator",
    "terminatornet",
    "bactermfinder",
    "ostir",
}
FARMSHARE_CPUS_MAX, FARMSHARE_MEM_GB_MAX, MEM_PER_CPU_GB, R14_CPUS = 16, 64, 4, 128


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)^\s*//.*$", "", text)


def _main(tool: str) -> str:
    path = MODULES / tool / "main.nf"
    assert path.is_file(), f"workflow/modules/local/{tool}/main.nf missing"
    return _strip_comments(path.read_text(encoding="utf-8"))


def _directives(text: str) -> str:
    """The process body before ``input:`` (where tag/label/conda/memory live)."""
    m = re.search(r"\bprocess\s+\w+\s*\{(.*?)\binput\s*:", text, flags=re.S)
    assert m, "process block with an input: section not found"
    return m.group(1)


def _block(text: str, name: str) -> str:
    """``script:`` (up to ``stub:``) or ``stub:`` (to the end of the process) block."""
    if name == "script":
        m = re.search(r"\bscript\s*:(.*?)\bstub\s*:", text, flags=re.S)
    else:
        m = re.search(r"\bstub\s*:(.*)", text, flags=re.S)
    assert m, f"no {name}: block"
    return m.group(1)


def _process_names() -> dict[str, Path]:
    out: dict[str, Path] = {}
    for path in sorted(WF.rglob("*.nf")):
        if "tests" in path.relative_to(WF).parts:
            continue
        for name in re.findall(r"(?m)^\s*process\s+(\w+)\s*\{", _strip_comments(path.read_text())):
            out[name] = path
    return out


# ---------------------------------------------------------------- files exist


@pytest.mark.parametrize("tool", TOOL_NAMES)
def test_module_files_exist(tool: str) -> None:
    assert (MODULES / tool / "main.nf").is_file(), f"modules/local/{tool}/main.nf missing"
    assert (MODULES / tool / "meta.yml").is_file(), f"modules/local/{tool}/meta.yml missing"


def test_annotate_subworkflow_exists_and_main_calls_it() -> None:
    assert SUBWF.is_file(), "workflow/subworkflows/local/annotate.nf missing"
    sub = _strip_comments(SUBWF.read_text())
    assert re.search(r"(?m)^\s*workflow\s+ANNOTATE\s*\{", sub), (
        "no `workflow ANNOTATE` in annotate.nf"
    )
    main = _strip_comments((WF / "main.nf").read_text())
    assert re.search(
        r"include\s*\{[^}]*\bANNOTATE\b[^}]*\}\s*from\s*['\"]\./subworkflows/local/annotate", main
    )
    assert re.search(r"\bANNOTATE\s*\(", main), "main.nf does not call ANNOTATE"


def test_annotate_scatter_gather_is_order_stable() -> None:
    """Promotech chunks / Promoter Calculator regions: flatten() then collect(sort: true)."""
    sub = _strip_comments(SUBWF.read_text())
    assert "flatten()" in sub
    assert len(re.findall(r"collect\s*\(\s*sort\s*:\s*true", sub)) >= 2, (
        "both scatters (promotech chunks, promoter_calculator regions) must gather sorted"
    )


# ---------------------------------------------------------------- example 7: conventions


@pytest.mark.parametrize("tool", TOOL_NAMES)
def test_process_is_named_after_the_tool(tool: str) -> None:
    assert re.search(rf"(?m)^\s*process\s+{tool.upper()}\s*\{{", _main(tool)), (
        f"modules/local/{tool}/main.nf must define process {tool.upper()}"
    )


@pytest.mark.parametrize("tool", TOOL_NAMES)
def test_process_has_tag_label_conda_stub_versions(tool: str) -> None:
    text = _main(tool)
    head = _directives(text)
    assert re.search(r"(?m)^\s*tag\s+[\"']\$\{meta_id\}[\"']", head), f'{tool}: tag "${{meta_id}}"'
    labels = re.findall(r"(?m)^\s*label\s+['\"](process_\w+)['\"]", head)
    assert labels == [LABELS[tool]], f"{tool}: label {labels} != {LABELS[tool]}"
    assert re.search(rf"(?m)^\s*conda\s+[\"']\$\{{projectDir\}}/envs/{tool}\.yml[\"']", head), (
        f'{tool}: conda "${{projectDir}}/envs/{tool}.yml"'
    )
    assert "stub:" in text, f"{tool}: no stub: block"
    for block in ("script", "stub"):
        body = _block(text, block)
        assert "versions.yml" in body, f"{tool}: {block} block writes no versions.yml"
        assert "${task.process}" in body, f"{tool}: versions.yml must be keyed by task.process"
        assert re.search(rf"(?m)^\s*{tool}\s*:", body), (
            f"{tool}: versions.yml has no `{tool}:` entry"
        )
    assert re.search(r"path\s*\(?\s*['\"]versions\.yml['\"]", text), (
        f"{tool}: versions.yml not emitted"
    )


@pytest.mark.parametrize("tool", TOOL_NAMES)
def test_conda_env_file_is_the_registry_one(tool: str) -> None:
    assert TOOLS[tool].env_file == f"workflow/envs/{tool}.yml"
    assert (REPO / TOOLS[tool].env_file).is_file()


@pytest.mark.parametrize("tool", TOOL_NAMES)
def test_module_does_not_hardcode_an_unpinned_container(tool: str) -> None:
    for image in re.findall(r"(?m)^\s*container\s+['\"]([^'\"]+)['\"]", _main(tool)):
        assert re.search(r"@sha256:[0-9a-f]{64}$", image), (
            f"{tool}: container {image} not by digest"
        )


@pytest.mark.parametrize("tool", TOOL_NAMES)
def test_meta_yml_records_name_and_licence(tool: str) -> None:
    meta = yaml.safe_load((MODULES / tool / "meta.yml").read_text(encoding="utf-8"))
    assert isinstance(meta, dict) and meta.get("name") == tool, f"{tool}: meta.yml name"

    def licences(node: Any) -> list[str]:
        if isinstance(node, dict):
            out = [str(x) for k, v in node.items() if k in {"licence", "license"}
                   for x in (v if isinstance(v, list) else [v])]  # fmt: skip
            return out + [s for v in node.values() for s in licences(v)]
        if isinstance(node, list):
            return [s for v in node for s in licences(v)]
        return []

    found = licences(meta)
    assert TOOLS[tool].licence in found, (
        f"{tool}: meta.yml licence {found} != {TOOLS[tool].licence}"
    )


# ---------------------------------------------------------------- registry and containers


# Process names that differ from their TOOLS key (BLAST+ is one tool; M3-13 runs its tblastn).
PROCESS_TOOL = {"tblastn": "blast"}


def _tool_key(process: str) -> str:
    return PROCESS_TOOL.get(process.lower(), process.lower())


def test_every_process_maps_to_a_tools_key() -> None:
    procs = _process_names()
    assert set(LABELS) <= {p.lower() for p in procs}, sorted(procs)
    for name, path in procs.items():
        if name == "GIAP_STEP":
            continue
        assert _tool_key(name) in TOOLS, f"process {name} ({path}) is not a TOOLS key"
        assert not TOOLS[_tool_key(name)].in_process, f"{name}: in-process tools get no process"


def test_every_with_name_in_containers_config_matches_a_process() -> None:
    procs = {p.upper() for p in _process_names()}
    for name in read_containers_config(CONTAINERS):
        assert _tool_key(name) in TOOLS, f"containers.config withName {name} is not a TOOLS key"
        if TOOLS[_tool_key(name)].module == "mapping":
            continue  # miniprot / BLAST processes belong to §7.2 (M3); pins land in M2-01
        assert name in procs, f"containers.config withName {name} matches no process"


@pytest.mark.parametrize("tool", TOOL_NAMES)
def test_every_tool_process_has_a_container_pin_or_is_conda_only(tool: str) -> None:
    pins = read_containers_config(CONTAINERS)
    if tool in CONDA_ONLY:
        assert tool.upper() not in pins
        header = CONTAINERS.read_text(encoding="utf-8")
        assert re.search(rf"(?m)^//.*\b{tool}\b", header), f"{tool}: not listed as conda-only"
    else:
        assert tool.upper() in pins, f"{tool}: no withName pin in containers.config"
        assert re.fullmatch(r"docker://\S+@sha256:[0-9a-f]{64}", pins[tool.upper()])


# ---------------------------------------------------------------- resources (§6.3, R14)


def _label_resources() -> dict[str, tuple[int, int]]:
    text = _strip_comments(BASE.read_text())
    out = {}
    for label, cpus, mem in re.findall(
        r"withLabel\s*:\s*['\"](process_\w+)['\"]\s*\{\s*cpus\s*=\s*(\d+)\s*;\s*memory\s*=\s*\{?\s*(\d+)\s*\.GB",
        text,
    ):
        out[label] = (int(cpus), int(mem))
    assert set(out) >= set(LABELS.values()), out
    return out


def _explicit_memory_gb(tool: str) -> int | None:
    m = re.search(r"(?m)^\s*memory\s*[=]?\s*\{?\s*['\"]?(\d+)\s*\.?\s*GB", _directives(_main(tool)))
    return int(m.group(1)) if m else None


@pytest.mark.parametrize("tool", sorted(EXPLICIT_MEMORY_GB))
def test_memory_heavy_processes_declare_memory(tool: str) -> None:
    assert _explicit_memory_gb(tool) == EXPLICIT_MEMORY_GB[tool], (
        f"{tool}: memory must be set explicitly to {EXPLICIT_MEMORY_GB[tool]}.GB"
    )


@pytest.mark.parametrize("tool", TOOL_NAMES)
def test_farmshare_task_resources_respect_r14(tool: str) -> None:
    """cpus = max(label cpus, ceil(mem / 4 GB)) with mem clamped to 64 GB, after a x2 retry too."""
    head = _directives(_main(tool))
    m = re.search(r"(?m)^\s*cpus\s*[=]?\s*(\d+)", head)
    label_cpus, label_mem = _label_resources()[LABELS[tool]]
    cpus = int(m.group(1)) if m else label_cpus
    assert cpus <= FARMSHARE_CPUS_MAX, f"{tool}: cpus {cpus} > {FARMSHARE_CPUS_MAX}"
    mem = _explicit_memory_gb(tool) or label_mem
    assert mem <= FARMSHARE_MEM_GB_MAX, f"{tool}: memory {mem} GB > {FARMSHARE_MEM_GB_MAX} GB"
    for attempt in (1, 2):
        clamped = min(mem * attempt, FARMSHARE_MEM_GB_MAX)
        need = max(cpus, math.ceil(clamped / MEM_PER_CPU_GB))
        assert need <= FARMSHARE_CPUS_MAX, f"{tool} attempt {attempt}: {need} CPUs"
    assert "maxForks" not in head or "bigmem" in head, f"{tool}: maxForks only for a bigmem label"


def test_farmshare_concurrency_cap_still_holds() -> None:
    text = _strip_comments(FARMSHARE.read_text())
    queue = int(re.search(r"queueSize\s*=\s*(\d+)", text).group(1))  # type: ignore[union-attr]
    assert queue * FARMSHARE_CPUS_MAX <= R14_CPUS


# ---------------------------------------------------------------- stub outputs


def test_stub_bakta_json_asset() -> None:
    """Q2 default: a small simulator-derived Bakta JSON that the real parser accepts."""
    from giap.annotate.bakta import parse_bakta_json

    assert STUB_BAKTA.is_file(), "workflow/assets/stub/bakta.json missing"
    assert STUB_BAKTA.stat().st_size <= 5 * 1024, "stub Bakta JSON must stay <= 5 kB"
    json.loads(STUB_BAKTA.read_text(encoding="utf-8"))
    ann = parse_bakta_json(STUB_BAKTA)
    assert len(ann.model.genes) >= 1, "stub Bakta JSON has no genes"


@pytest.mark.parametrize(
    ("tool", "needle"),
    [
        ("bakta", "assets/stub/bakta.json"),
        ("terminatornet", "NO TERMINATORS FOUND"),
        ("antismash", '"records"'),
        ("uniop", "uniop.pred"),
        ("ostir", "ostir.csv"),
        ("transtermhp", ".tt"),
    ],
)
def test_stub_block_writes_the_documented_minimal_output(tool: str, needle: str) -> None:
    assert needle in _block(_main(tool), "stub"), f"{tool} stub: {needle!r} missing"
