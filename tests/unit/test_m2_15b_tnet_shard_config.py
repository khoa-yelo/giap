"""M2-15b: Nextflow params, resources, modules and the run-review mapping of the sharded
TerminatorNet scan (docs/specs/M2-15b-terminatornet-sharding.md items 5, 6; lead decisions 1, 4).

Static pins, as tests/unit/test_tool_resources_pins.py does for the other ``withName`` blocks:
``params.terminatornet_sharded = true`` and ``params.terminatornet_chunk_cpus = 4`` (not
defaults.yaml keys); ``params.stub.terminatornet_chunks``; CHUNK asks for
``params.terminatornet_chunk_cpus`` CPUs and 4 GB x attempt, and on FarmShare
``max(chunk cpus, ceil(memory / 4000 MB))`` CPUs and 2 h x attempt; annotate.nf picks
``TERMINATORNET_SHARDED`` or the old ``TERMINATORNET`` from ``params.terminatornet_sharded``.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import pytest

from gippa.runner import nextflow as nf

WF = Path(__file__).resolve().parents[2] / "workflow"
CONF = WF / "conf"
MODS = WF / "modules" / "local"


def _block(text: str, selector: str) -> str:
    m = re.search(rf"withName:\s*'{re.escape(selector)}'\s*\{{(.*)\}}\s*$", text, re.M)
    assert m, f"no withName: '{selector}' block"
    return m[1]


def _params() -> str:
    text = (WF / "nextflow.config").read_text()
    m = re.search(r"^params\s*\{(.*?)^\}", text, re.M | re.S)
    assert m
    return m[1]


def test_params_defaults() -> None:
    p = _params()
    assert re.search(r"^\s*terminatornet_sharded\s*=\s*true\b", p, re.M), p
    assert re.search(r"^\s*terminatornet_chunk_cpus\s*=\s*4\b", p, re.M), p
    stub = re.search(r"^\s*stub\s*=\s*\[(.*)\]", p, re.M)
    assert stub and re.search(r"\bterminatornet_chunks\s*:\s*1\b", stub[1]), p


def test_params_are_not_defaults_yaml_keys() -> None:
    root = Path(__file__).resolve().parents[2] / "src" / "gippa"
    for y in root.rglob("defaults.yaml"):
        assert (
            "terminatornet_sharded" not in y.read_text()
            and "terminatornet_chunk_cpus" not in y.read_text()
        )


def test_chunk_request_base_config() -> None:
    body = _block((CONF / "base.config").read_text(), "TERMINATORNET_CHUNK")
    assert re.search(r"\bcpus\s*=\s*\{\s*params\.terminatornet_chunk_cpus\s*\}", body), body
    assert re.search(r"memory\s*=\s*\{\s*4\.GB\s*\*\s*task\.attempt\s*\}", body), body


def test_chunk_request_farmshare() -> None:
    body = _block((CONF / "farmshare.config").read_text(), "TERMINATORNET_CHUNK")
    assert re.search(
        r"Math\.max\(\s*params\.terminatornet_chunk_cpus\s*,\s*Math\.ceil\(task\.memory\.toBytes\(\)\s*/\s*4000\.MB",
        body,
    ), body
    assert re.search(r"time\s*=\s*\{\s*2\.h\s*\*\s*task\.attempt\s*\}", body), body
    # 4 GB -> 2 derived CPUs, the 8 GB retry -> 3: a 4-CPU chunk keeps asking for 4
    assert [max(4, math.ceil(4 * 1024 * a / 4000)) for a in (1, 2)] == [4, 4]


@pytest.mark.parametrize(
    ("name", "label"),
    [
        ("terminatornet_plan", "process_single"),
        ("terminatornet_merge", "process_single"),
        ("terminatornet_chunk", None),
    ],
)
def test_modules(name: str, label: str | None) -> None:
    text = (MODS / name / "main.nf").read_text()
    assert re.search(rf"^process\s+{name.upper()}\s*\{{", text, re.M), text
    assert 'conda "${projectDir}/envs/terminatornet.yml"' in text
    assert "gippa_terminatornet_shard" in text
    assert re.search(r"^\s*stub:", text, re.M)
    if label:
        assert f"label '{label}'" in text
    if name == "terminatornet_plan":
        assert "path 'plan.tsv'" in text and "params.stub.terminatornet_chunks" in text
    if name == "terminatornet_chunk":
        assert 'tag "${chunk_id}"' in text
        assert re.search(r"chunk\s+\$\{chunk_id\}\s+--cpus\s+\$\{task\.cpus\}", text), text
        assert 'path "calls_${chunk_id}.pkl"' in text
    if name == "terminatornet_merge":
        assert re.search(r"path\s+'terminatornet\.tsv',\s*emit:\s*tsv", text), text
        assert "path 'versions.yml'" in text and "GIPPA_COMMIT" in text
        assert "NO TERMINATORS FOUND" in text


def test_subworkflow_and_annotate_switch() -> None:
    sub = (WF / "subworkflows" / "local" / "terminatornet_sharded.nf").read_text()
    assert re.search(r"^workflow\s+TERMINATORNET_SHARDED\s*\{", sub, re.M)
    assert "splitCsv(header: true, sep: '\\t')" in sub and ".collect()" in sub
    for p in ("TERMINATORNET_PLAN", "TERMINATORNET_CHUNK", "TERMINATORNET_MERGE"):
        assert p in sub
    ann = (WF / "subworkflows" / "local" / "annotate.nf").read_text()
    assert re.search(
        r"include\s*\{\s*TERMINATORNET_SHARDED\s*\}\s*from\s*'\./terminatornet_sharded'", ann
    ), ann
    assert "params.terminatornet_sharded" in ann
    assert "include { TERMINATORNET } from '../../modules/local/terminatornet/main'" in ann


def test_prebuilt_env_selector_widened() -> None:
    text = (CONF / "prebuilt_envs.config").read_text()
    assert re.search(r"withName:\s*'TERMINATORNET\.\*'.*gippa-terminatornet", text), text


def _trace(tmp: Path, name: str) -> Path:
    work = tmp / "OUT" / "work" / "ab" / "cdef0123"
    work.mkdir(parents=True)
    t = tmp / "OUT" / "logs" / "trace.tsv"
    t.parent.mkdir(parents=True)
    t.write_text(
        f"task_id\thash\tname\tstatus\texit\tworkdir\n1\tab/cdef01\t{name}\tFAILED\t140\t{work}\n"
    )
    return t


@pytest.mark.parametrize(
    "name",
    [
        "ANNOTATE:TERMINATORNET_SHARDED:TERMINATORNET_PLAN (genome)",
        "ANNOTATE:TERMINATORNET_SHARDED:TERMINATORNET_CHUNK (c00000p0001)",
        "ANNOTATE:TERMINATORNET_SHARDED:TERMINATORNET_MERGE (genome)",
    ],
)
def test_run_review_maps_the_three_processes_to_terminatornet(tmp_path: Path, name: str) -> None:
    got = nf.failed_step(_trace(tmp_path, name))
    assert got is not None and got[0] == "terminatornet", got
