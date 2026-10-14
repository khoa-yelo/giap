"""Static checks on the Nextflow skeleton files (M2-06 public interface and edge cases).

No Nextflow needed: these read the committed files only. Behaviour of the configuration is tested
in ``tests/integration/test_nf_config.py`` (``needs_nextflow``).
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest
import yaml

REPO = Path(__file__).resolve().parents[2]
WF = REPO / "workflow"
GIAP_STEP = WF / "modules" / "local" / "giap_step"
SBATCH = REPO / "scripts" / "slurm" / "giap_run.sbatch"

FILES = [
    "main.nf",
    "nextflow.config",
    "conf/base.config",
    "conf/test.config",
    "conf/farmshare.config",
    "conf/farmshare_local.config",
    "conf/standard.config",
    "modules/local/giap_step/main.nf",
    "modules/local/giap_step/meta.yml",
    "tests/main.nf.test",
]


@pytest.mark.parametrize("rel", FILES)
def test_skeleton_file_exists(rel: str) -> None:
    assert (WF / rel).is_file(), f"workflow/{rel} missing"


def _strip_comments(text: str) -> str:
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    return re.sub(r"(?m)//.*$", "", text)


def test_no_plugins_block_anywhere_in_workflow() -> None:
    """§6.3: no Nextflow plugins in version 1."""
    files = [p for p in WF.rglob("*") if p.suffix in {".nf", ".config"} and p.is_file()]
    assert files, "no Nextflow files under workflow/"
    for path in files:
        code = _strip_comments(path.read_text())
        assert not re.search(r"\bplugins\s*\{", code), f"plugins block in {path}"
        assert "nf-schema" not in code and "nf-validation" not in code, path


def test_nextflow_config_includes_every_conf() -> None:
    text = _strip_comments((WF / "nextflow.config").read_text())
    for name in ("base", "test", "farmshare", "farmshare_local", "standard"):
        assert re.search(rf"includeConfig\s+['\"][^'\"]*conf/{name}\.config['\"]", text), name


def test_nextflow_config_declares_params_and_profiles() -> None:
    text = _strip_comments((WF / "nextflow.config").read_text())
    for p in ("threads", "memory", "outdir", "modules", "annotation"):
        assert re.search(rf"\b{p}\b", text), f"params.{p} not declared"
    for prof in (
        "apptainer",
        "conda",
        "docker",
        "farmshare",
        "farmshare_local",
        "test",
        "standard",
    ):
        assert re.search(rf"\b{prof}\s*\{{", text), f"profile {prof} not defined"


def test_base_config_defines_the_four_labels() -> None:
    text = _strip_comments((WF / "conf" / "base.config").read_text())
    for label in ("process_single", "process_low", "process_medium", "process_high"):
        assert re.search(rf"withLabel\s*:\s*['\"]?{label}\b", text), label
    assert "resourceLimits" in text + (WF / "nextflow.config").read_text()


def test_giap_step_module_calls_giap_step_and_has_stub_and_versions() -> None:
    text = (GIAP_STEP / "main.nf").read_text()
    assert re.search(r"process\s+GIAP_STEP\s*\{", text)
    assert re.search(r"\btag\b", text)
    assert re.search(r"\blabel\s+['\"]process_(single|low|medium|high)['\"]", text)
    assert "giap _step" in text
    assert "--config" in text and "--in" in text and "--out" in text
    assert re.search(r"(?m)^\s*stub:", text), "GIAP_STEP has no stub: block"
    assert "versions.yml" in text


def test_giap_step_meta_yml() -> None:
    meta = yaml.safe_load((GIAP_STEP / "meta.yml").read_text())
    assert isinstance(meta, dict)
    assert str(meta.get("name", "")).lower() == "giap_step"
    assert "input" in meta and "output" in meta
    assert "versions" in yaml.safe_dump(meta["output"])


def test_main_nf_includes_giap_step_or_annotate() -> None:
    """main.nf calls ANNOTATE (M2-07) when present, else only the GIAP_STEP smoke."""
    text = (WF / "main.nf").read_text()
    assert re.search(
        r"include\s*\{[^}]*\b(GIAP_STEP|ANNOTATE)\b[^}]*\}\s*from\s*['\"]\./(modules|subworkflows)/local/",
        text,
    )
    assert re.search(r"(?m)^\s*workflow\s*\{", text), "no entry workflow"


# ---------------------------------------------------------------- giap_run.sbatch (§3.8)


def _directives(text: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for m in re.finditer(r"(?m)^#SBATCH\s+--([\w-]+)=(\S+)", text):
        out[m.group(1)] = m.group(2)
    return out


def test_giap_run_sbatch_header() -> None:
    text = SBATCH.read_text()
    assert text.startswith("#!/bin/bash\n")
    d = _directives(text)
    assert d["job-name"].startswith("giap-")
    assert d["partition"] == "normal"
    assert d["qos"] == "normal"
    assert d["cpus-per-task"] == "2"
    assert d["mem"].upper() in {"8G", "8GB", "8192M"}
    assert d["time"] == "2-00:00:00"
    assert d["output"].endswith("/local/logs/%x.%j.out")
    assert "set -euo pipefail" in text
    assert re.search(r"source\s+\"?\$\{GIAP_REPO[^}]*\}/scripts/farmshare_env\.sh", text)
    assert "mktemp -d /tmp/giap." in text


def test_giap_run_sbatch_uses_farmshare_apptainer() -> None:
    text = _strip_comments(SBATCH.read_text())
    assert "farmshare,apptainer" in text
    assert "giap run" in text or "nextflow run" in text
