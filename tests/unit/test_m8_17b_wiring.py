"""M8-17b static wiring: schema, Nextflow configs, the ``requires_r`` marker and the R5 audit
docstring (docs/specs/M8-17-r-fgsea-engine.md §9, §10, A7; lead decisions on the M8-17a tests 7;
lead decision 9).

- ``enrichment_gsea.schema.json``: the flags pattern gains ``log2err_unestimated``; the ``method``
  description drops "reserved for M8-11"; the enum stays ``fgsea_simple``, ``fgsea_multilevel``.
- ``workflow/conf/prebuilt_envs.config`` sets ``env { GIPPA_FGSEA_PREFIX = "${GIPPA_ROOT ?:
  projectDir/../local}/micromamba/envs/gippa-fgsea" }``. Lead decision 9 (conda-only): no
  ``GIPPA_FGSEA_IMAGE`` anywhere under ``workflow/``.
- ``workflow/conf/base.config``: ``withName: 'FUNCTION_STEP'`` 1 CPU, 6 GB, 4 h; no ``maxForks``.
- ``pyproject.toml`` registers ``requires_r``.
- ``scripts/ci/check_no_rscript.py`` (from #192) documents the one R entry point.
"""

from __future__ import annotations

import json
import math
import re
import tomllib
from pathlib import Path

import numpy as np
import pytest
from helpers.fgsea_switch_data import RecordingEngine, validate_gsea
from helpers.gsea_data import gene_sets, singleton_frame, spread

REPO = Path(__file__).resolve().parents[2]
SCHEMA = REPO / "src" / "gippa" / "schemas" / "tables" / "enrichment_gsea.schema.json"
CONF = REPO / "workflow" / "conf"
FLAGS = (
    "ci_unestimated",
    "unestimated",
    "log2err_unestimated",
    "size_filtered",
    "genomic_cluster",
    "detection_confounded",
)


def _field(name: str) -> dict:
    d = json.loads(SCHEMA.read_text(encoding="utf-8"))
    (f,) = [f for f in d["fields"] if f["name"] == name]
    return f


# ---------------------------------------------------------------- schema (§10, A7)


@pytest.mark.parametrize(
    "cell",
    [
        "ci_unestimated;log2err_unestimated",
        "log2err_unestimated",
        "ci_unestimated;log2err_unestimated;genomic_cluster;detection_confounded",
    ],
)
def test_flags_pattern_accepts_log2err_unestimated(cell: str) -> None:
    pattern = _field("flags")["constraints"]["pattern"]
    assert re.fullmatch(pattern, cell), (pattern, cell)


@pytest.mark.parametrize("flag", FLAGS)
def test_flags_pattern_keeps_every_flag(flag: str) -> None:
    pattern = _field("flags")["constraints"]["pattern"]
    assert re.fullmatch(pattern, f"ci_unestimated;{flag}")


@pytest.mark.parametrize("cell", ["log2err_estimated", "log2err_unestimated;", "bogus"])
def test_flags_pattern_rejects_other_tokens(cell: str) -> None:
    assert not re.fullmatch(_field("flags")["constraints"]["pattern"], cell)


def test_method_description_and_enum() -> None:
    f = _field("method")
    assert "reserved for M8-11" not in f.get("description", "")
    assert f["constraints"]["enum"] == ["fgsea_simple", "fgsea_multilevel"]


def test_example_4_table_validates(tmp_path: Path) -> None:
    """A7: a written table with a ``log2err_unestimated`` row validates."""
    from gippa.config.model import GseaConfig
    from gippa.function.gsea import gsea_table
    from gippa.io.function_tables import write_enrichment_gsea_tsv
    from gippa.rng import for_step

    genes = singleton_frame(spread(60, np.random.default_rng(1)))
    tags = [f"g{i}" for i in range(1, 61)]
    coll = gene_sets(
        "kegg", {f"K0000{k}": [tags[(k - 1) + 6 * i] for i in range(10)] for k in range(1, 5)}
    )
    eng = RecordingEngine(rows={("main:kegg", "K00002"): {"p": 3.2e-15, "log2err": math.nan}})
    df = gsea_table(
        genes, [coll], (10, 500), GseaConfig(), 0.05, for_step("function.gsea", 7), engine=eng
    )
    out = tmp_path / "enrichment_gsea.tsv"
    write_enrichment_gsea_tsv(out, df)
    assert "log2err_unestimated" in out.read_text(encoding="utf-8")
    validate_gsea(out)


# ---------------------------------------------------------------- Nextflow configs (§9)


def _strip_comments(text: str) -> str:
    return "\n".join(line for line in text.splitlines() if not line.lstrip().startswith("//"))


def _block(text: str, opener: str) -> str:
    """The body of the first ``opener { ... }`` block (braces balanced)."""
    m = re.search(opener + r"\s*\{", text)
    assert m is not None, f"no {opener!r} block"
    depth, i = 1, m.end()
    while depth:
        c = text[i]
        depth += c == "{"
        depth -= c == "}"
        i += 1
    return text[m.end() : i - 1]


def test_prebuilt_envs_sets_gippa_fgsea_prefix() -> None:
    text = _strip_comments((CONF / "prebuilt_envs.config").read_text(encoding="utf-8"))
    env = _block(text, r"(?m)^\s*env")
    m = re.search(r"GIPPA_FGSEA_PREFIX\s*=\s*(\"[^\"]*\"|'[^']*')", env)
    assert m is not None, env
    value = m.group(1)
    assert value.rstrip("\"'").endswith("/micromamba/envs/gippa-fgsea"), value
    assert "GIPPA_ROOT" in value and "projectDir" in value and "/../local" in value, value


def test_no_image_launcher_anywhere_in_workflow() -> None:
    """Lead decision 9: fgsea is conda-only; the step never calls ``apptainer exec`` itself."""
    hits = [
        str(p.relative_to(REPO))
        for p in (REPO / "workflow").rglob("*")
        if p.is_file() and p.suffix in {".config", ".nf"} and "GIPPA_FGSEA_IMAGE" in p.read_text()
    ]
    assert hits == []


def test_base_config_function_step_resources() -> None:
    text = _strip_comments((CONF / "base.config").read_text(encoding="utf-8"))
    body = _block(text, r"withName\s*:\s*['\"]FUNCTION_STEP['\"]")
    assert re.search(r"\bcpus\s*=\s*1\b", body), body
    assert re.search(r"\bmemory\s*=\s*\{?\s*6\.GB\b", body), body
    assert re.search(r"\btime\s*=\s*\{?\s*4\.h\b", body), body
    assert "maxForks" not in body


def test_farmshare_gives_function_step_no_max_forks() -> None:
    for name in ("farmshare.config", "farmshare_local.config", "prebuilt_envs.config"):
        text = _strip_comments((CONF / name).read_text(encoding="utf-8"))
        if re.search(r"withName\s*:\s*['\"]FUNCTION_STEP['\"]", text):
            body = _block(text, r"withName\s*:\s*['\"]FUNCTION_STEP['\"]")
            assert "maxForks" not in body, name


# ---------------------------------------------------------------- marker and audit


def test_requires_r_marker_is_registered_in_pyproject() -> None:
    data = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
    markers = data["tool"]["pytest"]["ini_options"]["markers"]
    assert any(m.split(":", 1)[0].strip() == "requires_r" for m in markers), markers


def test_check_no_rscript_docstring_names_the_one_entry_point() -> None:
    path = REPO / "scripts" / "ci" / "check_no_rscript.py"
    assert path.is_file(), "scripts/ci/check_no_rscript.py arrives with #192; merge it first"
    text = path.read_text(encoding="utf-8")
    assert "GIPPA never shells out to R" not in text
    assert "src/ starts R only through workflow/bin/gippa_fgsea.R (ADR-0038)" in " ".join(
        text.split()
    )
