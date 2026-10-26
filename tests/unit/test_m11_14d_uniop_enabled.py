"""M11-14d: ``annotation.uniop.enabled``, its preflight resolution and ``params.run_uniop``.

docs/specs/M11-14-container-install.md item 1 option (c) "in detail (14d)": Config, Preflight
resolution (the six-row table, ``--publication`` included), Run; examples 8, 8b, 9 and 9b; AC 5
(the item 1 resolution table) and AC 8; edge case "Precomputed operon pairs"; Q-K13 (Khoa
2026-10-25: under ``--publication`` ``auto`` must resolve to ``true`` or exit 2, an explicit
``false`` is accepted). The step reasons, the manifest and the report text are in
``test_m11_14d_operon_source.py``.

Interface pinned here (implementer; the spec leaves these names open):

- ``UniopConfig.enabled: bool | Literal["auto"] = True`` (config model, ``defaults.yaml`` and
  ``src/giap/schemas/config.schema.json``).
- ``giap.annotate.tools.required_tools(cfg)``: ``uniop`` is needed under ``enabled: true`` and
  ``auto`` (preflight must probe it), never under ``false`` or a precomputed operon layer.
- ``giap.preflight.environment.resolve_uniop(cfg, *, installed, publication) -> UniopResolution``,
  pure. ``UniopResolution`` has:
  - ``run_uniop: bool``: whether UNIOP runs;
  - ``level``: ``"ok"`` | ``"warn"`` | ``"error"`` (error: ``environment.tools`` exit 2);
  - ``message: str``: the warning or error text (``""`` when ok);
  - ``cfg: GiapConfig``: the configuration written to ``resolved_config.yaml``, with
    ``annotation.uniop.enabled`` a bool (never ``"auto"``). How it records that ``auto``
    resolved to false (so ``operon_pairs`` can say "UniOP not installed") is the
    implementer's choice; the step tests round-trip it through ``write_params``.
- ``run_preflight(..., profile=..., publication=False)``: a new keyword; ``uniop`` installed
  means not in ``missing_tools`` (14f1). The ``auto`` warning and both errors live on the
  ``environment.tools`` check.
- ``giap run``: ``OUT/work/params.json`` gains ``run_uniop`` (bool), and
  ``resolved_config.yaml`` carries the resolved ``enabled``.
  ``write_params(outdir, cfg, inputs)`` sets ``run_uniop`` from ``cfg.annotation.uniop.enabled``
  of an already resolved ``cfg``.

Fakes only: ``$GIAP_ROOT`` is ``helpers.fake_env.fake_giap_root`` (``skip=["uniop"]``: not
installed). The profile is ``farmshare,apptainer``, where 14f1 takes the images from
``$GIAP_ROOT/containers`` and the envs from ``$GIAP_ROOT/micromamba/envs``. A fake ``nextflow``
on ``PATH`` records its argv. No network; never the real ``$GIAP_DB_DIR`` or ``$GIAP_ROOT``.
"""

from __future__ import annotations

import itertools
import json
import os
import random
import stat
import tomllib
from pathlib import Path
from typing import Any

import pytest
import yaml
from click.testing import CliRunner
from helpers.fake_env import fake_giap_root, write_setup_manifest

from giap.config import GiapConfig, json_schema, resolve_config
from giap.exitcodes import ExitCode

REPO = Path(__file__).resolve().parents[2]
PKG = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]["name"]
DEFAULTS_YAML = REPO / "src" / "giap" / "config" / "defaults.yaml"
CONFIG_SCHEMA = REPO / "src" / "giap" / "schemas" / "config.schema.json"
PROFILE = "farmshare,apptainer"

#: spec item 1 (c), row ``auto | no | no`` (the prefix up to the item 3 command)
AUTO_WARNING = (
    "uniop not installed (unlicensed, not distributed); operons use the distance model; "
    f"run: {PKG} setup --uniop"
)
#: spec item 1 (c), row ``auto | no | yes``; example 8b
AUTO_PUBLICATION_ERROR = (
    "uniop not installed; set annotation.uniop.enabled: false to publish with the "
    "distance-model operons (see Methods)"
)
ENABLED = [True, False, "auto"]
QUIET = ["preflight.min_n50_warn=1", "preflight.min_largest_contig_warn=1"]


def _set(enabled: bool | str) -> str:
    return f"annotation.uniop.enabled={str(enabled).lower()}"


def _cfg(enabled: bool | str, *extra: str) -> GiapConfig:
    return resolve_config(overrides=[*QUIET, _set(enabled), *extra])


def _env_mod() -> Any:
    from giap.preflight import environment

    return environment


def _resolve(cfg: GiapConfig, installed: bool, publication: bool) -> Any:
    fn = getattr(_env_mod(), "resolve_uniop", None)
    assert fn is not None, "giap.preflight.environment.resolve_uniop is not implemented"
    return fn(cfg, installed=installed, publication=publication)


# ---------------------------------------------------------------- config


def test_default_is_true_in_the_model() -> None:
    assert GiapConfig().annotation.uniop.enabled is True


def test_default_is_true_in_defaults_yaml() -> None:
    data = yaml.safe_load(DEFAULTS_YAML.read_text(encoding="utf-8"))
    assert data["annotation"]["uniop"]["enabled"] is True, data["annotation"]["uniop"]


def test_default_resolved_config_is_true() -> None:
    assert resolve_config().annotation.uniop.enabled is True


@pytest.mark.parametrize(("raw", "want"), [("true", True), ("false", False), ("auto", "auto")])
def test_set_accepts_the_three_values(raw: str, want: bool | str) -> None:
    cfg = resolve_config(overrides=[f"annotation.uniop.enabled={raw}"])
    assert cfg.annotation.uniop.enabled == want
    assert type(cfg.annotation.uniop.enabled) is type(want)


def test_config_file_accepts_auto(tmp_path: Path) -> None:
    path = tmp_path / "c.yaml"
    path.write_text("annotation:\n  uniop:\n    enabled: auto\n", encoding="utf-8")
    cfg = resolve_config(path)
    assert cfg.annotation.uniop.enabled == "auto"
    assert cfg.annotation.uniop.pair_threshold == 0.5  # deep merge keeps the siblings


@pytest.mark.parametrize("raw", ["maybe", "AUTO_", "[true]", "{a: 1}"])
def test_other_values_are_rejected(raw: str) -> None:
    with pytest.raises(ValueError):
        resolve_config(overrides=[f"annotation.uniop.enabled={raw}"])


def test_publication_does_not_change_enabled() -> None:
    assert resolve_config(publication=True).annotation.uniop.enabled is True
    cfg = resolve_config(overrides=[_set("auto")], publication=True)
    assert cfg.annotation.uniop.enabled == "auto"


def _uniop_schema(schema: dict[str, Any]) -> dict[str, Any]:
    """The ``UniopConfig`` object of the config JSON Schema (``$defs`` or inline)."""
    defs = schema.get("$defs", {})
    if "UniopConfig" in defs:
        node = defs["UniopConfig"]
        assert isinstance(node, dict)
        return node
    ann = schema["properties"]["annotation"]
    ann = defs[ann["$ref"].rsplit("/", 1)[1]] if "$ref" in ann else ann
    uni = ann["properties"]["uniop"]
    node = defs[uni["$ref"].rsplit("/", 1)[1]] if "$ref" in uni else uni
    assert isinstance(node, dict)
    return node


def _accepts(node: dict[str, Any], value: object) -> bool:
    import jsonschema

    try:
        jsonschema.Draft202012Validator(node).validate(value)
    except jsonschema.ValidationError:
        return False
    return True


@pytest.mark.parametrize("source", ["model", "file"])
def test_json_schema_has_enabled_true_false_auto(source: str) -> None:
    schema = (
        json_schema()
        if source == "model"
        else json.loads(CONFIG_SCHEMA.read_text(encoding="utf-8"))
    )
    props = _uniop_schema(schema)["properties"]
    assert "enabled" in props, sorted(props)
    node = props["enabled"]
    assert node.get("default") is True, node
    for ok in (True, False, "auto"):
        assert _accepts(node, ok), (ok, node)
    for bad in ("maybe", 1.5, None):
        assert not _accepts(node, bad), (bad, node)


# ---------------------------------------------------------------- required tools


def _required(cfg: GiapConfig) -> set[str]:
    from giap.annotate.tools import required_tools

    return {t.name for t in required_tools(cfg)}


@pytest.mark.parametrize(("enabled", "needed"), [(True, True), ("auto", True), (False, False)])
def test_required_tools_follow_enabled(enabled: bool | str, needed: bool) -> None:
    assert ("uniop" in _required(_cfg(enabled))) is needed


@pytest.mark.parametrize("enabled", ENABLED)
def test_required_tools_precomputed_operon_pairs_ignores_enabled(
    enabled: bool | str, tmp_path: Path
) -> None:
    layer = tmp_path / "operon_pairs.tsv"
    layer.write_text("gene_a\tgene_b\tuniop_p\tsource\n", encoding="utf-8")
    cfg = _cfg(enabled, f"annotation.precomputed.operon_pairs={layer}")
    assert "uniop" not in _required(cfg)


def test_required_tools_enabled_changes_only_uniop() -> None:
    on, off = _required(_cfg(True)), _required(_cfg(False))
    assert on - off == {"uniop"} and not off - on, (on, off)


# ---------------------------------------------------------------- resolution (pure; item 1 table)

#: (enabled, installed, publication) -> (level, run_uniop or None when not pinned)
TABLE: dict[tuple[bool | str, bool, bool], tuple[str, bool | None]] = {}
for _pub in (False, True):
    TABLE[(True, True, _pub)] = ("ok", True)
    TABLE[(True, False, _pub)] = ("error", None)
    TABLE[("auto", True, _pub)] = ("ok", True)
    TABLE[(False, True, _pub)] = ("ok", False)
    TABLE[(False, False, _pub)] = ("ok", False)
TABLE[("auto", False, False)] = ("warn", False)
TABLE[("auto", False, True)] = ("error", False)
CELLS = sorted(TABLE, key=str)


def test_table_covers_every_combination() -> None:
    assert set(TABLE) == set(itertools.product(ENABLED, (True, False), (True, False)))


@pytest.mark.parametrize(("enabled", "installed", "publication"), CELLS)
def test_resolution_table(enabled: bool | str, installed: bool, publication: bool) -> None:
    level, run = TABLE[(enabled, installed, publication)]
    res = _resolve(_cfg(enabled), installed, publication)
    assert res.level == level, (enabled, installed, publication, res)
    if run is not None:
        assert res.run_uniop is run, res
    if level == "ok":
        assert res.message == "", res


@pytest.mark.parametrize(("enabled", "installed", "publication"), CELLS)
def test_resolved_cfg_enabled_is_a_bool_equal_to_run_uniop(
    enabled: bool | str, installed: bool, publication: bool
) -> None:
    res = _resolve(_cfg(enabled), installed, publication)
    if res.level == "error":
        return
    got = res.cfg.annotation.uniop.enabled
    assert isinstance(got, bool), got
    assert got is res.run_uniop


@pytest.mark.parametrize(("enabled", "installed", "publication"), CELLS)
def test_uniop_never_runs_when_not_installed(
    enabled: bool | str, installed: bool, publication: bool
) -> None:
    res = _resolve(_cfg(enabled), installed, publication)
    if not installed and res.level != "error":
        assert res.run_uniop is False, res


def test_auto_not_installed_warning_text() -> None:
    res = _resolve(_cfg("auto"), installed=False, publication=False)
    assert AUTO_WARNING in res.message, res.message


def test_auto_not_installed_publication_error_text() -> None:
    res = _resolve(_cfg("auto"), installed=False, publication=True)
    assert AUTO_PUBLICATION_ERROR in res.message, res.message


def test_true_not_installed_error_names_uniop() -> None:
    res = _resolve(_cfg(True), installed=False, publication=False)
    assert "uniop" in res.message, res.message


def test_resolution_keeps_every_other_key() -> None:
    cfg = _cfg("auto", "annotation.uniop.pair_threshold=0.6", "seed=11")
    res = _resolve(cfg, installed=False, publication=False)
    a = res.cfg.model_dump(mode="json")
    b = cfg.model_dump(mode="json")
    a["annotation"]["uniop"].pop("enabled")
    b["annotation"]["uniop"].pop("enabled")
    extra = set(a["annotation"]["uniop"]) - set(b["annotation"]["uniop"])
    for key in extra:  # the implementer's record of the resolution, if any
        a["annotation"]["uniop"].pop(key)
    assert a == b


def test_resolution_is_idempotent() -> None:
    once = _resolve(_cfg("auto"), installed=False, publication=False)
    twice = _resolve(once.cfg, installed=False, publication=False)
    assert twice.level in {"ok", "warn"}
    assert twice.run_uniop is False
    assert twice.cfg.model_dump(mode="json") == once.cfg.model_dump(mode="json")


@pytest.mark.parametrize("enabled", ENABLED)
@pytest.mark.parametrize("publication", [False, True])
def test_precomputed_operon_pairs_needs_no_uniop(
    enabled: bool | str, publication: bool, tmp_path: Path
) -> None:
    """Edge case "Precomputed operon pairs": UniOP is not required and ``enabled`` is ignored."""
    layer = tmp_path / "operon_pairs.tsv"
    layer.write_text("gene_a\tgene_b\tuniop_p\tsource\n", encoding="utf-8")
    cfg = _cfg(enabled, f"annotation.precomputed.operon_pairs={layer}")
    res = _resolve(cfg, installed=False, publication=publication)
    assert res.level == "ok", res
    assert res.run_uniop is False


# ---------------------------------------------------------------- run_preflight (examples 8-9b)

PROTEIN = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDA"


@pytest.fixture
def inputs(tmp_path: Path) -> Any:
    from giap.preflight.run import PreflightInputs

    rng = random.Random(0)
    g = tmp_path / "in" / "g.fa"
    g.parent.mkdir(parents=True, exist_ok=True)
    g.write_text(
        f">chr1 [topology=circular]\n{''.join(rng.choice('ACGT') for _ in range(120_000))}\n"
    )
    p = tmp_path / "in" / "proteins.faa"
    p.write_text(f">P1\n{PROTEIN}\n")
    return PreflightInputs(genome=g, proteins=p)


@pytest.fixture
def db(tmp_path: Path) -> Path:
    return write_setup_manifest(tmp_path / "db")


def _root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch, *, uniop: bool) -> Path:
    root = fake_giap_root(tmp_path / "giap_root", skip=() if uniop else ("uniop",))
    monkeypatch.setenv("GIAP_ROOT", str(root))
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR", raising=False)
    return root


def _preflight(
    inputs: Any, db: Path, tmp_path: Path, enabled: bool | str, publication: bool, *extra: str
) -> tuple[Any, Any]:
    from giap.preflight.run import run_preflight

    report, _ = run_preflight(inputs, _cfg(enabled, *extra), None, "run", tmp_path / "OUT", db,
                              profile=PROFILE, publication=publication)  # fmt: skip
    tools = [c for c in report.checks if c.rule_id == "environment.tools"]
    assert len(tools) == 1, tools
    return report, tools[0]


def _missing(check: Any) -> set[str]:
    value = check.value if isinstance(check.value, dict) else {}
    return set(value.get("missing", []))


def _code(report: Any) -> int:
    from giap.preflight.run import exit_code

    return int(exit_code(report, False))


@pytest.mark.parametrize("enabled", [True, "auto"])
@pytest.mark.parametrize("publication", [False, True])
def test_preflight_installed_passes(
    inputs: Any, db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    enabled: bool | str, publication: bool,
) -> None:  # fmt: skip
    _root(tmp_path, monkeypatch, uniop=True)
    report, check = _preflight(inputs, db, tmp_path, enabled, publication)
    assert check.level == "ok", check
    assert _code(report) == ExitCode.SUCCESS


def test_example_8_auto_not_installed_warns_and_passes(
    inputs: Any, db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _root(tmp_path, monkeypatch, uniop=False)
    report, check = _preflight(inputs, db, tmp_path, "auto", False)
    assert check.level == "warn", check
    assert AUTO_WARNING in check.message, check.message
    assert "uniop" not in _missing(check), check.value
    assert _code(report) == ExitCode.SUCCESS
    from giap.preflight.run import exit_code

    assert exit_code(report, True) == ExitCode.WARNINGS_STRICT  # --strict turns it into exit 3


def test_example_8b_auto_not_installed_publication_exits_2(
    inputs: Any, db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _root(tmp_path, monkeypatch, uniop=False)
    report, check = _preflight(inputs, db, tmp_path, "auto", True)
    assert check.level == "error" and check.exit_code == 2, check
    assert AUTO_PUBLICATION_ERROR in check.message, check.message
    assert _code(report) == ExitCode.ENVIRONMENT


@pytest.mark.parametrize("publication", [False, True])
def test_example_9_true_not_installed_exits_2(
    inputs: Any, db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch, publication: bool
) -> None:
    _root(tmp_path, monkeypatch, uniop=False)
    report, check = _preflight(inputs, db, tmp_path, True, publication)
    assert check.level == "error" and check.exit_code == 2, check
    assert "uniop" in _missing(check), check.value
    assert _code(report) == ExitCode.ENVIRONMENT


@pytest.mark.parametrize("installed", [False, True])
@pytest.mark.parametrize("publication", [False, True])
def test_example_9b_false_passes_with_or_without_publication(
    inputs: Any, db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch,
    installed: bool, publication: bool,
) -> None:  # fmt: skip
    _root(tmp_path, monkeypatch, uniop=installed)
    report, check = _preflight(inputs, db, tmp_path, False, publication)
    assert check.level == "ok", check
    assert "uniop" not in _missing(check)
    assert _code(report) == ExitCode.SUCCESS


def test_auto_warning_does_not_hide_other_missing_tools(
    inputs: Any, db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Another tool missing is still an error; uniop is not among the missing under ``auto``."""
    root = fake_giap_root(tmp_path / "giap_root", skip=("uniop", "ostir"))
    monkeypatch.setenv("GIAP_ROOT", str(root))
    monkeypatch.delenv("NXF_APPTAINER_CACHEDIR", raising=False)
    report, check = _preflight(inputs, db, tmp_path, "auto", False)
    assert check.level == "error", check
    assert "ostir" in _missing(check) and "uniop" not in _missing(check), check.value
    assert _code(report) == ExitCode.ENVIRONMENT


def test_publication_keyword_defaults_to_false(
    inputs: Any, db: Path, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    from giap.preflight.run import run_preflight

    _root(tmp_path, monkeypatch, uniop=False)
    report, _ = run_preflight(inputs, _cfg("auto"), None, "run", tmp_path / "OUT", db,
                              profile=PROFILE)  # fmt: skip
    assert _code(report) == ExitCode.SUCCESS


# ---------------------------------------------------------------- write_params


@pytest.mark.parametrize("enabled", [True, False])
def test_write_params_run_uniop_from_the_resolved_value(tmp_path: Path, enabled: bool) -> None:
    from giap.runner.nextflow import write_params

    params, resolved = write_params(tmp_path / "OUT", _cfg(enabled), {})
    data = json.loads(params.read_text(encoding="utf-8"))
    assert data.get("run_uniop") is enabled, data.get("run_uniop")
    back = yaml.safe_load(resolved.read_text(encoding="utf-8"))
    assert back["annotation"]["uniop"]["enabled"] is enabled


@pytest.mark.parametrize(("installed", "want"), [(True, True), (False, False)])
def test_write_params_of_a_resolved_auto(tmp_path: Path, installed: bool, want: bool) -> None:
    from giap.runner.nextflow import write_params

    res = _resolve(_cfg("auto"), installed=installed, publication=False)
    params, resolved = write_params(tmp_path / "OUT", res.cfg, {})
    assert json.loads(params.read_text(encoding="utf-8"))["run_uniop"] is want
    back = yaml.safe_load(resolved.read_text(encoding="utf-8"))
    assert back["annotation"]["uniop"]["enabled"] is want
    assert resolve_config(resolved).annotation.uniop.enabled is want


def test_nextflow_config_declares_run_uniop_default_true() -> None:
    """``params.run_uniop`` is a declared param (``nextflow.config`` mirrors the defaults)."""
    import re

    text = (REPO / "workflow" / "nextflow.config").read_text(encoding="utf-8")
    block = text.split("params {", 1)[1].split("\n}", 1)[0]
    assert re.search(r"^\s*run_uniop\s*=\s*true\b", block, re.M), block


def test_annotate_subworkflow_reads_params_run_uniop() -> None:
    text = (REPO / "workflow" / "subworkflows" / "local" / "annotate.nf").read_text(
        encoding="utf-8"
    )
    assert "params.run_uniop" in text


# ---------------------------------------------------------------- giap run (fake nextflow)

FAKE_NF = """#!/bin/sh
printf '%s\\n' "$@" > "$FAKE_NF_ARGS"
exit 0
"""


@pytest.fixture
def genome(tmp_path: Path) -> Path:
    rng = random.Random(3)
    path = tmp_path / "in" / "genome.fa"
    path.parent.mkdir(parents=True, exist_ok=True)
    seq = "".join(rng.choice("ACGT") for _ in range(6000))
    path.write_text(f">chr1 [topology=circular]\n{seq}\n")
    return path


@pytest.fixture
def proteins(tmp_path: Path) -> Path:
    rng = random.Random(7)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    path = tmp_path / "in" / "proteins.faa"
    path.parent.mkdir(parents=True, exist_ok=True)
    recs = [f">prot{i}\nM{''.join(rng.choice(aa) for _ in range(150))}\n" for i in range(3)]
    path.write_text("".join(recs))
    return path


@pytest.fixture
def fake_nf(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    bindir = tmp_path / "fakebin"
    bindir.mkdir()
    exe = bindir / "nextflow"
    exe.write_text(FAKE_NF)
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    args = tmp_path / "nf_args.txt"
    monkeypatch.setenv("PATH", f"{bindir}{os.pathsep}{os.environ.get('PATH', '')}")
    monkeypatch.setenv("FAKE_NF_ARGS", str(args))
    return args


def _giap_run(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path,
    *, uniop: bool, sets: tuple[str, ...] = (), publication: bool = False,
) -> tuple[int, str, Path]:  # fmt: skip
    from giap.cli import cli

    _root(tmp_path, monkeypatch, uniop=uniop)
    db = write_setup_manifest(tmp_path / "db_env")
    monkeypatch.setenv("GIAP_DB_DIR", str(db))
    out = tmp_path / "OUT"
    argv = ["run", "--genome", str(genome), "--proteins", str(proteins), "--outdir", str(out),
            "--profile", PROFILE]  # fmt: skip
    for s in (*QUIET, *sets):
        argv += ["--set", s]
    if publication:
        argv.append("--publication")
    res = CliRunner().invoke(cli, argv)
    if res.exception is not None and not isinstance(res.exception, SystemExit):
        raise res.exception
    text = res.output + (res.stderr if res.stderr_bytes is not None else "")
    return res.exit_code, text, out


def _params(out: Path) -> dict[str, Any]:
    path = out / "work" / "params.json"
    assert path.is_file(), "giap run wrote no params.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def _resolved_enabled(out: Path) -> object:
    doc = yaml.safe_load((out / "resolved_config.yaml").read_text(encoding="utf-8"))
    return doc["annotation"]["uniop"]["enabled"]


def _preflight_messages(out: Path) -> str:
    doc = json.loads((out / "preflight" / "preflight.json").read_text(encoding="utf-8"))
    return "\n".join(str(c.get("message", "")) for c in doc["checks"])


@pytest.mark.usefixtures("fake_nf")
def test_run_default_installed_passes_run_uniop_true(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    code, text, out = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=True)
    assert code == ExitCode.SUCCESS, text
    assert _params(out)["run_uniop"] is True
    assert _resolved_enabled(out) is True


@pytest.mark.usefixtures("fake_nf")
def test_run_example_8_auto_not_installed(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    code, text, out = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=False,
                                sets=(_set("auto"),))  # fmt: skip
    assert code == ExitCode.SUCCESS, text
    assert _params(out)["run_uniop"] is False
    assert _resolved_enabled(out) is False
    assert AUTO_WARNING in _preflight_messages(out)


@pytest.mark.usefixtures("fake_nf")
def test_run_auto_installed_resolves_to_true(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    code, text, out = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=True,
                                sets=(_set("auto"),))  # fmt: skip
    assert code == ExitCode.SUCCESS, text
    assert _params(out)["run_uniop"] is True
    assert _resolved_enabled(out) is True


def test_run_example_8b_auto_publication_exits_2_before_nextflow(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path,
    fake_nf: Path,
) -> None:  # fmt: skip
    code, text, out = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=False,
                                sets=(_set("auto"),), publication=True)  # fmt: skip
    assert code == ExitCode.ENVIRONMENT, text
    assert AUTO_PUBLICATION_ERROR in text + _preflight_messages(out)
    assert not fake_nf.exists(), "nextflow was launched"


def test_run_example_9_true_not_installed_exits_2(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path,
    fake_nf: Path,
) -> None:  # fmt: skip
    code, text, _ = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=False)
    assert code == ExitCode.ENVIRONMENT, text
    assert "uniop" in text
    assert not fake_nf.exists(), "nextflow was launched"


@pytest.mark.usefixtures("fake_nf")
def test_run_example_9b_false_with_publication(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, genome: Path, proteins: Path
) -> None:
    code, text, out = _giap_run(tmp_path, monkeypatch, genome, proteins, uniop=False,
                                sets=(_set(False),), publication=True)  # fmt: skip
    assert code == ExitCode.SUCCESS, text
    assert _params(out)["run_uniop"] is False
    assert _resolved_enabled(out) is False
