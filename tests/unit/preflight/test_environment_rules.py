"""Environment preflight rules T-PRE-23 (tools and databases) and T-PRE-24 (resources)
(issue M1-11; GIPPA_SPEC.md §7.0 Environment rows; ADR-0005 items 3-4).

All probes are fakes (``EnvProbe`` is injected); no real tool is run.

Lead decisions used here: the only resource warning is RAM < 24 GB when ``promoters: genome``
(the caller passes that need as ``need_ram_bytes``); free disk never warns, whatever
``need_disk_bytes`` says (ADR-0005 item 4). Per the M1-12 review the check holds the verdict only,
no raw RAM/disk numbers, so ``preflight.json`` is identical across machines. Memory sizes are GiB
(2**30 bytes), as ``scripts/probe_env.py`` reports Slurm memory.
"""

from __future__ import annotations

import dataclasses
import json
from collections.abc import Mapping, Sequence

import pytest

from gippa.preflight.checks import Check
from gippa.preflight.environment import EnvProbe, Requirement, environment_checks, requirements

GIB = 1024**3
LEVELS = ("ok", "info", "warn", "error")
KINDS = ("binary", "env", "db", "model")
ENV_RULES = ("environment.tools", "environment.resources")


# ---------------------------------------------------------------- helpers


def found_all(reqs: Sequence[Requirement], **override: str | None) -> dict[str, str | None]:
    """Every requirement found at its minimum version (or "1.0"), then ``override`` applied."""
    found: dict[str, str | None] = {r.name: r.min_version or "1.0" for r in reqs}
    found.update(override)
    return found


def probe(found: Mapping[str, str | None], ram_gib: float = 64, disk_gib: float = 500) -> EnvProbe:
    return EnvProbe(
        found=found, free_disk_bytes=int(disk_gib * GIB), total_ram_bytes=int(ram_gib * GIB)
    )


def env_checks(
    reqs: Sequence[Requirement],
    p: EnvProbe,
    need_disk_bytes: int = 0,
    need_ram_bytes: int = 0,
) -> dict[str, Check]:
    out = environment_checks(reqs, p, need_disk_bytes, need_ram_bytes)
    ids = [c.rule_id for c in out]
    assert len(ids) == len(set(ids)), f"rule reported twice: {ids}"
    assert set(ids) == set(ENV_RULES), ids
    for c in out:
        assert c.input == "environment", c
        assert c.level in LEVELS
        assert c.rule, f"{c.rule_id}: empty human rule"
        json.dumps(c.value)
    return {c.rule_id: c for c in out}


def core_reqs() -> list[Requirement]:
    return requirements(["core"], "targeted")


# ---------------------------------------------------------------- requirements


def test_requirement_is_frozen_dataclass() -> None:
    r = Requirement(name="bakta", module="core", kind="binary", min_version="1.12")
    with pytest.raises(dataclasses.FrozenInstanceError):
        r.name = "x"  # type: ignore[misc]


def test_core_requires_bakta_db_v6() -> None:
    # run-wiring decision 4a: the bakta binary comes from the tool registry, not requirements()
    reqs = {r.name: r for r in core_reqs()}
    assert "bakta_db" in reqs
    assert reqs["bakta_db"].kind == "db"
    assert reqs["bakta_db"].module == "core"
    assert reqs["bakta_db"].min_version is not None
    assert reqs["bakta_db"].min_version.split(".")[0] == "6"


def test_requirements_well_formed() -> None:
    reqs = requirements(["core", "mobile", "bgc", "function"], "genome")
    assert reqs
    names = [r.name for r in reqs]
    assert len(names) == len(set(names)), names
    for r in reqs:
        assert r.kind in KINDS, r
        assert r.module in ("core", "mobile", "bgc", "function"), r


def test_optional_modules_off_not_required() -> None:
    reqs = core_reqs()
    assert {r.module for r in reqs} == {"core"}
    lowered = [r.name.lower() for r in reqs]
    assert not any("genomad" in n for n in lowered), lowered
    assert not any("antismash" in n for n in lowered), lowered


def test_mobile_and_bgc_modules_add_their_tools() -> None:
    reqs = requirements(["core", "mobile", "bgc"], "targeted")
    by_module: dict[str, list[str]] = {}
    for r in reqs:
        by_module.setdefault(r.module, []).append(r.name.lower())
    assert any("genomad" in n for n in by_module.get("mobile", [])), by_module
    assert any("antismash" in n for n in by_module.get("bgc", [])), by_module


def test_function_module_requires_go_kegg_cog_databases() -> None:
    assert not [r for r in core_reqs() if r.module == "function"]
    func = [r for r in requirements(["core", "function"], "targeted") if r.module == "function"]
    dbs = [r.name.lower() for r in func if r.kind == "db"]
    assert any(n == "go" or n.startswith("go") for n in dbs), dbs
    assert any("kegg" in n for n in dbs), dbs
    assert any("cog" in n for n in dbs), dbs


# ---------------------------------------------------------------- T-PRE-23 tools and databases


def test_T_PRE_23_pass() -> None:
    reqs = core_reqs()
    c = env_checks(reqs, probe(found_all(reqs, bakta="1.12.2", bakta_db="6.0")))
    tools = c["environment.tools"]
    assert tools.level == "ok", tools
    assert tools.exit_code == 0


def test_T_PRE_23_pass_optional_tools_absent_when_modules_off() -> None:
    reqs = core_reqs()
    found = found_all(reqs, bakta="1.12.2", bakta_db="6.0", genomad=None, antismash=None)
    assert env_checks(reqs, probe(found))["environment.tools"].level == "ok"


def test_T_PRE_23_fail() -> None:
    reqs = core_reqs()
    c = env_checks(reqs, probe(found_all(reqs, bakta="1.12.2", bakta_db=None)))
    tools = c["environment.tools"]
    assert tools.level == "error"
    assert tools.exit_code == 2
    assert tools.value == {"missing": ["bakta_db"]}
    assert "bakta_db" in tools.message


def test_T_PRE_23_fail_wrong_db_version() -> None:
    reqs = core_reqs()
    c = env_checks(reqs, probe(found_all(reqs, bakta="1.12.2", bakta_db="5.1")))
    tools = c["environment.tools"]
    assert tools.level == "error"
    assert tools.exit_code == 2
    assert tools.value == {"wrong_version": {"bakta_db": "5.1 (need 6.x)"}}


def test_T_PRE_23_fail_requirement_absent_from_probe_counts_as_missing() -> None:
    reqs = core_reqs()
    found = found_all(reqs, bakta_db="6.0")
    del found["promotech"]  # the Promotech model (decision 4b); tools come from the registry
    tools = env_checks(reqs, probe(found))["environment.tools"]
    assert tools.level == "error"
    assert tools.exit_code == 2
    assert tools.value == {"missing": ["promotech"]}


def test_T_PRE_23_fail_enabled_optional_module_tool_missing() -> None:
    reqs = requirements(["core", "mobile"], "targeted")
    missing = sorted(r.name for r in reqs if r.module == "mobile")
    assert missing
    found = found_all(reqs, bakta="1.12.2", bakta_db="6.0", **dict.fromkeys(missing))
    tools = env_checks(reqs, probe(found))["environment.tools"]
    assert tools.level == "error"
    assert tools.exit_code == 2
    assert isinstance(tools.value, dict)
    assert sorted(tools.value["missing"]) == missing


# ---------------------------------------------------------------- T-PRE-24 resources


def test_T_PRE_24_pass() -> None:
    reqs = requirements(["core"], "genome")
    found = found_all(reqs, bakta="1.12.2", bakta_db="6.0")
    res = env_checks(reqs, probe(found, ram_gib=64), need_ram_bytes=24 * GIB)[
        "environment.resources"
    ]
    assert res.level == "ok", res
    assert res.exit_code == 0


def test_T_PRE_24_fail() -> None:
    reqs = requirements(["core"], "genome")
    found = found_all(reqs, bakta="1.12.2", bakta_db="6.0")
    c = env_checks(reqs, probe(found, ram_gib=16), need_ram_bytes=24 * GIB)
    res = c["environment.resources"]
    assert res.level == "warn"
    assert res.exit_code == 0
    assert c["environment.tools"].level == "ok"
    # M1-12 lead decision: preflight.json is machine-independent, so the check holds the verdict
    # only; another machine that also lacks RAM yields the identical check (raw numbers → logs).
    other = env_checks(reqs, probe(found, ram_gib=20, disk_gib=7), need_ram_bytes=24 * GIB)
    assert other["environment.resources"] == res
    assert "16" not in json.dumps(res.value) + res.message


@pytest.mark.parametrize("ram_gib", [1, 16, 64])
def test_T_PRE_24_targeted_mode_never_warns_and_omits_raw_numbers(ram_gib: int) -> None:
    reqs = core_reqs()
    found = found_all(reqs, bakta="1.12.2", bakta_db="6.0")
    res = env_checks(reqs, probe(found, ram_gib=ram_gib, disk_gib=1))["environment.resources"]
    assert res.level in ("ok", "info"), res
    # M1-12 lead decision: no raw RAM/disk numbers in the check (machine-independent report).
    assert res == env_checks(reqs, probe(found, ram_gib=128, disk_gib=900))["environment.resources"]


def test_T_PRE_24_low_free_disk_never_warns() -> None:
    reqs = requirements(["core"], "genome")
    found = found_all(reqs, bakta="1.12.2", bakta_db="6.0")
    res = env_checks(
        reqs,
        probe(found, ram_gib=64, disk_gib=1),
        need_disk_bytes=100 * GIB,
        need_ram_bytes=24 * GIB,
    )["environment.resources"]
    assert res.level in ("ok", "info"), res


def test_T_PRE_24_exactly_24_gib_passes() -> None:
    reqs = requirements(["core"], "genome")
    found = found_all(reqs, bakta="1.12.2", bakta_db="6.0")
    res = env_checks(reqs, probe(found, ram_gib=24), need_ram_bytes=24 * GIB)[
        "environment.resources"
    ]
    assert res.level == "ok"
