"""Registry of external annotation tools, their pinned environments and containers (§7.1.1).

Pins, digests and licences are recorded in docs/adr/ADR-0008-tool-pins.md. GIAP never imports
these tools (R5); they run as separate processes from ``workflow/envs/<name>.yml`` or by digest.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Literal

import yaml

from giap.config import GiapConfig

Module = Literal["core", "mobile", "bgc", "function", "opt_in", "mapping", "pilot"]


@dataclass(frozen=True, slots=True)
class ToolSpec:
    name: str
    version: str
    licence: str
    env_file: str
    env_name: str
    version_cmd: tuple[str, ...]
    module: Module
    source_repo: str | None = None
    in_process: bool = False


def _pkg_version(dist: str) -> tuple[str, ...]:
    return ("python", "-c", f"import importlib.metadata as m; print(m.version('{dist}'))")


def _git_commit(name: str) -> tuple[str, ...]:
    # build_envs.sh writes the checked-out SHA next to the cloned source.
    return ("sh", "-c", f'cat "$CONDA_PREFIX/share/{name}/GIAP_COMMIT"')


def _tool(
    name: str,
    version: str,
    licence: str,
    version_cmd: tuple[str, ...],
    module: Module,
    source_repo: str | None = None,
) -> ToolSpec:
    return ToolSpec(
        name=name,
        version=version,
        licence=licence,
        env_file=f"workflow/envs/{name}.yml",
        env_name=f"giap-{name}",
        version_cmd=version_cmd,
        module=module,
        source_repo=source_repo,
    )


_GIT = {
    "uniop": ("3ed14dfc2a5d69b993d99a26e50c1f8e9723c087", "https://github.com/hongsua/UniOP"),
    "promotech": (
        "56251ad9b883ef831b4753fc623d5ec970fe65e0",
        "https://github.com/BioinformaticsLabAtMUN/PromoTech",
    ),
    "terminatornet": (
        "72180309545988dba63b8bd5db5d0b09bf4dba70",
        "https://github.com/btjaden/TerminatorNet",
    ),
    "bactermfinder": (
        "024870447f6288ebc9d29b314e957a6dfdfeb5f6",
        "https://github.com/BioinformaticsLabAtMUN/BacTermFinder",
    ),
}


def _git_tool(name: str, licence: str, module: Module) -> ToolSpec:
    sha, repo = _GIT[name]
    return _tool(name, f"git:{sha}", licence, _git_commit(name), module, repo)


_R_FGSEA_VERSION = ("R", "--vanilla", "--slave", "-e", "cat(as.character(packageVersion('fgsea')))")

TOOLS: Mapping[str, ToolSpec] = MappingProxyType(
    {
        t.name: t
        for t in (
            _tool("bakta", "1.12.2", "GPL-3.0", ("bakta", "--version"), "core"),
            _git_tool("uniop", "unlicensed", "core"),
            _git_tool("promotech", "GPL-3.0", "core"),
            _tool(
                "promoter_calculator",
                "1.2.4",
                "GPL-3.0",
                _pkg_version("promotercalculator"),
                "core",
            ),
            _git_tool("terminatornet", "MIT", "core"),
            _tool(
                "transtermhp",
                "2.09",
                "GPL-3.0",
                ("sh", "-c", 'ls "$CONDA_PREFIX/conda-meta" | grep "^transtermhp-"'),
                "core",
            ),
            _git_tool("bactermfinder", "GPL-3.0", "opt_in"),
            _tool("ostir", "1.1.3", "GPL-3.0", _pkg_version("ostir"), "core"),
            _tool("genomad", "1.12.0", "non-commercial", ("genomad", "--version"), "mobile"),
            _tool("antismash", "8.0.4", "AGPL-3.0", ("antismash", "--version"), "bgc"),
            _tool("miniprot", "0.18", "MIT", ("miniprot", "--version"), "mapping"),
            _tool("blast", "2.16.0", "public-domain", ("tblastn", "-version"), "mapping"),
            # R fgsea for the function step's GSEA, started by giap.io.fgsea_r (ADR-0038)
            _tool("fgsea", "1.36.2", "MIT", _R_FGSEA_VERSION, "function"),
            ToolSpec(
                name="codon_bias",
                version="0.5.0",
                licence="MIT",
                env_file="",
                env_name="",
                version_cmd=_pkg_version("codon-bias"),
                module="core",
                in_process=True,
            ),
        )
    }
)


# M10-08 (ADR-0054): envs of the A3 de novo pilot only. Never in TOOLS, so required_tools,
# `build_envs.sh all` and a run's tool versions never see them; `build_envs.sh <name>` builds them.
_TIMSCONVERT_COMMIT = "3cb96aef35bc0a15ce1b6cfe547386acc26a8b8b"  # tag v2.0.0
_PIP_GIT_COMMIT = (  # pip records the commit of a `git+https` install in direct_url.json
    "python",
    "-c",
    "import importlib.metadata as m, json; "
    "print(json.loads(m.distribution('timsconvert').read_text('direct_url.json'))"
    "['vcs_info']['commit_id'])",
)

PILOT_TOOLS: Mapping[str, ToolSpec] = MappingProxyType(
    {
        t.name: t
        for t in (
            _tool("cascadia", "0.0.7", "Apache-2.0", _pkg_version("cascadia"), "pilot"),
            _tool(
                "timsconvert", f"git:{_TIMSCONVERT_COMMIT}", "Apache-2.0", _PIP_GIT_COMMIT, "pilot"
            ),
        )
    }
)


def required_tools(cfg: GiapConfig, abundance: bool = False) -> list[ToolSpec]:
    """External tools a run with ``cfg`` needs (registry order); precomputed layers skip theirs.
    fgsea only when ``function`` runs GSEA, i.e. ``--abundance`` is given (M8-17 lead
    decision 8)."""
    ann = cfg.annotation
    pre = ann.precomputed
    term = ann.terminators
    mods = set(cfg.modules)
    no_terms = pre.terminators is None
    need = {
        "bakta": pre.bakta_json is None,
        "uniop": pre.operon_pairs is None,
        "promotech": pre.promoters is None,
        "promoter_calculator": ann.promoter_calculator.enabled and pre.tx_rates is None,
        "terminatornet": no_terms and term.tool == "terminatornet",
        "transtermhp": no_terms and (term.tool == "transtermhp" or term.cross_check_transtermhp),
        "bactermfinder": no_terms and term.bactermfinder,
        "ostir": pre.rbs is None,
        "genomad": "mobile" in mods and pre.regions is None,
        "antismash": "bgc" in mods and pre.regions is None,
        "miniprot": True,
        "blast": True,
        "fgsea": "function" in mods and abundance,
    }
    return [spec for name, spec in TOOLS.items() if need.get(name, False)]


# ---------------------------------------------------------------- env files and containers


@dataclass(frozen=True, slots=True)
class EnvFile:
    name: str
    channels: tuple[str, ...]
    conda: tuple[str, ...]
    pip: tuple[str, ...]
    stem: str = field(default="", compare=False)  # file stem; the name must be giap-<stem>


ALLOWED_CHANNELS = frozenset({"conda-forge", "bioconda"})
_CONDA_PIN = re.compile(r"^[A-Za-z0-9_.\-]+=[^=\s<>!~*]+(=[^=\s<>!~*]+)?$")
_PIP_PIN = re.compile(r"^[A-Za-z0-9_.\-\[\],]+==[^=\s<>!~*,;]+$")


def read_env_file(path: Path) -> EnvFile:
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    conda: list[str] = []
    pip: list[str] = []
    for dep in data.get("dependencies") or []:
        if isinstance(dep, dict):
            pip.extend(str(p) for p in dep.get("pip") or [])
        else:
            conda.append(str(dep))
    return EnvFile(
        str(data.get("name", "")),
        tuple(str(c) for c in data.get("channels") or []),
        tuple(conda),
        tuple(pip),
        Path(path).stem,
    )


def env_problems(env: EnvFile) -> list[str]:
    """Problems with ``env`` (§3.2, §4: exact pins, conda-forge + bioconda only); [] if clean."""
    problems = [f"channel '{c}' is not allowed" for c in env.channels if c not in ALLOWED_CHANNELS]
    if set(env.channels) != ALLOWED_CHANNELS and not problems:
        problems.append("channels must be exactly conda-forge and bioconda")
    problems += [
        f"{d}: conda dependency without exact version" for d in env.conda if not _CONDA_PIN.match(d)
    ]
    problems += [
        f"{d}: pip dependency must be pinned with ==" for d in env.pip if not _PIP_PIN.match(d)
    ]
    if env.stem and env.name != f"giap-{env.stem}":
        problems.append(f"env name {env.name!r} must be 'giap-{env.stem}'")
    return problems


_WITH_NAME = re.compile(
    r"withName\s*:\s*['\"]([^'\"]+)['\"]\s*\{[^}]*?\bcontainer\s*=\s*['\"]([^'\"]+)['\"]"
)


def read_containers_config(path: Path) -> dict[str, str]:
    """``withName: 'P' { container = '…' }`` entries of a Nextflow config; ``//`` lines skipped."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    text = "\n".join(line for line in lines if not line.lstrip().startswith("//"))
    return {m.group(1): m.group(2) for m in _WITH_NAME.finditer(text)}
