"""Fake tool installs and a real-format ``giap setup`` manifest for preflight tests.

Run-wiring lead decision 4 (docs/specs/M2-07-nextflow-tool-processes.md): ``environment.tools``
takes binaries/environments from the M2-01 tool registry (``giap.annotate.tools.TOOLS``) and
databases/models from the ``giap setup`` manifest. Tests never read the real ``$GIAP_DB_DIR`` or
``$GIAP_ROOT``; they build fakes with the helpers here:

- :func:`write_setup_manifest` copies ``tests/data/db_manifest/db_manifest.json`` (the structure of
  a real ``giap setup`` run, versions/URLs/paths faked) into a directory, optionally dropping
  items or changing their status.
- :func:`fake_giap_root` lays out ``<root>/containers/<image>.img`` for every tool pinned in
  ``workflow/conf/containers.config`` (file name as Nextflow's SingularityCache/
  ``scripts/slurm/build_envs.sbatch`` derive it: scheme stripped, ``/`` and ``:`` become ``-``) and
  ``<root>/micromamba/envs/giap-<tool>/.giap-built`` holding the sha256 of
  ``workflow/envs/<tool>.yml`` for every conda-only tool (``scripts/build_envs.sh``).
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable, Mapping
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
CONTAINERS_CONFIG = REPO / "workflow" / "conf" / "containers.config"
ENVS = REPO / "workflow" / "envs"
SETUP_MANIFEST = REPO / "tests" / "data" / "db_manifest" / "db_manifest.json"
MARKER = ".giap-built"


def image_file_name(ref: str) -> str:
    """``docker://a/b@sha256:<d>`` → ``a-b@sha256-<d>.img`` (build_envs.sbatch / Nextflow cache)."""
    return ref.removeprefix("docker://").replace("/", "-").replace(":", "-") + ".img"


def pinned_images() -> dict[str, str]:
    """Tool name (registry key) → pinned image reference, from containers.config."""
    from giap.annotate.tools import TOOLS, read_containers_config

    out: dict[str, str] = {}
    for process, ref in read_containers_config(CONTAINERS_CONFIG).items():
        name = process.lower()
        if name in TOOLS:  # TBLASTN reuses the blast image; the registry key is "blast"
            out[name] = ref
    return out


def conda_only_tools() -> list[str]:
    """Registry tools with an env file and no pinned container (prebuilt micromamba env)."""
    from giap.annotate.tools import TOOLS

    pinned = pinned_images()
    return sorted(n for n, t in TOOLS.items() if not t.in_process and n not in pinned)


def env_yml_sha256(tool: str) -> str:
    return hashlib.sha256((ENVS / f"{tool}.yml").read_bytes()).hexdigest()


def fake_giap_root(
    root: Path,
    *,
    skip: Iterable[str] = (),
    stale: Iterable[str] = (),
    no_marker: Iterable[str] = (),
) -> Path:
    """A fake ``$GIAP_ROOT`` with every pinned image and every conda-only env built.

    ``skip``: tools left out entirely; ``stale``: env marker holds a wrong hash; ``no_marker``:
    env directory present without ``.giap-built``.
    """
    skip, stale, no_marker = set(skip), set(stale), set(no_marker)
    (root / "containers").mkdir(parents=True, exist_ok=True)
    (root / "micromamba" / "envs").mkdir(parents=True, exist_ok=True)
    for tool, ref in pinned_images().items():
        if tool not in skip:
            (root / "containers" / image_file_name(ref)).write_bytes(b"fake sif image\n")
    for tool in conda_only_tools():
        if tool in skip:
            continue
        prefix = root / "micromamba" / "envs" / f"giap-{tool}"
        (prefix / "bin").mkdir(parents=True, exist_ok=True)
        if tool in no_marker:
            continue
        want = "0" * 64 if tool in stale else env_yml_sha256(tool)
        (prefix / MARKER).write_text(want + "\n")
    return root


def setup_manifest() -> dict[str, object]:
    data = json.loads(SETUP_MANIFEST.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    return data


def write_setup_manifest(
    db_dir: Path, *, drop: Iterable[str] = (), status: Mapping[str, str] | None = None
) -> Path:
    """Copy the real-format fixture to ``db_dir/db_manifest.json``; drop or re-status items."""
    data = setup_manifest()
    drop, status = set(drop), dict(status or {})
    items = [dict(i) for i in data["items"] if i["name"] not in drop]  # type: ignore[attr-defined]
    for item in items:
        if item["name"] in status:
            item["status"] = status[item["name"]]
    db_dir.mkdir(parents=True, exist_ok=True)
    (db_dir / "db_manifest.json").write_text(
        json.dumps({**data, "items": items}, indent=2), encoding="utf-8"
    )
    return db_dir
