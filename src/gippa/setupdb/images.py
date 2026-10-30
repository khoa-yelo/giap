"""``gippa setup --images``: pull the pinned container images (M11-14 item 3; ADR-0059).

The pins are the bundled ``workflow/conf/containers.config`` (``withName`` lines whose process,
lower-cased, is a ``TOOLS`` name) plus ``workflow/conf/images.json`` (which adds the
process-less fgsea). Apptainer pulls ``<cache>/<simpleName>.img``, the file Nextflow and
preflight look for; docker pulls into its own store and the pulled digest must equal the pin.
Each image is recorded as an ``image:<tool>`` item of ``db_manifest.json``.
"""

from __future__ import annotations

import dataclasses
import json
import os
import re
import subprocess
import time
from collections.abc import Mapping
from pathlib import Path
from typing import Literal

from gippa.annotate.tools import TOOLS, read_containers_config, required_tools
from gippa.config import GippaConfig
from gippa.exitcodes import ExitCode
from gippa.runner.nextflow import workflow_main
from gippa.setupdb.items import DbRecord, Status, today
from gippa.setupdb.manifest import read_records, write_manifest
from gippa.setupdb.space import disk_full, is_disk_full

__all__ = [
    "Runtime",
    "failed_images",
    "image_record",
    "image_set",
    "image_sizes",
    "image_versions",
    "install_images",
    "pinned_images",
    "pulled_digest",
    "simple_name",
]

Runtime = Literal["apptainer", "docker"]
_UPSTREAM_TAG = re.compile(r"//\s+([a-z0-9./_-]+):([A-Za-z0-9_.-]+)\s*")
_SIF_GIB = re.compile(r"^//\s+sif_gib((?:\s+[a-z0-9_]+=[0-9.]+)+)\s*$", re.MULTILINE)
_DIGEST_LINE = re.compile(r"^Digest:\s*(sha256:[0-9a-f]{64})\s*$", re.MULTILINE)
PULL_ATTEMPTS = 2  # one retry: a pull can fail transiently; it helps most when layers are cached
RETRY_WAIT_S = 5.0
TAIL_LINES = 5  # stderr lines kept in a failed record's message (#504)
NO_CACHE_HINT = "; set APPTAINER_CACHEDIR to keep the layers between attempts"


def _conf(conf_dir: Path | None) -> Path:
    return conf_dir if conf_dir is not None else workflow_main().parent / "conf"


def _images_json(conf_dir: Path) -> dict[str, dict[str, object]]:
    """``images.json`` entries by tool; ``{}`` when the file is missing or not an object."""
    try:
        data = json.loads((conf_dir / "images.json").read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return {}
    if not isinstance(data, dict):
        return {}
    return {str(k): v for k, v in data.items() if isinstance(v, dict)}


def simple_name(ref: str) -> str:
    """Nextflow's apptainer cache file name of ``ref`` (ADR-0059)."""
    return ref.removeprefix("docker://").replace("/", "-").replace(":", "-") + ".img"


def pinned_images(conf_dir: Path | None = None) -> dict[str, str]:
    """Registry tool → bare pinned ref (``<image>@sha256:<d>``). ``containers.config`` wins over
    ``images.json``, which adds process-less tools (fgsea)."""
    conf = _conf(conf_dir)
    out: dict[str, str] = {}
    for tool, entry in _images_json(conf).items():
        ref = entry.get("ref")
        if tool in TOOLS and isinstance(ref, str) and ref:
            out[tool] = ref.removeprefix("docker://")
    path = conf / "containers.config"
    if path.is_file():
        for process, ref in read_containers_config(path).items():
            if process.lower() in TOOLS:
                out[process.lower()] = ref.removeprefix("docker://")
    return dict(sorted(out.items()))


def image_versions(conf_dir: Path | None = None) -> dict[str, str]:
    """Tool → recorded version: the ``images.json`` tag, else the upstream tag the digest was
    resolved from (the ``//   <image>:<tag>`` comment lines of ``containers.config``)."""
    conf = _conf(conf_dir)
    path = conf / "containers.config"
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    upstream = {
        m.group(1): m.group(2)
        for line in text.splitlines()
        if (m := _UPSTREAM_TAG.fullmatch(line.strip()))
    }
    tags = {t: e.get("tag") for t, e in _images_json(conf).items()}
    out: dict[str, str] = {}
    for tool, ref in pinned_images(conf).items():
        tag = tags.get(tool)
        out[tool] = tag if isinstance(tag, str) and tag else upstream.get(ref.split("@", 1)[0], "")
    return out


def image_sizes(conf_dir: Path | None = None) -> dict[str, float]:
    """Tool -> SIF size in GiB from the ``//   sif_gib tool=x.y ...`` lines of
    ``containers.config`` (measured; for the setup plan)."""
    path = _conf(conf_dir) / "containers.config"
    text = path.read_text(encoding="utf-8") if path.is_file() else ""
    pairs = (kv.split("=") for m in _SIF_GIB.finditer(text) for kv in m.group(1).split())
    return {k: float(v) for k, v in pairs}


def image_set(cfg: GippaConfig, all_images: bool, conf_dir: Path | None = None) -> list[str]:
    """Sorted tool names to pull: every pin, or the pinned tools ``cfg`` needs (with GSEA)."""
    pins = pinned_images(conf_dir)
    if all_images:
        return sorted(pins)
    return sorted(t.name for t in required_tools(cfg, abundance=True) if t.name in pins)


def pulled_digest(output: str) -> str | None:
    """The ``Digest: sha256:<hex>`` line ``docker pull`` prints; ``None`` without one."""
    m = _DIGEST_LINE.search(output)
    return m.group(1) if m else None


def image_record(
    tool: str, ref: str, version: str, status: Status, message: str = "", path: Path | None = None
) -> DbRecord:
    """The ``image:<tool>`` manifest item; ``checksum`` is the pinned digest."""
    ok = status == "ok"
    size = path.stat().st_size if ok and path is not None and path.is_file() else 0
    return DbRecord(f"image:{tool}", version, ref, ref.split("@", 1)[1] if "@" in ref else "",
                    size, today() if ok else "", str(path) if path is not None else "", status,
                    message)  # fmt: skip


def _tail(text: str, n: int = TAIL_LINES) -> str:
    """The last ``n`` non-empty lines of ``text``, joined into one line."""
    lines = [ln.strip() for ln in text.splitlines() if ln.strip()]
    return " | ".join(lines[-n:]) if lines else "no output"


def _pull(runtime: Runtime, tool: str, ref: str, version: str, cache: Path) -> DbRecord:
    """Pull one image (apptainer: into ``cache``, skipped when already there). A failed pull is
    retried once (#504), except a full disk or quota, a digest mismatch or a missing runtime. An
    image another ``gippa setup`` finished during the wait is kept, not pulled again."""
    if runtime == "apptainer":
        dest = cache / simple_name(ref)
        if dest.is_file():
            return image_record(tool, ref, version, "ok", "", dest)
        cache.mkdir(parents=True, exist_ok=True)
        # the layer blobs would also land in ~/.apptainer/cache unless the user chose a cache
        own = ["--disable-cache"] if not os.environ.get("APPTAINER_CACHEDIR") else []
        cmd = ["apptainer", "pull", *own, str(dest), f"docker://{ref}"]
    else:
        dest = None
        cmd = ["docker", "pull", ref]
    rec = image_record(tool, ref, version, "absent", f"{runtime} pull {ref} was not run")
    for attempt in range(PULL_ATTEMPTS):
        if attempt:
            time.sleep(RETRY_WAIT_S)
            if dest is not None and dest.is_file():
                return image_record(tool, ref, version, "ok", "", dest)
        rec, retry = _attempt(runtime, cmd, dest, tool, ref, version)
        if not retry:
            break
    return rec


def _attempt(
    runtime: Runtime, cmd: list[str], dest: Path | None, tool: str, ref: str, version: str
) -> tuple[DbRecord, bool]:
    """One run of the pull: its record and whether a failure is worth one more run."""
    hint = NO_CACHE_HINT if "--disable-cache" in cmd else ""
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except OSError as exc:
        return image_record(tool, ref, version, "absent", f"{runtime} pull failed: {exc}"), False
    if proc.returncode != 0:
        if dest is not None:
            dest.unlink(missing_ok=True)  # no partial image is left for Nextflow (or a retry)
        if is_disk_full(f"{proc.stderr}\n{proc.stdout}"):
            where = dest or "the docker storage"
            full = disk_full(where, " or --cache-dir")
            return image_record(tool, ref, version, "failed", full), False
        why = _tail(proc.stderr or proc.stdout)
        msg = f"{runtime} pull {ref} failed: {why}{hint}"
        return image_record(tool, ref, version, "absent", msg), True
    if runtime == "docker":
        pin = ref.split("@", 1)[1]
        got = pulled_digest(proc.stdout + "\n" + proc.stderr)
        if got != pin:
            return image_record(tool, ref, version, "mismatch",
                                f"pulled {got or 'no digest'}, pinned {pin}"), False  # fmt: skip
        return image_record(tool, ref, version, "ok", f"pulled {ref}"), False
    if dest is None or not dest.is_file():
        msg = f"apptainer pull wrote no {dest}{hint}"
        return image_record(tool, ref, version, "absent", msg), True
    return image_record(tool, ref, version, "ok", f"pulled {ref}", dest), False


def failed_images(tools: list[str], db_dir: Path) -> tuple[list[str], list[str]]:
    """The ``tools`` not ``ok`` in ``db_dir``, for the summary at the end of ``gippa setup`` (#504):
    those a rerun can help, then those it cannot (a full disk, a digest mismatch)."""
    recs = read_records(db_dir) if tools else {}
    status = {t: r.status if (r := recs.get(f"image:{t}")) is not None else "absent" for t in tools}
    return ([t for t, s in status.items() if s == "absent"],
            [t for t, s in status.items() if s in ("failed", "mismatch")])  # fmt: skip


def _kept(rec: DbRecord, old: Mapping[str, DbRecord]) -> DbRecord:
    """An unchanged ``ok`` image keeps its previous record (idempotent reruns)."""
    prev = old.get(rec.name)
    if rec.status == "ok" and prev is not None and prev.status == "ok" and (
        prev.version, prev.url, prev.checksum, prev.path) == (
        rec.version, rec.url, rec.checksum, rec.path):  # fmt: skip
        return prev
    return dataclasses.replace(rec, message="") if rec.status == "ok" else rec


def install_images(
    tools: list[str], runtime: Runtime, cache: Path, db_dir: Path, conf_dir: Path | None = None
) -> ExitCode:
    """Pull ``tools``' pinned images (one pull per distinct ref) and record them in
    ``db_dir/db_manifest.json``. A failed pull or a digest mismatch is exit 2; the others are
    still pulled."""
    from gippa.setupdb.run import LOCK, _locked, _report, _say

    if not tools:
        return ExitCode.SUCCESS
    pins = pinned_images(conf_dir)
    versions = image_versions(conf_dir)
    db_dir.mkdir(parents=True, exist_ok=True)
    with _locked(db_dir) as acquired:
        if not acquired:
            _say(f"another gippa setup is running on {db_dir} ({LOCK}); try again later")
            return ExitCode.ENVIRONMENT
        old = read_records(db_dir)
        done: dict[str, DbRecord] = {}
        recs: list[DbRecord] = []
        for tool in tools:
            ref = pins[tool]
            first = done.get(ref)
            if first is None:
                rec = _pull(runtime, tool, ref, versions.get(tool, ""), cache)
                done[ref] = rec
            else:
                rec = dataclasses.replace(first, name=f"image:{tool}",
                                          version=versions.get(tool, ""))  # fmt: skip
            _report(rec)
            recs.append(_kept(rec, old))
        write_manifest(db_dir, recs)
    return ExitCode.SUCCESS if all(r.status == "ok" for r in recs) else ExitCode.ENVIRONMENT
