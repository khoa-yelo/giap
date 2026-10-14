#!/usr/bin/env python3
"""Probe the FarmShare environment and write .giap/capabilities.json (GIAP_SPEC.md §3.6).

Standard library only, so it runs before the project environment exists.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shutil
import socket
import subprocess
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime
from pathlib import Path

REPO = Path(os.environ.get("GIAP_REPO", Path(__file__).resolve().parents[1]))
ROOT = Path(os.environ.get("GIAP_ROOT", REPO / "local"))
DISK_BUDGET_GB = 400

# §3.2: one URL per pinned download source
SOURCES = {
    "pypi.org": "https://pypi.org/simple/",
    "conda.anaconda.org": "https://conda.anaconda.org/conda-forge/noarch/repodata.json.zst",
    "github.com": "https://github.com",
    "ghcr.io": "https://ghcr.io/v2/",
    "registry-1.docker.io": "https://registry-1.docker.io/v2/",
    "zenodo.org": "https://zenodo.org/records/14916843",
    "ftp.ncbi.nlm.nih.gov": "https://ftp.ncbi.nlm.nih.gov/",
    "www.cs.mun.ca": "https://www.cs.mun.ca/~lourdes/public/PromoTech_models/",
    "dl.secondarymetabolites.org": "https://dl.secondarymetabolites.org/",
    "rest.uniprot.org": "https://rest.uniprot.org/proteomes/UP000000625.json",
    "purl.obolibrary.org": "https://purl.obolibrary.org/obo/go/go-basic.obo",
    "rest.kegg.jp": "https://rest.kegg.jp/info/kegg",
    "regulondb.ccg.unam.mx": "https://regulondb.ccg.unam.mx/",
    "proteomecentral.proteomexchange.org": "https://proteomecentral.proteomexchange.org/",
}

TOOLS = {
    "python": ["python3", "--version"],
    "java": ["java", "-version"],
    "nextflow": ["nextflow", "-version"],
    "nf-test": ["nf-test", "version"],
    "micromamba": ["micromamba", "--version"],
    "uv": ["uv", "--version"],
    "gh": ["gh", "--version"],
    "apptainer": ["apptainer", "--version"],
}


def run(cmd: list[str], timeout: float = 30) -> tuple[int, str]:
    if shutil.which(cmd[0]) is None:
        return 127, ""
    try:
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
    except (subprocess.TimeoutExpired, OSError) as exc:
        return 124, str(exc)
    return p.returncode, (p.stdout + p.stderr).strip()


def first_version_line(text: str) -> str:
    for line in text.splitlines():
        line = line.strip(" >")
        if any(ch.isdigit() for ch in line):
            return line
    return text.splitlines()[0] if text else ""


def reachable(url: str, timeout: float = 5) -> dict[str, object]:
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "giap-probe"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return {"ok": True, "status": r.status}
    except urllib.error.HTTPError as exc:
        # Any HTTP answer means the host is reachable; 405/403/401 on HEAD are common
        return {"ok": exc.code < 500, "status": exc.code}
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        return {"ok": False, "error": str(exc)[:120]}


def du_gb(path: Path) -> float:
    code, out = run(["du", "-sk", str(path)], timeout=120)
    if code != 0 or not out:
        return 0.0
    return round(int(out.split()[0]) / 1024**2, 2)


def allocation() -> dict[str, object]:
    job = os.environ.get("SLURM_JOB_ID")
    info: dict[str, object] = {
        "slurm_job_id": job,
        "partition": os.environ.get("SLURM_JOB_PARTITION"),
        "cpus": int(os.environ.get("SLURM_CPUS_ON_NODE", os.cpu_count() or 1)),
    }
    mem_per_cpu = os.environ.get("SLURM_MEM_PER_CPU")
    mem_node = os.environ.get("SLURM_MEM_PER_NODE")
    if mem_node:
        info["mem_gb"] = round(int(mem_node) / 1024, 1)
    elif mem_per_cpu:
        info["mem_gb"] = round(int(mem_per_cpu) * int(info["cpus"]) / 1024, 1)  # type: ignore[arg-type]
    if job:
        _, left = run(["squeue", "-h", "-j", job, "-o", "%L"])
        info["time_left"] = left or None
    return info


def slurm() -> dict[str, object]:
    code, _ = run(
        ["sbatch", "--test-only", "--qos=dev", "-p", "normal", "-t", "1", "--wrap", "true"]
    )
    _, queue = run(["squeue", "--me", "-h", "-o", "%i %j %T %M %l"])
    giap_jobs = [line for line in queue.splitlines() if " giap-" in f" {line.split(' ', 1)[-1]}"]
    _, parts = run(["sinfo", "-h", "-o", "%P"])
    return {
        "sbatch_ok": code == 0,
        "gpu_partition": "gpu" in parts.split(),
        "giap_jobs": giap_jobs,
    }


def databases() -> dict[str, object]:
    db = Path(os.environ.get("GIAP_DB_DIR", ROOT / "db"))
    tools = Path(os.environ.get("GIAP_TOOLS_DIR", ROOT / "tools"))
    manifest = db / "db_manifest.json"
    out: dict[str, object] = {
        "db_dir": str(db),
        "bakta_light": (db / "bakta" / "db-light").is_dir(),
        "bakta_full": (db / "bakta" / "db").is_dir(),
        "genomad": (db / "genomad_db").is_dir(),
        "promotech_rf_hot": any(db.glob("promotech/**/RF-HOT.model")),
        "diann": any(tools.glob("diann-2.7.0*/**/diann*")),
    }
    if manifest.is_file():
        out["manifest_sha256"] = hashlib.sha256(manifest.read_bytes()).hexdigest()
    return out


def probe(check_net: bool) -> dict[str, object]:
    tools = {}
    for name, cmd in TOOLS.items():
        code, out = run(cmd)
        tools[name] = first_version_line(out) if code == 0 else None
    root_usage = {}
    if ROOT.is_dir():
        for sub in sorted(p for p in ROOT.iterdir() if p.is_dir()):
            root_usage[sub.name] = du_gb(sub)
    used = round(sum(root_usage.values()), 2)
    caps: dict[str, object] = {
        "probed_at": datetime.now(UTC).isoformat(timespec="seconds"),
        "host": socket.gethostname(),
        "os": platform.platform(),
        "allocation": allocation(),
        "giap_root": str(ROOT),
        "giap_root_used_gb": used,
        "giap_root_usage_gb": root_usage,
        "giap_root_over_budget": used > DISK_BUDGET_GB,
        "giap_root_free_gb": round(shutil.disk_usage(ROOT if ROOT.exists() else REPO).free / 1e9),
        "tmp_free_gb": round(shutil.disk_usage("/tmp").free / 1e9),
        "tools": tools,
        "slurm": slurm(),
        "databases": databases(),
    }
    if check_net:
        caps["sources"] = {host: reachable(url) for host, url in SOURCES.items()}
    return caps


def summary(caps: dict[str, object]) -> str:
    a = caps["allocation"]
    s = caps["slurm"]
    lines = [
        f"GIAP probe on {caps['host']} at {caps['probed_at']}",
        f"  allocation: job={a.get('slurm_job_id')} cpus={a.get('cpus')} "  # type: ignore[union-attr]
        f"mem={a.get('mem_gb')}GB time_left={a.get('time_left')}",  # type: ignore[union-attr]
        f"  local/: used={caps['giap_root_used_gb']} GB (budget {DISK_BUDGET_GB}) "
        f"free={caps['giap_root_free_gb']} GB; /tmp free={caps['tmp_free_gb']} GB",
        "  tools: "
        + ", ".join(f"{k}={'yes' if v else 'MISSING'}" for k, v in caps["tools"].items()),  # type: ignore[union-attr]
        f"  slurm: sbatch={s['sbatch_ok']} gpu_partition={s['gpu_partition']} "  # type: ignore[index]
        f"giap_jobs={len(s['giap_jobs'])}",  # type: ignore[index]
        "  databases: "
        + ", ".join(f"{k}={v}" for k, v in caps["databases"].items() if k != "db_dir"),  # type: ignore[union-attr]
    ]
    if "sources" in caps:
        bad = [h for h, r in caps["sources"].items() if not r["ok"]]  # type: ignore[union-attr]
        lines.append("  sources: all reachable" if not bad else f"  sources UNREACHABLE: {bad}")
    if not a.get("slurm_job_id"):  # type: ignore[union-attr]
        lines.append("  WARNING: not inside a Slurm allocation (R14)")
    if caps["giap_root_over_budget"]:
        lines.append("  WARNING: local/ is over the disk budget (§3.9)")
    return "\n".join(lines)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=REPO / ".giap" / "capabilities.json")
    ap.add_argument("--quiet", action="store_true", help="skip network checks; one-screen summary")
    ap.add_argument("--json", action="store_true", help="print JSON instead of the summary")
    args = ap.parse_args()
    caps = probe(check_net=not args.quiet)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(caps, indent=2, sort_keys=True) + "\n")
    print(json.dumps(caps, indent=2) if args.json else summary(caps))
    return 0


if __name__ == "__main__":
    sys.exit(main())
