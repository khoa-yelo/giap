"""Test data for the R fgsea engine wrapper (M8-17a; docs/specs/M8-17-r-fgsea-engine.md).

- :func:`load_call` rebuilds the ``EngineCall`` of a recorded directory
  (``tests/data/recorded/fgsea/1.36.2/<input>/``) from its ``stats.tsv``, ``jobs.tsv``,
  ``sets.tsv`` and the ``--seed`` of ``command.txt``.
- :func:`fake_prefix` lays out a fake micromamba prefix whose ``bin/`` holds a fake R front end.
  The wrapper starts ``workflow/bin/giap_fgsea.R`` whose shebang (``/usr/bin/env -S <R front end>
  --vanilla``) finds that fake first on the child ``PATH``; the fake is a small Python program
  driven by a JSON behaviour file. It logs its argv, cwd, environment and the input files it saw
  to ``<prefix>/log.json``. CI never runs R.
"""

from __future__ import annotations

import json
import os
import re
import stat
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parents[2]
RECORDED = REPO / "tests" / "data" / "recorded" / "fgsea" / "1.36.2"
SCRIPT = REPO / "workflow" / "bin" / "giap_fgsea.R"
GOLDEN = REPO / "tests" / "data" / "golden" / "r"
RESULTS_HEADER = "job_id\tpathway\tsize\tES\tNES\tpval\tlog2err\tleadingEdge"
# R's script front end, spelled in pieces only so that grep-based audits of the tree stay simple.
R_FRONT_END = "R" + "script"


def _rows(path: Path) -> list[list[str]]:
    lines = path.read_text(encoding="utf-8").split("\n")
    assert lines[-1] == "", f"{path} must end with a newline"
    return [line.split("\t") for line in lines[1:-1]]


def load_call(d: Path) -> Any:
    """The ``EngineCall`` whose ``write_call`` produced the recorded inputs in ``d``."""
    from giap.config.model import GseaConfig
    from giap.function.gsea import EngineCall, GseaJob

    stats = _rows(d / "stats.tsv")
    sets: dict[str, dict[str, list[str]]] = {}
    for job_id, set_id, member in _rows(d / "sets.tsv"):
        sets.setdefault(job_id, {}).setdefault(set_id, []).append(member)
    jobs = tuple(
        GseaJob(
            job_id=job_id,
            sets={k: tuple(v) for k, v in sets.get(job_id, {}).items()},
            min_size=int(lo),
            max_size=int(hi),
        )
        for job_id, lo, hi in _rows(d / "jobs.tsv")
    )
    command = (d / "command.txt").read_text(encoding="utf-8")
    seed = re.search(r"--seed (\d+)", command)
    assert seed is not None
    return EngineCall(
        tu_ids=tuple(r[0] for r in stats),
        scores=np.array([float(r[1]) for r in stats], dtype=np.float64),
        jobs=jobs,
        params=GseaConfig(),
        seed=int(seed.group(1)),
    )


_FAKE = r'''#!{python} -I
"""Fake R front end for giap_fgsea.R (tests only)."""
import json, os, subprocess, sys, time
from pathlib import Path

PREFIX = Path({prefix!r})
cfg = json.loads((PREFIX / "behaviour.json").read_text())
argv = sys.argv[1:]
opts = dict(zip(argv[2::2], argv[3::2]))
cwd = Path.cwd()
seen = {{p.name: p.read_text() for p in sorted(cwd.iterdir()) if p.is_file()}}
log = {{"argv": argv, "cwd": str(cwd), "env": dict(os.environ), "files": seen}}
(PREFIX / "log.json").write_text(json.dumps(log))
with (PREFIX / "calls.txt").open("a") as fh:
    fh.write(str(cwd) + "\n")
mode = cfg.get("mode", "ok")
if mode == "fail":
    sys.stderr.write(cfg.get("stderr", ""))
    sys.exit(cfg.get("code", 1))
if mode == "sleep":
    child = subprocess.Popen(["sleep", str(cfg.get("seconds", 20))])
    (PREFIX / "grandchild.pid").write_text(str(child.pid))
    time.sleep(cfg.get("seconds", 20))
    sys.exit(0)
if mode == "no_output":
    sys.exit(0)
jobs = [l.split("\t") for l in (cwd / opts["--jobs"]).read_text().splitlines()[1:]]
sets = {{}}
for line in (cwd / opts["--sets"]).read_text().splitlines()[1:]:
    job, set_id, member = line.split("\t")
    sets.setdefault(job, {{}}).setdefault(set_id, []).append(member)
if "results" in cfg:
    text = cfg["results"]
else:
    lines = [{header!r}]
    for job, lo, hi in jobs:
        for set_id, members in sets.get(job, {{}}).items():
            if int(lo) <= len(members) <= int(hi):
                size = len(members) + cfg.get("size_offset", 0)
                edge = ";".join(members[:2])
                lines.append(f"{{job}}\t{{set_id}}\t{{size}}\t0.5\t1.5\t0.01\t0.25\t{{edge}}")
    text = "\n".join(lines) + "\n"
(cwd / opts["--out"]).write_text(text)
versions = [
    "R\tR version 4.5.3 (2026-03-11)",
    "fgsea\t" + cfg.get("fgsea", "1.36.2"),
    "data.table\t1.18.6.1",
    "BiocParallel\t1.44.0",
    "seed\t" + opts.get("--seed", ""),
    "jobs\t" + str(len(jobs)),
] + ["warning\t" + j + "\t" + m for j, m in cfg.get("warnings", [])]
(cwd / opts["--versions"]).write_text("\n".join(versions) + "\n")
sys.exit(0)
'''


def fake_prefix(root: Path, behaviour: Mapping[str, Any] | None = None) -> Path:
    """``root`` as a fake env prefix: ``bin/<R front end>`` plus ``behaviour.json``."""
    (root / "bin").mkdir(parents=True, exist_ok=True)
    exe = root / "bin" / R_FRONT_END
    exe.write_text(
        _FAKE.format(python=sys.executable, prefix=str(root), header=RESULTS_HEADER),
        encoding="utf-8",
    )
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    set_behaviour(root, behaviour or {})
    return root


def set_behaviour(prefix: Path, behaviour: Mapping[str, Any]) -> None:
    (prefix / "behaviour.json").write_text(json.dumps(dict(behaviour)), encoding="utf-8")


def fake_log(prefix: Path) -> dict[str, Any]:
    return json.loads((prefix / "log.json").read_text(encoding="utf-8"))


def options(argv: list[str]) -> dict[str, str]:
    """``--key value`` pairs after the two leading arguments (``--vanilla``, the script)."""
    rest = argv[2:]
    assert len(rest) % 2 == 0, rest
    keys = rest[0::2]
    assert all(k.startswith("--") for k in keys), rest
    assert len(set(keys)) == len(keys), rest
    return dict(zip(keys, rest[1::2], strict=True))


def real_prefix() -> Path | None:
    """The built ``giap-fgsea`` env (``$GIAP_FGSEA_PREFIX``, else ``$GIAP_ROOT``'s), or None."""
    env = os.environ.get("GIAP_FGSEA_PREFIX")
    candidates = [Path(env)] if env else []
    root = os.environ.get("GIAP_ROOT")
    if root:
        candidates.append(Path(root) / "micromamba" / "envs" / "giap-fgsea")
    for c in candidates:
        if (c / "bin" / R_FRONT_END).is_file():
            return c
    return None
