"""Render every registered plot of one ``OUT/`` tree in a fresh interpreter, for the M9-26
byte-identity tests (docs/specs/M9-26-scale-and-determinism.md work item 4, example 2; M9-01
item 7; GIPPA_SPEC.md R7/TE7).

``render_in_subprocess`` runs ``python -m helpers.byte_identity OUT DEST`` with a chosen
``PYTHONHASHSEED``: two renders in two interpreters with different hash seeds catch any output
that depends on ``set``/``dict`` iteration of strings or on process state, which two calls in one
process would not. Each plot goes through its registry wrapper (declared skips skip), after
``gippa.viz.style.apply()`` as the ``plots`` step does, and every instance is written with
``gippa.io.plot_files.write_instance`` in ``svg``, ``pdf`` and ``png``.

``golden_out`` builds a small rich ``OUT/`` from the committed gippasim sim50k golden tables
(``tests/data/golden/sim50k``) plus a one-replicon ``preflight.json`` (``chr1``, 50 kb,
circular), so many plots draw something; nothing touches ``$GIPPA_DB_DIR`` or the network.
"""

from __future__ import annotations

import difflib
import json
import os
import subprocess
import sys
from collections.abc import Sequence
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
GOLDEN = REPO / "tests" / "data" / "golden" / "sim50k"
FORMATS = ("svg", "pdf", "png")
#: golden file -> its ``OUT/``-relative path (``gippa.io.run_view.PATHS``)
GOLDEN_FILES = {
    "m4/genes.tsv": "tables/genes.tsv",
    "m4/detection_summary.tsv": "tables/detection_summary.tsv",
    "m4/regions_stats.tsv": "tables/regions_stats.tsv",
    "m4/sensitivity.tsv": "tables/sensitivity.tsv",
    "m5/autocorrelation.tsv": "tables/autocorrelation.tsv",
    "m5/islands.tsv": "tables/islands.tsv",
    "mapping/mapping_fdr.tsv": "mapping/mapping_fdr.tsv",
    "mapping/peptide_hits.tsv": "mapping/peptide_hits.tsv",
    "mapping/protein_hits.tsv": "mapping/protein_hits.tsv",
    "mapping/unmapped.tsv": "mapping/unmapped.tsv",
}
SIM50K_LENGTH = 50_000


def golden_out(tmp: Path) -> Path:
    """An ``OUT/`` of the sim50k golden tables at their published paths, with a preflight."""
    from helpers.viz_landscape import write_preflight

    out = tmp / "OUT"
    for src, dst in GOLDEN_FILES.items():
        path = out / dst
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes((GOLDEN / src).read_bytes())
    write_preflight(out, [("chr1", SIM50K_LENGTH, "circular", "chromosome")])
    return out


def render_all(out: Path, dest: Path, pids: Sequence[str] = (), dpi: int = 100) -> dict[str, str]:
    """Render every registered plot (or only ``pids``) of ``out`` into ``dest``; returns
    ``{pid: "rendered" | "skip:<reason>"}``."""
    from gippa.config import resolve_config
    from gippa.io.plot_files import write_instance
    from gippa.io.run_view import RunView
    from gippa.steps.plots import load_plot_modules
    from gippa.viz import base, style

    load_plot_modules()
    style.apply()
    runs = (out / "resolved_config.yaml", out / "work" / "resolved_config.yaml")
    cfg = resolve_config(next((p for p in runs if p.is_file()), None))
    view = RunView.from_outdir(out)
    status: dict[str, str] = {}
    for pid, spec in base.registered().items():
        if pids and pid not in pids:
            continue
        result = spec.fn(view, cfg)
        if isinstance(result, base.Skip):
            status[pid] = f"skip:{result.reason}"
            continue
        for inst in result:
            write_instance(inst, pid, spec.name, dest, list(FORMATS), dpi)
        status[pid] = "rendered"
    return status


def render_in_subprocess(
    out: Path, dest: Path, hashseed: int, pids: Sequence[str] = (), timeout: float = 1800
) -> dict[str, str]:
    """``render_all`` in a fresh interpreter with ``PYTHONHASHSEED=hashseed``."""
    env = {**os.environ, "PYTHONHASHSEED": str(hashseed), "MPLBACKEND": "Agg"}
    paths = [str(REPO / "src"), str(REPO / "tests" / "sim"), str(REPO / "tests")]
    env["PYTHONPATH"] = os.pathsep.join([*paths, env.get("PYTHONPATH", "")])
    status = dest.parent / f"{dest.name}.status.json"
    cmd = [sys.executable, "-m", "helpers.byte_identity", str(out), str(dest), str(status)]
    cmd += list(pids)
    res = subprocess.run(cmd, capture_output=True, text=True, env=env, timeout=timeout)
    assert res.returncode == 0, f"{' '.join(cmd)}\n{res.stdout}\n{res.stderr}"
    data: dict[str, str] = json.loads(status.read_text(encoding="utf-8"))
    return data


def tree_bytes(root: Path) -> dict[str, bytes]:
    """Every file under ``root`` (``root``-relative POSIX path -> bytes)."""
    files = (p for p in sorted(root.rglob("*")) if p.is_file())
    return {p.relative_to(root).as_posix(): p.read_bytes() for p in files}


def differing(a: dict[str, bytes], b: dict[str, bytes]) -> list[str]:
    """Paths present in only one tree or with different bytes."""
    return sorted(k for k in a.keys() | b.keys() if a.get(k) != b.get(k))


def svg_diff_report(a: dict[str, bytes], b: dict[str, bytes], max_lines: int = 40) -> str:
    """A failure message: the differing paths, then for each differing ``.svg`` the first
    ``max_lines`` lines of a ``difflib.unified_diff`` of its two renders (e.g. a clip-path id
    that changed with a 1-ulp layout difference)."""
    paths = differing(a, b)
    parts = [f"differing: {paths}"]
    for path in paths:
        if not path.endswith(".svg") or path not in a or path not in b:
            continue
        la = a[path].decode("utf-8", "replace").splitlines()
        lb = b[path].decode("utf-8", "replace").splitlines()
        diff = list(difflib.unified_diff(la, lb, f"a/{path}", f"b/{path}", n=1, lineterm=""))
        parts.append("\n".join(diff[:max_lines]))
        if len(diff) > max_lines:
            parts.append(f"... ({len(diff) - max_lines} more diff lines)")
    return "\n".join(parts)


def _main(argv: Sequence[str]) -> int:
    out, dest, status, *pids = argv
    result = render_all(Path(out), Path(dest), pids)
    Path(status).write_text(json.dumps(result, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(_main(sys.argv[1:]))
