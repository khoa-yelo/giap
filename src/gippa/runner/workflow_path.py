"""Where the bundled Nextflow workflow lives; import-light, so ``gippa.cli`` can read the profile
names at load time without pulling in the step registry (cold ``gippa _step`` start)."""

from __future__ import annotations

from pathlib import Path


def workflow_main() -> Path:
    """``workflow/main.nf`` of this checkout (src layout) or of an installed package."""
    here = Path(__file__).resolve()
    candidates = [
        here.parents[3] / "workflow" / "main.nf",
        here.parents[1] / "workflow" / "main.nf",
    ]
    return next((p for p in candidates if p.is_file()), candidates[0])
