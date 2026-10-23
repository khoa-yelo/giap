"""M7-10 lead decision 1 in the Nextflow wiring (static, no Nextflow needed): the VARIANTS
subworkflow hands the optional C17 ``peptide_spectra`` to the ``variants`` step, and main.nf
passes the preflight-written file to it (``giap run`` puts it in params like the other
normalised preflight files). The key names follow the step registry, as in
tests/unit/test_nf_variants_registry.py.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
SUBWF = REPO / "workflow" / "subworkflows" / "local" / "variants.nf"
MAIN = REPO / "workflow" / "main.nf"


def test_variants_subworkflow_stages_peptide_spectra() -> None:
    text = SUBWF.read_text(encoding="utf-8")
    assert re.search(r"kv\(\s*'peptide_spectra'\s*,", text), "peptide_spectra is not staged"
    need = re.search(r"def needed = \[([^\]]*)\]", text)
    assert need is not None and "peptide_spectra" not in need.group(1), "C17 must stay optional"


def test_main_passes_peptide_spectra_to_variants() -> None:
    text = MAIN.read_text(encoding="utf-8")
    call = re.search(r"VARIANTS\((.*?)\)\s*\n", text, re.S)
    assert call is not None
    assert "peptide_spectra" in call.group(1), call.group(1)
